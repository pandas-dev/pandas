import re
import asyncio
import tokenize
from io import StringIO
from typing import Callable, List, Optional, Union, Generator, Tuple, ClassVar, Any
import warnings

import prompt_toolkit
from prompt_toolkit.buffer import Buffer
from prompt_toolkit.key_binding import KeyPressEvent
from prompt_toolkit.key_binding.bindings import named_commands as nc
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory, Suggestion, AutoSuggest
from prompt_toolkit.document import Document
from prompt_toolkit.history import History
from prompt_toolkit.shortcuts import PromptSession
from prompt_toolkit.layout.processors import (
    Processor,
    Transformation,
    TransformationInput,
)

from IPython.core.getipython import get_ipython
from IPython.utils.tokenutil import generate_tokens

from .filters import pass_through


def _get_query(document: Document):
    return document.lines[document.cursor_position_row]


class AppendAutoSuggestionInAnyLine(Processor):
    """
    Append the auto suggestion to lines other than the last (appending to the
    last line is natively supported by the prompt toolkit).

    This has a private `_debug` attribute that can be set to True to display
    debug information as virtual suggestion on the end of any line. You can do
    so with:

        >>> from IPython.terminal.shortcuts.auto_suggest import AppendAutoSuggestionInAnyLine
        >>> AppendAutoSuggestionInAnyLine._debug = True

    """

    _debug: ClassVar[bool] = False

    def __init__(self, style: str = "class:auto-suggestion") -> None:
        self.style = style

    def apply_transformation(self, ti: TransformationInput) -> Transformation:
        """
         Apply transformation to the line that is currently being edited.

         This is a variation of the original implementation in prompt toolkit
         that allows to not only append suggestions to any line, but also to show
         multi-line suggestions.

         As transformation are applied on a line-by-line basis; we need to trick
         a bit, and elide any line that is after the line we are currently
         editing, until we run out of completions. We cannot shift the existing
         lines

         There are multiple cases to handle:

         The completions ends before the end of the buffer:
             We can resume showing the normal line, and say that some code may
             be hidden.

        The completions ends at the end of the buffer
             We can just say that some code may be hidden.

         And separately:

         The completions ends beyond the end of the buffer
             We need to both say that some code may be hidden, and that some
             lines are not shown.

        """
        last_line_number = ti.document.line_count - 1
        is_last_line = ti.lineno == last_line_number

        noop = lambda text: Transformation(
            fragments=ti.fragments + [(self.style, " " + text if self._debug else "")]
        )
        if ti.document.line_count == 1:
            return noop("noop:oneline")
        if ti.document.cursor_position_row == last_line_number and is_last_line:
            # prompt toolkit already appends something; just leave it be
            return noop("noop:last line and cursor")

        # first everything before the current line is unchanged.
        if ti.lineno < ti.document.cursor_position_row:
            return noop("noop:before cursor")

        buffer = ti.buffer_control.buffer
        if not buffer.suggestion or not ti.document.is_cursor_at_the_end_of_line:
            return noop("noop:not eol")

        delta = ti.lineno - ti.document.cursor_position_row
        suggestions = buffer.suggestion.text.splitlines()

        if len(suggestions) == 0:
            return noop("noop: no suggestions")

        suggestions_longer_than_buffer: bool = (
            len(suggestions) + ti.document.cursor_position_row > ti.document.line_count
        )
        if prompt_toolkit.VERSION < (3, 0, 49):
            if len(suggestions) > 1 and prompt_toolkit.VERSION < (3, 0, 49):
                if ti.lineno == ti.document.cursor_position_row:
                    return Transformation(
                        fragments=ti.fragments
                        + [
                            (
                                "red",
                                "(Cannot show multiline suggestion; requires prompt_toolkit > 3.0.49)",
                            )
                        ]
                    )
                else:
                    return Transformation(fragments=ti.fragments)
            elif len(suggestions) == 1:
                if ti.lineno == ti.document.cursor_position_row:
                    return Transformation(
                        fragments=ti.fragments + [(self.style, suggestions[0])]
                    )
            return Transformation(fragments=ti.fragments)

        if delta == 0:
            suggestion = suggestions[0]
            return Transformation(fragments=ti.fragments + [(self.style, suggestion)])
        if is_last_line:
            if delta < len(suggestions):
                extra = f"; {len(suggestions) - delta} line(s) hidden"
                suggestion = f"… rest of suggestion ({len(suggestions) - delta} lines) and code hidden"
                return Transformation([(self.style, suggestion)])

            n_elided = len(suggestions)
            for i in range(len(suggestions)):
                ll = ti.get_line(last_line_number - i)
                el = "".join(l[1] for l in ll).strip()
                if el:
                    break
                else:
                    n_elided -= 1
            if n_elided:
                return Transformation([(self.style, f"… {n_elided} line(s) hidden")])
            else:
                return Transformation(
                    ti.get_line(last_line_number - len(suggestions) + 1)
                    + ([(self.style, "shift-last-line")] if self._debug else [])
                )

        elif delta < len(suggestions):
            suggestion = suggestions[delta]
            return Transformation([(self.style, suggestion)])
        else:
            shift = ti.lineno - len(suggestions) + 1
            return Transformation(ti.get_line(shift))


class NavigableAutoSuggestFromHistory(AutoSuggestFromHistory):
    """
    A subclass of AutoSuggestFromHistory that allow navigation to next/previous
    suggestion from history. To do so it remembers the current position, but it
    state need to carefully be cleared on the right events.
    """

    skip_lines: int
    _connected_apps: list[PromptSession]

    # handle to the currently running llm task that appends suggestions to the
    # current buffer; we keep a handle to it in order to cancell it when there is a cursor movement, or
    # another request.
    _llm_task: asyncio.Task | None = None

    # This is the instance of the LLM provider from jupyter-ai to which we forward the request
    # to generate inline completions.
    _llm_provider: Any | None
    _llm_prefixer: callable = lambda self, x: "wrong"

    def __init__(self):
        super().__init__()
        self.skip_lines = 0
        self._connected_apps = []
        self._llm_provider = None

    def reset_history_position(self, _: Buffer) -> None:
        self.skip_lines = 0

    def disconnect(self) -> None:
        self._cancel_running_llm_task()
        for pt_app in self._connected_apps:
            text_insert_event = pt_app.default_buffer.on_text_insert
            text_insert_event.remove_handler(self.reset_history_position)

    def connect(self, pt_app: PromptSession) -> None:
        self._connected_apps.append(pt_app)
        # note: `on_text_changed` could be used for a bit different behaviour
        # on character deletion (i.e. resetting history position on backspace)
        pt_app.default_buffer.on_text_insert.add_handler(self.reset_history_position)
        pt_app.default_buffer.on_cursor_position_changed.add_handler(self._dismiss)

    def get_suggestion(
        self, buffer: Buffer, document: Document
    ) -> Optional[Suggestion]:
        text = _get_query(document)

        if text.strip():
            for suggestion, _ in self._find_next_match(
                text, self.skip_lines, buffer.history
            ):
                return Suggestion(suggestion)

        return None

    def _dismiss(self, buffer, *args, **kwargs) -> None:
        self._cancel_running_llm_task()
        buffer.suggestion = None

    def _find_match(
        self, text: str, skip_lines: float, history: History, previous: bool
    ) -> Generator[Tuple[str, float], None, None]:
        """
        text : str
            Text content to find a match for, the user cursor is most of the
            time at the end of this text.
        skip_lines : float
            number of items to skip in the search, this is used to indicate how
            far in the list the user has navigated by pressing up or down.
            The float type is used as the base value is +inf
        history : History
            prompt_toolkit History instance to fetch previous entries from.
        previous : bool
            Direction of the search, whether we are looking previous match
            (True), or next match (False).

        Yields
        ------
        Tuple with:
        str:
            current suggestion.
        float:
            will actually yield only ints, which is passed back via skip_lines,
            which may be a +inf (float)


        """
        line_number = -1
        for string in reversed(list(history.get_strings())):
            for line in reversed(string.splitlines()):
                line_number += 1
                if not previous and line_number < skip_lines:
                    continue
                # do not return empty suggestions as these
                # close the auto-suggestion overlay (and are useless)
                if line.startswith(text) and len(line) > len(text):
                    yield line[len(text) :], line_number
                if previous and line_number >= skip_lines:
                    return

    def _find_next_match(
        self, text: str, skip_lines: float, history: History
    ) -> Generator[Tuple[str, float], None, None]:
        return self._find_match(text, skip_lines, history, previous=False)

    def _find_previous_match(self, text: str, skip_lines: float, history: History):
        return reversed(
            list(self._find_match(text, skip_lines, history, previous=True))
        )

    def up(self, query: str, other_than: str, history: History) -> None:
        self._cancel_running_llm_task()
        for suggestion, line_number in self._find_next_match(
            query, self.skip_lines, history
        ):
            # if user has history ['very.a', 'very', 'very.b'] and typed 'very'
            # we want to switch from 'very.b' to 'very.a' because a) if the
            # suggestion equals current text, prompt-toolkit aborts suggesting
            # b) user likely would not be interested in 'very' anyways (they
            # already typed it).
            if query + suggestion != other_than:
                self.skip_lines = line_number
                break
        else:
            # no matches found, cycle back to beginning
            self.skip_lines = 0

    def down(self, query: str, other_than: str, history: History) -> None:
        self._cancel_running_llm_task()
        for suggestion, line_number in self._find_previous_match(
            query, self.skip_lines, history
        ):
            if query + suggestion != other_than:
                self.skip_lines = line_number
                break
        else:
            # no matches found, cycle to end
            for suggestion, line_number in self._find_previous_match(
                query, float("Inf"), history
            ):
                if query + suggestion != other_than:
                    self.skip_lines = line_number
                    break

    def _cancel_running_llm_task(self) -> None:
        """
        Try to cancell the currently running llm_task if exists, and set it to None.
        """
        if self._llm_task is not None:
            if self._llm_task.done():
                self._llm_task = None
                return
            cancelled = self._llm_task.cancel()
            if cancelled:
                self._llm_task = None
            if not cancelled:
                warnings.warn(
                    "LLM task not cancelled, does your provider support cancellation?"
                )

    async def _trigger_llm(self, buffer) -> None:
        """
        This will ask the current llm provider a suggestion for the current buffer.

        If there is a currently running llm task, it will cancel it.
        """
        # we likely want to store the current cursor position, and cancel if the cursor has moved.
        try:
            import jupyter_ai_magics
            import jupyter_ai.completions.models as jai_models
        except ModuleNotFoundError:
            jai_models = None
        if not self._llm_provider:
            warnings.warn("No LLM provider found, cannot trigger LLM completions")
            return
        if jai_models is None:
            warnings.warn(
                "LLM Completion requires `jupyter_ai_magics` and `jupyter_ai` to be installed"
            )

        self._cancel_running_llm_task()

        async def error_catcher(buffer):
            """
            This catches and log any errors, as otherwise this is just
            lost in the void of the future running task.
            """
            try:
                await self._trigger_llm_core(buffer)
            except Exception as e:
                get_ipython().log.error("error")
                raise

        # here we need a cancellable task so we can't just await the error catched
        self._llm_task = asyncio.create_task(error_catcher(buffer))
        await self._llm_task

    async def _trigger_llm_core(self, buffer: Buffer):
        """
        This is the core of the current llm request.

        Here we build a compatible `InlineCompletionRequest` and ask the llm
        provider to stream it's response back to us iteratively setting it as
        the suggestion on the current buffer.

        Unlike with JupyterAi, as we do not have multiple cell, the cell number
        is always set to `0`, note that we _could_ set it to a new number each
        time and ignore threply from past numbers.

        We set the prefix to the current cell content, but could also inset the
        rest of the history or even just the non-fail history.

        In the same way, we do not have cell id.

        LLM provider may return multiple suggestion stream, but for the time
        being we only support one.

        Here we make the assumption that the provider will have
        stream_inline_completions, I'm not sure it is the case for all
        providers.
        """
        try:
            import jupyter_ai_magics
            import jupyter_ai.completions.models as jai_models
        except ModuleNotFoundError:
            jai_models = None

        hm = buffer.history.shell.history_manager
        prefix = self._llm_prefixer(hm)
        print(prefix)

        request = jai_models.InlineCompletionRequest(
            number=0,
            prefix=prefix + buffer.document.text,
            suffix="",
            mime="text/x-python",
            stream=True,
            path=None,
            language="python",
            cell_id=None,
        )

        async for reply_and_chunks in self._llm_provider.stream_inline_completions(
            request
        ):
            if isinstance(reply_and_chunks, jai_models.InlineCompletionReply):
                if len(reply_and_chunks.list.items) > 1:
                    raise ValueError(
                        "Terminal IPython cannot deal with multiple LLM suggestions at once"
                    )
                buffer.suggestion = Suggestion(
                    reply_and_chunks.list.items[0].insertText
                )
                buffer.on_suggestion_set.fire()
            elif isinstance(reply_and_chunks, jai_models.InlineCompletionStreamChunk):
                buffer.suggestion = Suggestion(reply_and_chunks.response.insertText)
                buffer.on_suggestion_set.fire()
        return


async def llm_autosuggestion(event: KeyPressEvent):
    """
    Ask the AutoSuggester from history to delegate to ask an LLM for completion

    This will first make sure that the current buffer have _MIN_LINES (7)
    available lines to insert the LLM completion

    Provisional as of 8.32, may change without warnigns

    """
    _MIN_LINES = 5
    provider = get_ipython().auto_suggest
    if not isinstance(provider, NavigableAutoSuggestFromHistory):
        return
    doc = event.current_buffer.document
    lines_to_insert = max(0, _MIN_LINES - doc.line_count + doc.cursor_position_row)
    for _ in range(lines_to_insert):
        event.current_buffer.insert_text("\n", move_cursor=False)

    await provider._trigger_llm(event.current_buffer)


def accept_or_jump_to_end(event: KeyPressEvent):
    """Apply autosuggestion or jump to end of line."""
    buffer = event.current_buffer
    d = buffer.document
    after_cursor = d.text[d.cursor_position :]
    lines = after_cursor.split("\n")
    end_of_current_line = lines[0].strip()
    suggestion = buffer.suggestion
    if (suggestion is not None) and (suggestion.text) and (end_of_current_line == ""):
        buffer.insert_text(suggestion.text)
    else:
        nc.end_of_line(event)


def accept(event: KeyPressEvent):
    """Accept autosuggestion"""
    buffer = event.current_buffer
    suggestion = buffer.suggestion
    if suggestion:
        buffer.insert_text(suggestion.text)
    else:
        nc.forward_char(event)


def discard(event: KeyPressEvent):
    """Discard autosuggestion"""
    buffer = event.current_buffer
    buffer.suggestion = None


def accept_word(event: KeyPressEvent):
    """Fill partial autosuggestion by word"""
    buffer = event.current_buffer
    suggestion = buffer.suggestion
    if suggestion:
        t = re.split(r"(\S+\s+)", suggestion.text)
        buffer.insert_text(next((x for x in t if x), ""))
    else:
        nc.forward_word(event)


def accept_character(event: KeyPressEvent):
    """Fill partial autosuggestion by character"""
    b = event.current_buffer
    suggestion = b.suggestion
    if suggestion and suggestion.text:
        b.insert_text(suggestion.text[0])


def accept_and_keep_cursor(event: KeyPressEvent):
    """Accept autosuggestion and keep cursor in place"""
    buffer = event.current_buffer
    old_position = buffer.cursor_position
    suggestion = buffer.suggestion
    if suggestion:
        buffer.insert_text(suggestion.text)
        buffer.cursor_position = old_position


def accept_and_move_cursor_left(event: KeyPressEvent):
    """Accept autosuggestion and move cursor left in place"""
    accept_and_keep_cursor(event)
    nc.backward_char(event)


def _update_hint(buffer: Buffer):
    if buffer.auto_suggest:
        suggestion = buffer.auto_suggest.get_suggestion(buffer, buffer.document)
        buffer.suggestion = suggestion


def backspace_and_resume_hint(event: KeyPressEvent):
    """Resume autosuggestions after deleting last character"""
    nc.backward_delete_char(event)
    _update_hint(event.current_buffer)


def resume_hinting(event: KeyPressEvent):
    """Resume autosuggestions"""
    pass_through.reply(event)
    # Order matters: if update happened first and event reply second, the
    # suggestion would be auto-accepted if both actions are bound to same key.
    _update_hint(event.current_buffer)


def up_and_update_hint(event: KeyPressEvent):
    """Go up and update hint"""
    current_buffer = event.current_buffer

    current_buffer.auto_up(count=event.arg)
    _update_hint(current_buffer)


def down_and_update_hint(event: KeyPressEvent):
    """Go down and update hint"""
    current_buffer = event.current_buffer

    current_buffer.auto_down(count=event.arg)
    _update_hint(current_buffer)


def accept_token(event: KeyPressEvent):
    """Fill partial autosuggestion by token"""
    b = event.current_buffer
    suggestion = b.suggestion

    if suggestion:
        prefix = _get_query(b.document)
        text = prefix + suggestion.text

        tokens: List[Optional[str]] = [None, None, None]
        substrings = [""]
        i = 0

        for token in generate_tokens(StringIO(text).readline):
            if token.type == tokenize.NEWLINE:
                index = len(text)
            else:
                index = text.index(token[1], len(substrings[-1]))
            substrings.append(text[:index])
            tokenized_so_far = substrings[-1]
            if tokenized_so_far.startswith(prefix):
                if i == 0 and len(tokenized_so_far) > len(prefix):
                    tokens[0] = tokenized_so_far[len(prefix) :]
                    substrings.append(tokenized_so_far)
                    i += 1
                tokens[i] = token[1]
                if i == 2:
                    break
                i += 1

        if tokens[0]:
            to_insert: str
            insert_text = substrings[-2]
            if tokens[1] and len(tokens[1]) == 1:
                insert_text = substrings[-1]
            to_insert = insert_text[len(prefix) :]
            b.insert_text(to_insert)
            return

    nc.forward_word(event)


Provider = Union[AutoSuggestFromHistory, NavigableAutoSuggestFromHistory, None]


def _swap_autosuggestion(
    buffer: Buffer,
    provider: NavigableAutoSuggestFromHistory,
    direction_method: Callable,
):
    """
    We skip most recent history entry (in either direction) if it equals the
    current autosuggestion because if user cycles when auto-suggestion is shown
    they most likely want something else than what was suggested (otherwise
    they would have accepted the suggestion).
    """
    suggestion = buffer.suggestion
    if not suggestion:
        return

    query = _get_query(buffer.document)
    current = query + suggestion.text

    direction_method(query=query, other_than=current, history=buffer.history)

    new_suggestion = provider.get_suggestion(buffer, buffer.document)
    buffer.suggestion = new_suggestion


def swap_autosuggestion_up(event: KeyPressEvent):
    """Get next autosuggestion from history."""
    shell = get_ipython()
    provider = shell.auto_suggest

    if not isinstance(provider, NavigableAutoSuggestFromHistory):
        return

    return _swap_autosuggestion(
        buffer=event.current_buffer, provider=provider, direction_method=provider.up
    )


def swap_autosuggestion_down(event: KeyPressEvent):
    """Get previous autosuggestion from history."""
    shell = get_ipython()
    provider = shell.auto_suggest

    if not isinstance(provider, NavigableAutoSuggestFromHistory):
        return

    return _swap_autosuggestion(
        buffer=event.current_buffer,
        provider=provider,
        direction_method=provider.down,
    )
