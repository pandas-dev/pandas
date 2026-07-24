"""Textual application for b2view."""

from __future__ import annotations

import contextlib
import io
import os
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np
from rich.markup import escape as markup_escape
from textual import work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.content import Content
from textual.css.query import NoMatches
from textual.screen import ModalScreen
from textual.theme import Theme
from textual.widgets import (
    Checkbox,
    DataTable,
    Footer,
    Header,
    Input,
    OptionList,
    ProgressBar,
    SelectionList,
    Static,
    Tree,
)
from textual.widgets._header import HeaderTitle
from textual.widgets.option_list import Option
from textual.widgets.selection_list import Selection

try:
    from textual_plotext import PlotextPlot
except ImportError:  # plotting is optional
    PlotextPlot = None

try:
    # Auto-selects the best terminal image protocol (kitty/iTerm2/sixel),
    # degrading to colored half-cells; used by the high-res 'h' plot view.
    from textual_image.widget import Image as TextualImage
except ImportError:  # high-res view is optional
    TextualImage = None

import blosc2
from blosc2.b2view.model import DataSliceLayout, StoreBrowser
from blosc2.b2view.render import (
    column_float_decimals,
    format_cell,
    make_metadata_renderable,
    make_preview_renderables,
)

if TYPE_CHECKING:
    from textual import events

_KIND_ICONS = {
    "group": "📁",
    "ndarray": "▦",
    "c2array": "▦",
    "ctable": "▤",
    "schunk": "▣",
    "unknown": "?",
}

# Source kinds whose data grid supports horizontal (column) paging.
_COL_PAGED_KINDS = frozenset({"ndarray2d", "ndarray_slice", "ctable"})

# Blosc2-branded palette layered over Textual's default dark canvas: only the
# logo colors are overridden (background/surface/panel stay None so they derive
# the same near-black as textual-dark).  Turquoise is used for all borders and
# scrollbars (deep blue is too dark to read on the dark canvas), with yellow as
# the accent for the focused pane's border.
BLOSC2_THEME = Theme(
    name="blosc2",
    primary="#007a86",  # turquoise
    secondary="#007a86",  # turquoise (deep blue reads poorly on a dark canvas)
    accent="#df9e00",  # yellow
    foreground="#e0e0e0",  # match textual-dark's foreground
    dark=True,
)


def _accent_chip(text: str) -> str:
    """A reverse-video status chip in the brand accent (dark text on yellow)."""
    return f"[$background on $accent] {text} [/]"


class B2ViewPanel(Vertical):
    """Pane container that can be maximized."""

    ALLOW_MAXIMIZE = True


class BufferedDataTable(DataTable):
    """DataTable with app-controlled page changes at row boundaries."""

    def action_cursor_down(self) -> None:
        app = self.app
        if getattr(app, "_dim_mode", False):
            getattr(app, "_dim_adjust", lambda _: None)(-1)
            return
        if self.cursor_row >= self.row_count - 1 and getattr(app, "page_table", lambda _: False)(1):
            return
        super().action_cursor_down()

    def action_cursor_up(self) -> None:
        app = self.app
        if getattr(app, "_dim_mode", False):
            getattr(app, "_dim_adjust", lambda _: None)(1)
            return
        if self.cursor_row <= 0 and getattr(app, "page_table", lambda _: False)(-1):
            return
        super().action_cursor_up()

    def action_cursor_right(self) -> None:
        app = self.app
        if getattr(app, "_dim_mode", False):
            getattr(app, "_dim_cursor", lambda _: None)(1)
            return
        if self.cursor_column >= len(self.columns) - 1 and getattr(
            app, "page_grid_columns", lambda _: False
        )(1):
            return
        super().action_cursor_right()

    def action_cursor_left(self) -> None:
        app = self.app
        if getattr(app, "_dim_mode", False):
            getattr(app, "_dim_cursor", lambda _: None)(-1)
            return
        if self.cursor_column <= 0 and getattr(app, "page_grid_columns", lambda _: False)(-1):
            return
        super().action_cursor_left()

    def action_page_down(self) -> None:
        if getattr(self.app, "page_table", lambda *a, **k: False)(1, align=True):
            return
        super().action_page_down()

    def action_page_up(self) -> None:
        if getattr(self.app, "page_table", lambda *a, **k: False)(-1, align=True):
            return
        super().action_page_up()

    def action_page_right(self) -> None:
        if getattr(self.app, "page_grid_columns", lambda _: False)(1):
            return
        super().action_page_right()

    def action_page_left(self) -> None:
        if getattr(self.app, "page_grid_columns", lambda _: False)(-1):
            return
        super().action_page_left()

    def action_select_cursor(self) -> None:
        app = self.app
        if getattr(app, "_dim_mode", False):
            getattr(app, "action_dim_toggle_nav", lambda: None)()
            return
        if getattr(app, "_drilldown_arg_cell", lambda: False)():
            return
        if getattr(app, "_inspect_cursor_cell", lambda: False)():
            return
        super().action_select_cursor()

    def _wheel_step(self) -> int:
        # Half the visible rows per tick; arrow keys remain the
        # single-step path (also for dim-mode index changes).
        return max(1, self.row_count // 2)

    def on_mouse_scroll_down(self, event: events.MouseScrollDown) -> None:
        # The grid holds exactly one viewport-sized page, so the default
        # scroll handler has nothing to scroll; move the cursor instead,
        # which pages at the edges just like the arrow keys.
        event.stop()
        event.prevent_default()
        for _ in range(self._wheel_step()):
            self.action_cursor_down()

    def on_mouse_scroll_up(self, event: events.MouseScrollUp) -> None:
        event.stop()
        event.prevent_default()
        for _ in range(self._wheel_step()):
            self.action_cursor_up()

    def on_resize(self, event) -> None:
        # The column/row windows are fitted to this table's size; re-check
        # whenever it changes (terminal resize, panel maximize, ...).
        getattr(self.app, "_on_data_table_resized", lambda: None)()

    def action_scroll_home(self) -> None:
        if getattr(self.app, "_grid_col_home", lambda: False)():
            pass
        else:
            super().action_scroll_home()

    def action_scroll_end(self) -> None:
        if getattr(self.app, "_grid_col_end", lambda: False)():
            pass
        else:
            super().action_scroll_end()


class HelpScreen(ModalScreen[None]):
    """Modal listing all key bindings, grouped by area."""

    CSS = """
    HelpScreen {
        align: center middle;
    }
    #help-dialog {
        width: 62;
        height: auto;
        max-height: 90%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #help-title {
        text-style: bold;
        margin-bottom: 1;
    }
    #help-body {
        height: auto;
    }
    """

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
    ]

    _SECTIONS: ClassVar = [
        (
            "Panels",
            [
                ("tab / shift+tab", "next / previous panel"),
                ("m", "maximize the focused panel"),
                ("r", "restore panel (or refresh the tree)"),
                ("q", "quit"),
            ],
        ),
        (
            "Tree",
            [
                ("up / down", "move between nodes"),
                ("enter", "select node (and expand groups)"),
            ],
        ),
        (
            "Data grid — rows",
            [
                ("up / down", "move cursor; pages at the edges"),
                ("pageup / pagedown", "previous / next page"),
                ("t / b", "first / last row"),
                ("g", "go to row..."),
                ("f", "filter rows (CTable)"),
                ("S", "sort by an indexed column, or the grouped result (CTable; R reverses)"),
                ("R", "reverse the current sort order (when sorted)"),
                ("G", "group by a dictionary/numeric column (CTable; p shows a bar chart)"),
                ("escape", "unlock a row window / clear the active filter, sort or group"),
            ],
        ),
        (
            "Data grid — columns",
            [
                ("left / right", "move cursor; pages at the edges"),
                ("s / e  (home / end)", "first / last column window"),
                ("c", "go to column (searchable name list; CTable, else index)"),
                ("/", "pick which columns to show (searchable multi-select; CTable)"),
                ("p", "plot a whole-column overview (needs textual-plotext)"),
                ("enter", "decode a skipped cell; or jump to an argmin/argmax row (grouped)"),
            ],
        ),
        (
            "Plot modal (after 'p')",
            [
                ("+ / -", "zoom in / out about the left edge"),
                ("left / right", "pan the zoomed window"),
                ("0", "reset to the whole series"),
                ("g", "type an exact start:stop row range"),
                ("v", "lock the data grid to the current range (esc unlocks)"),
                ("h", "high-res matplotlib image of the current range"),
                ("escape", "close the plot (q quits b2view)"),
            ],
        ),
        (
            "Dim mode (N-D arrays)",
            [
                ("d", "toggle dim mode"),
                ("left / right", "select the active dimension"),
                ("up / down", "change fixed index / scroll viewport"),
                ("enter", "toggle fixed <-> navigable"),
                ("escape", "exit dim mode"),
            ],
        ),
    ]

    def compose(self) -> ComposeResult:
        from rich.table import Table

        body = Table(show_header=False, box=None, padding=(0, 1))
        body.add_column("key", style="bold cyan", no_wrap=True)
        body.add_column("action")
        for i, (section, entries) in enumerate(self._SECTIONS):
            if i:
                body.add_row("", "")
            body.add_row(f"[bold]{section}[/bold]", "")
            for key, action in entries:
                body.add_row(key, action)
        with Vertical(id="help-dialog"):
            yield Static("b2view keys  (esc to close)", id="help-title")
            with VerticalScroll(id="help-body"):
                yield Static(body)

    def action_close(self) -> None:
        self.dismiss(None)


class GoToRowScreen(ModalScreen[int | None]):
    """Small modal asking for a global row number."""

    CSS = """
    GoToRowScreen {
        align: center middle;
    }
    #goto-dialog {
        width: 50;
        height: auto;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #goto-title {
        text-style: bold;
        margin-bottom: 1;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(self, *, nrows: int, current: int):
        super().__init__()
        self.nrows = nrows
        self.current = current

    def compose(self) -> ComposeResult:
        with Vertical(id="goto-dialog"):
            yield Static(f"Go to row 0..{self.nrows - 1} (current: {self.current})", id="goto-title")
            yield Input(placeholder="row number", id="goto-input")

    def on_mount(self) -> None:
        input_widget = self.query_one("#goto-input", Input)
        input_widget.value = str(self.current)
        input_widget.focus()
        # Pre-select the current value so the first keystroke replaces it (typing
        # a fresh number is the common case); arrows/edits still work as usual.
        input_widget.select_all()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        value = event.value.strip().replace("_", "")
        try:
            row = int(value)
        except ValueError:
            self.query_one("#goto-title", Static).update("Please enter an integer row number")
            return
        if not 0 <= row < self.nrows:
            self.query_one("#goto-title", Static).update(f"Row must be in range 0..{self.nrows - 1}")
            return
        self.dismiss(row)

    def action_cancel(self) -> None:
        self.dismiss(None)


class GoToColumnScreen(ModalScreen[int | None]):
    """Small modal asking for a column index or (for CTables) a column name."""

    CSS = """
    GoToColumnScreen {
        align: center middle;
    }
    #gotocol-dialog {
        width: 50;
        height: auto;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #gotocol-title {
        text-style: bold;
        margin-bottom: 1;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(self, *, ncols: int, current: int, names: list[str] | None = None):
        super().__init__()
        self.ncols = ncols
        self.current = current
        self.names = names

    def compose(self) -> ComposeResult:
        what = f"column 0..{self.ncols - 1}"
        if self.names:
            what += " or name"
        with Vertical(id="gotocol-dialog"):
            yield Static(f"Go to {what} (current: {self.current})", id="gotocol-title")
            yield Input(placeholder="column index or name", id="gotocol-input")

    def on_mount(self) -> None:
        input_widget = self.query_one("#gotocol-input", Input)
        input_widget.value = str(self.current)
        input_widget.focus()
        # Pre-select the current index so typing a column name (or a new index)
        # replaces it instead of appending (e.g. "0" + "payment.fare").
        input_widget.select_all()

    def _fail(self, message: str) -> None:
        self.query_one("#gotocol-title", Static).update(message)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        value = event.value.strip().replace("_", "")
        try:
            col = int(value)
        except ValueError:
            col = self._match_name(event.value.strip())
            if col is None:
                return
        if not 0 <= col < self.ncols:
            self._fail(f"Column must be in range 0..{self.ncols - 1}")
            return
        self.dismiss(col)

    def _match_name(self, value: str) -> int | None:
        """Resolve a column name (exact, or unique prefix) to its index."""
        if not self.names:
            self._fail("Please enter an integer column index")
            return None
        if value in self.names:
            return self.names.index(value)
        matches = [i for i, name in enumerate(self.names) if name.startswith(value)] if value else []
        if len(matches) == 1:
            return matches[0]
        self._fail(f"{'Ambiguous' if matches else 'Unknown'} column name {value!r}")
        return None

    def action_cancel(self) -> None:
        self.dismiss(None)


class ColumnSelectScreen(ModalScreen["int | None"]):
    """Searchable column picker: type to filter, ↑/↓ to move, Enter to choose.

    Dismisses with the chosen column's index into *names* (or None on cancel),
    a drop-in for :class:`GoToColumnScreen`'s result contract.  Used by the
    ``c`` go-to-column key (CTables, where columns have names) and by the
    scatter ``s`` key to pick the Y column from the visible-column universe.
    """

    CSS = """
    ColumnSelectScreen {
        align: center middle;
    }
    #colselect-dialog {
        width: 50;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #colselect-title {
        text-style: bold;
        margin-bottom: 1;
    }
    #colselect-list {
        height: auto;
        max-height: 16;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(self, *, names: list[str], title: str = "Select column"):
        super().__init__()
        self.names = names
        self._title = title

    def compose(self) -> ComposeResult:
        with Vertical(id="colselect-dialog"):
            yield Static(self._title, id="colselect-title")
            yield Input(placeholder="type to filter…", id="colselect-input")
            yield OptionList(id="colselect-list")

    def on_mount(self) -> None:
        self._populate("")
        self.query_one("#colselect-input", Input).focus()

    def _populate(self, query: str) -> None:
        """Refill the list with names matching *query* (case-insensitive substring).

        Each option's id is the column's original index into *names*, so a match
        resolves to the right column regardless of the current filtering.
        """
        q = query.strip().lower()
        option_list = self.query_one("#colselect-list", OptionList)
        option_list.clear_options()
        matches = [Option(name, id=str(i)) for i, name in enumerate(self.names) if q in name.lower()]
        option_list.add_options(matches)
        if matches:
            option_list.highlighted = 0

    def on_input_changed(self, event: Input.Changed) -> None:
        self._populate(event.value)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        # Enter from the filter input accepts the currently highlighted match.
        self._accept_highlighted()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        if event.option.id is not None:
            self.dismiss(int(event.option.id))

    def _accept_highlighted(self) -> None:
        option_list = self.query_one("#colselect-list", OptionList)
        idx = option_list.highlighted
        if idx is None:
            return
        option = option_list.get_option_at_index(idx)
        if option.id is not None:
            self.dismiss(int(option.id))

    def on_key(self, event: events.Key) -> None:
        # Focus stays in the filter Input; ↑/↓ drive the list highlight so the
        # user can type then arrow to a match without leaving the keyboard home.
        if event.key in ("down", "up"):
            option_list = self.query_one("#colselect-list", OptionList)
            (option_list.action_cursor_down if event.key == "down" else option_list.action_cursor_up)()
            event.stop()

    def action_cancel(self) -> None:
        self.dismiss(None)


class _ApplySelectionList(SelectionList):
    """A SelectionList whose Enter *applies* (dismisses) instead of toggling.

    SelectionList toggles the highlighted row on Space; Enter is inherited from
    OptionList as another toggle.  Here Enter is rebound so it bubbles up as
    "apply the chosen set", matching the Enter-applies convention of the other
    b2view modals (Space still toggles).
    """

    BINDINGS: ClassVar = [Binding("enter", "apply", "Apply", show=False)]

    def action_apply(self) -> None:
        self.screen.action_apply_filter()


class ColumnFilterScreen(ModalScreen["list[str] | None"]):
    """Searchable multi-select for which CTable columns to show.

    Type in the filter box to narrow the candidate list; Tab (or ↓) moves into
    the checkbox list where Space toggles a column; Enter applies the checked
    set and Escape cancels.  Opens with the currently-visible columns checked;
    applying an empty set (or all of them) shows every column.

    Dismisses with the chosen column names in table order, or ``None`` on cancel.
    """

    CSS = """
    ColumnFilterScreen {
        align: center middle;
    }
    #colfilter-dialog {
        width: 60;
        height: auto;
        max-height: 90%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #colfilter-title {
        text-style: bold;
        margin-bottom: 1;
    }
    #colfilter-list {
        height: auto;
        max-height: 18;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(self, *, names: list[str], selected: list[str]):
        super().__init__()
        self.names = names
        self._checked: set[str] = {n for n in selected if n in set(names)}
        self._visible: list[str] = list(names)
        self._populating = False

    def compose(self) -> ComposeResult:
        with Vertical(id="colfilter-dialog"):
            yield Static(
                "Show columns — type to filter · Tab/↓ to list · Space toggles · Enter applies",
                id="colfilter-title",
            )
            yield Input(placeholder="type to filter…", id="colfilter-input")
            yield _ApplySelectionList(id="colfilter-list")

    def on_mount(self) -> None:
        self._populate("")
        self.query_one("#colfilter-input", Input).focus()

    def _populate(self, query: str) -> None:
        """Refill the checkbox list with names matching *query*, preserving checks.

        ``_checked`` is the source of truth across re-filters (a checked column
        that scrolls out of the filtered view stays checked); each option's box
        is seeded from it.
        """
        q = query.strip().lower()
        sel = self.query_one("#colfilter-list", SelectionList)
        self._visible = [name for name in self.names if q in name.lower()]
        self._populating = True
        try:
            sel.clear_options()
            sel.add_options([Selection(name, name, name in self._checked) for name in self._visible])
            if self._visible:
                sel.highlighted = 0
        finally:
            self._populating = False

    def on_input_changed(self, event: Input.Changed) -> None:
        self._populate(event.value)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.action_apply_filter()

    def on_selection_list_selected_changed(self, event: SelectionList.SelectedChanged) -> None:
        # Merge the visible options' checkbox states into _checked, leaving any
        # checked-but-filtered-out columns untouched.  Skipped while _populate
        # is rebuilding the list (those toggles are not user actions).
        if self._populating:
            return
        selected = set(event.selection_list.selected)
        self._checked = (self._checked - set(self._visible)) | selected

    def on_key(self, event: events.Key) -> None:
        # ↓ from the filter box drops focus into the checkbox list (Tab also works).
        if event.key == "down" and self.query_one("#colfilter-input", Input).has_focus:
            self.query_one("#colfilter-list", SelectionList).focus()
            event.stop()

    def action_apply_filter(self) -> None:
        self.dismiss([name for name in self.names if name in self._checked])

    def action_cancel(self) -> None:
        self.dismiss(None)


class FilterScreen(ModalScreen[str | None]):
    """Small modal asking for a CTable filter (row expression or column pattern)."""

    CSS = """
    FilterScreen {
        align: center middle;
    }
    #filter-dialog {
        width: 70;
        height: auto;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #filter-title {
        text-style: bold;
        margin-bottom: 1;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(
        self,
        *,
        current: str | None = None,
        title: str = "Filter rows (empty clears)",
        placeholder: str = "e.g. payment.tips > 100 and trip.km > 0",
    ):
        super().__init__()
        self.current = current or ""
        self.title_text = title
        self.placeholder = placeholder

    def compose(self) -> ComposeResult:
        with Vertical(id="filter-dialog"):
            yield Static(self.title_text, id="filter-title")
            yield Input(placeholder=self.placeholder, id="filter-input")

    def on_mount(self) -> None:
        input_widget = self.query_one("#filter-input", Input)
        input_widget.value = self.current
        input_widget.focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.dismiss(event.value.strip())

    def action_cancel(self) -> None:
        self.dismiss(None)


class SortByScreen(ModalScreen["tuple[str, bool] | None"]):
    """Dropdown to sort a CTable by one of its columns.

    ↑/↓ to pick a column, ``r`` (or click) toggles reverse/descending, Enter
    applies.  ``labels`` may decorate the displayed names (e.g. mark indexed
    columns) while ``columns`` carries the real names returned on selection.
    Dismisses with ``(column, reverse)`` or None on cancel.
    """

    CSS = """
    SortByScreen {
        align: center middle;
    }
    #sortby-dialog {
        width: 60;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #sortby-title {
        text-style: bold;
        margin-bottom: 1;
    }
    #sortby-list {
        height: auto;
        max-height: 16;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel"), ("R", "toggle_reverse", "Reverse")]

    def __init__(
        self,
        *,
        columns: list[str],
        labels: list[str] | None = None,
        current: tuple[str, bool] | None = None,
        title: str = "Sort by column (Enter applies, R reverses)",
    ):
        super().__init__()
        self.columns = columns
        self._labels = labels or columns
        self._current = current
        self._title = title

    def compose(self) -> ComposeResult:
        cur_col, cur_rev = self._current or (None, False)
        with Vertical(id="sortby-dialog"):
            yield Static(self._title, id="sortby-title")
            yield OptionList(
                *(Option(name, id=str(i)) for i, name in enumerate(self._labels)), id="sortby-list"
            )
            yield Checkbox("Reverse (descending)", value=cur_rev, id="sortby-reverse")

    def on_mount(self) -> None:
        option_list = self.query_one("#sortby-list", OptionList)
        cur_col = (self._current or (None, False))[0]
        option_list.highlighted = self.columns.index(cur_col) if cur_col in self.columns else 0
        option_list.focus()

    def action_toggle_reverse(self) -> None:
        checkbox = self.query_one("#sortby-reverse", Checkbox)
        checkbox.value = not checkbox.value

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        if event.option.id is not None:
            reverse = self.query_one("#sortby-reverse", Checkbox).value
            self.dismiss((self.columns[int(event.option.id)], reverse))

    def action_cancel(self) -> None:
        self.dismiss(None)


class GroupByScreen(ModalScreen["tuple[str, str, str | None] | None"]):
    """Group a CTable by a key column and one aggregation.

    Three lists, Enter advances then applies: pick the key, the operation, and —
    for every operation but ``count rows`` — the value column.  ``count rows``
    (the keyless ``size`` aggregation) applies straight from the operation list.
    Dismisses with ``(key, op, value_col|None)`` or None on cancel.
    """

    CSS = """
    GroupByScreen {
        align: center middle;
    }
    #groupby-dialog {
        width: 80;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #groupby-title {
        text-style: bold;
        margin-bottom: 1;
    }
    .groupby-label {
        color: $text-muted;
    }
    #groupby-cols {
        height: auto;
    }
    .groupby-col {
        width: 1fr;
        height: auto;
    }
    #groupby-col-left {
        margin-right: 2;
    }
    #groupby-key, #groupby-op {
        height: auto;
        max-height: 10;
    }
    #groupby-value {
        height: auto;
        max-height: 20;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    # (label, library op); "count rows" is size (no value column needed).
    _ALL_OPS: ClassVar = [
        ("count rows", "size"),
        ("count (non-null)", "count"),
        ("sum", "sum"),
        ("mean", "mean"),
        ("min", "min"),
        ("max", "max"),
        ("argmin", "argmin"),
        ("argmax", "argmax"),
    ]

    def __init__(
        self,
        *,
        keys: list[str],
        values: list[str],
        current: tuple[str, str, str | None] | None = None,
    ):
        super().__init__()
        self.keys = keys
        self.values = values
        # Every op but "count rows" needs a value column; offer only it if none.
        self.ops = self._ALL_OPS if values else self._ALL_OPS[:1]
        self._current = current

    def compose(self) -> ComposeResult:
        with Vertical(id="groupby-dialog"):
            yield Static("Group by (Enter advances / applies)", id="groupby-title")
            with Horizontal(id="groupby-cols"):
                with Vertical(id="groupby-col-left", classes="groupby-col"):
                    yield Static("Key column", classes="groupby-label")
                    yield OptionList(
                        *(Option(name, id=str(i)) for i, name in enumerate(self.keys)),
                        id="groupby-key",
                    )
                    yield Static("Operation", classes="groupby-label")
                    yield OptionList(
                        *(Option(label, id=str(i)) for i, (label, _op) in enumerate(self.ops)),
                        id="groupby-op",
                    )
                with Vertical(classes="groupby-col"):
                    yield Static("Value column", classes="groupby-label")
                    yield OptionList(
                        *(Option(name, id=str(i)) for i, name in enumerate(self.values)),
                        id="groupby-value",
                    )

    def on_mount(self) -> None:
        key_list = self.query_one("#groupby-key", OptionList)
        op_list = self.query_one("#groupby-op", OptionList)
        value_list = self.query_one("#groupby-value", OptionList)
        cur_key, cur_op, cur_val = self._current or (None, None, None)
        key_list.highlighted = self.keys.index(cur_key) if cur_key in self.keys else 0
        op_labels = [op for _label, op in self.ops]
        op_idx = op_labels.index(cur_op) if cur_op in op_labels else 0
        op_list.highlighted = op_idx
        value_list.highlighted = self.values.index(cur_val) if cur_val in self.values else 0
        key_list.focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        # Enter advances key -> operation -> value; "count rows" applies straight
        # from the operation list (no value column needed).
        list_id = event.option_list.id
        if list_id == "groupby-key":
            self.query_one("#groupby-op", OptionList).focus()
            return
        if list_id == "groupby-op":
            op = self.ops[int(event.option.id)][1]
            if op == "size":
                self._apply(op, None)
            else:
                self.query_one("#groupby-value", OptionList).focus()
            return
        op = self.ops[self.query_one("#groupby-op", OptionList).highlighted or 0][1]
        self._apply(op, self.values[int(event.option.id)])

    def _apply(self, op: str, value_col: str | None) -> None:
        key_idx = self.query_one("#groupby-key", OptionList).highlighted
        if key_idx is None:
            return
        self.dismiss((self.keys[key_idx], op, value_col))

    def action_cancel(self) -> None:
        self.dismiss(None)


class GroupBarScreen(ModalScreen[None]):
    """Plot of a grouped view: bars for a categorical key (capped to the top
    groups), or a stem/impulse plot over all groups for a numeric key (discrete
    groups, so no connecting line between adjacent key values)."""

    CSS = """
    GroupBarScreen {
        align: center middle;
    }
    #groupbar-dialog {
        width: 90%;
        height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #groupbar-title {
        text-style: bold;
        height: 1;
    }
    #groupbar-widget {
        height: 1fr;
    }
    #groupbar-keys {
        height: 1;
        color: $text-muted;
    }
    """

    _KEYS_HINT = "h hi-res · esc close"

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
        ("h", "hires", "High-res"),
    ]

    def __init__(self, *, title_prefix: str, bars: dict):
        super().__init__()
        self.title_prefix = title_prefix
        self.numeric = bars.get("numeric", False)
        self.labels = [str(label) for label in bars.get("labels", [])]
        self.x = list(bars.get("x", []))
        self.values = list(bars.get("values", []))
        self.key = bars.get("key", "")
        self.agg = bars.get("agg", "")
        self.xlabel = bars.get("xlabel", self.key)
        total = bars.get("total", len(self.values))
        shown = len(self.values)
        cap = f" · top {shown} of {total} groups" if shown < total else f" · {total} groups"
        self.plot_title = f"{title_prefix}{cap}"

    def compose(self) -> ComposeResult:
        with Vertical(id="groupbar-dialog"):
            yield Static(markup_escape(self.plot_title), id="groupbar-title")
            yield PlotextPlot(id="groupbar-widget")
            yield Static(self._KEYS_HINT, id="groupbar-keys")

    def on_mount(self) -> None:
        widget = self.query_one(PlotextPlot)
        plt = widget.plt
        plt.clear_figure()
        if self.numeric:
            # Stem/impulse, not a connected line: groups are discrete points, so
            # bars-to-baseline show each group's magnitude without inventing a
            # curve between adjacent key values (which, on spiky/quantised data,
            # reads as a phantom baseline).
            if self.values:
                plt.bar(self.x, self.values)
                plt.xlabel(self.xlabel)
                plt.ylabel(self.agg)
        elif self.labels:
            plt.bar(self.labels, self.values)
        widget.refresh()

    def action_hires(self) -> None:
        """h key — open a high-res matplotlib plot over the plotext braille one."""
        if TextualImage is None or not _matplotlib_available():
            self.app.notify(
                "High-res view needs the 'textual-image' and 'matplotlib' packages",
                severity="warning",
            )
            return
        if not self.values:
            self.app.notify("No groups to plot", severity="warning")
            return
        if self.numeric:
            screen = HiResPlotScreen(
                mode="stem",
                title=self.plot_title,
                stem_data=(self.x, self.values),
                xlabel=self.xlabel,
                ylabel=self.agg,
            )
        else:
            screen = HiResPlotScreen(
                mode="bar",
                title=self.plot_title,
                bar_data=(self.labels, self.values),
                xlabel=self.key,
                ylabel=self.agg,
            )
        self.app.push_screen(screen)

    def action_close(self) -> None:
        self.dismiss(None)


def _plot_view(series: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Turn a ``plot_series`` result into drawable arrays + a method label.

    Drops all-NaN buckets (no finite extremes) and maps the read method to a
    human description shown in the title.
    """
    x = np.asarray(series["x"])
    ymin = np.asarray(series["ymin"], dtype=np.float64)
    ymax = np.asarray(series["ymax"], dtype=np.float64)
    finite = np.isfinite(ymin) & np.isfinite(ymax)
    x, ymin, ymax = x[finite], ymin[finite], ymax[finite]
    method = series.get("method")
    descr = {
        "summary": "min/max envelope",
        "reduce": "min/max envelope",
        "sorted": "min/max envelope",
    }.get(method, "sampled — may miss extremes")
    return x, ymin, ymax, descr


class PlotRangeScreen(ModalScreen["tuple[int, int] | None"]):
    """Small modal asking for an explicit ``start:stop`` row range."""

    CSS = """
    PlotRangeScreen {
        align: center middle;
    }
    #range-dialog {
        width: 50;
        height: auto;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #range-title {
        text-style: bold;
        margin-bottom: 1;
    }
    """

    BINDINGS: ClassVar = [("escape", "cancel", "Cancel")]

    def __init__(self, *, n: int, start: int, stop: int):
        super().__init__()
        self.n = n
        self.start = start
        self.stop = stop

    def compose(self) -> ComposeResult:
        with Vertical(id="range-dialog"):
            yield Static(
                f"Row range start:stop within 0..{self.n} (current {self.start}:{self.stop})",
                id="range-title",
            )
            yield Input(placeholder="start:stop", id="range-input")

    def on_mount(self) -> None:
        widget = self.query_one("#range-input", Input)
        widget.value = f"{self.start}:{self.stop}"
        widget.focus()

    def _parse(self, text: str) -> tuple[int, int] | None:
        if ":" not in text:
            return None
        lo, hi = text.split(":", 1)
        try:
            start = int(lo) if lo.strip() else 0
            stop = int(hi) if hi.strip() else self.n
        except ValueError:
            return None
        start = max(0, min(start, self.n))
        stop = max(0, min(stop, self.n))
        return None if stop <= start else (start, stop)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        parsed = self._parse(event.value.strip().replace("_", ""))
        if parsed is None:
            self.query_one("#range-title", Static).update("Enter a range as start:stop")
            return
        self.dismiss(parsed)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PlotScreen(ModalScreen["tuple[int, int] | None"]):
    """Modal plotting one numeric column; zoomable into a row sub-range.

    Keys: ``+``/``-`` zoom about the view's left edge, ``←``/``→`` pan, ``0`` reset to
    the whole series, ``g`` type an exact ``start:stop`` range.  Each change
    re-fetches the envelope for the new range (exact for sub-ranges) via the
    *fetch* closure, so zooming reveals detail the whole-series buckets hide.

    ``v`` dismisses with the current ``(row_start, row_stop)`` so the caller can
    jump the data grid to the range you navigated to; closing dismisses ``None``.
    """

    CSS = """
    PlotScreen {
        align: center middle;
    }
    #plot-dialog {
        width: 90%;
        height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #plot-title {
        text-style: bold;
        height: 1;
    }
    #plot-widget {
        height: 1fr;
    }
    #plot-keys {
        height: 1;
        color: $text-muted;
    }
    """

    _KEYS_HINT = "+/- zoom · ←/→ pan · 0 reset · g range · v view rows · h hi-res · s scatter · esc close"
    _MIN_WIDTH = 16  # smallest zoom window (rows), so the envelope still reads
    _HIRES_MAX_POINTS = 50_000  # above this, the hi-res raw view is strided-sampled
    _SCATTER_MAX_POINTS = 50_000  # above this, the col-vs-col scatter is strided-sampled

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
        ("plus", "zoom_in", "Zoom in"),
        ("equals_sign", "zoom_in", "Zoom in"),
        ("minus", "zoom_out", "Zoom out"),
        ("left", "pan_left", "Pan left"),
        ("right", "pan_right", "Pan right"),
        ("0", "reset_range", "Reset"),
        ("g", "goto_range", "Range"),
        ("v", "view_range", "View rows"),
        ("h", "hires", "High-res"),
        ("s", "scatter", "Scatter"),
    ]

    def __init__(
        self,
        *,
        title_prefix: str,
        fetch,
        n: int,
        row_start: int,
        row_stop: int,
        series: dict,
        raw_fetch=None,
        xcol: str | None = None,
        scatter_fetch=None,
        scatter_columns: list[str] | None = None,
    ):
        super().__init__()
        self.title_prefix = title_prefix
        self._fetch = fetch
        self._raw_fetch = raw_fetch  # (start, stop) -> {"x", "y", ...} raw read
        self._xcol = xcol  # X column name, for scatter labels
        self._scatter_fetch = scatter_fetch  # (ycol, start, stop) -> series, or None
        self._scatter_columns = scatter_columns or []
        self.n = n
        self.row_start = row_start
        self.row_stop = row_stop
        self._apply(series)

    def _apply(self, series: dict) -> None:
        x, ymin, ymax, descr = _plot_view(series)
        self.x = list(x)
        self.ymin = list(ymin)
        self.ymax = list(ymax)
        full = self.row_start == 0 and self.row_stop == self.n
        rng = "" if full else f" · rows {self.row_start}:{self.row_stop}"
        note = "" if self.x else " · (no finite values in range)"
        self.plot_title = f"{self.title_prefix} · {self.n} rows{rng} · {descr}{note}"

    def compose(self) -> ComposeResult:
        with Vertical(id="plot-dialog"):
            yield Static(markup_escape(self.plot_title), id="plot-title")
            yield PlotextPlot(id="plot-widget")
            yield Static(self._KEYS_HINT, id="plot-keys")

    def on_mount(self) -> None:
        self._redraw()

    def _redraw(self) -> None:
        widget = self.query_one(PlotextPlot)
        plt = widget.plt
        plt.clear_figure()
        if self.x:
            # Upper (max) and lower (min) envelope; a single line when they
            # coincide (a sampled series).
            plt.plot(self.x, self.ymax, marker="braille")
            if self.ymin != self.ymax:
                plt.plot(self.x, self.ymin, marker="braille")
        plt.xlabel("row")
        widget.refresh()
        self.query_one("#plot-title", Static).update(markup_escape(self.plot_title))

    def _set_range(self, start: int, stop: int) -> None:
        start = max(0, min(int(start), self.n))
        stop = max(0, min(int(stop), self.n))
        if stop <= start or (start, stop) == (self.row_start, self.row_stop):
            return
        self.row_start, self.row_stop = start, stop
        self._apply(self._fetch(start, stop))
        self._redraw()

    def _zoom(self, factor: float) -> None:
        width = self.row_stop - self.row_start
        new_w = width // 2 if factor < 1 else width * 2
        new_w = max(min(self._MIN_WIDTH, self.n), min(self.n, new_w))
        # Anchor on the left edge so the zoomed plot starts where it did before.
        start = max(0, min(self.row_start, self.n - new_w))
        self._set_range(start, start + new_w)

    def _pan(self, direction: int) -> None:
        width = self.row_stop - self.row_start
        delta = max(1, width // 4) * direction
        start = max(0, min(self.row_start + delta, self.n - width))
        self._set_range(start, start + width)

    def action_zoom_in(self) -> None:
        self._zoom(0.5)

    def action_zoom_out(self) -> None:
        self._zoom(2.0)

    def action_pan_left(self) -> None:
        self._pan(-1)

    def action_pan_right(self) -> None:
        self._pan(1)

    def action_reset_range(self) -> None:
        self._set_range(0, self.n)

    def action_goto_range(self) -> None:
        def _on_range(result: tuple[int, int] | None) -> None:
            if result is not None:
                self._set_range(*result)

        self.app.push_screen(PlotRangeScreen(n=self.n, start=self.row_start, stop=self.row_stop), _on_range)

    def action_view_range(self) -> None:
        """v key — close the plot and jump the data grid to the current range."""
        self.dismiss((self.row_start, self.row_stop))

    def action_hires(self) -> None:
        """h key — open a high-res matplotlib min/max envelope of the current range.

        Pushed on top of this screen so ``q`` returns to the braille view with
        the zoom intact.  It reuses the on-screen envelope data (so it always
        renders, no zoom gate); inside it ``r`` toggles to raw values, sampled
        when the range is wide (see ``read_series``' ``max_points``).
        """
        if self._raw_fetch is None or TextualImage is None or not _matplotlib_available():
            self.app.notify(
                "High-res view needs the 'textual-image' and 'matplotlib' packages",
                severity="warning",
            )
            return
        if not self.x:
            self.app.notify("No finite values in range", severity="warning")
            return
        self.app.push_screen(
            HiResPlotScreen(
                mode="envelope",
                title_prefix=self.title_prefix,
                n=self.n,
                row_start=self.row_start,
                row_stop=self.row_stop,
                envelope={"x": self.x, "ymin": self.ymin, "ymax": self.ymax},
                raw_fetch=self._raw_fetch,
            )
        )

    def action_scatter(self) -> None:
        """s key — scatter the current X column against a chosen Y column.

        The current zoom/position fixes the row range first, so the scatter read
        is bounded; only CTable sources (which carry a *scatter_fetch* closure)
        support it.  The Y column is picked from the visible-column universe.
        """
        if self._scatter_fetch is None or not self._scatter_columns:
            self.app.notify("Scatter needs a CTable column", severity="warning")
            return

        def _on_ycol(index: int | None) -> None:
            if index is None:
                return
            ycol = self._scatter_columns[index]
            try:
                series = self._scatter_fetch(ycol, self.row_start, self.row_stop)
            except ValueError as exc:
                self.app.notify(str(exc), severity="warning")
                return
            x = np.asarray(series["x"], dtype=np.float64)
            y = np.asarray(series["y"], dtype=np.float64)
            finite = np.isfinite(x) & np.isfinite(y)
            if not finite.any():
                self.app.notify("No finite points in range", severity="warning")
                return
            if series["sampled"]:
                self.app.notify(
                    f"Sampled {series['shown']} of "
                    f"{series['row_stop'] - series['row_start']} rows (stride {series['stride']})",
                    severity="warning",
                )
            self.app.push_screen(
                ScatterPlotScreen(
                    xcol=self._xcol or "x",
                    ycol=ycol,
                    series=series,
                )
            )

        self.app.push_screen(
            ColumnSelectScreen(names=self._scatter_columns, title="Scatter Y column"),
            _on_ycol,
        )

    def action_close(self) -> None:
        self.dismiss(None)


class ScatterPlotScreen(ModalScreen[None]):
    """A col-vs-col scatter (X column vs a chosen Y column) over a row range.

    Pushed on top of :class:`PlotScreen` by the ``s`` key once the user has
    framed a row range; both columns are read row-aligned over that range (see
    :meth:`B2View.read_xy`).  ``q``/``esc`` return to the braille plot
    underneath.  No zoom in v1.
    """

    CSS = """
    ScatterPlotScreen {
        align: center middle;
    }
    #scatter-dialog {
        width: 90%;
        height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #scatter-title {
        text-style: bold;
        height: 1;
    }
    #scatter-widget {
        height: 1fr;
    }
    #scatter-keys {
        height: 1;
        color: $text-muted;
    }
    """

    _KEYS_HINT = "h hi-res · esc back to plot"

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
        ("h", "hires", "High-res"),
    ]

    def __init__(self, *, xcol: str, ycol: str, series: dict):
        super().__init__()
        self.xcol = xcol
        self.ycol = ycol
        x = np.asarray(series["x"], dtype=np.float64)
        y = np.asarray(series["y"], dtype=np.float64)
        finite = np.isfinite(x) & np.isfinite(y)
        self.x = list(x[finite])
        self.y = list(y[finite])
        title = f"{xcol} vs {ycol} · rows {series['row_start']}:{series['row_stop']}"
        if series["sampled"]:
            width = series["row_stop"] - series["row_start"]
            title += f" · sampled {series['shown']}/{width} (stride {series['stride']})"
        self.plot_title = title

    def compose(self) -> ComposeResult:
        with Vertical(id="scatter-dialog"):
            yield Static(markup_escape(self.plot_title), id="scatter-title")
            yield PlotextPlot(id="scatter-widget")
            yield Static(self._KEYS_HINT, id="scatter-keys")

    def on_mount(self) -> None:
        widget = self.query_one(PlotextPlot)
        plt = widget.plt
        plt.clear_figure()
        if self.x:
            # braille gives 2x4 sub-cell resolution, matching PlotScreen's lines
            # (the default scatter marker is one full cell per point).
            plt.scatter(self.x, self.y, marker="braille")
        plt.xlabel(self.xcol)
        plt.ylabel(self.ycol)
        widget.refresh()

    def action_hires(self) -> None:
        """h key — open a high-res matplotlib scatter over the braille scatter.

        The point set is already bounded (read_xy strided it to fit), so no zoom
        gate is needed here, unlike PlotScreen's hi-res line view.
        """
        if TextualImage is None or not _matplotlib_available():
            self.app.notify(
                "High-res view needs the 'textual-image' and 'matplotlib' packages",
                severity="warning",
            )
            return
        if not self.x:
            self.app.notify("No finite points to plot", severity="warning")
            return
        self.app.push_screen(
            HiResPlotScreen(
                mode="scatter",
                title=self.plot_title,
                scatter_xy=(self.x, self.y),
                xlabel=self.xcol,
                ylabel=self.ycol,
            )
        )

    def action_close(self) -> None:
        # pop_screen (not dismiss): pushed without a result callback, so this
        # returns to the braille PlotScreen with its zoom intact.
        self.app.pop_screen()


def _matplotlib_available() -> bool:
    """Whether matplotlib can be imported (the high-res view needs it)."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


class HiResPlotScreen(ModalScreen[None]):
    """A high-res matplotlib image over the braille plot, in one of three modes.

    Rendered with matplotlib (Agg) to a PNG and shown via ``textual-image``,
    which auto-selects the best terminal protocol (kitty/iTerm2/sixel) and
    degrades to colored half-cells elsewhere.  ``q``/``esc``/``h`` return to the
    braille view underneath with its zoom intact.

    Modes:

    - ``"scatter"`` — a static col-vs-col scatter (from ``ScatterPlotScreen``).
    - ``"bar"`` / ``"stem"`` — a grouped result: bars for a categorical key,
      impulses (vertical lines to a 0 baseline) for a numeric key.
    - ``"envelope"`` — the min/max envelope of a column (from ``PlotScreen``'s
      ``h``), reusing the on-screen braille envelope data.
    - ``"raw"`` — the column's raw values, reached by toggling with ``r`` from
      envelope mode; fetched lazily (and strided-sampled when wide).

    Envelope and raw share one screen: ``r`` toggles between them in place.
    """

    CSS = """
    HiResPlotScreen {
        align: center middle;
    }
    #hires-dialog {
        width: 95%;
        height: 90%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #hires-title {
        text-style: bold;
        height: 1;
    }
    #hires-body {
        height: 1fr;
        width: 1fr;
        align: center middle;
    }
    #hires-image {
        width: 100%;
        height: 100%;
    }
    #hires-keys {
        height: 1;
        color: $text-muted;
    }
    """

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
        ("r", "toggle_raw", "Raw/envelope"),
    ]

    def __init__(
        self,
        *,
        mode: str = "raw",
        xlabel: str = "row",
        ylabel: str | None = None,
        # scatter / bar / stem modes:
        title: str | None = None,
        scatter_xy: tuple | None = None,
        bar_data: tuple | None = None,
        stem_data: tuple | None = None,
        # column (envelope/raw) modes:
        title_prefix: str | None = None,
        n: int | None = None,
        row_start: int | None = None,
        row_stop: int | None = None,
        envelope: dict | None = None,
        raw_fetch=None,
    ):
        super().__init__()
        self._mode = mode
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._static_title = title
        self._scatter_xy = scatter_xy
        self._bar_data = bar_data
        self._stem_data = stem_data
        self._title_prefix = title_prefix
        self._n = n
        self._row_start = row_start
        self._row_stop = row_stop
        self._envelope = envelope
        self._raw_fetch = raw_fetch
        self._raw_series: dict | None = None  # lazily fetched on first 'r'
        # Toggling is offered only for the column path (envelope + a raw fetch).
        self._can_toggle = envelope is not None and raw_fetch is not None

    @property
    def _keys_hint(self) -> str:
        if self._can_toggle:
            return "r raw/envelope · esc back to braille"
        if self._mode in ("bar", "stem"):
            return "esc · back to chart"
        return "esc · back to braille"

    def _current_title(self) -> str:
        if self._mode in ("scatter", "bar", "stem"):
            return self._static_title or ""
        full = self._row_start == 0 and self._row_stop == self._n
        rng = "" if full else f" · rows {self._row_start}:{self._row_stop}"
        base = f"{self._title_prefix} · {self._n} rows{rng}"
        if self._mode == "envelope":
            return f"{base} · min/max envelope"
        s = self._raw_series
        if s is not None and s.get("sampled"):
            width = s["row_stop"] - s["row_start"]
            return f"{base} · raw values · sampled {s['shown']}/{width} (stride {s['stride']})"
        return f"{base} · raw values"

    def compose(self) -> ComposeResult:
        with Vertical(id="hires-dialog"):
            yield Static(markup_escape(self._current_title()), id="hires-title")
            # A VerticalScroll is focusable, so the screen's key bindings fire
            # (the image widget itself is not focusable).
            yield VerticalScroll(id="hires-body")
            yield Static(self._keys_hint, id="hires-keys")

    async def on_mount(self) -> None:
        self.query_one("#hires-body", VerticalScroll).focus()
        await self._show()

    async def _show(self) -> None:
        """(Re)render the current mode into the image body and refresh the title."""
        body = self.query_one("#hires-body", VerticalScroll)
        # Await the removal so a re-render (the 'r' toggle) doesn't collide with
        # the still-present image widget on its #hires-image id.
        await body.remove_children()
        self.query_one("#hires-title", Static).update(markup_escape(self._current_title()))
        try:
            png = self._render_png()
        except Exception as exc:  # pragma: no cover - defensive
            await body.mount(Static(f"Could not render: {exc}"))
            return
        await body.mount(TextualImage(io.BytesIO(png), id="hires-image"))

    def _render_png(self) -> bytes:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 6), dpi=110)
        if self._mode == "scatter":
            x, y = self._scatter_xy
            ax.scatter(x, y, s=6, color="#1f77b4", alpha=0.6)
        elif self._mode == "bar":
            labels, values = self._bar_data
            positions = range(len(labels))
            ax.bar(positions, values, color="#1f77b4")
            ax.set_xticks(list(positions))
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        elif self._mode == "stem":
            # Impulses from a 0 baseline — discrete groups, not a continuous line.
            x, y = self._stem_data
            ax.vlines(x, 0, y, linewidth=0.8, color="#1f77b4")
            ax.set_ylim(bottom=0)
            ax.margins(x=0)
        elif self._mode == "envelope":
            env = self._envelope
            x, ymin, ymax = env["x"], env["ymin"], env["ymax"]
            ax.fill_between(x, ymin, ymax, color="#1f77b4", alpha=0.3)
            ax.plot(x, ymax, linewidth=0.6, color="#1f77b4")
            if list(ymin) != list(ymax):  # a single line when the band collapses
                ax.plot(x, ymin, linewidth=0.6, color="#1f77b4")
            ax.margins(x=0)
        else:  # raw
            s = self._raw_series
            ax.plot(s["x"], s["y"], linewidth=0.8, color="#1f77b4")
            ax.margins(x=0)
        ax.set_title(self._current_title(), fontsize=10)
        ax.set_xlabel(self._xlabel)
        if self._ylabel is not None:
            ax.set_ylabel(self._ylabel)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        plt.close(fig)
        return buf.getvalue()

    async def action_toggle_raw(self) -> None:
        """r key — toggle the column view between min/max envelope and raw values.

        No-op in scatter mode.  The raw series is fetched lazily on the first
        switch (and strided-sampled when wide, see read_series' max_points).
        """
        if not self._can_toggle:
            return
        if self._mode == "envelope":
            if self._raw_series is None:
                try:
                    self._raw_series = self._raw_fetch(self._row_start, self._row_stop)
                except Exception as exc:  # pragma: no cover - defensive
                    self.app.notify(f"Could not read raw values: {exc}", severity="warning")
                    return
            self._mode = "raw"
        else:
            self._mode = "envelope"
        await self._show()

    def action_close(self) -> None:
        # pop_screen (not dismiss): this screen is pushed without a result
        # callback, so it returns to the braille PlotScreen with its zoom intact.
        self.app.pop_screen()


class CellDetailScreen(ModalScreen[None]):
    """Pretty-printed view of a single decoded CTable cell.

    Reached with Return on an expensive (list/struct/object/ndarray) column
    whose grid cell shows a ``<...; skipped>`` placeholder; the value is decoded
    on demand.  The table stays underneath with its position intact (esc
    returns).
    """

    CSS = """
    CellDetailScreen {
        align: center middle;
    }
    #cell-dialog {
        width: 80%;
        height: auto;
        max-height: 90%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    #cell-title {
        text-style: bold;
        height: 1;
    }
    #cell-body {
        height: auto;
        max-height: 1fr;
        margin-top: 1;
    }
    #cell-keys {
        height: 1;
        color: $text-muted;
        margin-top: 1;
    }
    """

    BINDINGS: ClassVar = [
        ("escape", "close", "Close"),
        ("q", "app.quit", "Quit b2view"),
    ]

    def __init__(self, *, row: int, name: str, label: str, value: Any):
        super().__init__()
        self._row = row
        self._name = name
        self._label = label
        self._value = value

    def compose(self) -> ComposeResult:
        import pprint

        title = f"row {self._row} · {self._name} ({self._label})"
        text = pprint.pformat(self._value, width=100, sort_dicts=False)
        with Vertical(id="cell-dialog"):
            yield Static(markup_escape(title), id="cell-title")
            # A VerticalScroll is focusable, so the screen's key bindings fire.
            with VerticalScroll(id="cell-body"):
                yield Static(markup_escape(text))
            yield Static("esc · close", id="cell-keys")

    def on_mount(self) -> None:
        self.query_one("#cell-body", VerticalScroll).focus()

    def action_close(self) -> None:
        self.app.pop_screen()


def _http_download(url: str, dest: str, on_progress) -> None:
    """Stream *url* to *dest*, reporting ``on_progress(downloaded, total)``.

    Module-level (and free of Textual) so it can be monkeypatched in tests.
    Writes to a temp file beside *dest* and atomically renames on success, so an
    interrupted download never leaves a corrupt file at the final name.  *total*
    is ``None`` when the server sends no ``Content-Length``.
    """
    import httpx

    tmp = dest + ".part"
    with httpx.stream("GET", url, timeout=30) as resp:
        resp.raise_for_status()
        length = resp.headers.get("Content-Length")
        total = int(length) if length is not None else None
        downloaded = 0
        on_progress(0, total)
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_bytes(chunk_size=1 << 16):
                if not chunk:
                    continue
                fh.write(chunk)
                downloaded += len(chunk)
                on_progress(downloaded, total)
    os.replace(tmp, dest)


def _fetch_remote_size(info_url: str) -> int | None:
    """Return the stored size (``cbytes``) from the bundle's info endpoint, or None.

    The download stream itself sends no ``Content-Length``, so this companion
    metadata call is what makes the progress bar determinate.  Any failure
    (network, missing key) just yields None -> an indeterminate bar.
    """
    import httpx

    try:
        resp = httpx.get(info_url, timeout=15)
        resp.raise_for_status()
        return int(resp.json()["cbytes"])
    except Exception:
        return None


class DownloadScreen(ModalScreen["bool | str"]):
    """Centered "Downloading … Please wait…" message with a progress bar.

    Pushed at startup when ``--download`` needs to fetch the bundle.  Runs the
    download on a worker thread and dismisses with ``True`` on success or the
    error string on failure; the app opens the file (or exits) from there.
    """

    CSS = """
    DownloadScreen {
        align: center middle;
    }
    #download-dialog {
        width: auto;
        height: auto;
        border: thick $accent;
        background: $surface;
        padding: 2 4;
    }
    #download-message {
        text-style: bold;
        margin-bottom: 1;
        width: 100%;
        content-align: center middle;
    }
    """

    def __init__(
        self,
        url: str,
        *,
        dest: str,
        name: str,
        info_url: str | None = None,
        source_url: str | None = None,
    ):
        super().__init__()
        self._url = url
        self._dest = dest
        self._name = name
        self._info_url = info_url
        self._source_url = source_url

    def compose(self) -> ComposeResult:
        message = f"Downloading {self._name} file"
        if self._source_url:
            message += f"\nfrom {self._source_url}"
        with Vertical(id="download-dialog"):
            yield Static(message, id="download-message")
            yield ProgressBar(id="download-bar", show_eta=True)

    def on_mount(self) -> None:
        self._run_download()

    @work(thread=True)
    def _run_download(self) -> None:
        bar = self.query_one("#download-bar", ProgressBar)
        # Prefer the info endpoint's size (the download stream omits
        # Content-Length); fall back to it if info is unavailable.
        size = _fetch_remote_size(self._info_url) if self._info_url else None

        def on_progress(downloaded: int, content_total: int | None) -> None:
            total = size or content_total
            # An unknown total leaves the bar indeterminate (pulsing).  Marshal
            # the update onto the UI thread.
            if total:
                self.app.call_from_thread(bar.update, total=total, progress=downloaded)
            else:
                self.app.call_from_thread(bar.update, progress=downloaded)

        try:
            _http_download(self._url, self._dest, on_progress)
        except Exception as exc:  # network/IO failure -> surface to the app
            self.app.call_from_thread(self.dismiss, str(exc))
            return
        self.app.call_from_thread(self.dismiss, True)


class B2ViewHeader(Header):
    """App header that also shows the open bundle's filename, left of the title.

    The filename is rendered *into the stock ``HeaderTitle`` widget* (in the
    space left of the centered "b2view — Python-Blosc2 X" title) rather than as
    extra docked child widgets.  Adding docked children to the Header was found
    to break Tab focus cycling between the panels under the Windows test driver,
    so this keeps the Header's widget tree exactly as Textual builds it and only
    overrides what the title renders.  The filename takes only the room left over
    once the centered title is reserved, truncating (with an ellipsis) as the
    terminal narrows.  Set the name with :meth:`set_filename`.
    """

    _GAP = 2  # cells kept between the filename and the centered title

    DEFAULT_CSS = """
    B2ViewHeader {
        background: $primary;  /* Blosc2 turquoise brand bar */
        color: $foreground;
    }
    B2ViewHeader HeaderIcon {
        color: $accent;  /* yellow command-palette glyph */
        width: 4;  /* tighten the gap to the filename (stock is 8) */
    }
    B2ViewHeader HeaderTitle {
        content-align: left middle;  /* we place the title ourselves */
    }
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._label = ""

    def set_filename(self, label: str) -> None:
        self._label = label
        self._refresh_title()

    def on_resize(self) -> None:
        self._refresh_title()

    def _refresh_title(self) -> None:
        with contextlib.suppress(NoMatches):  # HeaderTitle may not be composed yet
            self.query_one(HeaderTitle).update(self.format_title())

    def format_title(self) -> Content:
        """Render the centered title with the filename in the left gutter.

        With no filename (or no room for one) this is just the stock centered
        title.  Otherwise the title is left-padded so it stays centered across
        the ``HeaderTitle`` region, and the filename fills the left gutter.
        """
        base = super().format_title()
        if not self._label:
            return base
        try:
            width = self.query_one(HeaderTitle).content_size.width
        except NoMatches:
            return base
        title_len = base.cell_length
        if width <= 0 or title_len >= width:
            return base  # not even room for the title alone
        # Left padding that centers the title across the full HeaderTitle width.
        left_pad = (width - title_len) // 2
        avail = left_pad - self._GAP  # cells the filename may use
        if avail < 2:  # below this there is not even room for " x"
            return base
        label = f" {self._label}"
        if len(label) > avail:
            label = label[: avail - 1] + "…"
        gutter = " " * (left_pad - len(label))
        return Content.assemble((label, "italic dim"), gutter, base)


class B2ViewApp(App):
    """Browse TreeStore hierarchy and preview objects."""

    TITLE = "b2view"  # header title (defaults to the class name otherwise)

    # Keep escape on our own layered exit (action_dim_exit: dim mode -> locked
    # row window -> row filter -> column filter).  Textual's default would have
    # escape *restore* a maximized panel first, silently shadowing that ladder
    # (e.g. escape after a plot's `v` would un-maximize instead of unlocking the
    # window).  Restoring a maximized panel stays on its dedicated `r` key.
    ESCAPE_TO_MINIMIZE = False

    CSS = """
    #main { height: 1fr; }
    #tree-pane { width: 35%; border: solid $primary; }
    #right-pane { width: 65%; }
    #top-row { height: 40%; }
    #meta-pane, #vlmeta-pane { width: 50%; border: solid $secondary; }
    #data-pane { height: 60%; border: solid $secondary; }
    #tree { height: 1fr; }
    #data-header { height: auto; padding: 0 1; }
    #data-table-row { height: 1fr; }
    #data-table { width: 1fr; height: 1fr; }
    #row-scrollbar { width: 1; height: 1fr; color: $primary; }
    #col-scrollbar { height: 1; width: 1fr; color: $primary; }
    #meta-scroll, #vlmeta-scroll, #data-scroll { height: 1fr; padding: 0 1; }
    #tree-pane:focus-within, #meta-pane:focus-within, #vlmeta-pane:focus-within, #data-pane:focus-within { border: heavy $accent; }
    B2ViewPanel.-maximized,
    #tree-pane.-maximized,
    #meta-pane.-maximized,
    #data-pane.-maximized { width: 1fr; height: 1fr; }
    """

    BINDINGS: ClassVar = [
        ("q", "quit", "Quit"),
        ("question_mark", "show_help", "Help"),
        ("tab", "focus_next_panel", "Next panel"),
        ("shift+tab", "focus_previous_panel", "Previous panel"),
        Binding("g", "go_to_row", "Go to row", show=False),
        ("m", "maximize_panel", "Maximize"),
        ("r", "restore_or_refresh", "Restore/Refresh"),
        Binding("t", "grid_row_top", "Top", show=False),
        Binding("b", "grid_row_bottom", "Bottom", show=False),
        Binding("s", "grid_col_start", "Row start", show=False),
        Binding("e", "grid_col_end", "Row end", show=False),
        Binding("c", "go_to_column", "Go to column", show=False),
        Binding("f", "filter_rows", "Filter rows", show=False),
        Binding("S", "sort_rows", "Sort by", show=False),
        Binding("R", "reverse_sort", "Reverse sort", show=False),
        Binding("G", "group_rows", "Group by", show=False),
        Binding("slash", "filter_columns", "Filter columns", show=False),
        Binding("p", "plot_column", "Plot column", show=False),
        Binding("d", "dim_cycle", "Dim mode", show=False),
        Binding("enter", "dim_toggle_nav", "Toggle nav", show=False),
        Binding("escape", "dim_exit", "Exit dim mode", show=False),
    ]

    def __init__(
        self,
        urlpath: str,
        *,
        start_path: str = "/",
        start_panel: str = "tree",
        start_maximized: bool = False,
        preview_rows: int = 20,
        preview_cols: int = 10,
        download_url: str | None = None,
        info_url: str | None = None,
    ):
        super().__init__()
        self.sub_title = f"Python-Blosc2 {blosc2.__version__}"  # shown beside the title in the header
        self.urlpath = urlpath
        self.download_url = download_url  # when set, fetch urlpath before browsing
        self.info_url = info_url  # optional: metadata endpoint giving the size
        # Header label: the path as given on the CLI, or the @public-relative
        # path for a download (set in on_mount once that is known).
        self._header_label = urlpath
        self.start_path = start_path
        self.start_panel = start_panel
        self.start_maximized = start_maximized
        self.preview_rows = preview_rows
        self.preview_cols = preview_cols
        self.browser: StoreBrowser | None = None
        self.loaded_paths: set[str] = set()
        self.selected_path = "/"
        self.table_page: dict | None = None
        self.table_buffer: dict | None = None
        self.grid_col_start = 0
        # Sticky visible-column count: (layout key, count).  Keeps the column
        # set stable across vertical scroll / sort reverse (see _load_table_page).
        self._col_fit: tuple[tuple, int] | None = None
        self._data_layout: DataSliceLayout | None = None
        self._active_dim = 0
        self._dim_mode = False
        self.loading_table_page = False
        # One-shot: apply the --panel start focus after the first update_panels,
        # once the data panel's display/contents have settled (see update_panels).
        self._apply_focus_on_next_update = False
        # Absolute (start, stop) of a locked row window from the plot's 'v' key.
        self.row_window: tuple[int, int] | None = None
        # Last applied group-by config (key, op, value_col), reused to pre-fill
        # the 'G' modal on any table — survives navigating to other nodes.
        self._last_group: tuple[str, str, str | None] | None = None

    def compose(self) -> ComposeResult:
        yield B2ViewHeader()
        with Horizontal(id="main"):
            with B2ViewPanel(id="tree-pane") as tree_pane:
                tree_pane.border_title = "tree"
                yield Tree("/", id="tree")
            with Vertical(id="right-pane"):
                with Horizontal(id="top-row"):
                    with B2ViewPanel(id="meta-pane") as meta_pane:
                        meta_pane.border_title = "meta"
                        with VerticalScroll(id="meta-scroll", can_focus=True):
                            yield Static("Select a node", id="metadata")
                    with B2ViewPanel(id="vlmeta-pane") as vlmeta_pane:
                        vlmeta_pane.border_title = "vlmeta"
                        with VerticalScroll(id="vlmeta-scroll", can_focus=True):
                            yield Static("", id="vlmetadata")
                with B2ViewPanel(id="data-pane") as data_pane:
                    data_pane.border_title = "data"
                    data_pane.border_subtitle = (
                        "?(help) | d(im mode) | filter: f(rows) /(cols) | S(ort) | G(roup) | "
                        "rows: t/b/g(oto) | cols: s/e/c(goto) | p(lot)"
                    )
                    yield Static("", id="data-header")
                    with Horizontal(id="data-table-row"):
                        yield BufferedDataTable(id="data-table", show_row_labels=True, zebra_stripes=True)
                        yield Static("", id="row-scrollbar")
                    yield Static("", id="col-scrollbar")
                    with VerticalScroll(id="data-scroll", can_focus=True):
                        yield Static("", id="preview")
        yield Footer()

    def on_mount(self) -> None:
        self.register_theme(BLOSC2_THEME)
        self.theme = "blosc2"
        if self.download_url:
            # Fetch the bundle first, then open it from _after_download.  The
            # message shows the @public-relative path (e.g. "large/foo.b2z"),
            # falling back to the local basename for any other URL shape.
            name = os.path.basename(self.urlpath)
            if "/@public/" in self.download_url:
                name = self.download_url.split("/@public/", 1)[1]
            self._header_label = name  # @public-relative path for the header
            # A browsable URL for the source root, so the user can see where the
            # file comes from: e.g. https://cat2.cloud/demo/?roots=@public
            source_url = None
            if "/api/" in self.download_url:
                source_url = self.download_url.split("/api/", 1)[0] + "/?roots=@public"
            self.push_screen(
                DownloadScreen(
                    self.download_url,
                    dest=self.urlpath,
                    name=name,
                    info_url=self.info_url,
                    source_url=source_url,
                ),
                self._after_download,
            )
        else:
            self._start_browsing()

    def _after_download(self, result: bool | str) -> None:
        """Resume browsing once DownloadScreen finishes (True ok, else error)."""
        if result is True:
            self._start_browsing()
        else:
            self.exit(message=f"Download failed: {result}")

    def _start_browsing(self) -> None:
        """Open the bundle and populate the tree (the normal startup path)."""
        self.browser = StoreBrowser(self.urlpath)
        self.query_one(B2ViewHeader).set_filename(self._header_label)
        tree = self.query_one("#tree", Tree)
        tree.root.data = "/"
        self.load_children(tree.root)
        tree.root.expand()
        self.query_one("#data-table-row", Horizontal).display = False
        self.query_one("#col-scrollbar", Static).display = False

        # Focus the requested start panel after the first update_panels has set
        # up the data panel (its display and contents), so 'data' lands on the
        # populated grid instead of racing the node selection.
        self._apply_focus_on_next_update = True
        if self.start_path and self.start_path != "/":
            self._navigate_to_path(self.start_path)
        else:
            self.call_after_refresh(self.update_panels, "/")

    def _apply_start_focus(self) -> None:
        """Focus the panel requested on startup (the --panel option), and
        maximize it too when --max was given."""
        self._focus_panel_by_name(self.start_panel)
        if self.start_maximized:
            # Defer a frame: .focus() above is scheduled, so the new focus (which
            # action_maximize_panel reads) isn't applied yet this tick.
            self.call_after_refresh(self.action_maximize_panel)

    def _focus_panel_by_name(self, name: str) -> None:
        """Focus a panel by its user-facing name."""
        panel_map = {
            "tree": lambda: self.query_one("#tree", Tree),
            "meta": lambda: self.query_one("#meta-scroll", VerticalScroll),
            "vlmeta": lambda: self.query_one("#vlmeta-scroll", VerticalScroll),
            "data": lambda: (
                self.query_one("#data-table", DataTable)
                if self.query_one("#data-table-row", Horizontal).display
                else self.query_one("#data-scroll", VerticalScroll)
            ),
        }
        getter = panel_map.get(name)
        if getter is not None:
            getter().focus()

    def _navigate_to_path(self, path: str) -> None:
        """Expand the tree and select the node at *path*."""
        tree = self.query_one("#tree", Tree)
        parts = [p for p in path.split("/") if p]
        node = tree.root
        # Walk down the tree expanding each level
        for part in parts:
            self.load_children(node)
            found = None
            for child in node.children:
                if child.label and child.label.plain.endswith(f" {part}"):
                    found = child
                    break
            if found is None:
                # Path not found — fall back to root
                self.call_after_refresh(self.update_panels, "/")
                tree.focus()
                return
            if found.allow_expand:
                self.load_children(found)
            found.expand()
            node = found

        # Selecting the node fires NodeSelected → on_tree_node_selected →
        # update_panels, which applies the one-shot start-panel focus once the
        # data panel is populated (see _apply_focus_on_next_update).
        def _do_select():
            tree.select_node(node)
            tree.scroll_to_node(node)

        self.call_after_refresh(_do_select)

    def on_unmount(self) -> None:
        if self.browser is not None:
            self.browser.close()

    def load_children(self, node) -> None:
        path = node.data or "/"
        if self.browser is None or path in self.loaded_paths:
            return
        for child in self.browser.list_children(path):
            icon = _KIND_ICONS.get(child.kind, "?")
            node.add(f"{icon} {child.name}", data=child.path, allow_expand=child.has_children)
        self.loaded_paths.add(path)

    def on_tree_node_expanded(self, event: Tree.NodeExpanded) -> None:
        self.load_children(event.node)

    def on_tree_node_selected(self, event: Tree.NodeSelected) -> None:
        path = event.node.data or "/"
        self.selected_path = path
        self.update_panels(path)
        if event.node.allow_expand:
            self.load_children(event.node)

    def update_panels(self, path: str) -> None:
        if self.browser is None:
            return
        metadata = self.query_one("#metadata", Static)
        data_header = self.query_one("#data-header", Static)
        data_table_row = self.query_one("#data-table-row", Horizontal)
        data_scroll = self.query_one("#data-scroll", VerticalScroll)
        preview = self.query_one("#preview", Static)
        vlmeta_pane = self.query_one("#vlmeta-pane", B2ViewPanel)
        vlmeta_widget = self.query_one("#vlmetadata", Static)
        try:
            info = self.browser.get_info(path)
            metadata.update(make_metadata_renderable(info))
            self.table_buffer = None
            self.grid_col_start = 0
            self._data_layout = None
            self._active_dim = 0
            self._dim_mode = False
            # A locked row window does not survive navigating to a node.
            self.row_window = None
            self.browser.clear_row_window(path)
            if info.kind == "group":
                data_header.display = False
                data_table_row.display = False
                data_scroll.display = True
                self.query_one("#col-scrollbar", Static).display = False
                data_header.update("")
                preview.update("Group node; select an array or table to preview.")
                self._update_vlmeta(vlmeta_pane, vlmeta_widget, path)
            else:
                if self._uses_grid_preview(info):
                    data_header.display = True
                    data_table_row.display = True
                    data_scroll.display = False
                    preview.update("")
                    shape = tuple(info.metadata.get("shape", ()) or ())
                    ndim = len(shape)
                    if ndim >= 1 and self._data_layout is None:
                        self._data_layout = DataSliceLayout.from_shape(shape)
                        self._active_dim = 0
                    data = self._load_table_page(path, 0)
                else:
                    data = self.browser.preview(path, max_rows=self.preview_rows, max_cols=self.preview_cols)
                if self._is_table_preview(data):
                    # A freshly selected node starts at the first column
                    self._update_data_table(data, cursor_col=0)
                    self._update_data_header(data)
                    self.call_after_refresh(self._ensure_viewport_consistent)
                else:
                    header, body = make_preview_renderables(data)
                    data_header.display = header is not None
                    data_table_row.display = False
                    data_scroll.display = True
                    self.query_one("#col-scrollbar", Static).display = False
                    data_header.update("" if header is None else header)
                    preview.update(body)
            self._update_vlmeta(vlmeta_pane, vlmeta_widget, path)
            self._reset_panel_scroll()
        except Exception as exc:
            metadata.update(f"Error reading {path}: {exc}")
            data_header.display = False
            data_table_row.display = False
            data_scroll.display = True
            self.query_one("#col-scrollbar", Static).display = False
            data_header.update("")
            preview.update("")
            self._update_vlmeta(vlmeta_pane, vlmeta_widget, None)
            self._reset_panel_scroll()

        # The data panel's display/contents are now settled; apply the one-shot
        # startup focus (deferred one frame so the target widget is rendered).
        if self._apply_focus_on_next_update:
            self._apply_focus_on_next_update = False
            self.call_after_refresh(self._apply_start_focus)

    @staticmethod
    def _format_vlmeta_value(value: Any) -> str:
        """Format a vlmeta value for display."""
        if isinstance(value, bool):
            return str(value)
        if isinstance(value, (int, float)):
            return str(value)
        if isinstance(value, (list, tuple)):
            return ", ".join(str(v) for v in value)
        if isinstance(value, dict):
            return ", ".join(f"{k}: {v}" for k, v in value.items())
        return str(value)

    def _update_vlmeta(self, pane, widget: Static, path: str | None) -> None:
        """Populate the vlmeta pane with variable-length metadata."""
        pane.display = True
        if path is None or self.browser is None:
            widget.update("<not available>")
            return
        try:
            info = self.browser.get_info(path)
            if info.user_attrs is None:
                widget.update("<not available>")
            elif not info.user_attrs:
                widget.update("")
            else:
                from rich.table import Table

                table = Table(show_header=False, box=None, expand=True)
                table.add_column("key", style="bold cyan", no_wrap=True)
                table.add_column("value")
                for k, v in info.user_attrs.items():
                    table.add_row(str(k), self._format_vlmeta_value(v))
                widget.update(table)
        except Exception:
            widget.update("<not available>")

    @staticmethod
    def _is_table_preview(data) -> bool:
        return isinstance(data, dict) and "data" in data and "columns" in data

    @staticmethod
    def _uses_grid_preview(info) -> bool:
        # 1D, 2D, 3D+ NDArray/C2Array all use grid preview; SChunk uses it for
        # the paged hex dump (rows of 16 bytes).
        return info.kind in {"ctable", "schunk"} or (
            info.kind in {"ndarray", "c2array"} and info.metadata.get("ndim", 0) >= 1
        )

    def _col_page_size(self) -> int:
        """Return the number of columns that fit in the current data table width."""
        table = self.query_one("#data-table", DataTable)
        width = table.size.width
        if width <= 1:
            return self.preview_cols
        # Each column uses roughly 9 characters (float format width) + 2 padding.
        # Row labels take about 6 characters.
        col_width = 11
        # Subtract row-label column space
        usable = max(1, width - 6)
        return max(1, usable // col_width)

    # DataTable pads each cell with one space on both sides (cell_padding=1)
    _CELL_PAD = 2

    def _data_table_width(self) -> int:
        return self.query_one("#data-table", DataTable).size.width

    def _col_avail_width(self, nrows: int) -> int:
        """Width available for data columns, or 0 before layout has settled."""
        width = self._data_table_width()
        if width <= 1:
            return 0
        label_width = len(str(max(0, int(nrows) - 1))) + self._CELL_PAD
        return max(1, width - label_width)

    def _candidate_max_cols(self) -> int:
        """Upper bound of columns worth fetching before the width-based trim."""
        width = self._data_table_width()
        if width <= 1:
            return self.preview_cols
        # The narrowest possible column is one character plus padding.
        return max(1, width // (1 + self._CELL_PAD))

    @classmethod
    def _measure_column_widths(cls, data: dict) -> list[int]:
        """Rendered width (content + padding) of every column in *data*."""
        widths = []
        for name in data["columns"]:
            cells = data["data"][name]
            decimals = column_float_decimals(cells)
            content = max(
                len(str(name)),
                max((len(format_cell(value, float_decimals=decimals)) for value in cells), default=1),
            )
            widths.append(content + cls._CELL_PAD)
        return widths

    def _trim_columns_to_fit(self, data: dict) -> dict:
        """Drop trailing columns of *data* that do not fit the table width.

        The preview fetches a generous candidate window of columns; this
        second pass measures the actual rendered cell widths and keeps only
        as many whole columns as truly fit the data table.
        """
        if data.get("source_kind") not in _COL_PAGED_KINDS:
            return data
        avail = self._col_avail_width(data["nrows"])
        if avail <= 0:
            return data  # layout not settled; keep the estimate-based window
        widths = self._measure_column_widths(data)
        keep = 0
        total = 0
        for width in widths:
            if keep >= 1 and total + width > avail:
                break
            total += width
            keep += 1
        return self._take_n_columns(data, keep)

    def _take_n_columns(self, data: dict, n: int) -> dict:
        """Keep the first *n* columns of a paged *data* window (clamped to range).

        Width-based fitting (:meth:`_trim_columns_to_fit`) is recomputed from the
        currently visible rows, so it can vary as you scroll or reverse a sort.
        Pinning the count keeps the visible column set stable across those
        re-renders (the sticky fit in :meth:`_load_table_page`).
        """
        if data.get("source_kind") not in _COL_PAGED_KINDS:
            return data
        keep = max(1, min(n, len(data["columns"])))
        if keep >= len(data["columns"]):
            return data
        kept = data["columns"][:keep]
        data = dict(data)
        data["data"] = {name: data["data"][name] for name in kept}
        data["columns"] = kept
        data["col_stop"] = data["col_start"] + keep
        data["hidden_columns"] = max(0, data["ncols"] - keep)
        return data

    def _fetch_columns_for_measure(self, col_start: int, count: int) -> dict:
        """Fetch the current page rows for columns [col_start, col_start+count)."""
        page = self.table_page
        max_rows = max(1, page["stop"] - page["start"])
        layout = self._data_layout
        if layout is not None and len(layout.shape) >= 1:
            probe = layout.copy_with(row_start=page["start"], col_start=col_start)
            return self.browser.preview(self.selected_path, max_rows=max_rows, max_cols=count, layout=probe)
        return self.browser.preview(
            self.selected_path,
            start=page["start"],
            stop=page["stop"],
            max_cols=count,
            col_start=col_start,
        )

    def _fit_col_start_backward(self, end: int) -> int:
        """Start of the widest whole-column window ending just before *end*."""
        page = self.table_page
        avail = self._col_avail_width(page["nrows"])
        if avail <= 0:
            return max(0, end - max(1, self._col_page_size()))
        candidate = min(end, max(1, avail // (1 + self._CELL_PAD)))
        cand_start = end - candidate
        widths = self._measure_column_widths(self._fetch_columns_for_measure(cand_start, candidate))
        start = end - 1  # always keep at least one column
        total = widths[-1]
        for i in range(len(widths) - 2, -1, -1):
            if total + widths[i] > avail:
                break
            total += widths[i]
            start = cand_start + i
        return max(0, start)

    def _table_page_size(self) -> int:
        table = self.query_one("#data-table", DataTable)
        # Keep only rows likely to be visible.  The DataTable header consumes one
        # line; fall back to the CLI limit before layout has assigned sizes.
        height = table.size.height
        if height <= 1:
            height = self.query_one("#data-pane", Vertical).size.height - 2
        return max(1, height - 1) if height > 1 else max(1, self.preview_rows)

    def _col_fit_key(self) -> tuple:
        """Identity of the current column layout for the sticky column-count fit.

        Changes (forcing a width re-fit) on a new node, a horizontal scroll, an
        ndarray dim/fixed-value change, a column filter, or a terminal resize —
        but not on vertical scroll, sort reverse or row filter, which keep the
        same columns.
        """
        layout = self._data_layout
        layout_sig = None
        if layout is not None:
            layout_sig = (
                tuple(layout.navigable_dims),
                tuple(sorted(layout.fixed_values.items())),
                tuple(layout.shape),
            )
        col_filter = self.browser.get_column_filter(self.selected_path) if self.browser else None
        return (self.selected_path, self.grid_col_start, layout_sig, col_filter, self._data_table_width())

    def _load_table_page(self, path: str, start: int) -> dict:
        if self.browser is None:
            raise RuntimeError("Store browser is not open")
        page_size = self._table_page_size()
        start = max(0, start)
        layout = self._data_layout

        if self.table_buffer is not None:
            buffer_start = self.table_buffer["start"]
            buffer_stop = self.table_buffer["stop"]
            buffer_kind = self.table_buffer.get("source_kind")
            if buffer_kind in {"ndarray2d", "ndarray_slice"}:
                same_columns = self.table_buffer.get(
                    "col_start"
                ) == self.grid_col_start and self.table_buffer.get("slice_indices") == (
                    [
                        layout.fixed_values.get(i, 0)
                        for i in range(len(layout.shape))
                        if i in layout.fixed_values
                    ]
                    if layout is not None
                    else []
                )
            elif buffer_kind == "ctable":
                same_columns = self.table_buffer.get("col_start") == self.grid_col_start
            else:
                same_columns = True
            if same_columns and buffer_start <= start and start + page_size <= buffer_stop:
                data = self._slice_table_buffer(start, page_size)
                self.table_page = data
                return data

        buffer_size = page_size * 10
        buffer_start = max(0, start - page_size * 4)

        if layout is not None and len(layout.shape) >= 1:
            # Use the layout-based preview for all array types (1D+)
            # Scalar view (0 navigable dims) always starts at 0
            if not layout.navigable_dims:
                start = 0
            self._sync_layout_scroll(start, layout)
            data = self.browser.preview(
                path,
                max_rows=buffer_size,
                max_cols=self._candidate_max_cols(),
                layout=layout,
            )
        else:
            # CTable or non-array objects — use legacy preview
            data = self.browser.preview(
                path,
                start=buffer_start,
                stop=buffer_start + buffer_size,
                max_rows=buffer_size,
                max_cols=self._candidate_max_cols(),
                col_start=self.grid_col_start,
            )
        # The visible column count is sticky for a given column layout: recompute
        # the width-based fit only when the layout key changes (node, horizontal
        # position, ndarray dims, column filter).  Vertical scrolling, reversing a
        # sort and row filtering keep the same columns instead of dropping one
        # when the freshly visible rows happen to measure wider.
        fit_key = self._col_fit_key()
        if self._col_fit is not None and self._col_fit[0] == fit_key:
            data = self._take_n_columns(data, self._col_fit[1])
        else:
            data = self._trim_columns_to_fit(data)
            # Only remember the count once the layout has settled (a real
            # width-based fit); before that the trim is a no-op and would pin a
            # bloated count that overflows the table on later renders.
            if data.get("source_kind") in _COL_PAGED_KINDS and self._data_table_width() > 1:
                self._col_fit = (fit_key, len(data["columns"]))
        data["viewport_width"] = self._data_table_width()
        self.table_buffer = data
        data = self._slice_table_buffer(start, page_size)
        self.table_page = data
        return data

    def _sync_layout_scroll(self, start: int, layout: DataSliceLayout) -> None:
        """Update the layout's row/col scroll positions to match the page start."""
        if layout is None:
            return
        navigable = layout.navigable_dims
        if len(navigable) >= 1:
            row_dim = navigable[0]
            # A locked window shortens the navigable row dim; scroll is
            # window-relative, so clamp to the window length.
            win_lo, win_hi = layout.row_window_bounds(row_dim)
            total = win_hi - win_lo
            layout.row_start = max(0, min(start, total))
            layout.row_stop = min(layout.row_start + self._table_page_size() * 10, total)
        if len(navigable) >= 2:
            col_dim = navigable[1]
            total = layout.shape[col_dim]
            layout.col_start = max(0, min(self.grid_col_start, total))
            layout.col_stop = min(layout.col_start + self._col_page_size(), total)

    def _slice_table_buffer(self, start: int, page_size: int) -> dict:
        if self.table_buffer is None:
            raise RuntimeError("No table buffer loaded")
        buffer = self.table_buffer
        offset = start - buffer["start"]
        available = max(0, buffer["stop"] - start)
        count = min(page_size, available)
        stop = start + count
        return {
            "start": start,
            "stop": stop,
            "nrows": buffer["nrows"],
            "columns": buffer["columns"],
            "hidden_columns": buffer["hidden_columns"],
            "data": {name: values[offset : offset + count] for name, values in buffer["data"].items()},
            **(
                {"row_labels": buffer["row_labels"][offset : offset + count]}
                if "row_labels" in buffer
                else {}
            ),
            **{
                key: buffer[key]
                for key in (
                    "source_kind",
                    "shape",
                    "col_start",
                    "col_stop",
                    "ncols",
                    "slice_indices",
                    "n_slices_per_dim",
                    "viewport_width",
                    "nbytes",
                    "typesize",
                )
                if key in buffer
            },
        }

    def _update_data_table(self, data: dict, *, cursor_row: int = 0, cursor_col: int | None = None) -> None:
        """Refresh the data grid; *cursor_col* None keeps the current column."""
        table = self.query_one("#data-table", DataTable)
        if cursor_col is None:
            cursor_col = table.cursor_column
        self.loading_table_page = True
        try:
            table.clear(columns=True)
            for name in data["columns"]:
                table.add_column(name, key=name)
            # Uniform decimals per float column, taken from the whole buffer
            # when available so the format is stable while paging rows.
            buffer = self.table_buffer
            source = buffer if buffer is not None and buffer["columns"] == data["columns"] else data
            decimals = {name: column_float_decimals(source["data"][name]) for name in data["columns"]}
            nrows = data["stop"] - data["start"]
            # SChunk hex dumps carry explicit (hex byte-offset) row labels;
            # everything else labels the gutter with the logical row number.
            row_labels = data.get("row_labels")
            for i in range(nrows):
                table.add_row(
                    *[
                        format_cell(data["data"][name][i], float_decimals=decimals[name])
                        for name in data["columns"]
                    ],
                    label=row_labels[i] if row_labels is not None else str(data["start"] + i),
                )
            nrows = data["stop"] - data["start"]
            cursor_row = min(max(0, cursor_row), max(0, nrows - 1))
            cursor_col = min(max(0, cursor_col), max(0, len(data["columns"]) - 1))
            table.cursor_coordinate = (cursor_row, cursor_col)
            table.scroll_home(animate=False)
            self._update_global_row_scrollbar(data)
            self._update_global_col_scrollbar(data)
        finally:
            self.call_after_refresh(self._finish_table_page_load)

    def _finish_table_page_load(self) -> None:
        self.loading_table_page = False

    def page_table(self, direction: int, *, align: bool = False) -> bool:
        if self.loading_table_page or self.table_page is None:
            return False
        page = self.table_page
        page_size = self._table_page_size()
        if direction > 0:
            if page["stop"] >= page["nrows"]:
                return False
            # An explicit page down re-aligns to the page grid: dim-mode
            # single-row scrolls (_scroll_navigable_viewport) can leave `start`
            # off a page_size boundary, and contiguous paging from `stop` would
            # carry that offset forever.  Snapping to the next page_size
            # multiple mirrors how column paging re-fits on each page.  For an
            # already-aligned page this equals `stop`, so cursor-edge paging
            # (align=False) is unchanged.
            new_start = (page["start"] // page_size + 1) * page_size if align else page["stop"]
            data = self._load_table_page(self.selected_path, new_start)
            cursor_row = 0
        else:
            if page["start"] <= 0:
                return False
            if align:
                # Previous grid line: floor for an off-grid start, start-page
                # for an aligned one (ceil-div keeps aligned pages contiguous).
                new_start = (-(-page["start"] // page_size) - 1) * page_size
            else:
                new_start = page["start"] - page_size
            new_start = max(0, new_start)
            data = self._load_table_page(self.selected_path, new_start)
            cursor_row = data["stop"] - data["start"] - 1
        self._update_data_table(data, cursor_row=cursor_row)
        self._update_data_header(data)
        return True

    def page_grid_columns(self, direction: int) -> bool:
        if self.loading_table_page or self.table_page is None:
            return False
        page = self.table_page
        if page.get("source_kind") not in _COL_PAGED_KINDS:
            return False
        # Whole-column windows of data-dependent size: paging right starts at
        # the first hidden column; paging left fits as many whole columns as
        # possible ending just before the current first one (no skips, no gaps).
        if direction > 0:
            if page["col_stop"] >= page["ncols"]:
                return False
            self.grid_col_start = page["col_stop"]
        else:
            if page["col_start"] <= 0:
                return False
            self.grid_col_start = self._fit_col_start_backward(page["col_start"])
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        cursor_col = 0 if direction > 0 else len(data["columns"]) - 1
        self._update_data_table(data, cursor_row=cursor_row, cursor_col=cursor_col)
        self._update_data_header(data)
        return True

    def _grid_col_home(self) -> bool:
        if self.table_page is None or self.table_page.get("source_kind") not in _COL_PAGED_KINDS:
            return False
        self.grid_col_start = 0
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, self.table_page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row, cursor_col=0)
        self._update_data_header(data)
        return True

    def _grid_col_end(self) -> bool:
        if self.table_page is None or self.table_page.get("source_kind") not in _COL_PAGED_KINDS:
            return False
        page = self.table_page
        # Jump to the widest whole-column window ending at the last column
        self.grid_col_start = self._fit_col_start_backward(page["ncols"])
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row, cursor_col=len(data["columns"]) - 1)
        self._update_data_header(data)
        return True

    def _update_data_header(self, data: dict) -> None:
        layout = self._data_layout
        header_parts: list[str] = []

        if layout is not None and len(layout.shape) >= 1:
            ndim = len(layout.shape)
            for i in range(ndim):
                is_active = i == self._active_dim

                if i in layout.fixed_values:
                    idx = layout.fixed_values[i]
                    part = f"d{i} [{idx}]"
                elif i in layout.navigable_dims:
                    pos = layout.navigable_dims.index(i)
                    if pos == 0:
                        s, e = data["start"], data["stop"]
                    else:
                        s, e = data.get("col_start", 0), data.get("col_stop", 0)
                    part = f"d{i}[{s}:{e}]"
                else:
                    part = f"d{i} ?"

                if is_active and self._dim_mode:
                    part = f"[bold]{part}[/bold]"
                header_parts.append(part)

            if self._dim_mode:
                # The whole line is accent-reversed below; this chip is a cutout
                # (accent text on the dark canvas) so it stands out against it.
                header_parts.append("[$accent on $background] DIM MODE [/]")
                header_parts.append("←→dim  ↑↓val  <Enter>fix/nav  <Esc>exit")
            elif self.row_window is not None:
                ws, we = self.row_window
                header_parts.append(_accent_chip(f"WINDOW {ws}:{we}"))
                header_parts.append("<Esc>unlock")
        elif data.get("source_kind") == "schunk":
            # The hex dump is paged in 16-byte rows; report it in bytes.
            header_parts.append(f"hex dump · {data.get('nbytes', 0)} bytes")
            if data.get("typesize", 1) > 1:
                header_parts.append(f"typesize {data['typesize']}")
        else:
            header_parts.append(f"rows {data['start']}:{data['stop']} of {data['nrows']}")
            if "col_start" in data:
                header_parts.append(f"cols {data['col_start']}:{data['col_stop']} of {data['ncols']}")
            header_parts.extend(self._window_and_filter_chips(data))

        line = ", ".join(header_parts)
        if self._dim_mode and layout is not None:
            line = f"[$background on $accent]{line}[/]"
        self.query_one("#data-header", Static).update(line)

    def _window_and_filter_chips(self, data: dict) -> list[str]:
        """Header chips for a locked row window and any active CTable filters."""
        chips: list[str] = []
        if self.row_window is not None:
            ws, we = self.row_window
            chips.append(_accent_chip(f"WINDOW {ws}:{we}"))
        if data.get("source_kind") == "ctable" and self.browser is not None:
            flt = self.browser.get_filter(self.selected_path)
            col_flt = self.browser.get_column_filter(self.selected_path)
            sort = self.browser.get_sort(self.selected_path)
            group = self.browser.get_group(self.selected_path)
            gsort = self.browser.get_group_sort(self.selected_path)
            if group:
                key, op, value_col = group
                label = "count" if op == "size" else f"{op}({markup_escape(value_col)})"
                chips.append(_accent_chip(f"GROUPED {markup_escape(key)} · {label}"))
                if gsort:
                    col, reverse = gsort
                    arrow = "▼" if reverse else "▲"
                    chips.append(_accent_chip(f"SORTED {arrow} {markup_escape(col)}"))
            if flt:
                total = self.browser.base_nrows(self.selected_path)
                chips.append(f"filter: [bold]{markup_escape(flt)}[/bold] ({total} total)")
            if sort:
                col, reverse = sort
                arrow = "▼" if reverse else "▲"
                chips.append(_accent_chip(f"SORTED {arrow} {markup_escape(col)}"))
            if col_flt:
                chips.append(f"cols: [bold]{markup_escape(col_flt)}[/bold]")
            if group:
                if gsort:
                    chips.append("<R>everse")
                if group[1] in ("argmin", "argmax"):
                    chips.append("<Enter>go to row")
                chips.append("<Esc>ungroup" if not gsort else "<Esc>unsort")
            else:
                if sort:
                    chips.append("<R>everse")
                if flt or col_flt or sort or self.row_window is not None:
                    chips.append("<Esc>unlock/clear")
        return chips

    def _make_global_scrollbar(self, *, start: int, stop: int, total: int, size: int, track: str) -> str:
        size = max(1, size)
        total = max(1, total)
        start = min(max(0, start), total)
        stop = min(max(start, stop), total)
        visible = max(1, stop - start)
        thumb_size = max(1, round(size * min(1.0, visible / total)))
        if total <= visible:
            thumb_start = 0
            thumb_size = size
        else:
            thumb_start = round((size - thumb_size) * (start / (total - visible)))
        thumb_stop = min(size, thumb_start + thumb_size)
        return "".join("█" if thumb_start <= i < thumb_stop else track for i in range(size))

    def _update_global_row_scrollbar(self, data: dict) -> None:
        scrollbar = self.query_one("#row-scrollbar", Static)
        height = max(1, self.query_one("#data-table", DataTable).size.height)
        bar = self._make_global_scrollbar(
            start=int(data["start"]),
            stop=int(data["stop"]),
            total=int(data["nrows"]),
            size=height,
            track="│",
        )
        scrollbar.update("\n".join(bar))

    def _update_global_col_scrollbar(self, data: dict) -> None:
        scrollbar = self.query_one("#col-scrollbar", Static)
        if data.get("source_kind") not in _COL_PAGED_KINDS:
            scrollbar.display = False
            scrollbar.update("")
            return
        scrollbar.display = True
        width = max(1, self.query_one("#data-table", DataTable).size.width)
        scrollbar.update(
            self._make_global_scrollbar(
                start=int(data["col_start"]),
                stop=int(data["col_stop"]),
                total=int(data["ncols"]),
                size=width,
                track="─",
            )
        )

    def _reset_panel_scroll(self) -> None:
        for selector in ("#meta-scroll", "#data-scroll"):
            self.query_one(selector, VerticalScroll).scroll_home(animate=False)
        data_table_row = self.query_one("#data-table-row", Horizontal)
        if data_table_row.display:
            self.query_one("#data-table", DataTable).scroll_home(animate=False)
            if self.table_page is not None:
                self._update_global_row_scrollbar(self.table_page)
                self._update_global_col_scrollbar(self.table_page)

    def _focusable_panels(self):
        data_table_row = self.query_one("#data-table-row", Horizontal)
        data_panel = (
            self.query_one("#data-table", DataTable)
            if data_table_row.display
            else self.query_one("#data-scroll", VerticalScroll)
        )
        return [
            self.query_one("#tree", Tree),
            self.query_one("#meta-scroll", VerticalScroll),
            self.query_one("#vlmeta-scroll", VerticalScroll),
            data_panel,
        ]

    def _focus_panel(self, step: int) -> None:
        panels = self._focusable_panels()
        focused = self.focused
        try:
            index = panels.index(focused)
        except ValueError:
            index = 0 if step > 0 else len(panels) - 1
        panels[(index + step) % len(panels)].focus()

    def action_focus_next_panel(self) -> None:
        self._focus_panel(1)

    def action_focus_previous_panel(self) -> None:
        self._focus_panel(-1)

    def _in_data_grid(self) -> bool:
        """Return True if focus is inside the data pane and a grid is active."""
        if self.table_page is None:
            return False
        if not self.query_one("#data-table-row", Horizontal).display:
            return False
        focused = self.focused
        if focused is None:
            return False
        pane = self.query_one("#data-pane", Vertical)
        return focused is pane or pane in focused.ancestors

    def action_go_to_row(self) -> None:
        if not self._in_data_grid():
            return
        current = self.table_page["start"] + self.query_one("#data-table", DataTable).cursor_row
        screen = GoToRowScreen(nrows=self.table_page["nrows"], current=current)
        self.push_screen(screen, self._go_to_row)

    _PLOT_MAX_POINTS = 2000

    def _inspect_cursor_cell(self) -> bool:
        """Return on a skipped CTable cell: decode just that cell into a modal.

        Returns True when the key was consumed (the cursor sat on an expensive
        ``<...; skipped>`` cell), so the data-table's default select handler is
        skipped.  Anything else (numeric/text cells, non-CTable grids) returns
        False and falls through.
        """
        if not self._in_data_grid() or self.table_page is None or self.browser is None:
            return False
        if self.table_page.get("source_kind") != "ctable":
            return False
        # skipped_columns lives on the buffer; _slice_table_buffer drops it.
        skipped = (self.table_buffer or {}).get("skipped_columns") or {}
        if not skipped:
            return False
        columns = self.table_page["columns"]
        table = self.query_one("#data-table", DataTable)
        cursor_col = table.cursor_column
        if not (0 <= cursor_col < len(columns)):
            return False
        name = columns[cursor_col]
        if name not in skipped:
            return False
        row = self.table_page["start"] + table.cursor_row
        try:
            value = self.browser.read_cell(self.selected_path, name, row)
        except Exception as exc:  # pragma: no cover - defensive
            self.notify(f"Could not decode cell: {exc}", severity="error")
            return True  # we owned the key; surface the failure
        self.push_screen(CellDetailScreen(row=row, name=name, label=skipped[name], value=value))
        return True

    def _drilldown_arg_cell(self) -> bool:
        """On an argmin/argmax cell of a grouped view, jump to that base-table row.

        The cell holds the logical row position of the group's extreme value, so
        we ungroup and land the cursor on that row, on the column the extreme was
        computed over.  Returns True when the key was consumed.
        """
        if not self._in_data_grid() or self.table_page is None or self.browser is None:
            return False
        group = self.browser.get_group(self.selected_path)
        if group is None:
            return False
        _key, op, value_col = group
        if op not in ("argmin", "argmax"):
            return False
        agg_col = self.browser.group_agg_column(self.selected_path)
        columns = self.table_page["columns"]
        table = self.query_one("#data-table", DataTable)
        cursor_col = table.cursor_column
        if not (0 <= cursor_col < len(columns)) or columns[cursor_col] != agg_col:
            return False
        cells = self.table_page["data"][agg_col]
        if not (0 <= table.cursor_row < len(cells)):
            return False
        pos = int(cells[table.cursor_row])
        if pos < 0:
            self.notify(f"No {op} row for this group (no non-null values)", severity="warning")
            return True
        self._clear_group()  # back to the base table...
        self._go_to_row_col(pos, value_col)  # ...landing on the extreme row/column
        return True

    def action_plot_column(self) -> None:
        """p key — plot a downsampled overview of the whole cursor column."""
        if not self._in_data_grid():
            return
        if PlotextPlot is None:
            self.notify("Plotting needs the 'textual-plotext' package", severity="warning")
            return
        if self.browser is not None and self.browser.get_group(self.selected_path):
            bars = self.browser.group_bars(self.selected_path)
            if not bars["values"]:
                self.notify("Nothing to plot for this group-by", severity="warning")
                return
            title = f"{self.selected_path} · {bars['agg']} by {bars['key']}"
            self.push_screen(GroupBarScreen(title_prefix=title, bars=bars))
            return
        buffer = self.table_buffer or self.table_page
        columns = buffer["columns"]
        if not columns:
            return
        cursor_col = self.query_one("#data-table", DataTable).cursor_column
        name = columns[min(max(0, cursor_col), len(columns) - 1)]
        # Cheap numeric check on the already-loaded buffer; this also rejects
        # expensive object columns before any whole-column strided read.
        sample = np.asarray(buffer["data"][name])
        if sample.dtype.kind not in "iufb":
            self.notify(f"Column {name!r} is not numeric", severity="warning")
            return

        column: str | int | None
        if buffer.get("source_kind") == "ctable":
            column = name
        elif name.isdigit():  # array grids label columns with global indices
            column = int(name)
        else:  # 1-D arrays (single navigable dim) have one "value" column
            column = None

        layout = self._data_layout

        def fetch(start: int, stop: int | None) -> dict:
            return self.browser.plot_series(
                self.selected_path,
                column=column,
                layout=layout,
                max_points=self._PLOT_MAX_POINTS,
                row_start=start,
                row_stop=stop,
            )

        def raw_fetch(start: int, stop: int | None) -> dict:
            # max_points caps the raw read: a wide range is strided-sampled
            # rather than refused, for the hi-res 'r' (raw values) view.
            return self.browser.read_series(
                self.selected_path,
                column=column,
                layout=layout,
                row_start=start,
                row_stop=stop,
                max_points=PlotScreen._HIRES_MAX_POINTS,
            )

        # Col-vs-col scatter ('s' key) is CTable-only; build the read closure and
        # the visible-column universe for the Y picker just for that source kind.
        scatter_fetch = None
        scatter_columns = None
        if buffer.get("source_kind") == "ctable":

            def scatter_fetch(ycol: str, start: int, stop: int) -> dict:
                return self.browser.read_xy(
                    self.selected_path,
                    xcol=name,
                    ycol=ycol,
                    layout=layout,
                    row_start=start,
                    row_stop=stop,
                    max_points=PlotScreen._SCATTER_MAX_POINTS,
                )

            scatter_columns = self.browser.column_names(self.selected_path) or columns

        series = fetch(0, None)  # whole series (uses the fast SUMMARY tier if any)
        x, _ymin, _ymax, _descr = _plot_view(series)
        if x.size == 0:
            self.notify(f"Column {name!r} has no finite values to plot", severity="warning")
            return
        self.push_screen(
            PlotScreen(
                title_prefix=f"{self.selected_path} · {name}",
                fetch=fetch,
                n=series["n"],
                row_start=series["row_start"],
                row_stop=series["row_stop"],
                series=series,
                raw_fetch=raw_fetch,
                xcol=name,
                scatter_fetch=scatter_fetch,
                scatter_columns=scatter_columns,
            ),
            self._view_plot_range,
        )

    def _view_plot_range(self, span: tuple[int, int] | None) -> None:
        """Lock the data grid to a row range chosen with 'v' in the plot modal.

        The grid is replaced in place with just the range, so paging cannot
        leave it (``esc`` unlocks).  CTable nodes use a zero-copy ``slice`` view;
        NDArray nodes narrow the layout's navigable row dim.  Other source kinds
        fall back to a cursor jump.
        """
        if span is None or self.table_page is None:
            return
        start, stop = span
        kind = self.table_page.get("source_kind")
        if kind == "ctable" and self.browser is not None:
            self._enter_row_window(start, stop, backend="ctable")
        elif kind in {"ndarray_slice", "ndarray2d"} and self._data_layout is not None:
            self._enter_row_window(start, stop, backend="ndarray")
        else:
            self._go_to_row(start)
            self.notify(f"Viewing rows {start}:{stop}")

    def _enter_row_window(self, start: int, stop: int, *, backend: str) -> None:
        """Replace the grid with a locked [start:stop] window (in place)."""
        if backend == "ctable":
            try:
                self.browser.set_row_window(self.selected_path, start, stop)
            except Exception as exc:  # pragma: no cover - defensive
                self.notify(f"Could not lock rows: {exc}", severity="error")
                return
        else:  # ndarray: narrow the navigable row dim, scroll back to its top
            self._data_layout.row_window = (start, stop)
            self._data_layout.row_start = 0
            self._data_layout.row_stop = 0
        self.row_window = (start, stop)
        self._reload_row_window(0)
        self.notify(f"Locked to rows {start}:{stop} · esc to unlock")

    def _exit_row_window(self) -> None:
        """Unlock the row window and restore the full grid."""
        if self.row_window is None:
            return
        if self.browser is not None:
            self.browser.clear_row_window(self.selected_path)
        if self._data_layout is not None and self._data_layout.row_window is not None:
            self._data_layout.row_window = None
            self._data_layout.row_start = 0
            self._data_layout.row_stop = 0
        self.row_window = None
        self._reload_row_window(0)

    def _reload_row_window(self, start: int) -> None:
        """Rebuild the data grid from scratch after a window change."""
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, start)
        self._update_data_table(data, cursor_row=0, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def action_go_to_column(self) -> None:
        if not self._in_data_grid():
            return
        page = self.table_page
        if page.get("source_kind") not in _COL_PAGED_KINDS:
            return
        if page["source_kind"] == "ctable":
            # Named columns -> pick from a searchable list (type to filter, ↑/↓,
            # Enter); the option ids index into the visible-column universe, which
            # is exactly what _go_to_column expects.
            names = self.browser.column_names(self.selected_path)
            screen: ModalScreen[int | None] = ColumnSelectScreen(names=names, title="Go to column")
        else:
            # N-D arrays have no column names; fall back to a numeric index entry.
            current = page["col_start"] + self.query_one("#data-table", DataTable).cursor_column
            screen = GoToColumnScreen(ncols=page["ncols"], current=current, names=None)
        self.push_screen(screen, self._go_to_column)

    def action_filter_rows(self) -> None:
        if not self._in_data_grid():
            return
        if self.table_page.get("source_kind") != "ctable":
            self.notify("Filtering is only supported for CTable nodes", severity="warning")
            return
        if self.browser.get_group(self.selected_path):
            self.notify("Ungroup (Esc) before filtering", severity="warning")
            return
        screen = FilterScreen(current=self.browser.get_filter(self.selected_path))
        self.push_screen(screen, self._apply_filter)

    def action_sort_rows(self) -> None:
        if not self._in_data_grid():
            return
        if self.table_page.get("source_kind") != "ctable":
            self.notify("Sorting is only supported for CTable nodes", severity="warning")
            return
        if self.browser.get_group(self.selected_path):
            # Sort the (tiny) grouped result by any of its columns — key or aggregate.
            columns = self.browser.column_names(self.selected_path) or []
            screen = SortByScreen(
                columns=columns,
                current=self.browser.get_group_sort(self.selected_path),
                title="Sort grouped result (Enter applies, R reverses)",
            )
            self.push_screen(screen, self._apply_group_sort)
            return
        # Offer every column. FULL-indexed columns reuse their pre-sorted
        # positions (instant); the rest materialise the sort key and lexsort on
        # demand — slower on a big table, but no whole-table copy.
        columns = self.browser.column_names(self.selected_path) or []
        if not columns:
            self.notify("No columns to sort by", severity="warning")
            return
        indexed = set(self.browser.full_index_columns(self.selected_path))
        title = (
            "Sort by column (Enter applies, R reverses)\n◆ = indexed (fast)"
            if indexed
            else "Sort by column (Enter applies, R reverses)\nnon-indexed columns are scanned"
        )
        labels = [f"◆ {c}" if c in indexed else c for c in columns]
        screen = SortByScreen(
            columns=columns,
            labels=labels,
            current=self.browser.get_sort(self.selected_path),
            title=title,
        )
        self.push_screen(screen, self._apply_sort)

    def _apply_sort(self, choice: tuple[str, bool] | None, *, reposition: bool = True) -> None:
        if choice is None or self.browser is None or self.table_page is None:
            return  # cancelled
        column, reverse = choice
        # A FULL-indexed column reuses its pre-sorted positions instantly; a
        # non-indexed column must materialise the key and lexsort, which can take
        # a while on a big table — run that off the UI thread with a spinner so
        # the app stays responsive.  An active filter forces the scan path too:
        # the index can't be reused over a filtered (where) view.
        indexed = column in set(self.browser.full_index_columns(self.selected_path))
        if indexed and not self.browser.get_filter(self.selected_path):
            try:
                self.browser.set_sort(self.selected_path, column, reverse)
            except Exception as exc:
                self.notify(f"Cannot sort: {exc}", severity="error")
                return
            self._finish_sort(column, reverse, reposition)
        else:
            self._sort_in_background(column, reverse, reposition)

    @work(thread=True, exclusive=True)
    def _sort_in_background(self, column: str, reverse: bool, reposition: bool) -> None:
        """Build the sort permutation off the UI thread, with a loading spinner."""
        table = self.query_one("#data-table", DataTable)
        self.app.call_from_thread(setattr, table, "loading", True)
        try:
            self.browser.set_sort(self.selected_path, column, reverse)
        except Exception as exc:
            self.app.call_from_thread(setattr, table, "loading", False)
            self.app.call_from_thread(self.notify, f"Cannot sort: {exc}", severity="error")
            return
        self.app.call_from_thread(self._finish_sort, column, reverse, reposition)

    def _finish_sort(self, column: str, reverse: bool, reposition: bool) -> None:
        """Repaint the grid after the sort view is in place (UI thread)."""
        self.query_one("#data-table", DataTable).loading = False
        self.row_window = None  # set_sort drops any window/filter; keep the chip in sync
        self.table_buffer = None
        # Park the cursor on the sorted column's first row.  On the initial sort
        # ('S') bring the column into view, clamping the window start to the tail
        # ('End') position so the natural left-to-right column order is preserved
        # and the last column shows a full window, not a lone column.  On reverse
        # ('R') leave the horizontal scroll where it is — only the order flips.
        names = self.browser.column_names(self.selected_path) or []
        col_idx = names.index(column) if column in names else 0
        if reposition:
            self.grid_col_start = min(col_idx, self._fit_col_start_backward(self.table_page["ncols"]))
        data = self._load_table_page(self.selected_path, 0)
        cursor_col = max(0, col_idx - data.get("col_start", 0))
        self._update_data_table(data, cursor_row=0, cursor_col=cursor_col)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def action_group_rows(self) -> None:
        if not self._in_data_grid():
            return
        if self.table_page.get("source_kind") != "ctable":
            self.notify("Grouping is only supported for CTable nodes", severity="warning")
            return
        keys = self.browser.group_key_columns(self.selected_path)
        if not keys:
            self.notify("No dictionary/numeric columns to group by", severity="warning")
            return
        # Pre-fill with this table's active group, else the last one used
        # anywhere; GroupByScreen ignores any field whose column is absent here.
        screen = GroupByScreen(
            keys=keys,
            values=self.browser.group_value_columns(self.selected_path),
            current=self.browser.get_group(self.selected_path) or self._last_group,
        )
        self.push_screen(screen, self._apply_group)

    def _apply_group(self, choice: tuple[str, str, str | None] | None) -> None:
        if choice is None or self.browser is None or self.table_page is None:
            return
        key, op, value_col = choice
        try:
            self.browser.set_group(self.selected_path, key, op, value_col)
        except Exception as exc:
            self.notify(f"Cannot group: {exc}", severity="error")
            return
        self._last_group = (key, op, value_col)  # remember for the next 'G' anywhere
        self.row_window = None  # set_group drops any window/filter; keep chips in sync
        self.table_buffer = None
        self.grid_col_start = 0
        data = self._load_table_page(self.selected_path, 0)
        # Park the cursor on the aggregate column (last column of the result).
        agg = self.browser.group_agg_column(self.selected_path)
        cols = data.get("columns", [])
        cursor_col = cols.index(agg) if agg in cols else 0
        self._update_data_table(data, cursor_row=0, cursor_col=cursor_col)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _apply_group_sort(self, choice: tuple[str, bool] | None) -> None:
        if choice is None or self.browser is None or self.table_page is None:
            return  # cancelled
        column, reverse = choice
        try:
            self.browser.set_group_sort(self.selected_path, column, reverse)
        except Exception as exc:
            self.notify(f"Cannot sort: {exc}", severity="error")
            return
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, 0)
        cols = data.get("columns", [])
        cursor_col = cols.index(column) if column in cols else 0
        self._update_data_table(data, cursor_row=0, cursor_col=cursor_col)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def action_reverse_sort(self) -> None:
        """Flip ascending/descending on the currently sorted column."""
        if (
            not self._in_data_grid()
            or self.table_page.get("source_kind") != "ctable"
            or self.browser is None
        ):
            return
        # While grouped, 'R' flips the sort on the grouped result.
        gsort = self.browser.get_group_sort(self.selected_path)
        if self.browser.get_group(self.selected_path) is not None:
            if gsort is not None:
                self._apply_group_sort((gsort[0], not gsort[1]))
            return
        sort = self.browser.get_sort(self.selected_path)
        if sort is None:
            return
        column, reverse = sort
        self._apply_sort((column, not reverse), reposition=False)

    def _clear_sort(self) -> None:
        """Escape out of a sort view, restoring original row order."""
        self.browser.clear_sort(self.selected_path)
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, 0)
        self._update_data_table(data, cursor_row=0, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _clear_group(self) -> None:
        """Escape out of a group-by view, restoring the base table."""
        self.browser.clear_group(self.selected_path)
        self.table_buffer = None
        self.grid_col_start = 0
        data = self._load_table_page(self.selected_path, 0)
        self._update_data_table(data, cursor_row=0, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _clear_group_sort(self) -> None:
        """Drop the sort on a grouped result, keeping the group itself."""
        self.browser.clear_group_sort(self.selected_path)
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, 0)
        self._update_data_table(data, cursor_row=0, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _apply_filter(self, expr: str | None) -> None:
        if expr is None or self.browser is None or self.table_page is None:
            return
        if expr == (self.browser.get_filter(self.selected_path) or ""):
            return
        try:
            self.browser.set_filter(self.selected_path, expr)
        except Exception as exc:
            self.notify(f"Invalid filter: {exc}", severity="error")
            return
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, 0)
        self._update_data_table(data, cursor_row=0, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def action_filter_columns(self) -> None:
        if not self._in_data_grid():
            return
        if self.table_page.get("source_kind") != "ctable":
            self.notify("Column filtering is only supported for CTable nodes", severity="warning")
            return
        if self.browser.get_group(self.selected_path):
            self.notify("Ungroup (Esc) before filtering columns", severity="warning")
            return
        # Preselect the currently-shown columns (column_names honors any active
        # selection); the picker universe is the full, unfiltered column set.
        screen = ColumnFilterScreen(
            names=self.browser.base_column_names(self.selected_path),
            selected=self.browser.column_names(self.selected_path) or [],
        )
        self.push_screen(screen, self._apply_column_selection)

    def _apply_column_selection(self, names: list[str] | None) -> None:
        if names is None or self.browser is None or self.table_page is None:
            return  # cancelled
        self.browser.set_column_selection(self.selected_path, names)
        self._reload_columns()

    def _reload_columns(self) -> None:
        """Reload the grid after the visible-column set changed."""
        self.grid_col_start = 0
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, self.table_page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _go_to_column(self, col: int | None) -> None:
        if col is None or self.table_page is None:
            return
        self.grid_col_start = col
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, self.table_page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row, cursor_col=0)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _focused_pane(self):
        focused = self.focused
        if focused is None:
            return None
        for selector in ("#tree-pane", "#meta-pane", "#vlmeta-pane", "#data-pane"):
            pane = self.query_one(selector, Vertical)
            if focused is pane or pane in focused.ancestors:
                return pane
        return None

    def action_maximize_panel(self) -> None:
        pane = self._focused_pane()
        if pane is None:
            self.notify("Focus a pane before maximizing", severity="warning")
            return
        if self.screen.maximize(pane, container=False):
            self.call_after_refresh(self._reload_table_for_current_viewport)

    def action_restore_or_refresh(self) -> None:
        if self.screen.maximized is not None:
            self.screen.maximized = None
            self.call_after_refresh(self._reload_table_for_current_viewport)
            return
        self.action_refresh()

    def _ensure_viewport_consistent(self) -> None:
        """Reload the page if it was sized before the layout had settled.

        The first page of a node may be loaded while the DataTable still has
        no size, in which case the CLI fallbacks (preview_rows/preview_cols)
        determine the window.  Later paging then uses the settled viewport
        sizes, so the windows would drift unless we reload once here.
        """
        page = self.table_page
        if page is None or not self.query_one("#data-table-row", Horizontal).display:
            return
        rows_loaded = page["stop"] - page["start"]
        rows_want = min(self._table_page_size(), page["nrows"] - page["start"])
        cols_ok = True
        if page.get("source_kind") in _COL_PAGED_KINDS:
            # The column window is fitted to the width current at load time
            cols_ok = page.get("viewport_width") == self._data_table_width()
        if rows_loaded == rows_want and cols_ok:
            return
        self._reload_table_for_current_viewport()

    def _on_data_table_resized(self) -> None:
        self.call_after_refresh(self._ensure_viewport_consistent)

    def _reload_table_for_current_viewport(self) -> None:
        """Reload the table page after layout changes such as maximize/restore."""
        if self.table_page is None or not self.query_one("#data-table-row", Horizontal).display:
            return
        current = self.table_page["start"] + self.query_one("#data-table", DataTable).cursor_row
        page_size = self._table_page_size()
        start = (current // page_size) * page_size
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, start)
        self._update_data_table(data, cursor_row=current - data["start"])
        self._update_data_header(data)

    def _go_to_row(self, row: int | None) -> None:
        if row is None or self.table_page is None:
            return
        page_size = self._table_page_size()
        start = (row // page_size) * page_size
        data = self._load_table_page(self.selected_path, start)
        self._update_data_table(data, cursor_row=row - data["start"])
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def _go_to_row_col(self, row: int, column: str) -> None:
        """Jump to *row* with the cursor on *column*, scrolling it into view."""
        if self.table_page is None:
            return
        names = self.browser.column_names(self.selected_path) or []
        col_idx = names.index(column) if column in names else 0
        self.grid_col_start = min(col_idx, self._fit_col_start_backward(self.table_page["ncols"]))
        page_size = self._table_page_size()
        start = (row // page_size) * page_size
        data = self._load_table_page(self.selected_path, start)
        cursor_col = max(0, col_idx - data.get("col_start", 0))
        self._update_data_table(data, cursor_row=row - data["start"], cursor_col=cursor_col)
        self._update_data_header(data)
        self.query_one("#data-table", DataTable).focus()

    def action_refresh(self) -> None:
        tree = self.query_one("#tree", Tree)
        node = tree.cursor_node or tree.root
        self.loaded_paths.discard(node.data or "/")
        node.remove_children()
        self.load_children(node)
        self.update_panels(node.data or "/")

    def _adjust_fixed_value(self, direction: int) -> None:
        """Adjust the fixed value of the active dimension (if it is fixed).

        The value clamps at the boundaries (no wrap-around).
        """
        layout = self._data_layout
        if layout is None or self.table_page is None:
            return
        dim = self._active_dim
        if dim not in layout.fixed_values:
            return
        total = layout.shape[dim]
        if total <= 0:
            return
        current = layout.fixed_values[dim]
        new_val = min(max(current + direction, 0), total - 1)
        if new_val == current:
            return
        new_fixed = dict(layout.fixed_values)
        new_fixed[dim] = new_val
        self._data_layout = layout.copy_with(fixed_values=new_fixed)
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, self.table_page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row)
        self._update_data_header(data)

    def _rebuild_layout(self, navigable: list[int]) -> DataSliceLayout:
        """Return a copy of the current layout with the given *navigable* dims.

        All non-navigable dimensions are fixed at their previous value (or 0).
        """
        layout = self._data_layout
        if layout is None:
            raise RuntimeError("No layout available")
        new_fixed: dict[int, int] = {}
        for d in range(len(layout.shape)):
            if d in navigable:
                continue
            if d in layout.fixed_values:
                new_fixed[d] = layout.fixed_values[d]
            else:
                new_fixed[d] = 0
        return layout.copy_with(fixed_values=new_fixed, navigable_dims=navigable)

    def _dim_toggle(self) -> None:
        """: key — toggle active dim between fixed and navigable."""
        layout = self._data_layout
        if layout is None or self.table_page is None:
            return
        dim = self._active_dim
        if dim not in range(len(layout.shape)):
            return

        if dim in layout.navigable_dims:
            # Navigable → fixed (at index 0)
            new_nav = [d for d in layout.navigable_dims if d != dim]
            self._data_layout = self._rebuild_layout(new_nav)
        elif dim in layout.fixed_values:
            # Fixed → navigable (if room)
            if len(layout.navigable_dims) >= 2:
                self.notify("At most 2 navigable dimensions are allowed")
                return
            new_nav = sorted(layout.navigable_dims + [dim])
            self._data_layout = self._rebuild_layout(new_nav)
        else:
            return

        # Refresh the display (DataTable for 1-2 nav dims, same path for 0)
        self.table_buffer = None
        data = self._load_table_page(self.selected_path, self.table_page["start"])
        cursor_row = self.query_one("#data-table", DataTable).cursor_row
        self._update_data_table(data, cursor_row=cursor_row)
        self._update_data_header(data)

    def _dim_cursor(self, direction: int) -> None:
        """In dim mode, move the active dimension up (+1) or down (-1)."""
        layout = self._data_layout
        if layout is None or len(layout.shape) < 1:
            return
        ndim = len(layout.shape)
        self._active_dim = (self._active_dim + direction) % ndim
        if self.table_page is not None:
            self._update_data_header(self.table_page)

    def _dim_adjust(self, direction: int) -> None:
        """In DIM mode, adjust the active dim: fixed value or navigable viewport."""
        layout = self._data_layout
        if layout is None or self.table_page is None:
            return
        dim = self._active_dim
        if dim in layout.fixed_values:
            self._adjust_fixed_value(direction)
        elif dim in layout.navigable_dims:
            self._scroll_navigable_viewport(direction)

    def _scroll_navigable_viewport(self, direction: int) -> None:
        """Shift the viewport of a navigable dimension by one step (clamps)."""
        layout = self._data_layout
        if layout is None or self.table_page is None:
            return
        dim = self._active_dim
        if dim not in layout.navigable_dims:
            return

        pos = layout.navigable_dims.index(dim)
        page = self.table_page
        total = layout.shape[dim]

        if pos == 0:
            # Row navigable dim — shift start by one row, keeping a full page
            max_start = max(0, total - self._table_page_size())
            new_start = min(max(page["start"] + direction, 0), max_start)
            if new_start == page["start"]:
                return
            self.table_buffer = None
            data = self._load_table_page(self.selected_path, new_start)
        else:
            # Column navigable dim — shift col_start by one whole column
            max_col = self._fit_col_start_backward(total)
            new_col = min(max(page["col_start"] + direction, 0), max_col)
            if new_col == page["col_start"]:
                return
            self.grid_col_start = new_col
            self.table_buffer = None
            data = self._load_table_page(self.selected_path, page["start"])

        self._update_data_table(data)
        self._update_data_header(data)

    def action_dim_cycle(self) -> None:
        """d key — toggle DIM mode on/off."""
        if not self._in_data_grid():
            return
        layout = self._data_layout
        if layout is None or len(layout.shape) < 1:
            self.notify("No dimensions to navigate")
            return

        self._dim_mode = not self._dim_mode
        if self.table_page is not None:
            self._update_data_header(self.table_page)

    def action_dim_toggle_nav(self) -> None:
        """Enter — toggle active dim between fixed and navigable (in dim mode)."""
        if not self._in_data_grid() or not self._dim_mode:
            return
        self._dim_toggle()

    def action_dim_exit(self) -> None:
        """Escape: exit dim mode, unlock a row window, or clear a CTable filter.

        One layer per press, peeling derived views before their base: dim mode,
        then the locked row window, group-sort, group, then the sort (built on
        the filter), then the row filter, then the column filter.
        """
        if self._dim_mode:
            self._dim_mode = False
            if self.table_page is not None:
                self._update_data_header(self.table_page)
            return
        if self.row_window is not None:
            self._exit_row_window()
            return
        if (
            not self._in_data_grid()
            or self.table_page.get("source_kind") != "ctable"
            or self.browser is None
        ):
            return
        if self.browser.get_group_sort(self.selected_path):
            self._clear_group_sort()
        elif self.browser.get_group(self.selected_path):
            self._clear_group()
        elif self.browser.get_sort(self.selected_path):
            self._clear_sort()
        elif self.browser.get_filter(self.selected_path):
            self._apply_filter("")
        elif self.browser.get_column_filter(self.selected_path):
            self.browser.set_column_selection(self.selected_path, None)
            self._reload_columns()

    def action_grid_row_top(self) -> None:
        """Jump to the first row of the table."""
        if not self._in_data_grid():
            return
        self._go_to_row(0)

    def action_grid_row_bottom(self) -> None:
        """Jump to the last row of the table."""
        if not self._in_data_grid():
            return
        self._go_to_row(self.table_page["nrows"] - 1)

    def action_show_help(self) -> None:
        self.push_screen(HelpScreen())

    def action_grid_col_start(self) -> None:
        """Jump to the first column window (alias of Home)."""
        if self._in_data_grid():
            self._grid_col_home()

    def action_grid_col_end(self) -> None:
        """Jump to the last column window (alias of End)."""
        if self._in_data_grid():
            self._grid_col_end()
