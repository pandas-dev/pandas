"""
Progress bar
------------
"""

from __future__ import annotations

__lazy_modules__ = {"plumbum.cli.termsize", "warnings"}

import abc
import sys
import time
import warnings
from typing import TYPE_CHECKING, Any, Generic

from plumbum._compat.typing import TypeVar
from plumbum.cli.termsize import get_terminal_size

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from plumbum._compat.typing import Self


T = TypeVar("T", default=int)


class ProgressBase(Generic[T], metaclass=abc.ABCMeta):
    """Base class for progress bars. Customize for types of progress bars.

    :param iterator: The iterator to wrap with a progress bar
    :param length: The length of the iterator (will use ``__len__`` if None)
    :param timer: Try to time the completion status of the iterator
    :param body: True if the slow portion occurs outside the iterator (in a loop, for example)
    :param has_output: True if the iteration body produces output to the screen (forces rewrite off)
    :param clear: Clear the progress bar afterwards, if applicable.
    """

    __slots__ = (
        "_start_time",
        "_value",
        "body",
        "clear",
        "has_output",
        "iter",
        "iterator",
        "length",
        "timer",
    )

    def __init__(
        self,
        iterator: Iterable[T] | None = None,
        length: int | None = None,
        timer: bool = True,
        body: bool = False,
        has_output: bool = False,
        clear: bool = True,
    ):
        if length is None:
            length = len(iterator)  # type: ignore[arg-type]
        elif iterator is None:
            iterator = range(length)  # type: ignore[assignment]
        elif length is None and iterator is None:
            raise TypeError("Expected either an iterator or a length")

        assert iterator is not None

        self.length = length
        self.iterator = iterator
        self.timer = timer
        self.body = body
        self.has_output = has_output
        self.clear = clear
        self.iter: Iterator[T] | None = None
        self._start_time: float | None = None
        self._value: int | None = None

    def __len__(self) -> int:
        return self.length

    def __iter__(self) -> Self:
        self.start()
        return self

    @abc.abstractmethod
    def start(self) -> None:
        """This should initialize the progress bar and the iterator"""
        self.iter = iter(self.iterator)
        self.value = -1 if self.body else 0
        self._start_time = time.perf_counter()

    def __next__(self) -> T:
        if self.iter is None:
            raise ValueError("Iteration not started")
        try:
            rval = next(self.iter)
            self.increment()
        except StopIteration:
            self.done()
            raise
        return rval

    def next(self) -> T:
        return next(self)

    @property
    def value(self) -> int | None:
        """This is the current value, as a property so setting it can be customized"""
        return self._value

    @value.setter
    def value(self, val: int) -> None:
        self._value = val

    @abc.abstractmethod
    def display(self) -> None:
        """Called to update the progress bar"""

    def increment(self) -> None:
        """Sets next value and displays the bar"""
        if self.value is None:
            raise ValueError("Iteration not started")
        self.value += 1
        self.display()

    def time_remaining(
        self,
    ) -> tuple[float, float] | tuple[None, None]:
        """Get the time remaining for the progress bar, guesses"""
        if self.value is None or self.value < 1:
            return None, None
        if self._start_time is None:
            raise ValueError("Iteration not started")
        elapsed_time = time.perf_counter() - self._start_time
        time_each = elapsed_time / self.value
        time_remaining = time_each * (self.length - self.value)
        return elapsed_time, time_remaining

    def str_time_remaining(self) -> str:
        """Returns a string version of time remaining"""
        if self.value is None or self.value < 1:
            return "Starting..."

        completed, remaining = self.time_remaining()
        return f"{completed:.0f} completed, {remaining:.0f} remaining"

    @abc.abstractmethod
    def done(self) -> None:
        """Is called when the iterator is done."""

    @classmethod
    def range(cls, *value: int, **kargs: Any) -> Self:
        """Fast shortcut to create a range based progress bar, assumes work done in body"""
        return cls(range(*value), body=True, **kargs)  # type: ignore[arg-type]

    @classmethod
    def wrap(
        cls, iterator: Iterable[Any], length: int | None = None, **kargs: Any
    ) -> Self:
        """Shortcut to wrap an iterator that does not do all the work internally"""
        return cls(iterator, length, body=True, **kargs)


class Progress(ProgressBase[T]):
    __slots__ = ("width",)

    def start(self) -> None:
        super().start()
        self.display()

    def done(self) -> None:
        self.value = self.length
        self.display()
        if self.clear and not self.has_output:
            sys.stdout.write("\r" + len(str(self)) * " " + "\r")
        else:
            sys.stdout.write("\n")
        sys.stdout.flush()

    def __str__(self) -> str:
        width = get_terminal_size(default=(0, 0))[0]
        if not self.length:
            self.width = 0
            return "0/0 complete"

        percent = max(self.value or 0, 0) / self.length
        ending = " " + (
            self.str_time_remaining()
            if self.timer
            else f"{self.value} of {self.length} complete"
        )
        if width - len(ending) < 10 or self.has_output:
            self.width = 0

            if self.timer:
                return f"{percent:.0%} complete: {self.str_time_remaining()}"

            return f"{percent:.0%} complete"

        self.width = width - len(ending) - 2 - 1
        nstars = int(percent * self.width)
        pbar = "[" + "*" * nstars + " " * (self.width - nstars) + "]" + ending

        str_percent = f" {percent:.0%} "

        return (
            pbar[: self.width // 2 - 2]
            + str_percent
            + pbar[self.width // 2 + len(str_percent) - 2 :]
        )

    def display(self) -> None:
        disptxt = str(self)
        if self.width == 0 or self.has_output:
            sys.stdout.write(disptxt + "\n")
        else:
            sys.stdout.write("\r")
            sys.stdout.write(disptxt)
        sys.stdout.flush()


class ProgressIPy(ProgressBase[T]):  # pragma: no cover
    __slots__ = ("_box", "_label", "prog")
    HTMLBOX = '<div class="widget-hbox widget-progress"><div class="widget-label" style="display:block;">{0}</div></div>'

    # pylint: disable-next=keyword-arg-before-vararg
    def __init__(self, iterator: Iterable[T] | None = None, *args: Any, **kargs: Any):
        # Ipython gives warnings when using widgets about the API potentially changing
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                from ipywidgets import HTML, HBox, IntProgress
            except ImportError:  # Support IPython < 4.0
                from IPython.html.widgets import HTML, HBox, IntProgress

        super().__init__(iterator, *args, **kargs)
        self.prog = IntProgress(max=self.length)
        self._label = HTML()
        self._box = HBox((self.prog, self._label))

    def start(self) -> None:
        from IPython.display import display

        display(self._box)  # type: ignore[no-untyped-call]
        super().start()

    @property
    def value(self) -> int:
        """This is the current value, -1 allowed (automatically fixed for display)"""
        if self._value is None:
            return -1
        return self._value

    @value.setter
    def value(self, val: int) -> None:
        # Pylint false positive
        self._value = val  # pylint: disable=assigning-non-slot
        self.prog.value = max(val, 0)
        self.prog.description = f"{self.value / self.length:.2%}"
        if self.timer and val > 0:
            self._label.value = self.HTMLBOX.format(self.str_time_remaining())

    def display(self) -> None:
        pass

    def done(self) -> None:
        if self.clear:
            self._box.close()


class ProgressAuto(ProgressBase[T]):  # pylint: disable=abstract-method
    """Automatically selects the best progress bar (IPython HTML or text). Does not work with qtconsole
    (as that is correctly identified as identical to notebook, since the kernel is the same); it will still
    iterate, but no graphical indication will be displayed.

    :param iterator: The iterator to wrap with a progress bar
    :param length: The length of the iterator (will use ``__len__`` if None)
    :param timer: Try to time the completion status of the iterator
    :param body: True if the slow portion occurs outside the iterator (in a loop, for example)
    """

    __slots__ = ()

    def __new__(cls, *args: Any, **kargs: Any) -> ProgressIPy[T] | Progress[T]:  # type: ignore[misc]
        """Uses the generator trick that if a cls instance is returned, the __init__ method is not called."""
        try:  # pragma: no cover
            __IPYTHON__  # type: ignore[name-defined] # noqa: B018
            from traitlets import TraitError

            try:
                return ProgressIPy(*args, **kargs)
            except TraitError:
                raise NameError() from None
        except (NameError, ImportError):
            return Progress(*args, **kargs)


ProgressAuto.register(ProgressIPy)
ProgressAuto.register(Progress)


def main() -> None:
    tst = Progress.range(20)
    for _ in tst:
        time.sleep(1)


__all__ = [
    "Progress",
    "ProgressAuto",
    "ProgressBase",
    "ProgressIPy",
    "main",
]


def __dir__() -> list[str]:
    return list(__all__)


if __name__ == "__main__":
    main()
