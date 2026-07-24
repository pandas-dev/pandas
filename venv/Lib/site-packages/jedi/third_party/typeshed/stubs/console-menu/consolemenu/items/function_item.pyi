from _typeshed import Incomplete
from collections.abc import Callable, Mapping, Sequence
from typing import Any

from consolemenu.console_menu import ConsoleMenu
from consolemenu.items import ExternalItem as ExternalItem

class FunctionItem(ExternalItem):
    function: Callable[..., Any]
    args: Sequence[Any]
    kwargs: Mapping[str, Any]
    return_value: Incomplete | None
    def __init__(
        self,
        text: str,
        function: Callable[..., Any],
        args: Sequence[Any] | None = None,
        kwargs: Mapping[str, Any] | None = None,
        menu: ConsoleMenu | None = None,
        should_exit: bool = False,
        menu_char: str | None = None,
    ) -> None: ...
    def action(self) -> None: ...
    def clean_up(self) -> None: ...
    def get_return(self) -> Any | None: ...
