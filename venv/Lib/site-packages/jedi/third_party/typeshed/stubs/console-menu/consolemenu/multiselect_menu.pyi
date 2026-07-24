from consolemenu import ConsoleMenu as ConsoleMenu
from consolemenu.console_menu import MenuItem
from consolemenu.items import SubmenuItem as SubmenuItem
from consolemenu.menu_formatter import MenuFormatBuilder

class MultiSelectMenu(ConsoleMenu):
    ack_item_completion: bool
    def __init__(
        self,
        title: str | None = None,
        subtitle: str | None = None,
        formatter: MenuFormatBuilder | None = None,
        prologue_text: str | None = None,
        epilogue_text: str | None = None,
        ack_item_completion: bool = True,
        show_exit_option: bool = True,
        exit_option_text: str = "Exit",
        clear_screen: bool = True,
    ) -> None: ...
    def append_item(self, item: MenuItem) -> None: ...
    current_option: int
    def process_user_input(self) -> None: ...
