CSI: str
OSC: str
BEL: str

def code_to_chars(code: int) -> str: ...
def set_title(title: str) -> str: ...
def clear_screen(mode: int = 2) -> str: ...
def clear_line(mode: int = 2) -> str: ...

class AnsiCodes:
    def __init__(self) -> None: ...

class AnsiCursor:
    def UP(self, n: int = 1) -> str: ...
    def DOWN(self, n: int = 1) -> str: ...
    def FORWARD(self, n: int = 1) -> str: ...
    def BACK(self, n: int = 1) -> str: ...
    def POS(self, x: int = 1, y: int = 1) -> str: ...

# All attributes in the following classes are string in instances and int in the class.
# We use str since that is the common case for users.
class AnsiFore(AnsiCodes):
    BLACK: str
    RED: str
    GREEN: str
    YELLOW: str
    BLUE: str
    MAGENTA: str
    CYAN: str
    WHITE: str
    RESET: str
    LIGHTBLACK_EX: str
    LIGHTRED_EX: str
    LIGHTGREEN_EX: str
    LIGHTYELLOW_EX: str
    LIGHTBLUE_EX: str
    LIGHTMAGENTA_EX: str
    LIGHTCYAN_EX: str
    LIGHTWHITE_EX: str

class AnsiBack(AnsiCodes):
    BLACK: str
    RED: str
    GREEN: str
    YELLOW: str
    BLUE: str
    MAGENTA: str
    CYAN: str
    WHITE: str
    RESET: str
    LIGHTBLACK_EX: str
    LIGHTRED_EX: str
    LIGHTGREEN_EX: str
    LIGHTYELLOW_EX: str
    LIGHTBLUE_EX: str
    LIGHTMAGENTA_EX: str
    LIGHTCYAN_EX: str
    LIGHTWHITE_EX: str

class AnsiStyle(AnsiCodes):
    BRIGHT: str
    DIM: str
    NORMAL: str
    RESET_ALL: str

Fore: AnsiFore
Back: AnsiBack
Style: AnsiStyle
Cursor: AnsiCursor
