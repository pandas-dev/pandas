import re
from collections.abc import Iterable, KeysView
from typing import Literal, overload

def is_ascii(s: str) -> bool: ...

class IniFile:
    defaultGroup: str
    fileExtension: str
    filename: str
    tainted: bool
    content: dict[str, dict[str, str]]
    warnings: list[str]
    errors: list[str]
    def __init__(self, filename: str | None = None) -> None: ...
    def __cmp__(self, other: IniFile) -> bool: ...
    def parse(self, filename: str, headers: Iterable[str] | None = None) -> None: ...
    @overload
    def get(
        self,
        key: str,
        group: str | None = None,
        locale: bool = False,
        type: Literal["string"] = "string",
        list: Literal[False] = False,
        strict: bool = False,
    ) -> str: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["string"], list: Literal[True], strict: bool = False
    ) -> list[str]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        type: Literal["string"] = "string",
        strict: bool = False,
    ) -> list[str]: ...
    @overload
    def get(
        self,
        key: str,
        group: str | None,
        locale: bool,
        type: Literal["boolean"],
        list: Literal[False] = False,
        strict: bool = False,
    ) -> bool: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["boolean"],
        group: str | None = None,
        locale: bool = False,
        list: Literal[False] = False,
        strict: bool = False,
    ) -> bool: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["boolean"], list: Literal[True], strict: bool = False
    ) -> list[bool]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["boolean"],
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        strict: bool = False,
    ) -> list[bool]: ...
    @overload
    def get(
        self,
        key: str,
        group: str | None,
        locale: bool,
        type: Literal["integer"],
        list: Literal[False] = False,
        strict: bool = False,
    ) -> int: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["integer"],
        group: str | None = None,
        locale: bool = False,
        list: Literal[False] = False,
        strict: bool = False,
    ) -> int: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["integer"], list: Literal[True], strict: bool = False
    ) -> list[int]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["integer"],
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        strict: bool = False,
    ) -> list[int]: ...

    # Float
    @overload
    def get(
        self,
        key: str,
        group: str | None,
        locale: bool,
        type: Literal["numeric"],
        list: Literal[False] = False,
        strict: bool = False,
    ) -> float: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["numeric"],
        group: str | None = None,
        locale: bool = False,
        list: Literal[False] = False,
        strict: bool = False,
    ) -> float: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["numeric"], list: Literal[True], strict: bool = False
    ) -> list[float]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["numeric"],
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        strict: bool = False,
    ) -> list[float]: ...

    # Regex
    @overload
    def get(
        self,
        key: str,
        group: str | None,
        locale: bool,
        type: Literal["regex"],
        list: Literal[False] = False,
        strict: bool = False,
    ) -> re.Pattern[str]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["regex"],
        group: str | None = None,
        locale: bool = False,
        list: Literal[False] = False,
        strict: bool = False,
    ) -> re.Pattern[str]: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["regex"], list: Literal[True], strict: bool = False
    ) -> list[re.Pattern[str]]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["regex"],
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        strict: bool = False,
    ) -> list[re.Pattern[str]]: ...
    # point
    @overload
    def get(
        self,
        key: str,
        group: str | None,
        locale: bool,
        type: Literal["point"],
        list: Literal[False] = False,
        strict: bool = False,
    ) -> tuple[int, int]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["point"],
        group: str | None = None,
        locale: bool = False,
        list: Literal[False] = False,
        strict: bool = False,
    ) -> tuple[int, int]: ...
    @overload
    def get(
        self, key: str, group: str | None, locale: bool, type: Literal["point"], list: Literal[True], strict: bool = False
    ) -> list[tuple[int, int]]: ...
    @overload
    def get(
        self,
        key: str,
        *,
        type: Literal["point"],
        list: Literal[True],
        group: str | None = None,
        locale: bool = False,
        strict: bool = False,
    ) -> list[tuple[int, int]]: ...
    def getList(self, string: str) -> list[str]: ...
    def validate(self, report: Literal["All", "Warnings", "Errors"] = "All") -> None: ...
    def checkGroup(self, group: str) -> None: ...
    def checkKey(self, key: str, value: str, group: str) -> None: ...
    def checkValue(
        self,
        key: str,
        value: str,
        type: Literal["string", "localestring", "boolean", "numeric", "integer", "regex", "point"] = "string",
        list: bool = False,
    ) -> None: ...
    def checkExtras(self) -> None: ...
    def checkBoolean(self, value: str) -> Literal[1, 2] | None: ...
    def checkNumber(self, value: str) -> Literal[1, 2] | None: ...
    def checkInteger(self, value: str) -> Literal[1] | None: ...
    def checkPoint(self, value: str) -> Literal[1] | None: ...
    def checkString(self, value: str) -> Literal[0, 1]: ...
    def checkRegex(self, value: str) -> Literal[1] | None: ...
    def write(self, filename: str | None = None, trusted: bool = False) -> None: ...
    def set(self, key: str, value: str, group: str | None = None, locale: bool = False) -> None: ...
    def addGroup(self, group: str) -> None: ...
    def removeGroup(self, group: str) -> bool: ...
    def removeKey(self, key: str, group: str | None = None, locales: bool = True) -> str: ...
    def groups(self) -> KeysView[str]: ...
    def hasGroup(self, group: str) -> bool: ...
    def hasKey(self, key: str, group: str | None = None) -> bool: ...
    def getFileName(self) -> str: ...
