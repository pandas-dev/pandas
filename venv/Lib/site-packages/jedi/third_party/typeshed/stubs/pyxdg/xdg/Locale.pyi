from collections.abc import Iterable

regex: str

def expand_languages(languages: Iterable[str] | None = None) -> list[str]: ...
def update(language: str | None = None) -> None: ...

langs: list[str]
