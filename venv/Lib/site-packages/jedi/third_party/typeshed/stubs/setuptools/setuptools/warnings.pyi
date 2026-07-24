from typing_extensions import TypeAlias

_DueDate: TypeAlias = tuple[int, int, int]  # time tuple

class SetuptoolsWarning(UserWarning):
    @classmethod
    def emit(
        cls,
        summary: str | None = None,
        details: str | None = None,
        due_date: _DueDate | None = None,
        see_docs: str | None = None,
        see_url: str | None = None,
        stacklevel: int = 2,
        **kwargs,
    ) -> None: ...

class InformationOnly(SetuptoolsWarning): ...
class SetuptoolsDeprecationWarning(SetuptoolsWarning): ...
