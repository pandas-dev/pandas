from collections.abc import Callable
from gettext import NullTranslations
from typing import Any, Protocol, type_check_only
from typing_extensions import Self

@type_check_only
class _TranslationStringFactory(Protocol):
    def __call__(
        self,
        msgid: str | TranslationString,
        mapping: dict[str, Any] | None = ...,
        default: str | None = ...,
        context: str | None = ...,
    ) -> TranslationString: ...

@type_check_only
class _ChameleonTranslate(Protocol):
    def __call__(
        self,
        msgid: str | TranslationString,
        mapping: dict[str, Any] | None = ...,
        context: str | None = ...,
        target_language: str | None = ...,
        default: str | None = ...,
    ) -> TranslationString: ...

@type_check_only
class _TranslatorPolicy(Protocol):
    def __call__(self, translations: NullTranslations, tstring: str, domain: str | None, context: str | None) -> str: ...

@type_check_only
class _Translator(Protocol):
    def __call__(
        self,
        tstring: str | TranslationString,
        domain: str | None = None,
        mapping: dict[str, Any] | None = None,
        context: str | None = None,
    ) -> str: ...

@type_check_only
class _PluralizerPolicy(Protocol):
    def __call__(
        self, translations: NullTranslations, singular: str, plural: str, n: int, domain: str | None, context: str | None
    ) -> str: ...

@type_check_only
class _Pluralizer(Protocol):
    def __call__(
        self,
        singular: str | TranslationString,
        plural: str | TranslationString,
        n: int,
        domain: str | None = None,
        mapping: dict[str, Any] | None = None,
        context: str | None = None,
    ) -> str: ...

class TranslationString(str):
    __slots__ = ("domain", "context", "default", "mapping")
    domain: str | None
    context: str | None
    default: str
    mapping: dict[str, Any] | None
    def __new__(
        self,
        msgid: str | Self,
        domain: str | None = None,
        default: str | None = None,
        mapping: dict[str, Any] | None = None,
        context: str | None = None,
    ) -> Self: ...
    def __mod__(self, options: dict[str, Any]) -> TranslationString: ...  # type: ignore[override]
    def interpolate(self, translated: str | None = None) -> str: ...
    def __reduce__(self) -> tuple[type[Self], tuple[str, str | None, str, dict[str, Any], str | None]]: ...

def TranslationStringFactory(factory_domain: str) -> _TranslationStringFactory: ...
def ChameleonTranslate(translator: Callable[[TranslationString], str] | None) -> _ChameleonTranslate: ...
def ugettext_policy(translations: NullTranslations, tstring: str, domain: str | None, context: str | None) -> str: ...
def dugettext_policy(translations: NullTranslations, tstring: str, domain: str | None, context: str | None) -> str: ...
def Translator(translations: NullTranslations | None = None, policy: _TranslatorPolicy | None = None) -> _Translator: ...
def ungettext_policy(
    translations: NullTranslations, singular: str, plural: str, n: int, domain: str | None, context: str | None
) -> str: ...
def dungettext_policy(
    translations: NullTranslations, singular: str, plural: str, n: int, domain: str | None, context: str | None
) -> str: ...
def Pluralizer(translations: NullTranslations | None = None, policy: _PluralizerPolicy | None = None) -> _Pluralizer: ...
