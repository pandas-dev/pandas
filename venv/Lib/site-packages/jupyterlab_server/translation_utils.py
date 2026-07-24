# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""
Localization utilities to find available language packs and packages with
localization data.
"""

from __future__ import annotations

import gettext
import importlib
import json
import locale
import os
import re
import sys
import traceback
from functools import lru_cache
from re import Pattern
from typing import Any

import babel
from packaging.version import parse as parse_version

# See compatibility note on `group` keyword in https://docs.python.org/3/library/importlib.metadata.html#entry-points
if sys.version_info < (3, 10):  # pragma: no cover
    from importlib_metadata import entry_points
else:  # pragma: no cover
    from importlib.metadata import entry_points

# Entry points
JUPYTERLAB_LANGUAGEPACK_ENTRY = "jupyterlab.languagepack"
JUPYTERLAB_LOCALE_ENTRY = "jupyterlab.locale"

# Constants
DEFAULT_LOCALE = "en"
SYS_LOCALE = locale.getlocale()[0] or DEFAULT_LOCALE
LOCALE_DIR = "locale"
LC_MESSAGES_DIR = "LC_MESSAGES"
DEFAULT_DOMAIN = "jupyterlab"
L10N_SCHEMA_NAME = "@jupyterlab/translation-extension:plugin"
PY37_OR_LOWER = sys.version_info[:2] <= (3, 7)

# Pseudo language locale for in-context translation
PSEUDO_LANGUAGE = "ach_UG"

_default_schema_context = "schema"
_default_settings_context = "settings"
_lab_i18n_config = "jupyter.lab.internationalization"

# mapping of schema translatable string selectors to translation context
DEFAULT_SCHEMA_SELECTORS = {
    "properties/.*/title": _default_settings_context,
    "properties/.*/description": _default_settings_context,
    "definitions/.*/properties/.*/title": _default_settings_context,
    "definitions/.*/properties/.*/description": _default_settings_context,
    "title": _default_schema_context,
    "description": _default_schema_context,
    # JupyterLab-specific
    r"jupyter\.lab\.setting-icon-label": _default_settings_context,
    r"jupyter\.lab\.menus/.*/label": "menu",
    r"jupyter\.lab\.toolbars/.*/label": "toolbar",
}


@lru_cache
def _get_default_schema_selectors() -> dict[Pattern, str]:
    return {
        re.compile("^/" + pattern + "$"): context
        for pattern, context in DEFAULT_SCHEMA_SELECTORS.items()
    }


def _prepare_schema_patterns(schema: dict) -> dict[Pattern, str]:
    return {
        **_get_default_schema_selectors(),
        **{
            re.compile("^/" + selector + "$"): _default_schema_context
            for selector in schema.get(_lab_i18n_config, {}).get("selectors", [])
        },
    }


# --- Private process helpers
# ----------------------------------------------------------------------------
def _get_installed_language_pack_locales() -> tuple[dict[str, Any], str]:
    """
    Get available installed language pack locales.

    Returns
    -------
    tuple
        A tuple, where the first item is the result and the second item any
        error messages.
    """
    data = {}
    messages = []
    for entry_point in entry_points(group=JUPYTERLAB_LANGUAGEPACK_ENTRY):
        try:
            data[entry_point.name] = os.path.dirname(entry_point.load().__file__)
        except Exception:  # pragma: no cover
            messages.append(traceback.format_exc())

    message = "\n".join(messages)
    return data, message


def _get_installed_package_locales() -> tuple[dict[str, Any], str]:
    """
    Get available installed packages containing locale information.

    Returns
    -------
    tuple
        A tuple, where the first item is the result and the second item any
        error messages. The value for the key points to the root location
        the package.
    """
    data = {}
    messages = []
    for entry_point in entry_points(group=JUPYTERLAB_LOCALE_ENTRY):
        try:
            data[entry_point.name] = os.path.dirname(entry_point.load().__file__)
        except Exception:
            messages.append(traceback.format_exc())

    message = "\n".join(messages)
    return data, message


# --- Helpers
# ----------------------------------------------------------------------------
def is_valid_locale(locale_: str) -> bool:
    """
    Check if a `locale_` value is valid.

    Parameters
    ----------
    locale_: str
        Language locale code.

    Notes
    -----
    A valid locale is in the form language (See ISO-639 standard) and an
    optional territory (See ISO-3166 standard).

    Examples of valid locales:
    - English: DEFAULT_LOCALE
    - Australian English: "en_AU"
    - Portuguese: "pt"
    - Brazilian Portuguese: "pt_BR"

    Examples of invalid locales:
    - Australian Spanish: "es_AU"
    - Brazilian German: "de_BR"
    """
    # Add exception for Norwegian
    if locale_ in {
        "no_NO",
    }:
        return True

    valid = False
    try:
        babel.Locale.parse(locale_)
        valid = True
    except (babel.core.UnknownLocaleError, ValueError):
        # Expected error if the locale is unknown
        pass

    return valid


def get_display_name(locale_: str, display_locale: str = DEFAULT_LOCALE) -> str:
    """
    Return the language name to use with a `display_locale` for a given language locale.

    Parameters
    ----------
    locale_: str
        The language name to use.
    display_locale: str, optional
        The language to display the `locale_`.

    Returns
    -------
    str
        Localized `locale_` and capitalized language name using `display_locale` as language.
    """
    locale_ = locale_ if is_valid_locale(locale_) else DEFAULT_LOCALE
    display_locale = display_locale if is_valid_locale(display_locale) else DEFAULT_LOCALE
    try:
        loc = babel.Locale.parse(locale_)
        display_name = loc.get_display_name(display_locale)
    except babel.UnknownLocaleError:
        display_name = display_locale
    if display_name:
        display_name = display_name[0].upper() + display_name[1:]
    return display_name  # type:ignore[return-value]


def merge_locale_data(
    language_pack_locale_data: dict[str, Any], package_locale_data: dict[str, Any]
) -> dict[str, Any]:
    """
    Merge language pack data with locale data bundled in packages.

    Parameters
    ----------
    language_pack_locale_data: dict
        The dictionary with language pack locale data.
    package_locale_data: dict
        The dictionary with package locale data.

    Returns
    -------
    dict
        Merged locale data.
    """
    result = language_pack_locale_data
    package_lp_metadata = language_pack_locale_data.get("", {})
    package_lp_version = package_lp_metadata.get("version", None)
    package_lp_domain = package_lp_metadata.get("domain", None)

    package_metadata = package_locale_data.get("", {})
    package_version = package_metadata.get("version", None)
    package_domain = package_metadata.get("domain", "None")

    if package_lp_version and package_version and package_domain == package_lp_domain:
        package_version = parse_version(package_version)
        package_lp_version = parse_version(package_lp_version)

        if package_version > package_lp_version:
            # If package version is more recent, then update keys of the language pack
            result = language_pack_locale_data.copy()
            result.update(package_locale_data)

    return result


def get_installed_packages_locale(locale_: str) -> tuple[dict, str]:
    """
    Get all jupyterlab extensions installed that contain locale data.

    Returns
    -------
    tuple
        A tuple in the form `(locale_data_dict, message)`,
        where the `locale_data_dict` is an ordered list
        of available language packs:
            >>> {"package-name": locale_data, ...}

    Examples
    --------
    - `entry_points={"jupyterlab.locale": "package-name = package_module"}`
    - `entry_points={"jupyterlab.locale": "jupyterlab-git = jupyterlab_git"}`
    """
    found_package_locales, message = _get_installed_package_locales()
    packages_locale_data = {}
    messages = message.split("\n")
    if not message:
        for package_name, package_root_path in found_package_locales.items():
            locales = {}
            try:
                locale_path = os.path.join(package_root_path, LOCALE_DIR)
                # Handle letter casing
                locales = {
                    loc.lower(): loc
                    for loc in os.listdir(locale_path)
                    if os.path.isdir(os.path.join(locale_path, loc))
                }
            except Exception:
                messages.append(traceback.format_exc())

            if locale_.lower() in locales:
                locale_json_path = os.path.join(
                    locale_path,
                    locales[locale_.lower()],
                    LC_MESSAGES_DIR,
                    f"{package_name}.json",
                )
                if os.path.isfile(locale_json_path):
                    try:
                        with open(locale_json_path, encoding="utf-8") as fh:
                            packages_locale_data[package_name] = json.load(fh)
                    except Exception:
                        messages.append(traceback.format_exc())

    return packages_locale_data, "\n".join(messages)


# --- API
# ----------------------------------------------------------------------------
def get_language_packs(display_locale: str = DEFAULT_LOCALE) -> tuple[dict, str]:
    """
    Return the available language packs installed in the system.

    The returned information contains the languages displayed in the current
    locale.

    Parameters
    ----------
    display_locale: str, optional
        Default is DEFAULT_LOCALE.

    Returns
    -------
    tuple
        A tuple in the form `(locale_data_dict, message)`.
    """
    found_locales, message = _get_installed_language_pack_locales()
    locales = {}
    messages = message.split("\n")
    if not message:
        invalid_locales = []
        valid_locales = []
        messages = []
        for locale_ in found_locales:
            if is_valid_locale(locale_):
                valid_locales.append(locale_)
            else:
                invalid_locales.append(locale_)

        display_locale_ = display_locale if display_locale in valid_locales else DEFAULT_LOCALE
        locales = {
            DEFAULT_LOCALE: {
                "displayName": (
                    get_display_name(DEFAULT_LOCALE, display_locale_)
                    if display_locale != PSEUDO_LANGUAGE
                    else "Default"
                ),
                "nativeName": get_display_name(DEFAULT_LOCALE, DEFAULT_LOCALE),
            }
        }
        for locale_ in valid_locales:
            locales[locale_] = {
                "displayName": get_display_name(locale_, display_locale_),
                "nativeName": get_display_name(locale_, locale_),
            }

        if invalid_locales:
            if PSEUDO_LANGUAGE in invalid_locales:
                invalid_locales.remove(PSEUDO_LANGUAGE)
                locales[PSEUDO_LANGUAGE] = {
                    "displayName": "Pseudo-language",
                    # Trick to ensure the proper language is selected in the language menu
                    "nativeName": (
                        "to translate the UI"
                        if display_locale != PSEUDO_LANGUAGE
                        else "Pseudo-language"
                    ),
                }
            # Check again as the pseudo-language was maybe the only invalid locale
            if invalid_locales:
                messages.append(f"The following locales are invalid: {invalid_locales}!")

    return locales, "\n".join(messages)


def get_language_pack(locale_: str) -> tuple:
    """
    Get a language pack for a given `locale_` and update with any installed
    package locales.

    Returns
    -------
    tuple
        A tuple in the form `(locale_data_dict, message)`.

    Notes
    -----
    We call `_get_installed_language_pack_locales` via a subprocess to
    guarantee the results represent the most up-to-date entry point
    information, which seems to be defined on interpreter startup.
    """
    found_locales, message = _get_installed_language_pack_locales()
    found_packages_locales, message = get_installed_packages_locale(locale_)
    locale_data = {}
    messages = message.split("\n")
    if (
        not message
        and (locale_ == PSEUDO_LANGUAGE or is_valid_locale(locale_))
        and locale_ in found_locales
    ):
        path = found_locales[locale_]
        for root, __, files in os.walk(path, topdown=False):
            for name in files:
                if name.endswith(".json"):
                    pkg_name = name.replace(".json", "")
                    json_path = os.path.join(root, name)
                    try:
                        with open(json_path, encoding="utf-8") as fh:
                            merged_data = json.load(fh)
                    except Exception:
                        messages.append(traceback.format_exc())

                    # Load packages with locale data and merge them
                    if pkg_name in found_packages_locales:
                        pkg_data = found_packages_locales[pkg_name]
                        merged_data = merge_locale_data(merged_data, pkg_data)

                    locale_data[pkg_name] = merged_data

        # Check if package locales exist that do not exists in language pack
        for pkg_name, data in found_packages_locales.items():
            if pkg_name not in locale_data:
                locale_data[pkg_name] = data

    return locale_data, "\n".join(messages)


# --- Translators
# ----------------------------------------------------------------------------
class TranslationBundle:
    """
    Translation bundle providing gettext translation functionality.
    """

    def __init__(self, domain: str, locale_: str):
        """Initialize the bundle."""
        self._domain = domain
        self._locale = locale_
        self._translator = gettext.NullTranslations()

        self.update_locale(locale_)

    def update_locale(self, locale_: str) -> None:
        """
        Update the locale.

        Parameters
        ----------
        locale_: str
            The language name to use.
        """
        # TODO: Need to handle packages that provide their own .mo files
        self._locale = locale_
        localedir = None
        if locale_ != DEFAULT_LOCALE:
            language_pack_module = f"jupyterlab_language_pack_{locale_}"
            try:
                mod = importlib.import_module(language_pack_module)
                assert mod.__file__ is not None
                localedir = os.path.join(os.path.dirname(mod.__file__), LOCALE_DIR)
            except Exception:  # noqa: S110
                # no-op
                pass

        self._translator = gettext.translation(
            self._domain, localedir=localedir, languages=(self._locale,), fallback=True
        )

    def gettext(self, msgid: str) -> str:
        """
        Translate a singular string.

        Parameters
        ----------
        msgid: str
            The singular string to translate.

        Returns
        -------
        str
            The translated string.
        """
        return self._translator.gettext(msgid)

    def ngettext(self, msgid: str, msgid_plural: str, n: int) -> str:
        """
        Translate a singular string with pluralization.

        Parameters
        ----------
        msgid: str
            The singular string to translate.
        msgid_plural: str
            The plural string to translate.
        n: int
            The number for pluralization.

        Returns
        -------
        str
            The translated string.
        """
        return self._translator.ngettext(msgid, msgid_plural, n)

    def pgettext(self, msgctxt: str, msgid: str) -> str:
        """
        Translate a singular string with context.

        Parameters
        ----------
        msgctxt: str
            The message context.
        msgid: str
            The singular string to translate.

        Returns
        -------
        str
            The translated string.
        """
        # Python 3.7 or lower does not offer translations based on context.
        # On these versions `pgettext` falls back to `gettext`
        if PY37_OR_LOWER:
            translation = self._translator.gettext(msgid)
        else:
            translation = self._translator.pgettext(msgctxt, msgid)

        return translation

    def npgettext(self, msgctxt: str, msgid: str, msgid_plural: str, n: int) -> str:
        """
        Translate a singular string with context and pluralization.

        Parameters
        ----------
        msgctxt: str
            The message context.
        msgid: str
            The singular string to translate.
        msgid_plural: str
            The plural string to translate.
        n: int
            The number for pluralization.

        Returns
        -------
        str
            The translated string.
        """
        # Python 3.7 or lower does not offer translations based on context.
        # On these versions `npgettext` falls back to `ngettext`
        if PY37_OR_LOWER:
            translation = self._translator.ngettext(msgid, msgid_plural, n)
        else:
            translation = self._translator.npgettext(msgctxt, msgid, msgid_plural, n)

        return translation

    # Shorthands
    def __(self, msgid: str) -> str:
        """
        Shorthand for gettext.

        Parameters
        ----------
        msgid: str
            The singular string to translate.

        Returns
        -------
        str
            The translated string.
        """
        return self.gettext(msgid)

    def _n(self, msgid: str, msgid_plural: str, n: int) -> str:
        """
        Shorthand for ngettext.

        Parameters
        ----------
        msgid: str
            The singular string to translate.
        msgid_plural: str
            The plural string to translate.
        n: int
            The number for pluralization.

        Returns
        -------
        str
            The translated string.
        """
        return self.ngettext(msgid, msgid_plural, n)

    def _p(self, msgctxt: str, msgid: str) -> str:
        """
        Shorthand for pgettext.

        Parameters
        ----------
        msgctxt: str
            The message context.
        msgid: str
            The singular string to translate.

        Returns
        -------
        str
            The translated string.
        """
        return self.pgettext(msgctxt, msgid)

    def _np(self, msgctxt: str, msgid: str, msgid_plural: str, n: int) -> str:
        """
        Shorthand for npgettext.

        Parameters
        ----------
        msgctxt: str
            The message context.
        msgid: str
            The singular string to translate.
        msgid_plural: str
            The plural string to translate.
        n: int
            The number for pluralization.

        Returns
        -------
        str
            The translated string.
        """
        return self.npgettext(msgctxt, msgid, msgid_plural, n)


class translator:
    """
    Translations manager.
    """

    _TRANSLATORS: dict[str, TranslationBundle] = {}
    _LOCALE = SYS_LOCALE

    @staticmethod
    def normalize_domain(domain: str) -> str:
        """Normalize a domain name.

        Parameters
        ----------
        domain: str
            Domain to normalize

        Returns
        -------
        str
            Normalized domain
        """
        return domain.replace("-", "_")

    @classmethod
    def set_locale(cls, locale_: str) -> None:
        """
        Set locale for the translation bundles based on the settings.

        Parameters
        ----------
        locale_: str
            The language name to use.
        """
        if locale_ == cls._LOCALE:
            # Nothing to do bail early
            return

        if is_valid_locale(locale_):
            cls._LOCALE = locale_
            for _, bundle in cls._TRANSLATORS.items():
                bundle.update_locale(locale_)

    @classmethod
    def load(cls, domain: str) -> TranslationBundle:
        """
        Load translation domain.

        The domain is usually the normalized ``package_name``.

        Parameters
        ----------
        domain: str
            The translations domain. The normalized python package name.

        Returns
        -------
        Translator
            A translator instance bound to the domain.
        """
        norm_domain = translator.normalize_domain(domain)
        if norm_domain in cls._TRANSLATORS:
            trans = cls._TRANSLATORS[norm_domain]
        else:
            trans = TranslationBundle(norm_domain, cls._LOCALE)
            cls._TRANSLATORS[norm_domain] = trans

        return trans

    @staticmethod
    def _translate_schema_strings(
        translations: Any,
        schema: dict,
        prefix: str = "",
        to_translate: dict[Pattern, str] | None = None,
    ) -> None:
        """Translate a schema in-place."""
        if to_translate is None:
            to_translate = _prepare_schema_patterns(schema)

        for key, value in schema.items():
            path = prefix + "/" + key

            if isinstance(value, str):
                matched = False
                for pattern, context in to_translate.items():  # noqa: B007
                    if pattern.fullmatch(path):
                        matched = True
                        break
                if matched:
                    schema[key] = translations.pgettext(context, value)
            elif isinstance(value, dict):
                translator._translate_schema_strings(
                    translations,
                    value,
                    prefix=path,
                    to_translate=to_translate,
                )
            elif isinstance(value, list):
                for i, element in enumerate(value):
                    if not isinstance(element, dict):
                        continue
                    translator._translate_schema_strings(
                        translations,
                        element,
                        prefix=path + "[" + str(i) + "]",
                        to_translate=to_translate,
                    )

    @staticmethod
    def translate_schema(schema: dict) -> dict:
        """Translate a schema.

        Parameters
        ----------
        schema: dict
            The schema to be translated

        Returns
        -------
        Dict
            The translated schema
        """
        if translator._LOCALE == DEFAULT_LOCALE:
            return schema

        translations = translator.load(
            schema.get(_lab_i18n_config, {}).get("domain", DEFAULT_DOMAIN)
        )

        new_schema = schema.copy()
        translator._translate_schema_strings(translations, new_schema)

        return new_schema
