from __future__ import annotations

import locale

# High performance method for English (no translation needed)
loc = locale.getlocale()[0]
if loc is None or loc.startswith("en") or loc == "C":

    class NullTranslation:
        def gettext(self, str1: str) -> str:  # pylint: disable=no-self-use
            return str1

        def ngettext(self, str1, strN, n):  # pylint: disable=no-self-use
            if n == 1:
                return str1.replace("{0}", str(n))

            return strN.replace("{0}", str(n))

    def get_translation_for(
        package_name: str,  # noqa: ARG001
    ) -> NullTranslation:
        return NullTranslation()

else:
    import gettext
    import sys

    if sys.version_info < (3, 9):
        from importlib_resources import as_file, files
    else:
        from importlib.resources import as_file, files

    def get_translation_for(package_name: str) -> gettext.NullTranslations:  # type: ignore[misc]
        """
        Find and return gettext translation for package
        """
        assert loc is not None

        if "." in package_name:
            package_name = ".".join(package_name.split(".")[:-1])

        localedir = None

        with as_file(files(package_name) / "i18n") as mydir:
            for localedir in mydir, None:
                localefile = gettext.find(package_name, localedir, languages=[loc])
                if localefile:
                    break

            return gettext.translation(
                package_name, localedir=localedir, languages=[loc], fallback=True
            )
