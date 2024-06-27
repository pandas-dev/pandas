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
    import os

    # If not installed with setuptools, this might not be available
    try:
        import pkg_resources
    except ImportError:
        pkg_resources = None  # type: ignore[assignment]

    local_dir = os.path.basename(__file__)

    def get_translation_for(package_name: str) -> gettext.NullTranslations:  # type: ignore[misc]
        """Find and return gettext translation for package
        (Try to find folder manually if setuptools does not exist)
        """

        if "." in package_name:
            package_name = ".".join(package_name.split(".")[:-1])
        localedir = None

        if pkg_resources is None:
            mydir = os.path.join(local_dir, "i18n")
        else:
            mydir = pkg_resources.resource_filename(package_name, "i18n")

        for localedir in mydir, None:
            localefile = gettext.find(package_name, localedir)
            if localefile:
                break

        return gettext.translation(package_name, localedir=localedir, fallback=True)
