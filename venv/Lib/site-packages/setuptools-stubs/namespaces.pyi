import itertools

flatten = itertools.chain.from_iterable

class Installer:
    nspkg_ext: str
    def install_namespaces(self) -> None: ...
    def uninstall_namespaces(self) -> None: ...

class DevelopInstaller(Installer): ...
