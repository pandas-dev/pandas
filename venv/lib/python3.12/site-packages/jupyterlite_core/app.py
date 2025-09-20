"""the JupyterLite CLI App(s)"""

import typing
from pathlib import Path

from jupyter_core.application import JupyterApp, base_aliases, base_flags
from jupyter_core.paths import jupyter_config_path
from traitlets import Bool, Instance, List, Unicode, default
from traitlets.utils.text import indent

from . import __version__
from .addons import get_addon_implementations, merge_addon_aliases, merge_addon_flags
from .config import LiteBuildConfig
from .constants import PHASES
from .manager import LiteManager

#: the total set of flags from discovered addons
lite_flags = {
    **{
        flag: value
        for flag, value in base_flags.items()
        if flag not in ["show-config", "show-config-json", "generate-config"]
    },
    "ignore-sys-prefix": (
        {"LiteBuildConfig": {"ignore_sys_prefix": True}},
        "Do not copy anything from sys.prefix",
    ),
    "no-sourcemaps": (
        {"LiteBuildConfig": {"no_sourcemaps": True}},
        "Strip all sourcemaps from applications and extensions",
    ),
    "no-unused-shared-packages": (
        {"LiteBuildConfig": {"no_unused_shared_packages": True}},
        "Remove shared packages not used by --apps",
    ),
    "no-libarchive": (
        {"LiteBuildConfig": {"no_libarchive": True}},
        "Do not try to use libarchive-c for archive operations",
    ),
}

lite_aliases = dict(
    **base_aliases,
    **{
        # meta options
        "disable-addons": "LiteBuildConfig.disable_addons",
        # input options
        "app-archive": "LiteBuildConfig.app_archive",
        "apps": "LiteBuildConfig.apps",
        # top-level
        "lite-dir": "LiteBuildConfig.lite_dir",
        # contents
        "contents": "LiteBuildConfig.contents",
        "ignore-contents": "LiteBuildConfig.ignore_contents",
        "extra-ignore-contents": "LiteBuildConfig.extra_ignore_contents",
        # settings
        "settings-overrides": "LiteBuildConfig.settings_overrides",
        # output options
        "output-dir": "LiteBuildConfig.output_dir",
        "output-archive": "LiteBuildConfig.output_archive",
        "source-date-epoch": "LiteBuildConfig.source_date_epoch",
        # server-specific things
        "port": "LiteBuildConfig.port",
        "base-url": "LiteBuildConfig.base_url",
    },
)


class DescribedMixin:
    """a self-describing mixin"""

    @property
    def description(self):  # pragma: no cover
        return self.__doc__.splitlines()[0].strip()


class BaseLiteApp(JupyterApp, LiteBuildConfig, DescribedMixin):
    """TODO: An undescribed app"""

    version = __version__

    config_file_name = Unicode("jupyter_lite_config").tag(config=True)

    config_file_paths = List(Unicode(help="Paths to search for jupyter_lite.(py|json)")).tag(
        config=True
    )

    @property
    def aliases(self):
        """Get CLI aliases, including ones provided by addons."""
        return merge_addon_aliases(lite_aliases)

    @property
    def flags(self):
        """Get CLI flags, including ones provided by addons."""
        return merge_addon_flags(lite_flags)

    @default("config_file_paths")
    def _config_file_paths_default(self):
        return [str(Path.cwd()), *jupyter_config_path()]

    def emit_alias_help(self):  # pragma: no cover
        """Yield the lines for alias part of the help.

        copied from:

        https://github.com/ipython/traitlets/blob/v5.8.0/traitlets/config/application.py#L517

        Unlike the upstream, this also takes Addon classes (and their parents)
        into consideration.
        """
        if not self.aliases:
            return

        classdict = {}
        # discover more classes
        for cls in [*self.classes, *get_addon_implementations().values()]:
            # include all parents (up to, but excluding Configurable) in available names
            for c in cls.mro()[:-3]:
                classdict[c.__name__] = c

        fhelp: typing.Optional[str]
        for alias, longname in self.aliases.items():
            try:
                if isinstance(longname, tuple):
                    name, fhelp = longname
                else:
                    name, fhelp = longname, None
                classname, traitname = name.split(".")[-2:]
                name = classname + "." + traitname

                cls = classdict[classname]

                trait = cls.class_traits(config=True)[traitname]
                fhelp = cls.class_get_trait_help(trait, helptext=fhelp).splitlines()

                aliases = (alias,) if not isinstance(alias, tuple) else alias
                aliases = sorted(aliases, key=len)  # type:ignore[assignment]
                aliases = ", ".join(("--%s" if len(m) > 1 else "-%s") % m for m in aliases)

                # reformat first line
                assert fhelp is not None
                fhelp[0] = fhelp[0].replace(f"--{name}", f"--{alias}")  # type:ignore
                yield from fhelp
                yield indent(f"Equivalent to: [--{name}]")
            except Exception as ex:
                self.log.error("Failed collecting help-message for alias %r, due to: %s", alias, ex)
                raise


class ManagedApp(BaseLiteApp):
    """An app with a LiteManager that can do some config fixing"""

    lite_manager = Instance(LiteManager)

    @default("lite_manager")
    def _default_manager(self):  # noqa: C901, PLR0912
        kwargs = dict(
            parent=self,
        )
        if self.lite_dir:
            kwargs["lite_dir"] = self.lite_dir
        if self.app_archive:
            kwargs["app_archive"] = self.app_archive
        if self.no_libarchive:
            kwargs["no_libarchive"] = self.no_libarchive
        if self.output_dir:
            kwargs["output_dir"] = self.output_dir
        if self.file_types:
            kwargs["file_types"] = self.file_types
        if self.extra_file_types:
            kwargs["extra_file_types"] = self.extra_file_types
        if self.contents:
            kwargs["contents"] = [Path(p) for p in self.contents]
        if self.ignore_contents:
            kwargs["ignore_contents"] = self.ignore_contents
        if self.extra_ignore_contents:
            kwargs["extra_ignore_contents"] = self.extra_ignore_contents
        if self.settings_overrides:
            kwargs["settings_overrides"] = [Path(p) for p in self.settings_overrides]
        if self.apps:
            kwargs["apps"] = self.apps
        if self.no_sourcemaps is not None:
            kwargs["no_sourcemaps"] = self.no_sourcemaps
        if self.no_unused_shared_packages is not None:
            kwargs["no_unused_shared_packages"] = self.no_unused_shared_packages
        if self.output_archive:
            kwargs["output_archive"] = Path(self.output_archive)
        if self.disable_addons:
            kwargs["disable_addons"] = self.disable_addons
        if self.source_date_epoch is not None:
            kwargs["source_date_epoch"] = self.source_date_epoch
        if self.port is not None:
            kwargs["port"] = self.port
        if self.base_url is not None:
            kwargs["base_url"] = self.base_url
        if self.federated_extensions is not None:
            kwargs["federated_extensions"] = self.federated_extensions
        if self.ignore_sys_prefix is not None:
            kwargs["ignore_sys_prefix"] = self.ignore_sys_prefix

        return LiteManager(**kwargs)

    def start(self):
        self.lite_manager.initialize()


# doit stuff
class LiteDoitApp(ManagedApp):
    """Run the doit command"""

    _doit_cmd = None

    def start(self):
        super().start()
        self.exit(self.lite_manager.doit_run(*self._doit_cmd))


# special non-task apps
class LiteRawDoitApp(LiteDoitApp):
    """use the full doit CLI, see https://pydoit.org/contents.html

    tell jupyter to not parse the arguments with --, e.g.

        jupyter-lite doit -- --help
    """

    def parse_command_line(self, argv=None):
        super().parse_command_line(argv)
        self._doit_cmd = [*(self.extra_args or [])]


class LiteListApp(LiteDoitApp):
    """describe a JupyterLite site"""

    _doit_cmd = ["list", "--all", "--status"]


# task app base class
class LiteTaskApp(LiteDoitApp):
    """run a doit task, optionally with --force"""

    force = Bool(False, help="forget previous runs of task and re-run from the beginning").tag(
        config=True
    )

    @property
    def flags(self):
        """CLI flags, including some custom ones."""
        return merge_addon_flags(
            {
                **lite_flags,
                "force": ({"LiteTaskApp": {"force": True}}, LiteTaskApp.force.help),
            }
        )

    _doit_task = None

    @property
    def _doit_cmd(self):
        return [f"{phase}{self._doit_task}" for phase in PHASES]

    def start(self):
        super().start()
        if self.force:
            for phase in PHASES:
                self.lite_manager.doit_run("forget", f"{phase}{self._doit_task}")


# the tasks
class LiteStatusApp(LiteTaskApp):
    """report about what a JupyterLite build _might_ do"""

    _doit_task = "status"


class LiteInitApp(LiteTaskApp):
    """initialize a JupyterLite site from an app archive baseline"""

    _doit_task = "init"


class LiteBuildApp(LiteTaskApp):
    """build a JupyterLite site, including user content"""

    _doit_task = "build"


class LiteCheckApp(LiteTaskApp):
    """verify a JupyterLite site, using available schema and rules"""

    _doit_task = "check"


class LiteServeApp(LiteTaskApp):
    """serve a JupyterLite site, using best available HTTP server"""

    _doit_task = "serve"


class LiteArchiveApp(LiteTaskApp):
    """build a JupyterLite app archive, which can be used as a baseline"""

    _doit_task = "archive"


class LiteApp(BaseLiteApp):
    """build ready-to-serve (or -publish) JupyterLite sites"""

    subcommands = {
        k: (v, v.__doc__.splitlines()[0].strip())
        for k, v in dict(
            # special apps
            list=LiteListApp,
            # task apps
            status=LiteStatusApp,
            init=LiteInitApp,
            build=LiteBuildApp,
            check=LiteCheckApp,
            serve=LiteServeApp,
            archive=LiteArchiveApp,
            # more special apps
            doit=LiteRawDoitApp,
        ).items()
    }


main = launch_new_instance = LiteApp.launch_instance

if __name__ == "__main__":  # pragma: nocover
    main()
