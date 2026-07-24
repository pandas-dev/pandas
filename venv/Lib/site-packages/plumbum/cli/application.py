from __future__ import annotations

__lazy_modules__ = {
    "collections",
    "contextlib",
    f"{__spec__.parent}.terminal",
    "functools",
    "inspect",
    "plumbum.colorlib",
    "plumbum.colorlib.styles",
    "plumbum.lib",
    "textwrap",
}

import contextlib
import errno
import functools
import inspect
import os
import sys
import typing
from collections import defaultdict
from textwrap import TextWrapper
from typing import Any, ClassVar, Literal, NoReturn, TypeVar

from plumbum import colors, local
from plumbum.cli.i18n import get_translation_for
from plumbum.colorlib.styles import Style
from plumbum.lib import getdoc

from .switches import (
    CountOf,
    Flag,
    MissingArgument,
    MissingMandatorySwitch,
    PositionalArgumentsError,
    Set,
    SubcommandError,
    SwitchCombinationError,
    SwitchError,
    UnknownSwitch,
    WrongArgumentType,
    switch,
)
from .terminal import get_terminal_size

if typing.TYPE_CHECKING:
    from collections.abc import Callable, Generator, MutableMapping, Sequence

    from plumbum._compat.typing import Self
    from plumbum.cli.switches import SwitchInfo

_translation = get_translation_for(__name__)
T_, ngettext = _translation.gettext, _translation.ngettext


T = TypeVar("T")


class ShowHelp(SwitchError):
    pass


class ShowHelpAll(SwitchError):
    pass


class ShowVersion(SwitchError):
    pass


class ShowCompletion(SwitchError):
    pass


class SwitchParseInfo:
    __slots__ = ("__weakref__", "index", "swname", "val")

    def __init__(self, swname: str, val: tuple[Any, ...], index: int):
        self.swname = swname
        self.val = val
        self.index = index


class Subcommand:
    __slots__ = ("name", "subapplication")

    def __init__(self, name: str, subapplication: type[Application] | str):
        self.name = name
        self.subapplication = subapplication

    def get(self) -> type[Application]:
        if isinstance(self.subapplication, str):
            modname, clsname = self.subapplication.rsplit(".", 1)
            mod = __import__(modname, None, None, "*")
            try:
                cls = getattr(mod, clsname)
            except AttributeError:
                raise ImportError(f"cannot import name {clsname}") from None
            self.subapplication = cls
        assert not isinstance(self.subapplication, str)
        return self.subapplication

    def __repr__(self) -> str:
        return T_("Subcommand({self.name}, {self.subapplication})").format(self=self)


_switch_groups = ["Switches", "Meta-switches"]
_switch_groups_l10n = [T_("Switches"), T_("Meta-switches")]


# ===================================================================================================
# CLI Application base class
# ===================================================================================================

X = TypeVar("X", "type[Application]", str)


class Application:
    """The base class for CLI applications; your "entry point" class should derive from it,
    define the relevant switch functions and attributes, and the ``main()`` function.
    The class defines two overridable "meta switches" for version (``-v``, ``--version``)
    help (``-h``, ``--help``), and help-all (``--help-all``).

    The signature of the main function matters: any positional arguments (e.g., non-switch
    arguments) given on the command line are passed to the ``main()`` function; if you wish
    to allow unlimited number of positional arguments, use varargs (``*args``). The names
    of the arguments will be shown in the help message.

    The classmethod ``run`` serves as the entry point of the class. It parses the command-line
    arguments, invokes switch functions and enters ``main``. You should **not override** this
    method.

    Usage::

        class FileCopier(Application):
            stat = Flag("p", "copy stat info as well")

            def main(self, src, dst):
                if self.stat:
                    shutil.copy2(src, dst)
                else:
                    shutil.copy(src, dst)

        if __name__ == "__main__":
            FileCopier.run()

    There are several class-level attributes you may set:

    * ``PROGNAME`` - the name of the program; if ``None`` (the default), it is set to the
      name of the executable (``argv[0]``); can be in color. If only a color, will be applied to the name.

    * ``VERSION`` - the program's version (defaults to ``1.0``, can be in color)

    * ``DESCRIPTION`` - a short description of your program (shown in help). If not set,
      the class' ``__doc__`` will be used. Can be in color.

    * ``DESCRIPTION_MORE`` - a detailed description of your program (shown in help). The text will be printed
      by paragraphs (specified by empty lines between them). The indentation of each paragraph will be the
      indentation of its first line. List items are identified by their first non-whitespace character being
      one of '-', '*', and '/'; so that they are not combined with preceding paragraphs. Bullet '/' is
      "invisible", meaning that the bullet itself will not be printed to the output.

    * ``USAGE`` - the usage line (shown in help).

    * ``COLOR_USAGE_TITLE`` - The color of the usage line's header.

    * ``COLOR_USAGE`` - The color of the usage line.

    * ``COLOR_GROUPS`` - A dictionary that sets colors for the groups, like Meta-switches, Switches,
      and Subcommands.

    * ``COLOR_GROUP_TITLES`` - A dictionary that sets colors for the group titles. If the dictionary is empty,
      it defaults to ``COLOR_GROUPS``.

    * ``SUBCOMMAND_HELPMSG`` - Controls the printing of extra "see subcommand -h" help message.
      Default is a message, set to False to remove.

    * ``ALLOW_ABBREV`` - Controls whether partial switch names are supported, for example '--ver' will match
      '--verbose'. Default is False for backward consistency with previous plumbum releases. Note that ambiguous
      abbreviations will not match, for example if --foothis and --foothat are defined, then --foo will not match.

    A note on sub-commands: when an application is the root, its ``parent`` attribute is set to
    ``None``. When it is used as a nested-command, ``parent`` will point to its direct ancestor.
    Likewise, when an application is invoked with a sub-command, its ``nested_command`` attribute
    will hold the chosen sub-application and its command-line arguments (a tuple); otherwise, it
    will be set to ``None``

    """

    # Some that are not typed None will always be set in __init__
    PROGNAME: str | Style = None  # type: ignore[assignment]
    DESCRIPTION: str | None = None
    DESCRIPTION_MORE: str | None = None
    VERSION: str | None = None
    USAGE: str | None = None
    COLOR_USAGE: Style = None  # type: ignore[assignment]
    COLOR_USAGE_TITLE: Style | None = None
    COLOR_GROUPS: MutableMapping[str, Style] = None  # type: ignore[assignment]
    COLOR_GROUP_TITLES: MutableMapping[str, Style] = None  # type: ignore[assignment]
    CALL_MAIN_IF_NESTED_COMMAND: bool = True
    SUBCOMMAND_HELPMSG: str | Literal[False] = T_(
        "see '{parent} {sub} --help' for more info"
    )
    ALLOW_ABBREV: bool = False

    parent: Self | None = None
    nested_command: tuple[type[Application], list[str]] | None = None
    _unbound_switches: ClassVar[tuple[str, ...]] = ()
    # Set transiently by ``helpall`` so that ``help`` demotes meta-switches for
    # the nested render only, without mutating the shared SwitchInfo objects.
    _hide_meta_switches: bool = False

    def __new__(cls, executable: object | None = None) -> Self:
        """Allows running the class directly as a shortcut for main.
        This is necessary for some setup scripts that want a single function,
        instead of an expression with a dot in it."""

        if executable is None:
            # This return value was not a class instance, so __init__ is never called
            cls.run()

        return super().__new__(cls)

    def __init__(self, executable: str | None = None):
        assert executable is not None

        # Filter colors

        if self.PROGNAME is None:
            self.PROGNAME = os.path.basename(executable)
        elif isinstance(self.PROGNAME, Style):
            self.PROGNAME = self.PROGNAME | os.path.basename(executable)
        elif colors.filter(self.PROGNAME) == "":
            self.PROGNAME = colors.extract(self.PROGNAME) | os.path.basename(executable)
        if self.DESCRIPTION is None:
            self.DESCRIPTION = getdoc(self)

        # Allow None for the colors
        self.COLOR_GROUPS = defaultdict(
            lambda: colors.do_nothing,
            {} if self.COLOR_GROUPS is None else self.COLOR_GROUPS,  # type: ignore[redundant-expr]
        )

        self.COLOR_GROUP_TITLES = defaultdict(
            lambda: colors.do_nothing,
            self.COLOR_GROUPS
            if self.COLOR_GROUP_TITLES is None  # type: ignore[redundant-expr]
            else self.COLOR_GROUP_TITLES,
        )
        if type(self).COLOR_USAGE is None:
            self.COLOR_USAGE = colors.do_nothing

        self.executable = executable
        self._switches_by_name: dict[str, SwitchInfo] = {}
        self._switches_by_func = {}
        self._switches_by_envar = {}
        self._subcommands = {}

        for cls in reversed(type(self).mro()):
            for obj in cls.__dict__.values():
                if isinstance(obj, Subcommand):
                    name = colors.filter(obj.name)
                    if name.startswith("-"):
                        raise SubcommandError(
                            T_("Sub-command names cannot start with '-'")
                        )
                    # it's okay for child classes to override sub-commands set by their parents
                    self._subcommands[name] = obj
                    continue

                swinfo = getattr(obj, "_switch_info", None)
                if not swinfo:
                    continue
                for name in swinfo.names:
                    if name in self._unbound_switches:
                        continue
                    if (
                        name in self._switches_by_name
                        and not self._switches_by_name[name].overridable
                    ):
                        raise SwitchError(
                            T_(
                                "Switch {name} already defined and is not overridable"
                            ).format(name=name)
                        )
                    self._switches_by_name[name] = swinfo
                    self._switches_by_func[swinfo.func] = swinfo
                    if swinfo.envname:
                        self._switches_by_envar[swinfo.envname] = swinfo

    @property
    def root_app(self) -> Application:
        return self.parent.root_app if self.parent else self

    @classmethod
    def unbind_switches(cls, *switch_names: str) -> None:
        """Unbinds the given switch names from this application. For example

        ::

            class MyApp(cli.Application):
                pass
            MyApp.unbind_switches("--version")

        """
        cls._unbound_switches += tuple(
            name.lstrip("-") for name in switch_names if name
        )

    @typing.overload
    @classmethod
    def subcommand(
        cls, name: str, subapp: None = ...
    ) -> Callable[[type[Application]], type[Application]]: ...

    @typing.overload
    @classmethod
    def subcommand(cls, name: str, subapp: type[Application]) -> type[Application]: ...

    @typing.overload
    @classmethod
    def subcommand(cls, name: str, subapp: str) -> str: ...

    @classmethod
    def subcommand(
        cls, name: str, subapp: type[Application] | str | None = None
    ) -> Callable[[type[Application]], type[Application]] | type[Application] | str:
        """Registers the given sub-application as a sub-command of this one. This method can be
        used both as a decorator and as a normal ``classmethod``::

            @MyApp.subcommand("foo")
            class FooApp(cli.Application):
                pass

        Or ::

            MyApp.subcommand("foo", FooApp)

        .. versionadded:: 1.1

        .. versionadded:: 1.3
            The sub-command can also be a string, in which case it is treated as a
            fully-qualified class name and is imported on demand. For example,

            MyApp.subcommand("foo", "fully.qualified.package.FooApp")

        """

        def wrapper(subapp: X) -> X:
            # Use the subcommand name (not subapp name) to ensure uniqueness
            # This allows the same subapp to be registered under multiple names
            attrname = f"_subcommand_{name}"
            setattr(cls, attrname, Subcommand(name, subapp))
            return subapp

        if subapp is None:
            return wrapper
        if isinstance(subapp, str):
            return wrapper(subapp)
        return wrapper(subapp)

    def _get_partial_matches(self, partialname: str) -> list[str]:
        return [
            switch_
            for switch_ in self._switches_by_name
            if switch_.startswith(partialname)
        ]

    def _lookup_switch(self, name: str, swinfo: SwitchInfo) -> SwitchInfo:
        try:
            return self._switches_by_name[name]
        except KeyError:
            raise SwitchError(
                T_(
                    "Switch {0} refers to an unknown switch {1} "
                    "in its requires/excludes"
                ).format(
                    ("-" if len(swinfo.names[0]) == 1 else "--") + swinfo.names[0],
                    name,
                )
            ) from None

    def _parse_args(
        self, argv: list[str]
    ) -> tuple[dict[Callable[..., None], Any], list[str]]:
        tailargs = []
        swfuncs: dict[Callable[..., None], Any] = {}
        index = 0

        while argv:
            index += 1
            a = argv.pop(0)
            val = None
            if a == "--":
                # end of options, treat the rest as tailargs
                tailargs.extend(argv)
                break

            if a in self._subcommands:
                subcmd = self._subcommands[a].get()
                self.nested_command = (
                    subcmd,
                    [self.PROGNAME + " " + self._subcommands[a].name, *argv],
                )
                break

            if a.startswith("--") and len(a) >= 3:
                # [--name], [--name=XXX], [--name, XXX], [--name, ==, XXX],
                # [--name=, XXX], [--name, =XXX]
                eqsign = a.find("=")
                has_eq = eqsign >= 0
                if has_eq:
                    name = a[2:eqsign]
                    argv.insert(0, a[eqsign:])
                else:
                    name = a[2:]

                # An exact match always wins over abbreviation (argparse-style).
                if self.ALLOW_ABBREV and name not in self._switches_by_name:
                    partials = self._get_partial_matches(name)
                    if len(partials) == 1:
                        name = partials[0]
                    elif len(partials) > 1:
                        raise UnknownSwitch(
                            T_("Ambiguous partial switch {0}").format("--" + name)
                        )

                swname = "--" + name
                if name not in self._switches_by_name:
                    raise UnknownSwitch(T_("Unknown switch {0}").format(swname))
                swinfo = self._switches_by_name[name]
                if not swinfo.argtype and has_eq:
                    raise SwitchError(
                        T_("Switch {0} does not take an argument").format(swname)
                    )
                if swinfo.argtype:
                    if not argv:
                        raise MissingArgument(
                            T_("Switch {0} requires an argument").format(swname)
                        )
                    a = argv.pop(0)
                    if a and a[0] == "=":
                        if len(a) >= 2:
                            val = a[1:]
                        else:
                            if not argv:
                                raise MissingArgument(
                                    T_("Switch {0} requires an argument").format(swname)
                                )
                            val = argv.pop(0)
                    else:
                        val = a

            elif a.startswith("-") and len(a) >= 2:
                # [-a], [-a, XXX], [-aXXX], [-abc]
                name = a[1]
                swname = "-" + name
                if name not in self._switches_by_name:
                    raise UnknownSwitch(T_("Unknown switch {0}").format(swname))
                swinfo = self._switches_by_name[name]
                if swinfo.argtype:
                    if len(a) >= 3:
                        val = a[2:]
                    else:
                        if not argv:
                            raise MissingArgument(
                                T_("Switch {0} requires an argument").format(swname)
                            )
                        val = argv.pop(0)
                elif len(a) >= 3:
                    argv.insert(0, "-" + a[2:])

            else:
                if a.startswith("-"):
                    raise UnknownSwitch(T_("Unknown switch {0}").format(a))
                tailargs.append(a)
                continue

            # handle argument
            if typing.TYPE_CHECKING:
                assert val is not None
            val = self._handle_argument(val, swinfo.argtype, name)

            if swinfo.func in swfuncs:
                if swinfo.list:
                    swfuncs[swinfo.func].val[0].append(val)
                else:
                    if swfuncs[swinfo.func].swname == swname:
                        raise SwitchError(T_("Switch {0} already given").format(swname))
                    raise SwitchError(
                        T_("Switch {0} already given ({1} is equivalent)").format(
                            swfuncs[swinfo.func].swname, swname
                        )
                    )
            elif swinfo.list:
                swfuncs[swinfo.func] = SwitchParseInfo(swname, ([val],), index)
            elif val is NotImplemented:
                swfuncs[swinfo.func] = SwitchParseInfo(swname, (), index)
            else:
                swfuncs[swinfo.func] = SwitchParseInfo(swname, (val,), index)

        # Extracting arguments from environment variables
        envindex = 0
        for env, swinfo in self._switches_by_envar.items():
            envindex -= 1
            envval = local.env.get(env)
            if envval is None:
                continue

            if swinfo.func in swfuncs:
                continue  # skip if overridden by command line arguments

            val = self._handle_argument(envval, swinfo.argtype, env)
            envname = f"${env}"
            if swinfo.list:
                # multiple values over environment variables are not supported,
                # this will require some sort of escaping and separator convention
                swfuncs[swinfo.func] = SwitchParseInfo(envname, ([val],), envindex)
            elif val is NotImplemented:
                swfuncs[swinfo.func] = SwitchParseInfo(envname, (), envindex)
            else:
                swfuncs[swinfo.func] = SwitchParseInfo(envname, (val,), envindex)

        return swfuncs, tailargs

    @classmethod
    def autocomplete(cls, argv: list[str]) -> None:
        """Handle shell completion requests.

        For bash, this uses COMP_WORDS and COMP_CWORD environment variables.
        For fish, this uses the completions switch.
        """
        comp_words = os.environ.get("COMP_WORDS", "")
        comp_cword = os.environ.get("COMP_CWORD", "")

        if comp_words and comp_cword:
            cls._handle_bash_completion(argv, comp_words, int(comp_cword))

    @classmethod
    def _handle_bash_completion(
        cls, argv: list[str], comp_words: str, comp_cword: int
    ) -> None:
        """Handle bash completion request."""
        words = comp_words.split()
        if comp_cword >= len(words):
            comp_cword = len(words) - 1
        partial = words[comp_cword] if comp_cword < len(words) else ""
        previous_words = words[1:comp_cword] if comp_cword > 0 else []

        inst = cls(argv[0] if argv else words[0])
        completions = inst.get_completions([*previous_words, partial])
        cls._print_completions_and_exit(completions)

    @classmethod
    def _print_completions_and_exit(cls, completions: list[str]) -> None:
        """Print completions to stdout and exit."""
        for completion in completions:
            print(completion)
        sys.exit(0)

    def get_completions(self, argv: list[str]) -> list[str]:
        """Get completion suggestions based on partial command line arguments.

        :param argv: List of command line arguments (not including the program name)
        :return: List of completion suggestions
        """
        if not argv:
            return self._get_all_completions()

        partial = argv[-1]
        previous = argv[:-1]

        if previous and previous[-1].startswith("-"):
            last_switch = previous[-1]
            swinfo = self._get_switch_info_for_arg(last_switch)
            if swinfo and swinfo.argtype:
                return self._get_switch_arg_completions(swinfo, partial)

        if partial.startswith("-"):
            swinfo = self._get_switch_info_for_arg(partial)
            if swinfo is not None and not swinfo.argtype:
                return self._get_all_completions()
            return self._get_switch_completions(partial)

        nested_result = self._try_nested_completion(argv)
        if nested_result is not None:
            return nested_result

        return self._get_all_completions(partial)

    def get_nested_completions(
        self, argv: list[str]
    ) -> tuple[Application | None, list[str]]:
        """Get the nested application and completions for subcommand-aware completion.

        :param argv: List of command line arguments
        :return: Tuple of (nested application or None, list of completions)
        """
        if not argv:
            return None, []

        for i, arg in enumerate(argv):
            if arg in self._subcommands:
                subapp_cls = self._subcommands[arg].get()
                subapp = subapp_cls(f"{self.PROGNAME} {arg}")
                subapp.parent = self
                remaining_args = argv[i + 1 :]
                return subapp, subapp.get_completions(remaining_args)

        return None, []

    def _get_switch_info_for_arg(self, switch_name: str) -> SwitchInfo | None:
        """Get switch info for a switch name (with or without dashes)."""
        name = switch_name.lstrip("-")
        if name in self._switches_by_name:
            return self._switches_by_name[name]
        if self.ALLOW_ABBREV:
            for sw_name, sw_info in self._switches_by_name.items():
                if sw_name.startswith(name):
                    return sw_info
        return None

    def _get_all_completions(self, partial: str = "") -> list[str]:
        """Get all possible completions (switches and subcommands)."""
        completions: list[str] = []

        switch_completions = self._get_switch_completions(partial)
        completions.extend(
            switch_completions
            + [
                subcmd_name
                for subcmd_name in self._subcommands
                if not partial or subcmd_name.startswith(partial)
            ]
        )

        return completions

    def _get_switch_completions(self, partial: str) -> list[str]:
        """Get completions for switches matching the partial string."""
        completions = []
        partial_name = partial.lstrip("-")

        for name in self._switches_by_name:
            if partial_name and not name.startswith(partial_name):
                continue
            prefix = "-" if len(name) == 1 else "--"
            completion = prefix + name
            if completion not in completions:
                completions.append(completion)

        return completions

    @staticmethod
    def _get_switch_arg_completions(swinfo: SwitchInfo, partial: str) -> list[str]:
        """Get completions for a switch argument based on its type/validator."""
        if swinfo.argtype is None:
            return []

        argtype = swinfo.argtype

        if hasattr(argtype, "choices"):
            try:
                choices = argtype.choices(partial)
                if isinstance(choices, set):
                    return sorted(choices)
                return list(choices)
            except (TypeError, AttributeError):
                pass

        if hasattr(argtype, "__name__"):
            type_name = argtype.__name__
            if type_name == "bool":
                return ["true", "false"]

        return []

    def _try_nested_completion(self, argv: list[str]) -> list[str] | None:
        """Try to get completions from nested subcommands."""
        nested_app, nested_completions = self.get_nested_completions(argv)
        if nested_app is not None:
            return nested_completions
        return None

    @classmethod
    def completion_script_bash(cls, prog_name: str) -> str:
        """Generate a bash completion script for this application.

        :param prog_name: The program name to use in the completion script
        :return: Bash completion script as a string
        """
        inst = cls(prog_name)

        switch_names = []
        for name in inst._switches_by_name:
            prefix = "-" if len(name) == 1 else "--"
            switch_names.append(prefix + name)

        subcommand_names = list(inst._subcommands.keys())

        script = f'''# Bash completion for {prog_name}
# Generated by plumbum.cli.Application

_{prog_name}_completion() {{
    local cur="${{COMP_WORDS[COMP_CWORD]}}"
    local prev="${{COMP_WORDS[COMP_CWORD-1]}}"
    local words="${{COMP_WORDS[@]}}"

    COMPREPLY=()

    # Static completions for switches and subcommands
    local switches="{" ".join(switch_names)}"
    local subcommands="{" ".join(subcommand_names)}"

    # If previous word is a switch that takes an argument, let the app handle it
    case "$prev" in'''

        for name, swinfo in inst._switches_by_name.items():
            if swinfo.argtype:
                prefix = "-" if len(name) == 1 else "--"
                script += f"""
        {prefix}{name})
            COMPREPLY=( $(compgen -W "$({prog_name} --complete $cur 2>/dev/null)" -- "$cur") )
            return
            ;;"""

        script += f"""
        *)
            ;;
    esac

    # Complete switches and subcommands
    if [[ "$cur" == -* ]]; then
        COMPREPLY=( $(compgen -W "$switches" -- "$cur") )
    else
        COMPREPLY=( $(compgen -W "$subcommands" -- "$cur") )
    fi
}}

complete -F _{prog_name}_completion {prog_name}
"""
        return script

    @classmethod
    def completion_script_fish(cls, prog_name: str) -> str:
        """Generate a fish completion script for this application.

        :param prog_name: The program name to use in the completion script
        :return: Fish completion script as a string
        """
        lines = [
            f"# Fish completion for {prog_name}",
            "# Generated by plumbum.cli.Application",
            "",
        ]
        inst = cls(prog_name)

        for name, swinfo in inst._switches_by_name.items():
            opt_flag = "-s" if len(name) == 1 else "-l"
            args = ""
            if swinfo.argtype:
                if hasattr(swinfo.argtype, "choices"):
                    try:
                        choices = swinfo.argtype.choices("")
                        choices_str = " ".join(
                            sorted(choices) if isinstance(choices, set) else choices
                        )
                        args = " -r -a '" + choices_str + "'"
                    except (TypeError, AttributeError):
                        args = " -r"
                else:
                    args = " -r"
            desc = ""
            if swinfo.help:
                escaped = swinfo.help.replace("'", "\\'")
                desc = " -d '" + escaped + "'"
            lines.append(f"complete -c {prog_name} {opt_flag} {name}{args}{desc}")

        for subcmd_name, subapp in inst._subcommands.items():
            subapp_cls = subapp.get()
            desc = ""
            if subapp_cls.DESCRIPTION:
                escaped = subapp_cls.DESCRIPTION.replace("'", "\\'")
                desc = " -d '" + escaped + "'"
            lines.append(
                f"complete -c {prog_name} -n '__fish_use_subcommand' -a '{subcmd_name}'{desc}"
            )

        return "\n".join(lines)

    @staticmethod
    def _handle_argument(val: str, argtype: Callable[[str], T] | None, name: str) -> T:
        if argtype is not None:
            try:
                return argtype(val)
            except (TypeError, ValueError) as ex:
                raise WrongArgumentType(
                    T_(
                        "Argument of {name} expected to be {argtype}, not {val!r}:\n    {ex!r}"
                    ).format(name=name, argtype=argtype, val=val, ex=ex)
                ) from None
        else:
            # TODO: This is required to handle (correctly) None, but probably could be done better
            return NotImplemented  # type: ignore[no-any-return]

    def _validate_args(
        self,
        swfuncs: dict[Callable[..., Any], SwitchParseInfo],
        tailargs: Sequence[str],
    ) -> tuple[list[tuple[Callable[..., Any], tuple[Any, ...]]], list[str]]:
        if self.help.__func__ in swfuncs:  # type: ignore[attr-defined]
            raise ShowHelp()
        if self.helpall.__func__ in swfuncs:  # type: ignore[attr-defined]
            raise ShowHelpAll()
        if self.version.__func__ in swfuncs:  # type: ignore[attr-defined]
            raise ShowVersion()
        if self.completions.__func__ in swfuncs:  # type: ignore[attr-defined]
            raise ShowCompletion()

        requirements = {}
        exclusions = {}
        for swinfo in self._switches_by_func.values():
            if swinfo.mandatory and swinfo.func not in swfuncs:
                raise MissingMandatorySwitch(
                    T_("Switch {0} is mandatory").format(
                        "/".join(
                            ("-" if len(n) == 1 else "--") + n for n in swinfo.names
                        )
                    )
                )
            requirements[swinfo.func] = {
                self._lookup_switch(req, swinfo) for req in swinfo.requires
            }
            exclusions[swinfo.func] = {
                self._lookup_switch(exc, swinfo) for exc in swinfo.excludes
            }

        # TODO: compute topological order

        gotten = set(swfuncs.keys())
        for func in gotten:
            missing = {f.func for f in requirements[func]} - gotten
            if missing:
                raise SwitchCombinationError(
                    T_("Given {0}, the following are missing {1}").format(
                        swfuncs[func].swname,
                        [self._switches_by_func[f].names[0] for f in missing],
                    )
                )
            invalid = {f.func for f in exclusions[func]} & gotten
            if invalid:
                raise SwitchCombinationError(
                    T_("Given {0}, the following are invalid {1}").format(
                        swfuncs[func].swname, [swfuncs[f].swname for f in invalid]
                    )
                )

        m = inspect.getfullargspec(self.main)

        max_args = sys.maxsize if m.varargs else len(m.args) - 1
        min_args = len(m.args) - 1 - (len(m.defaults) if m.defaults else 0)
        if len(tailargs) < min_args:
            raise PositionalArgumentsError(
                ngettext(
                    "Expected at least {0} positional argument, got {1}",
                    "Expected at least {0} positional arguments, got {1}",
                    min_args,
                ).format(min_args, tailargs)
            )
        if len(tailargs) > max_args:
            raise PositionalArgumentsError(
                ngettext(
                    "Expected at most {0} positional argument, got {1}",
                    "Expected at most {0} positional arguments, got {1}",
                    max_args,
                ).format(max_args, tailargs)
            )

        # Positional argument validation
        if hasattr(self.main, "positional"):
            new_tailargs = self._positional_validate(
                tailargs,
                self.main.positional,
                self.main.positional_varargs,
                m.args[1:],
                m.varargs,
            )

        elif hasattr(m, "annotations") and m.annotations:
            annotations = typing.get_type_hints(self.main)
            args_names = list(m.args[1:])
            positional: list[Any] = [None] * len(args_names)
            varargs = None

            # All args are positional, so convert kargs to positional
            for item, annotation in annotations.items():
                if item == m.varargs:
                    varargs = annotation
                elif item != "return":
                    positional[args_names.index(item)] = annotation

            new_tailargs = self._positional_validate(
                tailargs, positional, varargs, m.args[1:], m.varargs
            )

        else:
            new_tailargs = list(tailargs)

        ordered = [
            (f, a)
            for _, f, a in sorted((sf.index, f, sf.val) for f, sf in swfuncs.items())
        ]
        return ordered, new_tailargs

    def _positional_validate(
        self,
        args: Sequence[str],
        validator_list: list[Callable[[str], Any]],
        varargs: Callable[[str], Any] | None,
        argnames: list[str],
        varargname: str | None,
    ) -> list[str]:
        """Makes sure args follows the validation given input"""
        out_args = list(args)

        for i in range(min(len(args), len(validator_list))):
            if validator_list[i] is not None:
                out_args[i] = self._handle_argument(
                    args[i], validator_list[i], argnames[i]
                )

        if len(args) > len(validator_list):
            if varargs is not None:
                assert varargname is not None
                out_args[len(validator_list) :] = [
                    self._handle_argument(a, varargs, varargname)
                    for a in args[len(validator_list) :]
                ]
            else:
                out_args[len(validator_list) :] = args[len(validator_list) :]

        return out_args

    def _parse_and_dispatch(self, argv: list[str]) -> tuple[Application, int]:
        """Parses ``argv``, runs switches and ``main()``; returns the final
        instance (the nested subcommand's, if one ran) and the return code."""
        inst: Application = self
        retcode: int | None = 0
        try:
            swfuncs, tailargs = self._parse_args(argv)
            ordered, tailargs = self._validate_args(swfuncs, tailargs)
        except ShowHelp:
            self.help()
        except ShowHelpAll:
            self.helpall()
        except ShowVersion:
            self.version()
        except ShowCompletion:
            info = swfuncs[self.completions.__func__]  # type: ignore[attr-defined]
            self._print_completion(info.val[0])
        except SwitchError as ex:
            print(T_("Error: {0}").format(ex))
            print(T_("------"))
            self.help()
            retcode = 2
        else:
            for f, a in ordered:
                f(self, *a)

            cleanup = None
            if not self.nested_command or self.CALL_MAIN_IF_NESTED_COMMAND:
                retcode = self.main(*tailargs)
                cleanup = functools.partial(self.cleanup, retcode)
            if not retcode and self.nested_command:
                subapp, argv = self.nested_command
                subapp.parent = self
                inst, retcode = subapp.run(argv, exit=False)

            if cleanup:
                cleanup()

        return inst, retcode or 0

    @typing.overload
    @classmethod
    def run(
        cls,
        argv: list[str] | None = ...,
        *,
        exit: Literal[True] = ...,  # pylint: disable=redefined-builtin
    ) -> NoReturn: ...

    @typing.overload
    @classmethod
    def run(
        cls,
        argv: list[str] | None = ...,
        *,
        exit: Literal[False],  # pylint: disable=redefined-builtin
    ) -> tuple[Self, int]: ...

    @classmethod
    def run(
        cls,
        argv: list[str] | None = None,
        *,
        exit: bool = True,  # pylint: disable=redefined-builtin
    ) -> tuple[Self, int]:
        """
        Runs the application, taking the arguments from ``sys.argv`` by default if
        nothing is passed. If ``exit`` is
        ``True`` (the default), the function will exit with the appropriate return code;
        otherwise it will return a tuple of ``(inst, retcode)``, where ``inst`` is the
        application instance created internally by this function and ``retcode`` is the
        exit code of the application.

        .. note::
           Setting ``exit`` to ``False`` is intended for testing/debugging purposes only -- do
           not override it in other situations.
        """
        if argv is None:
            argv = sys.argv
        cls.autocomplete(argv)
        argv = list(argv)
        inst = cls(argv.pop(0))
        try:
            inst, retcode = inst._parse_and_dispatch(argv)  # type: ignore[assignment]
            if exit:
                # surface an EPIPE now, while we can still handle it below
                sys.stdout.flush()
        except OSError as exc:
            # The reader closed the pipe (e.g. output piped to ``head``). On
            # POSIX this is BrokenPipeError (EPIPE); on Windows the flush
            # raises OSError EINVAL instead. Re-raise anything else.
            # Never change the SIGPIPE disposition instead: that would make a
            # socket send() to a closed peer kill the whole process.
            if not isinstance(exc, BrokenPipeError) and exc.errno != errno.EINVAL:
                raise
            retcode = 1
            if exit:
                # Point stdout at devnull so the interpreter's final flush
                # doesn't raise on whatever is still buffered.
                with contextlib.suppress(OSError, ValueError):
                    devnull = os.open(os.devnull, os.O_WRONLY)
                    os.dup2(devnull, sys.stdout.fileno())

        if exit:
            sys.exit(retcode)
        else:
            return inst, retcode

    @classmethod
    def invoke(cls, *args: str, **switches: Any) -> tuple[Application, int]:
        """Invoke this application programmatically (as a function), in the same way ``run()``
        would. There are two key differences: the return value of ``main()`` is not converted to
        an integer (returned as-is), and exceptions are not swallowed either.

        :param args: any positional arguments for ``main()``
        :param switches: command-line switches are passed as keyword arguments,
                         e.g., ``foo=5`` for ``--foo=5``
        """

        inst = cls("")

        swfuncs = inst._parse_kwd_args(switches)
        ordered, tailargs = inst._validate_args(swfuncs, args)
        for f, a in ordered:
            f(inst, *a)

        cleanup = None
        retcode = 0
        if not inst.nested_command or inst.CALL_MAIN_IF_NESTED_COMMAND:
            retcode = inst.main(*tailargs)
            cleanup = functools.partial(inst.cleanup, retcode)
        if not retcode and inst.nested_command:
            subapp, argv = inst.nested_command
            subapp.parent = inst
            inst, retcode = subapp.run(argv, exit=False)

        if cleanup:
            cleanup()

        return inst, retcode

    def _parse_kwd_args(
        self, switches: dict[str, Any]
    ) -> dict[Callable[..., Any], SwitchParseInfo]:
        """Parses keywords (positional arguments), used by invoke."""
        swfuncs = {}
        p: tuple[Any, ...]
        for index, (swname, val) in enumerate(switches.items(), 1):
            switch_local = getattr(type(self), swname)
            swinfo = self._switches_by_func[switch_local._switch_info.func]
            if isinstance(switch_local, CountOf):
                p = (range(val),)
            elif swinfo.list and not hasattr(val, "__iter__"):
                raise SwitchError(
                    T_("Switch {0} must be a sequence (iterable)").format(swname)
                )
            elif not swinfo.argtype:
                # a flag
                # If val is False or None, skip adding the flag (treat as not present)
                if val in (False, None):
                    continue
                if val not in (True, Flag):
                    raise SwitchError(T_("Switch {0} is a boolean flag").format(swname))
                p = ()
            else:
                p = (val,)
            swfuncs[swinfo.func] = SwitchParseInfo(swname, p, index)
        return swfuncs

    @typing.no_type_check
    def main(self, *args: str) -> int:
        """Implement me (no need to call super)"""
        if self._subcommands:
            if args:
                print(T_("Unknown sub-command '{0}'").format(args[0]))
                print(T_("------"))
                self.help()
                return 1
            if not self.nested_command:
                print(T_("No sub-command given"))
                print(T_("------"))
                self.help()
                return 1
            return 0

        print(T_("main() not implemented"))
        return 1

    def cleanup(self, retcode: int) -> None:
        """Called after ``main()`` and all sub-applications have executed, to perform any necessary cleanup.

        :param retcode: the return code of ``main()``
        """

    @switch(
        ["--help-all"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints help messages of all sub-commands and quits"""),
    )
    def helpall(self) -> None:
        """Prints help messages of all sub-commands and quits"""
        self.help()
        print()

        if self._subcommands:
            for name, subcls in sorted(self._subcommands.items()):
                subapp = (subcls.get())(f"{self.PROGNAME} {name}")
                subapp.parent = self
                # Demote meta-switches in the nested help output. This must NOT
                # mutate the shared SwitchInfo objects (which live on the
                # function objects and are reused by every Application in the
                # process); instead flag this instance so help() regroups them
                # locally for this render only.
                subapp._hide_meta_switches = True
                subapp.helpall()

    @switch(
        ["-h", "--help"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints this help message and quits"""),
    )
    def help(self) -> None:  # @ReservedAssignment
        """Prints this help message and quits"""
        if self._get_prog_version():
            self.version()
            print()
        if self.DESCRIPTION:
            print(self.DESCRIPTION.strip() + "\n")

        def split_indentation(s: str) -> tuple[str, str]:
            """Identifies the initial indentation (all spaces) of the string and returns the indentation as well
            as the remainder of the line.
            """
            i = 0
            while i < len(s) and s[i] == " ":
                i += 1
            return s[:i], s[i:]

        def paragraphs(text: str) -> Generator[tuple[str, str, str], None, None]:
            """Yields each paragraph of text along with its initial and subsequent indentations to be used by
            textwrap.TextWrapper.

            Identifies list items from their first non-space character being one of bullets '-', '*', and '/'.
            However, bullet '/' is invisible and is removed from the list item.

            :param text: The text to separate into paragraphs
            """

            paragraph = None
            initial_indent = ""
            subsequent_indent = ""

            def current() -> Generator[tuple[str, str, str], None, None]:
                """Yields the current result if present."""
                if paragraph:
                    yield paragraph, initial_indent, subsequent_indent

            for part in text.lstrip("\n").split("\n"):
                indent, line = split_indentation(part)

                if len(line) == 0:
                    # Starting a new paragraph
                    yield from current()
                    yield "", "", ""

                    paragraph = None
                    initial_indent = ""
                    subsequent_indent = ""
                else:
                    # Adding to current paragraph
                    def is_list_item(line: str) -> bool:
                        """Returns true if the first element of 'line' is a bullet character."""
                        bullets = ["-", "*", "/"]
                        return line[0] in bullets

                    def has_invisible_bullet(line: str) -> bool:
                        """Returns true if the first element of 'line' is the invisible bullet ('/')."""
                        return line[0] == "/"

                    if is_list_item(line):
                        # Done with current paragraph
                        yield from current()

                        if has_invisible_bullet(line):
                            line = line[1:]

                        paragraph = line
                        initial_indent = indent

                        # Calculate extra indentation for subsequent lines of this list item
                        i = 1
                        while i < len(line) and line[i] == " ":
                            i += 1
                        subsequent_indent = indent + " " * i
                    elif not paragraph:
                        # Start a new paragraph
                        paragraph = line
                        initial_indent = indent
                        subsequent_indent = indent
                    else:
                        # Add to current paragraph
                        paragraph = paragraph + " " + line

            yield from current()

        def wrapped_paragraphs(text: str, width: int) -> Generator[str, None, None]:
            """Yields each line of each paragraph of text after wrapping them on 'width' number of columns.

            :param text: The text to yield wrapped lines of
            :param width: The width of the wrapped output
            """
            if not text:
                return

            width = max(width, 1)

            for paragraph, initial_indent, subsequent_indent in paragraphs(text):
                wrapper = TextWrapper(
                    width,
                    initial_indent=initial_indent,
                    subsequent_indent=subsequent_indent,
                )
                w = wrapper.wrap(paragraph)
                yield from w
                if len(w) == 0:
                    yield ""

        cols, _ = get_terminal_size()
        if self.DESCRIPTION_MORE is not None:
            for line in wrapped_paragraphs(self.DESCRIPTION_MORE, cols):
                print(line)

        m = inspect.getfullargspec(self.main)
        tailargs_str = m.args[1:]  # skip self
        if m.defaults:
            for i, d in enumerate(reversed(m.defaults)):
                tailargs_str[-i - 1] = f"[{tailargs_str[-i - 1]}={d}]"
        if m.varargs:
            tailargs_str.append(f"{m.varargs}...")
        tailargs = " ".join(tailargs_str)

        utc = self.COLOR_USAGE_TITLE or self.COLOR_USAGE
        print(utc | T_("Usage:"))

        with self.COLOR_USAGE:
            if not self.USAGE:
                if self._subcommands:
                    self.USAGE = T_(
                        "    {progname} [SWITCHES] [SUBCOMMAND [SWITCHES]] {tailargs}\n"
                    )
                else:
                    self.USAGE = T_("    {progname} [SWITCHES] {tailargs}\n")
            assert not isinstance(self.PROGNAME, Style)
            print(
                self.USAGE.format(
                    progname=colors.filter(self.PROGNAME), tailargs=tailargs
                )
            )

        by_groups: dict[str, list[SwitchInfo]] = {}
        for si in self._switches_by_func.values():
            group = (
                "Hidden-switches"
                if self._hide_meta_switches and si.group == "Meta-switches"
                else si.group
            )
            if group not in by_groups:
                by_groups[group] = []
            by_groups[group].append(si)

        def switchs(
            by_groups: dict[str, list[SwitchInfo]], show_groups: bool
        ) -> Generator[tuple[SwitchInfo, str, Style], None, None]:
            for grp, swinfos in sorted(by_groups.items(), key=lambda item: item[0]):
                if show_groups:
                    lgrp = T_(grp) if grp in _switch_groups else grp
                    print(self.COLOR_GROUP_TITLES[grp] | lgrp + ":")

                for si in sorted(swinfos, key=lambda x: x.names):
                    swnames = ", ".join(
                        ("-" if len(n) == 1 else "--") + n
                        for n in si.names
                        if n in self._switches_by_name
                        and self._switches_by_name[n] == si
                    )
                    if si.argtype:
                        if hasattr(si.argtype, "__name__"):
                            typename = si.argtype.__name__
                        else:
                            typename = str(si.argtype)
                        argtype = f" {si.argname.upper()}:{typename}"
                    else:
                        argtype = ""
                    prefix = swnames + argtype
                    yield si, prefix, self.COLOR_GROUPS[grp]

                if show_groups:
                    print()

        sw_width = (
            max(len(prefix) for si, prefix, color in switchs(by_groups, False)) + 4
        )
        description_indent = "    {0}{1}{2}"
        wrapper = TextWrapper(width=max(cols - min(sw_width, 60), 50) - 6)
        indentation = "\n" + " " * (cols - wrapper.width)

        for switch_info, prefix, color in switchs(by_groups, True):
            help_txt = switch_info.help or ""
            if switch_info.list:
                help_txt += T_("; may be given multiple times")
            if switch_info.mandatory:
                help_txt += T_("; required")
            if switch_info.requires:
                help_txt += T_("; requires {0}").format(
                    ", ".join(
                        (("-" if len(switch) == 1 else "--") + switch)
                        for switch in switch_info.requires
                    )
                )
            if switch_info.excludes:
                help_txt += T_("; excludes {0}").format(
                    ", ".join(
                        (("-" if len(switch) == 1 else "--") + switch)
                        for switch in switch_info.excludes
                    )
                )

            msg = indentation.join(
                wrapper.wrap(" ".join(ln.strip() for ln in help_txt.splitlines()))
            )

            if len(prefix) + wrapper.width >= cols:
                padding = indentation
            else:
                padding = " " * max(cols - wrapper.width - len(prefix) - 4, 1)
            print(description_indent.format(color | prefix, padding, color | msg))

        if self._subcommands:
            gc = self.COLOR_GROUP_TITLES["Sub-commands"]
            print(gc | T_("Sub-commands:"))
            for name, subcls in sorted(self._subcommands.items()):
                with gc:
                    subapp = subcls.get()
                    doc = subapp.DESCRIPTION or getdoc(subapp)
                    if self.SUBCOMMAND_HELPMSG:
                        help_str = doc + "; " if doc else ""
                        help_str += self.SUBCOMMAND_HELPMSG.format(
                            parent=self.PROGNAME, sub=name
                        )
                    else:
                        help_str = doc or ""

                    msg = indentation.join(
                        wrapper.wrap(
                            " ".join(ln.strip() for ln in help_str.splitlines())
                        )
                    )

                    if len(name) + wrapper.width >= cols:
                        padding = indentation
                    else:
                        padding = " " * max(cols - wrapper.width - len(name) - 4, 1)

                    bodycolor = (
                        colors.extract(subcls.name)
                        if colors.contains_colors(subcls.name)
                        else gc
                    )

                    print(
                        description_indent.format(
                            subcls.name, padding, bodycolor | colors.filter(msg)
                        )
                    )

    def _get_prog_version(self) -> str | None:
        ver: str | None = None
        curr: Application | None = self
        while curr is not None:
            ver = getattr(curr, "VERSION", None)
            if ver is not None:
                return ver  # type: ignore[no-any-return]
            curr = curr.parent
        return ver

    @switch(
        ["-v", "--version"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints the program's version and quits"""),
    )
    def version(self) -> None:
        """Prints the program's version and quits"""
        ver = self._get_prog_version()
        ver_name = ver if ver is not None else T_("(version not set)")
        print(f"{self.PROGNAME} {ver_name}")

    def _get_prog_name_for_completion(self) -> str | None:
        """Get program name for completion scripts, or None if not a root app."""
        if self.parent is not None:
            return None
        return (
            os.path.basename(self.executable) if self.executable else str(self.PROGNAME)
        )

    def _print_completion(self, shell: str) -> None:
        """Print completion script for the specified shell."""
        prog_name = self._get_prog_name_for_completion()
        if prog_name is None:
            print(T_("Shell completion must be generated from the root application"))
            return

        if shell == "bash":
            print(self.completion_script_bash(prog_name))
        elif shell == "fish":
            print(self.completion_script_fish(prog_name))
        else:
            print(T_("Unsupported shell: {shell}").format(shell=shell))
            print(T_("Supported shells are: bash, fish"))

    @switch(
        ["--completions"],
        Set("bash", "fish"),
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints shell completion script and quits"""),
    )
    def completions(self, shell: str) -> None:
        """Prints shell completion script and quits"""
        # Handled in _validate_args via ShowCompletion


__all__ = [
    "Application",
    "ShowCompletion",
    "ShowHelp",
    "ShowHelpAll",
    "ShowVersion",
]


def __dir__() -> list[str]:
    return list(__all__)
