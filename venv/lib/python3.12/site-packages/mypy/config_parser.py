from __future__ import annotations

import argparse
import configparser
import glob as fileglob
import os
import re
import sys
from io import StringIO

from mypy.errorcodes import error_codes

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from collections.abc import Mapping, MutableMapping, Sequence
from typing import Any, Callable, Final, TextIO, Union
from typing_extensions import TypeAlias as _TypeAlias

from mypy import defaults
from mypy.options import PER_MODULE_OPTIONS, Options

_CONFIG_VALUE_TYPES: _TypeAlias = Union[
    str, bool, int, float, dict[str, str], list[str], tuple[int, int]
]
_INI_PARSER_CALLABLE: _TypeAlias = Callable[[Any], _CONFIG_VALUE_TYPES]


class VersionTypeError(argparse.ArgumentTypeError):
    """Provide a fallback value if the Python version is unsupported."""

    def __init__(self, *args: Any, fallback: tuple[int, int]) -> None:
        self.fallback = fallback
        super().__init__(*args)


def parse_version(v: str | float) -> tuple[int, int]:
    m = re.match(r"\A(\d)\.(\d+)\Z", str(v))
    if not m:
        raise argparse.ArgumentTypeError(f"Invalid python version '{v}' (expected format: 'x.y')")
    major, minor = int(m.group(1)), int(m.group(2))
    if major == 2 and minor == 7:
        pass  # Error raised elsewhere
    elif major == 3:
        if minor < defaults.PYTHON3_VERSION_MIN[1]:
            msg = "Python 3.{} is not supported (must be {}.{} or higher)".format(
                minor, *defaults.PYTHON3_VERSION_MIN
            )

            if isinstance(v, float):
                msg += ". You may need to put quotes around your Python version"

            raise VersionTypeError(msg, fallback=defaults.PYTHON3_VERSION_MIN)
    else:
        raise argparse.ArgumentTypeError(
            f"Python major version '{major}' out of range (must be 3)"
        )
    return major, minor


def try_split(v: str | Sequence[str], split_regex: str = "[,]") -> list[str]:
    """Split and trim a str or list of str into a list of str"""
    if isinstance(v, str):
        items = [p.strip() for p in re.split(split_regex, v)]
        if items and items[-1] == "":
            items.pop(-1)
        return items
    return [p.strip() for p in v]


def validate_codes(codes: list[str]) -> list[str]:
    invalid_codes = set(codes) - set(error_codes.keys())
    if invalid_codes:
        raise argparse.ArgumentTypeError(
            f"Invalid error code(s): {', '.join(sorted(invalid_codes))}"
        )
    return codes


def validate_package_allow_list(allow_list: list[str]) -> list[str]:
    for p in allow_list:
        msg = f"Invalid allow list entry: {p}"
        if "*" in p:
            raise argparse.ArgumentTypeError(
                f"{msg} (entries are already prefixes so must not contain *)"
            )
        if "\\" in p or "/" in p:
            raise argparse.ArgumentTypeError(
                f"{msg} (entries must be packages like foo.bar not directories or files)"
            )
    return allow_list


def expand_path(path: str) -> str:
    """Expand the user home directory and any environment variables contained within
    the provided path.
    """

    return os.path.expandvars(os.path.expanduser(path))


def str_or_array_as_list(v: str | Sequence[str]) -> list[str]:
    if isinstance(v, str):
        return [v.strip()] if v.strip() else []
    return [p.strip() for p in v if p.strip()]


def split_and_match_files_list(paths: Sequence[str]) -> list[str]:
    """Take a list of files/directories (with support for globbing through the glob library).

    Where a path/glob matches no file, we still include the raw path in the resulting list.

    Returns a list of file paths
    """
    expanded_paths = []

    for path in paths:
        path = expand_path(path.strip())
        globbed_files = fileglob.glob(path, recursive=True)
        if globbed_files:
            expanded_paths.extend(globbed_files)
        else:
            expanded_paths.append(path)

    return expanded_paths


def split_and_match_files(paths: str) -> list[str]:
    """Take a string representing a list of files/directories (with support for globbing
    through the glob library).

    Where a path/glob matches no file, we still include the raw path in the resulting list.

    Returns a list of file paths
    """

    return split_and_match_files_list(split_commas(paths))


def check_follow_imports(choice: str) -> str:
    choices = ["normal", "silent", "skip", "error"]
    if choice not in choices:
        raise argparse.ArgumentTypeError(
            "invalid choice '{}' (choose from {})".format(
                choice, ", ".join(f"'{x}'" for x in choices)
            )
        )
    return choice


def check_junit_format(choice: str) -> str:
    choices = ["global", "per_file"]
    if choice not in choices:
        raise argparse.ArgumentTypeError(
            "invalid choice '{}' (choose from {})".format(
                choice, ", ".join(f"'{x}'" for x in choices)
            )
        )
    return choice


def split_commas(value: str) -> list[str]:
    # Uses a bit smarter technique to allow last trailing comma
    # and to remove last `""` item from the split.
    items = value.split(",")
    if items and items[-1] == "":
        items.pop(-1)
    return items


# For most options, the type of the default value set in options.py is
# sufficient, and we don't have to do anything here.  This table
# exists to specify types for values initialized to None or container
# types.
ini_config_types: Final[dict[str, _INI_PARSER_CALLABLE]] = {
    "python_version": parse_version,
    "custom_typing_module": str,
    "custom_typeshed_dir": expand_path,
    "mypy_path": lambda s: [expand_path(p.strip()) for p in re.split("[,:]", s)],
    "files": split_and_match_files,
    "quickstart_file": expand_path,
    "junit_xml": expand_path,
    "junit_format": check_junit_format,
    "follow_imports": check_follow_imports,
    "no_site_packages": bool,
    "plugins": lambda s: [p.strip() for p in split_commas(s)],
    "always_true": lambda s: [p.strip() for p in split_commas(s)],
    "always_false": lambda s: [p.strip() for p in split_commas(s)],
    "untyped_calls_exclude": lambda s: validate_package_allow_list(
        [p.strip() for p in split_commas(s)]
    ),
    "enable_incomplete_feature": lambda s: [p.strip() for p in split_commas(s)],
    "disable_error_code": lambda s: validate_codes([p.strip() for p in split_commas(s)]),
    "enable_error_code": lambda s: validate_codes([p.strip() for p in split_commas(s)]),
    "package_root": lambda s: [p.strip() for p in split_commas(s)],
    "cache_dir": expand_path,
    "python_executable": expand_path,
    "strict": bool,
    "exclude": lambda s: [s.strip()],
    "packages": try_split,
    "modules": try_split,
}

# Reuse the ini_config_types and overwrite the diff
toml_config_types: Final[dict[str, _INI_PARSER_CALLABLE]] = ini_config_types.copy()
toml_config_types.update(
    {
        "python_version": parse_version,
        "mypy_path": lambda s: [expand_path(p) for p in try_split(s, "[,:]")],
        "files": lambda s: split_and_match_files_list(try_split(s)),
        "junit_format": lambda s: check_junit_format(str(s)),
        "follow_imports": lambda s: check_follow_imports(str(s)),
        "plugins": try_split,
        "always_true": try_split,
        "always_false": try_split,
        "untyped_calls_exclude": lambda s: validate_package_allow_list(try_split(s)),
        "enable_incomplete_feature": try_split,
        "disable_error_code": lambda s: validate_codes(try_split(s)),
        "enable_error_code": lambda s: validate_codes(try_split(s)),
        "package_root": try_split,
        "exclude": str_or_array_as_list,
        "packages": try_split,
        "modules": try_split,
    }
)


def _parse_individual_file(
    config_file: str, stderr: TextIO | None = None
) -> tuple[MutableMapping[str, Any], dict[str, _INI_PARSER_CALLABLE], str] | None:

    if not os.path.exists(config_file):
        return None

    parser: MutableMapping[str, Any]
    try:
        if is_toml(config_file):
            with open(config_file, "rb") as f:
                toml_data = tomllib.load(f)
            # Filter down to just mypy relevant toml keys
            toml_data = toml_data.get("tool", {})
            if "mypy" not in toml_data:
                return None
            toml_data = {"mypy": toml_data["mypy"]}
            parser = destructure_overrides(toml_data)
            config_types = toml_config_types
        else:
            parser = configparser.RawConfigParser()
            parser.read(config_file)
            config_types = ini_config_types

    except (tomllib.TOMLDecodeError, configparser.Error, ConfigTOMLValueError) as err:
        print(f"{config_file}: {err}", file=stderr)
        return None

    if os.path.basename(config_file) in defaults.SHARED_CONFIG_NAMES and "mypy" not in parser:
        return None

    return parser, config_types, config_file


def _find_config_file(
    stderr: TextIO | None = None,
) -> tuple[MutableMapping[str, Any], dict[str, _INI_PARSER_CALLABLE], str] | None:

    current_dir = os.path.abspath(os.getcwd())

    while True:
        for name in defaults.CONFIG_NAMES + defaults.SHARED_CONFIG_NAMES:
            config_file = os.path.relpath(os.path.join(current_dir, name))
            ret = _parse_individual_file(config_file, stderr)
            if ret is None:
                continue
            return ret

        if any(
            os.path.exists(os.path.join(current_dir, cvs_root)) for cvs_root in (".git", ".hg")
        ):
            break
        parent_dir = os.path.dirname(current_dir)
        if parent_dir == current_dir:
            break
        current_dir = parent_dir

    for config_file in defaults.USER_CONFIG_FILES:
        ret = _parse_individual_file(config_file, stderr)
        if ret is None:
            continue
        return ret

    return None


def parse_config_file(
    options: Options,
    set_strict_flags: Callable[[], None],
    filename: str | None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> None:
    """Parse a config file into an Options object.

    Errors are written to stderr but are not fatal.

    If filename is None, fall back to default config files.
    """
    stdout = stdout or sys.stdout
    stderr = stderr or sys.stderr

    ret = (
        _parse_individual_file(filename, stderr)
        if filename is not None
        else _find_config_file(stderr)
    )
    if ret is None:
        return
    parser, config_types, file_read = ret

    options.config_file = file_read
    os.environ["MYPY_CONFIG_FILE_DIR"] = os.path.dirname(os.path.abspath(file_read))

    if "mypy" not in parser:
        if filename or os.path.basename(file_read) not in defaults.SHARED_CONFIG_NAMES:
            print(f"{file_read}: No [mypy] section in config file", file=stderr)
    else:
        section = parser["mypy"]
        prefix = f"{file_read}: [mypy]: "
        updates, report_dirs = parse_section(
            prefix, options, set_strict_flags, section, config_types, stderr
        )
        for k, v in updates.items():
            setattr(options, k, v)
        options.report_dirs.update(report_dirs)

    for name, section in parser.items():
        if name.startswith("mypy-"):
            prefix = get_prefix(file_read, name)
            updates, report_dirs = parse_section(
                prefix, options, set_strict_flags, section, config_types, stderr
            )
            if report_dirs:
                print(
                    prefix,
                    "Per-module sections should not specify reports ({})".format(
                        ", ".join(s + "_report" for s in sorted(report_dirs))
                    ),
                    file=stderr,
                )
            if set(updates) - PER_MODULE_OPTIONS:
                print(
                    prefix,
                    "Per-module sections should only specify per-module flags ({})".format(
                        ", ".join(sorted(set(updates) - PER_MODULE_OPTIONS))
                    ),
                    file=stderr,
                )
                updates = {k: v for k, v in updates.items() if k in PER_MODULE_OPTIONS}

            globs = name[5:]
            for glob in globs.split(","):
                # For backwards compatibility, replace (back)slashes with dots.
                glob = glob.replace(os.sep, ".")
                if os.altsep:
                    glob = glob.replace(os.altsep, ".")

                if any(c in glob for c in "?[]!") or any(
                    "*" in x and x != "*" for x in glob.split(".")
                ):
                    print(
                        prefix,
                        "Patterns must be fully-qualified module names, optionally "
                        "with '*' in some components (e.g spam.*.eggs.*)",
                        file=stderr,
                    )
                else:
                    options.per_module_options[glob] = updates


def get_prefix(file_read: str, name: str) -> str:
    if is_toml(file_read):
        module_name_str = 'module = "%s"' % "-".join(name.split("-")[1:])
    else:
        module_name_str = name

    return f"{file_read}: [{module_name_str}]:"


def is_toml(filename: str) -> bool:
    return filename.lower().endswith(".toml")


def destructure_overrides(toml_data: dict[str, Any]) -> dict[str, Any]:
    """Take the new [[tool.mypy.overrides]] section array in the pyproject.toml file,
    and convert it back to a flatter structure that the existing config_parser can handle.

    E.g. the following pyproject.toml file:

        [[tool.mypy.overrides]]
        module = [
            "a.b",
            "b.*"
        ]
        disallow_untyped_defs = true

        [[tool.mypy.overrides]]
        module = 'c'
        disallow_untyped_defs = false

    Would map to the following config dict that it would have gotten from parsing an equivalent
    ini file:

        {
            "mypy-a.b": {
                disallow_untyped_defs = true,
            },
            "mypy-b.*": {
                disallow_untyped_defs = true,
            },
            "mypy-c": {
                disallow_untyped_defs: false,
            },
        }
    """
    if "overrides" not in toml_data["mypy"]:
        return toml_data

    if not isinstance(toml_data["mypy"]["overrides"], list):
        raise ConfigTOMLValueError(
            "tool.mypy.overrides sections must be an array. Please make "
            "sure you are using double brackets like so: [[tool.mypy.overrides]]"
        )

    result = toml_data.copy()
    for override in result["mypy"]["overrides"]:
        if "module" not in override:
            raise ConfigTOMLValueError(
                "toml config file contains a [[tool.mypy.overrides]] "
                "section, but no module to override was specified."
            )

        if isinstance(override["module"], str):
            modules = [override["module"]]
        elif isinstance(override["module"], list):
            modules = override["module"]
        else:
            raise ConfigTOMLValueError(
                "toml config file contains a [[tool.mypy.overrides]] "
                "section with a module value that is not a string or a list of "
                "strings"
            )

        for module in modules:
            module_overrides = override.copy()
            del module_overrides["module"]
            old_config_name = f"mypy-{module}"
            if old_config_name not in result:
                result[old_config_name] = module_overrides
            else:
                for new_key, new_value in module_overrides.items():
                    if (
                        new_key in result[old_config_name]
                        and result[old_config_name][new_key] != new_value
                    ):
                        raise ConfigTOMLValueError(
                            "toml config file contains "
                            "[[tool.mypy.overrides]] sections with conflicting "
                            f"values. Module '{module}' has two different values for '{new_key}'"
                        )
                    result[old_config_name][new_key] = new_value

    del result["mypy"]["overrides"]
    return result


def parse_section(
    prefix: str,
    template: Options,
    set_strict_flags: Callable[[], None],
    section: Mapping[str, Any],
    config_types: dict[str, Any],
    stderr: TextIO = sys.stderr,
) -> tuple[dict[str, object], dict[str, str]]:
    """Parse one section of a config file.

    Returns a dict of option values encountered, and a dict of report directories.
    """
    results: dict[str, object] = {}
    report_dirs: dict[str, str] = {}

    # Because these fields exist on Options, without proactive checking, we would accept them
    # and crash later
    invalid_options = {
        "enabled_error_codes": "enable_error_code",
        "disabled_error_codes": "disable_error_code",
    }

    for key in section:
        invert = False
        options_key = key
        if key in config_types:
            ct = config_types[key]
        elif key in invalid_options:
            print(
                f"{prefix}Unrecognized option: {key} = {section[key]}"
                f" (did you mean {invalid_options[key]}?)",
                file=stderr,
            )
            continue
        else:
            dv = getattr(template, key, None)
            if dv is None:
                if key.endswith("_report"):
                    report_type = key[:-7].replace("_", "-")
                    if report_type in defaults.REPORTER_NAMES:
                        report_dirs[report_type] = str(section[key])
                    else:
                        print(f"{prefix}Unrecognized report type: {key}", file=stderr)
                    continue
                if key.startswith("x_"):
                    pass  # Don't complain about `x_blah` flags
                elif key.startswith("no_") and hasattr(template, key[3:]):
                    options_key = key[3:]
                    invert = True
                elif key.startswith("allow") and hasattr(template, "dis" + key):
                    options_key = "dis" + key
                    invert = True
                elif key.startswith("disallow") and hasattr(template, key[3:]):
                    options_key = key[3:]
                    invert = True
                elif key.startswith("show_") and hasattr(template, "hide_" + key[5:]):
                    options_key = "hide_" + key[5:]
                    invert = True
                elif key == "strict":
                    pass  # Special handling below
                else:
                    print(f"{prefix}Unrecognized option: {key} = {section[key]}", file=stderr)
                if invert:
                    dv = getattr(template, options_key, None)
                else:
                    continue
            ct = type(dv)
        v: Any = None
        try:
            if ct is bool:
                if isinstance(section, dict):
                    v = convert_to_boolean(section.get(key))
                else:
                    v = section.getboolean(key)  # type: ignore[attr-defined]  # Until better stub
                if invert:
                    v = not v
            elif callable(ct):
                if invert:
                    print(f"{prefix}Can not invert non-boolean key {options_key}", file=stderr)
                    continue
                try:
                    v = ct(section.get(key))
                except VersionTypeError as err_version:
                    print(f"{prefix}{key}: {err_version}", file=stderr)
                    v = err_version.fallback
                except argparse.ArgumentTypeError as err:
                    print(f"{prefix}{key}: {err}", file=stderr)
                    continue
            else:
                print(f"{prefix}Don't know what type {key} should have", file=stderr)
                continue
        except ValueError as err:
            print(f"{prefix}{key}: {err}", file=stderr)
            continue
        if key == "strict":
            if v:
                set_strict_flags()
            continue
        results[options_key] = v

    # These two flags act as per-module overrides, so store the empty defaults.
    if "disable_error_code" not in results:
        results["disable_error_code"] = []
    if "enable_error_code" not in results:
        results["enable_error_code"] = []

    return results, report_dirs


def convert_to_boolean(value: Any | None) -> bool:
    """Return a boolean value translating from other types if necessary."""
    if isinstance(value, bool):
        return value
    if not isinstance(value, str):
        value = str(value)
    if value.lower() not in configparser.RawConfigParser.BOOLEAN_STATES:
        raise ValueError(f"Not a boolean: {value}")
    return configparser.RawConfigParser.BOOLEAN_STATES[value.lower()]


def split_directive(s: str) -> tuple[list[str], list[str]]:
    """Split s on commas, except during quoted sections.

    Returns the parts and a list of error messages."""
    parts = []
    cur: list[str] = []
    errors = []
    i = 0
    while i < len(s):
        if s[i] == ",":
            parts.append("".join(cur).strip())
            cur = []
        elif s[i] == '"':
            i += 1
            while i < len(s) and s[i] != '"':
                cur.append(s[i])
                i += 1
            if i == len(s):
                errors.append("Unterminated quote in configuration comment")
                cur.clear()
        else:
            cur.append(s[i])
        i += 1
    if cur:
        parts.append("".join(cur).strip())

    return parts, errors


def mypy_comments_to_config_map(line: str, template: Options) -> tuple[dict[str, str], list[str]]:
    """Rewrite the mypy comment syntax into ini file syntax."""
    options = {}
    entries, errors = split_directive(line)
    for entry in entries:
        if "=" not in entry:
            name = entry
            value = None
        else:
            name, value = (x.strip() for x in entry.split("=", 1))

        name = name.replace("-", "_")
        if value is None:
            value = "True"
        options[name] = value

    return options, errors


def parse_mypy_comments(
    args: list[tuple[int, str]], template: Options
) -> tuple[dict[str, object], list[tuple[int, str]]]:
    """Parse a collection of inline mypy: configuration comments.

    Returns a dictionary of options to be applied and a list of error messages
    generated.
    """

    errors: list[tuple[int, str]] = []
    sections = {}

    for lineno, line in args:
        # In order to easily match the behavior for bools, we abuse configparser.
        # Oddly, the only way to get the SectionProxy object with the getboolean
        # method is to create a config parser.
        parser = configparser.RawConfigParser()
        options, parse_errors = mypy_comments_to_config_map(line, template)

        if "python_version" in options:
            errors.append((lineno, "python_version not supported in inline configuration"))
            del options["python_version"]

        parser["dummy"] = options
        errors.extend((lineno, x) for x in parse_errors)

        stderr = StringIO()
        strict_found = False

        def set_strict_flags() -> None:
            nonlocal strict_found
            strict_found = True

        new_sections, reports = parse_section(
            "", template, set_strict_flags, parser["dummy"], ini_config_types, stderr=stderr
        )
        errors.extend((lineno, x) for x in stderr.getvalue().strip().split("\n") if x)
        if reports:
            errors.append((lineno, "Reports not supported in inline configuration"))
        if strict_found:
            errors.append(
                (
                    lineno,
                    'Setting "strict" not supported in inline configuration: specify it in '
                    "a configuration file instead, or set individual inline flags "
                    '(see "mypy -h" for the list of flags enabled in strict mode)',
                )
            )

        sections.update(new_sections)

    return sections, errors


def get_config_module_names(filename: str | None, modules: list[str]) -> str:
    if not filename or not modules:
        return ""

    if not is_toml(filename):
        return ", ".join(f"[mypy-{module}]" for module in modules)

    return "module = ['%s']" % ("', '".join(sorted(modules)))


class ConfigTOMLValueError(ValueError):
    pass
