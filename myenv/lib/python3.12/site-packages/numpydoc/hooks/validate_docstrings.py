"""Run numpydoc validation on contents of a file."""

import argparse
import ast
import configparser
import os
import re
import sys

try:
    import tomllib
except ImportError:
    import tomli as tomllib

from pathlib import Path
from typing import Sequence, Tuple, Union

from tabulate import tabulate

from .. import docscrape, validate
from .utils import find_project_root


class AstValidator(validate.Validator):
    """
    Overrides the :class:`Validator` to work entirely with the AST.

    Parameters
    ----------
    ast_node : ast.AST
        The node under inspection.
    filename : os.PathLike
        The file where the node is defined.
    obj_name : str
        A name for the node to use in the listing of issues for the file as a whole.
    """

    def __init__(
        self, *, ast_node: ast.AST, filename: os.PathLike, obj_name: str
    ) -> None:
        self.node: ast.AST = ast_node
        self.raw_doc: str = ast.get_docstring(self.node, clean=False) or ""
        self.clean_doc: str = ast.get_docstring(self.node, clean=True)
        self.doc: docscrape.NumpyDocString = docscrape.NumpyDocString(self.raw_doc)

        self._source_file: os.PathLike = Path(filename).resolve()
        self._name: str = obj_name

        self.is_class: bool = isinstance(ast_node, ast.ClassDef)
        self.is_module: bool = isinstance(ast_node, ast.Module)

    @staticmethod
    def _load_obj(name):
        raise NotImplementedError("AstValidator does not support this method.")

    @property
    def name(self) -> str:
        return self._name

    @property
    def is_function_or_method(self) -> bool:
        return isinstance(self.node, (ast.FunctionDef, ast.AsyncFunctionDef))

    @property
    def is_generator_function(self) -> bool:
        if not self.is_function_or_method:
            return False
        for child in ast.iter_child_nodes(self.node):
            if isinstance(child, ast.Expr) and isinstance(child.value, ast.Yield):
                return True
        return False

    @property
    def type(self) -> str:
        if self.is_function_or_method:
            return "function"
        if self.is_class:
            return "type"
        if self.is_module:
            return "module"
        raise ValueError("Unknown type.")

    @property
    def source_file_name(self) -> str:
        return self._source_file

    @property
    def source_file_def_line(self) -> int:
        return self.node.lineno if not self.is_module else 1

    @property
    def signature_parameters(self) -> Tuple[str]:
        def extract_signature(node):
            args_node = node.args
            params = []
            for arg_type in ["posonlyargs", "args", "vararg", "kwonlyargs", "kwarg"]:
                entries = getattr(args_node, arg_type)
                if arg_type in ["vararg", "kwarg"]:
                    if entries and arg_type == "vararg":
                        params.append(f"*{entries.arg}")
                    if entries and arg_type == "kwarg":
                        params.append(f"**{entries.arg}")
                else:
                    params.extend([arg.arg for arg in entries])
            params = tuple(params)
            if params and params[0] in {"self", "cls"}:
                return params[1:]
            return params

        params = tuple()
        if self.is_function_or_method:
            params = extract_signature(self.node)
        elif self.is_class:
            for child in self.node.body:
                if isinstance(child, ast.FunctionDef) and child.name == "__init__":
                    params = extract_signature(child)
        return params

    @property
    def method_source(self) -> str:
        with open(self.source_file_name) as file:
            source = ast.get_source_segment(file.read(), self.node)
        return source


class DocstringVisitor(ast.NodeVisitor):
    """
    Visits nodes in the AST from a given module and reporting numpydoc issues.

    Parameters
    ----------
    filepath : str
        The absolute or relative path to the file to inspect.
    config : dict
        Configuration options for reviewing flagged issues.
    """

    def __init__(
        self,
        filepath: str,
        config: dict,
    ) -> None:
        self.config: dict = config
        self.filepath: str = filepath
        self.module_name: str = Path(self.filepath).stem
        self.stack: list[str] = []
        self.findings: list = []

    def _ignore_issue(self, node: ast.AST, check: str) -> bool:
        """
        Check whether the issue should be ignored.

        Parameters
        ----------
        node : ast.AST
            The node under inspection.
        check : str
            The code for the check being evaluated.

        Return
        ------
        bool
            Whether the issue should be excluded from the report.
        """
        if check not in self.config["checks"]:
            return True

        if self.config["overrides"]:
            try:
                pattern = self.config["overrides"][check]
                if re.search(pattern, ast.get_docstring(node)) is not None:
                    return True
            except KeyError:
                pass

        return False

    def _get_numpydoc_issues(self, node: ast.AST) -> None:
        """
        Get numpydoc validation issues.

        Parameters
        ----------
        node : ast.AST
            The node under inspection.
        """
        name = ".".join(self.stack)
        report = validate.validate(
            name, AstValidator, ast_node=node, filename=self.filepath
        )
        self.findings.extend(
            [
                [f'{self.filepath}:{report["file_line"]}', name, check, description]
                for check, description in report["errors"]
                if not self._ignore_issue(node, check)
            ]
        )

    def visit(self, node: ast.AST) -> None:
        """
        Visit a node in the AST and report on numpydoc validation issues.

        Parameters
        ----------
        node : ast.AST
            The node to visit.
        """
        if isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            self.stack.append(
                self.module_name if isinstance(node, ast.Module) else node.name
            )

            if not (
                self.config["exclude"]
                and re.search(self.config["exclude"], ".".join(self.stack))
            ):
                self._get_numpydoc_issues(node)

            self.generic_visit(node)
            _ = self.stack.pop()


def parse_config(dir_path: os.PathLike = None) -> dict:
    """
    Parse config information from a pyproject.toml or setup.cfg file.

    This function looks in the provided directory path first for a
    pyproject.toml file. If it finds that, it won't look for a setup.cfg
    file.

    Parameters
    ----------
    dir_path : os.PathLike
        An absolute or relative path to a directory containing
        either a pyproject.toml file specifying a
        [tool.numpydoc_validation] section or a setup.cfg file
        specifying a [tool:numpydoc_validation] section.
        For example, ``~/my_project``. If not provided, the hook
        will try to find the project root directory.

    Returns
    -------
    dict
        Config options for the numpydoc validation hook.
    """
    options = {"checks": {"all"}, "exclude": set(), "overrides": {}}
    dir_path = Path(dir_path).expanduser().resolve()

    toml_path = dir_path / "pyproject.toml"
    cfg_path = dir_path / "setup.cfg"

    def compile_regex(expressions):
        return (
            re.compile(r"|".join(exp for exp in expressions if exp))
            if expressions
            else None
        )

    def extract_check_overrides(options, config_items):
        for option, value in config_items:
            if option.startswith("override_"):
                _, check = option.split("_")
                if value:
                    options["overrides"][check.upper()] = compile_regex(value)

    if toml_path.is_file():
        with open(toml_path, "rb") as toml_file:
            pyproject_toml = tomllib.load(toml_file)
            config = pyproject_toml.get("tool", {}).get("numpydoc_validation", {})
            options["checks"] = set(config.get("checks", options["checks"]))

            global_exclusions = config.get("exclude", options["exclude"])
            options["exclude"] = set(
                global_exclusions
                if not isinstance(global_exclusions, str)
                else [global_exclusions]
            )

            extract_check_overrides(options, config.items())

    elif cfg_path.is_file():
        config = configparser.ConfigParser()
        config.read(cfg_path)
        numpydoc_validation_config_section = "tool:numpydoc_validation"
        try:
            try:
                options["checks"] = set(
                    config.get(numpydoc_validation_config_section, "checks")
                    .rstrip(",")
                    .split(",")
                    or options["checks"]
                )
            except configparser.NoOptionError:
                pass
            try:
                options["exclude"] = set(
                    config.get(numpydoc_validation_config_section, "exclude")
                    .rstrip(",")
                    .split(",")
                    or options["exclude"]
                )
            except configparser.NoOptionError:
                pass

            extract_check_overrides(
                options, config.items(numpydoc_validation_config_section)
            )

        except configparser.NoSectionError:
            pass

    options["checks"] = validate.get_validation_checks(options["checks"])
    options["exclude"] = compile_regex(options["exclude"])
    return options


def process_file(filepath: os.PathLike, config: dict) -> "list[list[str]]":
    """
    Run numpydoc validation on a file.

    Parameters
    ----------
    filepath : path-like
        The absolute or relative path to the file to inspect.
    config : dict
        Configuration options for reviewing flagged issues.

    Returns
    -------
    list[list[str]]
        A list of [name, check, description] lists for flagged issues.
    """
    with open(filepath) as file:
        module_node = ast.parse(file.read(), filepath)

    docstring_visitor = DocstringVisitor(filepath=str(filepath), config=config)
    docstring_visitor.visit(module_node)

    return docstring_visitor.findings


def main(argv: Union[Sequence[str], None] = None) -> int:
    """Run the numpydoc validation hook."""

    project_root_from_cwd, config_file = find_project_root(["."])
    config_options = parse_config(project_root_from_cwd)
    ignored_checks = (
        "\n  "
        + "\n  ".join(
            [
                f"- {check}: {validate.ERROR_MSGS[check]}"
                for check in set(validate.ERROR_MSGS.keys()) - config_options["checks"]
            ]
        )
        + "\n"
    )

    parser = argparse.ArgumentParser(
        description="Run numpydoc validation on files with option to ignore individual checks.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "files", type=str, nargs="+", help="File(s) to run numpydoc validation on."
    )
    parser.add_argument(
        "--config",
        type=str,
        help=(
            "Path to a directory containing a pyproject.toml or setup.cfg file.\n"
            "The hook will look for it in the root project directory.\n"
            "If both are present, only pyproject.toml will be used.\n"
            "Options must be placed under\n"
            "    - [tool:numpydoc_validation] for setup.cfg files and\n"
            "    - [tool.numpydoc_validation] for pyproject.toml files."
        ),
    )
    parser.add_argument(
        "--ignore",
        type=str,
        nargs="*",
        help=(
            f"""Check codes to ignore.{
                ' Currently ignoring the following from '
                f'{Path(project_root_from_cwd) / config_file}: {ignored_checks}'
                'Values provided here will be in addition to the above, unless an alternate config is provided.'
                if config_options["checks"] else ''
            }"""
        ),
    )

    args = parser.parse_args(argv)
    project_root, _ = find_project_root(args.files)
    config_options = parse_config(args.config or project_root)
    config_options["checks"] -= set(args.ignore or [])

    findings = []
    for file in args.files:
        findings.extend(process_file(file, config_options))

    if findings:
        print(
            tabulate(
                findings,
                headers=["file", "item", "check", "description"],
                tablefmt="grid",
                maxcolwidths=50,
            ),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
