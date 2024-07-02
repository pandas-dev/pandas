"""

Enforce that all usages of tm.assert_produces_warning use
the "match" argument. This will help ensure that users always see
the correct warning message.

tm.assert_produces_warning(None) is excluded as no warning is
produced.

"""


import argparse
import ast
from collections.abc import Sequence
import sys

ERROR_MESSAGE = (
    "{path}:{lineno}:{col_offset}: "
    '"match" argument missing in tm.assert_produces_warning'
    "\n"
)


class MatchArgForWarningsChecker(ast.NodeVisitor):
    def __init__(self) -> None:
        self.error_set = []

    def visit_Call(self, node) -> None:
        if ( isinstance(node.func, ast.Attribute)  and
            node.func.attr == "assert_produces_warning"):
            # only check for attribute function of class/module tm
            if ( isinstance(node.func.value, ast.Name) and
                node.func.value.id == "tm" ):
                # ignore tm.assert_produces_warning(None)/tm.assert_produces_warning()
                if ( len(node.args) == 0 or
                    (isinstance(node.args[0], ast.Constant) and
                    node.args[0].value is None) ):
                    return
                if not any(keyword.arg == "match" for keyword in node.keywords):
                    self.error_set.append((node.lineno, node.col_offset))


# Returns true if a file fails the check
def check_for_match_arg(content: str, filename: str) -> bool:
    tree = ast.parse(content)
    visitor = MatchArgForWarningsChecker()
    visitor.visit(tree)

    if len(visitor.error_set) == 0:
        return False

    for error in visitor.error_set:
        msg = ERROR_MESSAGE.format(
                lineno=error[0],
                col_offset=error[1],
                path=filename,
            )
        sys.stdout.write(msg)

    return True


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")

    args = parser.parse_args(argv)
    isMatchMissing = False

    for filename in args.paths:
        with open(filename, encoding="utf-8") as fd:
            content = fd.read()
            isMatchMissing = check_for_match_arg(content, filename) | isMatchMissing

    if isMatchMissing:
        sys.exit(1)


if __name__ == "__main__":
    main()
