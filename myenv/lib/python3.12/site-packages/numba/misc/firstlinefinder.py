"""
This module provides helper functions to find the first line of a function
body.
"""

import ast


class FindDefFirstLine(ast.NodeVisitor):
    """
    Attributes
    ----------
    first_stmt_line : int or None
        This stores the first statement line number if the definition is found.
        Or, ``None`` if the definition is not found.
    """

    def __init__(self, code):
        """
        Parameters
        ----------
        code :
            The function's code object.
        """
        self._co_name = code.co_name
        self._co_firstlineno = code.co_firstlineno
        self.first_stmt_line = None

    def _visit_children(self, node):
        for child in ast.iter_child_nodes(node):
            super().visit(child)

    def visit_FunctionDef(self, node: ast.FunctionDef):
        if node.name == self._co_name:
            # Name of function matches.

            # The `def` line may match co_firstlineno.
            possible_start_lines = set([node.lineno])
            if node.decorator_list:
                # Has decorators.
                # The first decorator line may match co_firstlineno.
                first_decor = node.decorator_list[0]
                possible_start_lines.add(first_decor.lineno)
            # Does the first lineno match?
            if self._co_firstlineno in possible_start_lines:
                # Yes, we found the function.
                # So, use the first statement line as the first line.
                if node.body:
                    first_stmt = node.body[0]
                    if _is_docstring(first_stmt):
                        # Skip docstring
                        first_stmt = node.body[1]
                    self.first_stmt_line = first_stmt.lineno
                    return
                else:
                    # This is probably unreachable.
                    # Function body cannot be bare. It must at least have
                    # A const string for docstring or a `pass`.
                    pass
        self._visit_children(node)


def _is_docstring(node):
    if isinstance(node, ast.Expr):
        if (isinstance(node.value, ast.Constant)
                and isinstance(node.value.value, str)):
            return True
    return False


def get_func_body_first_lineno(pyfunc):
    """
    Look up the first line of function body using the file in
    ``pyfunc.__code__.co_filename``.

    Returns
    -------
    lineno : int; or None
        The first line number of the function body; or ``None`` if the first
        line cannot be determined.
    """
    co = pyfunc.__code__
    try:
        with open(co.co_filename) as fin:
            file_content = fin.read()
    except (FileNotFoundError, OSError):
        return
    else:
        tree = ast.parse(file_content)
        finder = FindDefFirstLine(co)
        finder.visit(tree)
        return finder.first_stmt_line
