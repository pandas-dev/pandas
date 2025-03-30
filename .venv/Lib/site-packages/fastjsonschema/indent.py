def indent(func):
    """
    Decorator for allowing to use method as normal method or with
    context manager for auto-indenting code blocks.
    """
    def wrapper(self, line, *args, optimize=True, **kwds):
        last_line = self._indent_last_line
        line = func(self, line, *args, **kwds)
        # When two blocks have the same condition (such as value has to be dict),
        # do the check only once and keep it under one block.
        if optimize and last_line == line:
            self._code.pop()
        self._indent_last_line = line
        return Indent(self, line)
    return wrapper


class Indent:
    def __init__(self, instance, line):
        self.instance = instance
        self.line = line

    def __enter__(self):
        self.instance._indent += 1

    def __exit__(self, type_, value, traceback):
        self.instance._indent -= 1
        self.instance._indent_last_line = self.line
