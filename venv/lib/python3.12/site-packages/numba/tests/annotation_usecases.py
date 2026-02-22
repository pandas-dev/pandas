"""
Usecases with Python 3 function annotations.  This is a separate module
in order to avoid syntax errors with Python 2.
"""


class AnnotatedClass:
    """
    A class with annotated methods.
    """

    def __init__(self, v: int):
        self.x = v

    def add(self, v: int) -> int:
        return self.x + v
