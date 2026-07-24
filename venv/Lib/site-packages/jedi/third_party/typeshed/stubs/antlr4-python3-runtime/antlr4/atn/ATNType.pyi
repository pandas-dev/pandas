from enum import IntEnum

class ATNType(IntEnum):
    LEXER = 0
    PARSER = 1
    @classmethod
    def fromOrdinal(cls, i: int): ...
