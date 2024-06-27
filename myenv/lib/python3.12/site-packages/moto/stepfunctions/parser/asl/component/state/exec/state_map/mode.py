from enum import Enum

from moto.stepfunctions.parser.asl.antlr.runtime.ASLLexer import ASLLexer


class Mode(Enum):
    Inline = ASLLexer.INLINE
    Distributed = ASLLexer.DISTRIBUTED
