from enum import Enum

from moto.stepfunctions.parser.asl.antlr.runtime.ASLLexer import ASLLexer


class ExecutionType(Enum):
    Standard = ASLLexer.STANDARD
