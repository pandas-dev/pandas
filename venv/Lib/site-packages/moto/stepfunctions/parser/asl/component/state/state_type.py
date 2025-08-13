from enum import Enum

from moto.stepfunctions.parser.asl.antlr.runtime.ASLLexer import ASLLexer


class StateType(Enum):
    Task = ASLLexer.TASK
    Pass = ASLLexer.PASS
    Choice = ASLLexer.CHOICE
    Fail = ASLLexer.FAIL
    Succeed = ASLLexer.SUCCEED
    Wait = ASLLexer.WAIT
    Map = ASLLexer.MAP
    Parallel = ASLLexer.PARALLEL
