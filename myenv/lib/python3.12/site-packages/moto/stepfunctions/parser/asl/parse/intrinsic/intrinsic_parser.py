import abc

from antlr4 import CommonTokenStream, InputStream

from moto.stepfunctions.parser.asl.antlr.runtime.ASLIntrinsicLexer import (
    ASLIntrinsicLexer,
)
from moto.stepfunctions.parser.asl.antlr.runtime.ASLIntrinsicParser import (
    ASLIntrinsicParser,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.function import Function
from moto.stepfunctions.parser.asl.parse.intrinsic.preprocessor import Preprocessor


class IntrinsicParser(abc.ABC):
    @staticmethod
    def parse(src: str) -> Function:
        input_stream = InputStream(src)
        lexer = ASLIntrinsicLexer(input_stream)
        stream = CommonTokenStream(lexer)
        parser = ASLIntrinsicParser(stream)
        tree = parser.func_decl()
        preprocessor = Preprocessor()
        function: Function = preprocessor.visit(tree)
        return function
