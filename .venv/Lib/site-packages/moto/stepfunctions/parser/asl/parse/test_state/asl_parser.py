from typing import Tuple

from antlr4 import CommonTokenStream, InputStream, ParserRuleContext

from moto.stepfunctions.parser.asl.antlr.runtime.ASLLexer import ASLLexer
from moto.stepfunctions.parser.asl.antlr.runtime.ASLParser import ASLParser
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.parse.asl_parser import (
    AmazonStateLanguageParser,
    ASLParserException,
    SyntaxErrorListener,
)
from moto.stepfunctions.parser.asl.parse.test_state.preprocessor import (
    TestStatePreprocessor,
)


class TestStateAmazonStateLanguageParser(AmazonStateLanguageParser):
    @staticmethod
    def parse(definition: str) -> Tuple[EvalComponent, ParserRuleContext]:
        # Attempt to build the AST and look out for syntax errors.
        syntax_error_listener = SyntaxErrorListener()

        input_stream = InputStream(definition)
        lexer = ASLLexer(input_stream)
        stream = CommonTokenStream(lexer)
        parser = ASLParser(stream)
        parser.removeErrorListeners()
        parser.addErrorListener(syntax_error_listener)
        # Unlike the main Program parser, TestState parsing occurs at a state declaration level.
        tree = parser.state_decl_body()

        errors = syntax_error_listener.errors
        if errors:
            raise ASLParserException(errors=errors)

        # Attempt to preprocess the AST into evaluation components.
        preprocessor = TestStatePreprocessor()
        test_state_program = preprocessor.visit(tree)

        return test_state_program, tree
