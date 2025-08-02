import abc

from moto.stepfunctions.parser.asl.antlr.runtime.ASLParserVisitor import (
    ASLParserVisitor,
)
from moto.stepfunctions.parser.asl.parse.asl_parser import AmazonStateLanguageParser


class StaticAnalyser(ASLParserVisitor, abc.ABC):
    def analyse(self, definition: str) -> None:
        _, parser_rule_context = AmazonStateLanguageParser.parse(definition)
        self.visit(parser_rule_context)
