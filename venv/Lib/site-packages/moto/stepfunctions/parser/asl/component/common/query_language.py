from __future__ import annotations

import enum

from moto.stepfunctions.parser.asl.antlr.runtime.ASLLexer import ASLLexer
from moto.stepfunctions.parser.asl.component.component import Component


class QueryLanguageMode(enum.Enum):
    JSONPath = ASLLexer.JSONPATH
    JSONata = ASLLexer.JSONATA

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"QueryLanguageMode.{self}({self.value})"


DEFAULT_QUERY_LANGUAGE_MODE: QueryLanguageMode = QueryLanguageMode.JSONPath


class QueryLanguage(Component):
    query_language_mode: QueryLanguageMode

    def __init__(
        self, query_language_mode: QueryLanguageMode = DEFAULT_QUERY_LANGUAGE_MODE
    ):
        self.query_language_mode = query_language_mode
