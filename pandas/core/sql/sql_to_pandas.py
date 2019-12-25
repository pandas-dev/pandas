"""
Convert sql statement to run on pandas dataframes
"""
from pathlib import Path
import os
from lark import Lark, UnexpectedToken
from pandas.core.sql.parsers import SQLTransformer
from pandas.core.sql.sql_exception import InvalidQueryException

SHOW_TREE = False
SHOW_DF = False

with open(file=os.path.join(Path(__file__).parent, "sql.grammar")) as sql_grammar_file:
    GRAMMAR_TEXT = sql_grammar_file.read()


class SqlToPandas:
    """
    Class that handles conversion from sql to pandas data frame methods.
    """
    def __init__(self, sql: str, all_global_vars):
        self.sql = sql
        if SHOW_TREE:
            self.parser = Lark(GRAMMAR_TEXT, parser='lalr')
        else:
            self.parser = Lark(GRAMMAR_TEXT, parser='lalr', transformer=SQLTransformer(all_global_vars))
        self.ast = self.parse_sql()
        if SHOW_TREE or SHOW_DF:
            print("Result:")
        if SHOW_TREE:
            print(self.ast)
            print(self.ast.pretty())
        if SHOW_DF:
            print(self.ast)
        self.data_frame = self.ast

    def parse_sql(self):
        """
        Splits the sql into tokens
        :return:
        """
        try:
            return self.parser.parse(self.sql)
        except UnexpectedToken as err:
            message = f"Expected one of the following input(s): {err.expected}\n" \
                      f"Unexpected input at line {err.line}, column {err.column}\n" \
                      f"{err.get_context(self.sql)}"
            raise InvalidQueryException(message)
