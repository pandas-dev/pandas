"""
Convert dataframe_sql statement to run on pandas dataframes
"""
import os
from pathlib import Path

from lark import Lark, UnexpectedToken

from pandas.core.frame import DataFrame

from pandas.io.dataframe_sql.parsers import SQLTransformer
from pandas.io.dataframe_sql.sql_exception import InvalidQueryException
from pandas.io.dataframe_sql.sql_objects import AmbiguousColumn

SHOW_TREE = False
SHOW_DF = False


with open(file=os.path.join(Path(__file__).parent, "sql.grammar")) as sql_grammar_file:
    GRAMMAR_TEXT = sql_grammar_file.read()


def register_temp_table(frame: DataFrame, table_name: str):
    table_info = TableInfo()
    table_info.register_temporary_table(frame, table_name)


def remove_temp_table(table_name: str):
    table_info = TableInfo()
    table_info.remove_temp_table(table_name)


def query(sql: str):
    return SqlToPandas(sql).data_frame


class SqlToPandas:
    """
    Class that handles conversion from dataframe_sql to pandas data frame methods.
    """

    def __init__(self, sql: str):
        self.sql = sql
        table_info = TableInfo()
        if SHOW_TREE:
            self.parser = Lark(GRAMMAR_TEXT, parser="lalr")
        else:
            self.parser = Lark(
                GRAMMAR_TEXT,
                parser="lalr",
                transformer=SQLTransformer(
                    table_info.dataframe_name_map,
                    table_info.dataframe_map,
                    table_info.column_name_map,
                    table_info.column_to_dataframe_name,
                ),
            )
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
        Splits the dataframe_sql into tokens
        :return:
        """
        try:
            return self.parser.parse(self.sql)
        except UnexpectedToken as err:
            message = (
                f"Expected one of the following input(s): {err.expected}\n"
                f"Unexpected input at line {err.line}, column {err.column}\n"
                f"{err.get_context(self.sql)}"
            )
            raise InvalidQueryException(message)


class TableInfo:
    column_to_dataframe_name = {}
    column_name_map = {}
    dataframe_name_map = {}
    dataframe_map = {}

    def add_column_to_column_to_dataframe_name_map(self, column, table):
        """
        Adds a column to the column_to_dataframe_name_map
        :param column:
        :param table:
        :return:
        """
        if self.column_to_dataframe_name.get(column) is None:
            self.column_to_dataframe_name[column] = table
        elif isinstance(self.column_to_dataframe_name[column], AmbiguousColumn):
            self.column_to_dataframe_name[column].tables.append(table)
        else:
            original_table = self.column_to_dataframe_name[column]
            self.column_to_dataframe_name[column] = AmbiguousColumn(
                [original_table, table]
            )

    def register_temporary_table(self, frame: DataFrame, table_name: str):
        """
        Registers dataframe info so that the data frame is prepared to be used
        with dataframe_sql
        :param frame:
        :param table_name:
        :return:
        """
        if table_name.lower() in self.dataframe_name_map:
            raise Exception(
                f"A table {table_name.lower()} has already been registered. Keep in "
                f"mind that table names are case insensitive"
            )

        self.dataframe_name_map[table_name.lower()] = table_name
        self.dataframe_map[table_name] = frame
        self.column_name_map[table_name] = {}
        for column in frame.columns:
            lower_column = column.lower()
            self.column_name_map[table_name][lower_column] = column
            self.add_column_to_column_to_dataframe_name_map(lower_column, table_name)

    def remove_temp_table(self, table_name: str):
        """
        Unregisters a dataframe from table info
        :param table_name:
        :return:
        """
        if table_name.lower() not in self.dataframe_name_map:
            raise Exception(f"Table {table_name.lower()} is not registered")
        real_table_name = self.dataframe_name_map[table_name.lower()]

        columns = self.dataframe_map[real_table_name].columns.to_list()
        for column in columns:
            lower_column = column.lower()
            value = self.column_to_dataframe_name[lower_column]
            if isinstance(value, AmbiguousColumn):
                value.tables.remove(real_table_name)
                if len(value.tables) == 1:
                    self.column_to_dataframe_name[lower_column] = value.tables[0]
            else:
                del self.column_to_dataframe_name[lower_column]

        del self.dataframe_name_map[table_name.lower()]
        del self.dataframe_map[real_table_name]
        del self.column_name_map[real_table_name]
