"""
Module containing all lark transformer classes
"""
import re
from datetime import datetime, date
from typing import Tuple, List
from lark import Transformer, v_args
from lark.lexer import Token
from lark.tree import Tree
from pandas.core.frame import DataFrame
from pandas.core.reshape.merge import merge
from pandas.core.reshape.concat import concat
from pandas.io.dataframe_sql.sql_objects import (
    AmbiguousColumn,
    Column,
    Subquery,
    Literal,
    Number,
    String,
    Date,
    Bool,
    Expression,
    Value,
)
from pandas.io.dataframe_sql.sql_exception import DataFrameDoesNotExist

DEBUG = False
PRINT = False

ORDER_TYPES = ["asc", "desc", "ascending", "descending"]
ORDER_TYPES_MAPPING = {
    "asc": "asc",
    "desc": "desc",
    "ascending": "asc",
    "descending": "desc",
}
GET_TABLE_REGEX = re.compile(
    r"^(?P<table>[a-z_]\w*)\.(?P<column>[a-z_]\w*)$", re.IGNORECASE
)
FUNCTION_MAPPING = {
    "average": "mean",
    "avg": "mean",
    "mean": "mean",
    "maximum": "max",
    "max": "max",
    "minimum": "min",
    "min": "min",
}
PANDAS_TYPE_PYTHON_TYPE_FUNCTION = {
    "object": str,
    "int64": int,
    "float64": float,
    "bool": bool,
}
PANDAS_TYPE_TO_SQL_TYPE = {
    "object": String,
    "int64": Number,
    "float64": Number,
    "bool": Bool,
    "datetime64": Date,
}


def num_eval(arg):
    """
    Takes an argument that may be a string or number and outputs a number
    :param arg:
    :return:
    """
    assert isinstance(arg, (Token, float, int))
    if isinstance(arg, str):
        # pylint: disable=eval-used
        return eval(arg)
    return arg


def get_wrapper_value(value):
    """
    If the value is a literal return it's value
    :param value:
    :return:
    """
    if isinstance(value, Value):
        return value.value
    return value


class TransformerBaseClass(Transformer):
    """
    Base class for transformers
    """

    def __init__(
        self,
        dataframe_name_map=None,
        dataframe_map=None,
        column_name_map=None,
        column_to_dataframe_name=None,
        _temp_dataframes_dict=None,
    ):
        Transformer.__init__(self, visit_tokens=False)
        self.dataframe_name_map = dataframe_name_map
        self.dataframe_map = dataframe_map
        self.column_name_map = column_name_map
        self.column_to_dataframe_name = column_to_dataframe_name
        self._temp_dataframes_dict = _temp_dataframes_dict

    def get_frame(self, frame_name) -> DataFrame:
        """
        Returns the dataframe with the name given
        :param frame_name:
        :return:
        """
        if isinstance(frame_name, Token):
            frame_name = frame_name.value
        if isinstance(frame_name, Subquery):
            frame_name = frame_name.name
        return self.dataframe_map[frame_name]

    def set_column_value(self, column: Column):
        """
        Sets the column value based on what it is in the dataframe
        :param column:
        :return:
        """
        if column.name != "*":
            dataframe_name = self.column_to_dataframe_name[column.name.lower()]
            if isinstance(dataframe_name, AmbiguousColumn):
                raise Exception(f"Ambiguous column reference: {column.name}")
            dataframe = self.get_frame(dataframe_name)
            column_true_name = self.column_name_map[dataframe_name][column.name.lower()]
            column.value = dataframe[column_true_name]
            column.table = dataframe_name

    def column_name(self, name_list_format: List[str]):
        """
        Returns a column token with the name extracted
        :param name_list_format: List formatted name
        :return: Tree with column token
        """
        name = "".join(name_list_format)
        column = Column(name="".join(name))
        self.set_column_value(column)
        return column


# pylint: disable=no-self-use, too-many-public-methods, too-many-instance-attributes
class InternalTransformer(TransformerBaseClass):
    """
    Evaluates subtrees with knowledge of provided tables that are in the proper scope
    """

    def __init__(
        self, tables, dataframe_map, column_name_map, column_to_dataframe_name
    ):
        TransformerBaseClass.__init__(
            self, dataframe_map=dataframe_map, column_name_map=column_name_map
        )
        self.tables = [
            table.name if isinstance(table, Subquery) else table for table in tables
        ]
        self.column_to_dataframe_name = {}
        for column in column_to_dataframe_name:
            table = column_to_dataframe_name.get(column)
            if isinstance(table, AmbiguousColumn):
                table_name = self.tables[0]
                if table_name in table.tables:
                    self.column_to_dataframe_name[column] = table_name
            if table in self.tables:
                self.column_to_dataframe_name[column] = table

        # These variables need instance scope for rank instance method
        self.partition_func_dict = {
            self.set_rank_regular: self.set_rank_regular_partition,
            self.set_rank_dense: self.set_rank_dense_partition,
        }
        self.partition_rank_counter = {}
        self.partition_rank_offset = {}
        self.rank_counter = 1
        self.rank_offset = 0
        self.rank_map = {}
        self.last_key = None

    def mul(self, args):
        """
        Returns the product two numbers
        """
        arg1 = args[0]
        arg2 = args[1]
        return num_eval(arg1) * num_eval(arg2)

    def expression_mul(self, args):
        """
        Returns the product of two expressions
        :param args:
        :return:
        """
        arg1 = args[0]
        arg2 = args[1]
        return arg1 * arg2

    def add(self, args):
        """
        Returns the sum two numbers
        """
        arg1 = args[0]
        arg2 = args[1]
        return num_eval(arg1) + num_eval(arg2)

    def expression_add(self, args):
        """
        Returns the sum of two expressions
        :param args:
        :return:
        """
        arg1 = args[0]
        arg2 = args[1]
        return arg1 + arg2

    def sub(self, args):
        """
        Returns the difference between two numbers
        """
        arg1 = args[0]
        arg2 = args[1]
        return num_eval(arg1) - num_eval(arg2)

    def expression_sub(self, args):
        """
        Returns the difference between two expressions
        :param args:
        :return:
        """
        arg1 = args[0]
        arg2 = args[1]
        return arg1 - arg2

    def div(self, args):
        """
        Returns the division of two numbers
        """
        arg1 = args[0]
        arg2 = args[1]
        return num_eval(arg1) / num_eval(arg2)

    def expression_div(self, args):
        """
        Returns the difference between two expressions
        :param args:
        :return:
        """
        arg1 = args[0]
        arg2 = args[1]
        return arg1 / arg2

    def number(self, numerical_value):
        """
        Return a number token with a numeric value as a child
        :param numerical_value:
        :return:
        """
        return Number(num_eval(numerical_value[0]))

    def string(self, string_token):
        """
        Return value of the token associated with the string
        :param string_token:
        :return:
        """
        return String(string_token[0].value)

    def timestamp_expression(self, date_list):
        """
        Return a timestamp object
        :param date_list:
        :return:
        """
        return date_list[0]

    @staticmethod
    def int_token_list(token_list):
        """
        Returns a list of integer from a list of tokens
        :param self:
        :return:
        """
        return [int(token.value) for token in token_list]

    def date(self, date_list):
        """
        Returns list with correct date integers
        :param date_list:
        :return:
        """
        return self.int_token_list(date_list)

    def time(self, time_list):
        """
        Returns list with correct time integers
        :param time_list:
        :return:
        """
        return self.int_token_list(time_list)

    def custom_timestamp(self, datetime_list):
        """
        Return a custom time stamp based on user input
        :param datetime_list:
        :return:
        """
        return Date(datetime(*(datetime_list[0] + datetime_list[1])))

    def datetime_now(self, *extra_args):
        """
        Return current date and time
        :param extra_args: Arguments that lark parser must pass in
        :return:
        """
        date_value = Date(datetime.now())
        date_value.set_alias("now()")
        return date_value

    def date_today(self, *extra_args):
        """
        Return current date
        :param extra_args: Arguments that lark parser must pass in
        :return:
        """
        date_value = Date(date.today())
        date_value.set_alias("today()")
        return date_value

    def equals(self, expressions):
        """
        Compares two expressions for equality
        :param expressions:
        :return:
        """
        return expressions[0] == expressions[1]

    def greater_than(self, expressions):
        """
        Performs a greater than expression
        :param expressions:
        :return:
        """
        return expressions[0] > expressions[1]

    def less_than(self, expressions):
        """
        Performs a less than expression
        :param expressions:
        :return:
        """
        return expressions[0] < expressions[1]

    def between(self, expressions):
        """
        Performs a less than or equal and greater than or equal
        :param expressions:
        :return:
        """
        return (expressions[0] >= expressions[1]) & (expressions[0] <= expressions[2])

    def in_expr(self, expressions):
        """
        Evaluate in expression
        :param expressions:
        :return:
        """
        in_list = [
            expression.value if isinstance(expression, Literal) else expression
            for expression in expressions[1:]
        ]
        return expressions[0].value.isin(in_list)

    def not_in_expr(self, expressions):
        """
        Negate in expr
        :param expressions:
        :return:
        """
        return ~self.in_expr(expressions)

    def bool_expression(self, bool_expression):
        """
        Return the bool expression
        :param bool_expression:
        :return: boolean expression
        """
        return bool_expression[0]

    def negated_bool_expression(self, bool_expression):
        """
        Returns a negated boolean expression
        :param bool_expression:
        :return:
        """
        print("Negated bool", bool_expression)

    def where_expr(self, truth_value_dataframe):
        """
        Return a where token
        :param truth_value_dataframe:
        :return: Token
        """
        return Token("where_expr", truth_value_dataframe[0])

    def function_name(self, function_name):
        """
        Returns the function name tree
        :param function_name:
        :return:
        """
        function_name = function_name[0].lower()
        true_function_name = FUNCTION_MAPPING.get(function_name)
        if true_function_name:
            return Tree("aggregate", true_function_name)
        return Tree("function", function_name)

    def alias_string(self, name: str):
        """
        Returns an alias token with the name extracted
        :param name:
        :return:
        """
        return Tree("alias", str(name[0]))

    def from_expression(self, expression):
        """
        Return a from expression token
        :param expression:
        :return: Token from expression
        """
        expression = expression[0]
        if isinstance(expression, Subquery):
            value = expression
        else:
            value = expression.value
        return Token("from_expression", value)

    def when_then(self, when_then):
        """
        When / then expression
        :param when_then:
        :return:
        """
        then_value = get_wrapper_value(when_then[1])
        return when_then[0], then_value

    def case_expression(self, when_expressions):
        """
        Handles dataframe_sql case expressions
        :param when_expressions:
        :return:
        """
        # TODO Possibly a problem when dealing with booleans
        new_column = when_expressions[0][0]
        for when_expression in when_expressions:
            if isinstance(when_expression, Tuple):
                new_column[when_expression[0]] = when_expression[1]
            else:
                new_column[new_column == False] = get_wrapper_value(when_expression) # noqa
        return Expression(value=new_column)

    def rank_form(self, form):
        """
        Returns the rank form
        :param form:
        :return:
        """
        return form

    def order_asc(self, column):
        """
        Return expression in asc order
        :param expression:
        :return:
        """
        column = column[0]
        column.value = column.value.copy().sort_values(ascending=True)
        return Token("order", (column, True))

    def order_desc(self, column):
        """
        Return expression in asc order
        :param expression:
        :return:
        """
        column = column[0]
        column.value = column.value.copy().sort_values(ascending=False)
        return Token("order", (column, False))

    def set_rank_dense(self, row):
        """
        Set rank dense rank without gaps between consecutive ranks
        :param row: A pandas row object
        :return: The dense rank of the current row
        """
        key = str(list(row))
        if self.rank_map.get(key):
            return self.rank_map[key]
        self.rank_map[key] = self.rank_counter
        rank = self.rank_counter
        self.rank_counter += 1
        return rank

    def set_rank_regular(self, row):
        """
        Set rank in traditional SQL database fashion
        :param row: A pandas row object
        :return: The rank of the current row
        """
        key = str(list(row))
        if self.rank_map.get(key):
            self.rank_offset += 1
            return self.rank_map[key]
        rank = self.rank_counter + self.rank_offset
        self.rank_map[key] = rank
        self.rank_counter += 1
        return rank

    def set_rank_regular_partition(self, row, partition_slice):
        """
        Set rank with a partition
        :param row: A pandas row object
        :param partition_slice: Integer to slice the values in the row based on
        partition
        :return: The partitioned rank of the current row
        """
        row_list = list(row)[:partition_slice]
        partition_list = list(row)[partition_slice:]
        key = str(row_list)
        partition_key = str(partition_list)
        if self.rank_map.get(partition_key):
            if self.rank_map[partition_key].get(key):
                self.partition_rank_offset[partition_key] += 1
                return self.rank_map[partition_key][key]
            self.partition_rank_counter[partition_key] += 1
            rank = (
                self.partition_rank_counter[partition_key]
                + self.partition_rank_offset[partition_key]
            )
            self.rank_map[partition_key][key] = rank
        else:
            rank = 1
            self.rank_map[partition_key] = {}
            self.partition_rank_counter[partition_key] = rank
            self.partition_rank_offset[partition_key] = 0
            self.rank_map[partition_key][key] = rank
        return rank

    def set_rank_dense_partition(self, row, partition_slice):
        """
        Set rank with a partition
        :param row: A pandas row object
        :param partition_slice: Integer to slice the values in the row based on
        partition
        :return: The partitioned rank of the current row
        """
        row_list = list(row)[:partition_slice]
        partition_list = list(row)[partition_slice:]
        key = str(row_list)
        partition_key = str(partition_list)
        if self.rank_map.get(partition_key):
            if self.rank_map[partition_key].get(key):
                return self.rank_map[partition_key][key]
            self.partition_rank_counter[partition_key] += 1
            rank = self.partition_rank_counter[partition_key]
            self.rank_map[partition_key][key] = rank
        else:
            rank = 1
            self.rank_map[partition_key] = {}
            self.partition_rank_counter[partition_key] = 1
            self.rank_map[partition_key][key] = rank
        return rank

    def partition_by(self, column_list):
        """
        Returns a partition token containing the corresponding column
        :param column_list: List containing only one column
        :return:
        """
        column = column_list[0]
        return Token("partition", column)

    def rank(self, tokens, rank_function):
        """
        Returns the evaluated rank expressions
        :param tokens: Tokens making up the rank expression
        :param rank_function: Function to be used in rank evaluation
        :return:
        """
        expressions = tokens[0]
        series_list = []
        order_list = []
        partition_list = []
        for token in expressions:
            if token.type == "order":
                column = token.value[0]
                ascending = token.value[1]
                series_list.append(column.value)
                order_list.append(ascending)
            elif token.type == "partition":
                column = token.value
                partition_list.append(column.value)
        rank_df = concat(series_list, axis=1)
        column_names = rank_df.columns.to_list()
        if partition_list:
            rank_df = concat(series_list + partition_list, axis=1)
            rank_function = self.partition_func_dict[rank_function]
            column_slice = len(column_names)
        rank_df.sort_values(by=column_names, ascending=order_list, inplace=True)
        if partition_list:
            rank = rank_df.apply(rank_function, args=(column_slice,), axis=1)
            self.partition_rank_counter = {}
            self.partition_rank_offset = {}
        else:
            rank = rank_df.apply(rank_function, axis=1)
            self.rank_counter = 1
            self.rank_offset = 0
        self.rank_map = {}
        rank_df["rank"] = rank
        return Expression(value=rank)

    def rank_expression(self, tokens):
        """
        Handles rank expressions
        :param tokens:
        :return:
        """
        return self.rank(tokens, self.set_rank_regular)

    def dense_rank_expression(self, tokens):
        """
        Handles dense_rank_expressions
        :param tokens:
        :return:
        """
        return self.rank(tokens, self.set_rank_dense)

    def select_expression(self, expression_and_alias):
        """
        Returns the appropriate object for the given expression
        :param expression_and_alias: An expression token and A token containing the
        name to be assigned
        :return:
        """
        expression = expression_and_alias[0]
        alias = ""
        if len(expression_and_alias) == 2:
            alias = expression_and_alias[1]
        if isinstance(expression, Tree):
            value = expression.children
            if expression.data == "sql_function":
                function = value[0].children
                value = value[1]
                expression = Expression(value=value, function=function)

        if alias:
            expression.set_alias(alias.children)
        return expression

    def join(self, *args):
        """
        Extracts the join expression
        :param args: Arguments that are passed to the join
        :return: join expression
        """
        return args[0]

    def group_by(self, column):
        """
        Returns a group token
        :param column: Column to group by
        :return: group token
        """
        column = column[0]
        return Token("group", str(column.name))

    def as_type(self, column_and_type):
        """
        Extracts token type and returns tree object with expression and type
        :param expression: Expression to be evaluated / the name of a column
        :param typename: Data type
        :return:
        """
        column = column_and_type[0]
        typename = column_and_type[1]
        column.typename = typename.value
        return column

    def literal_cast(self, value_and_type: list):
        """
        Cast variable as the given pandas_type for a literal
        :param value_and_type: Value and pandas dtype to be cast as
        :return:
        """
        value_wrapper = value_and_type[0]
        pandas_type = value_and_type[1]
        if pandas_type == "datetime64":
            date_value = datetime.strptime(value_wrapper.value, "%Y-%m-%d")
            return Date(date_value)
        conversion_func = PANDAS_TYPE_PYTHON_TYPE_FUNCTION[pandas_type]
        new_type = PANDAS_TYPE_TO_SQL_TYPE[pandas_type]
        new_value = new_type(conversion_func(value_wrapper.value))
        return new_value


# pylint: disable=no-self-use
class HavingTransformer(TransformerBaseClass):
    """
    Transformer for having clauses since group by needs to be applied first
    """

    # pylint: disable=too-many-arguments
    def __init__(
        self, tables, group_by, dataframe_map, column_name_map, column_to_dataframe_name
    ):
        self.tables = tables
        self.group_by = group_by
        TransformerBaseClass.__init__(
            self,
            dataframe_map=dataframe_map,
            column_name_map=column_name_map,
            column_to_dataframe_name=column_to_dataframe_name,
        )

    def aggregate(self, function_name_list_form):
        """
        Return the string representation fo aggregate function name instead of list
        :param function_name_list_form:
        :return:
        """
        return "".join(function_name_list_form)

    def function_name(self, tokens):
        """
        Extracts function name from token
        :param fucntion_name:
        :return:
        """
        return tokens[0].value

    def sql_function(self, function_expr):
        """
        Handles presence of functions in an expression
        :param function_expr: Function expression
        :return:
        """
        aggregate_name = function_expr[0]
        column = function_expr[1]
        table = self.dataframe_map[column.table]
        aggregates = {column.name: aggregate_name}
        if self.group_by:
            new_series = (
                table.groupby(self.group_by).aggregate(aggregates).reset_index()
            )
        else:
            new_series = table.aggregate(aggregates).to_frame().transpose()
        return new_series[column.name]

    def having_expr(self, having_expr):
        """
        Handles having expressions
        :param having_expr:
        :return:
        """
        internal_transformer = InternalTransformer(
            self.tables,
            self.dataframe_map,
            self.column_name_map,
            self.column_to_dataframe_name,
        )
        having_expr = Tree("having_expr", having_expr)
        return internal_transformer.transform(having_expr)


# pylint: disable=no-self-use, super-init-not-called
@v_args(inline=True)
class SQLTransformer(TransformerBaseClass):
    """
    Transformer for the lark dataframe_sql parser
    """

    def __init__(
        self,
        dataframe_name_map={},
        dataframe_map={},
        column_name_map={},
        column_to_dataframe_name={},
    ):
        TransformerBaseClass.__init__(
            self,
            dataframe_name_map,
            dataframe_map,
            column_name_map,
            column_to_dataframe_name,
            _temp_dataframes_dict={},
        )

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

    def table(self, table_name, alias=""):
        """
        Check for existence of pandas dataframe with same name
        If not exists raise DataFrameDoesNotExist
        Otherwise return the name of the actual DataFrame
        :return:
        """
        table_name = table_name.lower()
        if table_name not in self.dataframe_name_map:
            raise DataFrameDoesNotExist(table_name)
        if alias:
            self.dataframe_name_map[alias] = self.dataframe_name_map[table_name]
        return Token("table", self.dataframe_name_map[table_name])

    def order_by_expression(self, rank_tree):
        """
        Returns the column name for the order expression
        :param rank_tree: Tree containing order info
        :return:
        """
        order_type = rank_tree.data
        ascending = order_type == "order_asc"
        return Token("order_by", (rank_tree.children[0].children, ascending))

    def integer(self, integer_token):
        """
        Returns the integer value
        :param integer_token:
        :return:
        """
        integer_value = int(integer_token.value)
        return integer_value

    def limit_count(self, limit_count_value):
        """
        Returns a limit token
        :param limit_count_value:
        :return:
        """
        return Token("limit", limit_count_value)

    def query_expr(self, query_info, *args):
        """
        Handles the full query, including order and set operations such as union
        :param query_info: Map of all query information
        :param args: Additional arguments aside from query info
        :return: Query info
        """
        order_by = []
        limit = None
        for token in args:
            if isinstance(token, Token):
                if token.type == "order_by":
                    order_by.append(token.value)
                elif token.type == "limit":
                    limit = token.value
        query_info["order_by"] = order_by
        query_info["limit"] = limit
        return query_info

    def subquery(self, query_info, alias):
        """
        Handle subqueries amd return a subquery object
        :param query_info:
        :param alias:
        :return:
        """
        alias_name = alias.children[0].value
        self.dataframe_map[alias_name] = self.to_dataframe(query_info)
        subquery = Subquery(name=alias_name, query_info=query_info)
        self.column_name_map[alias_name] = {}
        for column in self.dataframe_map[alias_name].columns:
            self.add_column_to_column_to_dataframe_name_map(column.lower(), alias_name)
            self.column_name_map[alias_name][column.lower()] = column
        return subquery

    # pylint: disable=arguments-differ
    def column_name(self, *names):
        full_name = ".".join([str(name) for name in names])
        return Tree("column_name", full_name)

    def join(self, join_expression):
        """
        Handle join tree
        :param join_expression:
        :return:
        """
        return join_expression

    def get_lower_columns(self, table_name):
        """
        Returns a list of lower case column names for a given table name
        :param column_list:
        :return:
        """
        return [column.lower() for column in list(self.get_frame(table_name).columns)]

    def determine_column_side(self, column, left_table, right_table):
        """
        Check if column table prefix is one of the two tables (if there is one) AND
        the column has to be in one of the two tables
        """
        column_match = GET_TABLE_REGEX.match(column)
        column_table = ""
        if column_match:
            column = column_match.group("column").lower()
            column_table = column_match.group("table").lower()

        left_columns = self.column_name_map[left_table]
        right_columns = self.column_name_map[right_table]
        if column not in left_columns and column not in right_columns:
            raise Exception("Column not found")

        left_table = left_table.lower()
        right_table = right_table.lower()
        if column_table:
            if column_table == left_table and column in left_columns:
                return "left", column
            if column_table == right_table and column in right_columns:
                return "right", column
            raise Exception("Table specified in join columns not present in join")
        if column in left_columns and column in right_columns:
            raise Exception(
                f"Ambiguous column: {column}\nSpecify table name with table_name"
                f".{column}"
            )
        if column in left_columns:
            return "left", column
        if column in right_columns:
            return "right", column
        raise Exception("Column does not exist in either table")

    def join_expression(self, *args):
        """
        Evaluate a join into one dataframe using a merge method
        :return:
        """
        # There will only ever be four args if a join is specified and three if a
        # join isn't specified
        if len(args) == 3:
            join_type = "inner"
            table1 = args[0]
            table2 = args[1]
            join_condition = args[2]
        else:
            table1 = args[0]
            join_type = args[1]
            table2 = args[2]
            join_condition = args[3]
            if "outer" in join_type:
                match = re.match(r"(?P<type>.*)\souter", join_type)
                join_type = match.group("type")
            if join_type in ("full", "cross"):
                join_type = "outer"

        # Check that there is a column from both sides
        column_comparison = join_condition.children[0].children
        column1 = str(column_comparison[0].children)
        column2 = str(column_comparison[1].children)

        column1_side, column1 = self.determine_column_side(column1, table1, table2)
        column2_side, column2 = self.determine_column_side(column2, table1, table2)
        if column1_side == column2_side:
            raise Exception("Join columns must be one column from each join table!")
        column1 = self.column_name_map[table1][column1]
        column2 = self.column_name_map[table2][column2]
        if column1_side == "left":
            left_on = column1
            right_on = column2
        else:
            left_on = column2
            right_on = column1

        dictionary_name = f"{table1}x{table2}"
        self.dataframe_map[dictionary_name] = self.get_frame(table1).merge(
            right=self.get_frame(table2),
            how=join_type,
            left_on=left_on,
            right_on=right_on,
        )
        return Subquery(dictionary_name, query_info="")

    @staticmethod
    def has_star(column_list: List[str]):
        """
        Returns true if any columns have a star
        :param column_list:
        :return:
        """
        for column_name in column_list:
            if re.match(r"\*", column_name):
                return True
        return False

    @staticmethod
    def handle_non_token_non_tree(query_info: dict, token, token_pos):
        """
        Handles non token non tree items and extracts necessary query information
        from it
        :param query_info: Dictionary of all info about the query
        :param token: Item being handled
        :param token_pos: Ordinal position of the item
        :return:
        """
        query_info["all_names"].append(token.final_name)
        query_info["name_order"][token.final_name] = token_pos

        if token.typename:
            query_info["conversions"][token.final_name] = token.typename

        if isinstance(token, Column):
            query_info["columns"].append(token)
            if token.alias:
                query_info["aliases"][token.name] = token.alias

        if isinstance(token, Expression):
            query_info["expressions"].append(token)
            if token.alias:
                query_info["aliases"][str(token.value)] = token.alias
            if token.function:
                query_info["aggregates"][token.alias] = token.function

        if isinstance(token, Literal):
            query_info["literals"].append(token)

    def handle_token(self, query_info: dict, token, token_pos):
        """
        Handles token and extracts necessary query information from it
        :param query_info: Dictionary of all info about the query
        :param token: Item being handled
        :param token_pos: Ordinal position of the token
        :return:
        """
        if isinstance(token, Token):
            if token.type == "from_expression":
                query_info["frame_names"].append(token.value)
            elif token.type == "group":
                query_info["group_columns"].append(token.value)
            elif token.type == "where_expr":
                query_info["where_expr"] = token.value
        elif isinstance(token, Tree):
            if token.data == "having_expr":
                query_info["having_expr"] = token
        else:
            self.handle_non_token_non_tree(query_info, token, token_pos)

    def select(self, *select_expressions: Tuple[Tree]):
        """
        Forms the final sequence of methods that will be executed
        :param select_expressions:
        :return:
        """
        if DEBUG:
            print("Select Expressions:", select_expressions)

        tables = []
        query_info = {
            "columns": [],
            "expressions": [],
            "literals": [],
            "frame_names": [],
            "aliases": {},
            "all_names": [],
            "name_order": {},
            "conversions": {},
            "aggregates": {},
            "group_columns": [],
            "where_expr": None,
            "distinct": False,
            "having_expr": None,
        }
        for select_expression in select_expressions:
            if isinstance(select_expression, Tree):
                if select_expression.data == "from_expression":
                    tables.append(select_expression.children[0])
                elif select_expression.data == "having_expr":
                    query_info["having_expr"] = select_expression

        select_expressions_no_having = tuple(
            select_expression
            for select_expression in select_expressions
            if isinstance(select_expression, Tree)
            and select_expression.data != "having_expr"
            or not isinstance(select_expression, Tree)
        )

        internal_transformer = InternalTransformer(
            tables,
            self.dataframe_map,
            self.column_name_map,
            self.column_to_dataframe_name,
        )
        select_expressions = internal_transformer.transform(
            Tree("select", select_expressions_no_having)
        ).children

        if isinstance(select_expressions[0], Token):
            if str(select_expressions[0]) == "distinct":
                query_info["distinct"] = True
            select_expressions = select_expressions[1:]

        for token_pos, token in enumerate(select_expressions):
            self.handle_token(query_info, token, token_pos)

        having_expr = query_info["having_expr"]
        if having_expr is not None:
            having_expr = (
                HavingTransformer(
                    tables,
                    query_info["group_columns"],
                    self.dataframe_map,
                    self.column_name_map,
                    self.column_to_dataframe_name,
                )
                .transform(having_expr)
                .children[0]
            )

        query_info["having_expr"] = having_expr
        return query_info

    def cross_join(self, df1, df2):
        """
        Returns the crossjoin between two dataframes
        :param df1: Dataframe1
        :param df2: Dataframe2
        :return: Crossjoined dataframe
        """
        temp_key_name = "_cross_join_tempkey"
        df1[temp_key_name] = 1
        df2[temp_key_name] = 1
        new_frame = merge(df1, df2, on=temp_key_name).drop(columns=[temp_key_name])
        df1.drop(columns=[temp_key_name], inplace=True)
        if df1 is not df2:
            df2.drop(columns=[temp_key_name], inplace=True)
        return new_frame

    @staticmethod
    def handle_aggregation(aggregates, group_columns, dataframe: DataFrame):
        """
        Handles all aggregation operations when translating from dictionary info
        to dataframe
        :param aggregates:
        :param group_columns:
        :return:
        """
        if group_columns and not aggregates:
            dataframe = (
                dataframe.groupby(group_columns)
                .size()
                .to_frame("size")
                .reset_index()
                .drop(columns=["size"])
            )
        elif aggregates and not group_columns:
            dataframe = dataframe.aggregate(aggregates).to_frame().transpose()
        elif aggregates and group_columns:
            dataframe = (
                dataframe.groupby(group_columns).aggregate(aggregates).reset_index()
            )

        return dataframe

    def handle_naming(self, columns: list, aliases: dict, first_frame: DataFrame):
        """
        Returns frame with appropriately selected and named columns
        :param columns:
        :param aliases:
        :param first_frame:
        :return:
        """
        column_names = [column.name for column in columns]
        if self.has_star(column_names):
            new_frame: DataFrame = first_frame.copy()
        else:
            column_names = []
            for column in columns:
                true_column_name = self.column_name_map[column.table][
                    column.name.lower()
                ]
                column_names.append(true_column_name)
                if aliases.get(true_column_name) is None:
                    aliases[true_column_name] = column.name
            new_frame: DataFrame = first_frame[column_names].rename(columns=aliases)

        return new_frame

    def to_dataframe(self, query_info):
        """
        Returns the dataframe resulting from the SQL query
        :return:
        """
        if DEBUG:
            print("Query info:", query_info)
        having_expr = query_info["having_expr"]

        frame_names = query_info["frame_names"]
        if not frame_names:
            raise Exception("No table specified")

        first_frame = self.get_frame(frame_names[0])
        for frame_name in frame_names[1:]:
            next_frame = self.get_frame(frame_name)
            first_frame = self.cross_join(first_frame, next_frame)

        new_frame = self.handle_naming(
            query_info["columns"], query_info["aliases"], first_frame
        )

        for expression in query_info["expressions"]:
            new_frame[expression.alias] = expression.evaluate()

        for literal in query_info["literals"]:
            new_frame[literal.alias] = literal.value

        if query_info["conversions"]:
            return new_frame.astype(query_info["conversions"])

        if query_info["where_expr"] is not None:
            new_frame = new_frame[query_info["where_expr"]]

        new_frame = self.handle_aggregation(
            query_info["aggregates"], query_info["group_columns"], new_frame
        )

        if having_expr is not None:
            new_frame = new_frame[having_expr]

        if query_info["distinct"]:
            new_frame.drop_duplicates(keep="first", inplace=True)

        order_by = query_info["order_by"]
        if order_by:
            new_frame.sort_values(
                by=[pair[0] for pair in order_by],
                ascending=[pair[1] for pair in order_by],
                inplace=True,
            )

        if query_info["limit"] is not None:
            new_frame = new_frame.head(query_info["limit"])

        return new_frame

    def set_expr(self, query_info):
        """
        Return different expression with set relational operations performed
        :param query_info:
        :return:
        """
        return self.to_dataframe(query_info).reset_index(drop=True)

    def union_all(self, frame1: DataFrame, frame2: DataFrame):
        """
        Return union all of two dataframes
        :param frame1: Left dataframe
        :param frame2: Right dataframe
        :return:
        """
        return concat([frame1, frame2], ignore_index=True).reset_index(drop=True)

    def union_distinct(self, frame1: DataFrame, frame2: DataFrame):
        """
        Return union all of two dataframes
        :param frame1: Left dataframe
        :param frame2: Right dataframe
        :return:
        """
        return (
            concat([frame1, frame2], ignore_index=True)
            .drop_duplicates()
            .reset_index(drop=True)
        )

    def intersect_distinct(self, frame1: DataFrame, frame2: DataFrame):
        """
        Return intersection of two dataframes
        :param frame1: Left dataframe
        :param frame2: Right dataframe
        :return:
        """
        return merge(
            left=frame1, right=frame2, how="inner", on=frame1.columns.to_list()
        ).reset_index(drop=True)

    def except_distinct(self, frame1: DataFrame, frame2: DataFrame):
        """
        Return first dataframe excluding everything that's also in the second dataframe,
        no duplicates
        :param frame1: Left dataframe
        :param frame2: Right dataframe
        :return:
        """
        return (
            frame1[~frame1.isin(frame2).all(axis=1)]
            .drop_duplicates()
            .reset_index(drop=True)
        )

    def except_all(self, frame1: DataFrame, frame2: DataFrame):
        """
        Return first dataframe excluding everything that's also in the second dataframe,
        with duplicates
        :param frame1: Left dataframe
        :param frame2: Right dataframe
        :return:
        """
        return frame1[~frame1.isin(frame2).all(axis=1)].reset_index(drop=True)

    def final(self, dataframe):
        """
        Returns the final dataframe
        :param dataframe:
        :return:
        """
        return dataframe
