"""
Module containing all sql objects
"""
from pandas import Series


# pylint: disable=too-few-public-methods
class AmbiguousColumn:
    """
    Class for identifying ambiguous table names
    """

    def __init__(self, tables):
        self.tables = tables

    def __repr__(self):
        return f"AmbiguousColumn({','.join(self.tables)})"


class Value:
    """
    Parent class for expressions and columns
    """

    def __init__(self, value, alias='', typename=''):
        self.value = value
        self.alias = alias
        self.typename = typename
        self.final_name = alias

    def __repr__(self):
        if isinstance(self.value, Series):
            print_value = "SeriesObject"
        else:
            print_value = self.value

        display = f"{type(self).__name__}(final_name={self.final_name}, value={print_value}"
        if self.alias:
            display += f", alias={self.alias}"
        if self.typename:
            display += f", type={self.typename}"
        return display

    def __add__(self, other):
        other_name = self.get_other_name(other)
        other = self.get_other_value(other)
        return Expression(value=self.value + other, alias=f'{self.final_name}_add_{other_name}')

    def __sub__(self, other):
        other_name = self.get_other_name(other)
        other = self.get_other_value(other)
        return Expression(value=self.value - other, alias=f'{self.final_name}_sub_{other_name}')

    def __mul__(self, other):
        other_name = self.get_other_name(other)
        other = self.get_other_value(other)
        return Expression(value=self.value * other, alias=f'{self.final_name}_mul_{other_name}')

    def __truediv__(self, other):
        other_name = self.get_other_name(other)
        other = self.get_other_value(other)
        return Expression(value=self.value / other, alias=f"{self.final_name}_div_{other_name}")

    @staticmethod
    def get_other_name(other):
        """
        Gets the name representation for the other value
        :param other:
        :return:
        """
        if isinstance(other, Column):
            return other.name
        if isinstance(other, Expression):
            return other.alias
        if isinstance(other, Literal):
            return str(other.value)
        return str(other)

    @staticmethod
    def get_other_value(other):
        """
        Return the appropriate value based on the type of other
        :param other:
        :return:
        """
        if isinstance(other, (Literal, Column, Expression)):
            return other.value
        return other

    def set_alias(self, alias):
        """
        Sets the alias and final name for the value object
        :param alias:
        :return:
        """
        self.alias = alias
        self.final_name = alias


class Literal(Value):
    """
    Stores literal data
    """

    literal_count = 0

    def __init__(self, value, alias=''):
        super(Literal, self).__init__(value, alias)
        if not alias:
            self.alias = f"_literal{self.literal_count}"
            self.literal_count += 1

    def __gt__(self, other):
        if isinstance(other, Literal):
            return self.value > other.value
        return self.value > other

    def __lt__(self, other):
        if isinstance(other, Literal):
            return self.value < other.value
        return self.value < other

    def __ge__(self, other):
        if isinstance(other, Literal):
            return self.value >= other.value
        return self.value >= other

    def __le__(self, other):
        if isinstance(other, Literal):
            return self.value <= other.value
        return self.value <= other

    def __repr__(self):
        return super(Literal, self).__repr__() + ")"


class Number(Literal):
    """
    Stores numerical data
    """

    def __init__(self, value):
        super(Number, self).__init__(value)


class String(Literal):
    """
    Store information about a string literal
    """

    def __init__(self, value):
        super(String, self).__init__(value)


class Date(Literal):
    """
    Store information about a date literal
    """

    def __init__(self, value):
        super(Date, self).__init__(value)


class Bool(Literal):
    """
    Store information about a date literal
    """

    def __init__(self, value):
        super(Bool, self).__init__(value)


class Expression(Value):
    """
    Store information about an expression
    """
    expressions = 0

    def __init__(self, value, alias='', typename='', function=''):
        super(Expression, self).__init__(value, alias, typename)
        self.function = function
        if self.alias:
            self.final_name = self.alias
        else:
            if isinstance(self.value, Series):
                self.final_name = f"_expression{self.expressions}"
                self.alias = self.final_name
                self.expressions += 1
            else:
                self.final_name = str(self.value)
            if self.function:
                if isinstance(self.value, Column):
                    expression_name = self.value.name
                else:
                    expression_name = str(self.value)
                self.alias = self.function + "_" + expression_name
                self.final_name = self.alias
        self.has_columns = True

    def __repr__(self):
        display = super(Expression, self).__repr__()
        if self.function:
            display += f", function={self.function}"
        return display + ")"

    def evaluate(self):
        """
        Returns the value from the expression
        :return:
        """
        if isinstance(self.value, Column):
            return self.value.value
        return self.value


class Column(Value):
    """
    Store information about columns
    """

    def __init__(self, name: str, alias='', typename='', value=None):
        super(Column, self).__init__(value, alias, typename)
        self.name = name
        if self.alias:
            self.final_name = self.alias
        else:
            self.final_name = self.name
        self.table = None

    def __repr__(self):
        display = super(Column, self).__repr__()
        display += f", name={self.name}"
        display += f", table={self.table}"
        return display + ")"

    def __eq__(self, other):
        other = self.get_other_value(other)
        return self.value == other

    def __gt__(self, other):
        other = self.get_other_value(other)
        return self.value > other

    def __lt__(self, other):
        other = self.get_other_value(other)
        return self.value < other

    def __ge__(self, other):
        other = self.get_other_value(other)
        return self.value >= other

    def __le__(self, other):
        other = self.get_other_value(other)
        return self.value <= other

    def set_value(self, new_value: Series):
        """
        Set the value of the column to value
        :param new_value:
        :return:
        """
        self.value = new_value


class Subquery:
    """
    Wrapper for subqueries
    """

    def __init__(self, name: str, query_info: dict):
        self.name = name
        self.query_info = query_info

    def __repr__(self):
        return f"Subquery(name={self.name}, query_info={self.query_info})"
