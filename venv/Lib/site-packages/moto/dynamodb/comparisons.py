import re
from collections import deque, namedtuple
from typing import Any, Deque, Dict, Iterable, List, Optional, Tuple, Union

from moto.dynamodb.exceptions import ConditionAttributeIsReservedKeyword
from moto.dynamodb.models.dynamo_type import Item
from moto.dynamodb.parsing.reserved_keywords import ReservedKeywords


def create_condition_expression_parser(
    expr: Optional[str],
    names: Optional[Dict[str, str]],
    values: Optional[Dict[str, Dict[str, str]]],
) -> "ConditionExpressionParser":
    return ConditionExpressionParser(expr, names, values)


def get_filter_expression(
    expr: Optional[str],
    names: Optional[Dict[str, str]],
    values: Optional[Dict[str, Dict[str, str]]],
) -> Union["Op", "Func"]:
    """
    Parse a filter expression into an Op.

    Examples
        expr = 'Id > 5 AND attribute_exists(test) AND Id BETWEEN 5 AND 6 OR length < 6 AND contains(test, 1) AND 5 IN (4,5, 6) OR (Id < 5 AND 5 > Id)'
        expr = 'Id > 5 AND Subs < 7'
    """
    parser = create_condition_expression_parser(expr, names, values)
    return parser.parse()


def get_expected(expected: Dict[str, Any]) -> Union["Op", "Func"]:
    """
    Parse a filter expression into an Op.

    Examples
        expr = 'Id > 5 AND attribute_exists(test) AND Id BETWEEN 5 AND 6 OR length < 6 AND contains(test, 1) AND 5 IN (4,5, 6) OR (Id < 5 AND 5 > Id)'
        expr = 'Id > 5 AND Subs < 7'
    """
    ops: Dict[str, Any] = {
        "EQ": OpEqual,
        "NE": OpNotEqual,
        "LE": OpLessThanOrEqual,
        "LT": OpLessThan,
        "GE": OpGreaterThanOrEqual,
        "GT": OpGreaterThan,
        "NOT_NULL": FuncAttrExists,
        "NULL": FuncAttrNotExists,
        "CONTAINS": FuncContains,
        "NOT_CONTAINS": FuncNotContains,
        "BEGINS_WITH": FuncBeginsWith,
        "IN": FuncIn,
        "BETWEEN": FuncBetween,
    }

    # NOTE: Always uses ConditionalOperator=AND
    conditions: List[Union["Op", "Func"]] = []
    for key, cond in expected.items():
        path = AttributePath([key])
        if "Exists" in cond:
            if cond["Exists"]:
                conditions.append(FuncAttrExists(path))
            else:
                conditions.append(FuncAttrNotExists(path))
        elif "Value" in cond:
            conditions.append(OpEqual(path, AttributeValue(cond["Value"])))
        elif "ComparisonOperator" in cond:
            operator_name = cond["ComparisonOperator"]
            values = [AttributeValue(v) for v in cond.get("AttributeValueList", [])]
            OpClass = ops[operator_name]
            conditions.append(OpClass(path, *values))

    # NOTE: Ignore ConditionalOperator
    ConditionalOp = OpAnd
    if conditions:
        output = conditions[0]
        for condition in conditions[1:]:
            output = ConditionalOp(output, condition)
    else:
        return OpDefault(None, None)  # type: ignore[arg-type]

    return output


class Op:
    """
    Base class for a FilterExpression operator
    """

    OP = ""

    def __init__(
        self, lhs: Union["Func", "Op", "Operand"], rhs: Union["Func", "Op", "Operand"]
    ):
        self.lhs = lhs
        self.rhs = rhs

    def expr(self, item: Optional[Item]) -> bool:
        raise NotImplementedError(f"Expr not defined for {type(self)}")

    def __repr__(self) -> str:
        return f"({self.lhs} {self.OP} {self.rhs})"


# TODO add tests for all of these

EQ_FUNCTION = lambda item_value, test_value: item_value == test_value  # noqa
NE_FUNCTION = lambda item_value, test_value: item_value != test_value  # noqa
LE_FUNCTION = lambda item_value, test_value: item_value <= test_value  # noqa
LT_FUNCTION = lambda item_value, test_value: item_value < test_value  # noqa
GE_FUNCTION = lambda item_value, test_value: item_value >= test_value  # noqa
GT_FUNCTION = lambda item_value, test_value: item_value > test_value  # noqa

COMPARISON_FUNCS = {
    "EQ": EQ_FUNCTION,
    "=": EQ_FUNCTION,
    "NE": NE_FUNCTION,
    "!=": NE_FUNCTION,
    "LE": LE_FUNCTION,
    "<=": LE_FUNCTION,
    "LT": LT_FUNCTION,
    "<": LT_FUNCTION,
    "GE": GE_FUNCTION,
    ">=": GE_FUNCTION,
    "GT": GT_FUNCTION,
    ">": GT_FUNCTION,
    # NULL means the value should not exist at all
    "NULL": lambda item_value: False,
    # NOT_NULL means the value merely has to exist, and values of None are valid
    "NOT_NULL": lambda item_value: True,
    "CONTAINS": lambda item_value, test_value: test_value in item_value,
    "NOT_CONTAINS": lambda item_value, test_value: test_value not in item_value,
    "BEGINS_WITH": lambda item_value, test_value: item_value.startswith(test_value),
    "IN": lambda item_value, *test_values: item_value in test_values,
    "BETWEEN": lambda item_value, lower_test_value, upper_test_value: lower_test_value
    <= item_value
    <= upper_test_value,
}


def get_comparison_func(range_comparison: str) -> Any:
    return COMPARISON_FUNCS.get(range_comparison)


class RecursionStopIteration(StopIteration):
    pass


class ConditionExpressionParser:
    def __init__(
        self,
        condition_expression: Optional[str],
        expression_attribute_names: Optional[Dict[str, str]],
        expression_attribute_values: Optional[Dict[str, Dict[str, str]]],
    ):
        self.condition_expression = condition_expression
        self.expression_attribute_names = expression_attribute_names
        self.expression_attribute_values = expression_attribute_values

        self.expr_attr_names_found: List[str] = []

    def parse(self) -> Union[Op, "Func"]:
        """Returns a syntax tree for the expression.

        The tree, and all of the nodes in the tree are a tuple of
        - kind: str
        - children/value:
            list of nodes for parent nodes
            value for leaf nodes

        Raises ValueError if the condition expression is invalid
        Raises KeyError if expression attribute names/values are invalid

        Here are the types of nodes that can be returned.
        The types of child nodes are denoted with a colon (:).
        An arbitrary number of children is denoted with ...

        Condition:
            ('OR', [lhs : Condition, rhs : Condition])
            ('AND', [lhs: Condition, rhs: Condition])
            ('NOT', [argument: Condition])
            ('PARENTHESES', [argument: Condition])
            ('FUNCTION', [('LITERAL', function_name: str), argument: Operand, ...])
            ('BETWEEN', [query: Operand, low: Operand, high: Operand])
            ('IN', [query: Operand, possible_value: Operand, ...])
            ('COMPARISON', [lhs: Operand, ('LITERAL', comparator: str), rhs: Operand])

        Operand:
            ('EXPRESSION_ATTRIBUTE_VALUE', value: dict, e.g. {'S': 'foobar'})
            ('PATH', [('LITERAL', path_element: str), ...])
            NOTE: Expression attribute names will be expanded
            ('FUNCTION', [('LITERAL', 'size'), argument: Operand])

        Literal:
            ('LITERAL', value: str)

        """
        if not self.condition_expression:
            return OpDefault(None, None)  # type: ignore[arg-type]
        nodes = self._lex_condition_expression()
        nodes = self._parse_paths(nodes)
        # NOTE: The docs say that functions should be parsed after
        # IN, BETWEEN, and comparisons like <=.
        # However, these expressions are invalid as function arguments,
        # so it is okay to parse functions first. This needs to be done
        # to interpret size() correctly as an operand.
        nodes = self._apply_functions(nodes)
        nodes = self._apply_comparator(nodes)
        nodes = self._apply_in(nodes)
        nodes = self._apply_between(nodes)
        nodes = self._apply_parens_and_booleans(nodes)
        node = nodes[0]

        self.expr_attr_names_found.extend(self._find_literals(node))

        op = self._make_op_condition(node)
        return op

    class Kind:
        """Enum defining types of nodes in the syntax tree."""

        # Condition nodes
        # ---------------
        OR = "OR"
        AND = "AND"
        NOT = "NOT"
        PARENTHESES = "PARENTHESES"
        FUNCTION = "FUNCTION"
        BETWEEN = "BETWEEN"
        IN = "IN"
        COMPARISON = "COMPARISON"

        # Operand nodes
        # -------------
        EXPRESSION_ATTRIBUTE_VALUE = "EXPRESSION_ATTRIBUTE_VALUE"
        PATH = "PATH"

        # Literal nodes
        # --------------
        LITERAL = "LITERAL"

    class Nonterminal:
        """Enum defining nonterminals for productions."""

        CONDITION = "CONDITION"
        OPERAND = "OPERAND"
        COMPARATOR = "COMPARATOR"
        FUNCTION_NAME = "FUNCTION_NAME"
        IDENTIFIER = "IDENTIFIER"
        AND = "AND"
        OR = "OR"
        NOT = "NOT"
        BETWEEN = "BETWEEN"
        IN = "IN"
        COMMA = "COMMA"
        LEFT_PAREN = "LEFT_PAREN"
        RIGHT_PAREN = "RIGHT_PAREN"
        WHITESPACE = "WHITESPACE"

    Node = namedtuple("Node", ["nonterminal", "kind", "text", "value", "children"])

    @classmethod
    def _find_literals(cls, parent: Node) -> List[str]:  # type: ignore
        literals: List[str] = []
        if parent.kind == "LITERAL" and parent.nonterminal == "IDENTIFIER":
            literals.append(parent.text)
        else:
            for child in parent.children:
                literals.extend(cls._find_literals(child))
        return literals

    @classmethod
    def raise_exception_if_keyword(cls, attribute: str) -> None:
        if attribute.upper() in ReservedKeywords.get_reserved_keywords():
            raise ConditionAttributeIsReservedKeyword(attribute)

    def _lex_condition_expression(self) -> Deque[Node]:
        nodes: Deque[ConditionExpressionParser.Node] = deque()
        remaining_expression = self.condition_expression
        while remaining_expression:
            node, remaining_expression = self._lex_one_node(remaining_expression)
            if node.nonterminal == self.Nonterminal.WHITESPACE:
                continue
            nodes.append(node)
        return nodes

    def _lex_one_node(self, remaining_expression: str) -> Tuple[Node, str]:
        # TODO: Handle indexing like [1]
        attribute_regex = r"(:|#)?[A-z0-9\-_]+"
        patterns = [
            (self.Nonterminal.WHITESPACE, re.compile(r"^ +")),
            (
                self.Nonterminal.COMPARATOR,
                re.compile(
                    "^("
                    # Put long expressions first for greedy matching
                    "<>|"
                    "<=|"
                    ">=|"
                    "=|"
                    "<|"
                    ">)"
                ),
            ),
            (
                self.Nonterminal.OPERAND,
                re.compile(rf"^{attribute_regex}(\.{attribute_regex}|\[[0-9]\])*"),
            ),
            (self.Nonterminal.COMMA, re.compile(r"^,")),
            (self.Nonterminal.LEFT_PAREN, re.compile(r"^\(")),
            (self.Nonterminal.RIGHT_PAREN, re.compile(r"^\)")),
        ]

        for nonterminal, pattern in patterns:
            match = pattern.match(remaining_expression)
            if match:
                match_text = match.group()
                break
        else:  # pragma: no cover
            raise ValueError(
                f"Cannot parse condition starting at:{remaining_expression}"
            )

        node = self.Node(
            nonterminal=nonterminal,
            kind=self.Kind.LITERAL,
            text=match_text,
            value=match_text,
            children=[],
        )

        remaining_expression = remaining_expression[len(match_text) :]

        return node, remaining_expression

    def _parse_paths(self, nodes: Deque[Node]) -> Deque[Node]:
        output: Deque[ConditionExpressionParser.Node] = deque()

        while nodes:
            node = nodes.popleft()

            if node.nonterminal == self.Nonterminal.OPERAND:
                path = node.value.replace("[", ".[").split(".")
                children = [self._parse_path_element(name) for name in path]
                if len(children) == 1:
                    child = children[0]
                    if child.nonterminal != self.Nonterminal.IDENTIFIER:
                        output.append(child)
                        continue
                else:
                    for child in children:
                        self._assert(
                            child.nonterminal == self.Nonterminal.IDENTIFIER,
                            f"Cannot use {child.text} in path",
                            [node],
                        )
                output.append(
                    self.Node(
                        nonterminal=self.Nonterminal.OPERAND,
                        kind=self.Kind.PATH,
                        text=node.text,
                        value=None,
                        children=children,
                    )
                )
            else:
                output.append(node)
        return output

    def _parse_path_element(self, name: str) -> Node:
        reserved = {
            "and": self.Nonterminal.AND,
            "or": self.Nonterminal.OR,
            "in": self.Nonterminal.IN,
            "between": self.Nonterminal.BETWEEN,
            "not": self.Nonterminal.NOT,
        }

        functions = {
            "attribute_exists",
            "attribute_not_exists",
            "attribute_type",
            "begins_with",
            "contains",
            "size",
        }

        if name.lower() in reserved:
            # e.g. AND
            nonterminal = reserved[name.lower()]
            return self.Node(
                nonterminal=nonterminal,
                kind=self.Kind.LITERAL,
                text=name,
                value=name,
                children=[],
            )
        elif name in functions:
            # e.g. attribute_exists
            return self.Node(
                nonterminal=self.Nonterminal.FUNCTION_NAME,
                kind=self.Kind.LITERAL,
                text=name,
                value=name,
                children=[],
            )
        elif name.startswith(":"):
            # e.g. :value0
            return self.Node(
                nonterminal=self.Nonterminal.OPERAND,
                kind=self.Kind.EXPRESSION_ATTRIBUTE_VALUE,
                text=name,
                value=self._lookup_expression_attribute_value(name),
                children=[],
            )
        elif name.startswith("#"):
            # e.g. #name0
            return self.Node(
                nonterminal=self.Nonterminal.IDENTIFIER,
                kind=self.Kind.LITERAL,
                text=name,
                value=self._lookup_expression_attribute_name(name),
                children=[],
            )
        elif name.startswith("["):
            # e.g. [123]
            if not name.endswith("]"):  # pragma: no cover
                raise ValueError(f"Bad path element {name}")
            return self.Node(
                nonterminal=self.Nonterminal.IDENTIFIER,
                kind=self.Kind.LITERAL,
                text=name,
                value=int(name[1:-1]),
                children=[],
            )
        else:
            # e.g. ItemId
            self.raise_exception_if_keyword(name)
            return self.Node(
                nonterminal=self.Nonterminal.IDENTIFIER,
                kind=self.Kind.LITERAL,
                text=name,
                value=name,
                children=[],
            )

    def _lookup_expression_attribute_value(self, name: str) -> Dict[str, str]:
        return self.expression_attribute_values[name]  # type: ignore[index]

    def _lookup_expression_attribute_name(self, name: str) -> str:
        return self.expression_attribute_names[name]  # type: ignore[index]

    # NOTE: The following constructions are ordered from high precedence to low precedence
    # according to
    # https://docs.aws.amazon.com/amazondynamodb/latest/developerguide/Expressions.OperatorsAndFunctions.html#Expressions.OperatorsAndFunctions.Precedence
    #
    # = <> < <= > >=
    # IN
    # BETWEEN
    # attribute_exists attribute_not_exists begins_with contains
    # Parentheses
    # NOT
    # AND
    # OR
    #
    # The grammar is taken from
    # https://docs.aws.amazon.com/amazondynamodb/latest/developerguide/Expressions.OperatorsAndFunctions.html#Expressions.OperatorsAndFunctions.Syntax
    #
    # condition-expression ::=
    #     operand comparator operand
    #     operand BETWEEN operand AND operand
    #     operand IN ( operand (',' operand (, ...) ))
    #     function
    #     condition AND condition
    #     condition OR condition
    #     NOT condition
    #     ( condition )
    #
    # comparator ::=
    #     =
    #     <>
    #     <
    #     <=
    #     >
    #     >=
    #
    # function ::=
    #     attribute_exists (path)
    #     attribute_not_exists (path)
    #     attribute_type (path, type)
    #     begins_with (path, substr)
    #     contains (path, operand)
    #     size (path)

    def _matches(self, nodes: Deque[Node], production: List[str]) -> bool:
        """Check if the nodes start with the given production.

        Parameters
        ----------
        nodes: list of Node
        production: list of str
            The name of a Nonterminal, or '*' for anything

        """
        if len(nodes) < len(production):
            return False
        for i in range(len(production)):
            if production[i] == "*":
                continue
            expected = getattr(self.Nonterminal, production[i])
            if nodes[i].nonterminal != expected:
                return False
        return True

    def _apply_comparator(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := operand comparator operand."""
        output: Deque[ConditionExpressionParser.Node] = deque()

        while nodes:
            if self._matches(nodes, ["*", "COMPARATOR"]):
                self._assert(
                    self._matches(nodes, ["OPERAND", "COMPARATOR", "OPERAND"]),
                    "Bad comparison",
                    list(nodes)[:3],
                )
                lhs = nodes.popleft()
                comparator = nodes.popleft()
                rhs = nodes.popleft()
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.COMPARISON,
                        text=" ".join([lhs.text, comparator.text, rhs.text]),
                        value=None,
                        children=[lhs, comparator, rhs],
                    )
                )
            else:
                output.append(nodes.popleft())
        return output

    def _apply_in(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := operand IN ( operand , ... )."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["*", "IN"]):
                self._assert(
                    self._matches(nodes, ["OPERAND", "IN", "LEFT_PAREN"]),
                    "Bad IN expression",
                    list(nodes)[:3],
                )
                lhs = nodes.popleft()
                in_node = nodes.popleft()
                left_paren = nodes.popleft()
                all_children = [lhs, in_node, left_paren]
                rhs = []
                while True:
                    if self._matches(nodes, ["OPERAND", "COMMA"]):
                        operand = nodes.popleft()
                        separator = nodes.popleft()
                        all_children += [operand, separator]
                        rhs.append(operand)
                    elif self._matches(nodes, ["OPERAND", "RIGHT_PAREN"]):
                        operand = nodes.popleft()
                        separator = nodes.popleft()
                        all_children += [operand, separator]
                        rhs.append(operand)
                        break  # Close
                    else:
                        self._assert(False, "Bad IN expression starting at", nodes)
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.IN,
                        text=" ".join([t.text for t in all_children]),
                        value=None,
                        children=[lhs] + rhs,
                    )
                )
            else:
                output.append(nodes.popleft())
        return output

    def _apply_between(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := operand BETWEEN operand AND operand."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["*", "BETWEEN"]):
                self._assert(
                    self._matches(
                        nodes, ["OPERAND", "BETWEEN", "OPERAND", "AND", "OPERAND"]
                    ),
                    "Bad BETWEEN expression",
                    list(nodes)[:5],
                )
                lhs = nodes.popleft()
                between_node = nodes.popleft()
                low = nodes.popleft()
                and_node = nodes.popleft()
                high = nodes.popleft()
                all_children = [lhs, between_node, low, and_node, high]
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.BETWEEN,
                        text=" ".join([t.text for t in all_children]),
                        value=None,
                        children=[lhs, low, high],
                    )
                )
            else:
                output.append(nodes.popleft())
        return output

    def _apply_functions(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := function_name (operand , ...)."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        either_kind = {self.Kind.PATH, self.Kind.EXPRESSION_ATTRIBUTE_VALUE}
        expected_argument_kind_map = {
            "attribute_exists": [{self.Kind.PATH}],
            "attribute_not_exists": [{self.Kind.PATH}],
            "attribute_type": [either_kind, {self.Kind.EXPRESSION_ATTRIBUTE_VALUE}],
            "begins_with": [either_kind, either_kind],
            "contains": [either_kind, either_kind],
            "size": [{self.Kind.PATH}],
        }
        while nodes:
            if self._matches(nodes, ["FUNCTION_NAME"]):
                self._assert(
                    self._matches(
                        nodes, ["FUNCTION_NAME", "LEFT_PAREN", "OPERAND", "*"]
                    ),
                    "Bad function expression at",
                    list(nodes)[:4],
                )
                function_name = nodes.popleft()
                left_paren = nodes.popleft()
                all_children = [function_name, left_paren]
                arguments = []
                while True:
                    if self._matches(nodes, ["OPERAND", "COMMA"]):
                        operand = nodes.popleft()
                        separator = nodes.popleft()
                        all_children += [operand, separator]
                        arguments.append(operand)
                    elif self._matches(nodes, ["OPERAND", "RIGHT_PAREN"]):
                        operand = nodes.popleft()
                        separator = nodes.popleft()
                        all_children += [operand, separator]
                        arguments.append(operand)
                        break  # Close paren
                    else:
                        self._assert(
                            False,
                            "Bad function expression",
                            all_children + list(nodes)[:2],
                        )
                expected_kinds = expected_argument_kind_map[function_name.value]
                self._assert(
                    len(arguments) == len(expected_kinds),
                    "Wrong number of arguments in",
                    all_children,
                )
                for i in range(len(expected_kinds)):
                    self._assert(
                        arguments[i].kind in expected_kinds[i],
                        f"Wrong type for argument {i} in",
                        all_children,
                    )
                if function_name.value == "size":
                    nonterminal = self.Nonterminal.OPERAND
                else:
                    nonterminal = self.Nonterminal.CONDITION
                nodes.appendleft(
                    self.Node(
                        nonterminal=nonterminal,
                        kind=self.Kind.FUNCTION,
                        text=" ".join([t.text for t in all_children]),
                        value=None,
                        children=[function_name] + arguments,
                    )
                )
            else:
                output.append(nodes.popleft())
        return output

    def _apply_parens_and_booleans(
        self, nodes: Deque[Node], left_paren: Any = None
    ) -> Deque[Node]:
        """Apply condition := ( condition ) and booleans."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["LEFT_PAREN"]):
                parsed = self._apply_parens_and_booleans(
                    nodes, left_paren=nodes.popleft()
                )
                self._assert(len(parsed) >= 1, "Failed to close parentheses at", nodes)
                parens = parsed.popleft()
                self._assert(
                    parens.kind == self.Kind.PARENTHESES,
                    "Failed to close parentheses at",
                    nodes,
                )
                output.append(parens)
                nodes = parsed
            elif self._matches(nodes, ["RIGHT_PAREN"]):
                self._assert(left_paren is not None, "Unmatched ) at", nodes)
                close_paren = nodes.popleft()
                children = self._apply_booleans(output)
                all_children = [left_paren] + list(children) + [close_paren]
                return deque(
                    [
                        self.Node(
                            nonterminal=self.Nonterminal.CONDITION,
                            kind=self.Kind.PARENTHESES,
                            text=" ".join([t.text for t in all_children]),
                            value=None,
                            children=list(children),
                        )
                    ]
                    + list(nodes)
                )
            else:
                output.append(nodes.popleft())

        self._assert(left_paren is None, "Unmatched ( at", list(output))
        return self._apply_booleans(output)

    def _apply_booleans(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply and, or, and not constructions."""
        nodes = self._apply_not(nodes)
        nodes = self._apply_and(nodes)
        nodes = self._apply_or(nodes)
        # The expression should reduce to a single condition
        self._assert(len(nodes) == 1, "Unexpected expression at", list(nodes)[1:])
        self._assert(
            nodes[0].nonterminal == self.Nonterminal.CONDITION,
            "Incomplete condition",
            nodes,
        )
        return nodes

    def _apply_not(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := NOT condition."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["NOT"]):
                self._assert(
                    self._matches(nodes, ["NOT", "CONDITION"]),
                    "Bad NOT expression",
                    list(nodes)[:2],
                )
                not_node = nodes.popleft()
                child = nodes.popleft()
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.NOT,
                        text=" ".join([not_node.text, child.text]),
                        value=None,
                        children=[child],
                    )
                )
            else:
                output.append(nodes.popleft())

        return output

    def _apply_and(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := condition AND condition."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["*", "AND"]):
                self._assert(
                    self._matches(nodes, ["CONDITION", "AND", "CONDITION"]),
                    "Bad AND expression",
                    list(nodes)[:3],
                )
                lhs = nodes.popleft()
                and_node = nodes.popleft()
                rhs = nodes.popleft()
                all_children = [lhs, and_node, rhs]
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.AND,
                        text=" ".join([t.text for t in all_children]),
                        value=None,
                        children=[lhs, rhs],
                    )
                )
            else:
                output.append(nodes.popleft())

        return output

    def _apply_or(self, nodes: Deque[Node]) -> Deque[Node]:
        """Apply condition := condition OR condition."""
        output: Deque[ConditionExpressionParser.Node] = deque()
        while nodes:
            if self._matches(nodes, ["*", "OR"]):
                self._assert(
                    self._matches(nodes, ["CONDITION", "OR", "CONDITION"]),
                    "Bad OR expression",
                    list(nodes)[:3],
                )
                lhs = nodes.popleft()
                or_node = nodes.popleft()
                rhs = nodes.popleft()
                all_children = [lhs, or_node, rhs]
                nodes.appendleft(
                    self.Node(
                        nonterminal=self.Nonterminal.CONDITION,
                        kind=self.Kind.OR,
                        text=" ".join([t.text for t in all_children]),
                        value=None,
                        children=[lhs, rhs],
                    )
                )
            else:
                output.append(nodes.popleft())

        return output

    def _make_operand(self, node: Node) -> "Operand":
        if node.kind == self.Kind.PATH:
            return AttributePath([child.value for child in node.children])
        elif node.kind == self.Kind.EXPRESSION_ATTRIBUTE_VALUE:
            return AttributeValue(node.value)
        elif node.kind == self.Kind.FUNCTION:
            # size()
            function_node = node.children[0]
            arguments = node.children[1:]
            function_name = function_node.value
            arguments = [self._make_operand(arg) for arg in arguments]
            return FUNC_CLASS[function_name](*arguments)
        else:  # pragma: no cover
            raise ValueError(f"Unknown operand: {node}")

    def _make_op_condition(self, node: Node) -> Union["Func", Op]:
        if node.kind == self.Kind.OR:
            lhs, rhs = node.children
            return OpOr(self._make_op_condition(lhs), self._make_op_condition(rhs))
        elif node.kind == self.Kind.AND:
            lhs, rhs = node.children
            return OpAnd(self._make_op_condition(lhs), self._make_op_condition(rhs))
        elif node.kind == self.Kind.NOT:
            (child,) = node.children
            return OpNot(self._make_op_condition(child))
        elif node.kind == self.Kind.PARENTHESES:
            (child,) = node.children
            return self._make_op_condition(child)
        elif node.kind == self.Kind.FUNCTION:
            function_node = node.children[0]
            arguments = node.children[1:]
            function_name = function_node.value
            arguments = [self._make_operand(arg) for arg in arguments]
            return FUNC_CLASS[function_name](*arguments)
        elif node.kind == self.Kind.BETWEEN:
            query, low, high = node.children
            return FuncBetween(
                self._make_operand(query),
                self._make_operand(low),
                self._make_operand(high),
            )
        elif node.kind == self.Kind.IN:
            query = node.children[0]
            possible_values = node.children[1:]
            query = self._make_operand(query)
            possible_values = [self._make_operand(v) for v in possible_values]
            return FuncIn(query, *possible_values)
        elif node.kind == self.Kind.COMPARISON:
            lhs, comparator, rhs = node.children
            return COMPARATOR_CLASS[comparator.value](
                self._make_operand(lhs), self._make_operand(rhs)
            )
        else:  # pragma: no cover
            raise ValueError(f"Unknown expression node kind {node.kind}")

    def _assert(self, condition: bool, message: str, nodes: Iterable[Node]) -> None:
        if not condition:
            raise ValueError(message + " " + " ".join([t.text for t in nodes]))


class Operand:
    def expr(self, item: Optional[Item]) -> Any:
        raise NotImplementedError

    def get_type(self, item: Optional[Item]) -> Optional[str]:
        raise NotImplementedError


class AttributePath(Operand):
    def __init__(self, path: List[Any]):
        """Initialize the AttributePath.

        Parameters
        ----------
        path: list of int/str

        """
        assert len(path) >= 1
        self.path = path

    def _get_attr(self, item: Optional[Item]) -> Any:
        if item is None:
            return None

        base = self.path[0]
        if base not in item.attrs:
            return None
        attr = item.attrs[base]

        for name in self.path[1:]:
            attr = attr.child_attr(name)
            if attr is None:
                return None

        return attr

    def expr(self, item: Optional[Item]) -> Any:
        attr = self._get_attr(item)
        if attr is None:
            return None
        else:
            return attr.cast_value

    def get_type(self, item: Optional[Item]) -> Optional[str]:
        attr = self._get_attr(item)
        if attr is None:
            return None
        else:
            return attr.type

    def __repr__(self) -> str:
        return ".".join(self.path)


class AttributeValue(Operand):
    def __init__(self, value: Dict[str, Any]):
        """Initialize the AttributePath.

        Parameters
        ----------
        value: dict
            e.g. {'N': '1.234'}

        """
        self.type = list(value.keys())[0]
        self.value = value[self.type]

    def expr(self, item: Optional[Item]) -> Any:
        # TODO: Reuse DynamoType code
        if self.type == "N":
            try:
                return int(self.value)
            except ValueError:
                return float(self.value)
        elif self.type in ["SS", "NS", "BS"]:
            sub_type = self.type[0]
            return set([AttributeValue({sub_type: v}).expr(item) for v in self.value])
        elif self.type == "L":
            return [AttributeValue(v).expr(item) for v in self.value]
        elif self.type == "M":
            return dict(
                [(k, AttributeValue(v).expr(item)) for k, v in self.value.items()]
            )
        else:
            return self.value
        return self.value

    def get_type(self, item: Optional[Item]) -> str:
        return self.type

    def __repr__(self) -> str:
        return repr(self.value)


class OpDefault(Op):
    OP = "NONE"

    def expr(self, item: Optional[Item]) -> bool:
        """If no condition is specified, always True."""
        return True


class OpNot(Op):
    OP = "NOT"

    def __init__(self, lhs: Union["Func", Op]):
        super().__init__(lhs, None)  # type: ignore[arg-type]

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        return not lhs

    def __str__(self) -> str:
        return f"({self.OP} {self.lhs})"


class OpAnd(Op):
    OP = "AND"

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        return lhs and self.rhs.expr(item)


class OpLessThan(Op):
    OP = "<"

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        # In python3 None is not a valid comparator when using < or > so must be handled specially
        if lhs is not None and rhs is not None:
            return lhs < rhs
        else:
            return False


class OpGreaterThan(Op):
    OP = ">"

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        # In python3 None is not a valid comparator when using < or > so must be handled specially
        if lhs is not None and rhs is not None:
            return lhs > rhs
        else:
            return False


class OpEqual(Op):
    OP = "="

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        return lhs == rhs


class OpNotEqual(Op):
    OP = "<>"

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        return lhs != rhs


class OpLessThanOrEqual(Op):
    OP = "<="

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        # In python3 None is not a valid comparator when using < or > so must be handled specially
        if lhs is not None and rhs is not None:
            return lhs <= rhs
        else:
            return False


class OpGreaterThanOrEqual(Op):
    OP = ">="

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        rhs = self.rhs.expr(item)
        # In python3 None is not a valid comparator when using < or > so must be handled specially
        if lhs is not None and rhs is not None:
            return lhs >= rhs
        else:
            return False


class OpOr(Op):
    OP = "OR"

    def expr(self, item: Optional[Item]) -> bool:
        lhs = self.lhs.expr(item)
        return lhs or self.rhs.expr(item)


class Func:
    """
    Base class for a FilterExpression function
    """

    FUNC = "Unknown"

    def __init__(self, *arguments: Any):
        self.arguments = arguments

    def expr(self, item: Optional[Item]) -> bool:
        raise NotImplementedError

    def __repr__(self) -> str:
        return f"{self.FUNC}({' '.join([repr(arg) for arg in self.arguments])})"


class FuncAttrExists(Func):
    FUNC = "attribute_exists"

    def __init__(self, attribute: Operand):
        self.attr = attribute
        super().__init__(attribute)

    def expr(self, item: Optional[Item]) -> bool:
        return self.attr.get_type(item) is not None


def FuncAttrNotExists(attribute: Operand) -> Any:
    return OpNot(FuncAttrExists(attribute))


class FuncAttrType(Func):
    FUNC = "attribute_type"

    def __init__(self, attribute: Operand, _type: Func):
        self.attr = attribute
        self.type = _type
        super().__init__(attribute, _type)

    def expr(self, item: Optional[Item]) -> bool:
        return self.attr.get_type(item) == self.type.expr(item)  # type: ignore[comparison-overlap]


class FuncBeginsWith(Func):
    FUNC = "begins_with"

    def __init__(self, attribute: Operand, substr: Operand):
        self.attr = attribute
        self.substr = substr
        super().__init__(attribute, substr)

    def expr(self, item: Optional[Item]) -> bool:
        if self.attr.get_type(item) != "S":
            return False
        if self.substr.get_type(item) != "S":
            return False
        return self.attr.expr(item).startswith(self.substr.expr(item))


class FuncContains(Func):
    FUNC = "contains"

    def __init__(self, attribute: Operand, operand: Operand):
        self.attr = attribute
        self.operand = operand
        super().__init__(attribute, operand)

    def expr(self, item: Optional[Item]) -> bool:
        if self.attr.get_type(item) in ("S", "SS", "NS", "BS", "L"):
            try:
                return self.operand.expr(item) in self.attr.expr(item)
            except TypeError:
                return False
        return False


def FuncNotContains(attribute: Operand, operand: Operand) -> OpNot:
    return OpNot(FuncContains(attribute, operand))


class FuncSize(Func):
    FUNC = "size"

    def __init__(self, attribute: Operand):
        self.attr = attribute
        super().__init__(attribute)

    def expr(self, item: Optional[Item]) -> int:  # type: ignore[override]
        if self.attr.get_type(item) is None:
            raise ValueError(f"Invalid attribute name {self.attr}")

        if self.attr.get_type(item) in ("S", "SS", "NS", "B", "BS", "L", "M"):
            return len(self.attr.expr(item))
        raise ValueError("Invalid filter expression")


class FuncBetween(Func):
    FUNC = "BETWEEN"

    def __init__(self, attribute: Operand, start: Operand, end: Operand):
        self.attr = attribute
        self.start = start
        self.end = end
        super().__init__(attribute, start, end)

    def expr(self, item: Optional[Item]) -> bool:
        # In python3 None is not a valid comparator when using < or > so must be handled specially
        start = self.start.expr(item)
        attr = self.attr.expr(item)
        end = self.end.expr(item)
        # Need to verify whether start has a valid value
        # Can't just check  'if start', because start could be 0, which is a valid integer
        start_has_value = start is not None and (isinstance(start, int) or start)
        end_has_value = end is not None and (isinstance(end, int) or end)
        if start_has_value and attr and end_has_value:
            return start <= attr <= end
        elif start is None and attr is None:
            # None is between None and None as well as None is between None and any number
            return True
        elif start is None and attr and end:
            return attr <= end
        else:
            return False


class FuncIn(Func):
    FUNC = "IN"

    def __init__(self, attribute: Operand, *possible_values: Any):
        self.attr = attribute
        self.possible_values = possible_values
        super().__init__(attribute, *possible_values)

    def expr(self, item: Optional[Item]) -> bool:
        for possible_value in self.possible_values:
            if self.attr.expr(item) == possible_value.expr(item):
                return True

        return False


COMPARATOR_CLASS = {
    "<": OpLessThan,
    ">": OpGreaterThan,
    "<=": OpLessThanOrEqual,
    ">=": OpGreaterThanOrEqual,
    "=": OpEqual,
    "<>": OpNotEqual,
}

FUNC_CLASS: Dict[str, Any] = {
    "attribute_exists": FuncAttrExists,
    "attribute_not_exists": FuncAttrNotExists,
    "attribute_type": FuncAttrType,
    "begins_with": FuncBeginsWith,
    "contains": FuncContains,
    "size": FuncSize,
    "between": FuncBetween,
}
