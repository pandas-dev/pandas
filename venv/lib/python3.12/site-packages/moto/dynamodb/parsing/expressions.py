# type: ignore
import abc
import logging
from abc import abstractmethod
from collections import deque

from moto.dynamodb.exceptions import InvalidTokenException, InvalidUpdateExpression
from moto.dynamodb.parsing.ast_nodes import (
    ExpressionAttribute,
    ExpressionAttributeName,
    ExpressionAttributeValue,
    ExpressionPathDescender,
    ExpressionSelector,
    ExpressionValueOperator,
    UpdateExpression,
    UpdateExpressionAddAction,
    UpdateExpressionAddActions,
    UpdateExpressionAddClause,
    UpdateExpressionDeleteAction,
    UpdateExpressionDeleteActions,
    UpdateExpressionDeleteClause,
    UpdateExpressionFunction,
    UpdateExpressionGroupedValue,
    UpdateExpressionPath,
    UpdateExpressionRemoveAction,
    UpdateExpressionRemoveActions,
    UpdateExpressionRemoveClause,
    UpdateExpressionSetAction,
    UpdateExpressionSetActions,
    UpdateExpressionSetClause,
    UpdateExpressionValue,
)
from moto.dynamodb.parsing.tokens import ExpressionTokenizer, Token

logger = logging.getLogger(__name__)


class NestableExpressionParserMixin:
    """
    For nodes that can be nested in themselves (recursive). Take for example UpdateExpression's grammar:

    UpdateExpression => UpdateExpressionClause*
    UpdateExpression => UpdateExpressionClause* UpdateExpression

    If we consider it of structure
    NestableExpression => TargetClause*
    NestableExpression => TargetClause* NestableExpression

    This pattern comes back multiple times. This Mixin adds re-usability for that type of pattern.

    This approach is taken since it allows to remain the ordering of the Nodes as how the corresponding tokens where
    in the originating expression.
    """

    def __init__(self):
        self.target_clauses = deque()

    def _parse_target_clause(self, factory_class):
        """

        Args:
            factory_class: The factory for the target clause e.g.  UpdateExpressionSetClauseParser

        Returns:

        """
        pos = self.token_pos
        fc = factory_class.__class__.__name__
        logger.debug(
            f"Move token pos {pos} to continue parsing with specific factory class {fc}"
        )
        # noinspection PyProtectedMember
        ast, token_pos = factory_class(**self._initializer_args())._parse_with_pos()
        self.target_clauses.append(ast)
        logger.debug(f"Continue where previous parsing ended {token_pos}")
        self.token_pos = token_pos

    @abstractmethod
    def _initializer_args(self):
        """
        Get the arguments of the initializer. This is implemented by the calling class. See ExpressionParser for an
        example.

        Returns:
            dict: A dictionary of the initializer arguments
        """

    @classmethod
    @abstractmethod
    def _nestable_class(cls):
        """
        Get the class of the Node that will be created that would be nested. For the example in the docstring this would
        be UpdateExpression

        Returns:
            class: The class of the Nodes that will be created.
        """

    def _create_node(self):
        """
        target_clauses has the nodes in order of encountering. Go through them backwards and build the tree bottom up.

        This way left-deep-descending traversal will process nodes in order.

        Continuing the example of an UpdateExpression:
            For example SET a=3 REMOVE b
                  UpdateExpression
                  /              \
             SET a=3        UpdateExpression
                                    |
                                 REMOVE b

        self.target_clauses looks like:  ( SET a=3 >> REMOVE b )
        Returns:
            moto.dynamodb.ast_nodes.Node: Node of an AST representing the Expression as produced by the factory.
        """
        cn = self.__class__.__name__
        assert len(self.target_clauses) > 0, f"No nodes for {cn}"
        target_node = self._nestable_class()(children=[self.target_clauses.pop()])
        while len(self.target_clauses) > 0:
            target_node = self._nestable_class()(
                children=[self.target_clauses.pop(), target_node]
            )
        return target_node


class ExpressionParser(metaclass=abc.ABCMeta):
    """Abstract class"""

    def __init__(self, expression_token_list, token_pos=0):
        """

        Args:
            expression_token_list:
            token_pos(int): Location where parsing is
        """
        self.token_list = expression_token_list
        self.token_pos = token_pos

    def _initializer_args(self):
        return {"expression_token_list": self.token_list, "token_pos": self.token_pos}

    @abstractmethod
    def _parse(self):
        """
        Start parsing the token_list from token_pos for the factory type.

        Returns:
            moto.dynamodb.ast_nodes.Node: AST which is root node of resulting abstract syntax tree
        """

    @classmethod
    def is_possible_start(cls, token):
        return token is not None and cls._is_possible_start(token)

    @classmethod
    @abstractmethod
    def _is_possible_start(cls, token):
        """

        Args:
            token(moto.dynamodb.tokens.Token):

        Returns:
            bool: True if token is a possible start for entries processed by `cls`
        """

    def _parse_with_pos(self):
        """
        Start parsing the token_list from token_pos for the factory type and also return the resulting token_pos.

        Returns:
            (ast, token_pos): tuple of AST which is root node of resulting abstract syntax tree and token_pos is the
                              position in the tokenlist.
        """
        return self._parse(), self.token_pos

    def parse(self):
        return self._parse()

    def get_next_token_type(self):
        """
        Get the type of the next token to be processed

        Returns:
            str: Token type or None if no more next token
        """
        try:
            return self.get_next_token().type
        except AttributeError:
            return None

    def get_next_token(self):
        """
        Get the next token to be processed

        Returns:
            moto.dynamodb.tokens.Token: or None if no more next token
        """
        try:
            return self.token_list[self.token_pos]
        except IndexError:
            return None

    def get_next_token_value(self):
        """
        Get the value of the next token to be processed

        Returns:
            str: value or None if no more next token
        """
        try:
            return self.get_next_token().value
        except AttributeError:
            return None

    def is_at_end(self):
        """Return boolean indicating whether we are at end of the parsing"""
        return self.token_pos == len(self.token_list)

    def is_at_start(self):
        """Return boolean indicating whether we are at start of the parsing"""
        return self.token_pos == 0

    def get_last_token_value(self):
        """Get the last token that was correctly parsed or return empty string"""
        if self.token_pos > 0:
            return self.token_list[self.token_pos - 1].value
        else:
            return ""

    def get_last_token_type(self):
        """Get the last token type that was correctly parsed or return None"""
        if self.token_pos > 0:
            return self.token_list[self.token_pos - 1].type
        else:
            return None

    def get_2nd_last_token_value_if_last_was_whitespace(self):
        """Get the 2nd last token that was correctly parsed if last one was whitespace or return empty string"""
        if self.token_pos > 1 and self.get_last_token_type() == Token.WHITESPACE:
            return self.token_list[self.token_pos - 2].value
        else:
            return ""

    def get_following_token_value(self):
        """Get the token value after the one that is being parsed or empty string if non existent."""
        try:
            return self.token_list[self.token_pos + 1].value
        except IndexError:
            return ""

    def get_following_token_type(self):
        """Get the token type after the one that is being parsed or None if non existent."""
        try:
            return self.token_list[self.token_pos + 1].type
        except IndexError:
            return None

    def get_2nd_following_token_value_if_following_was_whitespace(self):
        """Get the 2nd following token that was correctly parsed if 1st one was whitespace or return empty string"""
        if self.get_following_token_type() == Token.WHITESPACE:
            try:
                return self.token_list[self.token_pos + 2].value
            except IndexError:
                return ""
        else:
            return ""

    def skip_white_space(self):
        try:
            while self.get_next_token_type() == Token.WHITESPACE:
                self.token_pos += 1
        except IndexError:
            assert self.token_pos > 0, "We should always have positive indexes"
            logger.debug("We are out of range so end is reached")

    def process_token_of_type(self, token_type):
        """
        Maker sure the next token is of type `token_type` if not raise unexpected token
        Args:
            token_type: A token type

        Returns:
            str: The value if the token is of type `token_type`
        """
        if self.get_next_token_type() == token_type:
            token_value = self.get_next_token_value()
            self.goto_next_significant_token()
            return token_value
        else:
            self.raise_unexpected_token()

    def goto_next_significant_token(self):
        """Continue past current token and skip all whitespaces"""
        self.token_pos += 1
        self.skip_white_space()

    def raise_unexpected_token(self):
        if self.is_at_end():
            problematic_token = "<EOF>"
            problematic_token_in_near = ""
        else:
            problematic_token_in_near = problematic_token = self.get_next_token_value()

        near = "".join(
            [
                self.get_2nd_last_token_value_if_last_was_whitespace(),
                self.get_last_token_value(),
                problematic_token_in_near,
                self.get_following_token_value(),
                self.get_2nd_following_token_value_if_following_was_whitespace(),
            ]
        )

        raise InvalidTokenException(problematic_token, near)


class NestableBinExpressionParser(ExpressionParser):
    """
    For nodes that can be nested in themselves (recursive) but with an operation. Take for example
    UpdateExpressionValue's grammar:

    Value => Operand*
    Value => Operand* + Value
    Value => Operand* - Value

    If we consider it of structure
    NestableBinExpression => TargetClause*
    NestableBinExpression => TargetClause* BinOp NestableBinExpression

    This pattern comes back multiple times. This Mixin adds re-usability for that type of pattern.

    This approach is taken since it allows to remain the ordering of the Nodes as how the corresponding tokens where
    in the originating expression.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.target_nodes = deque()

    def _parse_target_clause(self, factory_class):
        """

        Args:
            factory_class: The factory for the target clause e.g.  UpdateExpressionSetClauseParser

        Returns:

        """
        # noinspection PyProtectedMember
        ast, self.token_pos = factory_class(
            **self._initializer_args()
        )._parse_with_pos()
        self.target_nodes.append(ast)
        logger.debug(f"Continue where previous parsing ended {self.token_pos}")

    def _parse(self):
        self._parse_target_clause(self._operand_factory_class())
        while self._binop_factory_class().is_possible_start(self.get_next_token()):
            self._parse_target_clause(self._binop_factory_class())
            if self._operand_factory_class().is_possible_start(self.get_next_token()):
                self._parse_target_clause(self._operand_factory_class())
            else:
                self.raise_unexpected_token()
        return self._create_node()

    @abstractmethod
    def _operand_factory_class(self):
        """
        Get the Parser class of the Operands for the Binary operations/actions.

        Returns:
            class:
        """

    @abstractmethod
    def _binop_factory_class(self):
        """
        Get a factory that gets the possible binary operation.

        Returns:
            class: A class extending ExpressionParser
        """

    def _create_node(self):
        """
        target_clauses has the nodes in order of encountering. Go through them forward and build the tree bottom up.
        For simplicity docstring will use Operand Node rather than the specific node

        This way left-deep-descending traversal will process nodes in order.

        Continuing the example of an UpdateExpressionValue:
            For example value => a + :val - :val2
                                  UpdateExpressionValue
                                  /            |      \
                    UpdateExpressionValue     BinOp   Operand
                           /    |      |       |        |
        UpdateExpressionValue  BinOp Operand   -       :val2
               /                |       |
            Operand             +      :val
              |
              a

        self.target_nodes looks like: (  a >> + >> :val >> - >> :val2 )
        Returns:
            moto.dynamodb.ast_nodes.Node: Node of an AST representing the Expression as produced by the factory.
        """
        if len(self.target_nodes) == 1:
            return UpdateExpressionValue(children=[self.target_nodes.popleft()])
        else:
            target_node = UpdateExpressionValue(
                children=[
                    self.target_nodes.popleft(),
                    self.target_nodes.popleft(),
                    self.target_nodes.popleft(),
                ]
            )
            while len(self.target_nodes) >= 2:
                target_node = UpdateExpressionValue(
                    children=[
                        target_node,
                        self.target_nodes.popleft(),
                        self.target_nodes.popleft(),
                    ]
                )
            assert len(self.target_nodes) == 0
            return target_node


class UpdateExpressionParser(ExpressionParser, NestableExpressionParserMixin):
    """
    Parser to create update expressions
    """

    @classmethod
    def _sub_factories(cls):
        return [
            UpdateExpressionSetClauseParser,
            UpdateExpressionAddClauseParser,
            UpdateExpressionDeleteClauseParser,
            UpdateExpressionRemoveClauseParser,
        ]

    @classmethod
    def _is_possible_start(cls, token):
        pass

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        NestableExpressionParserMixin.__init__(self)

    @classmethod
    def _nestable_class(cls):
        return UpdateExpression

    def _parse_expression_clause(self, factory_class):
        return self._parse_target_clause(factory_class)

    def _parse_by_a_subfactory(self):
        for sub_factory in self._sub_factories():
            if sub_factory.is_possible_start(self.get_next_token()):
                self._parse_expression_clause(sub_factory)
                return True
        return False

    def _parse(self):
        """
        Update Expression is the top-most node therefore it is expected to end up at the end of the expression.
        """
        while True:
            self.skip_white_space()
            if self.is_at_end():
                logger.debug("End reached")
                break
            elif self._parse_by_a_subfactory():
                continue
            else:
                self.raise_unexpected_token()

        return self._create_node()

    @classmethod
    def make(cls, expression_str) -> UpdateExpression:
        token_list = ExpressionTokenizer.make_list(expression_str)
        return cls(token_list).parse()


class UpdateExpressionSetClauseParser(ExpressionParser):
    """
    UpdateExpressionSetClause => SET SetActions
    """

    @classmethod
    def _is_possible_start(cls, token):
        return token.type == Token.ATTRIBUTE and token.value.upper() == "SET"

    def _parse(self):
        assert self.is_possible_start(self.get_next_token())
        self.goto_next_significant_token()
        ast, self.token_pos = UpdateExpressionSetActionsParser(
            **self._initializer_args()
        )._parse_with_pos()
        # noinspection PyProtectedMember
        return UpdateExpressionSetClause(children=[ast])


class UpdateExpressionActionsParser(ExpressionParser, NestableExpressionParserMixin):
    """
    UpdateExpressionSetActions
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        NestableExpressionParserMixin.__init__(self)

    @classmethod
    def _is_possible_start(cls, token):
        cn = cls._nestable_class().__name__
        raise RuntimeError(f"{cn} cannot be identified by the next token.")

    @classmethod
    @abstractmethod
    def _nestable_class(cls):
        return UpdateExpressionSetActions

    @classmethod
    @abstractmethod
    def _nested_expression_parser_class(cls):
        """Returns the parser for the query part that creates the nested nodes"""

    def _parse(self):
        """
        UpdateExpressionSetActions is inside the expression so it can be followed by others. Process SetActions one by
        one until no more SetAction.
        """
        self.skip_white_space()

        while self._nested_expression_parser_class().is_possible_start(
            self.get_next_token()
        ):
            self._parse_target_clause(self._nested_expression_parser_class())
            self.skip_white_space()
            if self.get_next_token_type() == Token.COMMA:
                self.goto_next_significant_token()
                if self.is_at_end():
                    # The expression should not end with a comma
                    self.raise_unexpected_token()
            else:
                break

        if len(self.target_clauses) == 0:
            nc = self._nestable_class().__name__
            nepc = self._nested_expression_parser_class().__name__
            logger.debug(f"Didn't encounter a single {nc} in {nepc}.")
            self.raise_unexpected_token()

        return self._create_node()


class UpdateExpressionSetActionsParser(UpdateExpressionActionsParser):
    """
    UpdateExpressionSetActions
    """

    @classmethod
    def _nested_expression_parser_class(cls):
        return UpdateExpressionSetActionParser

    @classmethod
    def _nestable_class(cls):
        return UpdateExpressionSetActions


class UpdateExpressionSetActionParser(ExpressionParser):
    """
    SetAction => Path = Value

    So we create an UpdateExpressionSetAction Node that has 2 children. Left child Path and right child Value.
    """

    @classmethod
    def _is_possible_start(cls, token):
        return UpdateExpressionPathParser.is_possible_start(token)

    def _parse(self):
        """
        UpdateExpressionSetActionParser only gets called when expecting a SetAction. So we should be aggressive on
        raising invalid Tokens.  We can thus do the following:
          1) Process path
          2) skip whitespace if there are any
          3) Process equal-sign token
          4) skip whitespace if there are any
          3) Process value

        """
        path, self.token_pos = UpdateExpressionPathParser(
            **self._initializer_args()
        )._parse_with_pos()
        self.skip_white_space()
        self.process_token_of_type(Token.EQUAL_SIGN)
        self.skip_white_space()
        value, self.token_pos = UpdateExpressionValueParser(
            **self._initializer_args()
        )._parse_with_pos()
        return UpdateExpressionSetAction(children=[path, value])


class UpdateExpressionPathParser(ExpressionParser):
    """
    Paths are selectors within items to specify a part within an Item. DynamoDB does not impose much restrictions on the
    data it stores but it does store more strict restrictions on how they are represented in UpdateExpression's.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.path_nodes = []

    @classmethod
    def _is_possible_start(cls, token):
        """
        Args:
            token(Token): the token to be checked

        Returns:
            bool: Whether the token could be the start of an UpdateExpressionPath
        """
        if token.type == Token.ATTRIBUTE_NAME:
            return True
        elif token.type == Token.ATTRIBUTE and token.value.upper() != "REMOVE":
            """We have to make sure remove is not passed"""
            return True
        return False

    def _parse(self):
        return self.process_path()

    def process_path(self):
        self.parse_path()
        return UpdateExpressionPath(children=self.path_nodes)

    def parse_path(self):
        """
        A path is comprised of:
          - Attribute: the name of an attribute as how it is stored which has no special characters
          - ATTRIBUTE_NAME: A placeholder that has no special characters except leading # to refer to attributes that
                            have a name that is not allowed in an UpdateExpression)
          - DOT's: These are used to decent in a nested structure. When a DOT is in a path expression it is never part
                   of an attribute name but always means to descent into a MAP. We will call each descend a patch
                   chain
          - SELECTORs: E.g.: [1] These are used to select an element in ordered datatypes like a list.

        Whitespaces can be between all these elements that build a path. For SELECTORs it is also allowed to have
        whitespaces between brackets and numbers but the number cannot be split up with spaces

        Attributes and attribute_names must be separated with DOT's.
        Returns:
            UpdateExpressionPath:
        """
        self.parse_path_chain()
        while self.is_next_token_start_of_patch_chain():
            self.process_dot()
            self.parse_path_chain()

    def is_next_token_start_of_patch_chain(self):
        return self.get_next_token_type() == Token.DOT

    def process_dot(self):
        self.path_nodes.append(ExpressionPathDescender())
        self.goto_next_significant_token()

    def parse_path_chain(self):
        self.process_attribute_identifying_token()
        self.skip_white_space()
        while self.is_next_token_start_of_selector():
            self.process_selector()
            self.skip_white_space()

    def process_attribute_identifying_token(self):
        if self.get_next_token_type() == Token.ATTRIBUTE:
            self.path_nodes.append(ExpressionAttribute(self.get_next_token_value()))
        elif self.get_next_token_type() == Token.ATTRIBUTE_NAME:
            self.path_nodes.append(ExpressionAttributeName(self.get_next_token_value()))
        else:
            self.raise_unexpected_token()

        self.goto_next_significant_token()

    def is_next_token_start_of_selector(self):
        return self.get_next_token_type() == Token.OPEN_SQUARE_BRACKET

    def process_selector(self):
        """
        Process the selector is only called when a selector must be processed. So do the following actions:
         - skip opening bracket
         - skip optional spaces
         - read numeric literal
         - skip optional spaces
         - pass closing bracket
        """
        self.process_token_of_type(Token.OPEN_SQUARE_BRACKET)
        selector_value = self.process_token_of_type(Token.NUMBER)
        self.process_token_of_type(Token.CLOSE_SQUARE_BRACKET)
        self.path_nodes.append(ExpressionSelector(selector_value))


class UpdateExpressionValueParser(NestableBinExpressionParser):
    @classmethod
    def _is_possible_start(cls, token):
        return UpdateExpressionOperandParser.is_possible_start(token)

    def _operand_factory_class(self):
        return UpdateExpressionOperandParser

    def _binop_factory_class(self):
        return UpdateExpressionValueOperatorParser


class UpdateExpressionGroupedValueParser(ExpressionParser):
    """
    A grouped value is an Update Expression value clause that is surrounded by round brackets. Each Operand can be
    a grouped value by itself.
    """

    def _parse(self):
        self.process_token_of_type(Token.OPEN_ROUND_BRACKET)
        value, self.token_pos = UpdateExpressionValueParser(
            **self._initializer_args()
        )._parse_with_pos()
        self.process_token_of_type(Token.CLOSE_ROUND_BRACKET)
        return UpdateExpressionGroupedValue(children=value)

    @classmethod
    def _is_possible_start(cls, token):
        return token.type == Token.OPEN_ROUND_BRACKET


class UpdateExpressionValueOperatorParser(ExpressionParser):
    OPERATION_TOKENS = [Token.PLUS_SIGN, Token.MINUS_SIGN]

    @classmethod
    def _is_possible_start(cls, token):
        return token.type in cls.OPERATION_TOKENS

    def _parse(self):
        operation_value = self.get_next_token_value()
        assert operation_value in self.OPERATION_TOKENS
        self.goto_next_significant_token()
        return ExpressionValueOperator(operation_value)


class UpdateExpressionOperandParser(ExpressionParser):
    """
    Grammar
    Operand* => AttributeValue
    Operand* => UpdateExpressionFunction
    Operand* => Path
    Operand* => GroupedValue
    """

    @classmethod
    def _sub_factories(cls):
        return [
            UpdateExpressionAttributeValueParser,
            UpdateExpressionFunctionParser,
            UpdateExpressionPathParser,
            UpdateExpressionGroupedValueParser,
        ]

    @classmethod
    def _is_possible_start(cls, token):
        return any(parser.is_possible_start(token) for parser in cls._sub_factories())

    def _parse(self):
        for factory in self._sub_factories():
            if factory.is_possible_start(self.get_next_token()):
                node, self.token_pos = factory(
                    **self._initializer_args()
                )._parse_with_pos()
                return node
        self.raise_unexpected_token()


class UpdateExpressionAttributeValueParser(ExpressionParser):
    def _parse(self):
        attr_value = ExpressionAttributeValue(
            self.process_token_of_type(Token.ATTRIBUTE_VALUE)
        )
        return attr_value

    @classmethod
    def _is_possible_start(cls, token):
        return token.type == Token.ATTRIBUTE_VALUE


class UpdateExpressionAttributeValueOrPathParser(ExpressionParser):
    def _parse(self):
        if UpdateExpressionAttributeValueParser.is_possible_start(
            self.get_next_token()
        ):
            token, self.token_pos = UpdateExpressionAttributeValueParser(
                **self._initializer_args()
            )._parse_with_pos()
        else:
            token, self.token_pos = UpdateExpressionPathParser(
                **self._initializer_args()
            )._parse_with_pos()
        return token

    @classmethod
    def _is_possible_start(cls, token):
        return any(
            [
                UpdateExpressionAttributeValueParser.is_possible_start(token),
                UpdateExpressionPathParser.is_possible_start(token),
            ]
        )


class UpdateExpressionFunctionParser(ExpressionParser):
    """
    A helper to process a function of an Update Expression
    """

    # Map function to the factories for its elements
    FUNCTIONS = {
        "if_not_exists": [
            UpdateExpressionPathParser,
            UpdateExpressionAttributeValueOrPathParser,
        ],
        "list_append": [UpdateExpressionOperandParser, UpdateExpressionOperandParser],
    }

    @classmethod
    def _is_possible_start(cls, token):
        """
        Check whether a token is supposed to be a function
        Args:
            token(Token): the token to check

        Returns:
            bool: True if token is the start of a function.
        """
        if token.type == Token.ATTRIBUTE:
            return token.value in cls.FUNCTIONS.keys()
        else:
            return False

    def _parse(self):
        function_name = self.get_next_token_value()
        if function_name not in self.FUNCTIONS.keys():
            # Function names are case sensitive
            raise InvalidUpdateExpression(function_name)
        self.goto_next_significant_token()
        self.process_token_of_type(Token.OPEN_ROUND_BRACKET)
        function_elements = [function_name]
        function_arguments = self.FUNCTIONS[function_name]
        for i, func_elem_factory in enumerate(function_arguments):
            func_elem, self.token_pos = func_elem_factory(
                **self._initializer_args()
            )._parse_with_pos()
            function_elements.append(func_elem)
            if i + 1 < len(function_arguments):
                self.skip_white_space()
                self.process_token_of_type(Token.COMMA)
        self.process_token_of_type(Token.CLOSE_ROUND_BRACKET)
        return UpdateExpressionFunction(children=function_elements)


class UpdateExpressionRemoveClauseParser(ExpressionParser):
    """
    UpdateExpressionRemoveClause => REMOVE RemoveActions
    """

    def _parse(self):
        assert self.is_possible_start(self.get_next_token())
        self.goto_next_significant_token()
        ast, self.token_pos = UpdateExpressionRemoveActionsParser(
            **self._initializer_args()
        )._parse_with_pos()
        # noinspection PyProtectedMember
        return UpdateExpressionRemoveClause(children=[ast])

    @classmethod
    def _is_possible_start(cls, token):
        """REMOVE is not a keyword"""
        return token.type == Token.ATTRIBUTE and token.value.upper() == "REMOVE"


class UpdateExpressionRemoveActionsParser(UpdateExpressionActionsParser):
    """
    UpdateExpressionSetActions
    """

    @classmethod
    def _nested_expression_parser_class(cls):
        return UpdateExpressionRemoveActionParser

    @classmethod
    def _nestable_class(cls):
        return UpdateExpressionRemoveActions


class UpdateExpressionRemoveActionParser(ExpressionParser):
    """
    RemoveAction => Path = Value

    So we create an UpdateExpressionSetAction Node that has 2 children. Left child Path and right child Value.
    """

    @classmethod
    def _is_possible_start(cls, token):
        return UpdateExpressionPathParser.is_possible_start(token)

    def _parse(self):
        """
        UpdateExpressionRemoveActionParser only gets called when expecting a RemoveAction. So we should be aggressive on
        raising invalid Tokens.  We can thus do the following:
          1) Process path
          2) skip whitespace if there are any

        """
        path, self.token_pos = UpdateExpressionPathParser(
            **self._initializer_args()
        )._parse_with_pos()
        self.skip_white_space()
        return UpdateExpressionRemoveAction(children=[path])


class UpdateExpressionAddClauseParser(ExpressionParser):
    def _parse(self):
        assert self.is_possible_start(self.get_next_token())
        self.goto_next_significant_token()
        ast, self.token_pos = UpdateExpressionAddActionsParser(
            **self._initializer_args()
        )._parse_with_pos()
        # noinspection PyProtectedMember
        return UpdateExpressionAddClause(children=[ast])

    @classmethod
    def _is_possible_start(cls, token):
        return token.type == Token.ATTRIBUTE and token.value.upper() == "ADD"


class UpdateExpressionAddActionsParser(UpdateExpressionActionsParser):
    """
    UpdateExpressionSetActions
    """

    @classmethod
    def _nested_expression_parser_class(cls):
        return UpdateExpressionAddActionParser

    @classmethod
    def _nestable_class(cls):
        return UpdateExpressionAddActions


class UpdateExpressionPathValueParser(ExpressionParser, metaclass=abc.ABCMeta):
    def _parse_path_and_value(self):
        """
        UpdateExpressionAddActionParser only gets called when expecting an AddAction. So we should be aggressive on
        raising invalid Tokens.  We can thus do the following:
          1) Process path
          2) skip whitespace if there are any
          3) Process a value
          4) skip whitespace if there are any

        Returns:
            [path, value]: A list containing the Path node and the AttributeValue nodes
        """
        path, self.token_pos = UpdateExpressionPathParser(
            **self._initializer_args()
        )._parse_with_pos()
        self.skip_white_space()
        value, self.token_pos = UpdateExpressionAttributeValueParser(
            **self._initializer_args()
        )._parse_with_pos()
        self.skip_white_space()
        return [path, value]


class UpdateExpressionAddActionParser(UpdateExpressionPathValueParser):
    @classmethod
    def _is_possible_start(cls, token):
        return UpdateExpressionPathParser.is_possible_start(token)

    def _parse(self):
        return UpdateExpressionAddAction(children=self._parse_path_and_value())


class UpdateExpressionDeleteClauseParser(ExpressionParser):
    def _parse(self):
        assert self.is_possible_start(self.get_next_token())
        self.goto_next_significant_token()
        ast, self.token_pos = UpdateExpressionDeleteActionsParser(
            **self._initializer_args()
        )._parse_with_pos()
        # noinspection PyProtectedMember
        return UpdateExpressionDeleteClause(children=[ast])

    @classmethod
    def _is_possible_start(cls, token):
        return token.type == Token.ATTRIBUTE and token.value.upper() == "DELETE"


class UpdateExpressionDeleteActionsParser(UpdateExpressionActionsParser):
    """
    UpdateExpressionSetActions
    """

    @classmethod
    def _nested_expression_parser_class(cls):
        return UpdateExpressionDeleteActionParser

    @classmethod
    def _nestable_class(cls):
        return UpdateExpressionDeleteActions


class UpdateExpressionDeleteActionParser(UpdateExpressionPathValueParser):
    @classmethod
    def _is_possible_start(cls, token):
        return UpdateExpressionPathParser.is_possible_start(token)

    def _parse(self):
        return UpdateExpressionDeleteAction(children=self._parse_path_and_value())
