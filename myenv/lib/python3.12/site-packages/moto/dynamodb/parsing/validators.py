"""
See docstring class Validator below for more details on validation
"""

from abc import abstractmethod
from copy import deepcopy
from typing import Any, Callable, Dict, List, Type, Union

from moto.dynamodb.exceptions import (
    AttributeDoesNotExist,
    AttributeIsReservedKeyword,
    EmptyKeyAttributeException,
    ExpressionAttributeNameNotDefined,
    ExpressionAttributeValueNotDefined,
    IncorrectOperandType,
    InvalidAttributeTypeError,
    InvalidUpdateExpressionInvalidDocumentPath,
    MockValidationException,
    ProvidedKeyDoesNotExist,
    UpdateHashRangeKeyException,
)
from moto.dynamodb.models.dynamo_type import DynamoType, Item
from moto.dynamodb.models.table import Table
from moto.dynamodb.parsing.ast_nodes import (  # type: ignore
    DDBTypedValue,
    DepthFirstTraverser,
    ExpressionAttribute,
    ExpressionAttributeName,
    ExpressionAttributeValue,
    ExpressionPathDescender,
    ExpressionSelector,
    ExpressionValueOperator,
    Node,
    NoneExistingPath,
    UpdateExpressionAddAction,
    UpdateExpressionClause,
    UpdateExpressionDeleteAction,
    UpdateExpressionFunction,
    UpdateExpressionPath,
    UpdateExpressionRemoveAction,
    UpdateExpressionSetAction,
    UpdateExpressionValue,
)
from moto.dynamodb.parsing.reserved_keywords import ReservedKeywords


class ExpressionAttributeValueProcessor(DepthFirstTraverser):  # type: ignore[misc]
    def __init__(self, expression_attribute_values: Dict[str, Dict[str, Any]]):
        self.expression_attribute_values = expression_attribute_values

    def _processing_map(
        self,
    ) -> Dict[
        Type[ExpressionAttributeValue],
        Callable[[ExpressionAttributeValue], DDBTypedValue],
    ]:
        return {
            ExpressionAttributeValue: self.replace_expression_attribute_value_with_value
        }

    def replace_expression_attribute_value_with_value(
        self, node: ExpressionAttributeValue
    ) -> DDBTypedValue:
        """A node representing an Expression Attribute Value. Resolve and replace value"""
        assert isinstance(node, ExpressionAttributeValue)
        attribute_value_name = node.get_value_name()
        try:
            target = self.expression_attribute_values[attribute_value_name]
        except KeyError:
            raise ExpressionAttributeValueNotDefined(
                attribute_value=attribute_value_name
            )
        return DDBTypedValue(DynamoType(target))


class ExpressionPathResolver:
    def __init__(self, expression_attribute_names: Dict[str, str]):
        self.expression_attribute_names = expression_attribute_names

    @classmethod
    def raise_exception_if_keyword(cls, attribute: Any) -> None:  # type: ignore[misc]
        if attribute.upper() in ReservedKeywords.get_reserved_keywords():
            raise AttributeIsReservedKeyword(attribute)

    def resolve_expression_path(
        self, item: Item, update_expression_path: UpdateExpressionPath
    ) -> Union[NoneExistingPath, DDBTypedValue]:
        assert isinstance(update_expression_path, UpdateExpressionPath)
        return self.resolve_expression_path_nodes(item, update_expression_path.children)

    def resolve_expression_path_nodes(
        self, item: Item, update_expression_path_nodes: List[Node]
    ) -> Union[NoneExistingPath, DDBTypedValue]:
        target = item.attrs

        for child in update_expression_path_nodes:
            # First replace placeholder with attribute_name
            attr_name = None
            if isinstance(child, ExpressionAttributeName):
                attr_placeholder = child.get_attribute_name_placeholder()
                try:
                    attr_name = self.expression_attribute_names[attr_placeholder]
                except KeyError:
                    raise ExpressionAttributeNameNotDefined(attr_placeholder)
            elif isinstance(child, ExpressionAttribute):
                attr_name = child.get_attribute_name()
                self.raise_exception_if_keyword(attr_name)
            if attr_name is not None:
                # Resolv attribute_name
                try:
                    target = target[attr_name]
                except (KeyError, TypeError):
                    if child == update_expression_path_nodes[-1]:
                        return NoneExistingPath(creatable=True)
                    return NoneExistingPath()
            else:
                if isinstance(child, ExpressionPathDescender):
                    continue
                elif isinstance(child, ExpressionSelector):
                    index = child.get_index()
                    if target.is_list():  # type: ignore
                        try:
                            target = target[index]
                        except IndexError:
                            # When a list goes out of bounds when assigning that is no problem when at the assignment
                            # side. It will just append to the list.
                            if child == update_expression_path_nodes[-1]:
                                return NoneExistingPath(creatable=True)
                            return NoneExistingPath()
                    else:
                        raise InvalidUpdateExpressionInvalidDocumentPath
                else:
                    raise NotImplementedError(f"Path resolution for {type(child)}")
        return DDBTypedValue(target)

    def resolve_expression_path_nodes_to_dynamo_type(
        self, item: Item, update_expression_path_nodes: List[Node]
    ) -> Any:
        node = self.resolve_expression_path_nodes(item, update_expression_path_nodes)
        if isinstance(node, NoneExistingPath):
            raise ProvidedKeyDoesNotExist()
        assert isinstance(node, DDBTypedValue)
        return node.get_value()


class ExpressionAttributeResolvingProcessor(DepthFirstTraverser):  # type: ignore[misc]
    def _processing_map(
        self,
    ) -> Dict[Type[UpdateExpressionClause], Callable[[DDBTypedValue], DDBTypedValue]]:
        return {
            UpdateExpressionSetAction: self.disable_resolving,
            UpdateExpressionPath: self.process_expression_path_node,
        }

    def __init__(self, expression_attribute_names: Dict[str, str], item: Item):
        self.expression_attribute_names = expression_attribute_names
        self.item = item
        self.resolving = False

    def pre_processing_of_child(
        self,
        parent_node: Union[
            UpdateExpressionSetAction,
            UpdateExpressionRemoveAction,
            UpdateExpressionDeleteAction,
            UpdateExpressionAddAction,
        ],
        child_id: int,
    ) -> None:
        """
        We have to enable resolving if we are processing a child of UpdateExpressionSetAction that is not first.
        Because first argument is path to be set, 2nd argument would be the value.
        """
        if isinstance(
            parent_node,
            (
                UpdateExpressionSetAction,
                UpdateExpressionRemoveAction,
                UpdateExpressionDeleteAction,
                UpdateExpressionAddAction,
            ),
        ):
            if child_id == 0:
                self.resolving = False
            else:
                self.resolving = True

    def disable_resolving(self, node: DDBTypedValue) -> DDBTypedValue:
        self.resolving = False
        return node

    def process_expression_path_node(self, node: DDBTypedValue) -> DDBTypedValue:
        """Resolve ExpressionAttribute if not part of a path and resolving is enabled."""
        if self.resolving:
            return self.resolve_expression_path(node)
        else:
            # Still resolve but return original note to make sure path is correct Just make sure nodes are creatable.
            result_node = self.resolve_expression_path(node)
            if (
                isinstance(result_node, NoneExistingPath)
                and not result_node.is_creatable()
            ):
                raise InvalidUpdateExpressionInvalidDocumentPath()

            return node

    def resolve_expression_path(
        self, node: DDBTypedValue
    ) -> Union[NoneExistingPath, DDBTypedValue]:
        return ExpressionPathResolver(
            self.expression_attribute_names
        ).resolve_expression_path(self.item, node)


class UpdateExpressionFunctionEvaluator(DepthFirstTraverser):  # type: ignore[misc]
    """
    At time of writing there are only 2 functions for DDB UpdateExpressions. They both are specific to the SET
    expression as per the official AWS docs:
        https://docs.aws.amazon.com/amazondynamodb/latest/developerguide/
        Expressions.UpdateExpressions.html#Expressions.UpdateExpressions.SET
    """

    def _processing_map(
        self,
    ) -> Dict[
        Type[UpdateExpressionFunction],
        Callable[[UpdateExpressionFunction], DDBTypedValue],
    ]:
        return {UpdateExpressionFunction: self.process_function}

    def process_function(self, node: UpdateExpressionFunction) -> DDBTypedValue:
        assert isinstance(node, UpdateExpressionFunction)
        function_name = node.get_function_name()
        first_arg = node.get_nth_argument(1)
        second_arg = node.get_nth_argument(2)

        if function_name == "if_not_exists":
            if isinstance(first_arg, NoneExistingPath):
                result = second_arg
            else:
                result = first_arg
            assert isinstance(result, (DDBTypedValue, NoneExistingPath))
            return result
        elif function_name == "list_append":
            if isinstance(first_arg, NoneExistingPath):
                raise MockValidationException(
                    "The provided expression refers to an attribute that does not exist in the item"
                )
            first_arg = deepcopy(
                self.get_list_from_ddb_typed_value(first_arg, function_name)
            )
            second_arg = self.get_list_from_ddb_typed_value(second_arg, function_name)
            for list_element in second_arg.value:
                first_arg.value.append(list_element)
            return DDBTypedValue(first_arg)
        else:
            raise NotImplementedError(f"Unsupported function for moto {function_name}")

    @classmethod
    def get_list_from_ddb_typed_value(  # type: ignore[misc]
        cls, node: DDBTypedValue, function_name: str
    ) -> DynamoType:
        assert isinstance(node, DDBTypedValue)
        dynamo_value = node.get_value()
        assert isinstance(dynamo_value, DynamoType)
        if not dynamo_value.is_list():
            raise IncorrectOperandType(function_name, dynamo_value.type)
        return dynamo_value


class NoneExistingPathChecker(DepthFirstTraverser):  # type: ignore[misc]
    """
    Pass through the AST and make sure there are no none-existing paths.
    """

    def _processing_map(self) -> Dict[Type[NoneExistingPath], Callable[[Node], None]]:
        return {NoneExistingPath: self.raise_none_existing_path}

    def raise_none_existing_path(self, node: Node) -> None:
        raise AttributeDoesNotExist


class ExecuteOperations(DepthFirstTraverser):  # type: ignore[misc]
    def _processing_map(
        self,
    ) -> Dict[Type[UpdateExpressionValue], Callable[[Node], DDBTypedValue]]:
        return {UpdateExpressionValue: self.process_update_expression_value}

    def process_update_expression_value(self, node: Node) -> DDBTypedValue:
        """
        If an UpdateExpressionValue only has a single child the node will be replaced with the childe.
        Otherwise it has 3 children and the middle one is an ExpressionValueOperator which details how to combine them
        Args:
            node(Node):

        Returns:
            Node: The resulting node of the operation if present or the child.
        """
        assert isinstance(node, UpdateExpressionValue)
        if len(node.children) == 1:
            return node.children[0]
        elif len(node.children) == 3:
            operator_node = node.children[1]
            assert isinstance(operator_node, ExpressionValueOperator)
            operator = operator_node.get_operator()
            left_operand = self.get_dynamo_value_from_ddb_typed_value(node.children[0])
            right_operand = self.get_dynamo_value_from_ddb_typed_value(node.children[2])
            if operator == "+":
                return self.get_sum(left_operand, right_operand)
            elif operator == "-":
                return self.get_subtraction(left_operand, right_operand)
            else:
                raise NotImplementedError(f"Moto does not support operator {operator}")
        else:
            raise NotImplementedError(
                "UpdateExpressionValue only has implementations for 1 or 3 children."
            )

    @classmethod
    def get_dynamo_value_from_ddb_typed_value(cls, node: DDBTypedValue) -> DynamoType:  # type: ignore[misc]
        assert isinstance(node, DDBTypedValue)
        dynamo_value = node.get_value()
        assert isinstance(dynamo_value, DynamoType)
        return dynamo_value

    @classmethod
    def get_sum(  # type: ignore[misc]
        cls, left_operand: DynamoType, right_operand: DynamoType
    ) -> DDBTypedValue:
        """
        Args:
            left_operand(DynamoType):
            right_operand(DynamoType):

        Returns:
            DDBTypedValue:
        """
        try:
            return DDBTypedValue(left_operand + right_operand)
        except TypeError:
            raise IncorrectOperandType("+", left_operand.type)

    @classmethod
    def get_subtraction(  # type: ignore[misc]
        cls, left_operand: DynamoType, right_operand: DynamoType
    ) -> DDBTypedValue:
        """
        Args:
            left_operand(DynamoType):
            right_operand(DynamoType):

        Returns:
            DDBTypedValue:
        """
        try:
            return DDBTypedValue(left_operand - right_operand)
        except TypeError:
            raise IncorrectOperandType("-", left_operand.type)


class EmptyStringKeyValueValidator(DepthFirstTraverser):  # type: ignore[misc]
    def __init__(self, key_attributes: List[str]):
        self.key_attributes = key_attributes

    def _processing_map(
        self,
    ) -> Dict[
        Type[UpdateExpressionSetAction],
        Callable[[UpdateExpressionSetAction], UpdateExpressionSetAction],
    ]:
        return {UpdateExpressionSetAction: self.check_for_empty_string_key_value}

    def check_for_empty_string_key_value(
        self, node: UpdateExpressionSetAction
    ) -> UpdateExpressionSetAction:
        """A node representing a SET action. Check that keys are not being assigned empty strings"""
        assert isinstance(node, UpdateExpressionSetAction)
        assert len(node.children) == 2
        key = node.children[0].children[0].children[0]
        val_node = node.children[1].children[0]
        if (
            not val_node.value
            and val_node.type in ["S", "B"]
            and key in self.key_attributes
        ):
            raise EmptyKeyAttributeException(key_in_index=True)
        return node


class TypeMismatchValidator(DepthFirstTraverser):  # type: ignore[misc]
    def __init__(self, key_attributes_type: List[Dict[str, str]]):
        self.key_attributes_type = key_attributes_type

    def _processing_map(
        self,
    ) -> Dict[
        Type[UpdateExpressionSetAction],
        Callable[[UpdateExpressionSetAction], UpdateExpressionSetAction],
    ]:
        return {UpdateExpressionSetAction: self.check_for_type_mismatch}

    def check_for_type_mismatch(
        self, node: UpdateExpressionSetAction
    ) -> UpdateExpressionSetAction:
        """A node representing a SET action. Check that type matches with the definition"""
        assert isinstance(node, UpdateExpressionSetAction)
        assert len(node.children) == 2
        key = node.children[0].children[0].children[0]
        val_node = node.children[1].children[0]
        for dct in self.key_attributes_type:
            if dct["AttributeName"] == key and dct["AttributeType"] != val_node.type:
                raise InvalidAttributeTypeError(
                    key, dct["AttributeType"], val_node.type
                )
        return node


class UpdateHashRangeKeyValidator(DepthFirstTraverser):  # type: ignore[misc]
    def __init__(
        self,
        table_key_attributes: List[str],
        expression_attribute_names: Dict[str, str],
    ):
        self.table_key_attributes = table_key_attributes
        self.expression_attribute_names = expression_attribute_names

    def _processing_map(
        self,
    ) -> Dict[
        Type[UpdateExpressionPath],
        Callable[[UpdateExpressionPath], UpdateExpressionPath],
    ]:
        return {UpdateExpressionPath: self.check_for_hash_or_range_key}

    def check_for_hash_or_range_key(
        self, node: UpdateExpressionPath
    ) -> UpdateExpressionPath:
        """Check that hash and range keys are not updated"""
        key_to_update = node.children[0].children[0]
        key_to_update = self.expression_attribute_names.get(
            key_to_update, key_to_update
        )
        if key_to_update in self.table_key_attributes:
            raise UpdateHashRangeKeyException(key_to_update)
        return node


class Validator:
    """
    A validator is used to validate expressions which are passed in as an AST.
    """

    def __init__(
        self,
        expression: Node,
        expression_attribute_names: Dict[str, str],
        expression_attribute_values: Dict[str, Dict[str, Any]],
        item: Item,
        table: Table,
    ):
        """
        Besides validation the Validator should also replace referenced parts of an item which is cheapest upon
        validation.

        Args:
            expression(Node): The root node of the AST representing the expression to be validated
            expression_attribute_names(ExpressionAttributeNames):
            expression_attribute_values(ExpressionAttributeValues):
            item(Item): The item which will be updated (pointed to by Key of update_item)
        """
        self.expression_attribute_names = expression_attribute_names
        self.expression_attribute_values = expression_attribute_values
        self.item = item
        self.table = table
        self.processors = self.get_ast_processors()
        self.node_to_validate = deepcopy(expression)

    @abstractmethod
    def get_ast_processors(self) -> List[DepthFirstTraverser]:  # type: ignore[misc]
        """Get the different processors that go through the AST tree and processes the nodes."""

    def validate(self) -> Node:
        n = self.node_to_validate
        for processor in self.processors:
            n = processor.traverse(n)
        return n


class UpdateExpressionValidator(Validator):
    def get_ast_processors(self) -> List[DepthFirstTraverser]:
        """Get the different processors that go through the AST tree and processes the nodes."""
        processors = [
            UpdateHashRangeKeyValidator(
                self.table.table_key_attrs, self.expression_attribute_names or {}
            ),
            ExpressionAttributeValueProcessor(self.expression_attribute_values),
            ExpressionAttributeResolvingProcessor(
                self.expression_attribute_names, self.item
            ),
            UpdateExpressionFunctionEvaluator(),
            NoneExistingPathChecker(),
            ExecuteOperations(),
            TypeMismatchValidator(self.table.attr),
            EmptyStringKeyValueValidator(self.table.attribute_keys),
        ]
        return processors
