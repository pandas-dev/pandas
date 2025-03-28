# type: ignore
import abc
from abc import abstractmethod
from collections import deque

from moto.dynamodb.models import DynamoType

from ..exceptions import (
    DuplicateUpdateExpression,
    MockValidationException,
    TooManyClauses,
)
from ..utils import extract_duplicates


class Node(metaclass=abc.ABCMeta):
    def __init__(self, children=None):
        self.type = self.__class__.__name__
        assert children is None or isinstance(children, list)
        self.children = children
        self.parent = None

        if isinstance(children, list):
            for child in children:
                if isinstance(child, Node):
                    child.set_parent(self)

    def set_parent(self, parent_node):
        self.parent = parent_node

    def validate(self, limit_set_actions: bool = False) -> None:
        if self.type == "UpdateExpression":
            add_clauses = self.find_clauses([UpdateExpressionAddClause])
            remove_clauses = self.find_clauses([UpdateExpressionRemoveClause])

            # Only allow single ADD/REMOVE clauses
            if len(add_clauses) > 1:
                raise TooManyClauses("ADD")
            if len(remove_clauses) > 1:
                raise TooManyClauses("REMOVE")

            add_actions = self.find_clauses([UpdateExpressionAddAction])
            delete_actions = self.find_clauses([UpdateExpressionDeleteAction])
            set_actions = self.find_clauses([UpdateExpressionSetAction])

            # Cannot ADD/DELETE to the same path
            add_paths = [a.get_value() for a in add_actions]
            delete_paths = [a.get_value() for a in delete_actions]
            if set(add_paths).intersection(delete_paths):
                raise DuplicateUpdateExpression(names=[*add_paths, *delete_paths])

            # Ensure SET has no duplicates
            # We currently only check for duplicates
            # We should also check for partial duplicates, i.e. [attr, attr.sub] is also invalid
            set_attributes = [s.children[0].to_str() for s in set_actions]
            duplicates = extract_duplicates(set_attributes)
            if duplicates:
                # There might be more than one attribute duplicated:
                # they may get mixed up in the Error Message which is inline with actual boto3 Error Messages
                raise DuplicateUpdateExpression(duplicates)

            set_clauses = self.find_clauses([UpdateExpressionSetClause])
            if limit_set_actions and len(set_clauses) > 1:
                raise MockValidationException(
                    'Invalid UpdateExpression: The "SET" section can only be used once in an update expression;'
                )

    def normalize(self):
        """
        Flatten the Add-/Delete-/Remove-/Set-Action children within this Node
        """
        if self.type == "UpdateExpression":
            # We can have multiple REMOVE attr[idx] expressions, such as attr[i] and attr[i+2]
            # If we remove attr[i] first, attr[i+2] suddenly refers to a different item
            # So we sort them in reverse order - we can remove attr[i+2] first, attr[i] still refers to the same item

            # Behaviour that is unknown, for now:
            # What happens if we SET and REMOVE on the same list - what takes precedence?
            # We're assuming this is executed in original order

            remove_actions = []
            sorted_actions = []
            possible_clauses = [
                UpdateExpressionAddAction,
                UpdateExpressionDeleteAction,
                UpdateExpressionRemoveAction,
                UpdateExpressionSetAction,
            ]
            for action in self.find_clauses(possible_clauses):
                if isinstance(action, UpdateExpressionRemoveAction):
                    # Keep these separate for now
                    remove_actions.append(action)
                else:
                    if len(remove_actions) > 0:
                        # Remove-actions were found earlier
                        # Now that we have other action-types, that means we've found all possible Remove-actions
                        # Sort them appropriately
                        sorted_actions.extend(sorted(remove_actions, reverse=True))
                        remove_actions.clear()
                    # Add other actions by insertion order
                    sorted_actions.append(action)
            # Remove actions were found last
            if len(remove_actions) > 0:
                sorted_actions.extend(sorted(remove_actions, reverse=True))

            self.children = sorted_actions

    def find_clauses(self, clause_types):
        clauses = []
        for child in self.children or []:
            if type(child) in clause_types:
                clauses.append(child)
            elif isinstance(child, Expression):
                clauses.extend(child.find_clauses(clause_types))
            elif isinstance(child, UpdateExpressionFunction):
                clauses.extend(child.find_clauses(clause_types))
        return clauses


class LeafNode(Node):
    """A LeafNode is a Node where none of the children are Nodes themselves."""

    def __init__(self, children=None):
        super().__init__(children)


class Expression(Node, metaclass=abc.ABCMeta):
    """
    Abstract Syntax Tree representing the expression

    For the Grammar start here and jump down into the classes at the righ-hand side to look further. Nodes marked with
    a star are abstract and won't appear in the final AST.

    Expression* => UpdateExpression
    Expression* => ConditionExpression
    """


class UpdateExpression(Expression):
    """
    UpdateExpression => UpdateExpressionClause*
    UpdateExpression => UpdateExpressionClause* UpdateExpression
    """


class UpdateExpressionClause(UpdateExpression, metaclass=abc.ABCMeta):
    """
    UpdateExpressionClause* => UpdateExpressionSetClause
    UpdateExpressionClause* => UpdateExpressionRemoveClause
    UpdateExpressionClause* => UpdateExpressionAddClause
    UpdateExpressionClause* => UpdateExpressionDeleteClause
    """


class UpdateExpressionSetClause(UpdateExpressionClause):
    """
    UpdateExpressionSetClause => SET SetActions
    """


class UpdateExpressionSetActions(UpdateExpressionClause):
    """
    UpdateExpressionSetClause => SET SetActions

    SetActions => SetAction
    SetActions => SetAction , SetActions

    """


class UpdateExpressionSetAction(UpdateExpressionClause):
    """
    SetAction => Path = Value
    """


class UpdateExpressionAction(UpdateExpressionClause):
    def get_value(self):
        expression_path = self.children[0]
        expression_selector = expression_path.children[-1]
        return expression_selector.children[0]


class UpdateExpressionRemoveActions(UpdateExpressionClause):
    """
    UpdateExpressionSetClause => REMOVE RemoveActions

    RemoveActions => RemoveAction
    RemoveActions => RemoveAction , RemoveActions
    """


class UpdateExpressionRemoveAction(UpdateExpressionAction):
    """
    RemoveAction => Path
    """

    def __lt__(self, other):
        self_value = self.get_value()
        other_value = other.get_value()
        if isinstance(self_value, int) and isinstance(other_value, int):
            return self_value < other_value
        else:
            return str(self_value) < str(other_value)


class UpdateExpressionAddActions(UpdateExpressionClause):
    """
    UpdateExpressionAddClause => ADD RemoveActions

    AddActions => AddAction
    AddActions => AddAction , AddActions
    """


class UpdateExpressionAddAction(UpdateExpressionAction):
    """
    AddAction => Path Value
    """


class UpdateExpressionDeleteActions(UpdateExpressionClause):
    """
    UpdateExpressionDeleteClause => DELETE RemoveActions

    DeleteActions => DeleteAction
    DeleteActions => DeleteAction , DeleteActions
    """


class UpdateExpressionDeleteAction(UpdateExpressionAction):
    """
    DeleteAction => Path Value
    """


class UpdateExpressionPath(UpdateExpressionClause):
    def to_str(self):
        return "".join(x.to_str() for x in self.children)


class UpdateExpressionValue(UpdateExpressionClause):
    """
    Value => Operand
    Value => Operand + Value
    Value => Operand - Value
    """


class UpdateExpressionGroupedValue(UpdateExpressionClause):
    """
    GroupedValue => ( Value )
    """


class UpdateExpressionRemoveClause(UpdateExpressionClause):
    """
    UpdateExpressionRemoveClause => REMOVE RemoveActions
    """


class UpdateExpressionAddClause(UpdateExpressionClause):
    """
    UpdateExpressionAddClause => ADD AddActions
    """


class UpdateExpressionDeleteClause(UpdateExpressionClause):
    """
    UpdateExpressionDeleteClause => DELETE DeleteActions
    """


class ExpressionPathDescender(Node):
    """Node identifying descender into nested structure (.) in expression"""

    def to_str(self):
        return "."


class ExpressionSelector(LeafNode):
    """Node identifying selector [selection_index] in expresion"""

    def __init__(self, selection_index):
        try:
            super().__init__(children=[int(selection_index)])
        except ValueError:
            assert (
                False
            ), "Expression selector must be an int, this is a bug in the moto library."

    def get_index(self):
        return self.children[0]

    def to_str(self):
        return f"[{self.get_index()}]"


class ExpressionAttribute(LeafNode):
    """An attribute identifier as used in the DDB item"""

    def __init__(self, attribute):
        super().__init__(children=[attribute])

    def get_attribute_name(self):
        return self.children[0]

    def to_str(self):
        return self.get_attribute_name()


class ExpressionAttributeName(LeafNode):
    """An ExpressionAttributeName is an alias for an attribute identifier"""

    def __init__(self, attribute_name):
        super().__init__(children=[attribute_name])

    def get_attribute_name_placeholder(self):
        return self.children[0]

    def to_str(self):
        return self.get_attribute_name_placeholder()


class ExpressionAttributeValue(LeafNode):
    """An ExpressionAttributeValue is an alias for an value"""

    def __init__(self, value):
        super().__init__(children=[value])

    def get_value_name(self):
        return self.children[0]


class ExpressionValueOperator(LeafNode):
    """An ExpressionValueOperator is an operation that works on 2 values"""

    def __init__(self, value):
        super().__init__(children=[value])

    def get_operator(self):
        return self.children[0]


class UpdateExpressionFunction(Node):
    """
    A Node representing a function of an Update Expression. The first child is the function name the others are the
    arguments.
    """

    def get_function_name(self):
        return self.children[0]

    def get_nth_argument(self, n=1):
        """Return nth element where n is a 1-based index."""
        assert n >= 1
        return self.children[n]


class DDBTypedValue(Node):
    """
    A node representing a DDBTyped value. This can be any structure as supported by DyanmoDB. The node only has 1 child
    which is the value of type `DynamoType`.
    """

    def __init__(self, value):
        assert isinstance(value, DynamoType), "DDBTypedValue must be of DynamoType"
        super().__init__(children=[value])

    def get_value(self):
        return self.children[0]


class NoneExistingPath(LeafNode):
    """A placeholder for Paths that did not exist in the Item."""

    def __init__(self, creatable=False):
        super().__init__(children=[creatable])

    def is_creatable(self):
        """Can this path be created if need be. For example path creating element in a dictionary or creating a new
        attribute under root level of an item."""
        return self.children[0]


class DepthFirstTraverser(object):
    """
    Helper class that allows depth first traversal and to implement custom processing for certain AST nodes. The
    processor of a node must return the new resulting node. This node will be placed in the tree. Processing of a
    node using this traverser should therefore only transform child nodes.  The returned node will get the same parent
    as the node before processing had.
    """

    @abstractmethod
    def _processing_map(self):
        """
        A map providing a processing function per node class type to a function that takes in a Node object and
        processes it. A Node can only be processed by a single function and they are considered in order. Therefore if
        multiple classes from a single class hierarchy strain are used the more specific classes have to be put before
        the less specific ones. That requires overriding `nodes_to_be_processed`. If no multiple classes form a single
        class hierarchy strain are used the default implementation of `nodes_to_be_processed` should be OK.
        Returns:
            dict: Mapping a Node Class to a processing function.
        """
        pass

    def nodes_to_be_processed(self):
        """Cached accessor for getting Node types that need to be processed."""
        return tuple(k for k in self._processing_map().keys())

    def process(self, node):
        """Process a Node"""
        for class_key, processor in self._processing_map().items():
            if isinstance(node, class_key):
                return processor(node)

    def pre_processing_of_child(self, parent_node, child_id):
        """Hook that is called pre-processing of the child at position `child_id`"""
        pass

    def traverse_node_recursively(self, node, child_id=-1):
        """
        Traverse nodes depth first processing nodes bottom up (if root node is considered the top).

        Args:
            node(Node): The node which is the last node to be processed but which allows to identify all the
                             work (which is in the children)
            child_id(int): The index in the list of children from the parent that this node corresponds to

        Returns:
            Node: The node of the new processed AST
        """
        if isinstance(node, Node):
            parent_node = node.parent
            if node.children is not None:
                for i, child_node in enumerate(node.children):
                    self.pre_processing_of_child(node, i)
                    self.traverse_node_recursively(child_node, i)
            # noinspection PyTypeChecker
            if isinstance(node, self.nodes_to_be_processed()):
                node = self.process(node)
                node.parent = parent_node
                parent_node.children[child_id] = node
        return node

    def traverse(self, node):
        return self.traverse_node_recursively(node)


class NodeDepthLeftTypeFetcher(object):
    """Helper class to fetch a node of a specific type. Depth left-first traversal"""

    def __init__(self, node_type, root_node):
        assert issubclass(node_type, Node)
        self.node_type = node_type
        self.root_node = root_node
        self.queue = deque()
        self.add_nodes_left_to_right_depth_first(self.root_node)

    def add_nodes_left_to_right_depth_first(self, node):
        if isinstance(node, Node) and node.children is not None:
            for child_node in node.children:
                self.add_nodes_left_to_right_depth_first(child_node)
                self.queue.append(child_node)
        self.queue.append(node)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        while len(self.queue) > 0:
            candidate = self.queue.popleft()
            if isinstance(candidate, self.node_type):
                return candidate
        raise StopIteration
