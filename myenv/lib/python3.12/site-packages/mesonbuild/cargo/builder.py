# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2023 Intel Corporation

"""Provides helpers for building AST

This is meant to make building Meson AST from foreign (largely declarative)
build descriptions easier.
"""

from __future__ import annotations
import builtins
import dataclasses
import typing as T

from .. import mparser


def _token(tid: str, filename: str, value: mparser.TV_TokenTypes) -> mparser.Token[mparser.TV_TokenTypes]:
    """Create a Token object, but with the line numbers stubbed out.

    :param tid: the token id (such as string, number, etc)
    :param filename: the filename that the token was generated from
    :param value: the value of the token
    :return: A Token object
    """
    return mparser.Token(tid, filename, -1, -1, -1, (-1, -1), value)


def string(value: str, filename: str) -> mparser.StringNode:
    """Build A StringNode

    :param value: the value of the string
    :param filename: the file that the value came from
    :return: A StringNode
    """
    return mparser.StringNode(_token('string', filename, value))


def number(value: int, filename: str) -> mparser.NumberNode:
    """Build A NumberNode

    :param value: the value of the number
    :param filename: the file that the value came from
    :return: A NumberNode
    """
    return mparser.NumberNode(_token('number', filename, value))


def bool(value: builtins.bool, filename: str) -> mparser.BooleanNode:
    """Build A BooleanNode

    :param value: the value of the boolean
    :param filename: the file that the value came from
    :return: A BooleanNode
    """
    return mparser.BooleanNode(_token('bool', filename, value))


def array(value: T.List[mparser.BaseNode], filename: str) -> mparser.ArrayNode:
    """Build an Array Node

    :param value: A list of nodes to insert into the array
    :param filename: The file the array is from
    :return: An ArrayNode built from the arguments
    """
    args = mparser.ArgumentNode(_token('array', filename, 'unused'))
    args.arguments = value
    return mparser.ArrayNode(args, -1, -1, -1, -1)


def identifier(value: str, filename: str) -> mparser.IdNode:
    """Build A IdNode

    :param value: the value of the boolean
    :param filename: the file that the value came from
    :return: A BooleanNode
    """
    return mparser.IdNode(_token('id', filename, value))


def method(name: str, id_: mparser.IdNode,
           pos: T.Optional[T.List[mparser.BaseNode]] = None,
           kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
           ) -> mparser.MethodNode:
    """Create a method call.

    :param name: the name of the method
    :param id_: the object to call the method of
    :param pos: a list of positional arguments, defaults to None
    :param kw: a dictionary of keyword arguments, defaults to None
    :return: a method call object
    """
    args = mparser.ArgumentNode(_token('array', id_.filename, 'unused'))
    if pos is not None:
        args.arguments = pos
    if kw is not None:
        args.kwargs = {identifier(k, id_.filename): v for k, v in kw.items()}
    return mparser.MethodNode(id_.filename, -1, -1, id_, name, args)


def function(name: str, filename: str,
             pos: T.Optional[T.List[mparser.BaseNode]] = None,
             kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
             ) -> mparser.FunctionNode:
    """Create a function call.

    :param name: the name of the function
    :param filename: The name of the current file being evaluated
    :param pos: a list of positional arguments, defaults to None
    :param kw: a dictionary of keyword arguments, defaults to None
    :return: a method call object
    """
    args = mparser.ArgumentNode(_token('array', filename, 'unused'))
    if pos is not None:
        args.arguments = pos
    if kw is not None:
        args.kwargs = {identifier(k, filename): v for k, v in kw.items()}
    return mparser.FunctionNode(filename, -1, -1, -1, -1, name, args)


def equal(lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
    """Create an equality operation

    :param lhs: The left hand side of the equal
    :param rhs: the right hand side of the equal
    :return: A compraison node
    """
    return mparser.ComparisonNode('==', lhs, rhs)


def or_(lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.OrNode:
    """Create and OrNode

    :param lhs: The Left of the Node
    :param rhs: The Right of the Node
    :return: The OrNode
    """
    return mparser.OrNode(lhs, rhs)


def and_(lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.AndNode:
    """Create an AndNode

    :param lhs: The left of the And
    :param rhs: The right of the And
    :return: The AndNode
    """
    return mparser.AndNode(lhs, rhs)


def not_(value: mparser.BaseNode, filename: str) -> mparser.NotNode:
    """Create a not node

    :param value: The value to negate
    :param filename: the string filename
    :return: The NotNode
    """
    return mparser.NotNode(_token('not', filename, ''), value)


def assign(value: mparser.BaseNode, varname: str, filename: str) -> mparser.AssignmentNode:
    """Create an AssignmentNode

    :param value: The rvalue
    :param varname: The lvalue
    :param filename: The filename
    :return: An AssignmentNode
    """
    return mparser.AssignmentNode(filename, -1, -1, varname, value)


def block(filename: str) -> mparser.CodeBlockNode:
    return mparser.CodeBlockNode(_token('node', filename, ''))


@dataclasses.dataclass
class Builder:

    filename: str

    def assign(self, value: mparser.BaseNode, varname: str) -> mparser.AssignmentNode:
        return assign(value, varname, self.filename)

    def string(self, value: str) -> mparser.StringNode:
        """Build A StringNode

        :param value: the value of the string
        :return: A StringNode
        """
        return string(value, self.filename)

    def number(self, value: int) -> mparser.NumberNode:
        """Build A NumberNode

        :param value: the value of the number
        :return: A NumberNode
        """
        return number(value, self.filename)

    def bool(self, value: builtins.bool) -> mparser.BooleanNode:
        """Build A BooleanNode

        :param value: the value of the boolean
        :return: A BooleanNode
        """
        return bool(value, self.filename)

    def array(self, value: T.List[mparser.BaseNode]) -> mparser.ArrayNode:
        """Build an Array Node

        :param value: A list of nodes to insert into the array
        :return: An ArrayNode built from the arguments
        """
        return array(value, self.filename)

    def identifier(self, value: str) -> mparser.IdNode:
        """Build A IdNode

        :param value: the value of the boolean
        :return: A BooleanNode
        """
        return identifier(value, self.filename)

    def method(self, name: str, id_: mparser.IdNode,
               pos: T.Optional[T.List[mparser.BaseNode]] = None,
               kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
               ) -> mparser.MethodNode:
        """Create a method call.

        :param name: the name of the method
        :param id_: the object to call the method of
        :param pos: a list of positional arguments, defaults to None
        :param kw: a dictionary of keyword arguments, defaults to None
        :return: a method call object
        """
        return method(name, id_, pos or [], kw or {})

    def function(self, name: str,
                 pos: T.Optional[T.List[mparser.BaseNode]] = None,
                 kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
                 ) -> mparser.FunctionNode:
        """Create a function call.

        :param name: the name of the function
        :param pos: a list of positional arguments, defaults to None
        :param kw: a dictionary of keyword arguments, defaults to None
        :return: a method call object
        """
        return function(name, self.filename, pos or [], kw or {})

    def equal(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
        """Create an equality operation

        :param lhs: The left hand side of the equal
        :param rhs: the right hand side of the equal
        :return: A compraison node
        """
        return equal(lhs, rhs)

    def or_(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.OrNode:
        """Create and OrNode

        :param lhs: The Left of the Node
        :param rhs: The Right of the Node
        :return: The OrNode
        """
        return or_(lhs, rhs)

    def and_(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.AndNode:
        """Create an AndNode

        :param lhs: The left of the And
        :param rhs: The right of the And
        :return: The AndNode
        """
        return and_(lhs, rhs)

    def not_(self, value: mparser.BaseNode, filename: str) -> mparser.NotNode:
        """Create a not node

        :param value: The value to negate
        :return: The NotNode
        """
        return not_(value, self.filename)
