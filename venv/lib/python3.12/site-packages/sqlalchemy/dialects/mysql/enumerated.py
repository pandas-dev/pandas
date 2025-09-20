# dialects/mysql/enumerated.py
# Copyright (C) 2005-2025 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php

from __future__ import annotations

import enum
import re
from typing import Any
from typing import Dict
from typing import Optional
from typing import Set
from typing import Type
from typing import TYPE_CHECKING
from typing import Union

from .types import _StringType
from ... import exc
from ... import sql
from ... import util
from ...sql import sqltypes
from ...sql import type_api

if TYPE_CHECKING:
    from ...engine.interfaces import Dialect
    from ...sql.elements import ColumnElement
    from ...sql.type_api import _BindProcessorType
    from ...sql.type_api import _ResultProcessorType
    from ...sql.type_api import TypeEngine
    from ...sql.type_api import TypeEngineMixin


class ENUM(type_api.NativeForEmulated, sqltypes.Enum, _StringType):
    """MySQL ENUM type."""

    __visit_name__ = "ENUM"

    native_enum = True

    def __init__(self, *enums: Union[str, Type[enum.Enum]], **kw: Any) -> None:
        """Construct an ENUM.

        E.g.::

          Column("myenum", ENUM("foo", "bar", "baz"))

        :param enums: The range of valid values for this ENUM.  Values in
          enums are not quoted, they will be escaped and surrounded by single
          quotes when generating the schema.  This object may also be a
          PEP-435-compliant enumerated type.

          .. versionadded: 1.1 added support for PEP-435-compliant enumerated
             types.

        :param strict: This flag has no effect.

         .. versionchanged:: The MySQL ENUM type as well as the base Enum
            type now validates all Python data values.

        :param charset: Optional, a column-level character set for this string
          value.  Takes precedence to 'ascii' or 'unicode' short-hand.

        :param collation: Optional, a column-level collation for this string
          value.  Takes precedence to 'binary' short-hand.

        :param ascii: Defaults to False: short-hand for the ``latin1``
          character set, generates ASCII in schema.

        :param unicode: Defaults to False: short-hand for the ``ucs2``
          character set, generates UNICODE in schema.

        :param binary: Defaults to False: short-hand, pick the binary
          collation type that matches the column's character set.  Generates
          BINARY in schema.  This does not affect the type of data stored,
          only the collation of character data.

        """
        kw.pop("strict", None)
        self._enum_init(enums, kw)  # type: ignore[arg-type]
        _StringType.__init__(self, length=self.length, **kw)

    @classmethod
    def adapt_emulated_to_native(
        cls,
        impl: Union[TypeEngine[Any], TypeEngineMixin],
        **kw: Any,
    ) -> ENUM:
        """Produce a MySQL native :class:`.mysql.ENUM` from plain
        :class:`.Enum`.

        """
        if TYPE_CHECKING:
            assert isinstance(impl, ENUM)
        kw.setdefault("validate_strings", impl.validate_strings)
        kw.setdefault("values_callable", impl.values_callable)
        kw.setdefault("omit_aliases", impl._omit_aliases)
        return cls(**kw)

    def _object_value_for_elem(self, elem: str) -> Union[str, enum.Enum]:
        # mysql sends back a blank string for any value that
        # was persisted that was not in the enums; that is, it does no
        # validation on the incoming data, it "truncates" it to be
        # the blank string.  Return it straight.
        if elem == "":
            return elem
        else:
            return super()._object_value_for_elem(elem)

    def __repr__(self) -> str:
        return util.generic_repr(
            self, to_inspect=[ENUM, _StringType, sqltypes.Enum]
        )


# TODO: SET is a string as far as configuration but does not act like
# a string at the python level.  We either need to make a py-type agnostic
# version of String as a base to be used for this, make this some kind of
# TypeDecorator, or just vendor it out as its own type.
class SET(_StringType):
    """MySQL SET type."""

    __visit_name__ = "SET"

    def __init__(self, *values: str, **kw: Any):
        """Construct a SET.

        E.g.::

          Column("myset", SET("foo", "bar", "baz"))

        The list of potential values is required in the case that this
        set will be used to generate DDL for a table, or if the
        :paramref:`.SET.retrieve_as_bitwise` flag is set to True.

        :param values: The range of valid values for this SET. The values
          are not quoted, they will be escaped and surrounded by single
          quotes when generating the schema.

        :param convert_unicode: Same flag as that of
         :paramref:`.String.convert_unicode`.

        :param collation: same as that of :paramref:`.String.collation`

        :param charset: same as that of :paramref:`.VARCHAR.charset`.

        :param ascii: same as that of :paramref:`.VARCHAR.ascii`.

        :param unicode: same as that of :paramref:`.VARCHAR.unicode`.

        :param binary: same as that of :paramref:`.VARCHAR.binary`.

        :param retrieve_as_bitwise: if True, the data for the set type will be
          persisted and selected using an integer value, where a set is coerced
          into a bitwise mask for persistence.  MySQL allows this mode which
          has the advantage of being able to store values unambiguously,
          such as the blank string ``''``.   The datatype will appear
          as the expression ``col + 0`` in a SELECT statement, so that the
          value is coerced into an integer value in result sets.
          This flag is required if one wishes
          to persist a set that can store the blank string ``''`` as a value.

          .. warning::

            When using :paramref:`.mysql.SET.retrieve_as_bitwise`, it is
            essential that the list of set values is expressed in the
            **exact same order** as exists on the MySQL database.

        """
        self.retrieve_as_bitwise = kw.pop("retrieve_as_bitwise", False)
        self.values = tuple(values)
        if not self.retrieve_as_bitwise and "" in values:
            raise exc.ArgumentError(
                "Can't use the blank value '' in a SET without "
                "setting retrieve_as_bitwise=True"
            )
        if self.retrieve_as_bitwise:
            self._inversed_bitmap: Dict[str, int] = {
                value: 2**idx for idx, value in enumerate(self.values)
            }
            self._bitmap: Dict[int, str] = {
                2**idx: value for idx, value in enumerate(self.values)
            }
        length = max([len(v) for v in values] + [0])
        kw.setdefault("length", length)
        super().__init__(**kw)

    def column_expression(
        self, colexpr: ColumnElement[Any]
    ) -> ColumnElement[Any]:
        if self.retrieve_as_bitwise:
            return sql.type_coerce(
                sql.type_coerce(colexpr, sqltypes.Integer) + 0, self
            )
        else:
            return colexpr

    def result_processor(
        self, dialect: Dialect, coltype: Any
    ) -> Optional[_ResultProcessorType[Any]]:
        if self.retrieve_as_bitwise:

            def process(value: Union[str, int, None]) -> Optional[Set[str]]:
                if value is not None:
                    value = int(value)

                    return set(util.map_bits(self._bitmap.__getitem__, value))
                else:
                    return None

        else:
            super_convert = super().result_processor(dialect, coltype)

            def process(value: Union[str, Set[str], None]) -> Optional[Set[str]]:  # type: ignore[misc]  # noqa: E501
                if isinstance(value, str):
                    # MySQLdb returns a string, let's parse
                    if super_convert:
                        value = super_convert(value)
                        assert value is not None
                    if TYPE_CHECKING:
                        assert isinstance(value, str)
                    return set(re.findall(r"[^,]+", value))
                else:
                    # mysql-connector-python does a naive
                    # split(",") which throws in an empty string
                    if value is not None:
                        value.discard("")
                    return value

        return process

    def bind_processor(
        self, dialect: Dialect
    ) -> _BindProcessorType[Union[str, int]]:
        super_convert = super().bind_processor(dialect)
        if self.retrieve_as_bitwise:

            def process(
                value: Union[str, int, set[str], None],
            ) -> Union[str, int, None]:
                if value is None:
                    return None
                elif isinstance(value, (int, str)):
                    if super_convert:
                        return super_convert(value)  # type: ignore[arg-type, no-any-return]  # noqa: E501
                    else:
                        return value
                else:
                    int_value = 0
                    for v in value:
                        int_value |= self._inversed_bitmap[v]
                    return int_value

        else:

            def process(
                value: Union[str, int, set[str], None],
            ) -> Union[str, int, None]:
                # accept strings and int (actually bitflag) values directly
                if value is not None and not isinstance(value, (int, str)):
                    value = ",".join(value)
                if super_convert:
                    return super_convert(value)  # type: ignore
                else:
                    return value

        return process

    def adapt(self, cls: type, **kw: Any) -> Any:
        kw["retrieve_as_bitwise"] = self.retrieve_as_bitwise
        return util.constructor_copy(self, cls, *self.values, **kw)

    def __repr__(self) -> str:
        return util.generic_repr(
            self,
            to_inspect=[SET, _StringType],
            additional_kw=[
                ("retrieve_as_bitwise", False),
            ],
        )
