# dialects/mysql/json.py
# Copyright (C) 2005-2025 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
from __future__ import annotations

from typing import Any
from typing import TYPE_CHECKING

from ... import types as sqltypes

if TYPE_CHECKING:
    from ...engine.interfaces import Dialect
    from ...sql.type_api import _BindProcessorType
    from ...sql.type_api import _LiteralProcessorType


class JSON(sqltypes.JSON):
    """MySQL JSON type.

    MySQL supports JSON as of version 5.7.
    MariaDB supports JSON (as an alias for LONGTEXT) as of version 10.2.

    :class:`_mysql.JSON` is used automatically whenever the base
    :class:`_types.JSON` datatype is used against a MySQL or MariaDB backend.

    .. seealso::

        :class:`_types.JSON` - main documentation for the generic
        cross-platform JSON datatype.

    The :class:`.mysql.JSON` type supports persistence of JSON values
    as well as the core index operations provided by :class:`_types.JSON`
    datatype, by adapting the operations to render the ``JSON_EXTRACT``
    function at the database level.

    """

    pass


class _FormatTypeMixin:
    def _format_value(self, value: Any) -> str:
        raise NotImplementedError()

    def bind_processor(self, dialect: Dialect) -> _BindProcessorType[Any]:
        super_proc = self.string_bind_processor(dialect)  # type: ignore[attr-defined]  # noqa: E501

        def process(value: Any) -> Any:
            value = self._format_value(value)
            if super_proc:
                value = super_proc(value)
            return value

        return process

    def literal_processor(
        self, dialect: Dialect
    ) -> _LiteralProcessorType[Any]:
        super_proc = self.string_literal_processor(dialect)  # type: ignore[attr-defined]  # noqa: E501

        def process(value: Any) -> str:
            value = self._format_value(value)
            if super_proc:
                value = super_proc(value)
            return value  # type: ignore[no-any-return]

        return process


class JSONIndexType(_FormatTypeMixin, sqltypes.JSON.JSONIndexType):
    def _format_value(self, value: Any) -> str:
        if isinstance(value, int):
            formatted_value = "$[%s]" % value
        else:
            formatted_value = '$."%s"' % value
        return formatted_value


class JSONPathType(_FormatTypeMixin, sqltypes.JSON.JSONPathType):
    def _format_value(self, value: Any) -> str:
        return "$%s" % (
            "".join(
                [
                    "[%s]" % elem if isinstance(elem, int) else '."%s"' % elem
                    for elem in value
                ]
            )
        )
