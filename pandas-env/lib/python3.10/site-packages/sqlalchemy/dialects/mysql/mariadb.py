# dialects/mysql/mariadb.py
# Copyright (C) 2005-2025 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors
from .base import MariaDBIdentifierPreparer
from .base import MySQLDialect
from .base import MySQLTypeCompiler
from ...sql import sqltypes


class INET4(sqltypes.TypeEngine[str]):
    """INET4 column type for MariaDB

    .. versionadded:: 2.0.37
    """

    __visit_name__ = "INET4"


class INET6(sqltypes.TypeEngine[str]):
    """INET6 column type for MariaDB

    .. versionadded:: 2.0.37
    """

    __visit_name__ = "INET6"


class MariaDBTypeCompiler(MySQLTypeCompiler):
    def visit_INET4(self, type_, **kwargs) -> str:
        return "INET4"

    def visit_INET6(self, type_, **kwargs) -> str:
        return "INET6"


class MariaDBDialect(MySQLDialect):
    is_mariadb = True
    supports_statement_cache = True
    name = "mariadb"
    preparer = MariaDBIdentifierPreparer
    type_compiler_cls = MariaDBTypeCompiler


def loader(driver):
    driver_mod = __import__(
        "sqlalchemy.dialects.mysql.%s" % driver
    ).dialects.mysql
    driver_cls = getattr(driver_mod, driver).dialect

    return type(
        "MariaDBDialect_%s" % driver,
        (
            MariaDBDialect,
            driver_cls,
        ),
        {"supports_statement_cache": True},
    )
