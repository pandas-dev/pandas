from __future__ import annotations

from typing import Any

import pytest

pytest.importorskip("sqlalchemy")


def test_modulo_operator(sql_con: Any) -> None:
    query = "SELECT 10 % 3"
    result = sql_con.execute(query)
    assert result.scalar() == 1


def test_like_pattern(sql_con: Any) -> None:
    query = "SELECT 'abc' LIKE 'a%'"
    result = sql_con.execute(query)
    assert result.scalar() == 1


def test_sqlalchemy_selectable(sql_con: Any) -> None:
    from sqlalchemy import (
        literal,
        select,
    )

    stmt = select(literal("hello"))
    result = sql_con.execute(stmt)
    assert result.scalar() == "hello"
