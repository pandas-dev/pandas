from __future__ import annotations

import pytest
from typing import Any

pytest.importorskip("sqlalchemy")


def test_modulo_operator(sql_con: Any) -> None:
    # Example test for modulo operator escaping
    query = "SELECT 10 % 3"
    result = sql_con.execute(query)
    assert result.scalar() == 1


def test_like_pattern(sql_con: Any) -> None:
    # Example test for LIKE pattern with percent signs
    query = "SELECT 'abc' LIKE 'a%'"
    result = sql_con.execute(query)
    assert result.scalar() == 1


def test_sqlalchemy_selectable(sql_con: Any) -> None:
    # Example test using a SQLAlchemy selectable
    from sqlalchemy import select, literal

    stmt = select(literal("hello"))
    result = sql_con.execute(stmt)
    assert result.scalar() == "hello"
