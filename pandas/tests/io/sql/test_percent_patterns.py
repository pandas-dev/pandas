# pandas/tests/io/sql/test_percent_patterns.py
import os
import pytest
from typing import TYPE_CHECKING

sa = pytest.importorskip("sqlalchemy")

if TYPE_CHECKING:
    from sqlalchemy.engine import Engine

PG = os.environ.get("PANDAS_TEST_POSTGRES_URI")
URL = PG or "sqlite+pysqlite:///:memory:"


def _eng() -> "Engine":
    return sa.create_engine(URL)


def test_text_modulo() -> None:
    import pandas as pd

    with _eng().connect() as c:
        df = pd.read_sql(sa.text("SELECT 5 % 2 AS r"), c)
    assert df.iloc[0, 0] == 1


def test_like_single_percent() -> None:
    import pandas as pd

    with _eng().connect() as c:
        df = pd.read_sql(
            sa.text("SELECT 'John' AS fullname WHERE 'John' LIKE 'John%'"),
            c,
        )
    assert len(df) == 1


def test_sqlalchemy_expr_percent_operator() -> None:
    from sqlalchemy import (
        literal,
        select,
    )

    import pandas as pd

    with _eng().connect() as c:
        df = pd.read_sql(select((literal(7) % literal(3)).label("r")), c)
    assert df.iloc[0, 0] == 1
