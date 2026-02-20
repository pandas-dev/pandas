# testing/suite/test_cte.py
# Copyright (C) 2005-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors

from .. import fixtures
from ..assertions import eq_
from ..schema import Column
from ..schema import Table
from ... import column
from ... import ForeignKey
from ... import Integer
from ... import select
from ... import String
from ... import testing
from ... import values


class CTETest(fixtures.TablesTest):
    __sparse_driver_backend__ = True
    __requires__ = ("ctes",)

    run_inserts = "each"
    run_deletes = "each"

    @classmethod
    def define_tables(cls, metadata):
        Table(
            "some_table",
            metadata,
            Column("id", Integer, primary_key=True),
            Column("data", String(50)),
            Column("parent_id", ForeignKey("some_table.id")),
        )

        Table(
            "some_other_table",
            metadata,
            Column("id", Integer, primary_key=True),
            Column("data", String(50)),
            Column("parent_id", Integer),
        )

    @classmethod
    def insert_data(cls, connection):
        connection.execute(
            cls.tables.some_table.insert(),
            [
                {"id": 1, "data": "d1", "parent_id": None},
                {"id": 2, "data": "d2", "parent_id": 1},
                {"id": 3, "data": "d3", "parent_id": 1},
                {"id": 4, "data": "d4", "parent_id": 3},
                {"id": 5, "data": "d5", "parent_id": 3},
            ],
        )

    def test_select_nonrecursive_round_trip(self, connection):
        some_table = self.tables.some_table

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte")
        )
        result = connection.execute(
            select(cte.c.data).where(cte.c.data.in_(["d4", "d5"]))
        )
        eq_(result.fetchall(), [("d4",)])

    def test_select_recursive_round_trip(self, connection):
        some_table = self.tables.some_table

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte", recursive=True)
        )

        cte_alias = cte.alias("c1")
        st1 = some_table.alias()
        # note that SQL Server requires this to be UNION ALL,
        # can't be UNION
        cte = cte.union_all(
            select(st1).where(st1.c.id == cte_alias.c.parent_id)
        )
        result = connection.execute(
            select(cte.c.data)
            .where(cte.c.data != "d2")
            .order_by(cte.c.data.desc())
        )
        eq_(
            result.fetchall(),
            [("d4",), ("d3",), ("d3",), ("d1",), ("d1",), ("d1",)],
        )

    def test_insert_from_select_round_trip(self, connection):
        some_table = self.tables.some_table
        some_other_table = self.tables.some_other_table

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte")
        )
        connection.execute(
            some_other_table.insert().from_select(
                ["id", "data", "parent_id"], select(cte)
            )
        )
        eq_(
            connection.execute(
                select(some_other_table).order_by(some_other_table.c.id)
            ).fetchall(),
            [(2, "d2", 1), (3, "d3", 1), (4, "d4", 3)],
        )

    @testing.requires.ctes_with_update_delete
    @testing.requires.update_from
    def test_update_from_round_trip(self, connection):
        some_table = self.tables.some_table
        some_other_table = self.tables.some_other_table

        connection.execute(
            some_other_table.insert().from_select(
                ["id", "data", "parent_id"], select(some_table)
            )
        )

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte")
        )
        connection.execute(
            some_other_table.update()
            .values(parent_id=5)
            .where(some_other_table.c.data == cte.c.data)
        )
        eq_(
            connection.execute(
                select(some_other_table).order_by(some_other_table.c.id)
            ).fetchall(),
            [
                (1, "d1", None),
                (2, "d2", 5),
                (3, "d3", 5),
                (4, "d4", 5),
                (5, "d5", 3),
            ],
        )

    @testing.requires.ctes_with_update_delete
    @testing.requires.delete_from
    def test_delete_from_round_trip(self, connection):
        some_table = self.tables.some_table
        some_other_table = self.tables.some_other_table

        connection.execute(
            some_other_table.insert().from_select(
                ["id", "data", "parent_id"], select(some_table)
            )
        )

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte")
        )
        connection.execute(
            some_other_table.delete().where(
                some_other_table.c.data == cte.c.data
            )
        )
        eq_(
            connection.execute(
                select(some_other_table).order_by(some_other_table.c.id)
            ).fetchall(),
            [(1, "d1", None), (5, "d5", 3)],
        )

    @testing.requires.ctes_with_update_delete
    def test_delete_scalar_subq_round_trip(self, connection):
        some_table = self.tables.some_table
        some_other_table = self.tables.some_other_table

        connection.execute(
            some_other_table.insert().from_select(
                ["id", "data", "parent_id"], select(some_table)
            )
        )

        cte = (
            select(some_table)
            .where(some_table.c.data.in_(["d2", "d3", "d4"]))
            .cte("some_cte")
        )
        connection.execute(
            some_other_table.delete().where(
                some_other_table.c.data
                == select(cte.c.data)
                .where(cte.c.id == some_other_table.c.id)
                .scalar_subquery()
            )
        )
        eq_(
            connection.execute(
                select(some_other_table).order_by(some_other_table.c.id)
            ).fetchall(),
            [(1, "d1", None), (5, "d5", 3)],
        )

    @testing.variation("values_named", [True, False])
    @testing.variation("cte_named", [True, False])
    @testing.variation("literal_binds", [True, False])
    @testing.requires.ctes_with_values
    def test_values_named_via_cte(
        self, connection, values_named, cte_named, literal_binds
    ):

        cte1 = (
            values(
                column("col1", String),
                column("col2", Integer),
                literal_binds=bool(literal_binds),
                name="some name" if values_named else None,
            )
            .data([("a", 2), ("b", 3)])
            .cte("cte1" if cte_named else None)
        )

        stmt = select(cte1)

        rows = connection.execute(stmt).all()
        eq_(rows, [("a", 2), ("b", 3)])
