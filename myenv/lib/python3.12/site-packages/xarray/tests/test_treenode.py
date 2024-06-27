from __future__ import annotations

from collections.abc import Iterator
from typing import cast

import pytest

from xarray.core.iterators import LevelOrderIter
from xarray.core.treenode import InvalidTreeError, NamedNode, NodePath, TreeNode


class TestFamilyTree:
    def test_lonely(self):
        root: TreeNode = TreeNode()
        assert root.parent is None
        assert root.children == {}

    def test_parenting(self):
        john: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        mary._set_parent(john, "Mary")

        assert mary.parent == john
        assert john.children["Mary"] is mary

    def test_no_time_traveller_loops(self):
        john: TreeNode = TreeNode()

        with pytest.raises(InvalidTreeError, match="cannot be a parent of itself"):
            john._set_parent(john, "John")

        with pytest.raises(InvalidTreeError, match="cannot be a parent of itself"):
            john.children = {"John": john}

        mary: TreeNode = TreeNode()
        rose: TreeNode = TreeNode()
        mary._set_parent(john, "Mary")
        rose._set_parent(mary, "Rose")

        with pytest.raises(InvalidTreeError, match="is already a descendant"):
            john._set_parent(rose, "John")

        with pytest.raises(InvalidTreeError, match="is already a descendant"):
            rose.children = {"John": john}

    def test_parent_swap(self):
        john: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        mary._set_parent(john, "Mary")

        steve: TreeNode = TreeNode()
        mary._set_parent(steve, "Mary")

        assert mary.parent == steve
        assert steve.children["Mary"] is mary
        assert "Mary" not in john.children

    def test_multi_child_family(self):
        mary: TreeNode = TreeNode()
        kate: TreeNode = TreeNode()
        john: TreeNode = TreeNode(children={"Mary": mary, "Kate": kate})
        assert john.children["Mary"] is mary
        assert john.children["Kate"] is kate
        assert mary.parent is john
        assert kate.parent is john

    def test_disown_child(self):
        mary: TreeNode = TreeNode()
        john: TreeNode = TreeNode(children={"Mary": mary})
        mary.orphan()
        assert mary.parent is None
        assert "Mary" not in john.children

    def test_doppelganger_child(self):
        kate: TreeNode = TreeNode()
        john: TreeNode = TreeNode()

        with pytest.raises(TypeError):
            john.children = {"Kate": 666}

        with pytest.raises(InvalidTreeError, match="Cannot add same node"):
            john.children = {"Kate": kate, "Evil_Kate": kate}

        john = TreeNode(children={"Kate": kate})
        evil_kate: TreeNode = TreeNode()
        evil_kate._set_parent(john, "Kate")
        assert john.children["Kate"] is evil_kate

    def test_sibling_relationships(self):
        mary: TreeNode = TreeNode()
        kate: TreeNode = TreeNode()
        ashley: TreeNode = TreeNode()
        TreeNode(children={"Mary": mary, "Kate": kate, "Ashley": ashley})
        assert kate.siblings["Mary"] is mary
        assert kate.siblings["Ashley"] is ashley
        assert "Kate" not in kate.siblings

    def test_ancestors(self):
        tony: TreeNode = TreeNode()
        michael: TreeNode = TreeNode(children={"Tony": tony})
        vito = TreeNode(children={"Michael": michael})
        assert tony.root is vito
        assert tony.parents == (michael, vito)
        assert tony.ancestors == (vito, michael, tony)


class TestGetNodes:
    def test_get_child(self):
        steven: TreeNode = TreeNode()
        sue = TreeNode(children={"Steven": steven})
        mary = TreeNode(children={"Sue": sue})
        john = TreeNode(children={"Mary": mary})

        # get child
        assert john._get_item("Mary") is mary
        assert mary._get_item("Sue") is sue

        # no child exists
        with pytest.raises(KeyError):
            john._get_item("Kate")

        # get grandchild
        assert john._get_item("Mary/Sue") is sue

        # get great-grandchild
        assert john._get_item("Mary/Sue/Steven") is steven

        # get from middle of tree
        assert mary._get_item("Sue/Steven") is steven

    def test_get_upwards(self):
        sue: TreeNode = TreeNode()
        kate: TreeNode = TreeNode()
        mary = TreeNode(children={"Sue": sue, "Kate": kate})
        john = TreeNode(children={"Mary": mary})

        assert sue._get_item("../") is mary
        assert sue._get_item("../../") is john

        # relative path
        assert sue._get_item("../Kate") is kate

    def test_get_from_root(self):
        sue: TreeNode = TreeNode()
        mary = TreeNode(children={"Sue": sue})
        john = TreeNode(children={"Mary": mary})  # noqa

        assert sue._get_item("/Mary") is mary


class TestSetNodes:
    def test_set_child_node(self):
        john: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        john._set_item("Mary", mary)

        assert john.children["Mary"] is mary
        assert isinstance(mary, TreeNode)
        assert mary.children == {}
        assert mary.parent is john

    def test_child_already_exists(self):
        mary: TreeNode = TreeNode()
        john: TreeNode = TreeNode(children={"Mary": mary})
        mary_2: TreeNode = TreeNode()
        with pytest.raises(KeyError):
            john._set_item("Mary", mary_2, allow_overwrite=False)

    def test_set_grandchild(self):
        rose: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        john: TreeNode = TreeNode()

        john._set_item("Mary", mary)
        john._set_item("Mary/Rose", rose)

        assert john.children["Mary"] is mary
        assert isinstance(mary, TreeNode)
        assert "Rose" in mary.children
        assert rose.parent is mary

    def test_create_intermediate_child(self):
        john: TreeNode = TreeNode()
        rose: TreeNode = TreeNode()

        # test intermediate children not allowed
        with pytest.raises(KeyError, match="Could not reach"):
            john._set_item(path="Mary/Rose", item=rose, new_nodes_along_path=False)

        # test intermediate children allowed
        john._set_item("Mary/Rose", rose, new_nodes_along_path=True)
        assert "Mary" in john.children
        mary = john.children["Mary"]
        assert isinstance(mary, TreeNode)
        assert mary.children == {"Rose": rose}
        assert rose.parent == mary
        assert rose.parent == mary

    def test_overwrite_child(self):
        john: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        john._set_item("Mary", mary)

        # test overwriting not allowed
        marys_evil_twin: TreeNode = TreeNode()
        with pytest.raises(KeyError, match="Already a node object"):
            john._set_item("Mary", marys_evil_twin, allow_overwrite=False)
        assert john.children["Mary"] is mary
        assert marys_evil_twin.parent is None

        # test overwriting allowed
        marys_evil_twin = TreeNode()
        john._set_item("Mary", marys_evil_twin, allow_overwrite=True)
        assert john.children["Mary"] is marys_evil_twin
        assert marys_evil_twin.parent is john


class TestPruning:
    def test_del_child(self):
        john: TreeNode = TreeNode()
        mary: TreeNode = TreeNode()
        john._set_item("Mary", mary)

        del john["Mary"]
        assert "Mary" not in john.children
        assert mary.parent is None

        with pytest.raises(KeyError):
            del john["Mary"]


def create_test_tree() -> tuple[NamedNode, NamedNode]:
    # a
    # ├── b
    # │   ├── d
    # │   └── e
    # │       ├── f
    # │       └── g
    # └── c
    #     └── h
    #         └── i
    a: NamedNode = NamedNode(name="a")
    b: NamedNode = NamedNode()
    c: NamedNode = NamedNode()
    d: NamedNode = NamedNode()
    e: NamedNode = NamedNode()
    f: NamedNode = NamedNode()
    g: NamedNode = NamedNode()
    h: NamedNode = NamedNode()
    i: NamedNode = NamedNode()

    a.children = {"b": b, "c": c}
    b.children = {"d": d, "e": e}
    e.children = {"f": f, "g": g}
    c.children = {"h": h}
    h.children = {"i": i}

    return a, f


class TestIterators:

    def test_levelorderiter(self):
        root, _ = create_test_tree()
        result: list[str | None] = [
            node.name for node in cast(Iterator[NamedNode], LevelOrderIter(root))
        ]
        expected = [
            "a",  # root Node is unnamed
            "b",
            "c",
            "d",
            "e",
            "h",
            "f",
            "g",
            "i",
        ]
        assert result == expected


class TestAncestry:

    def test_parents(self):
        _, leaf_f = create_test_tree()
        expected = ["e", "b", "a"]
        assert [node.name for node in leaf_f.parents] == expected

    def test_lineage(self):
        _, leaf_f = create_test_tree()
        expected = ["f", "e", "b", "a"]
        assert [node.name for node in leaf_f.lineage] == expected

    def test_ancestors(self):
        _, leaf_f = create_test_tree()
        ancestors = leaf_f.ancestors
        expected = ["a", "b", "e", "f"]
        for node, expected_name in zip(ancestors, expected):
            assert node.name == expected_name

    def test_subtree(self):
        root, _ = create_test_tree()
        subtree = root.subtree
        expected = [
            "a",
            "b",
            "c",
            "d",
            "e",
            "h",
            "f",
            "g",
            "i",
        ]
        for node, expected_name in zip(subtree, expected):
            assert node.name == expected_name

    def test_descendants(self):
        root, _ = create_test_tree()
        descendants = root.descendants
        expected = [
            "b",
            "c",
            "d",
            "e",
            "h",
            "f",
            "g",
            "i",
        ]
        for node, expected_name in zip(descendants, expected):
            assert node.name == expected_name

    def test_leaves(self):
        tree, _ = create_test_tree()
        leaves = tree.leaves
        expected = [
            "d",
            "f",
            "g",
            "i",
        ]
        for node, expected_name in zip(leaves, expected):
            assert node.name == expected_name

    def test_levels(self):
        a, f = create_test_tree()

        assert a.level == 0
        assert f.level == 3

        assert a.depth == 3
        assert f.depth == 3

        assert a.width == 1
        assert f.width == 3


class TestRenderTree:
    def test_render_nodetree(self):
        sam: NamedNode = NamedNode()
        ben: NamedNode = NamedNode()
        mary: NamedNode = NamedNode(children={"Sam": sam, "Ben": ben})
        kate: NamedNode = NamedNode()
        john: NamedNode = NamedNode(children={"Mary": mary, "Kate": kate})
        expected_nodes = [
            "NamedNode()",
            "\tNamedNode('Mary')",
            "\t\tNamedNode('Sam')",
            "\t\tNamedNode('Ben')",
            "\tNamedNode('Kate')",
        ]
        expected_str = "NamedNode('Mary')"
        john_repr = john.__repr__()
        mary_str = mary.__str__()

        assert mary_str == expected_str

        john_nodes = john_repr.splitlines()
        assert len(john_nodes) == len(expected_nodes)
        for expected_node, repr_node in zip(expected_nodes, john_nodes):
            assert expected_node == repr_node


def test_nodepath():
    path = NodePath("/Mary")
    assert path.root == "/"
    assert path.stem == "Mary"
