from unittest import mock

import pytest
from traitlets import Bool, Bunch

from jupyterlite_core import addons
from jupyterlite_core.addons.base import BaseAddon


@pytest.mark.parametrize(
    "expect_warn,aliases,flags",
    [
        ["--foo cannot be redefined", {"foo": "Foo.foo"}, {}],
        [
            "--foo cannot be redefined",
            {"foo": "Foo.foo"},
            {"bar": ({"LiteBuildConfig": {"log_level": "DEBUG"}}, "woo")},
        ],
        [
            "--bar cannot redefine",
            {},
            {"bar": ({"BadAddon": {"foo": False}}, "boo")},
        ],
    ],
)
def test_cli_no_redefine_warnings(aliases, flags, expect_warn, some_entry_point_addons):
    assert len(addons.get_addon_entry_points(1)) == 1
    with pytest.warns() as warned:
        aliases = addons.merge_addon_aliases(aliases, force=True)
        flags = addons.merge_addon_flags(flags, force=True)
    print("aliases", aliases)
    print("flags", flags)
    assert len(warned) == 1
    assert expect_warn in f"{warned[0].message}"


@pytest.fixture
def some_entry_point_addons():
    class BadAddon(BaseAddon):
        # traits
        foo = Bool(False, help="just foo").tag(config=True)
        baz = Bool(False, help="just baz").tag(config=True)

        # CLI
        aliases = {"foo": "BadAddon.foo"}
        flags = {"bar": ({"BadAddon": {"foo": True}}, "foo")}

    with mock.patch.multiple(
        addons, entry_points=lambda group: [Bunch(name="foo", load=lambda: BadAddon)]
    ):
        yield
