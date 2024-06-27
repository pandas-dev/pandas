from __future__ import annotations

import os
import pathlib
import site
import stat
import sys
from collections import OrderedDict
from contextlib import contextmanager

import pytest
import yaml

import dask.config
from dask.config import (
    _get_paths,
    canonical_name,
    collect,
    collect_env,
    collect_yaml,
    config,
    deserialize,
    ensure_file,
    expand_environment_variables,
    get,
    merge,
    paths_containing_key,
    pop,
    refresh,
    rename,
    serialize,
    update,
    update_defaults,
)


def test_canonical_name():
    c = {"foo-bar": 1, "fizz_buzz": 2}
    assert canonical_name("foo-bar", c) == "foo-bar"
    assert canonical_name("foo_bar", c) == "foo-bar"
    assert canonical_name("fizz-buzz", c) == "fizz_buzz"
    assert canonical_name("fizz_buzz", c) == "fizz_buzz"
    assert canonical_name("new-key", c) == "new-key"
    assert canonical_name("new_key", c) == "new_key"


def test_update():
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": OrderedDict({"b": 2})}
    update(b, a)
    assert b == {"x": 1, "y": {"a": 1, "b": 2}, "z": 3}

    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"a": 3, "b": 2}}
    update(b, a, priority="old")
    assert b == {"x": 2, "y": {"a": 3, "b": 2}, "z": 3}


def test_update_new_defaults():
    d = {"x": 1, "y": 1, "z": {"a": 1, "b": 1}}
    o = {"x": 1, "y": 2, "z": {"a": 1, "b": 2}, "c": 2, "c2": {"d": 2}}
    n = {"x": 3, "y": 3, "z": OrderedDict({"a": 3, "b": 3}), "c": 3, "c2": {"d": 3}}
    assert update(o, n, priority="new-defaults", defaults=d) == {
        "x": 3,
        "y": 2,
        "z": {"a": 3, "b": 2},
        "c": 2,
        "c2": {"d": 2},
    }
    assert update(o, n, priority="new-defaults", defaults=o) == update(
        o, n, priority="new"
    )
    assert update(o, n, priority="new-defaults", defaults=None) == update(
        o, n, priority="old"
    )


def test_update_defaults():
    defaults = [
        {"a": 1, "b": {"c": 1}},
        {"a": 2, "b": {"d": 2}},
    ]
    current = {"a": 2, "b": {"c": 1, "d": 3}, "extra": 0}
    new = {"a": 0, "b": {"c": 0, "d": 0}, "new-extra": 0}
    update_defaults(new, current, defaults=defaults)

    assert defaults == [
        {"a": 1, "b": {"c": 1}},
        {"a": 2, "b": {"d": 2}},
        {"a": 0, "b": {"c": 0, "d": 0}, "new-extra": 0},
    ]
    assert current == {"a": 0, "b": {"c": 0, "d": 3}, "extra": 0, "new-extra": 0}


def test_update_list_to_dict():
    a = {"x": {"y": 1, "z": 2}}
    b = {"x": [1, 2], "w": 3}
    update(b, a)
    assert b == {"x": {"y": 1, "z": 2}, "w": 3}


def test_update_dict_to_list():
    a = {"x": [1, 2]}
    b = {"x": {"y": 1, "z": 2}, "w": 3}
    update(b, a)
    assert b == {"x": [1, 2], "w": 3}


def test_merge():
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"b": 2}}

    expected = {"x": 2, "y": {"a": 1, "b": 2}, "z": 3}

    c = merge(a, b)
    assert c == expected


def test_collect_yaml_paths(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"b": 2}}

    expected = {"x": 2, "y": {"a": 1, "b": 2}, "z": 3}

    pa, pb = tmp_path / "a.yaml", tmp_path / "b.yaml"
    pa.write_text(yaml.dump(a))
    pb.write_text(yaml.dump(b))

    configs = list(collect_yaml(paths=[pa, pb], return_paths=True))
    assert configs == [(pa, a), (pb, b)]

    # `return_paths` defaults to False
    config = merge(*collect_yaml(paths=[pa, pb]))
    assert config == expected


def test_paths_containing_key(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"b": 2}}

    pa, pb = tmp_path / "a.yaml", tmp_path / "b.yaml"
    pa.write_text(yaml.dump(a))
    pb.write_text(yaml.dump(b))

    paths = list(paths_containing_key("y.a", paths=[pa, pb]))
    assert paths == [pa]

    assert not list(paths_containing_key("w", paths=[pa, pb]))
    assert not list(paths_containing_key("x", paths=[]))


def test_collect_yaml_dir(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"b": 2}}

    expected = {"x": 2, "y": {"a": 1, "b": 2}, "z": 3}

    pa, pb = tmp_path / "a.yaml", tmp_path / "b.yaml"
    pa.write_text(yaml.dump(a))
    pb.write_text(yaml.dump(b))

    config = merge(*collect_yaml(paths=[tmp_path]))
    assert config == expected


@contextmanager
def no_read_permissions(path):
    perm_orig = stat.S_IMODE(os.stat(path).st_mode)
    perm_new = perm_orig ^ stat.S_IREAD
    try:
        os.chmod(path, perm_new)
        yield
    finally:
        os.chmod(path, perm_orig)


@pytest.mark.skipif(
    sys.platform == "win32", reason="Can't make writeonly file on windows"
)
@pytest.mark.parametrize("kind", ["directory", "file"])
def test_collect_yaml_permission_errors(tmp_path, kind):
    a = {"x": 1, "y": 2}
    b = {"y": 3, "z": 4}

    pa, pb = tmp_path / "a.yaml", tmp_path / "b.yaml"
    pa.write_text(yaml.dump(a))
    pb.write_text(yaml.dump(b))

    if kind == "directory":
        cant_read = tmp_path
        expected = {}
    else:
        cant_read = pa
        expected = b

    with no_read_permissions(cant_read):
        config = merge(*collect_yaml(paths=[tmp_path]))
        assert config == expected


def test_collect_yaml_malformed_file(tmp_path):
    (tmp_path / "a.yaml").write_bytes(b"{")

    with pytest.raises(ValueError) as rec:
        list(collect_yaml(paths=[tmp_path]))
    assert "a.yaml" in str(rec.value)
    assert "is malformed" in str(rec.value)
    assert "original error message" in str(rec.value)


def test_collect_yaml_no_top_level_dict(tmp_path):
    (tmp_path / "a.yaml").write_bytes(b"[1234]")

    with pytest.raises(ValueError) as rec:
        list(collect_yaml(paths=[tmp_path]))
    assert "a.yaml" in str(rec.value)
    assert "is malformed" in str(rec.value)
    assert "must have a dict" in str(rec.value)


def test_env():
    env = {
        "DASK_A_B": "123",
        "DASK_C": "True",
        "DASK_D": "hello",
        "DASK_E__X": "123",
        "DASK_E__Y": "456",
        "DASK_F": '[1, 2, "3"]',
        "DASK_G": "/not/parsable/as/literal",
        "FOO": "not included",
    }

    expected = {
        "a_b": 123,
        "c": True,
        "d": "hello",
        "e": {"x": 123, "y": 456},
        "f": [1, 2, "3"],
        "g": "/not/parsable/as/literal",
    }

    res = collect_env(env)
    assert res == expected


@pytest.mark.parametrize(
    "preproc", [lambda x: x, lambda x: x.lower(), lambda x: x.upper()]
)
@pytest.mark.parametrize(
    "v,out", [("None", None), ("Null", None), ("False", False), ("True", True)]
)
def test_env_special_values(preproc, v, out):
    env = {"DASK_A": preproc(v)}
    res = collect_env(env)
    assert res == {"a": out}


def test_collect(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 2, "z": 3, "y": {"b": 2}}
    env = {"DASK_W": "4"}

    expected = {"w": 4, "x": 2, "y": {"a": 1, "b": 2}, "z": 3}

    pa, pb = tmp_path / "a.yaml", tmp_path / "b.yaml"
    pa.write_text(yaml.dump(a))
    pb.write_text(yaml.dump(b))

    config = collect([pa, pb], env=env)
    assert config == expected


def test_collect_env_none(monkeypatch):
    monkeypatch.setenv("DASK_FOO", "bar")
    config = collect([])
    assert config.get("foo") == "bar"


def test_get():
    d = {"x": 1, "y": {"a": 2}}

    assert get("x", config=d) == 1
    assert get("y.a", config=d) == 2
    assert get("y.b", 123, config=d) == 123
    with pytest.raises(KeyError):
        get("y.b", config=d)


def test_ensure_file(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    b = {"x": 123}

    source = tmp_path / "source.yaml"
    dest = tmp_path / "dest"
    destination = tmp_path / "dest" / "source.yaml"
    source.write_text(yaml.dump(a))

    ensure_file(source=source, destination=dest, comment=False)

    result = yaml.safe_load(destination.read_text())
    assert result == a

    # don't overwrite old config files
    source.write_text(yaml.dump(b))

    ensure_file(source=source, destination=dest, comment=False)

    result = yaml.safe_load(destination.read_text())
    assert result == a

    os.remove(destination)

    # Write again, now with comments
    ensure_file(source=source, destination=dest, comment=True)

    text = destination.read_text()
    assert "123" in text

    result = yaml.safe_load(destination.read_text())
    assert not result


def test_set():
    with dask.config.set(abc=123):
        assert config["abc"] == 123
        with dask.config.set(abc=456):
            assert config["abc"] == 456
        assert config["abc"] == 123

    assert "abc" not in config

    with dask.config.set({"abc": 123}):
        assert config["abc"] == 123
    assert "abc" not in config

    with dask.config.set({"abc.x": 1, "abc.y": 2, "abc.z.a": 3}):
        assert config["abc"] == {"x": 1, "y": 2, "z": {"a": 3}}
    assert "abc" not in config

    d = {}
    dask.config.set({"abc.x": 123}, config=d)
    assert d["abc"]["x"] == 123

    # Respect previous hyphenation, if any, or new hyphenation, if previous is not found
    d = {"e_f": 0, "g-h": 1}
    dask.config.set({"a_b": 2, "c-d": 3, "e-f": 4, "g_h": 5}, config=d)
    assert d == {"a_b": 2, "c-d": 3, "e_f": 4, "g-h": 5}


def test_set_kwargs():
    with dask.config.set(foo__bar=1, foo__baz=2):
        assert config["foo"] == {"bar": 1, "baz": 2}
    assert "foo" not in config

    # Mix kwargs and dict, kwargs override
    with dask.config.set({"foo.bar": 1, "foo.baz": 2}, foo__buzz=3, foo__bar=4):
        assert config["foo"] == {"bar": 4, "baz": 2, "buzz": 3}
    assert "foo" not in config

    # Mix kwargs and nested dict, kwargs override
    with dask.config.set({"foo": {"bar": 1, "baz": 2}}, foo__buzz=3, foo__bar=4):
        assert config["foo"] == {"bar": 4, "baz": 2, "buzz": 3}
    assert "foo" not in config


def test_set_nested():
    with dask.config.set({"abc": {"x": 123}}):
        assert config["abc"] == {"x": 123}
        with dask.config.set({"abc.y": 456}):
            assert config["abc"] == {"x": 123, "y": 456}
        assert config["abc"] == {"x": 123}
    assert "abc" not in config


def test_set_hard_to_copyables():
    import threading

    with dask.config.set(x=threading.Lock()):
        with dask.config.set(y=1):
            pass


@pytest.mark.parametrize("mkdir", [True, False])
def test_ensure_file_directory(mkdir, tmp_path):
    a = {"x": 1, "y": {"a": 1}}

    source = tmp_path / "source.yaml"
    dest = tmp_path / "dest"
    source.write_text(yaml.dump(a))

    if mkdir:
        dest.mkdir()

    ensure_file(source=source, destination=dest)

    assert dest.is_dir()
    assert (dest / "source.yaml").exists()


def test_ensure_file_defaults_to_DASK_CONFIG_directory(tmp_path):
    a = {"x": 1, "y": {"a": 1}}
    source = tmp_path / "source.yaml"
    source.write_text(yaml.dump(a))

    destination = tmp_path / "dask"
    PATH = dask.config.PATH
    try:
        dask.config.PATH = destination
        ensure_file(source=source)
    finally:
        dask.config.PATH = PATH

    assert destination.is_dir()
    [fn] = os.listdir(destination)
    assert os.path.split(fn)[1] == os.path.split(source)[1]


def test_pop():
    config = {"foo": {"ba-r": 1, "baz": 2}, "asd": 3}
    with pytest.raises(KeyError):
        pop("x", config=config)
    assert pop("x", default=4, config=config) == 4
    assert pop("foo.ba_r", config=config) == 1
    assert pop("asd", config=config) == 3
    assert config == {"foo": {"baz": 2}}


def test_refresh():
    defaults = []
    config = {}

    update_defaults({"a": 1}, config=config, defaults=defaults)
    assert config == {"a": 1}

    refresh(paths=[], env={"DASK_B": "2"}, config=config, defaults=defaults)
    assert config == {"a": 1, "b": 2}

    refresh(paths=[], env={"DASK_C": "3"}, config=config, defaults=defaults)
    assert config == {"a": 1, "c": 3}


@pytest.mark.parametrize(
    "inp,out",
    [
        ("1", "1"),
        (1, 1),
        ("$FOO", "foo"),
        ([1, "$FOO"], [1, "foo"]),
        ((1, "$FOO"), (1, "foo")),
        ({1, "$FOO"}, {1, "foo"}),
        ({"a": "$FOO"}, {"a": "foo"}),
        ({"a": "A", "b": [1, "2", "$FOO"]}, {"a": "A", "b": [1, "2", "foo"]}),
    ],
)
def test_expand_environment_variables(monkeypatch, inp, out):
    monkeypatch.setenv("FOO", "foo")
    assert expand_environment_variables(inp) == out


def test_env_var_canonical_name(monkeypatch):
    value = 3
    monkeypatch.setenv("DASK_A_B", str(value))
    d = {}
    dask.config.refresh(config=d)
    assert get("a_b", config=d) == value
    assert get("a-b", config=d) == value


def test_get_set_canonical_name():
    c = {"x-y": {"a_b": 123}}

    keys = ["x_y.a_b", "x-y.a-b", "x_y.a-b"]
    for k in keys:
        assert dask.config.get(k, config=c) == 123

    with dask.config.set({"x_y": {"a-b": 456}}, config=c):
        for k in keys:
            assert dask.config.get(k, config=c) == 456

    # No change to new keys in sub dicts
    with dask.config.set({"x_y": {"a-b": {"c_d": 1}, "e-f": 2}}, config=c):
        assert dask.config.get("x_y.a-b", config=c) == {"c_d": 1}
        assert dask.config.get("x_y.e_f", config=c) == 2


@pytest.mark.parametrize("key", ["custom_key", "custom-key"])
def test_get_set_roundtrip(key):
    value = 123
    with dask.config.set({key: value}):
        assert dask.config.get("custom_key") == value
        assert dask.config.get("custom-key") == value


def test_merge_None_to_dict():
    assert dask.config.merge({"a": None, "c": 0}, {"a": {"b": 1}}) == {
        "a": {"b": 1},
        "c": 0,
    }


def test_core_file():
    assert "temporary-directory" in dask.config.config
    assert "dataframe" in dask.config.config
    assert "compression" in dask.config.get("dataframe.shuffle")


def test_schema():
    jsonschema = pytest.importorskip("jsonschema")

    root_dir = pathlib.Path(__file__).parent.parent
    config = yaml.safe_load((root_dir / "dask.yaml").read_text())
    schema = yaml.safe_load((root_dir / "dask-schema.yaml").read_text())

    jsonschema.validate(config, schema)


def test_schema_is_complete():
    root_dir = pathlib.Path(__file__).parent.parent
    config = yaml.safe_load((root_dir / "dask.yaml").read_text())
    schema = yaml.safe_load((root_dir / "dask-schema.yaml").read_text())

    def test_matches(c, s):
        for k, v in c.items():
            if list(c) != list(s["properties"]):
                raise ValueError(
                    "\nThe dask.yaml and dask-schema.yaml files are not in sync.\n"
                    "This usually happens when we add a new configuration value,\n"
                    "but don't add the schema of that value to the dask-schema.yaml file\n"
                    "Please modify these files to include the missing values: \n\n"
                    "    dask.yaml:        {}\n"
                    "    dask-schema.yaml: {}\n\n"
                    "Examples in these files should be a good start, \n"
                    "even if you are not familiar with the jsonschema spec".format(
                        sorted(c), sorted(s["properties"])
                    )
                )
            if isinstance(v, dict):
                test_matches(c[k], s["properties"][k])

    test_matches(config, schema)


def test_rename():
    aliases = {
        "foo-bar": "foo.bar",
        "x.y": "foo.y",
        "a.b": "ab",
        "not-found": "not.found",
        "super.old": "super",
    }
    config = {
        "foo_bar": 1,
        "x": {"y": 2, "z": 3},
        "a": {"b": None},
        "super": {"old": 4},
    }
    with pytest.warns(FutureWarning) as w:
        rename(aliases, config=config)
    assert [str(wi.message) for wi in w.list] == [
        "Dask configuration key 'foo-bar' has been deprecated; please use 'foo.bar' instead",
        "Dask configuration key 'x.y' has been deprecated; please use 'foo.y' instead",
        "Dask configuration key 'a.b' has been deprecated; please use 'ab' instead",
        "Dask configuration key 'super.old' has been deprecated; please use 'super' instead",
    ]
    assert config == {
        "foo": {"bar": 1, "y": 2},
        "x": {"z": 3},
        "a": {},
        "ab": None,
        "super": 4,
    }


@pytest.mark.parametrize(
    "args,kwargs",
    [
        ((), {"fuse_ave_width": 123}),
        (({"fuse_ave_width": 123},), {}),
        (({"fuse-ave-width": 123},), {}),
    ],
)
def test_deprecations_on_set(args, kwargs):
    with pytest.warns(FutureWarning) as info:
        with dask.config.set(*args, **kwargs):
            assert dask.config.get("optimization.fuse.ave-width") == 123
    assert "optimization.fuse.ave-width" in str(info[0].message)


def test_deprecations_on_env_variables(monkeypatch):
    d = {}
    monkeypatch.setenv("DASK_FUSE_AVE_WIDTH", "123")
    with pytest.warns(FutureWarning) as info:
        dask.config.refresh(config=d)
    assert "optimization.fuse.ave-width" in str(info[0].message)
    assert get("optimization.fuse.ave-width", config=d) == 123


@pytest.mark.parametrize("key", ["fuse-ave-width", "fuse_ave_width"])
def test_deprecations_on_yaml(tmp_path, key):
    d = {}
    (tmp_path / "dask.yaml").write_text(yaml.dump({key: 123}))

    with pytest.warns(FutureWarning) as info:
        dask.config.refresh(config=d, paths=[tmp_path])
    assert "optimization.fuse.ave-width" in str(info[0].message)
    assert get("optimization.fuse.ave-width", config=d) == 123


def test_get_override_with():
    with dask.config.set({"foo": "bar"}):
        # If override_with is None get the config key
        assert dask.config.get("foo") == "bar"
        assert dask.config.get("foo", override_with=None) == "bar"

        # Otherwise pass the default straight through
        assert dask.config.get("foo", override_with="baz") == "baz"
        assert dask.config.get("foo", override_with=False) is False
        assert dask.config.get("foo", override_with=True) is True
        assert dask.config.get("foo", override_with=123) == 123
        assert dask.config.get("foo", override_with={"hello": "world"}) == {
            "hello": "world"
        }
        assert dask.config.get("foo", override_with=["one"]) == ["one"]


def test_config_serialization():
    # Use context manager without changing the value to ensure test side effects are restored
    with dask.config.set({"array.svg.size": dask.config.get("array.svg.size")}):
        # Take a round trip through the serialization
        serialized = serialize({"array": {"svg": {"size": 150}}})
        config = deserialize(serialized)

        dask.config.update(dask.config.global_config, config)
        assert dask.config.get("array.svg.size") == 150


def test_config_inheritance():
    config = collect_env(
        {"DASK_INTERNAL_INHERIT_CONFIG": serialize({"array": {"svg": {"size": 150}}})}
    )
    assert dask.config.get("array.svg.size", config=config) == 150


def test__get_paths(monkeypatch):
    # These settings are used by Dask's config system. We temporarily
    # remove them to avoid interference from the machine where tests
    # are being run.
    monkeypatch.delenv("DASK_CONFIG", raising=False)
    monkeypatch.delenv("DASK_ROOT_CONFIG", raising=False)
    monkeypatch.setattr(site, "PREFIXES", [])

    expected = [
        "/etc/dask",
        os.path.join(sys.prefix, "etc", "dask"),
        os.path.join(os.path.expanduser("~"), ".config", "dask"),
    ]
    paths = _get_paths()
    assert paths == expected
    assert len(paths) == len(set(paths))  # No duplicate paths

    with monkeypatch.context() as m:
        m.setenv("DASK_CONFIG", "foo-bar")
        paths = _get_paths()
        assert paths == expected + ["foo-bar"]
        assert len(paths) == len(set(paths))

    with monkeypatch.context() as m:
        m.setenv("DASK_ROOT_CONFIG", "foo-bar")
        paths = _get_paths()
        assert paths == ["foo-bar"] + expected[1:]
        assert len(paths) == len(set(paths))

    with monkeypatch.context() as m:
        prefix = os.path.join("include", "this", "path")
        m.setattr(site, "PREFIXES", site.PREFIXES + [prefix])
        paths = _get_paths()
        assert os.path.join(prefix, "etc", "dask") in paths
        assert len(paths) == len(set(paths))


def test_default_search_paths():
    # Ensure _get_paths() is used for default paths
    assert dask.config.paths == _get_paths()
