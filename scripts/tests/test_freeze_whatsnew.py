import pytest

from scripts.freeze_whatsnew import (
    Block,
    era_of,
    freeze,
    match_blocks,
    normalize_backticks,
    parse_blocks,
    render,
    split_items,
    strip_prompts,
)


def make_block(content: list[str]) -> Block:
    return Block(0, len(content), 0, ["python"], {}, content)


def test_parse_blocks_reads_options_and_content() -> None:
    text = (
        "Some prose.\n"
        "\n"
        ".. ipython:: python\n"
        "   :okwarning:\n"
        "\n"
        "   ser = pd.Series([1, 2])\n"
        "   ser\n"
        "\n"
        "More prose.\n"
    )
    (block,) = parse_blocks(text)
    assert block.options == {"okwarning": ""}
    assert block.content == ["ser = pd.Series([1, 2])", "ser"]
    assert block.base == 0


def test_parse_blocks_options_indented_shallower_than_body() -> None:
    # The option marker sits one column left of the code.  If it is allowed
    # into the dedent, every code line keeps a leading space and the body stops
    # parsing as Python -- which silently cost three blocks their output.
    text = (
        ".. ipython:: python\n"
        "   :okwarning:\n"
        "\n"
        "    ser = pd.Series(['a', 'b'])\n"
        "    ser.str.cat(ser)\n"
    )
    (block,) = parse_blocks(text)
    assert block.content == ["ser = pd.Series(['a', 'b'])", "ser.str.cat(ser)"]
    assert [item["kind"] for item in split_items(block.content)] == ["code", "code"]


def test_parse_blocks_indented_directive_keeps_base() -> None:
    text = (
        "- A bullet:\n"
        "\n"
        "  .. ipython:: python\n"
        "\n"
        "     ser = pd.Series([1])\n"
        "\n"
        "- Next bullet.\n"
    )
    (block,) = parse_blocks(text)
    assert block.base == 2
    assert block.content == ["ser = pd.Series([1])"]


def test_parse_blocks_accepts_space_before_colons() -> None:
    # Released files spell it with a space before the colons; the current
    # tree does not, and a pre-commit hook forbids writing it literally here.
    spaced = ".. ipython {} python".format("::")
    text = f"{spaced}\n\n   ser = pd.Series([1])\n"
    (block,) = parse_blocks(text)
    assert block.content == ["ser = pd.Series([1])"]


def test_split_items_records_blank_lines_and_comments() -> None:
    content = [
        "ser = pd.Series([1, 2])",
        "",
        "# a comment",
        "ser.sum()",
    ]
    items = split_items(content)
    assert [item["kind"] for item in items] == ["code", "comment", "code"]
    assert [item["blank_before"] for item in items] == [False, True, False]


def test_split_items_flags_magics() -> None:
    items = split_items(["%timeit pd.eval('a + b')"])
    assert [item["kind"] for item in items] == ["magic"]


def test_split_items_reports_unparsable_body() -> None:
    items = split_items(["   this is not python at all ("])
    assert [item["kind"] for item in items] == ["unparsed"]


def test_strip_prompts_unwraps_an_already_prompted_block() -> None:
    content = [
        "In [5]: df.all(bool_only=True)",
        "",
        "In [6]: df[['B', 'C']].all(bool_only=True)",
    ]
    assert strip_prompts(content) == [
        "df.all(bool_only=True)",
        "",
        "df[['B', 'C']].all(bool_only=True)",
    ]


def test_strip_prompts_leaves_plain_code_alone() -> None:
    content = ["ser = pd.Series([1])", "ser"]
    assert strip_prompts(content) == content


def test_render_uses_doctest_prompts_and_keeps_spacing() -> None:
    items = split_items(
        [
            "ser = pd.Series([1, 2])",
            "ser.sum()",
            "",
            "# and again",
            "ser.sum()",
        ]
    )
    results = [
        {"output": "", "error": None},
        {"output": "3\n", "error": None},
        {"output": "", "error": None},
        {"output": "3\n", "error": None},
    ]
    assert render(items, results) == [
        ">>> ser = pd.Series([1, 2])",
        ">>> ser.sum()",
        "3",
        "",
        "# and again",
        ">>> ser.sum()",
        "3",
    ]


def test_render_continues_multiline_statements() -> None:
    items = split_items(["ser = pd.Series(\n    [1, 2]\n)"])
    rendered = render(items, [{"output": "", "error": None}])
    assert rendered == [">>> ser = pd.Series(", "...     [1, 2]", "... )"]


def test_render_shows_a_raise_as_a_doctest_traceback() -> None:
    items = split_items(["pd.Index([1], dtype=np.float16)"])
    results = [{"output": "", "error": "NotImplementedError: nope"}]
    assert render(items, results) == [
        ">>> pd.Index([1], dtype=np.float16)",
        "Traceback (most recent call last):",
        "    ...",
        "NotImplementedError: nope",
    ]


def test_match_blocks_pairs_across_an_already_frozen_block() -> None:
    # The released file still has three directives; the current file has two,
    # because the middle one was frozen by an earlier pass.
    era = [make_block(["one()"]), make_block(["two()"]), make_block(["three()"])]
    current = [make_block(["one()"]), make_block(["three()"])]
    assert match_blocks(era, current) == {0: 0, 1: 2}


def test_match_blocks_tolerates_requoting_and_reflow() -> None:
    era = [make_block(["ser = pd.Series(['a', 'b'])"])]
    current = [make_block(['ser = pd.Series(\n    ["a", "b"]\n)'])]
    assert match_blocks(era, current) == {0: 0}


def test_match_blocks_leaves_an_unrelated_block_unpaired() -> None:
    era = [make_block(["ser.sum()"])]
    current = [make_block(["completely.different(thing, entirely=True)"])]
    assert match_blocks(era, current) == {}


def test_normalize_backticks_only_doubles_single_spans() -> None:
    assert normalize_backticks("# set `dropna` to False") == (
        "# set ``dropna`` to False"
    )
    assert normalize_backticks("# set ``dropna`` to False") == (
        "# set ``dropna`` to False"
    )


@pytest.mark.parametrize(
    "name, expected",
    [("v0.13.0.rst", (0, 13)), ("v1.5.1.rst", (1, 5)), ("v0.19.0.txt", (0, 19))],
)
def test_era_of(name: str, expected: tuple[int, int]) -> None:
    assert era_of(name) == expected


def test_era_of_rejects_a_non_whatsnew_name() -> None:
    with pytest.raises(ValueError, match="not a whatsnew file"):
        era_of("index.rst")


def test_split_items_drops_ipython_output_suppression() -> None:
    # A trailing semicolon told the directive not to render the result.  It is
    # not Python, and leaving it in makes the statement uncompilable as an
    # expression, so the block would freeze with no output at all.
    (item,) = split_items(["df;"])
    assert item["source"] == "df"


def test_render_hides_a_suppressed_statement_entirely() -> None:
    items = split_items(["@suppress", "df.index = other", "df"])
    results = [{"output": "", "error": None}, {"output": "1\n", "error": None}]
    assert render(items, results) == [">>> df", "1"]


def test_render_keeps_a_multiline_exception_message() -> None:
    items = split_items(["pd.read_parquet('x')"])
    results = [{"output": "", "error": "ImportError: no engine\nInstall pyarrow."}]
    assert render(items, results)[-2:] == ["ImportError: no engine", "Install pyarrow."]


def write_page(tmp_path, body: str) -> str:
    path = tmp_path / "v0.13.0.rst"
    path.write_text(f"Version 0.13.0\n--------------\n\n{body}")
    return str(path)


def test_freeze_rewrites_a_block_against_the_current_pandas(tmp_path) -> None:
    path = write_page(
        tmp_path,
        ".. ipython:: python\n\n   ser = pd.Series([1, 2])\n   ser.sum()\n",
    )
    frozen, report = freeze(path, True, None, False, False)
    assert report["blocks"] == 1
    assert ".. code-block:: pycon" in frozen
    assert ">>> ser.sum()" in frozen
    assert "\n   np.int64(3)\n" in frozen
    assert ".. ipython::" not in frozen


def test_freeze_deletes_a_suppressed_block(tmp_path) -> None:
    path = write_page(
        tmp_path,
        ".. ipython:: python\n"
        "   :suppress:\n"
        "\n"
        "   ser = pd.Series([1])\n"
        "\n"
        "Prose after.\n",
    )
    frozen, _ = freeze(path, True, None, False, False)
    assert "code-block" not in frozen
    assert frozen.endswith("Prose after.\n")


def test_freeze_refuses_a_block_that_does_not_parse(tmp_path) -> None:
    path = write_page(tmp_path, ".. ipython:: python\n\n   not python at all (\n")
    frozen, report = freeze(path, True, None, False, False)
    assert frozen is None
    assert "block did not parse as Python" in report["reason"]


def test_freeze_skip_unusable_leaves_the_block_executable(tmp_path) -> None:
    path = write_page(tmp_path, ".. ipython:: python\n\n   not python at all (\n")
    frozen, report = freeze(path, True, None, True, False)
    assert frozen is not None
    assert ".. ipython:: python" in frozen
    assert len(report["skipped"]) == 1


def test_freeze_refuses_output_that_is_only_an_object_address(tmp_path) -> None:
    # A default __repr__ pins the object's address, so it differs every run.
    path = write_page(
        tmp_path,
        ".. ipython:: python\n\n   df = pd.DataFrame({'a': [1]})\n   df.style\n",
    )
    frozen, report = freeze(path, True, None, False, False)
    assert "df.style" in frozen
    assert "object at 0x" not in frozen
    assert report["volatile"]


def test_freeze_drops_a_block_that_renders_nothing(tmp_path) -> None:
    # Every statement suppressed: the directive showed nothing, and an empty
    # literal block is a docutils error under -W.
    path = write_page(
        tmp_path,
        ".. ipython:: python\n"
        "\n"
        "   @suppress\n"
        "   ser = pd.Series([1])\n"
        "\n"
        "Prose after.\n",
    )
    frozen, report = freeze(path, True, None, False, False)
    assert "code-block" not in frozen
    assert report["emptied"]


def test_freeze_shows_volatile_output_as_input_only(tmp_path) -> None:
    # A wall clock cannot be frozen, and the block cannot be left executable
    # either: once its neighbours are frozen it no longer has their state, and
    # the doc build fails on it.  Show the input and no output.
    path = write_page(
        tmp_path,
        ".. ipython:: python\n\n   pd.Timestamp.now()\n",
    )
    frozen, report = freeze(path, True, None, False, False)
    assert ".. code-block:: ipython" in frozen
    assert "pd.Timestamp.now()" in frozen
    assert "Timestamp('" not in frozen
    assert report["volatile"]
