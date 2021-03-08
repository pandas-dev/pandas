import textwrap

import pytest

from pandas import DataFrame

jinja2 = pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


class TestStylerHTML:
    def setup(self):
        loader = jinja2.PackageLoader("pandas", "io/formats/templates")
        env = jinja2.Environment(loader=loader, trim_blocks=True)
        self.styles_template = env.get_template("html_all.tpl")
        self.no_styles_template = env.get_template("html_excl.tpl")

    def test_html_template_extends_options(self):
        # make sure if templates are edited tests are updated as are setup fixtures
        with open("pandas/io/formats/templates/html.tpl") as file:
            result = file.read()
        expected = (
            '{% if no_styles %}\n{% extends "html_excl.tpl" %}\n'
            '{% else %}\n{% extends "html_all.tpl" %}\n'
        )
        assert expected in result

    def test_w3_html_format(self):
        s = (
            Styler(
                DataFrame([[2.61], [2.69]], index=["a", "b"], columns=["A"]),
                uuid_len=0,
            )
            .set_table_styles([{"selector": "th", "props": "att2:v2;"}])
            .applymap(lambda x: "att1:v1;")
            .set_table_attributes('class="my-cls1" style="attr3:v3;"')
            .set_td_classes(DataFrame(["my-cls2"], index=["a"], columns=["A"]))
            .format("{:.1f}")
            .set_caption("A comprehensive test")
        )
        expected = """<style type="text/css">
#T__ th {
  att2: v2;
}
#T__row0_col0, #T__row1_col0 {
  att1: v1;
}
</style>
<table id="T__" class="my-cls1" style="attr3:v3;">
  <caption>A comprehensive test</caption>
  <thead>
    <tr>
      <th class="blank level0" ></th>
      <th class="col_heading level0 col0" >A</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0" >a</th>
      <td id="T__row0_col0" class="data row0 col0 my-cls2" >2.6</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1" >b</th>
      <td id="T__row1_col0" class="data row1 col0" >2.7</td>
    </tr>
  </tbody>
</table>
"""
        assert expected == s.render()

    def test_no_styles(self):
        s = Styler(DataFrame([[2.61], [2.69]], index=["a", "b"], columns=["A"]))
        result = s.to_html(no_styles=True)
        expected = """<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
</head>
<body>
<table>
  <thead>
    <tr>
      <th ></th>
      <th >A</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th >a</th>
      <td >2.610000</td>
    </tr>
    <tr>
      <th >b</th>
      <td >2.690000</td>
    </tr>
  </tbody>
</table>
</body>
</html>
"""
        assert result == expected

    def test_styles(self):
        s = Styler(
            DataFrame([[2.61], [2.69]], index=["a", "b"], columns=["A"]), uuid="abc"
        )
        result = s.to_html()
        expected = """<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <style type="text/css">
  </style>
</head>
<body>
<table id="T_abc_">
  <thead>
    <tr>
      <th class="blank level0" ></th>
      <th class="col_heading level0 col0" >A</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T_abc_level0_row0" class="row_heading level0 row0" >a</th>
      <td id="T_abc_row0_col0" class="data row0 col0" >2.610000</td>
    </tr>
    <tr>
      <th id="T_abc_level0_row1" class="row_heading level0 row1" >b</th>
      <td id="T_abc_row1_col0" class="data row1 col0" >2.690000</td>
    </tr>
  </tbody>
</table>
</body>
</html>
"""
        assert result == expected

    def test_block_names(self):
        # catch accidental removal of a block
        expected1 = {
            "before_style",
            "style",
            "table_styles",
            "before_cellstyle",
            "cellstyle",
            "before_table",
            "table",
            "caption",
            "thead",
            "tbody",
            "after_table",
            "before_head_rows",
            "head_tr",
            "after_head_rows",
            "before_rows",
            "tr",
            "after_rows",
        }
        result1 = set(self.styles_template.blocks)
        assert result1 == expected1

        expected2 = {v for v in expected1 if "style" not in v}
        result2 = set(self.no_styles_template.blocks)
        assert result2 == expected2


def test_from_custom_template(tmpdir):
    p = tmpdir.mkdir("templates").join("myhtml.tpl")
    p.write(
        textwrap.dedent(
            """\
        {% extends "html.tpl" %}
        {% block table %}
        <h1>{{ table_title|default("My Table") }}</h1>
        {{ super() }}
        {% endblock table %}"""
        )
    )
    result = Styler.from_custom_template(str(tmpdir.join("templates")), "myhtml.tpl")
    assert issubclass(result, Styler)
    assert result.env is not Styler.env
    assert result.template is not Styler.template
    styler = result(DataFrame({"A": [1, 2]}))
    assert styler.render()
