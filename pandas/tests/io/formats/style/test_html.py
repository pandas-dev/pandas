import pytest

from pandas import DataFrame

jinja2 = pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


class TestStylerHTML:
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
      <th class="blank level0" >&nbsp;</th>
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
        result = s.to_html
        return result
