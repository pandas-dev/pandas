import numpy as np

from pandas import (
    DataFrame,
    IndexSlice,
)


class Render:

    params = [[12, 24, 36], [12, 120]]
    param_names = ["cols", "rows"]

    def setup(self, cols, rows):
        self.df = DataFrame(
            np.random.randn(rows, cols),
            columns=[f"float_{i+1}" for i in range(cols)],
            index=[f"row_{i+1}" for i in range(rows)],
        )

    def time_apply_render(self, cols, rows):
        self._style_apply()
        self.st._render_html(True, True)

    def peakmem_apply_render(self, cols, rows):
        self._style_apply()
        self.st._render_html(True, True)

    def time_classes_render(self, cols, rows):
        self._style_classes()
        self.st._render_html(True, True)

    def peakmem_classes_render(self, cols, rows):
        self._style_classes()
        self.st._render_html(True, True)

    def time_tooltips_render(self, cols, rows):
        self._style_tooltips()
        self.st._render_html(True, True)

    def peakmem_tooltips_render(self, cols, rows):
        self._style_tooltips()
        self.st._render_html(True, True)

    def time_format_render(self, cols, rows):
        self._style_format()
        self.st._render_html(True, True)

    def peakmem_format_render(self, cols, rows):
        self._style_format()
        self.st._render_html(True, True)

    def time_apply_format_hide_render(self, cols, rows):
        self._style_apply_format_hide()
        self.st._render_html(True, True)

    def peakmem_apply_format_hide_render(self, cols, rows):
        self._style_apply_format_hide()
        self.st._render_html(True, True)

    def _style_apply(self):
        def _apply_func(s):
            return [
                "background-color: lightcyan" if s.name == "row_1" else "" for v in s
            ]

        self.st = self.df.style.apply(_apply_func, axis=1)

    def _style_classes(self):
        classes = self.df.applymap(lambda v: ("cls-1" if v > 0 else ""))
        classes.index, classes.columns = self.df.index, self.df.columns
        self.st = self.df.style.set_td_classes(classes)

    def _style_format(self):
        ic = int(len(self.df.columns) / 4 * 3)
        ir = int(len(self.df.index) / 4 * 3)
        # apply a formatting function
        # subset is flexible but hinders vectorised solutions
        self.st = self.df.style.format(
            "{:,.3f}", subset=IndexSlice["row_1":f"row_{ir}", "float_1":f"float_{ic}"]
        )

    def _style_apply_format_hide(self):
        self.st = self.df.style.applymap(lambda v: "color: red;")
        self.st.format("{:.3f}")
        self.st.hide(self.st.index[1:], axis=0)
        self.st.hide(self.st.columns[1:], axis=1)

    def _style_tooltips(self):
        ttips = DataFrame("abc", index=self.df.index[::2], columns=self.df.columns[::2])
        self.st = self.df.style.set_tooltips(ttips)
        self.st.hide(self.st.index[12:], axis=0)
        self.st.hide(self.st.columns[12:], axis=1)
