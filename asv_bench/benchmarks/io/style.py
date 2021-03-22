import numpy as np

from pandas import DataFrame


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
        self.st.render()

    def peakmem_apply_render(self, cols, rows):
        self._style_apply()
        self.st.render()

    def time_classes_render(self, cols, rows):
        self._style_classes()
        self.st.render()

    def peakmem_classes_render(self, cols, rows):
        self._style_classes()
        self.st.render()

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
