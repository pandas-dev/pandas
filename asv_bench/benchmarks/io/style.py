import numpy as np

from pandas import DataFrame


class RenderApply:

    params = [[12, 24, 36], [12, 120]]
    param_names = ["cols", "rows"]

    def setup(self, cols, rows):
        self.df = DataFrame(
            np.random.randn(rows, cols),
            columns=[f"float_{i+1}" for i in range(cols)],
            index=[f"row_{i+1}" for i in range(rows)],
        )
        self._style_apply()

    def time_render(self, cols, rows):
        self.st.render()

    def peakmem_apply(self, cols, rows):
        self._style_apply()

    def peakmem_render(self, cols, rows):
        self.st.render()

    def _style_apply(self):
        def _apply_func(s):
            return [
                "background-color: lightcyan" if s.name == "row_1" else "" for v in s
            ]

        self.st = self.df.style.apply(_apply_func, axis=1)
