import uuid
import os
import functools
from os.path import abspath, dirname

from pandas.io.templating import *


# TODO: caption

class HTMLStyler(Styler):
    def clone(self):
        """We always return a modified copy when extra styling is added

        `clone()` provides a copy to apply changes to
        """
        import copy

        c = HTMLStyler(self.df,
                       template=self.template,
                       engine_instance=self.engine_instance)
        c.style = copy.deepcopy(self.style)
        c.cell_context = copy.deepcopy(self.cell_context)
        return c

    def __init__(self, *args, **kwds):

        if not kwds.get('template'):
            tmpl_filename = os.path.join(TEMPLATE_DIR, "html")
            with open(tmpl_filename) as f:
                kwds['template'] = f.read()

        super(HTMLStyler, self).__init__(*args, **kwds)

        import html_helpers
        import types

        for n in dir(html_helpers):
            cand = getattr(html_helpers, n)
            if callable(cand):
                def f(cand):
                    @functools.wraps(cand)
                    def f(self, *args, **kwds):
                        self_copy = self.clone()
                        cand(self_copy, *args, **kwds)
                        # self.style.extend(sd.style)
                        # TODO, need to merge-with
                        # self.cell_context = sd.cell_context
                        return self_copy

                    return f

                setattr(self, n, types.MethodType(f(cand), self))

    def _repr_html_(self):
        return self.render()

        # def __dir__(self):
        #     import html_helpers
        #     return ["render","to_context"] + dir(html_helpers)

        # def __getattr__(self, key):
        #     import html_helpers
        #     return getattr(html_helpers,key)