#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pandas.io.templating import ITemplateEngine


class Jinja2Engine(ITemplateEngine):
    def __init__(self):
        from jinja2 import Template

    def load(self, s):
        from jinja2 import Template

        self.t = Template(s)

    def render(self, ctx):
        if self.t:
            return self.t.render(**ctx)
        else:
            raise AssertionError("Template not initialized, cannot render")