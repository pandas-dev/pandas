
from io import StringIO
from html.parser import HTMLParser
import json
import os
import re
import tempfile
import shutil

import traitlets

from ..widgets import IntSlider, IntText, Text, Widget, jslink, HBox, widget_serialization, widget as widget_module
from ..embed import embed_data, embed_snippet, embed_minimal_html, dependency_state


class CaseWidget(Widget):
    """Widget to test dependency traversal"""

    a = traitlets.Instance(Widget, allow_none=True).tag(sync=True, **widget_serialization)
    b = traitlets.Instance(Widget, allow_none=True).tag(sync=True, **widget_serialization)

    _model_name = traitlets.Unicode('CaseWidgetModel').tag(sync=True)

    other = traitlets.Dict().tag(sync=True, **widget_serialization)




class TestEmbed:

    def teardown(self):
        for w in tuple(widget_module._instances.values()):
            w.close()

    def test_embed_data_simple(self):
        w = IntText(4)
        state = dependency_state(w, drop_defaults=True)
        data = embed_data(views=w, drop_defaults=True, state=state)

        state = data['manager_state']['state']
        views = data['view_specs']

        assert len(state) == 3
        assert len(views) == 1

        model_names = [s['model_name'] for s in state.values()]
        assert 'IntTextModel' in model_names

    def test_cors(self):
        w = IntText(4)
        code = embed_snippet(w)
         # 1 is from the require
        assert len(re.findall(' crossorigin', code)) > 1
        f = StringIO()
        embed_minimal_html(f, w)
        assert len(re.findall(' crossorigin', f.getvalue())) > 1

        code = embed_snippet(w, cors=False, requirejs=False)
        assert ' crossorigin' not in code
        f = StringIO()
        embed_minimal_html(f, w, cors=False, requirejs=False)
        assert ' crossorigin' not in f.getvalue()

        code = embed_snippet(w, cors=False, requirejs=True)
        assert len(re.findall(' crossorigin', code)) == 1 # 1 is from the require, which is ok
        f = StringIO()
        embed_minimal_html(f, w, cors=False, requirejs=True)
        assert len(re.findall(' crossorigin', f.getvalue())) == 1 # 1 is from the require, which is ok

    def test_escape(self):
        w = Text('<script A> <ScRipt> </Script> <!-- --> <b>hi</b>')
        code = embed_snippet(w)
        assert code.find(r'<script A>') == -1
        assert code.find(r'\u003cscript A> \u003cScRipt> \u003c/Script> \u003c!-- --> <b>hi</b>') >= 0

        f = StringIO()
        embed_minimal_html(f, w)
        content = f.getvalue()
        assert content.find(r'<script A>') == -1
        assert content.find(r'\u003cscript A> \u003cScRipt> \u003c/Script> \u003c!-- --> <b>hi</b>') >= 0

    def test_embed_data_two_widgets(self):
        w1 = IntText(4)
        w2 = IntSlider(min=0, max=100)
        jslink((w1, 'value'), (w2, 'value'))
        state = dependency_state([w1, w2], drop_defaults=True)
        data = embed_data(views=[w1, w2], drop_defaults=True, state=state)

        state = data['manager_state']['state']
        views = data['view_specs']

        assert len(state) == 7
        assert len(views) == 2

        model_names = [s['model_name'] for s in state.values()]
        assert 'IntTextModel' in model_names
        assert 'IntSliderModel' in model_names

    def test_embed_data_complex(self):
        w1 = IntText(4)
        w2 = IntSlider(min=0, max=100)
        jslink((w1, 'value'), (w2, 'value'))

        w3 = CaseWidget()
        w3.a = w1

        w4 = CaseWidget()
        w4.a = w3
        w4.other['test'] = w2

        # Add a circular reference:
        w3.b = w4

        # Put it in an HBox
        HBox(children=[w4])

        state = dependency_state(w3)

        assert len(state) == 9

        model_names = [s['model_name'] for s in state.values()]
        assert 'IntTextModel' in model_names
        assert 'IntSliderModel' in model_names
        assert 'CaseWidgetModel' in model_names
        assert 'LinkModel' in model_names

        # Check that HBox is not collected
        assert 'HBoxModel' not in model_names

        # Check that views make sense:

        data = embed_data(views=w3, drop_defaults=True, state=state)
        assert state is data['manager_state']['state']
        views = data['view_specs']
        assert len(views) == 1


    def test_snippet(self):

        class Parser(HTMLParser):
            state = 'initial'
            states = []

            def handle_starttag(self, tag, attrs):
                attrs = dict(attrs)
                if tag == 'script' and attrs.get('type', '') == "application/vnd.jupyter.widget-state+json":
                    self.state = 'widget-state'
                    self.states.append(self.state)
                elif tag == 'script' and attrs.get('type', '') == "application/vnd.jupyter.widget-view+json":
                    self.state = 'widget-view'
                    self.states.append(self.state)

            def handle_endtag(self, tag):
                self.state = 'initial'

            def handle_data(self, data):
                if self.state == 'widget-state':
                    manager_state = json.loads(data)['state']
                    assert len(manager_state) == 3
                    self.states.append('check-widget-state')
                elif self.state == 'widget-view':
                    view = json.loads(data)
                    assert isinstance(view, dict)
                    self.states.append('check-widget-view')

        w = IntText(4)
        state = dependency_state(w, drop_defaults=True)
        snippet = embed_snippet(views=w, drop_defaults=True, state=state)
        parser = Parser()
        parser.feed(snippet)
        print(parser.states)
        assert parser.states == ['widget-state', 'check-widget-state', 'widget-view', 'check-widget-view']

    def test_minimal_html_filename(self):
        w = IntText(4)

        tmpd = tempfile.mkdtemp()

        try:
            output = os.path.join(tmpd, 'test.html')
            state = dependency_state(w, drop_defaults=True)
            embed_minimal_html(output, views=w, drop_defaults=True, state=state)
            # Check that the file is written to the intended destination:
            with open(output, 'r') as f:
                content = f.read()
            assert content.splitlines()[0] == '<!DOCTYPE html>'
        finally:
            shutil.rmtree(tmpd)

    def test_minimal_html_filehandle(self):
        w = IntText(4)
        output = StringIO()
        state = dependency_state(w, drop_defaults=True)
        embed_minimal_html(output, views=w, drop_defaults=True, state=state)
        content = output.getvalue()
        assert content.splitlines()[0] == '<!DOCTYPE html>'
