import os
from nose import SkipTest

# this is a mess. Getting failures on a python 2.7 build with
# whenever we try to import jinja, whether it's installed or not.
# so we're explicitly skipping that one *before* we try to import
# jinja. We still need to export the imports as globals,
# since importing Styler tries to import jinja2.
job_name = os.environ.get('JOB_NAME', None)
if job_name == '27_slow_nnet_LOCALE':
    raise SkipTest("No jinja")
try:
    from pandas.core.style import Styler
except ImportError:
    raise SkipTest("No Jinja2")

import copy

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.util.testing import TestCase
import pandas.util.testing as tm


class TestStyler(TestCase):

    def setUp(self):
        np.random.seed(24)
        self.s = DataFrame({'A': np.random.permutation(range(6))})
        self.df = DataFrame({'A': [0, 1], 'B': np.random.randn(2)})
        self.f = lambda x: x
        self.g = lambda x: x

        def h(x, foo='bar'):
            return pd.Series(['color: %s' % foo], index=x.index, name=x.name)

        self.h = h
        self.styler = Styler(self.df)
        self.attrs = pd.DataFrame({'A': ['color: red', 'color: blue']})
        self.dataframes = [
            self.df,
            pd.DataFrame({'f': [1., 2.], 'o': ['a', 'b'],
                          'c': pd.Categorical(['a', 'b'])})
        ]

    def test_update_ctx(self):
        self.styler._update_ctx(self.attrs)
        expected = {(0, 0): ['color: red'],
                    (1, 0): ['color: blue']}
        self.assertEqual(self.styler.ctx, expected)

    def test_update_ctx_flatten_multi(self):
        attrs = DataFrame({"A": ['color: red; foo: bar',
                                 'color: blue; foo: baz']})
        self.styler._update_ctx(attrs)
        expected = {(0, 0): ['color: red', ' foo: bar'],
                    (1, 0): ['color: blue', ' foo: baz']}
        self.assertEqual(self.styler.ctx, expected)

    def test_update_ctx_flatten_multi_traliing_semi(self):
        attrs = DataFrame({"A": ['color: red; foo: bar;',
                                 'color: blue; foo: baz;']})
        self.styler._update_ctx(attrs)
        expected = {(0, 0): ['color: red', ' foo: bar'],
                    (1, 0): ['color: blue', ' foo: baz']}
        self.assertEqual(self.styler.ctx, expected)

    def test_copy(self):
        s2 = copy.copy(self.styler)
        self.assertTrue(self.styler is not s2)
        self.assertTrue(self.styler.ctx is s2.ctx)  # shallow
        self.assertTrue(self.styler._todo is s2._todo)

        self.styler._update_ctx(self.attrs)
        self.styler.highlight_max()
        self.assertEqual(self.styler.ctx, s2.ctx)
        self.assertEqual(self.styler._todo, s2._todo)

    def test_deepcopy(self):
        s2 = copy.deepcopy(self.styler)
        self.assertTrue(self.styler is not s2)
        self.assertTrue(self.styler.ctx is not s2.ctx)
        self.assertTrue(self.styler._todo is not s2._todo)

        self.styler._update_ctx(self.attrs)
        self.styler.highlight_max()
        self.assertNotEqual(self.styler.ctx, s2.ctx)
        self.assertEqual(s2._todo, [])
        self.assertNotEqual(self.styler._todo, s2._todo)

    def test_clear(self):
        s = self.df.style.highlight_max()._compute()
        self.assertTrue(len(s.ctx) > 0)
        self.assertTrue(len(s._todo) > 0)
        s.clear()
        self.assertTrue(len(s.ctx) == 0)
        self.assertTrue(len(s._todo) == 0)

    def test_render(self):
        df = pd.DataFrame({"A": [0, 1]})
        style = lambda x: pd.Series(["color: red", "color: blue"], name=x.name)
        s = Styler(df, uuid='AB').apply(style).apply(style, axis=1)
        s.render()
        # it worked?

    def test_render_double(self):
        df = pd.DataFrame({"A": [0, 1]})
        style = lambda x: pd.Series(["color: red; border: 1px",
                                     "color: blue; border: 2px"], name=x.name)
        s = Styler(df, uuid='AB').apply(style)
        s.render()
        # it worked?

    def test_set_properties(self):
        df = pd.DataFrame({"A": [0, 1]})
        result = df.style.set_properties(color='white',
                                         size='10px')._compute().ctx
        # order is deterministic
        v = ["color: white", "size: 10px"]
        expected = {(0, 0): v, (1, 0): v}
        self.assertEqual(result.keys(), expected.keys())
        for v1, v2 in zip(result.values(), expected.values()):
            self.assertEqual(sorted(v1), sorted(v2))

    def test_set_properties_subset(self):
        df = pd.DataFrame({'A': [0, 1]})
        result = df.style.set_properties(subset=pd.IndexSlice[0, 'A'],
                                         color='white')._compute().ctx
        expected = {(0, 0): ['color: white']}
        self.assertEqual(result, expected)

    def test_apply_axis(self):
        df = pd.DataFrame({'A': [0, 0], 'B': [1, 1]})
        f = lambda x: ['val: %s' % x.max() for v in x]
        result = df.style.apply(f, axis=1)
        self.assertEqual(len(result._todo), 1)
        self.assertEqual(len(result.ctx), 0)
        result._compute()
        expected = {(0, 0): ['val: 1'], (0, 1): ['val: 1'],
                    (1, 0): ['val: 1'], (1, 1): ['val: 1']}
        self.assertEqual(result.ctx, expected)

        result = df.style.apply(f, axis=0)
        expected = {(0, 0): ['val: 0'], (0, 1): ['val: 1'],
                    (1, 0): ['val: 0'], (1, 1): ['val: 1']}
        result._compute()
        self.assertEqual(result.ctx, expected)
        result = df.style.apply(f)  # default
        result._compute()
        self.assertEqual(result.ctx, expected)

    def test_apply_subset(self):
        axes = [0, 1]
        slices = [pd.IndexSlice[:], pd.IndexSlice[:, ['A']],
                  pd.IndexSlice[[1], :], pd.IndexSlice[[1], ['A']],
                  pd.IndexSlice[:2, ['A', 'B']]]
        for ax in axes:
            for slice_ in slices:
                result = self.df.style.apply(self.h, axis=ax, subset=slice_,
                                             foo='baz')._compute().ctx
                expected = dict(((r, c), ['color: baz'])
                                for r, row in enumerate(self.df.index)
                                for c, col in enumerate(self.df.columns)
                                if row in self.df.loc[slice_].index
                                and col in self.df.loc[slice_].columns)
                self.assertEqual(result, expected)

    def test_applymap_subset(self):
        def f(x):
            return 'foo: bar'

        slices = [pd.IndexSlice[:], pd.IndexSlice[:, ['A']],
                  pd.IndexSlice[[1], :], pd.IndexSlice[[1], ['A']],
                  pd.IndexSlice[:2, ['A', 'B']]]

        for slice_ in slices:
            result = self.df.style.applymap(f, subset=slice_)._compute().ctx
            expected = dict(((r, c), ['foo: bar'])
                            for r, row in enumerate(self.df.index)
                            for c, col in enumerate(self.df.columns)
                            if row in self.df.loc[slice_].index
                            and col in self.df.loc[slice_].columns)
            self.assertEqual(result, expected)

    def test_empty(self):
        df = pd.DataFrame({'A': [1, 0]})
        s = df.style
        s.ctx = {(0, 0): ['color: red'],
                 (1, 0): ['']}

        result = s._translate()['cellstyle']
        expected = [{'props': [['color', ' red']], 'selector': 'row0_col0'},
                    {'props': [['', '']], 'selector': 'row1_col0'}]
        self.assertEqual(result, expected)

    def test_bar(self):
        df = pd.DataFrame({'A': [0, 1, 2]})
        result = df.style.bar()._compute().ctx
        expected = {
            (0, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,#d65f5f 0.0%, transparent 0%)'],
            (1, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,#d65f5f 50.0%, transparent 0%)'],
            (2, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,#d65f5f 100.0%, transparent 0%)']
        }
        self.assertEqual(result, expected)

        result = df.style.bar(color='red', width=50)._compute().ctx
        expected = {
            (0, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,red 0.0%, transparent 0%)'],
            (1, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,red 25.0%, transparent 0%)'],
            (2, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,red 50.0%, transparent 0%)']
        }
        self.assertEqual(result, expected)

        df['C'] = ['a'] * len(df)
        result = df.style.bar(color='red', width=50)._compute().ctx
        self.assertEqual(result, expected)
        df['C'] = df['C'].astype('category')
        result = df.style.bar(color='red', width=50)._compute().ctx
        self.assertEqual(result, expected)

    def test_highlight_null(self, null_color='red'):
        df = pd.DataFrame({'A': [0, np.nan]})
        result = df.style.highlight_null()._compute().ctx
        expected = {(0, 0): [''],
                    (1, 0): ['background-color: red']}
        self.assertEqual(result, expected)

    def test_nonunique_raises(self):
        df = pd.DataFrame([[1, 2]], columns=['A', 'A'])
        with tm.assertRaises(ValueError):
            df.style

        with tm.assertRaises(ValueError):
            Styler(df)

    def test_caption(self):
        styler = Styler(self.df, caption='foo')
        result = styler.render()
        self.assertTrue(all(['caption' in result, 'foo' in result]))

        styler = self.df.style
        result = styler.set_caption('baz')
        self.assertTrue(styler is result)
        self.assertEqual(styler.caption, 'baz')

    def test_uuid(self):
        styler = Styler(self.df, uuid='abc123')
        result = styler.render()
        self.assertTrue('abc123' in result)

        styler = self.df.style
        result = styler.set_uuid('aaa')
        self.assertTrue(result is styler)
        self.assertEqual(result.uuid, 'aaa')

    def test_table_styles(self):
        style = [{'selector': 'th', 'props': [('foo', 'bar')]}]
        styler = Styler(self.df, table_styles=style)
        result = ' '.join(styler.render().split())
        self.assertTrue('th { foo: bar; }' in result)

        styler = self.df.style
        result = styler.set_table_styles(style)
        self.assertTrue(styler is result)
        self.assertEqual(styler.table_styles, style)

    def test_table_attributes(self):
        attributes = 'class="foo" data-bar'
        styler = Styler(self.df, table_attributes=attributes)
        result = styler.render()
        self.assertTrue('class="foo" data-bar' in result)

        result = self.df.style.set_table_attributes(attributes).render()
        self.assertTrue('class="foo" data-bar' in result)

    def test_precision(self):
        with pd.option_context('display.precision', 10):
            s = Styler(self.df)
        self.assertEqual(s.precision, 10)
        s = Styler(self.df, precision=2)
        self.assertEqual(s.precision, 2)

        s2 = s.set_precision(4)
        self.assertTrue(s is s2)
        self.assertEqual(s.precision, 4)

    def test_apply_none(self):
        def f(x):
            return pd.DataFrame(np.where(x == x.max(), 'color: red', ''),
                                index=x.index, columns=x.columns)
        result = (pd.DataFrame([[1, 2], [3, 4]])
                  .style.apply(f, axis=None)._compute().ctx)
        self.assertEqual(result[(1, 1)], ['color: red'])

    def test_trim(self):
        result = self.df.style.render()  # trim=True
        self.assertEqual(result.count('#'), 0)

        result = self.df.style.highlight_max().render()
        self.assertEqual(result.count('#'), len(self.df.columns))

    def test_highlight_max(self):
        df = pd.DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        # max(df) = min(-df)
        for max_ in [True, False]:
            if max_:
                attr = 'highlight_max'
            else:
                df = -df
                attr = 'highlight_min'
            result = getattr(df.style, attr)()._compute().ctx
            self.assertEqual(result[(1, 1)], ['background-color: yellow'])

            result = getattr(df.style, attr)(color='green')._compute().ctx
            self.assertEqual(result[(1, 1)], ['background-color: green'])

            result = getattr(df.style, attr)(subset='A')._compute().ctx
            self.assertEqual(result[(1, 0)], ['background-color: yellow'])

            result = getattr(df.style, attr)(axis=0)._compute().ctx
            expected = {(1, 0): ['background-color: yellow'],
                        (1, 1): ['background-color: yellow'],
                        (0, 1): [''], (0, 0): ['']}
            self.assertEqual(result, expected)

            result = getattr(df.style, attr)(axis=1)._compute().ctx
            expected = {(0, 1): ['background-color: yellow'],
                        (1, 1): ['background-color: yellow'],
                        (0, 0): [''], (1, 0): ['']}
            self.assertEqual(result, expected)

        # separate since we cant negate the strs
        df['C'] = ['a', 'b']
        result = df.style.highlight_max()._compute().ctx
        expected = {(1, 1): ['background-color: yellow']}

        result = df.style.highlight_min()._compute().ctx
        expected = {(0, 0): ['background-color: yellow']}

    def test_export(self):
        f = lambda x: 'color: red' if x > 0 else 'color: blue'
        g = lambda x, y, z: 'color: %s' if x > 0 else 'color: %s' % z
        style1 = self.styler
        style1.applymap(f)\
            .applymap(g, y='a', z='b')\
            .highlight_max()
        result = style1.export()
        style2 = self.df.style
        style2.use(result)
        self.assertEqual(style1._todo, style2._todo)
        style2.render()


@tm.mplskip
class TestStylerMatplotlibDep(TestCase):

    def test_background_gradient(self):
        df = pd.DataFrame([[1, 2], [2, 4]], columns=['A', 'B'])
        for axis in [0, 1, 'index', 'columns']:
            for cmap in [None, 'YlOrRd']:
                result = df.style.background_gradient(cmap=cmap)._compute().ctx
                self.assertTrue(all("#" in x[0] for x in result.values()))
                self.assertEqual(result[(0, 0)], result[(0, 1)])
                self.assertEqual(result[(1, 0)], result[(1, 1)])

        result = (df.style.background_gradient(subset=pd.IndexSlice[1, 'A'])
                    ._compute().ctx)
        self.assertEqual(result[(1, 0)], ['background-color: #fff7fb'])
