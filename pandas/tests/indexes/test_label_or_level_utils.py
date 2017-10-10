import pytest
import pandas as pd
from pandas import compat
import pandas.util.testing as tm
from pandas.core.dtypes.missing import array_equivalent


class TestLabelOrLevelUtils(object):
    """
    Test NDFrame utility methods used by operations that allow users to
    specify a mixture of levels and labels
    """

    # Setup
    # =====
    def setup_method(self):
        self.df1 = pd.DataFrame({'L1': [1, 2, 3],
                                 'L2': [11, 12, 13],
                                 'L3': ['A', 'B', 'C']})

    # Data preparation helpers
    # ========================
    @staticmethod
    def prepare_df(df, levels=None, axis=0):
        if levels:
            if isinstance(levels, compat.string_types):
                levels = [levels]
            df = df.set_index(levels)

        if axis == 1:
            # Transpose so index levels become column levels
            df = df.T
        return df

    def prepare_df1(self, levels=None, axis=0):
        """Return DataFrame with specified levels (list of any of 'L1', 'L2',
        and 'L3'). Remaining keys are left as labels"""
        return self.prepare_df(self.df1, levels=levels, axis=axis)

    def prepare_df_ambig(self, axis=0):
        """Return DataFrame with levels 'L1' and 'L2' and
        labels 'L1' and 'L3' """
        df = self.df1.set_index(['L1', 'L2'])
        df['L1'] = df['L3']

        if axis == 1:
            df = df.T

        return df

    def prepare_df_duplabels(self, axis=0):
        """Return DataFrame with level 'L1' and labels 'L2', 'L3', and 'L2' """
        df = self.df1.set_index(['L1'])
        df = pd.concat([df, df['L2']], axis=1)

        if axis == 1:
            df = df.T

        return df

    # Test is label/level reference
    # =============================
    @staticmethod
    def check_level_reference(frame, levels, axis):
        for level in levels:
            assert frame._is_level_reference(level, axis=axis)
            assert not frame._is_label_reference(level, axis=axis)
            assert frame._is_label_or_level_reference(level, axis=axis)

    @staticmethod
    def check_label_reference(frame, labels, axis):
        for label in labels:
            assert frame._is_label_reference(label, axis=axis)
            assert not frame._is_level_reference(label, axis=axis)
            assert frame._is_label_or_level_reference(label, axis=axis)

    # DataFrame
    # ---------
    @pytest.mark.parametrize('axis', [0, 1])
    def test_is_level_or_label_reference_df_simple(self, axis):

        # df has no named levels on axis
        df = self.prepare_df1(axis=axis)
        self.check_label_reference(df, ['L1', 'L2', 'L3'], axis=axis)

        # Set L1 as level on axis
        df = self.prepare_df1('L1', axis=axis)
        self.check_level_reference(df, ['L1'], axis=axis)
        self.check_label_reference(df, ['L2', 'L3'], axis=axis)

        # Set L1 and L2 as levels on axis
        df = self.prepare_df1(['L1', 'L2'], axis=axis)
        self.check_level_reference(df, ['L1', 'L2'], axis=axis)
        self.check_label_reference(df, ['L3'], axis=axis)

        # Set L1, L2, and L3 as levels on axis
        df = self.prepare_df1(['L1', 'L2', 'L3'], axis=axis)
        self.check_level_reference(df, ['L1', 'L2', 'L3'], axis=axis)

    @pytest.mark.parametrize('axis', [0, 1])
    def test_is_level_reference_df_ambig(self, axis):

        df = self.prepare_df_ambig(axis=axis)

        # df has both an on-axis level and off-axis label named L1
        # Therefore L1 should reference the label, not the level
        self.check_label_reference(df, ['L1'], axis=axis)

        # df has an on-axis level named L2 and it is not ambiguous
        # Therefore L2 is an level reference
        self.check_level_reference(df, ['L2'], axis=axis)

        # df has a column named L3 and it not an level reference
        self.check_label_reference(df, ['L3'], axis=axis)

    # Series
    # ------
    def test_is_level_reference_series_simple_axis0(self):

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2
        self.check_level_reference(s, ['L1'], axis=0)
        assert not s._is_level_reference('L2')

        # Make series with L1 and L2 as index
        s = self.df1.set_index(['L1', 'L2']).L3
        self.check_level_reference(s, ['L1', 'L2'], axis=0)
        assert not s._is_level_reference('L3')

    def test_is_level_reference_series_axis1_error(self):

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2

        with tm.assert_raises_regex(ValueError, "No axis named 1"):
            s._is_level_reference('L1', axis=1)

    # Test _check_label_or_level_ambiguity_df
    # =======================================

    # DataFrame
    # ---------
    @pytest.mark.parametrize('axis', [0, 1])
    def test_check_label_or_level_ambiguity_df(self, axis):

        df = self.prepare_df_ambig(axis=axis)

        # df has both an on-axis level and off-axis label named L1
        # Therefore L1 is ambiguous
        with tm.assert_produces_warning(FutureWarning,
                                        clear=True,
                                        check_stacklevel=False) as w:

            assert df._check_label_or_level_ambiguity('L1', axis=axis)
            warning_msg = w[0].message.args[0]
            if axis == 0:
                assert warning_msg.startswith("'L1' is both an index level "
                                              "and a column label")
            else:
                assert warning_msg.startswith("'L1' is both a column level "
                                              "and an index label")

        # df has an on-axis level named L2 and it is not ambiguous
        # No warning should be raised
        with tm.assert_produces_warning(None):
            assert not df._check_label_or_level_ambiguity('L2', axis=axis)

        # df has an off-axis label named L3 and it is not ambiguous
        with tm.assert_produces_warning(None):
            assert not df._is_level_reference('L3', axis=axis)

    # Series
    # ------
    @pytest.mark.parametrize('axis', [0, 1])
    def test_check_label_or_level_ambiguity_series(self, axis):

        # A series has only one axis and references are never ambiguous,
        # regardless of what axis is considered

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2
        with tm.assert_produces_warning(None):
            assert not s._check_label_or_level_ambiguity('L1', axis=axis)
            assert not s._check_label_or_level_ambiguity('L2', axis=axis)

        # Make series with L1 and L2 as index
        s = self.df1.set_index(['L1', 'L2']).L3
        with tm.assert_produces_warning(None):
            assert not s._check_label_or_level_ambiguity('L1', axis=axis)
            assert not s._check_label_or_level_ambiguity('L2', axis=axis)
            assert not s._check_label_or_level_ambiguity('L3', axis=axis)

    # Test _get_label_or_level_values
    # ===============================

    # DataFrame
    # ---------
    @staticmethod
    def check_labels(frame, labels, axis):
        for label in labels:
            if axis == 0:
                expected = frame[label]._values
            else:
                expected = frame.loc[label]._values

            result = frame._get_label_or_level_values(label, axis=axis)
            assert array_equivalent(expected, result)

    @staticmethod
    def check_levels(frame, levels, axis):
        for level in levels:
            if axis == 0:
                expected = frame.index.get_level_values(level=level)._values
            else:
                expected = (frame.columns
                            .get_level_values(level=level)
                            ._values)

            result = frame._get_label_or_level_values(level, axis=axis)
            assert array_equivalent(expected, result)

    @pytest.mark.parametrize('axis', [0, 1])
    def test_get_label_or_level_values_df_simple(self, axis):

        # ### df has no named index levels ###
        df = self.prepare_df1(axis=axis)
        self.check_labels(df, ['L1', 'L2', 'L3'], axis=axis)

        # ### Set L1 as index level ###
        df = self.prepare_df1('L1', axis=axis)
        self.check_labels(df, ['L2', 'L3'], axis=axis)
        self.check_levels(df, ['L1'], axis=axis)

        # ### Set L1 and L2 as index levels ###
        df = self.prepare_df1(['L1', 'L2'], axis=axis)
        self.check_labels(df, ['L3'], axis=axis)
        self.check_levels(df, ['L1', 'L2'], axis=axis)

        # ### Set L1, L2, and L3 as index levels ###
        df = self.prepare_df1(['L1', 'L2', 'L3'], axis=axis)
        self.check_levels(df, ['L1', 'L2', 'L3'], axis=axis)

    @pytest.mark.parametrize('axis', [0, 1])
    def test_get_label_or_level_values_df_ambig(self, axis):
        df = self.prepare_df_ambig(axis=axis)

        # df has both an on-axis level and off-axis label named L1
        # Therefore L1 is ambiguous but will default to label
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            self.check_labels(df, ['L1'], axis=axis)

        # df has an on-axis level named L2 and it is not ambiguous
        with tm.assert_produces_warning(None):
            self.check_levels(df, ['L2'], axis=axis)

        # df has an off-axis label named L3 and it is not ambiguous
        with tm.assert_produces_warning(None):
            self.check_labels(df, ['L3'], axis=axis)

    @pytest.mark.parametrize('axis', [0, 1])
    def test_get_label_or_level_values_df_duplabels(self, axis):

        df = self.prepare_df_duplabels(axis=axis)

        # df has unambiguous level 'L1'
        self.check_levels(df, ['L1'], axis=axis)

        # df has unique label 'L3'
        self.check_labels(df, ['L3'], axis=axis)

        # df has duplicate labels 'L2'
        if axis == 0:
            expected_msg = "The column label 'L2' is not unique"
        else:
            expected_msg = "The index label 'L2' is not unique"

        with tm.assert_raises_regex(ValueError, expected_msg):
            self.check_labels(df, ['L2'], axis=axis)

    # Series
    # ------
    def test_get_label_or_level_values_series_axis0(self):

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2
        self.check_levels(s, ['L1'], axis=0)

        # Make series with L1 and L2 as index
        s = self.df1.set_index(['L1', 'L2']).L3
        self.check_levels(s, ['L1', 'L2'], axis=0)

    def test_get_label_or_level_values_series_axis1_error(self):

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2

        with tm.assert_raises_regex(ValueError, "No axis named 1"):
            s._get_label_or_level_values('L1', axis=1)

    # Test _drop_labels_or_levels
    # ===========================
    @staticmethod
    def check_labels_dropped(frame, labels, axis):
        for label in labels:
            df_dropped = frame._drop_labels_or_levels(label, axis=axis)

            if axis == 0:
                assert label in frame.columns
                assert label not in df_dropped.columns
            else:
                assert label in frame.index
                assert label not in df_dropped.index

    @staticmethod
    def check_levels_dropped(frame, levels, axis):
        for level in levels:
            df_dropped = frame._drop_labels_or_levels(level, axis=axis)

            if axis == 0:
                assert level in frame.index.names
                assert level not in df_dropped.index.names
            else:
                assert level in frame.columns.names
                assert level not in df_dropped.columns.names

    # DataFrame
    # ---------
    @pytest.mark.parametrize('axis', [0, 1])
    def test_drop_labels_or_levels_df(self, axis):

        # ### df has no named index levels ###
        df = self.prepare_df1(axis=axis)
        self.check_labels_dropped(df, ['L1', 'L2', 'L3'], axis=axis)

        with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
            df._drop_labels_or_levels('L4', axis=axis)

        # ### Set L1 as index level ###
        df = self.prepare_df1('L1', axis=axis)
        self.check_labels_dropped(df, ['L2', 'L3'], axis=axis)
        self.check_levels_dropped(df, ['L1'], axis=axis)

        with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
            df._drop_labels_or_levels('L4', axis=axis)

        # ### Set L1 and L2 as index levels ###
        df = self.prepare_df1(['L1', 'L2'], axis=axis)
        self.check_labels_dropped(df, ['L3'], axis=axis)
        self.check_levels_dropped(df, ['L1', 'L2'], axis=axis)

        with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
            df._drop_labels_or_levels('L4', axis=axis)

        # ### Set L1, L2, and L3 as index levels ###
        df = self.prepare_df1(['L1', 'L2', 'L3'], axis=axis)
        self.check_levels_dropped(df, ['L1', 'L2', 'L3'], axis=axis)

        with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
            df._drop_labels_or_levels('L4', axis=axis)

    # Series
    # ------
    def test_drop_labels_or_levels_series(self):

        # Make series with L1 as index
        s = self.df1.set_index('L1').L2
        self.check_levels_dropped(s, ['L1'], axis=0)

        with tm.assert_raises_regex(ValueError, "not valid index levels"):
            s._drop_labels_or_levels('L4', axis=0)

        # Make series with L1 and L2 as index
        s = self.df1.set_index(['L1', 'L2']).L3
        self.check_levels_dropped(s, ['L1', 'L2'], axis=0)

        with tm.assert_raises_regex(ValueError, "not valid index levels"):
            s._drop_labels_or_levels('L4', axis=0)
