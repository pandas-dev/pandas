from pandas.core.api import Series

from pandas.stats.ols import OLS, MovingOLS
from pandas.stats.plm import PanelOLS, MovingPanelOLS, NonPooledPanelOLS
import pandas.stats.common as common

def ols(**kwargs):
    """Returns the appropriate OLS object depending on whether you need
    simple or panel OLS, and a full-sample or rolling/expanding OLS.

    Parameters
    ----------
    y: Series for simple OLS.  DataFrame for panel OLS.
    x: Series, DataFrame, or dict of Series for simple OLS.
       Dict of DataFrame for panel OLS.
    intercept: bool
        True if you want an intercept.  Defaults to True.
    nw_lags: None or int
        Number of Newey-West lags.  Defaults to None.
    nw_overlap: bool
        Whether there are overlaps in the NW lags.  Defaults to False.
    window_type: {'full sample', 'rolling', 'expanding'}
        'full sample' by default
    window: int
        size of window (for rolling/expanding OLS). If window passed and no
        explicit window_type, 'rolling" will be used as the window_type

    Panel OLS options:
        pool: bool
            Whether to run pooled panel regression.  Defaults to true.
        weights: DataFrame
            Weight for each observation.  The weights are not normalized;
            they're multiplied directly by each observation.
        entity_effects: bool
            Whether to account for entity fixed effects.  Defaults to false.
        time_effects: bool
            Whether to account for time fixed effects.  Defaults to false.
        x_effects: list
            List of x's to account for fixed effects.  Defaults to none.
        dropped_dummies: dict
            Key is the name of the variable for the fixed effect.
            Value is the value of that variable for which we drop the dummy.

            For entity fixed effects, key equals 'entity'.

            By default, the first dummy is dropped if no dummy is specified.
        cluster: {'time', 'entity'}
            cluster variances

    Returns
    -------
    The appropriate OLS object, which allows you to obtain betas and various
    statistics, such as std err, t-stat, etc.

    Examples
    --------
    # Run simple OLS.
    result = ols(y=y, x=x)

    # Run rolling simple OLS with window of size 10.
    result = ols(y=y, x=x, window_type='rolling', window=10)
    print result.beta

    result = ols(y=y, x=x, nw_lags=1)

    # Set up LHS and RHS for data across all items
    y = A
    x = {'B' : B, 'C' : C}

    # Run panel OLS.
    result = ols(y=y, x=x)

    # Run expanding panel OLS with window 10 and entity clustering.
    result = ols(y=y, x=x, cluster='entity', window_type='expanding', window=10)
    """
    pool = kwargs.get('pool')
    if 'pool' in kwargs:
        del kwargs['pool']

    window_type = kwargs.get('window_type')
    window = kwargs.get('window')

    if window_type is None:
        if window is None:
            window_type = common.FULL_SAMPLE
        else:
            window_type = common.ROLLING
    else:
        window_type = common._get_window_type(window_type)

    y = kwargs.get('y')
    if window_type == common.FULL_SAMPLE:
        for rolling_field in ('window_type', 'window', 'min_periods'):
            if rolling_field in kwargs:
                del kwargs[rolling_field]

        if isinstance(y, Series):
            klass = OLS
        else:
            if pool == False:
                klass = NonPooledPanelOLS
            else:
                klass = PanelOLS
    else:
        if isinstance(y, Series):
            klass = MovingOLS
        else:
            if pool == False:
                klass = NonPooledPanelOLS
            else:
                klass = MovingPanelOLS

    return klass(**kwargs)
