import pandas as pd
import pytest

@pytest.mark.xfail(strict=True)
def test():
    index = pd.period_range(start='2018-01', periods=24, freq='M')
    periodSerie = pd.Series(range(24),index=index)
    periodSerie.index.name = 'Month'
    periodSerie.groupby(periodSerie.index.month).sum()
