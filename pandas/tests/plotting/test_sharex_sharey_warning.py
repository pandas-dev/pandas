import pytest
import pandas as pd

def test_plot_sharex_sharey_warning():
    # Create a simple DataFrame
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

    # Check for UserWarning when sharex or sharey is set with subplots=False
    with pytest.warns(UserWarning, match="parameters are ignored when 'subplots' is set to False") as record:
        df.plot(sharex=True, sharey=True, subplots=False)
    print(record[0].message)
