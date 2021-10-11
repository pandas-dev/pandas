import pandas as pd
import pytest

def test_basic():
    #gh43564
    df = pd.DataFrame({'gender': ['male', 'male', 'female', 'male', 'female', 'male'], 
                      'education': ['low', 'medium', 'high', 'low', 'high', 'low'],
                      'country': ['US', 'FR', 'US', 'FR', 'FR', 'US']})
    print("SeriesGroupBy:\n", df.groupby('country')['gender'].value_counts(normalize=True))


    def compute_proportions(df, var):
        return (df[var]
            .value_counts(normalize=True)
        )
    print("apply:\n", df.groupby('country')
      .apply(compute_proportions, var=['gender', 'education'])
    )

    # but fails on DataFrameGroupBy
    print("DataFrameGroupBy:\n", df.groupby('country')[['gender', 'education']].value_counts(normalize=True))
