import numpy as np
import pandas as pd
import pytest

def calculate_expected_output(df, window_size, min_periods):
    results = []
    for start in range(len(df) - window_size + 1):
        end = start + window_size
        window_df = df.iloc[start:end]
        if window_df[['value', 'weight']].notna().sum().min() >= min_periods:
            result = window_df['value'] * window_df['weight']
            results.append(result.iloc[-1])
        else:
            results.append(np.nan)

    results = [np.nan] * (window_size - 1) + results
    return pd.Series(results, index=df.index)

@pytest.mark.parametrize('window_size', [3, 5])
def test_rolling_apply_multi_column_simple(window_size):
    # Create a sample DataFrame
    data = {'value': [1, 2, 3, 4, 5],
            'weight': [0.1, 0.2, 0.3, 0.4, 0.5]}
    df = pd.DataFrame(data)

    # Define a simple function to multiply 'value' and 'weight' columns
    def multiply_columns(window_df):
        return window_df['value'] * window_df['weight']

    # Apply the function using multi-column rolling apply
    result = df.rolling(window=window_size).apply(multiply_columns, multi_column=True)

    # Calculate the expected output using the separate function
    expected = calculate_expected_output(df, window_size, min_periods=1)

    # Assert the result matches the expected output
    pd.testing.assert_series_equal(result, expected)

@pytest.mark.parametrize('min_periods', [1, 2])
def test_rolling_apply_multi_column_missing_values(min_periods):
    # Create a sample DataFrame with missing values
    data = {'value': [1, np.nan, 3, np.nan, 5],
            'weight': [0.1, 0.2, np.nan, 0.4, 0.5]}
    df = pd.DataFrame(data)

    # Define a function to multiply 'value' and 'weight' columns
    def multiply_columns(window_df):
        return window_df['value'] * window_df['weight']

    # Apply the function using multi-column rolling apply
    result = df.rolling(window=3, min_periods=min_periods).apply(multiply_columns, multi_column=True)

    # Calculate the expected output using the separate function
    expected = calculate_expected_output(df, window_size=3, min_periods=min_periods)

    # Assert the result matches the expected output
    pd.testing.assert_series_equal(result, expected)

def test_rolling_apply_multi_column_large_dataframe():
    # Create a larger sample DataFrame
    data = {'value': np.arange(1, 101),
            'weight': np.arange(0.1, 10.1, 0.1)}
    df = pd.DataFrame(data)

    # Define a function to multiply 'value' and 'weight' columns
    def multiply_columns(window_df):
        return window_df['value'] * window_df['weight']

    # Apply the function using multi-column rolling apply
    result = df.rolling(window=5).apply(multiply_columns, multi_column=True)

    # Calculate the expected output using the separate function
    expected = calculate_expected_output(df, window_size=5, min_periods=1)

    # Assert the result matches the expected output
    pd.testing.assert_series_equal(result, expected)