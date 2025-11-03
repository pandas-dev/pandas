"""
Unit Testing II - Mocking & Stubbing: DateTime Operations
Malikarjuna
"""

import pytest
import pandas as pd
import numpy as np
from datetime import datetime, timedelta


class TestDateTimeOperationsMocking:
    """Test time-series operations using mocks (FR-6)"""
    
    def test_timestamp_now_mocked(self, monkeypatch):
        """
        Test current timestamp creation with controlled time
        
        
        """
        # Setup: Fix current time to specific moment
        fixed_time = pd.Timestamp('2024-01-15 12:00:00')
        
        def mock_now(tz=None):
            return fixed_time
        
        monkeypatch.setattr(pd.Timestamp, 'now', staticmethod(mock_now))
        
        # Execute
        result = pd.Timestamp.now()
        
        # Verify: Time is exactly as mocked
        assert result == fixed_time
        assert result.year == 2024
        assert result.month == 1
        assert result.day == 15
    
    def test_date_range_generation(self, monkeypatch):
        """
        Test date range generation for time-series
        
        Test Oracle : Creating a date range for 365 days at daily frequency
        should produce exactly 365 timestamps
        
        
        """
        # Setup: Mock date range
        expected_dates = pd.date_range('2023-01-01', periods=365, freq='D')
        
        original_date_range = pd.date_range
        
        def mock_date_range(start=None, end=None, periods=None, freq=None, **kwargs):
            if periods == 365 and freq == 'D':
                return expected_dates
            return original_date_range(start, end, periods, freq, **kwargs)
        
        monkeypatch.setattr(pd, 'date_range', mock_date_range)
        
        # Execute
        result = pd.date_range('2023-01-01', periods=365, freq='D')
        
        # Verify Test Oracle: Exactly 365 dates
        assert len(result) == 365
        assert result[0] == pd.Timestamp('2023-01-01')
        assert result[-1] == pd.Timestamp('2023-12-31')
    
    def test_time_series_resampling(self, monkeypatch):
        """
        Test time-series resampling operation (FR-6)
        
        
        """
        # Setup: Create time-series data
        dates = pd.date_range('2023-01-01', periods=100, freq='h')
        df = pd.DataFrame({
            'value': np.random.rand(100)
        }, index=dates)
        
        # Mock the resample method to return a controlled result
        original_resample = pd.DataFrame.resample
        
        def mock_resample(self, rule, **kwargs):
            if rule == 'D':
                # Return a mock resampler that returns daily means
                class MockResampler:
                    def mean(inner_self):
                        return pd.DataFrame({
                            'value': [0.5, 0.6, 0.4, 0.7]
                        }, index=pd.date_range('2023-01-01', periods=4, freq='D'))
                return MockResampler()
            return original_resample(self, rule, **kwargs)
        
        monkeypatch.setattr(pd.DataFrame, 'resample', mock_resample)
        
        # Execute
        result = df.resample('D').mean()
        
        # Verify
        assert len(result) == 4
        assert 'value' in result.columns
        
    def test_rolling_window_operations(self, monkeypatch):
        """
        Test rolling window calculations (FR-6)
        
        Test Oracle: Rolling mean with window=7 on 30-day data should
        produce 30 values with first 6 as NaN
        
       
        """
        # Setup: Time-series data
        dates = pd.date_range('2023-01-01', periods=30, freq='D')
        df = pd.DataFrame({
            'price': range(30)
        }, index=dates)
        
        # Mock rolling method
        original_rolling = pd.DataFrame.rolling
        
        def mock_rolling(self, window, **kwargs):
            if window == 7:
                class MockRoller:
                    def mean(inner_self):
                        expected_result = pd.Series(
                            [np.nan]*6 + list(range(3, 27)),
                            index=dates
                        )
                        return expected_result
                return MockRoller()
            return original_rolling(self, window, **kwargs)
        
        monkeypatch.setattr(pd.DataFrame, 'rolling', mock_rolling)
        
        # Execute
        result = df['price'].rolling(window=7).mean()
        
        # Verify Test Oracle
        assert len(result) == 30
        assert pd.isna(result.iloc[:6]).all()  # First 6 are NaN
        
    def test_datetime_parsing_with_format(self, monkeypatch):
        """
        Test datetime string parsing with custom format
        
    
        """
        # Setup: Mock parsing of custom date format
        date_strings = ['2023-01-15', '2023-02-20', '2023-03-25']
        expected_dates = pd.DatetimeIndex([
            pd.Timestamp('2023-01-15'),
            pd.Timestamp('2023-02-20'),
            pd.Timestamp('2023-03-25')
        ])
        
        original_to_datetime = pd.to_datetime
        
        def mock_to_datetime(arg, format=None, **kwargs):
            if format == '%Y-%m-%d' and arg == date_strings:
                return expected_dates
            return original_to_datetime(arg, format=format, **kwargs)
        
        monkeypatch.setattr(pd, 'to_datetime', mock_to_datetime)
        
        # Execute
        result = pd.to_datetime(date_strings, format='%Y-%m-%d')
        
        # Verify
        assert len(result) == 3
        assert result[0] == pd.Timestamp('2023-01-15')
