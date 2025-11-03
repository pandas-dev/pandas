"""
Mutation Testing: Real Integration Tests for Pandas DataFrame Operations
These tests actually exercise pandas code without mocking core functionality
"""
import pytest
import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path


class TestDataFrameOperations:
    """Tests that exercise real pandas code for mutation testing"""
    
    def test_dataframe_concat_basic(self):
        """
        Test DataFrame concatenation with real pandas code
        Tests pandas.core.reshape.concat.concat()
        """
        # Create test data
        df1 = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        df2 = pd.DataFrame({'A': [7, 8, 9], 'B': [10, 11, 12]})
        
        # Execute concatenation
        result = pd.concat([df1, df2])
        
        # Verify
        assert result.shape == (6, 2)
        assert len(result) == 6
        assert list(result.columns) == ['A', 'B']
        assert result['A'].tolist() == [1, 2, 3, 7, 8, 9]
        
    def test_dataframe_concat_axis1(self):
        """
        Test DataFrame concatenation along columns
        Tests axis parameter handling
        """
        df1 = pd.DataFrame({'A': [1, 2, 3]})
        df2 = pd.DataFrame({'B': [4, 5, 6]})
        
        result = pd.concat([df1, df2], axis=1)
        
        assert result.shape == (3, 2)
        assert list(result.columns) == ['A', 'B']
        assert result['A'].tolist() == [1, 2, 3]
        assert result['B'].tolist() == [4, 5, 6]
        
    def test_dataframe_merge_inner(self):
        """
        Test DataFrame merge operation (inner join)
        Tests pandas.core.reshape.merge.merge()
        """
        df1 = pd.DataFrame({'key': ['A', 'B', 'C'], 'value1': [1, 2, 3]})
        df2 = pd.DataFrame({'key': ['B', 'C', 'D'], 'value2': [4, 5, 6]})
        
        result = pd.merge(df1, df2, on='key', how='inner')
        
        assert result.shape == (2, 3)
        assert len(result) == 2
        assert list(result['key']) == ['B', 'C']
        assert list(result['value1']) == [2, 3]
        assert list(result['value2']) == [4, 5]
        
    def test_dataframe_merge_left(self):
        """
        Test DataFrame merge with left join
        Tests join type handling
        """
        df1 = pd.DataFrame({'key': ['A', 'B', 'C'], 'value1': [1, 2, 3]})
        df2 = pd.DataFrame({'key': ['B', 'C', 'D'], 'value2': [4, 5, 6]})
        
        result = pd.merge(df1, df2, on='key', how='left')
        
        assert result.shape == (3, 3)
        assert len(result) == 3
        assert list(result['key']) == ['A', 'B', 'C']
        assert list(result['value1']) == [1, 2, 3]
        # Value for 'A' should be NaN since it's not in df2
        assert pd.isna(result.loc[result['key'] == 'A', 'value2'].iloc[0])


class TestDataFrameGroupBy:
    """Tests for GroupBy operations - real pandas functionality"""
    
    def test_groupby_sum(self):
        """
        Test GroupBy sum aggregation
        Tests pandas.core.groupby.GroupBy.sum()
        """
        df = pd.DataFrame({
            'category': ['A', 'B', 'A', 'B', 'A'],
            'values': [10, 20, 30, 40, 50]
        })
        
        result = df.groupby('category')['values'].sum()
        
        assert result['A'] == 90  # 10 + 30 + 50
        assert result['B'] == 60  # 20 + 40
        assert len(result) == 2
        
    def test_groupby_mean(self):
        """
        Test GroupBy mean aggregation
        Tests aggregation function handling
        """
        df = pd.DataFrame({
            'category': ['A', 'B', 'A', 'B', 'A'],
            'values': [10, 20, 30, 40, 50]
        })
        
        result = df.groupby('category')['values'].mean()
        
        assert result['A'] == 30.0  # (10 + 30 + 50) / 3
        assert result['B'] == 30.0  # (20 + 40) / 2
        
    def test_groupby_count(self):
        """
        Test GroupBy count aggregation
        Tests counting functionality
        """
        df = pd.DataFrame({
            'category': ['A', 'B', 'A', 'B', 'A', 'C'],
            'values': [10, 20, 30, 40, 50, 60]
        })
        
        result = df.groupby('category')['values'].count()
        
        assert result['A'] == 3
        assert result['B'] == 2
        assert result['C'] == 1
        assert len(result) == 3


class TestDataFrameFiltering:
    """Tests for DataFrame filtering and selection"""
    
    def test_filter_by_condition(self):
        """
        Test filtering DataFrame by boolean condition
        Tests pandas indexing and boolean operations
        """
        df = pd.DataFrame({
            'A': [1, 2, 3, 4, 5],
            'B': [10, 20, 30, 40, 50]
        })
        
        result = df[df['A'] > 2]
        
        assert len(result) == 3
        assert list(result['A']) == [3, 4, 5]
        assert list(result['B']) == [30, 40, 50]
        
    def test_filter_multiple_conditions(self):
        """
        Test filtering with multiple conditions (AND)
        Tests compound boolean operations
        """
        df = pd.DataFrame({
            'A': [1, 2, 3, 4, 5],
            'B': [10, 20, 30, 40, 50]
        })
        
        result = df[(df['A'] > 2) & (df['B'] < 50)]
        
        assert len(result) == 2
        assert list(result['A']) == [3, 4]
        assert list(result['B']) == [30, 40]
        
    def test_isin_filter(self):
        """
        Test filtering using isin() method
        Tests membership testing
        """
        df = pd.DataFrame({
            'category': ['A', 'B', 'C', 'D', 'E'],
            'value': [1, 2, 3, 4, 5]
        })
        
        result = df[df['category'].isin(['A', 'C', 'E'])]
        
        assert len(result) == 3
        assert list(result['category']) == ['A', 'C', 'E']
        assert list(result['value']) == [1, 3, 5]


class TestDataFrameCSVIO:
    """Tests for real CSV I/O operations (not mocked)"""
    
    def test_read_write_csv_basic(self, tmp_path):
        """
        Test writing and reading CSV files
        Tests real pandas I/O functionality
        """
        # Create test DataFrame
        df = pd.DataFrame({
            'col1': [1, 2, 3],
            'col2': ['a', 'b', 'c'],
            'col3': [1.1, 2.2, 3.3]
        })
        
        # Write to CSV
        csv_file = tmp_path / "test.csv"
        df.to_csv(csv_file, index=False)
        
        # Read back
        result = pd.read_csv(csv_file)
        
        # Verify
        assert result.shape == (3, 3)
        assert list(result.columns) == ['col1', 'col2', 'col3']
        assert list(result['col1']) == [1, 2, 3]
        assert list(result['col2']) == ['a', 'b', 'c']
        
    def test_read_csv_with_delimiter(self, tmp_path):
        """
        Test reading CSV with custom delimiter
        Tests delimiter parameter handling
        """
        # Create TSV file
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("col1\tcol2\tcol3\n1\ta\t1.1\n2\tb\t2.2\n")
        
        # Read with tab delimiter
        result = pd.read_csv(tsv_file, sep='\t')
        
        assert result.shape == (2, 3)
        assert list(result.columns) == ['col1', 'col2', 'col3']
        assert list(result['col1']) == [1, 2]
        
    def test_read_csv_with_header(self, tmp_path):
        """
        Test reading CSV with specific header row
        Tests header parameter
        """
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("# Comment line\ncol1,col2\n1,2\n3,4\n")
        
        # Read with header on row 1 (0-indexed)
        result = pd.read_csv(csv_file, header=1)
        
        assert result.shape == (2, 2)
        assert list(result.columns) == ['col1', 'col2']


class TestDataFrameSorting:
    """Tests for DataFrame sorting operations"""
    
    def test_sort_values_ascending(self):
        """
        Test sorting DataFrame by values (ascending)
        Tests pandas.DataFrame.sort_values()
        """
        df = pd.DataFrame({
            'A': [3, 1, 2],
            'B': [6, 4, 5]
        })
        
        result = df.sort_values('A')
        
        assert list(result['A']) == [1, 2, 3]
        assert list(result['B']) == [4, 5, 6]
        
    def test_sort_values_descending(self):
        """
        Test sorting DataFrame by values (descending)
        Tests ascending parameter
        """
        df = pd.DataFrame({
            'A': [1, 3, 2],
            'B': [4, 6, 5]
        })
        
        result = df.sort_values('A', ascending=False)
        
        assert list(result['A']) == [3, 2, 1]
        assert list(result['B']) == [6, 5, 4]
        
    def test_sort_values_multiple_columns(self):
        """
        Test sorting by multiple columns
        Tests multi-column sorting
        """
        df = pd.DataFrame({
            'A': [1, 2, 1, 2],
            'B': [4, 3, 2, 1]
        })
        
        result = df.sort_values(['A', 'B'])
        
        assert list(result['A']) == [1, 1, 2, 2]
        assert list(result['B']) == [2, 4, 1, 3]
