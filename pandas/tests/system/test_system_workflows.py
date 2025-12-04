"""
System-level black-box tests for pandas end-to-end workflows.

These tests validate complete user workflows through public APIs only,
treating the system as a black box without referencing internal implementation.

Test Categories:
1. Data Loading and Export Workflow (Sandeep Ramavath)
2. Data Cleaning and Transformation Workflow (Nithikesh Bobbili)
3. Aggregation and Analysis Workflow (Mallikarjuna)
"""
import os
import tempfile
import numpy as np
import pandas as pd
import pytest


class TestDataIOWorkflow:
    """
    System tests by Sandeep Ramavath.
    Validates end-to-end data import/export workflows through public API.
    """
    
    def test_csv_roundtrip_workflow(self, tmp_path):
        """
        Test Case: CSV Data Import-Export Workflow
        
        Pre-conditions:
        - Temporary directory available for file operations
        - pandas library installed and functional
        
        Test Steps:
        1. Create DataFrame with mixed data types using public API
        2. Export DataFrame to CSV file
        3. Import CSV file back into new DataFrame
        4. Verify data integrity and type preservation
        
        Expected Results:
        - CSV file created successfully
        - Data round-trips without loss
        - Numeric, string, and datetime types preserved
        - All values match original dataset
        """
        # Step 1: Create DataFrame with mixed types through public API
        original_data = pd.DataFrame({
            'id': [1, 2, 3, 4, 5],
            'name': ['Alice', 'Bob', 'Charlie', 'David', 'Eve'],
            'score': [95.5, 87.3, 92.1, 88.7, 91.4],
            'date': pd.date_range('2024-01-01', periods=5),
            'active': [True, False, True, True, False]
        })
        
        # Step 2: Export to CSV using public API
        csv_path = tmp_path / "test_data.csv"
        original_data.to_csv(csv_path, index=False)
        
        # Verify file exists
        assert csv_path.exists(), "CSV file should be created"
        
        # Step 3: Import from CSV using public API
        loaded_data = pd.read_csv(csv_path, parse_dates=['date'])
        
        # Step 4: Verify data integrity
        assert len(loaded_data) == 5, "Should load 5 rows"
        assert list(loaded_data.columns) == ['id', 'name', 'score', 'date', 'active']
        assert loaded_data['id'].tolist() == [1, 2, 3, 4, 5]
        assert loaded_data['name'].tolist() == ['Alice', 'Bob', 'Charlie', 'David', 'Eve']
        assert loaded_data['score'].tolist() == [95.5, 87.3, 92.1, 88.7, 91.4]
        assert loaded_data['active'].tolist() == [True, False, True, True, False]
        
        # Verify datetime parsing
        assert pd.api.types.is_datetime64_any_dtype(loaded_data['date'])


class TestDataCleaningWorkflow:
    """
    System tests by Nithikesh Bobbili.
    Validates end-to-end data cleaning and transformation workflows.
    """
    
    def test_missing_data_handling_workflow(self):
        """
        Test Case: Missing Data Cleaning Workflow
        
        Pre-conditions:
        - pandas library available
        - No external dependencies required
        
        Test Steps:
        1. Create DataFrame with missing values using public API
        2. Detect missing values using public methods
        3. Fill missing values using multiple strategies
        4. Verify all missing values handled correctly
        
        Expected Results:
        - Missing values correctly identified
        - Forward fill propagates last valid value
        - Backward fill propagates next valid value
        - Constant fill replaces with specified value
        - No missing values remain after filling
        """
        # Step 1: Create DataFrame with missing data
        data = pd.DataFrame({
            'A': [1, np.nan, 3, np.nan, 5],
            'B': [np.nan, 2, np.nan, 4, 5],
            'C': [1, 2, 3, 4, np.nan]
        })
        
        # Step 2: Detect missing values using public API
        missing_count = data.isnull().sum()
        assert missing_count['A'] == 2, "Column A should have 2 missing values"
        assert missing_count['B'] == 2, "Column B should have 2 missing values"
        assert missing_count['C'] == 1, "Column C should have 1 missing value"
        
        # Step 3a: Fill missing values with forward fill
        filled_ffill = data.ffill()
        assert filled_ffill.isnull().sum().sum() == 1, "Should have 1 remaining NaN at start"
        assert filled_ffill.loc[1, 'A'] == 1.0, "Should forward fill from previous value"
        
        # Step 3b: Fill missing values with backward fill
        filled_bfill = data.bfill()
        assert filled_bfill.isnull().sum().sum() == 1, "Should have 1 remaining NaN at end"
        assert filled_bfill.loc[0, 'B'] == 2.0, "Should backward fill from next value"
        
        # Step 3c: Fill with constant value
        filled_constant = data.fillna(0)
        assert filled_constant.isnull().sum().sum() == 0, "No missing values should remain"
        assert filled_constant.loc[1, 'A'] == 0.0, "Should fill with constant value"
        
        # Step 4: Verify complete workflow
        original_shape = data.shape
        assert filled_constant.shape == original_shape, "Shape should be preserved"
class TestAggregationWorkflow:
    """
    System tests by Mallikarjuna.
    Validates end-to-end data aggregation and analysis workflows.
    """
    
    def test_groupby_aggregation_workflow(self):
        """
        Test Case: Group-by Aggregation Analysis Workflow
        
        Pre-conditions:
        - pandas library functional
        - Sufficient memory for operations
        
        Test Steps:
        1. Create DataFrame with categorical and numeric data
        2. Group data by category using public API
        3. Apply multiple aggregation functions
        4. Verify aggregated results for each category
        5. Verify multiple aggregation functions work simultaneously
        
        Expected Results:
        - Data groups correctly by category
        - Mean aggregation produces correct averages
        - Sum aggregation produces correct totals
        - Count aggregation shows correct group sizes
        - Multiple aggregations work in single operation
        """
        # Step 1: Create DataFrame with categorical data
        data = pd.DataFrame({
            'category': ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B'],
            'value': [10, 20, 15, 25, 20, 30, 25, 35],
            'quantity': [1, 2, 3, 4, 5, 6, 7, 8]
        })
        
        # Step 2: Group by category using public API
        grouped = data.groupby('category')
        
        # Step 3a: Apply mean aggregation
        mean_result = grouped['value'].mean()
        assert mean_result['A'] == 17.5, "Category A mean should be 17.5"
        assert mean_result['B'] == 27.5, "Category B mean should be 27.5"
        
        # Step 3b: Apply sum aggregation
        sum_result = grouped['value'].sum()
        assert sum_result['A'] == 70, "Category A sum should be 70"
        assert sum_result['B'] == 110, "Category B sum should be 110"
        
        # Step 3c: Apply count aggregation
        count_result = grouped.size()
        assert count_result['A'] == 4, "Category A should have 4 items"
        assert count_result['B'] == 4, "Category B should have 4 items"
        
        # Step 4: Apply multiple aggregations simultaneously
        multi_agg = grouped['value'].agg(['mean', 'sum', 'count'])
        
        # Step 5: Verify multi-aggregation results
        assert multi_agg.loc['A', 'mean'] == 17.5
        assert multi_agg.loc['A', 'sum'] == 70
        assert multi_agg.loc['A', 'count'] == 4
        assert multi_agg.loc['B', 'mean'] == 27.5
        assert multi_agg.loc['B', 'sum'] == 110
        assert multi_agg.loc['B', 'count'] == 4
        
        # Verify shape of result
        assert multi_agg.shape == (2, 3), "Should have 2 categories and 3 aggregations"
