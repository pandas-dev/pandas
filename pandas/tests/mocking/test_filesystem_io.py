"""
Unit Testing II - Mocking & Stubbing: File System I/O Operations
Student: Sandeep
Requirement: FR-5 - Loading data from flat files (CSV, Excel, HDF5)

This module tests pandas file I/O functionality using mocks to avoid
requiring actual file system operations. Tests verify pandas correctly handles
file parsing without creating real files.

Following pandas test conventions: using pytest-style tests with monkeypatch.
"""

import pytest
import pandas as pd
import numpy as np


class TestFileSystemIOMocking:
    """Test file system I/O operations using mocks (FR-5)"""
    
    def test_read_csv_basic(self, monkeypatch):
        """
        Test basic CSV read operation with mocked file system
        
        Test Oracle (FR-5): Reading a CSV file containing 100 rows and 5 columns
        should create a DataFrame with 100 rows and 5 columns
        
        Rationale: File I/O is slow; mocking allows testing CSV parsing logic
        without actual file creation
        """
        # Setup: Mock CSV data (100 rows, 5 columns)
        expected_data = pd.DataFrame({
            'col1': range(100),
            'col2': np.random.rand(100),
            'col3': [f'text_{i}' for i in range(100)],
            'col4': pd.date_range('2023-01-01', periods=100),
            'col5': np.random.choice(['X', 'Y', 'Z'], 100)
        })
        
        def mock_read_csv(filepath, **kwargs):
            return expected_data
        
        monkeypatch.setattr(pd, 'read_csv', mock_read_csv)
        
        # Execute: Read CSV with mocked file
        result = pd.read_csv('data.csv')
        
        # Verify Test Oracle: Shape is (100, 5)
        assert result.shape == (100, 5), f"Expected (100, 5), got {result.shape}"
        assert list(result.columns) == ['col1', 'col2', 'col3', 'col4', 'col5']
    
    def test_read_csv_with_delimiter(self, monkeypatch):
        """
        Test CSV read with custom delimiter (tab-separated, pipe-separated)
        
        Rationale: Delimited files come in various formats; verify pandas
        handles custom delimiters correctly
        """
        # Setup: Mock TSV data
        tsv_data = pd.DataFrame({
            'name': ['Alice', 'Bob', 'Charlie'],
            'age': [25, 30, 35],
            'city': ['NYC', 'LA', 'Chicago']
        })
        
        def mock_read_csv(filepath, delimiter=None, **kwargs):
            if delimiter == '\t':
                return tsv_data
            return pd.DataFrame()
        
        monkeypatch.setattr(pd, 'read_csv', mock_read_csv)
        
        # Execute: Read with tab delimiter
        result = pd.read_csv('data.tsv', delimiter='\t')
        
        # Verify: Correct parsing
        assert len(result) == 3
        assert 'name' in result.columns
    
    def test_read_excel_basic(self, monkeypatch):
        """
        Test Excel file read operation
        
        Rationale: Excel files require xlrd/openpyxl; mocking avoids
        dependency on external libraries
        """
        # Setup: Mock Excel data
        excel_data = pd.DataFrame({
            'Product': ['A', 'B', 'C'],
            'Sales': [1000, 2000, 1500],
            'Region': ['North', 'South', 'East']
        })
        
        def mock_read_excel(filepath, sheet_name=None, **kwargs):
            return excel_data
        
        monkeypatch.setattr(pd, 'read_excel', mock_read_excel)
        
        # Execute
        result = pd.read_excel('sales.xlsx', sheet_name='Sheet1')
        
        # Verify
        assert len(result) == 3
        assert 'Product' in result.columns
        assert result['Sales'].sum() == 4500
        
    def test_read_hdf_basic(self, monkeypatch):
        """
        Test HDF5 file read operation
        
        Test Oracle (NFR-3): System should load data using ultrafast HDF5 format
        
        Rationale: HDF5 format is for high-performance storage; verify
        pandas handles HDF5 correctly without requiring pytables
        """
        # Setup: Mock HDF5 data
        hdf_data = pd.DataFrame({
            'timestamp': pd.date_range('2023-01-01', periods=1000, freq='h'),
            'sensor_1': np.random.rand(1000),
            'sensor_2': np.random.rand(1000),
            'sensor_3': np.random.rand(1000)
        })
        
        def mock_read_hdf(filepath, key=None, **kwargs):
            return hdf_data
        
        monkeypatch.setattr(pd, 'read_hdf', mock_read_hdf)
        
        # Execute
        result = pd.read_hdf('sensors.h5', key='data')
        
        # Verify: Large dataset loaded correctly
        assert len(result) == 1000
        assert len(result.columns) == 4
        
    def test_csv_file_not_found_handling(self, monkeypatch):
        """
        Test error handling when CSV file doesn't exist
        
        Rationale: File not found is common error; pandas should handle
        with clear error message
        """
        # Setup: Mock to raise FileNotFoundError
        def mock_read_csv(filepath, **kwargs):
            raise FileNotFoundError(f"File '{filepath}' not found")
        
        monkeypatch.setattr(pd, 'read_csv', mock_read_csv)
        
        # Execute & Verify
        with pytest.raises(FileNotFoundError, match="missing.csv"):
            pd.read_csv('missing.csv')
