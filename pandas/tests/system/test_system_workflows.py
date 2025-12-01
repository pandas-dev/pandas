"""
System-level black-box tests for pandas end-to-end workflows.

These tests validate complete user workflows through public APIs only,
treating the system as a black box without referencing internal implementation.

Test Categories:
Data Loading and Export Workflow (Sandeep Ramavath)
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


