"""
Unit Testing II - Mocking & Stubbing: Database I/O Operations
Nithikesh Reddy
"""
import pytest
import pandas as pd
import numpy as np


class TestDatabaseIOMocking:
    """Test database I/O operations using mocks (FR-5)"""
    
    def test_read_sql_basic(self, monkeypatch):
        """
        Test basic SQL read operation with mocked database connection
        
        Test Oracle (FR-5): Reading a SQL query that returns 100 rows and 3 columns
        should create a DataFrame with 100 rows and 3 columns
    
        """
        # Setup: Mock data that would come from database
        expected_data = pd.DataFrame({
            'id': range(100),
            'name': [f'user_{i}' for i in range(100)],
            'value': np.random.rand(100)
        })
        
        def mock_read_sql(query, con, **kwargs):
            return expected_data
        
        # Apply mock
        monkeypatch.setattr(pd, 'read_sql', mock_read_sql)
        
        # Execute: Read from "database"
        result = pd.read_sql("SELECT * FROM users", con=None)
        
        # Verify Test Oracle: Shape is (100, 3)
        assert result.shape == (100, 3), f"Expected (100, 3), got {result.shape}"
        assert list(result.columns) == ['id', 'name', 'value']
        assert len(result) == 100
    
    def test_read_sql_empty_result(self, monkeypatch):
        """
        Test SQL query returning empty result set

        """
        # Setup: Mock empty result
        empty_data = pd.DataFrame(columns=['id', 'name', 'value'])
        
        def mock_read_sql(query, con, **kwargs):
            return empty_data
        
        monkeypatch.setattr(pd, 'read_sql', mock_read_sql)
        
        # Execute
        result = pd.read_sql("SELECT * FROM empty_table", con=None)
        
        # Verify: Empty DataFrame with correct columns
        assert len(result) == 0
        assert list(result.columns) == ['id', 'name', 'value']
        assert isinstance(result, pd.DataFrame)
    
    def test_read_sql_with_parameters(self, monkeypatch):
        """
        Test parameterized SQL queries
        """
        # Setup: Mock filtered data
        filtered_data = pd.DataFrame({
            'id': [5],
            'name': ['user_5'],
            'value': [0.5]
        })
        
        def mock_read_sql(query, con, params=None, **kwargs):
            if params and params.get('user_id') == 5:
                return filtered_data
            return pd.DataFrame()
        
        monkeypatch.setattr(pd, 'read_sql', mock_read_sql)
        
        # Execute: Parameterized query
        result = pd.read_sql(
            "SELECT * FROM users WHERE id = :user_id",
            con=None,
            params={'user_id': 5}
        )
        
        # Verify: Filtered result
        assert len(result) == 1
        assert result['id'].iloc[0] == 5
        
    def test_read_sql_dtype_handling(self, monkeypatch):
        """
        Test SQL result data type conversion
        
        Test Oracle (FR-5): SQL INTEGER should convert to int64, VARCHAR to string,
        DECIMAL to float64 in the resulting DataFrame
        
        """
        # Setup: Mock with specific dtypes (using dict to avoid dtype conversion)
        typed_data = pd.DataFrame({
            'int_col': [1, 2, 3],
            'str_col': ['a', 'b', 'c'],
            'float_col': [1.1, 2.2, 3.3]
        })
        # Explicitly set dtypes to ensure consistency
        typed_data['int_col'] = typed_data['int_col'].astype('int64')
        typed_data['float_col'] = typed_data['float_col'].astype('float64')
        
        def mock_read_sql(query, con, **kwargs):
            return typed_data
        
        monkeypatch.setattr(pd, 'read_sql', mock_read_sql)
        
        # Execute
        result = pd.read_sql("SELECT * FROM typed_table", con=None)
        
        # Verify Test Oracle: Correct data types
        assert result['int_col'].dtype == np.int64
        # In pandas 3.0, strings may use string dtype instead of object
        assert result['str_col'].dtype in [object, 'string', pd.StringDtype()]
        assert result['float_col'].dtype == np.float64
        
    def test_read_sql_connection_error_handling(self, monkeypatch):
        """
        Test error handling when database connection fails
    
        """
        # Setup: Mock to raise connection error
        def mock_read_sql(query, con, **kwargs):
            raise ConnectionError("Unable to connect to database")
        
        monkeypatch.setattr(pd, 'read_sql', mock_read_sql)
        
        # Execute & Verify: Should raise ConnectionError
        with pytest.raises(ConnectionError, match="Unable to connect"):
            pd.read_sql("SELECT * FROM users", con=None)
