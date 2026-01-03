from __future__ import annotations
from pandas.core.dtypes.missing import isna
from pandas.core.dtypes.common import is_numeric_dtype
from typing import TYPE_CHECKING, overload
import numpy as np

if TYPE_CHECKING:
    from pandas import DataFrame, Series



# --- Array-level helpers ---
def detect_outliers_array(
        values: np.ndarray, 
        method: str = 'iqr',
        **kwargs
        ) -> np.ndarray:
   """
   Detect outliers in an array, using pandas' missing value system.
    
   This works with:
   - np.ndarray with NaN
   - np.ndarray with NaT (datetime)
   - ExtensionArrays (IntegerArray, etc.)
    
   Parameters
   ----------
   values : np.ndarray or ExtensionArray
       The underlying array from a Series
   method : str
        'iqr' or 'zscore'
   **kwargs
        threshold parameter
        
   Returns
   -------
   np.ndarray of bool
        Outlier mask
   """  
   
   # Convert ExtensionArray to Numpy
   if not isinstance(values, np.ndarray):
       values = np.asanyarray(values)
    
   # Option for cases when all are missing values
   result = np.zeros(len(values), dtype = bool)
   na_mask = isna(values)
   if na_mask.all():
       return result
   
   # Checking whether it's numeric or not
   if not is_numeric_dtype(values):
       return result
   
   # Extracting valid values
   valid_values = values[~na_mask]
   
   # Dispatch for IQR and Z-score methods   
   if method == 'iqr':
       threshold = kwargs.get('threshold', 1.5)
       q1, q3 = np.percentile(valid_values, [25, 75])
       iqr = q3 - q1
    
       # Edge cases: is the valid_values too small or are the values identical?
       if (len(valid_values) < 4) or  (iqr == 0):
            return result
    
       # Compute outliers
       lower_bound = q1 - threshold * iqr
       upper_bound = q3 + threshold * iqr
       valid_outliers = ((valid_values < lower_bound) | (valid_values > upper_bound))
   elif method == 'zscore':
       threshold = kwargs.get('threshold', 3)
       mean = valid_values.mean()
       std = valid_values.std(ddof = 1)
       if len(valid_values) < 3 or std == 0:
            return result
       z_scores = np.abs((valid_values - mean) / std)
       valid_outliers = z_scores > threshold       
   else:
       raise ValueError(f'Unknown method: {method}')
       
   # Map back so it includes outliers
   result[~na_mask] = valid_outliers    

   return result   


def _get_outlier_bounds(
        values: np.ndarray, 
        method: str = 'iqr',
        **kwargs
        ) -> tuple[float, float] | None:
    """
    Get the lower and upper bounds for outlier detection.
    
    Returns (lower_bound, upper_bound) or None if can't determine.
    """
    # Remove missing values
    na_mask = isna(values)
    if na_mask.all():
        return None
    
    valid_values = values[~na_mask]
    
    if method == 'iqr':
        threshold = kwargs.get('threshold', 1.5)
        q1, q3 = np.percentile(valid_values, [25, 75])
        iqr = q3 - q1
        
        if len(valid_values) < 4 or iqr == 0:
            return None        
                
        lower = q1 - threshold * iqr
        upper = q3 + threshold * iqr
        return lower, upper
    
    elif method == 'zscore':
        threshold = kwargs.get('threshold', 3.0)
        mean = valid_values.mean()
        std = valid_values.std(ddof = 1)
        
        if len(valid_values) < 3 or std == 0:
            return None
        
        lower = mean - threshold * std
        upper = mean + threshold * std
        return lower, upper
    
    else:
        raise ValueError(f'Unknown method: {method}')



# --- Series-level helpers ---
def _detect_outliers_series(
        series: Series,
        method: str = 'iqr',
        **kwargs
        ) -> Series:
    """
    Detect outliers in a pandas Series.
    
    This extracts the underlying array, detects outliers, and returns
    a Series with proper index.
    
    Parameters
    ----------
    series : pd.Series
        Input Series
    method : str
        Detection method
    **kwargs
        Method parameters
        
    Returns
    -------
    pd.Series of bool
        Outlier mask with same index as input
    """
    from pandas import Series
    # Extract the underlying array and masking it to be returned
    # as part of the outliers identification series
    values = series._values
    outlier_mask = detect_outliers_array(
        values,
        method = method,
        **kwargs
        )
    result = Series(
        outlier_mask,
        index = series.index,
        name = series.name,
        dtype = bool,
        copy = False
    )    
    return result


def _drop_outliers_series(
        series: Series, 
        method: str = 'iqr',
        **kwargs
        ) -> Series:
    """
    Drop outliers from a Series.
    
    Parameters
    ----------
    series : pd.Series
    method : str
    **kwargs
        
    Returns
    -------
    pd.Series
        Series with outliers removed
    """
    # Get outlier mask
    outlier_mask = _detect_outliers_series(
        series, 
        method = method,
        **kwargs
        )
    return series[~outlier_mask]


def _fill_outliers_series(
        series: Series,
        value = None,
        method = {'filling': None, 'detection': 'iqr'},
        **kwargs
        ) -> Series:
    """
    Fill outliers in a Series.
    
    Parameters
    ----------
    series : Series
        Input Series
    value : scalar, optional
        Value to fill outliers with
    method : dict 
        - filling: {'mean', 'median', 'clip', 'ffill', 'bfill'}, optional
        - detection: 'iqr' or 'zscore'
        Method to fill with and detect outlier

    **kwargs
        Parameters for outlier detection
        
    Returns
    -------
    Series
        Series with outliers filled
    """    
    if value is None and method['filling'] is None:
        raise ValueError('Must specify either "value" or "method["filling"]"')
    
    if value is not None and method['filling'] is not None:
        raise ValueError('Cannot specify both "value" and "method["filling"]"')
    
    # Detect outliers
    outlier_mask = _detect_outliers_series(series, method = method['detection'], **kwargs)
    
    # If no outliers, return copy
    if not outlier_mask.any():
        return series.copy()
    
    # Make a copy to modify
    result = series.copy() if method['filling'] not in ['mean', 'median', 'clip'] else \
        series.astype(float).copy()
    
    # Fill with constant value
    if value is not None:
        result[outlier_mask] = value
        return result
    
    # Fill with method
    if method['filling'] == 'mean':
        fill_value = series[~outlier_mask].mean()
        result[outlier_mask] = fill_value
    
    elif method['filling'] == 'median':
        fill_value = series[~outlier_mask].median()
        result[outlier_mask] = fill_value
    
    elif method['filling'] == 'clip':
        # Get bounds and clip
        values = series._values
        bounds = _get_outlier_bounds(values, method = method['detection'], **kwargs)
        if bounds is not None:
            lower, upper = bounds
            result = result.clip(lower = lower, upper = upper)
    
    elif method['filling'] == 'ffill':
        # Mark outliers as NaN, then forward fill
        result[outlier_mask] = np.nan
        result = result.ffill()
    
    elif method['filling'] == 'bfill':
        # Mark outliers as NaN, then backward fill
        result[outlier_mask] = np.nan
        result = result.bfill()
    
    else:
        raise ValueError(f'Invalid option for "filling_method": {method["filling"]}')
    
    return result


# --- DataFrame-level helpers ---
def _drop_outliers_dataframe(
        df: DataFrame, 
        method: str = 'iqr', 
        axis = 0,
        **kwargs
        ) -> DataFrame:
    """
    Drop outliers from a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
    method : str
    axis : {0, 1, None}
        0 or None: drop rows with outliers in any column
        1: drop columns with outliers in any row
    **kwargs
        
    Returns
    -------
    pd.DataFrame
        DataFrame with outliers removed
    """  
    from pandas import Series
    if axis is None:
        axis = 0    
        
    if axis == 0:        
        # Track which rows have outliers
        rows_with_outliers = Series(False, index = df.index)        
        for col in df.columns:
            if is_numeric_dtype(df[col]):
                col_outliers = _detect_outliers_series(
                    df[col], 
                    method = method,
                    **kwargs
                    )
                rows_with_outliers = rows_with_outliers | col_outliers
        return df[~rows_with_outliers]    
    
    elif axis == 1:
        # Drop columns that have outliers in any row        
        cols_to_keep = []        
        for col in df.columns:
            if is_numeric_dtype(df[col]):
                col_outliers = _detect_outliers_series(
                    df[col],
                    method = method,
                    **kwargs
                    )
                
                # Check if there's any outlier in the column
                if not col_outliers.any():
                    cols_to_keep.append(col)
            else:
                # Keep non-numeric columns
                cols_to_keep.append(col)
        
        return df[cols_to_keep]
    
    else:
        raise ValueError(f'axis must be 0 or 1, got {axis}')
        

def _fill_outliers_dataframe(
        df: DataFrame,
        value = None,
        method = {'filling': None, 'detection': 'iqr'},
        **kwargs
        ) -> DataFrame:
    """
    Fill outliers in a DataFrame.
    
    Parameters
    ----------
    df : DataFrame
        Input DataFrame
    value : scalar, dict, Series, or DataFrame, optional
        Value(s) to fill outliers with    
    method : dict 
        - filling: {'mean', 'median', 'clip', 'ffill', 'bfill'}, optional
        - detection: 'iqr' or 'zscore'
        Method to fill with and detect outlier
    **kwargs
        Parameters for outlier detection
        
    Returns
    -------
    DataFrame
        DataFrame with outliers filled
    """
    from pandas import Series, DataFrame
    
    result = df.copy()
    
    # Handle different value types
    if value is not None:
        if isinstance(value, dict):
            # Different value per column
            for col, col_value in value.items():
                if col in result.columns and is_numeric_dtype(result[col]):
                    result[col] = _fill_outliers_series(
                        result[col],
                        value = col_value,
                        method = method,
                        **kwargs
                    )
        
        elif isinstance(value, Series):
            # Use Series values for matching columns
            for col in result.columns:
                if col in value.index and is_numeric_dtype(result[col]):
                    result[col] = _fill_outliers_series(
                        result[col],
                        value = value[col],
                        method = method,
                        **kwargs
                    )
        
        elif isinstance(value, DataFrame):
            # Use corresponding DataFrame values
            for col in result.columns:
                if col in value.columns and is_numeric_dtype(result[col]):
                    outlier_mask = _detect_outliers_series(
                        result[col],
                        method = method['detection'],
                        **kwargs
                    )
                    result.loc[outlier_mask, col] = value.loc[outlier_mask, col]
        
        else:
            # Scalar - fill all numeric columns
            for col in result.columns:
                if is_numeric_dtype(result[col]):
                    result[col] = _fill_outliers_series(
                        result[col],
                        value = value,
                        method = method,
                        **kwargs
                    )
    
    else:
        # Method-based filling - apply to each numeric column
        for col in result.columns:
            if is_numeric_dtype(result[col]):
                result[col] = _fill_outliers_series(
                    result[col],
                    method = method,
                    **kwargs
                )
    
    return result


# --- Public API ---
@overload 
def drop_outliers(
        obj: Series, 
        method: str = 'iqr', 
        axis: int | None = None,
        **kwargs
        ) -> Series: 
    ... 

@overload 
def drop_outliers(
        obj: DataFrame, 
        method: str = 'iqr', 
        axis: int | None = 0, 
        **kwargs
        ) -> DataFrame: 
    ...
def drop_outliers(
        obj: Series | DataFrame, 
        method: str = 'iqr', 
        axis: int | None = None,
        **kwargs
        ) -> Series | DataFrame:
    """
    Remove outliers from Series or DataFrame.
    
    Parameters
    ----------
    obj : Series or DataFrame
        Data to remove outliers from
   method : {'iqr', 'zscore'}, default 'iqr'
        Outlier detection method
    axis : {0, 1, None}, default None
        For DataFrame:
        - 0 or None: drop rows with outliers
        - 1: drop columns with outliers
        For Series: ignored
    **kwargs
        threshold : float
            For 'iqr': IQR multiplier, default 1.5
            For 'zscore': z-score threshold, default 3.0
            
    Returns
    -------
    Same type as obj
        Object with outliers removed
        
    Examples
    --------
    >>> from pandas import DataFrame, Series  
    >>> from pandas.core.outliers import drop_outliers
    >>> s = Series([1, 2, 3, 4, 5, 100])
    >>> drop_outliers(s)
    0    1
    1    2
    2    3
    3    4
    4    5   
    dtype: int64
    
    >>> df = DataFrame({'A': [1, 2, 3, 4, 5, 100], 'B': [10, 20, 30, 40, 50, 60]})
    >>> drop_outliers(df)
         A   B
    0    1  10
    1    2  20
    2    3  30
    3    4  40
    4    5  50
    """
    from pandas import Series, DataFrame
    
    if isinstance(obj, Series):
        return _drop_outliers_series(obj, method = method, **kwargs)
    
    elif isinstance(obj, DataFrame):
        return _drop_outliers_dataframe(obj, method = method, axis = axis, **kwargs)
    
    else:
        raise TypeError(
            f'drop_outliers expects Series or DataFrame, got {type(obj).__name__}'
        )

@overload
def is_outlier(
        obj: Series | DataFrame, method: str = 'iqr', **kwargs
        ) -> Series[bool] | DataFrame[bool]:
    ...

def is_outlier(
        obj: Series | DataFrame, method: str = 'iqr', **kwargs
        ) -> Series[bool] | DataFrame[bool]:
    """
    Detect outliers element-wise.
    
    Parameters
    ----------
    obj : Series or DataFrame
        Data to check for outliers
    method : {'iqr', 'zscore'}, default 'iqr'
        Detection method
    **kwargs
        threshold parameter
        
    Returns
    -------
    Same type as obj, dtype bool
        Boolean mask indicating outlier positions
        
    Examples
    --------
    >>> from pandas import Series
    >>> s = Series([1, 2, 3, 4, 100])
    >>> is_outlier(s)
    0    False
    1    False
    2    False
    3    False
    4    True
    dtype: bool
    """
    from pandas import Series, DataFrame
    
    if isinstance(obj, Series):
        return _detect_outliers_series(
            obj, method = method, **kwargs
            )
    
    elif isinstance(obj, DataFrame):
        return obj.apply(
            lambda col: _detect_outliers_series(col, method = method, **kwargs)
            )
    else:
        raise TypeError(
            f'is_outlier expects Series or DataFrame, got {type(obj).__name__}'
        )


@overload
def has_outliers(obj: Series | DataFrame, method: str = 'iqr', **kwargs) -> bool:
    ...

def has_outliers(obj: Series | DataFrame, method: str = 'iqr', **kwargs) -> bool:
    """
    Check if object contains any outliers.
    
    Parameters
    ----------
    obj : Series or DataFrame
    method : str
    axis : {{0, 1, None}}
        For DataFrame:
        - None: check if ANY value is outlier
        - 0: check each column
        - 1: check each row
        
    Returns
    -------
    bool
        
    Examples
    --------
    >>> from pandas import Series, DataFrame
    >>> from pandas.core.outliers import has_outliers
    >>> import numpy as np
    >>> values = np.linspace(1, 5, num = 1000).tolist() + [100]
    >>> s = Series(values)
    >>> has_outliers(s)
    True
    
    >>> value_b = np.linspace(2, 3.5, num = 1000).tolist() + [55]
    >>> df = DataFrame({'A': values, 'B': value_b})
    >>> has_outliers(df)
    True
    """   
    from pandas import Series, DataFrame
    
    if isinstance(obj, Series):
        mask = is_outlier(obj, method = method, **kwargs)
        return bool(mask.any())
    
    elif isinstance(obj, DataFrame):
            for col in obj.columns:
                mask = is_outlier(obj[col], method = method, **kwargs)
                if mask.any():
                    return True
            return False
    else:
        raise TypeError(
            f'has_outliers expects Series or DataFrame, got {type(obj).__name__}'
        )

@overload
def fill_outliers(
    obj: Series,
    value = None,
    method: dict[str, str | None] = {'filling': None, 'detection': 'iqr'},
    **kwargs
) -> Series: ...

@overload
def fill_outliers(
    obj: DataFrame,
    value = None,
    method: dict[str | None, str] = {'filling': None, 'detection': 'iqr'},
    **kwargs
) -> DataFrame: ...

def fill_outliers(
        obj,
        value = None,
        method: dict[str | None, str] = {'filling': None, 'detection': 'iqr'},
        **kwargs
        ) -> Series | DataFrame:
    """
    Fill outliers with a specified value or method.
    
    Parameters
    ----------
    obj : Series or DataFrame
        Object to fill outliers in
    value : scalar, dict, Series, or DataFrame, optional
        Value to use to fill outliers. Ignored if method is specified.
        - scalar: Fill all outliers with this value
        - dict/Series: Fill each column with different values (DataFrame only)
        - DataFrame: Use corresponding values from this DataFrame
    method : dict
        Method to use for "filling": {'mean', 'median', 'clip', 'ffill', 'bfill'}, optional
        - 'mean': Replace with mean of non-outlier values
        - 'median': Replace with median of non-outlier values
        - 'clip': Clip to outlier detection boundaries
        - 'ffill': Forward fill from last non-outlier value
        - 'bfill': Backward fill from next non-outlier value
        If None, must provide value.
    Method to use for outlier "detection": {'iqr', 'zscore'}, default 'iqr'        
    **kwargs
        threshold : float
            Threshold for outlier detection
            
    Returns
    -------
    Same type as obj
        Object with outliers filled
        
    Examples
    --------
    Fill with a constant value:
    >>> from pandas import DataFrame, Series
    >>> s = Series([1, 2, 3, 4, 100])
    >>> fill_outliers(s, value = 5)
    0    1
    1    2
    2    3
    3    4
    4    5
    dtype: int64
    
    Fill with median of non-outlier values:
    
    >>> fill_outliers(s, method = {'filling': 'median', 'detection': 'iqr'})
    0    1.0
    1    2.0
    2    3.0
    3    4.0
    4    2.5
    dtype: float64
    
    Clip outliers to boundaries:
    
    >>> fill_outliers(s, method = {'filling': 'clip', 'detection': 'iqr'})
    0    1.0
    1    2.0
    2    3.0
    3    4.0
    4    7.0
    dtype: float64
    
    Fill DataFrame with different values per column:
    
    >>> df = DataFrame({'A': [1, 2, 3, 4, 5, 100], 'B': [10, 20, 30, 40, 50, 60]})
    >>> fill_outliers(df, value = {'A': 0, 'B': 999})
         A    B
    0    1   10
    1    2   20
    2    3   30
    3    4   40
    4    5   50
    5    0   60
    """
    from pandas import Series, DataFrame
    
    if isinstance(obj, Series):
        return _fill_outliers_series(
            obj,
            value = value,
            method = method,
            **kwargs
        )
    
    elif isinstance(obj, DataFrame):
        return _fill_outliers_dataframe(
            obj,
            value = value,
            method = method,
            **kwargs
        )
    
    else:
        raise TypeError(
            f'fill_outliers expects Series or DataFrame, got {type(obj).__name__}'
        )