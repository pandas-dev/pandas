# Polars Engine Implementation for pandas read_csv

This document summarizes the implementation of the polars engine for pandas' `read_csv` function.

## Files Modified/Created

### 1. Core Implementation
- **`pandas/io/parsers/polars_parser_wrapper.py`** - New file implementing the PolarsParserWrapper class
- **`pandas/_typing.py`** - Updated CSVEngine type to include "polars"
- **`pandas/io/parsers/readers.py`** - Updated to include polars engine support

### 2. Compatibility Support
- **`pandas/compat/polars.py`** - New file for polars compatibility checks
- **`pandas/compat/__init__.py`** - Updated to export HAS_POLARS

### 3. Test Infrastructure  
- **`pandas/tests/io/parser/conftest.py`** - Updated to include PolarsParser class and test fixtures
- **`pandas/tests/io/parser/test_polars_engine.py`** - New test file for polars engine specific tests

## Key Features Implemented

### Basic Functionality
- ✅ Reading CSV files with polars engine
- ✅ Converting polars DataFrame to pandas DataFrame
- ✅ Support for file paths and file-like objects
- ✅ Lazy evaluation using polars scan_csv when possible

### Supported Options
- ✅ `sep` - Field delimiter
- ✅ `header` - Row number(s) to use as column names
- ✅ `skiprows` - Lines to skip at start of file
- ✅ `na_values` - Additional strings to recognize as NA/NaN
- ✅ `names` - List of column names to use
- ✅ `usecols` - Return subset of columns (string names only)
- ✅ `nrows` - Number of rows to read
- ✅ `quotechar` - Character used to quote fields
- ✅ `comment` - Character(s) to treat as comment
- ✅ `encoding` - Encoding to use for UTF when reading
- ✅ `dtype` - Data type for data or columns (dict mapping)

### Unsupported Options (raises ValueError)
- ❌ `chunksize` - Not supported (similar to pyarrow)
- ❌ `iterator` - Not supported (similar to pyarrow)
- ❌ `skipfooter` - Not supported
- ❌ `float_precision` - Not supported
- ❌ `thousands` - Not supported
- ❌ `memory_map` - Not supported
- ❌ `dialect` - Not supported
- ❌ `quoting` - Not supported
- ❌ `lineterminator` - Not supported
- ❌ `converters` - Not supported
- ❌ `dayfirst` - Not supported
- ❌ `skipinitialspace` - Not supported
- ❌ `low_memory` - Not supported
- ❌ Callable `usecols` - Not supported
- ❌ Dict `na_values` - Not supported

## Performance Benefits

The polars engine is designed to provide:

1. **Fast CSV parsing** - Polars has state-of-the-art CSV parsing performance
2. **Memory efficiency** - Lazy evaluation where possible
3. **Parallel processing** - Polars can utilize multiple CPU cores
4. **Column pruning** - Only read requested columns when using `usecols`
5. **Predicate pushdown** - Future optimization for row filtering

## Usage Examples

```python
import pandas as pd

# Basic usage
df = pd.read_csv("data.csv", engine="polars")

# With options
df = pd.read_csv("data.csv", 
                 engine="polars",
                 usecols=["name", "age"],
                 nrows=1000,
                 na_values=["NULL", "N/A"])

# Custom column names
df = pd.read_csv("data.csv", 
                 engine="polars",
                 names=["col1", "col2", "col3"],
                 header=None)
```

## Error Handling

The implementation includes comprehensive error handling:

1. **Missing polars dependency** - Graceful ImportError with suggestion to install polars
2. **Unsupported options** - Clear ValueError messages listing unsupported parameters
3. **Polars parsing errors** - Wrapped in pandas ParserError with context
4. **File handling errors** - Proper cleanup and error propagation

## Testing

A comprehensive test suite has been implemented covering:

- Basic functionality tests
- Option validation tests
- Error condition tests
- Comparison with other engines
- Edge cases and compatibility

## Future Enhancements

Potential improvements for future versions:

1. **Enhanced dtype mapping** - Better support for pandas-specific dtypes
2. **Date parsing** - Leverage polars' built-in date parsing capabilities
3. **Index handling** - More sophisticated index column processing
4. **Streaming support** - Large file processing with minimal memory usage
5. **Schema inference** - Automatic optimal dtype detection

## Documentation Updates

The implementation includes updated documentation:

- Engine parameter documentation in `read_csv` docstring
- Version notes indicating experimental status
- Clear listing of supported and unsupported options

## Implementation Notes

### Design Decisions

1. **Lazy evaluation preferred** - Uses `scan_csv` for file paths when possible
2. **Pandas compatibility first** - All results converted to pandas DataFrame
3. **Error parity** - Similar error handling to existing engines
4. **Test infrastructure reuse** - Leverages existing parser test framework

### Limitations

1. **Experimental status** - Marked as experimental similar to pyarrow engine
2. **Option subset** - Only supports subset of pandas read_csv options
3. **Polars dependency** - Requires polars to be installed
4. **Performance trade-off** - Conversion to pandas may negate some performance benefits

This implementation provides a solid foundation for using polars as a high-performance CSV parsing engine within pandas while maintaining compatibility with the existing pandas API.
