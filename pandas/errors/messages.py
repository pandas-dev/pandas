"""
Error message templates shared across pandas internals. Not public API.
"""

# for a value setitem-like validation rejects; callers may append a dtype-specific hint
INVALID_SETITEM_VALUE = "Invalid value '{value!s}' for dtype '{dtype}'"
# for list-likes and NA-like scalars, whose text is long or misleading
INVALID_SETITEM_LISTLIKE = "Invalid value for dtype '{dtype}'"
