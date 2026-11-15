"""
Shared error message templates, kept in one place so the wording stays consistent.
"""

# For a reduction, accumulation, or groupby operation a dtype cannot perform.
# Use as REDUCTION_MESSAGE.format(dtype=..., name=...).
REDUCTION_MESSAGE = "dtype '{dtype}' does not support operation '{name}'"
