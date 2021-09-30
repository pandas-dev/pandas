"""
Common location for shared fused types
"""

from numpy cimport (
    float32_t,
    float64_t,
    int64_t,
    uint64_t,
)

ctypedef fused numeric_object_t:
    float64_t
    float32_t
    int64_t
    uint64_t
    object
