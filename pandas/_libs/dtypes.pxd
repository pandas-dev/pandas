"""
Common location for shared fused types
"""

from numpy cimport (
    float32_t,
    float64_t,
    int64_t,
    uint64_t,
)

ctypedef fused rank_t:
    float64_t
    float32_t
    int64_t
    uint64_t
    object
