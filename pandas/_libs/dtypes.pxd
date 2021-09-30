"""
Common location for shared fused types
"""

from numpy cimport (
    float32_t,
    float64_t,
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)

ctypedef fused numeric_t:
    int8_t
    int16_t
    int32_t
    int64_t

    uint8_t
    uint16_t
    uint32_t
    uint64_t

    float32_t
    float64_t

ctypedef fused numeric_object_t:
    numeric_t
    object

ctypedef fused iu_64_floating_t:
    float64_t
    float32_t
    int64_t
    uint64_t

ctypedef fused iu_64_floating_obj_t:
    iu_64_floating_t
    object
