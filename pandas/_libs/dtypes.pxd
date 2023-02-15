"""
Common location for shared fused types
"""

from numpy cimport float32_t
from numpy cimport float64_t
from numpy cimport int8_t
from numpy cimport int16_t
from numpy cimport int32_t
from numpy cimport int64_t
from numpy cimport uint8_t
from numpy cimport uint16_t
from numpy cimport uint32_t
from numpy cimport uint64_t

# All numeric types except complex
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

# All numeric types + object, doesn't include complex
ctypedef fused numeric_object_t:
    numeric_t
    object
