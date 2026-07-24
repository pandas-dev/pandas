// Primitives related to librt.strings.StringWriter that get linked statically
// with compiled modules, instead of being called via a capsule.

#include "stringwriter_extra_ops.h"

// All StringWriter operations are currently implemented as inline functions
// in stringwriter_extra_ops.h, or use the exported capsule API directly.
