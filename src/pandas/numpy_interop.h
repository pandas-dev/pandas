// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_NUMPY_INTEROP_H
#define PANDAS_NUMPY_INTEROP_H

#include "pandas/status.h"
#include "pandas/types.h"

namespace pandas {

Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type);

} // namespace pandas

#endif // PANDAS_NUMPY_INTEROP_H
