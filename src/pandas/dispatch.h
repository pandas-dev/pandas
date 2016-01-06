// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_DISPATCH_H
#define PANDAS_DISPATCH_H

#include "pandas/types.h"

namespace pandas {

class Status;

Status primitive_type_from_enum(TypeEnum tp_enum, DataType** out);

}

#endif // PANDAS_DISPATCH_H
