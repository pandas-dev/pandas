"""Test that C functions used in primitives are declared in a header such as CPy.h."""

from __future__ import annotations

import glob
import os
import re
import unittest

from mypyc.ir.deps import HeaderDep, SourceDep
from mypyc.ir.ops import PrimitiveDescription
from mypyc.primitives import (
    bytearray_ops,
    bytes_ops,
    dict_ops,
    exc_ops,
    float_ops,
    generic_ops,
    int_ops,
    librt_strings_ops,
    librt_vecs_ops,
    list_ops,
    misc_ops,
    registry,
    set_ops,
    str_ops,
    tuple_ops,
    weakref_ops,
)


class TestHeaderInclusion(unittest.TestCase):
    def test_primitives_included_in_header(self) -> None:
        base_dir = os.path.join(os.path.dirname(__file__), "..", "lib-rt")
        with open(os.path.join(base_dir, "CPy.h")) as f:
            header = f.read()
        with open(os.path.join(base_dir, "pythonsupport.h")) as f:
            header += f.read()

        def check_name(name: str) -> None:
            if name.startswith("CPy"):
                assert re.search(
                    rf"\b{name}\b", header
                ), f'"{name}" is used in mypyc.primitives but not declared in CPy.h'

        all_ops = []
        for values in [
            registry.method_call_ops.values(),
            registry.binary_ops.values(),
            registry.unary_ops.values(),
            registry.function_ops.values(),
        ]:
            for ops in values:
                all_ops.extend(ops)

        for module in [
            bytes_ops,
            str_ops,
            dict_ops,
            list_ops,
            bytearray_ops,
            generic_ops,
            int_ops,
            misc_ops,
            tuple_ops,
            exc_ops,
            float_ops,
            set_ops,
            weakref_ops,
            librt_vecs_ops,
            librt_strings_ops,
        ]:
            for name in dir(module):
                val = getattr(module, name, None)
                if isinstance(val, PrimitiveDescription):
                    all_ops.append(val)

        # Find additional headers via extra C source file dependencies.
        for op in all_ops:
            if op.dependencies:
                for dep in op.dependencies:
                    if isinstance(dep, (SourceDep, HeaderDep)):
                        header_fnam = os.path.join(base_dir, dep.get_header())
                        if os.path.isfile(header_fnam):
                            with open(os.path.join(base_dir, header_fnam)) as f:
                                header += f.read()

        for op in all_ops:
            if op.c_function_name is not None:
                check_name(op.c_function_name)

        primitives_path = os.path.join(os.path.dirname(__file__), "..", "primitives")
        for fnam in glob.glob(f"{primitives_path}/*.py"):
            with open(fnam) as f:
                content = f.read()
            for name in re.findall(r'c_function_name=["\'](CPy[A-Z_a-z0-9]+)', content):
                check_name(name)
