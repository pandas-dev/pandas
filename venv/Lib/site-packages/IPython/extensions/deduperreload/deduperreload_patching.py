from __future__ import annotations
import ctypes
import sys
from typing import Any

NOT_FOUND: object = object()
NULL: object = object()
_MAX_FIELD_SEARCH_OFFSET = 50

POINTER_N_BYTES = ctypes.sizeof(ctypes.c_void_p)


class DeduperReloaderPatchingMixin:
    @staticmethod
    def infer_field_offset(
        obj: object,
        field: str,
    ) -> int:
        field_value = getattr(obj, field, NOT_FOUND)
        if field_value is NOT_FOUND:
            return -1
        obj_addr = ctypes.c_void_p.from_buffer(ctypes.py_object(obj)).value
        field_addr = ctypes.c_void_p.from_buffer(ctypes.py_object(field_value)).value
        if obj_addr is None or field_addr is None:
            return -1
        ret = -1
        for offset in range(1, _MAX_FIELD_SEARCH_OFFSET):
            if (
                ctypes.cast(
                    obj_addr + POINTER_N_BYTES * offset, ctypes.POINTER(ctypes.c_void_p)
                ).contents.value
                == field_addr
            ):
                ret = offset
                break
        return ret

    @classmethod
    def try_write_readonly_attr(
        cls,
        obj: object,
        field: str,
        new_value: object,
        offset: int | None = None,
    ) -> None:
        prev_value = getattr(obj, field, NOT_FOUND)
        if prev_value is NOT_FOUND:
            return
        if offset is None:
            offset = cls.infer_field_offset(obj, field)
        if offset == -1:
            return
        obj_addr = ctypes.c_void_p.from_buffer(ctypes.py_object(obj)).value
        if new_value is NULL:
            new_value_addr: int | None = 0
        else:
            new_value_addr = ctypes.c_void_p.from_buffer(
                ctypes.py_object(new_value)
            ).value
        if obj_addr is None or new_value_addr is None:
            return
        if prev_value is not None:
            ctypes.pythonapi.Py_DecRef(ctypes.py_object(prev_value))
        if new_value not in (None, NULL):
            ctypes.pythonapi.Py_IncRef(ctypes.py_object(new_value))
        ctypes.cast(
            obj_addr + POINTER_N_BYTES * offset, ctypes.POINTER(ctypes.c_void_p)
        ).contents.value = new_value_addr

    @classmethod
    def try_patch_readonly_attr(
        cls,
        old: object,
        new: object,
        field: str,
        new_is_value: bool = False,
        offset: int = -1,
    ) -> None:
        old_value = getattr(old, field, NOT_FOUND)
        new_value = new if new_is_value else getattr(new, field, NOT_FOUND)
        if old_value is NOT_FOUND or new_value is NOT_FOUND:
            return
        elif old_value is new_value:
            return
        elif old_value is not None and offset < 0:
            offset = cls.infer_field_offset(old, field)
        elif offset < 0:
            assert not new_is_value
            assert new_value is not None
            offset = cls.infer_field_offset(new, field)
        cls.try_write_readonly_attr(old, field, new_value, offset=offset)

    @classmethod
    def try_patch_attr(
        cls,
        old: object,
        new: object,
        field: str,
        new_is_value: bool = False,
        offset: int = -1,
    ) -> None:
        try:
            setattr(old, field, new if new_is_value else getattr(new, field))
        except (AttributeError, TypeError, ValueError):
            cls.try_patch_readonly_attr(old, new, field, new_is_value, offset)

    @classmethod
    def patch_function(
        cls, to_patch_to: Any, to_patch_from: Any, is_method: bool
    ) -> None:
        new_closure = []
        for freevar, closure_val in zip(
            to_patch_from.__code__.co_freevars or [], to_patch_from.__closure__ or []
        ):
            if (
                callable(closure_val.cell_contents)
                and freevar in to_patch_to.__code__.co_freevars
            ):
                new_closure.append(
                    to_patch_to.__closure__[
                        to_patch_to.__code__.co_freevars.index(freevar)
                    ]
                )
            else:
                new_closure.append(closure_val)
        # lambdas may complain if there is more than one freevar
        cls.try_patch_attr(to_patch_to, to_patch_from, "__code__")
        offset = -1
        if to_patch_to.__closure__ is None and to_patch_from.__closure__ is not None:
            offset = cls.infer_field_offset(to_patch_from, "__closure__")
        if to_patch_to.__closure__ is not None or to_patch_from.__closure__ is not None:
            cls.try_patch_readonly_attr(
                to_patch_to,
                tuple(new_closure) or NULL,
                "__closure__",
                new_is_value=True,
                offset=offset,
            )
        for attr in ("__defaults__", "__kwdefaults__", "__doc__", "__dict__"):
            cls.try_patch_attr(to_patch_to, to_patch_from, attr)
        if is_method:
            cls.try_patch_readonly_attr(to_patch_to, to_patch_from, "__self__")
