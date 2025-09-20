import ctypes
import threading
from ctypes import CFUNCTYPE, c_int, c_int32
from ctypes.util import find_library
import gc
import locale
import os
import platform
import re
import subprocess
import sys
import unittest
from contextlib import contextmanager
from tempfile import mkstemp

from llvmlite import ir
from llvmlite import binding as llvm
from llvmlite.binding import ffi
from llvmlite.tests import TestCase

# arvm7l needs extra ABI symbols to link successfully
if platform.machine() == 'armv7l':
    llvm.load_library_permanently('libgcc_s.so.1')


is_conda_package = unittest.skipUnless(llvm.package_format == "conda",
                                       ("conda package test only, have "
                                        f"{llvm.package_format}"))
is_wheel_package = unittest.skipUnless(llvm.package_format == "wheel",
                                       ("wheel package test only, have "
                                        f"{llvm.package_format}"))

_HAVE_LIEF = False
try:
    import lief  # noqa: F401
    _HAVE_LIEF = True
except ImportError:
    pass


needs_lief = unittest.skipUnless(_HAVE_LIEF, "test needs py-lief package")


def no_de_locale():
    cur = locale.setlocale(locale.LC_ALL)
    try:
        locale.setlocale(locale.LC_ALL, 'de_DE')
    except locale.Error:
        return True
    else:
        return False
    finally:
        locale.setlocale(locale.LC_ALL, cur)


asm_sum = r"""
    ; ModuleID = '<string>'
    source_filename = "asm_sum.c"
    target triple = "{triple}"
    %struct.glob_type = type {{ i64, [2 x i64]}}
    %struct.glob_type_vec = type {{ i64, <2 x i64>}}

    @glob = global i32 0
    @glob_b = global i8 0
    @glob_f = global float 1.5
    @glob_struct = global %struct.glob_type {{i64 0, [2 x i64] [i64 0, i64 0]}}

    define i32 @sum(i32 %.1, i32 %.2) {{
      %.3 = add i32 %.1, %.2
      %.4 = add i32 0, %.3
      ret i32 %.4
    }}
    """

asm_sum2 = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define i32 @sum(i32 %.1, i32 %.2) {{
      %.3 = add i32 %.1, %.2
      ret i32 %.3
    }}
    """

asm_sum3 = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define i64 @sum(i64 %.1, i64 %.2) {{
      %.3 = add i64 %.1, %.2
      %.4 = add i64 5, %.3
      %.5 = add i64 -5, %.4
      ret i64 %.5
    }}
    """

asm_sum4 = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define i32 @sum(i32 %.1, i32 %.2) {{
        %.3 = add i32 %.1, %.2
        ret i32 0
    }}
    """

asm_mul = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"
    @mul_glob = global i32 0

    define i32 @mul(i32 %.1, i32 %.2) {{
      %.3 = mul i32 %.1, %.2
      ret i32 %.3
    }}
    """

asm_square_sum = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"
    @mul_glob = global i32 0

    declare i32 @sum(i32, i32)
    define i32 @square_sum(i32 %.1, i32 %.2) {{
      %.3 = call i32 @sum(i32 %.1, i32 %.2)
      %.4 = mul i32 %.3, %.3
      ret i32 %.4
    }}
    """

asm_getversion = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    declare i8* @Py_GetVersion()

    define void @getversion(i32 %.1, i32 %.2) {{
      %1 = call i8* @Py_GetVersion()
      ret void
    }}
    """

if platform.python_implementation() == 'PyPy':
    asm_getversion = asm_getversion.replace('Py_GetVersion', 'PyPy_GetVersion')

# `fadd` used on integer inputs
asm_parse_error = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define i32 @sum(i32 %.1, i32 %.2) {{
      %.3 = fadd i32 %.1, %.2
      ret i32 %.3
    }}
    """

# "%.bug" definition references itself
asm_verification_fail = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @sum() {{
      %.bug = add i32 1, %.bug
      ret void
    }}
    """

asm_sum_declare = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    declare i32 @sum(i32 %.1, i32 %.2)
    """

asm_vararg_declare = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    declare i32 @vararg(i32 %.1, ...)
    """

asm_double_inaccurate = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @foo() {{
      %const = fadd fp128 0xLF3CB1CCF26FBC178452FB4EC7F91DEAD, 0xL00000000000000000000000000000001
      ret void
    }}
    """  # noqa E501

asm_double_locale = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @foo() {{
      %const = fadd double 0.0, 3.14
      ret void
    }}
    """


asm_inlineasm = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @foo() {{
      call void asm sideeffect "nop", ""()
      ret void
    }}
    """

asm_inlineasm2 = """
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @inlineme() {{
        ret void
    }}

    define i32 @caller(i32 %.1, i32 %.2) {{
    entry:
      %stack = alloca i32
      store i32 %.1, i32* %stack
      br label %main
    main:
      %loaded = load i32, i32* %stack
      %.3 = add i32 %loaded, %.2
      %.4 = add i32 0, %.3
      call void @inlineme()
      ret i32 %.4
    }}
"""

asm_inlineasm3 = """
; ModuleID = 'test.c'
source_filename = "test.c"
target triple = "{triple}"

; Function Attrs: noinline nounwind optnone ssp uwtable
define void @inlineme() noinline !dbg !15 {{
  ret void, !dbg !18
}}

; Function Attrs: noinline nounwind optnone ssp uwtable
define i32 @foo(i32 %0, i32 %1) !dbg !19 {{
  %3 = alloca i32, align 4
  %4 = alloca i32, align 4
  store i32 %0, i32* %3, align 4
  call void @llvm.dbg.declare(metadata i32* %3, metadata !23, metadata !DIExpression()), !dbg !24
  store i32 %1, i32* %4, align 4
  call void @llvm.dbg.declare(metadata i32* %4, metadata !25, metadata !DIExpression()), !dbg !26
  call void @inlineme(), !dbg !27
  %5 = load i32, i32* %3, align 4, !dbg !28
  %6 = load i32, i32* %4, align 4, !dbg !29
  %7 = add nsw i32 %5, %6, !dbg !30
  ret i32 %7, !dbg !31
}}

; Function Attrs: nofree nosync nounwind readnone speculatable willreturn
declare void @llvm.dbg.declare(metadata, metadata, metadata) #1

attributes #1 = {{ nofree nosync nounwind readnone speculatable willreturn }}

!llvm.module.flags = !{{!1, !2, !3, !4, !5, !6, !7, !8, !9, !10}}
!llvm.dbg.cu = !{{!11}}
!llvm.ident = !{{!14}}

!0 = !{{i32 2, !"SDK Version", [2 x i32] [i32 12, i32 3]}}
!1 = !{{i32 7, !"Dwarf Version", i32 4}}
!2 = !{{i32 2, !"Debug Info Version", i32 3}}
!3 = !{{i32 1, !"wchar_size", i32 4}}
!4 = !{{i32 1, !"branch-target-enforcement", i32 0}}
!5 = !{{i32 1, !"sign-return-address", i32 0}}
!6 = !{{i32 1, !"sign-return-address-all", i32 0}}
!7 = !{{i32 1, !"sign-return-address-with-bkey", i32 0}}
!8 = !{{i32 7, !"PIC Level", i32 2}}
!9 = !{{i32 7, !"uwtable", i32 1}}
!10 = !{{i32 7, !"frame-pointer", i32 1}}
!11 = distinct !DICompileUnit(language: DW_LANG_C99, file: !12, producer: "Apple clang version 13.1.6 (clang-1316.0.21.2.3)", isOptimized: false, runtimeVersion: 0, emissionKind: FullDebug, enums: !13, splitDebugInlining: false, nameTableKind: None, sysroot: "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk", sdk: "MacOSX.sdk")
!12 = !DIFile(filename: "test.c", directory: "/")
!13 = !{{}}
!14 = !{{!"Apple clang version 13.1.6 (clang-1316.0.21.2.3)"}}
!15 = distinct !DISubprogram(name: "inlineme", scope: !12, file: !12, line: 1, type: !16, scopeLine: 1, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition, unit: !11, retainedNodes: !13)
!16 = !DISubroutineType(types: !17)
!17 = !{{null}}
!18 = !DILocation(line: 1, column: 22, scope: !15)
!19 = distinct !DISubprogram(name: "foo", scope: !12, file: !12, line: 3, type: !20, scopeLine: 3, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition, unit: !11, retainedNodes: !13)
!20 = !DISubroutineType(types: !21)
!21 = !{{!22, !22, !22}}
!22 = !DIBasicType(name: "int", size: 32, encoding: DW_ATE_signed)
!23 = !DILocalVariable(name: "a", arg: 1, scope: !19, file: !12, line: 3, type: !22)
!24 = !DILocation(line: 3, column: 13, scope: !19)
!25 = !DILocalVariable(name: "b", arg: 2, scope: !19, file: !12, line: 3, type: !22)
!26 = !DILocation(line: 3, column: 20, scope: !19)
!27 = !DILocation(line: 4, column: 5, scope: !19)
!28 = !DILocation(line: 5, column: 12, scope: !19)
!29 = !DILocation(line: 5, column: 16, scope: !19)
!30 = !DILocation(line: 5, column: 14, scope: !19)
!31 = !DILocation(line: 5, column: 5, scope: !19)
"""  # noqa E501

licm_asm = r"""
; ModuleID = "<string>"
target triple = "{triple}"

define double @licm(i32 %0) {{
  %2 = alloca i32, align 4
  %3 = alloca double, align 8
  %4 = alloca i32, align 4
  %5 = alloca double, align 8
  store i32 %0, i32* %2, align 4
  store double 0.000000e+00, double* %3, align 8
  store i32 0, i32* %4, align 4
  br label %6

6:                                                ; preds = %14, %1
  %7 = load i32, i32* %4, align 4
  %8 = load i32, i32* %2, align 4
  %9 = icmp slt i32 %7, %8
  br i1 %9, label %10, label %17

10:                                               ; preds = %6
  store double 7.000000e+00, double* %5, align 8
  %11 = load double, double* %5, align 8
  %12 = load double, double* %3, align 8
  %13 = fadd double %12, %11
  store double %13, double* %3, align 8
  br label %14

14:                                               ; preds = %10
  %15 = load i32, i32* %4, align 4
  %16 = add nsw i32 %15, 1
  store i32 %16, i32* %4, align 4
  br label %6

17:                                               ; preds = %6
  %18 = load double, double* %3, align 8
  ret double %18
}}
"""  # noqa E501

asm_global_ctors = r"""
    ; ModuleID = "<string>"
    target triple = "{triple}"

    @A = global i32 undef

    define void @ctor_A()
    {{
      store i32 10, i32* @A
      ret void
    }}

    define void @dtor_A()
    {{
      store i32 20, i32* @A
      ret void
    }}

    define i32 @foo()
    {{
      %.2 = load i32, i32* @A
      %.3 = add i32 %.2, 2
      ret i32 %.3
    }}

    @llvm.global_ctors = appending global [1 x {{i32, void ()*, i8*}}] [{{i32, void ()*, i8*}} {{i32 0, void ()* @ctor_A, i8* null}}]
    @llvm.global_dtors = appending global [1 x {{i32, void ()*, i8*}}] [{{i32, void ()*, i8*}} {{i32 0, void ()* @dtor_A, i8* null}}]
    """  # noqa E501

asm_ext_ctors = r"""
    ; ModuleID = "<string>"
    target triple = "{triple}"

    @A = external global i32

    define void @ctor_A()
    {{
      store i32 10, i32* @A
      ret void
    }}

    define void @dtor_A()
    {{
      store i32 20, i32* @A
      ret void
    }}

    define i32 @foo()
    {{
      %.2 = load i32, i32* @A
      %.3 = add i32 %.2, 2
      ret i32 %.3
    }}

    @llvm.global_ctors = appending global [1 x {{i32, void ()*, i8*}}] [{{i32, void ()*, i8*}} {{i32 0, void ()* @ctor_A, i8* null}}]
    @llvm.global_dtors = appending global [1 x {{i32, void ()*, i8*}}] [{{i32, void ()*, i8*}} {{i32 0, void ()* @dtor_A, i8* null}}]
    """  # noqa E501


asm_nonalphanum_blocklabel = """; ModuleID = ""
target triple = "unknown-unknown-unknown"
target datalayout = ""

define i32 @"foo"()
{
"<>!*''#":
  ret i32 12345
}
"""  # noqa W291 # trailing space needed for match later


asm_null_constant = r"""
    ; ModuleID = '<string>'
    target triple = "{triple}"

    define void @foo(i64* %.1) {{
      ret void
    }}

    define void @bar() {{
      call void @foo(i64* null)
      ret void
    }}
"""


riscv_asm_ilp32 = [
    'addi\tsp, sp, -16',
    'sw\ta1, 8(sp)',
    'sw\ta2, 12(sp)',
    'fld\tfa5, 8(sp)',
    'fmv.w.x\tfa4, a0',
    'fcvt.d.s\tfa4, fa4',
    'fadd.d\tfa5, fa4, fa5',
    'fsd\tfa5, 8(sp)',
    'lw\ta0, 8(sp)',
    'lw\ta1, 12(sp)',
    'addi\tsp, sp, 16',
    'ret'
]


riscv_asm_ilp32f = [
    'addi\tsp, sp, -16',
    'sw\ta0, 8(sp)',
    'sw\ta1, 12(sp)',
    'fld\tfa5, 8(sp)',
    'fcvt.d.s\tfa4, fa0',
    'fadd.d\tfa5, fa4, fa5',
    'fsd\tfa5, 8(sp)',
    'lw\ta0, 8(sp)',
    'lw\ta1, 12(sp)',
    'addi\tsp, sp, 16',
    'ret'
]


riscv_asm_ilp32d = [
    'fcvt.d.s\tfa5, fa0',
    'fadd.d\tfa0, fa5, fa1',
    'ret'
]


asm_attributes = r"""
declare void @a_readonly_func(i8 *) readonly

declare i8* @a_arg0_return_func(i8* returned, i32*)
"""

asm_alloca_optnone = r"""
define double @foo(i32 %i, double %j) optnone noinline {
    %I = alloca i32		; <i32*> [#uses=4]
    %J = alloca double		; <double*> [#uses=2]
    store i32 %i, i32* %I
    store double %j, double* %J
    %t1 = load i32, i32* %I		; <i32> [#uses=1]
    %t2 = add i32 %t1, 1		; <i32> [#uses=1]
    store i32 %t2, i32* %I
    %t3 = load i32, i32* %I		; <i32> [#uses=1]
    %t4 = sitofp i32 %t3 to double		; <double> [#uses=1]
    %t5 = load double, double* %J		; <double> [#uses=1]
    %t6 = fmul double %t4, %t5		; <double> [#uses=1]
    ret double %t6
}
"""

asm_declaration = r"""
declare void @test_declare(i32* )
"""

# This produces the following output from objdump:
#
# $ objdump -D 632.elf
#
# 632.elf:     file format elf64-x86-64
#
#
# Disassembly of section .text:
#
# 0000000000000000 <__arybo>:
#    0:	48 c1 e2 20          	shl    $0x20,%rdx
#    4:	48 09 c2             	or     %rax,%rdx
#    7:	48 89 d0             	mov    %rdx,%rax
#    a:	48 c1 c0 3d          	rol    $0x3d,%rax
#    e:	48 31 d0             	xor    %rdx,%rax
#   11:	48 b9 01 20 00 04 80 	movabs $0x7010008004002001,%rcx
#   18:	00 10 70
#   1b:	48 0f af c8          	imul   %rax,%rcx

issue_632_elf = \
    "7f454c4602010100000000000000000001003e00010000000000000000000000000000" \
    "0000000000e0000000000000000000000040000000000040000500010048c1e2204809" \
    "c24889d048c1c03d4831d048b90120000480001070480fafc800000000000000000000" \
    "0000000000000000000000000000002f0000000400f1ff000000000000000000000000" \
    "00000000070000001200020000000000000000001f00000000000000002e7465787400" \
    "5f5f617279626f002e6e6f74652e474e552d737461636b002e737472746162002e7379" \
    "6d746162003c737472696e673e00000000000000000000000000000000000000000000" \
    "0000000000000000000000000000000000000000000000000000000000000000000000" \
    "00000000000000001f0000000300000000000000000000000000000000000000a80000" \
    "0000000000380000000000000000000000000000000100000000000000000000000000" \
    "000001000000010000000600000000000000000000000000000040000000000000001f" \
    "000000000000000000000000000000100000000000000000000000000000000f000000" \
    "01000000000000000000000000000000000000005f0000000000000000000000000000" \
    "0000000000000000000100000000000000000000000000000027000000020000000000" \
    "0000000000000000000000000000600000000000000048000000000000000100000002" \
    "00000008000000000000001800000000000000"


issue_632_text = \
    "48c1e2204809c24889d048c1c03d4831d048b90120000480001070480fafc8"


asm_tli_exp2 = r"""
; ModuleID = '<lambda>'
source_filename = "<string>"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-pc-windows-msvc"

declare float @llvm.exp2.f32(float %casted)

define float @foo(i16 %arg) {
entry:
  %casted = sitofp i16 %arg to float
  %ret = call float @llvm.exp2.f32(float %casted)
  ret float %ret
}
"""  # noqa E501

asm_phi_blocks = r"""
; ModuleID = '<string>'
target triple = "{triple}"

define void @foo(i32 %N) {{
  ; unnamed block for testing
  %cmp4 = icmp sgt i32 %N, 0
  br i1 %cmp4, label %for.body, label %for.cond.cleanup

for.cond.cleanup:
  ret void

for.body:
  %i.05 = phi i32 [ %inc, %for.body ], [ 0, %0 ]
  %inc = add nuw nsw i32 %i.05, 1
  %exitcond.not = icmp eq i32 %inc, %N
  br i1 %exitcond.not, label %for.cond.cleanup, label %for.body
}}
"""

asm_cpp_class = r"""
; Source C++
;-----------------------------------------
; class MyClass;
;
; class MyClassDefined{
;     MyClass *member;
;     MyClass *m2;
;     MyClass *m3;
; };
;
; void foo(MyClass *c, MyClassDefined){ }
;-----------------------------------------
; LLVM-IR by: clang -arch arm64 -S -emit-llvm file.cpp
; ModuleID = 'file.cpp'
source_filename = "class.cpp"
target datalayout = "e-m:o-i64:64-i128:128-n32:64-S128"
target triple = "arm64-apple-macosx13.3.0"

%class.MyClass = type opaque
%class.MyClassDefined = type { %class.MyClass*, %class.MyClass*, %class.MyClass* }

; Function Attrs: noinline nounwind optnone ssp uwtable(sync)
define void @_Z3fooP7MyClass14MyClassDefined(%class.MyClass* noundef %0, %class.MyClassDefined* noundef %1) {
  %3 = alloca %class.MyClass*, align 8
  store %class.MyClass* %0, %class.MyClass** %3, align 8
  ret void
}

!llvm.module.flags = !{!0, !1, !2, !3, !4, !5, !6, !7, !8}
!llvm.ident = !{!9}

!0 = !{i32 2, !"SDK Version", [2 x i32] [i32 13, i32 3]}
!1 = !{i32 1, !"wchar_size", i32 4}
!2 = !{i32 8, !"branch-target-enforcement", i32 0}
!3 = !{i32 8, !"sign-return-address", i32 0}
!4 = !{i32 8, !"sign-return-address-all", i32 0}
!5 = !{i32 8, !"sign-return-address-with-bkey", i32 0}
!6 = !{i32 7, !"PIC Level", i32 2}
!7 = !{i32 7, !"uwtable", i32 1}
!8 = !{i32 7, !"frame-pointer", i32 1}
!9 = !{!"Apple clang version 14.0.3 (clang-1403.0.22.14.1)"}

""" # noqa

asm_cpp_vector = r"""; Source C++
;-----------------------------------------

; struct Vector2D{
;     float x, y;
; };
;
; void foo(Vector2D vec, Vector2D *out) {
;     *out = vec;
; }
;-----------------------------------------
; LLVM-IR by: clang -arch x86_64 -S -emit-llvm file.cpp
; ModuleID = 'file.cpp'
source_filename = "class.cpp"
target datalayout = "e-m:o-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx13.3.0"

%struct.Vector2D = type { float, float }

; Function Attrs: noinline nounwind optnone ssp uwtable
define void @_Z3foo8Vector2DPS_(<2 x float> %0, %struct.Vector2D* noundef %1) #0 {
  %3 = alloca %struct.Vector2D, align 4
  %4 = alloca %struct.Vector2D*, align 8
  %5 = bitcast %struct.Vector2D* %3 to <2 x float>*
  store <2 x float> %0, <2 x float>* %5, align 4
  store %struct.Vector2D* %1, %struct.Vector2D** %4, align 8
  %6 = load %struct.Vector2D*, %struct.Vector2D** %4, align 8
  %7 = bitcast %struct.Vector2D* %6 to i8*
  %8 = bitcast %struct.Vector2D* %3 to i8*
  call void @llvm.memcpy.p0i8.p0i8.i64(i8* align 4 %7, i8* align 4 %8, i64 8, i1 false)
  ret void
}

; Function Attrs: argmemonly nofree nounwind willreturn
declare void @llvm.memcpy.p0i8.p0i8.i64(i8* noalias nocapture writeonly, i8* noalias nocapture readonly, i64, i1 immarg) #1

attributes #0 = { noinline nounwind optnone ssp uwtable "darwin-stkchk-strong-link" "frame-pointer"="all" "min-legal-vector-width"="64" "no-trapping-math"="true" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "tune-cpu"="generic" }
attributes #1 = { argmemonly nofree nounwind willreturn }

!llvm.module.flags = !{!0, !1, !2, !3, !4}
!llvm.ident = !{!5}

!0 = !{i32 2, !"SDK Version", [2 x i32] [i32 13, i32 3]}
!1 = !{i32 1, !"wchar_size", i32 4}
!2 = !{i32 7, !"PIC Level", i32 2}
!3 = !{i32 7, !"uwtable", i32 2}
!4 = !{i32 7, !"frame-pointer", i32 2}
!5 = !{!"Apple clang version 14.0.3 (clang-1403.0.22.14.1)"}

""" # noqa


class BaseTest(TestCase):

    def setUp(self):
        llvm.initialize_native_target()
        llvm.initialize_native_asmprinter()
        gc.collect()
        self.old_garbage = gc.garbage[:]
        gc.garbage[:] = []

    def tearDown(self):
        # Test that no uncollectable objects were created
        # (llvmlite objects have a __del__ so a reference cycle could
        # create some).
        gc.collect()
        self.assertEqual(gc.garbage, [])
        # This will probably put any existing garbage in gc.garbage again
        del self.old_garbage

    def module(self, asm=asm_sum, context=None):
        asm = asm.format(triple=llvm.get_default_triple())
        mod = llvm.parse_assembly(asm, context)
        return mod

    def glob(self, name='glob', mod=None):
        if mod is None:
            mod = self.module()
        return mod.get_global_variable(name)

    def target_machine(self, *, jit):
        target = llvm.Target.from_default_triple()
        return target.create_target_machine(jit=jit)


class TestDependencies(BaseTest):
    """
    Test DLL dependencies are within a certain expected set.
    """

    @unittest.skipUnless(sys.platform.startswith('linux'),
                         "Linux-specific test")
    @unittest.skipUnless(os.environ.get('LLVMLITE_DIST_TEST'),
                         "Distribution-specific test")
    def test_linux(self):
        lib_path = ffi.lib._name
        env = os.environ.copy()
        env['LANG'] = 'C'
        p = subprocess.Popen(["objdump", "-p", lib_path],
                             stdout=subprocess.PIPE, env=env)
        out, _ = p.communicate()
        self.assertEqual(0, p.returncode)
        # Parse library dependencies
        lib_pat = re.compile(r'^([+-_a-zA-Z0-9]+)\.so(?:\.\d+){0,3}$')
        deps = set()
        for line in out.decode().splitlines():
            parts = line.split()
            if parts and parts[0] == 'NEEDED':
                dep = parts[1]
                m = lib_pat.match(dep)
                if len(parts) != 2 or not m:
                    self.fail("invalid NEEDED line: %r" % (line,))
                deps.add(m.group(1))
        # Sanity check that our dependencies were parsed ok
        if 'libc' not in deps or 'libpthread' not in deps:
            self.fail("failed parsing dependencies? got %r" % (deps,))
        # Ensure all dependencies are expected
        allowed = set(['librt', 'libdl', 'libpthread', 'libz', 'libm',
                       'libgcc_s', 'libc', 'ld-linux', 'ld64', 'libzstd',
                       'libstdc++'])
        if platform.python_implementation() == 'PyPy':
            allowed.add('libtinfo')

        fails = []
        for dep in deps:
            if not dep.startswith('ld-linux-') and dep not in allowed:
                fails.append(dep)
        if len(fails) == 1:
            self.fail("unexpected dependency %r in %r" % (fails[0], deps))
        elif len(fails) > 1:
            self.fail("unexpected dependencies %r in %r" % (fails, deps))
        else:
            pass  # test passes


class TestRISCVABI(BaseTest):
    """
    Test calling convention of floating point arguments of RISC-V
    using different ABI.
    """
    triple = "riscv32-unknown-linux"

    def setUp(self):
        super().setUp()
        llvm.initialize_all_targets()
        llvm.initialize_all_asmprinters()

    def check_riscv_target(self):
        try:
            llvm.Target.from_triple(self.triple)
        except RuntimeError as e:
            if "No available targets are compatible with triple" in str(e):
                self.skipTest("RISCV target unsupported by linked LLVM.")
            else:
                raise e

    def riscv_target_machine(self, **kwarg):
        lltarget = llvm.Target.from_triple(self.triple)
        return lltarget.create_target_machine(**kwarg)

    def fpadd_ll_module(self):
        f64 = ir.DoubleType()
        f32 = ir.FloatType()
        fnty = ir.FunctionType(f64, (f32, f64))
        module = ir.Module()
        func = ir.Function(module, fnty, name="fpadd")
        block = func.append_basic_block()
        builder = ir.IRBuilder(block)
        a, b = func.args
        arg0 = builder.fpext(a, f64)
        result = builder.fadd(arg0, b)
        builder.ret(result)

        llmod = llvm.parse_assembly(str(module))
        llmod.verify()
        return llmod

    def break_up_asm(self, asm):
        asm_list = []
        for line in asm.splitlines():
            s_line = line.strip()
            if not (s_line.startswith(".") or s_line.startswith("fpadd")
                    or s_line == ""):
                asm_list.append(s_line)
        return asm_list

    def test_rv32d_ilp32(self):
        self.check_riscv_target()
        llmod = self.fpadd_ll_module()
        target = self.riscv_target_machine(features="+f,+d", abiname="ilp32")
        self.assertEqual(self.break_up_asm(target.emit_assembly(llmod)),
                         riscv_asm_ilp32)

    def test_rv32d_ilp32f(self):
        self.check_riscv_target()
        llmod = self.fpadd_ll_module()
        target = self.riscv_target_machine(features="+f,+d", abiname="ilp32f")
        self.assertEqual(self.break_up_asm(target.emit_assembly(llmod)),
                         riscv_asm_ilp32f)

    def test_rv32d_ilp32d(self):
        self.check_riscv_target()
        llmod = self.fpadd_ll_module()
        target = self.riscv_target_machine(features="+f,+d", abiname="ilp32d")
        self.assertEqual(self.break_up_asm(target.emit_assembly(llmod)),
                         riscv_asm_ilp32d)


class TestMisc(BaseTest):
    """
    Test miscellaneous functions in llvm.binding.
    """

    def test_parse_assembly(self):
        self.module(asm_sum)

    def test_parse_assembly_error(self):
        with self.assertRaises(RuntimeError) as cm:
            self.module(asm_parse_error)
        s = str(cm.exception)
        self.assertIn("parsing error", s)
        self.assertIn("invalid operand type", s)

    def test_nonalphanum_block_name(self):
        mod = ir.Module()
        ft = ir.FunctionType(ir.IntType(32), [])
        fn = ir.Function(mod, ft, "foo")
        bd = ir.IRBuilder(fn.append_basic_block(name="<>!*''#"))
        bd.ret(ir.Constant(ir.IntType(32), 12345))
        asm = str(mod)
        self.assertEqual(asm, asm_nonalphanum_blocklabel)

    def test_global_context(self):
        gcontext1 = llvm.context.get_global_context()
        gcontext2 = llvm.context.get_global_context()
        assert gcontext1 == gcontext2

    def test_dylib_symbols(self):
        llvm.add_symbol("__xyzzy", 1234)
        llvm.add_symbol("__xyzzy", 5678)
        addr = llvm.address_of_symbol("__xyzzy")
        self.assertEqual(addr, 5678)
        addr = llvm.address_of_symbol("__foobar")
        self.assertIs(addr, None)

    def test_get_default_triple(self):
        triple = llvm.get_default_triple()
        self.assertIsInstance(triple, str)
        self.assertTrue(triple)

    def test_get_process_triple(self):
        # Sometimes we get synonyms for PPC
        def normalize_ppc(arch):
            if arch == 'powerpc64le':
                return 'ppc64le'
            else:
                return arch

        triple = llvm.get_process_triple()
        default = llvm.get_default_triple()
        self.assertIsInstance(triple, str)
        self.assertTrue(triple)

        default_arch = normalize_ppc(default.split('-')[0])
        triple_arch = normalize_ppc(triple.split('-')[0])
        # Arch must be equal
        self.assertEqual(default_arch, triple_arch)

    def test_get_host_cpu_features(self):
        features = llvm.get_host_cpu_features()
        # Check the content of `features`
        self.assertIsInstance(features, dict)
        self.assertIsInstance(features, llvm.FeatureMap)
        for k, v in features.items():
            self.assertIsInstance(k, str)
            self.assertTrue(k)  # single feature string cannot be empty
            self.assertIsInstance(v, bool)
        self.assertIsInstance(features.flatten(), str)

        re_term = r"[+\-][a-zA-Z0-9\._-]+"
        regex = r"^({0}|{0}(,{0})*)?$".format(re_term)
        # quick check for our regex
        self.assertIsNotNone(re.match(regex, ""))
        self.assertIsNotNone(re.match(regex, "+aa"))
        self.assertIsNotNone(re.match(regex, "+a,-bb"))
        # check CpuFeature.flatten()
        if len(features) == 0:
            self.assertEqual(features.flatten(), "")
        else:
            self.assertIsNotNone(re.match(regex, features.flatten()))

    def test_get_host_cpu_name(self):
        cpu = llvm.get_host_cpu_name()
        self.assertIsInstance(cpu, str)
        self.assertTrue(cpu)

    def test_initfini(self):
        code = """if 1:
            from llvmlite import binding as llvm

            llvm.initialize_native_target()
            llvm.initialize_native_asmprinter()
            llvm.initialize_all_targets()
            llvm.initialize_all_asmprinters()
            llvm.shutdown()
            """
        subprocess.check_call([sys.executable, "-c", code])

    def test_deprecated_init(self):
        regex = r"llvmlite.binding.initialize\(\) is deprecated"
        with self.assertRaisesRegex(RuntimeError, expected_regex=regex):
            llvm.initialize()

    def test_set_option(self):
        # We cannot set an option multiple times (LLVM would exit() the
        # process), so run the code in a subprocess.
        code = """if 1:
            from llvmlite import binding as llvm

            llvm.set_option("progname", "-debug-pass=Disabled")
            """
        subprocess.check_call([sys.executable, "-c", code])

    def test_version(self):
        major, minor, patch = llvm.llvm_version_info
        # one of these can be valid
        valid = (20,)
        self.assertIn(major, valid)
        self.assertIn(patch, range(9))

    def test_check_jit_execution(self):
        llvm.check_jit_execution()

    @unittest.skipIf(no_de_locale(), "Locale not available")
    def test_print_double_locale(self):
        m = self.module(asm_double_locale)
        expect = str(m)
        # Change the locale so that comma is used as decimal-point
        # to trigger the LLVM bug (llvmlite issue #80)
        locale.setlocale(locale.LC_ALL, 'de_DE')
        # The LLVM bug is trigged by print the module with double constant
        got = str(m)
        # Changing the locale should not affect the LLVM IR
        self.assertEqual(expect, got)

    def test_no_accidental_warnings(self):
        code = "from llvmlite import binding"
        flags = "-Werror"
        cmdargs = [sys.executable, flags, "-c", code]
        subprocess.check_call(cmdargs)


class TestModuleRef(BaseTest):

    def test_str(self):
        mod = self.module()
        s = str(mod).strip()
        self.assertTrue(s.startswith('; ModuleID ='), s)

    def test_close(self):
        mod = self.module()
        str(mod)
        mod.close()
        with self.assertRaises(ctypes.ArgumentError):
            str(mod)
        mod.close()

    def test_with(self):
        mod = self.module()
        str(mod)
        with mod:
            str(mod)
        with self.assertRaises(ctypes.ArgumentError):
            str(mod)
        with self.assertRaises(RuntimeError):
            with mod:
                pass

    def test_name(self):
        mod = self.module()
        mod.name = "foo"
        self.assertEqual(mod.name, "foo")
        mod.name = "bar"
        self.assertEqual(mod.name, "bar")

    def test_source_file(self):
        mod = self.module()
        self.assertEqual(mod.source_file, "asm_sum.c")

    def test_data_layout(self):
        mod = self.module()
        s = mod.data_layout
        self.assertIsInstance(s, str)
        mod.data_layout = s
        self.assertEqual(s, mod.data_layout)

    def test_triple(self):
        mod = self.module()
        s = mod.triple
        self.assertEqual(s, llvm.get_default_triple())
        mod.triple = ''
        self.assertEqual(mod.triple, '')

    def test_verify(self):
        # Verify successful
        mod = self.module()
        self.assertIs(mod.verify(), None)
        # Verify failed
        mod = self.module(asm_verification_fail)
        with self.assertRaises(RuntimeError) as cm:
            mod.verify()
        s = str(cm.exception)
        self.assertIn("%.bug = add i32 1, %.bug", s)

    def test_get_function(self):
        mod = self.module()
        fn = mod.get_function("sum")
        self.assertIsInstance(fn, llvm.ValueRef)
        self.assertEqual(fn.name, "sum")

        with self.assertRaises(NameError):
            mod.get_function("foo")

        # Check that fn keeps the module instance alive
        del mod
        str(fn.module)

    def test_get_struct_type(self):
        mod = self.module()
        st_ty = mod.get_struct_type("struct.glob_type")
        self.assertEqual(st_ty.name, "struct.glob_type")
        # also match struct names of form "%struct.glob_type.{some_index}"
        self.assertIsNotNone(re.match(
            r'%struct\.glob_type(\.[\d]+)? = type { i64, \[2 x i64\] }',
            str(st_ty)))

        with self.assertRaises(NameError):
            mod.get_struct_type("struct.doesnt_exist")

    def test_get_global_variable(self):
        mod = self.module()
        gv = mod.get_global_variable("glob")
        self.assertIsInstance(gv, llvm.ValueRef)
        self.assertEqual(gv.name, "glob")

        with self.assertRaises(NameError):
            mod.get_global_variable("bar")

        # Check that gv keeps the module instance alive
        del mod
        str(gv.module)

    def test_global_variables(self):
        mod = self.module()
        it = mod.global_variables
        del mod
        globs = sorted(it, key=lambda value: value.name)
        self.assertEqual(len(globs), 4)
        self.assertEqual([g.name for g in globs],
                         ["glob", "glob_b", "glob_f", "glob_struct"])

    def test_functions(self):
        mod = self.module()
        it = mod.functions
        del mod
        funcs = list(it)
        self.assertEqual(len(funcs), 1)
        self.assertEqual(funcs[0].name, "sum")

    def test_structs(self):
        mod = self.module()
        it = mod.struct_types
        del mod
        structs = list(it)
        self.assertEqual(len(structs), 1)
        self.assertIsNotNone(re.match(r'struct\.glob_type(\.[\d]+)?',
                                      structs[0].name))
        self.assertIsNotNone(re.match(
            r'%struct\.glob_type(\.[\d]+)? = type { i64, \[2 x i64\] }',
            str(structs[0])))

    def test_link_in(self):
        dest = self.module()
        src = self.module(asm_mul)
        dest.link_in(src)
        self.assertEqual(
            sorted(f.name for f in dest.functions), ["mul", "sum"])
        dest.get_function("mul")
        dest.close()
        with self.assertRaises(ctypes.ArgumentError):
            src.get_function("mul")

    def test_link_in_preserve(self):
        dest = self.module()
        src2 = self.module(asm_mul)
        dest.link_in(src2, preserve=True)
        self.assertEqual(
            sorted(f.name for f in dest.functions), ["mul", "sum"])
        dest.close()
        self.assertEqual(sorted(f.name for f in src2.functions), ["mul"])
        src2.get_function("mul")

    def test_link_in_error(self):
        # Raise an error by trying to link two modules with the same global
        # definition "sum".
        dest = self.module()
        src = self.module(asm_sum2)
        with self.assertRaises(RuntimeError) as cm:
            dest.link_in(src)
        self.assertIn("symbol multiply defined", str(cm.exception))

    def test_as_bitcode(self):
        mod = self.module()
        bc = mod.as_bitcode()
        # Refer to http://llvm.org/docs/doxygen/html/ReaderWriter_8h_source.html#l00064  # noqa E501
        # and http://llvm.org/docs/doxygen/html/ReaderWriter_8h_source.html#l00092  # noqa E501
        bitcode_wrapper_magic = b'\xde\xc0\x17\x0b'
        bitcode_magic = b'BC'
        self.assertTrue(bc.startswith(bitcode_magic) or
                        bc.startswith(bitcode_wrapper_magic))

    def test_parse_bitcode_error(self):
        with self.assertRaises(RuntimeError) as cm:
            llvm.parse_bitcode(b"")
        self.assertIn("LLVM bitcode parsing error", str(cm.exception))
        self.assertIn(
            "file too small to contain bitcode header", str(cm.exception),
        )

    def test_bitcode_roundtrip(self):
        # create a new context to avoid struct renaming
        context1 = llvm.create_context()
        bc = self.module(context=context1).as_bitcode()
        context2 = llvm.create_context()
        mod = llvm.parse_bitcode(bc, context2)
        self.assertEqual(mod.as_bitcode(), bc)

        mod.get_function("sum")
        mod.get_global_variable("glob")

    def test_cloning(self):
        m = self.module()
        cloned = m.clone()
        self.assertIsNot(cloned, m)
        self.assertEqual(cloned.as_bitcode(), m.as_bitcode())


class JITTestMixin(object):
    """
    Mixin for ExecutionEngine tests.
    """

    def get_sum(self, ee, func_name="sum"):
        ee.finalize_object()
        cfptr = ee.get_function_address(func_name)
        self.assertTrue(cfptr)
        return CFUNCTYPE(c_int, c_int, c_int)(cfptr)

    def test_run_code(self):
        mod = self.module()
        with self.jit(mod) as ee:
            cfunc = self.get_sum(ee)
            res = cfunc(2, -5)
            self.assertEqual(-3, res)

    def test_close(self):
        ee = self.jit(self.module())
        ee.close()
        ee.close()
        with self.assertRaises(ctypes.ArgumentError):
            ee.finalize_object()

    def test_with(self):
        ee = self.jit(self.module())
        with ee:
            pass
        with self.assertRaises(RuntimeError):
            with ee:
                pass
        with self.assertRaises(ctypes.ArgumentError):
            ee.finalize_object()

    def test_module_lifetime(self):
        mod = self.module()
        ee = self.jit(mod)
        ee.close()
        mod.close()

    def test_module_lifetime2(self):
        mod = self.module()
        ee = self.jit(mod)
        mod.close()
        ee.close()

    def test_add_module(self):
        ee = self.jit(self.module())
        mod = self.module(asm_mul)
        ee.add_module(mod)
        with self.assertRaises(KeyError):
            ee.add_module(mod)
        self.assertFalse(mod.closed)
        ee.close()
        self.assertTrue(mod.closed)

    def test_add_module_lifetime(self):
        ee = self.jit(self.module())
        mod = self.module(asm_mul)
        ee.add_module(mod)
        mod.close()
        ee.close()

    def test_add_module_lifetime2(self):
        ee = self.jit(self.module())
        mod = self.module(asm_mul)
        ee.add_module(mod)
        ee.close()
        mod.close()

    def test_remove_module(self):
        ee = self.jit(self.module())
        mod = self.module(asm_mul)
        ee.add_module(mod)
        ee.remove_module(mod)
        with self.assertRaises(KeyError):
            ee.remove_module(mod)
        self.assertFalse(mod.closed)
        ee.close()
        self.assertFalse(mod.closed)

    def test_target_data(self):
        mod = self.module()
        ee = self.jit(mod)
        td = ee.target_data
        # A singleton is returned
        self.assertIs(ee.target_data, td)
        str(td)
        del mod, ee
        str(td)

    def test_target_data_abi_enquiries(self):
        mod = self.module()
        ee = self.jit(mod)
        td = ee.target_data
        gv_i32 = mod.get_global_variable("glob")
        gv_i8 = mod.get_global_variable("glob_b")
        gv_struct = mod.get_global_variable("glob_struct")
        # A global is a pointer, it has the ABI size of a pointer
        pointer_size = 4 if sys.maxsize < 2 ** 32 else 8
        for g in (gv_i32, gv_i8, gv_struct):
            self.assertEqual(td.get_abi_size(g.type), pointer_size)

        self.assertEqual(td.get_abi_size(gv_i32.global_value_type), 4)
        self.assertEqual(td.get_abi_alignment(gv_i32.global_value_type), 4)

        self.assertEqual(td.get_abi_size(gv_i8.global_value_type), 1)
        self.assertIn(td.get_abi_alignment(gv_i8.global_value_type), (1, 2, 4))

        self.assertEqual(td.get_abi_size(gv_struct.global_value_type), 24)
        self.assertIn(td.get_abi_alignment(gv_struct.global_value_type), (4, 8))

    def test_object_cache_notify(self):
        notifies = []

        def notify(mod, buf):
            notifies.append((mod, buf))

        mod = self.module()
        ee = self.jit(mod)
        ee.set_object_cache(notify)

        self.assertEqual(len(notifies), 0)
        cfunc = self.get_sum(ee)
        cfunc(2, -5)
        self.assertEqual(len(notifies), 1)
        # The right module object was found
        self.assertIs(notifies[0][0], mod)
        self.assertIsInstance(notifies[0][1], bytes)

        notifies[:] = []
        mod2 = self.module(asm_mul)
        ee.add_module(mod2)
        cfunc = self.get_sum(ee, "mul")
        self.assertEqual(len(notifies), 1)
        # The right module object was found
        self.assertIs(notifies[0][0], mod2)
        self.assertIsInstance(notifies[0][1], bytes)

    def test_object_cache_getbuffer(self):
        notifies = []
        getbuffers = []

        def notify(mod, buf):
            notifies.append((mod, buf))

        def getbuffer(mod):
            getbuffers.append(mod)

        mod = self.module()
        ee = self.jit(mod)
        ee.set_object_cache(notify, getbuffer)

        # First return None from getbuffer(): the object is compiled normally
        self.assertEqual(len(notifies), 0)
        self.assertEqual(len(getbuffers), 0)
        cfunc = self.get_sum(ee)
        self.assertEqual(len(notifies), 1)
        self.assertEqual(len(getbuffers), 1)
        self.assertIs(getbuffers[0], mod)
        sum_buffer = notifies[0][1]

        # Recreate a new EE, and use getbuffer() to return the previously
        # compiled object.

        def getbuffer_successful(mod):
            getbuffers.append(mod)
            return sum_buffer

        notifies[:] = []
        getbuffers[:] = []
        # Use another source module to make sure it is ignored
        mod = self.module(asm_mul)
        ee = self.jit(mod)
        ee.set_object_cache(notify, getbuffer_successful)

        self.assertEqual(len(notifies), 0)
        self.assertEqual(len(getbuffers), 0)
        cfunc = self.get_sum(ee)
        self.assertEqual(cfunc(2, -5), -3)
        self.assertEqual(len(notifies), 0)
        self.assertEqual(len(getbuffers), 1)


class JITWithTMTestMixin(JITTestMixin):

    def test_emit_assembly(self):
        """Test TargetMachineRef.emit_assembly()"""
        target_machine = self.target_machine(jit=True)
        mod = self.module()
        ee = self.jit(mod, target_machine)  # noqa F841 # Keeps pointers alive
        raw_asm = target_machine.emit_assembly(mod)
        self.assertIn("sum", raw_asm)
        target_machine.set_asm_verbosity(True)
        raw_asm_verbose = target_machine.emit_assembly(mod)
        self.assertIn("sum", raw_asm)
        self.assertNotEqual(raw_asm, raw_asm_verbose)

    def test_emit_object(self):
        """Test TargetMachineRef.emit_object()"""
        target_machine = self.target_machine(jit=True)
        mod = self.module()
        ee = self.jit(mod, target_machine)  # noqa F841 # Keeps pointers alive
        code_object = target_machine.emit_object(mod)
        self.assertIsInstance(code_object, bytes)
        if sys.platform.startswith('linux'):
            # Sanity check
            self.assertIn(b"ELF", code_object[:10])


class TestMCJit(BaseTest, JITWithTMTestMixin):
    """
    Test JIT engines created with create_mcjit_compiler().
    """

    def jit(self, mod, target_machine=None):
        if target_machine is None:
            target_machine = self.target_machine(jit=True)
        return llvm.create_mcjit_compiler(mod, target_machine)


# There are some memory corruption issues with OrcJIT on AArch64 - see Issue
# #1000. Since OrcJIT is experimental, and we don't test regularly during
# llvmlite development on non-x86 platforms, it seems safest to skip these
# tests on non-x86 platforms.
# After LLVM20 upgrades skip on X86 too as ORCJit complains about missing
# JITDyLib symbol.
# TODO: Investigate this further.
@unittest.skip("OrcJIT support is experimental")
class TestOrcLLJIT(BaseTest):

    def jit(self, asm=asm_sum, func_name="sum", target_machine=None,
            add_process=False, func_type=CFUNCTYPE(c_int, c_int, c_int),
            suppress_errors=False):
        lljit = llvm.create_lljit_compiler(target_machine,
                                           use_jit_link=False,
                                           suppress_errors=suppress_errors)
        builder = llvm.JITLibraryBuilder()
        if add_process:
            builder.add_current_process()
        rt = builder\
            .add_ir(asm.format(triple=llvm.get_default_triple()))\
            .export_symbol(func_name)\
            .link(lljit, func_name)
        cfptr = rt[func_name]
        self.assertTrue(cfptr)
        self.assertEqual(func_name, rt.name)
        return lljit, rt, func_type(cfptr)

    # From test_dylib_symbols
    def test_define_symbol(self):
        lljit = llvm.create_lljit_compiler()
        rt = llvm.JITLibraryBuilder().import_symbol("__xyzzy", 1234)\
            .export_symbol("__xyzzy").link(lljit, "foo")
        self.assertEqual(rt["__xyzzy"], 1234)

    def test_lookup_undefined_symbol_fails(self):
        lljit = llvm.create_lljit_compiler()
        with self.assertRaisesRegex(RuntimeError, 'No such library'):
            lljit.lookup("foo", "__foobar")
        rt = llvm.JITLibraryBuilder().import_symbol("__xyzzy", 1234)\
            .export_symbol("__xyzzy").link(lljit, "foo")
        self.assertNotEqual(rt["__xyzzy"], 0)
        with self.assertRaisesRegex(RuntimeError,
                                    'Symbols not found.*__foobar'):
            lljit.lookup("foo", "__foobar")

    def test_jit_link(self):
        if sys.platform == "win32":
            with self.assertRaisesRegex(RuntimeError,
                                        'JITLink .* Windows'):
                llvm.create_lljit_compiler(use_jit_link=True)
        else:
            self.assertIsNotNone(llvm.create_lljit_compiler(use_jit_link=True))

    def test_run_code(self):
        (lljit, rt, cfunc) = self.jit()
        with lljit:
            res = cfunc(2, -5)
            self.assertEqual(-3, res)

    def test_close(self):
        (lljit, rt, cfunc) = self.jit()
        lljit.close()
        lljit.close()
        with self.assertRaises(AssertionError):
            lljit.lookup("foo", "fn")

    def test_with(self):
        (lljit, rt, cfunc) = self.jit()
        with lljit:
            pass
        with self.assertRaises(RuntimeError):
            with lljit:
                pass
        with self.assertRaises(AssertionError):
            lljit.lookup("foo", "fn")

    def test_add_ir_module(self):
        (lljit, rt_sum, cfunc_sum) = self.jit()
        rt_mul = llvm.JITLibraryBuilder() \
            .add_ir(asm_mul.format(triple=llvm.get_default_triple())) \
            .export_symbol("mul") \
            .link(lljit, "mul")
        res = CFUNCTYPE(c_int, c_int, c_int)(rt_mul["mul"])(2, -5)
        self.assertEqual(-10, res)
        self.assertNotEqual(lljit.lookup("sum", "sum")["sum"], 0)
        self.assertNotEqual(lljit.lookup("mul", "mul")["mul"], 0)
        with self.assertRaises(RuntimeError):
            lljit.lookup("sum", "mul")
        with self.assertRaises(RuntimeError):
            lljit.lookup("mul", "sum")

    def test_remove_module(self):
        (lljit, rt_sum, _) = self.jit()
        del rt_sum
        gc.collect()
        with self.assertRaises(RuntimeError):
            lljit.lookup("sum", "sum")
        lljit.close()

    def test_lib_depends(self):
        (lljit, rt_sum, cfunc_sum) = self.jit()
        rt_mul = llvm.JITLibraryBuilder() \
            .add_ir(asm_square_sum.format(triple=llvm.get_default_triple())) \
            .export_symbol("square_sum") \
            .add_jit_library("sum") \
            .link(lljit, "square_sum")
        res = CFUNCTYPE(c_int, c_int, c_int)(rt_mul["square_sum"])(2, -5)
        self.assertEqual(9, res)

    def test_target_data(self):
        (lljit, rt, _) = self.jit()
        td = lljit.target_data
        # A singleton is returned
        self.assertIs(lljit.target_data, td)
        str(td)
        del lljit
        str(td)

    def test_global_ctors_dtors(self):
        # test issue #303
        # (https://github.com/numba/llvmlite/issues/303)
        shared_value = c_int32(0)
        lljit = llvm.create_lljit_compiler()
        builder = llvm.JITLibraryBuilder()
        rt = builder \
            .add_ir(asm_ext_ctors.format(triple=llvm.get_default_triple())) \
            .import_symbol("A", ctypes.addressof(shared_value)) \
            .export_symbol("foo") \
            .link(lljit, "foo")
        foo = rt["foo"]
        self.assertTrue(foo)
        self.assertEqual(CFUNCTYPE(c_int)(foo)(), 12)
        del rt
        self.assertNotEqual(shared_value.value, 20)

    def test_lookup_current_process_symbol_fails(self):
        # An attempt to lookup a symbol in the current process (Py_GetVersion,
        # in this case) should fail with an appropriate error if we have not
        # enabled searching the current process for symbols.
        msg = 'Failed to materialize symbols:.*getversion'
        with self.assertRaisesRegex(RuntimeError, msg):
            self.jit(asm_getversion, "getversion", suppress_errors=True)

    def test_lookup_current_process_symbol(self):
        self.jit(asm_getversion, "getversion", None, True)

    def test_thread_safe(self):
        lljit = llvm.create_lljit_compiler()
        llvm_ir = asm_sum.format(triple=llvm.get_default_triple())

        def compile_many(i):
            def do_work():
                tracking = []
                for c in range(50):
                    tracking.append(llvm.JITLibraryBuilder()
                                    .add_ir(llvm_ir)
                                    .export_symbol("sum")
                                    .link(lljit, f"sum_{i}_{c}"))

            return do_work

        ths = [threading.Thread(target=compile_many(i))
               for i in range(os.cpu_count())]
        for th in ths:
            th.start()
        for th in ths:
            th.join()

    def test_add_object_file(self):
        target_machine = self.target_machine(jit=False)
        mod = self.module()
        lljit = llvm.create_lljit_compiler(target_machine)
        rt = llvm.JITLibraryBuilder()\
            .add_object_img(target_machine.emit_object(mod))\
            .export_symbol("sum")\
            .link(lljit, "sum")
        sum = CFUNCTYPE(c_int, c_int, c_int)(rt["sum"])
        self.assertEqual(sum(2, 3), 5)

    def test_add_object_file_from_filesystem(self):
        target_machine = self.target_machine(jit=False)
        mod = self.module()
        obj_bin = target_machine.emit_object(mod)
        temp_desc, temp_path = mkstemp()

        try:
            with os.fdopen(temp_desc, "wb") as f:
                f.write(obj_bin)
            lljit = llvm.create_lljit_compiler(target_machine)
            rt = llvm.JITLibraryBuilder() \
                .add_object_file(temp_path) \
                .export_symbol("sum") \
                .link(lljit, "sum")
            sum = CFUNCTYPE(c_int, c_int, c_int)(rt["sum"])
            self.assertEqual(sum(2, 3), 5)
        finally:
            os.unlink(temp_path)


class TestValueRef(BaseTest):

    def test_str(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        self.assertEqual(str(glob), "@glob = global i32 0")

    def test_name(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        self.assertEqual(glob.name, "glob")
        glob.name = "foobar"
        self.assertEqual(glob.name, "foobar")

    def test_linkage(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        linkage = glob.linkage
        self.assertIsInstance(glob.linkage, llvm.Linkage)
        glob.linkage = linkage
        self.assertEqual(glob.linkage, linkage)
        for linkage in ("internal", "external"):
            glob.linkage = linkage
            self.assertIsInstance(glob.linkage, llvm.Linkage)
            self.assertEqual(glob.linkage.name, linkage)

    def test_visibility(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        visibility = glob.visibility
        self.assertIsInstance(glob.visibility, llvm.Visibility)
        glob.visibility = visibility
        self.assertEqual(glob.visibility, visibility)
        for visibility in ("hidden", "protected", "default"):
            glob.visibility = visibility
            self.assertIsInstance(glob.visibility, llvm.Visibility)
            self.assertEqual(glob.visibility.name, visibility)

    def test_storage_class(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        storage_class = glob.storage_class
        self.assertIsInstance(glob.storage_class, llvm.StorageClass)
        glob.storage_class = storage_class
        self.assertEqual(glob.storage_class, storage_class)
        for storage_class in ("dllimport", "dllexport", "default"):
            glob.storage_class = storage_class
            self.assertIsInstance(glob.storage_class, llvm.StorageClass)
            self.assertEqual(glob.storage_class.name, storage_class)

    def test_add_function_attribute(self):
        mod = self.module()
        fn = mod.get_function("sum")
        fn.add_function_attribute("nocapture")
        with self.assertRaises(ValueError) as raises:
            fn.add_function_attribute("zext")
        self.assertEqual(str(raises.exception), "no such attribute 'zext'")

    def test_module(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        self.assertIs(glob.module, mod)

    def test_type(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        tp = glob.type
        self.assertIsInstance(tp, llvm.TypeRef)

    def test_type_name(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        tp = glob.type
        self.assertEqual(tp.name, "")
        st = mod.get_global_variable("glob_struct")
        self.assertIsNotNone(re.match(r"struct\.glob_type(\.[\d]+)?",
                                      st.global_value_type.name))

    def test_type_printing_variable(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        tp = glob.global_value_type
        self.assertEqual(str(tp), 'i32')

    def test_type_printing_function(self):
        mod = self.module()
        fn = mod.get_function("sum")
        self.assertEqual(str(fn.global_value_type), "i32 (i32, i32)")

    def test_type_printing_struct(self):
        mod = self.module()
        st = mod.get_global_variable("glob_struct")
        self.assertTrue(st.type.is_pointer)
        self.assertIsNotNone(re.match(r'ptr', str(st.type)))
        self.assertIsNotNone(re.match(
            r"%struct\.glob_type(\.[\d]+)? = type { i64, \[2 x i64\] }",
            str(st.global_value_type)))

    def test_close(self):
        glob = self.glob()
        glob.close()
        glob.close()

    def test_is_declaration(self):
        defined = self.module().get_function('sum')
        declared = self.module(asm_sum_declare).get_function('sum')
        self.assertFalse(defined.is_declaration)
        self.assertTrue(declared.is_declaration)

    def test_module_global_variables(self):
        mod = self.module(asm_sum)
        gvars = list(mod.global_variables)
        self.assertEqual(len(gvars), 4)
        for v in gvars:
            self.assertTrue(v.is_global)

    def test_module_functions(self):
        mod = self.module()
        funcs = list(mod.functions)
        self.assertEqual(len(funcs), 1)
        func = funcs[0]
        self.assertTrue(func.is_function)
        self.assertEqual(func.name, 'sum')

        with self.assertRaises(ValueError):
            func.instructions
        with self.assertRaises(ValueError):
            func.operands
        with self.assertRaises(ValueError):
            func.opcode

    def test_function_arguments(self):
        mod = self.module()
        func = mod.get_function('sum')
        self.assertTrue(func.is_function)
        args = list(func.arguments)
        self.assertEqual(len(args), 2)
        self.assertTrue(args[0].is_argument)
        self.assertTrue(args[1].is_argument)
        self.assertEqual(args[0].name, '.1')
        self.assertEqual(str(args[0].type), 'i32')
        self.assertEqual(args[1].name, '.2')
        self.assertEqual(str(args[1].type), 'i32')

        with self.assertRaises(ValueError):
            args[0].blocks
        with self.assertRaises(ValueError):
            args[0].arguments

    def test_function_blocks(self):
        func = self.module().get_function('sum')
        blocks = list(func.blocks)
        self.assertEqual(len(blocks), 1)
        block = blocks[0]
        self.assertTrue(block.is_block)

    def test_block_instructions(self):
        func = self.module().get_function('sum')
        insts = list(list(func.blocks)[0].instructions)
        self.assertEqual(len(insts), 3)
        self.assertTrue(insts[0].is_instruction)
        self.assertTrue(insts[1].is_instruction)
        self.assertTrue(insts[2].is_instruction)
        self.assertEqual(insts[0].opcode, 'add')
        self.assertEqual(insts[1].opcode, 'add')
        self.assertEqual(insts[2].opcode, 'ret')

    def test_instruction_operands(self):
        func = self.module().get_function('sum')
        add = list(list(func.blocks)[0].instructions)[0]
        self.assertEqual(add.opcode, 'add')
        operands = list(add.operands)
        self.assertEqual(len(operands), 2)
        self.assertTrue(operands[0].is_operand)
        self.assertTrue(operands[1].is_operand)
        self.assertEqual(operands[0].name, '.1')
        self.assertEqual(str(operands[0].type), 'i32')
        self.assertEqual(operands[1].name, '.2')
        self.assertEqual(str(operands[1].type), 'i32')

    def test_function_attributes(self):
        ver = llvm.llvm_version_info[0]
        readonly_attrs = [b'memory(read)' if ver > 15 else b'readonly']
        mod = self.module(asm_attributes)
        for func in mod.functions:
            attrs = list(func.attributes)
            if func.name == 'a_readonly_func':
                self.assertEqual(attrs, readonly_attrs)
            elif func.name == 'a_arg0_return_func':
                self.assertEqual(attrs, [])
                args = list(func.arguments)
                self.assertEqual(list(args[0].attributes), [b'returned'])
                self.assertEqual(list(args[1].attributes), [])

    def test_value_kind(self):
        mod = self.module()
        self.assertEqual(mod.get_global_variable('glob').value_kind,
                         llvm.ValueKind.global_variable)
        func = mod.get_function('sum')
        self.assertEqual(func.value_kind, llvm.ValueKind.function)
        block = list(func.blocks)[0]
        self.assertEqual(block.value_kind, llvm.ValueKind.basic_block)
        inst = list(block.instructions)[1]
        self.assertEqual(inst.value_kind, llvm.ValueKind.instruction)
        self.assertEqual(list(inst.operands)[0].value_kind,
                         llvm.ValueKind.constant_int)
        self.assertEqual(list(inst.operands)[1].value_kind,
                         llvm.ValueKind.instruction)

        iasm_func = self.module(asm_inlineasm).get_function('foo')
        iasm_inst = list(list(iasm_func.blocks)[0].instructions)[0]
        self.assertEqual(list(iasm_inst.operands)[0].value_kind,
                         llvm.ValueKind.inline_asm)

    def test_is_constant(self):
        mod = self.module()
        self.assertTrue(mod.get_global_variable('glob').is_constant)
        constant_operands = 0
        for func in mod.functions:
            self.assertTrue(func.is_constant)
            for block in func.blocks:
                self.assertFalse(block.is_constant)
                for inst in block.instructions:
                    self.assertFalse(inst.is_constant)
                    for op in inst.operands:
                        if op.is_constant:
                            constant_operands += 1

        self.assertEqual(constant_operands, 1)

    def test_constant_int(self):
        mod = self.module()
        func = mod.get_function('sum')
        insts = list(list(func.blocks)[0].instructions)
        self.assertEqual(insts[1].opcode, 'add')
        operands = list(insts[1].operands)
        self.assertTrue(operands[0].is_constant)
        self.assertFalse(operands[1].is_constant)
        self.assertEqual(operands[0].get_constant_value(), 0)
        with self.assertRaises(ValueError):
            operands[1].get_constant_value()

        mod = self.module(asm_sum3)
        func = mod.get_function('sum')
        insts = list(list(func.blocks)[0].instructions)
        posint64 = list(insts[1].operands)[0]
        negint64 = list(insts[2].operands)[0]
        self.assertEqual(posint64.get_constant_value(), 5)
        self.assertEqual(negint64.get_constant_value(signed_int=True), -5)

        # Convert from unsigned arbitrary-precision integer to signed i64
        as_u64 = negint64.get_constant_value(signed_int=False)
        as_i64 = int.from_bytes(as_u64.to_bytes(8, 'little'), 'little',
                                signed=True)
        self.assertEqual(as_i64, -5)

    def test_constant_fp(self):
        mod = self.module(asm_double_locale)
        func = mod.get_function('foo')
        insts = list(list(func.blocks)[0].instructions)
        self.assertEqual(len(insts), 2)
        self.assertEqual(insts[0].opcode, 'fadd')
        operands = list(insts[0].operands)
        self.assertTrue(operands[0].is_constant)
        self.assertAlmostEqual(operands[0].get_constant_value(), 0.0)
        self.assertTrue(operands[1].is_constant)
        self.assertAlmostEqual(operands[1].get_constant_value(), 3.14)

        mod = self.module(asm_double_inaccurate)
        func = mod.get_function('foo')
        inst = list(list(func.blocks)[0].instructions)[0]
        operands = list(inst.operands)
        with self.assertRaises(ValueError):
            operands[0].get_constant_value()
        self.assertAlmostEqual(operands[1].get_constant_value(round_fp=True), 0)

    def test_constant_as_string(self):
        mod = self.module(asm_null_constant)
        func = mod.get_function('bar')
        inst = list(list(func.blocks)[0].instructions)[0]
        arg = list(inst.operands)[0]
        self.assertTrue(arg.is_constant)
        self.assertEqual(arg.get_constant_value(), 'ptr null')

    def test_incoming_phi_blocks(self):
        mod = self.module(asm_phi_blocks)
        func = mod.get_function('foo')
        blocks = list(func.blocks)
        instructions = list(blocks[-1].instructions)
        self.assertTrue(instructions[0].is_instruction)
        self.assertEqual(instructions[0].opcode, 'phi')

        incoming_blocks = list(instructions[0].incoming_blocks)
        self.assertEqual(len(incoming_blocks), 2)
        self.assertTrue(incoming_blocks[0].is_block)
        self.assertTrue(incoming_blocks[1].is_block)
        # Test reference to blocks (named or unnamed)
        self.assertEqual(incoming_blocks[0], blocks[-1])
        self.assertEqual(incoming_blocks[1], blocks[0])

        # Test case that should fail
        self.assertNotEqual(instructions[1].opcode, 'phi')
        with self.assertRaises(ValueError):
            instructions[1].incoming_blocks


class TestTypeRef(BaseTest):

    def test_str(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        self.assertEqual(str(glob.global_value_type), "i32")
        glob_struct_type = mod.get_struct_type("struct.glob_type")
        self.assertEqual(str(glob_struct_type),
                         "%struct.glob_type = type { i64, [2 x i64] }")

        elements = list(glob_struct_type.elements)
        self.assertEqual(len(elements), 2)
        self.assertEqual(str(elements[0]), "i64")
        self.assertEqual(str(elements[1]), "[2 x i64]")

    def test_type_kind(self):
        mod = self.module()
        glob = mod.get_global_variable("glob")
        self.assertEqual(glob.type.type_kind, llvm.TypeKind.pointer)
        self.assertTrue(glob.type.is_pointer)

        glob_struct = mod.get_global_variable("glob_struct")
        self.assertEqual(glob_struct.type.type_kind, llvm.TypeKind.pointer)
        self.assertTrue(glob_struct.type.is_pointer)

        stype = glob_struct.global_value_type
        self.assertEqual(stype.type_kind, llvm.TypeKind.struct)
        self.assertTrue(stype.is_struct)

        stype_a, stype_b = stype.elements
        self.assertEqual(stype_a.type_kind, llvm.TypeKind.integer)
        self.assertEqual(stype_b.type_kind, llvm.TypeKind.array)
        self.assertTrue(stype_b.is_array)

        glob_vec_struct_type = mod.get_struct_type("struct.glob_type_vec")
        _, vector_type = glob_vec_struct_type.elements
        self.assertEqual(vector_type.type_kind, llvm.TypeKind.vector)
        self.assertTrue(vector_type.is_vector)

        funcptr = mod.get_function("sum").type
        self.assertEqual(funcptr.type_kind, llvm.TypeKind.pointer)
        functype = mod.get_function("sum").global_value_type
        self.assertEqual(functype.type_kind, llvm.TypeKind.function)

    def test_element_count(self):
        mod = self.module()
        glob_struct_type = mod.get_struct_type("struct.glob_type")
        _, array_type = glob_struct_type.elements
        self.assertEqual(array_type.element_count, 2)
        with self.assertRaises(ValueError):
            glob_struct_type.element_count

    def test_type_width(self):
        mod = self.module()
        glob_struct_type = mod.get_struct_type("struct.glob_type")
        glob_vec_struct_type = mod.get_struct_type("struct.glob_type_vec")
        integer_type, array_type = glob_struct_type.elements
        _, vector_type = glob_vec_struct_type.elements
        self.assertEqual(integer_type.type_width, 64)
        self.assertEqual(vector_type.type_width, 64 * 2)

        # Structs and arrays are not primitive types
        self.assertEqual(glob_struct_type.type_width, 0)
        self.assertEqual(array_type.type_width, 0)

    def test_vararg_function(self):
        # Variadic function
        mod = self.module(asm_vararg_declare)
        func = mod.get_function('vararg')
        decltype = func.global_value_type
        self.assertTrue(decltype.is_function_vararg)

        mod = self.module(asm_sum_declare)
        func = mod.get_function('sum')
        decltype = func.global_value_type
        self.assertFalse(decltype.is_function_vararg)

        # test that the function pointer type cannot use is_function_vararg
        self.assertTrue(func.type.is_pointer)
        with self.assertRaises(ValueError) as raises:
            func.type.is_function_vararg
        self.assertIn("Type ptr is not a function", str(raises.exception))

    def test_function_typeref_as_ir(self):
        mod = self.module()

        [fn] = list(mod.functions)
        # .type gives a pointer type, a problem if it's opaque (llvm15+)
        self.assertEqual(fn.type.type_kind, llvm.TypeKind.pointer)
        self.assertFalse(fn.type.is_function)
        # Use .global_value_type instead
        fnty = fn.global_value_type
        self.assertEqual(fnty.type_kind, llvm.TypeKind.function)
        self.assertTrue(fnty.is_function)
        # Run .as_ir() to get llvmlite.ir.FunctionType
        tyir = fnty.as_ir(ir.global_context)
        self.assertIsInstance(tyir, ir.FunctionType)
        self.assertEqual(tyir.args, (ir.IntType(32), ir.IntType(32)))
        self.assertEqual(tyir.return_type ,ir.IntType(32))

    def test_void_typeref_as_ir(self):
        # Void type can only be used as return-type of llvmlite.ir.FunctionType.
        fnty = ir.FunctionType(ir.VoidType(), ())
        irmod = ir.Module()
        fn = ir.Function(irmod, fnty, "foo")
        mod = self.module(str(irmod))
        fn = mod.get_function("foo")
        gvty = fn.global_value_type
        self.assertEqual(fnty.return_type,
                         gvty.as_ir(ir.global_context).return_type)

    def test_global_typeref_as_ir(self):
        from llvmlite.binding.typeref import _TypeKindToIRType
        ctx = ir.Context()

        skipped = {
            "function",     # tested in test_function_typeref_as_ir
            "void",         # tested in test_void_typeref_as_ir
        }

        makers = {}

        def maker_half():
            yield ir.HalfType()

        makers['half'] = maker_half

        def maker_float():
            yield ir.FloatType()

        makers['float'] = maker_float

        def maker_double():
            yield ir.DoubleType()

        makers['double'] = maker_double

        def maker_integer():
            yield ir.IntType(32)

        makers['integer'] = maker_integer

        def maker_pointer():
            yield ir.PointerType(ir.IntType(8))
            # opaque struct ptr
            yield ctx.get_identified_type("myclass").as_pointer()
            # named struct with defined body
            myclass2 = ctx.get_identified_type("myclass2")
            myclass2.set_body(ir.IntType(8))
            yield myclass2.as_pointer()

        makers['pointer'] = maker_pointer

        def maker_array():
            yield ir.ArrayType(ir.IntType(8), 123)

        makers['array'] = maker_array

        def maker_vector():
            yield ir.VectorType(ir.FloatType(), 2)

        makers['vector'] = maker_vector

        def maker_struct():
            yield ir.LiteralStructType([ir.FloatType(), ir.IntType(64)])
            yield ir.LiteralStructType([ir.FloatType(), ir.IntType(64)],
                                       packed=True)

        makers['struct'] = maker_struct

        # Ensure that number of supported TypeKind matches number of makers
        self.assertEqual({x.name for x in _TypeKindToIRType.keys()},
                         set(makers.keys()) | set(skipped))

        # Test each type-kind
        for type_kind, irtype in _TypeKindToIRType.items():
            if type_kind.name in skipped:
                continue
            for ty in makers[type_kind.name]():
                with self.subTest(f"{type_kind!s} -> {ty}"):
                    irmod = ir.Module(context=ctx)
                    ir.GlobalVariable(irmod, ty, name='gv')
                    asm = str(irmod)
                    mod = llvm.parse_assembly(asm)
                    gv = mod.get_global_variable("gv")
                    gvty = gv.global_value_type
                    got = gvty.as_ir(ir.Context())  # fresh context
                    self.assertEqual(got, ty)
                    self.assertIsInstance(got, irtype)

    def _check_typeref_as_ir_for_wrappers(self, asm, target_symbol):
        # Get a clang++ defined function from a llvm ir
        mod = llvm.parse_assembly(asm)
        cppfn = mod.get_function(target_symbol)
        cppfntype = cppfn.global_value_type

        # Get the function type into a new context
        my_context = ir.Context()  # don't populate global context
        ty = cppfntype.as_ir(ir_ctx=my_context)

        # Build a wrapper module for the cpp function
        wrapper_mod = ir.Module(context=my_context)
        # declare the original function
        declfn = ir.Function(wrapper_mod, ty, name=cppfn.name)
        # populate the wrapper function
        wrapfn = ir.Function(wrapper_mod, ty, name="wrapper")
        builder = ir.IRBuilder(wrapfn.append_basic_block())
        # just call the original function
        builder.call(declfn, wrapfn.args)
        builder.ret_void()
        # Create a new LLVM module with the wrapper
        new_mod = llvm.parse_assembly(str(wrapper_mod))
        self.assertTrue(new_mod.get_function(declfn.name).is_declaration,
                        msg="declfn must not have a body")
        # Merge/link the original module into the new module
        new_mod.link_in(mod, preserve=True)
        self.assertEqual(len(list(new_mod.functions)),
                         len(list(mod.functions)) + 1,
                         msg="the only new function is the wrapper")
        self.assertFalse(new_mod.get_function(declfn.name).is_declaration,
                         msg="declfn must have a body now")
        self.assertEqual(new_mod.get_function(declfn.name).global_value_type,
                         new_mod.get_function(wrapfn.name).global_value_type,
                         msg="declfn and wrapfn must have the same llvm Type")

    def test_typeref_as_ir_for_wrappers_of_cpp_class(self):
        """Exercise extracting C++ defined class types.
        Contains both opaque and non-opaque class definitions.
        """
        self._check_typeref_as_ir_for_wrappers(
            asm_cpp_class,
            "_Z3fooP7MyClass14MyClassDefined",
        )

    def test_typeref_as_ir_for_wrappers_of_cpp_vector_struct(self):
        """Exercise extracting C++ struct types that are passed as vectors.

        IA64 ABI on x86_64 will put struct with two floats as
        a vector of two floats.
        """
        self._check_typeref_as_ir_for_wrappers(
            asm_cpp_vector,
            "_Z3foo8Vector2DPS_",
        )


class TestTarget(BaseTest):

    def test_from_triple(self):
        f = llvm.Target.from_triple
        with self.assertRaises(RuntimeError) as cm:
            f("foobar")
        self.assertIn("No available targets are compatible with",
                      str(cm.exception))
        triple = llvm.get_default_triple()
        target = f(triple)
        self.assertEqual(target.triple, triple)
        target.close()

    def test_create_target_machine(self):
        target = llvm.Target.from_triple(llvm.get_default_triple())
        # With the default settings
        target.create_target_machine('', '', 1, 'default', 'default')
        # With the host's CPU
        cpu = llvm.get_host_cpu_name()
        target.create_target_machine(cpu, '', 1, 'default', 'default')

    def test_name(self):
        t = llvm.Target.from_triple(llvm.get_default_triple())
        u = llvm.Target.from_default_triple()
        self.assertIsInstance(t.name, str)
        self.assertEqual(t.name, u.name)

    def test_description(self):
        t = llvm.Target.from_triple(llvm.get_default_triple())
        u = llvm.Target.from_default_triple()
        self.assertIsInstance(t.description, str)
        self.assertEqual(t.description, u.description)

    def test_str(self):
        target = llvm.Target.from_triple(llvm.get_default_triple())
        s = str(target)
        self.assertIn(target.name, s)
        self.assertIn(target.description, s)

    def test_get_parts_from_triple(self):
        # Tests adapted from llvm-14::llvm/unittests/ADT/TripleTest.cpp
        cases = [
            ("x86_64-scei-ps4",
             llvm.targets.Triple(Arch="x86_64", SubArch='',
                                 Vendor="scei", OS="ps4",
                                 Env="unknown", ObjectFormat="ELF")),
            ("x86_64-sie-ps4",
             llvm.targets.Triple(Arch="x86_64", SubArch='',
                                 Vendor="scei", OS="ps4",
                                 Env="unknown", ObjectFormat="ELF")),
            ("powerpc-dunno-notsure",
             llvm.targets.Triple(Arch="powerpc", SubArch='',
                                 Vendor="unknown", OS="unknown",
                                 Env="unknown", ObjectFormat="ELF")),
            ("powerpcspe-unknown-freebsd",
             llvm.targets.Triple(Arch="powerpc", SubArch='spe',
                                 Vendor="unknown", OS="freebsd",
                                 Env="unknown", ObjectFormat="ELF")),
            ("armv6hl-none-linux-gnueabi",
             llvm.targets.Triple(Arch="arm", SubArch='v6hl',
                                 Vendor="unknown", OS="linux",
                                 Env="gnueabi", ObjectFormat="ELF")),
            ("i686-unknown-linux-gnu",
             llvm.targets.Triple(Arch="i386", SubArch='',
                                 Vendor="unknown", OS="linux",
                                 Env="gnu", ObjectFormat="ELF")),
            ("i686-apple-macosx",
             llvm.targets.Triple(Arch="i386", SubArch='',
                                 Vendor="apple", OS="macosx",
                                 Env="unknown", ObjectFormat="MachO")),
            ("i686-dunno-win32",
             llvm.targets.Triple(Arch="i386", SubArch='',
                                 Vendor="unknown", OS="windows",
                                 Env="msvc", ObjectFormat="COFF")),
            ("s390x-ibm-zos",
             llvm.targets.Triple(Arch="s390x", SubArch='',
                                 Vendor="ibm", OS="zos",
                                 Env="unknown", ObjectFormat="GOFF")),
            ("wasm64-wasi",
             llvm.targets.Triple(Arch="wasm64", SubArch='',
                                 Vendor="unknown", OS="wasi",
                                 Env="unknown", ObjectFormat="Wasm")),
        ]

        for case in cases:
            triple_str, triple_obj = case
            res = llvm.get_triple_parts(triple_str)

            self.assertEqual(res, triple_obj)


class TestTargetData(BaseTest):

    def target_data(self):
        return llvm.create_target_data("e-m:e-i64:64-f80:128-n8:16:32:64-S128")

    def test_get_abi_size(self):
        td = self.target_data()
        glob = self.glob()
        self.assertEqual(td.get_abi_size(glob.type), 8)

    def test_get_pointee_abi_size(self):
        td = self.target_data()

        glob = self.glob()
        self.assertEqual(td.get_abi_size(glob.global_value_type), 4)

        glob = self.glob("glob_struct")
        self.assertEqual(td.get_abi_size(glob.global_value_type), 24)

    def test_get_struct_element_offset(self):
        td = self.target_data()
        glob = self.glob("glob_struct")

        with self.assertRaises(ValueError):
            td.get_element_offset(glob.type, 0)

        struct_type = glob.global_value_type
        self.assertEqual(td.get_element_offset(struct_type, 0), 0)
        self.assertEqual(td.get_element_offset(struct_type, 1), 8)


class TestTargetMachine(BaseTest):

    def test_target_data_from_tm(self):
        tm = self.target_machine(jit=False)
        td = tm.target_data
        mod = self.module()
        gv_i32 = mod.get_global_variable("glob")
        # A global is a pointer, it has the ABI size of a pointer
        pointer_size = 4 if sys.maxsize < 2 ** 32 else 8
        self.assertEqual(td.get_abi_size(gv_i32.type), pointer_size)


class TestDylib(BaseTest):

    def test_bad_library(self):
        with self.assertRaises(RuntimeError):
            llvm.load_library_permanently("zzzasdkf;jasd;l")

    @unittest.skipUnless(platform.system() in ["Linux"],
                         "test only works on Linux")
    def test_libm(self):
        libm = find_library("m")
        llvm.load_library_permanently(libm)


class TestAnalysis(BaseTest):
    def build_ir_module(self):
        m = ir.Module()
        ft = ir.FunctionType(ir.IntType(32), [ir.IntType(32), ir.IntType(32)])
        fn = ir.Function(m, ft, "foo")
        bd = ir.IRBuilder(fn.append_basic_block())
        x, y = fn.args
        z = bd.add(x, y)
        bd.ret(z)
        return m

    def test_get_function_cfg_on_ir(self):
        mod = self.build_ir_module()
        foo = mod.get_global('foo')
        dot_showing_inst = llvm.get_function_cfg(foo)
        dot_without_inst = llvm.get_function_cfg(foo, show_inst=False)
        inst = "%.5 = add i32 %.1, %.2"
        self.assertIn(inst, dot_showing_inst)
        self.assertNotIn(inst, dot_without_inst)

    def test_function_cfg_on_llvm_value(self):
        defined = self.module().get_function('sum')
        dot_showing_inst = llvm.get_function_cfg(defined, show_inst=True)
        dot_without_inst = llvm.get_function_cfg(defined, show_inst=False)
        # Check "digraph"
        prefix = 'digraph'
        self.assertIn(prefix, dot_showing_inst)
        self.assertIn(prefix, dot_without_inst)
        # Check function name
        fname = "CFG for 'sum' function"
        self.assertIn(fname, dot_showing_inst)
        self.assertIn(fname, dot_without_inst)
        # Check instruction
        inst = "%.3 = add i32 %.1, %.2"
        self.assertIn(inst, dot_showing_inst)
        self.assertNotIn(inst, dot_without_inst)


class TestTypeParsing(BaseTest):
    @contextmanager
    def check_parsing(self):
        mod = ir.Module()
        # Yield to caller and provide the module for adding
        # new GV.
        yield mod
        # Caller yield back and continue with testing
        asm = str(mod)
        llvm.parse_assembly(asm)

    def test_literal_struct(self):
        # Natural layout
        with self.check_parsing() as mod:
            typ = ir.LiteralStructType([ir.IntType(32)])
            gv = ir.GlobalVariable(mod, typ, "foo")
            # Also test constant text repr
            gv.initializer = ir.Constant(typ, [1])

        # Packed layout
        with self.check_parsing() as mod:
            typ = ir.LiteralStructType([ir.IntType(32)],
                                       packed=True)
            gv = ir.GlobalVariable(mod, typ, "foo")
            # Also test constant text repr
            gv.initializer = ir.Constant(typ, [1])

        # Packed layout created from Constant.literal_struct
        with self.check_parsing() as mod:
            const = ir.Constant.literal_struct([ir.IntType(32)(1),
                                                ir.IntType(32)(2)],
                                               packed=True)
            gv = ir.GlobalVariable(mod, const.type, "foo")
            gv.initializer = const


class TestGlobalConstructors(TestMCJit):
    @unittest.skipIf(platform.system() == "Darwin",
                     "__cxa_atexit is broken on OSX in MCJIT")
    def test_global_ctors_dtors(self):
        # test issue #303
        # (https://github.com/numba/llvmlite/issues/303)
        mod = self.module(asm_global_ctors)
        ee = self.jit(mod)
        ee.finalize_object()

        ee.run_static_constructors()

        # global variable should have been initialized
        ptr_addr = ee.get_global_value_address("A")
        ptr_t = ctypes.POINTER(ctypes.c_int32)
        ptr = ctypes.cast(ptr_addr, ptr_t)
        self.assertEqual(ptr.contents.value, 10)

        foo_addr = ee.get_function_address("foo")
        foo = ctypes.CFUNCTYPE(ctypes.c_int32)(foo_addr)
        self.assertEqual(foo(), 12)

        ee.run_static_destructors()

        # destructor should have run
        self.assertEqual(ptr.contents.value, 20)


class TestGlobalVariables(BaseTest):
    def check_global_variable_linkage(self, linkage, has_undef=True):
        # This test default initializer on global variables with different
        # linkages.  Some linkages requires an initializer be present, while
        # it is optional for others.  This test uses ``parse_assembly()``
        # to verify that we are adding an `undef` automatically if user didn't
        # specific one for certain linkages.  It is a IR syntax error if the
        # initializer is not present for certain linkages e.g. "external".
        mod = ir.Module()
        typ = ir.IntType(32)
        gv = ir.GlobalVariable(mod, typ, "foo")
        gv.linkage = linkage
        asm = str(mod)
        # check if 'undef' is present
        if has_undef:
            self.assertIn('undef', asm)
        else:
            self.assertNotIn('undef', asm)
        # parse assembly to ensure correctness
        self.module(asm)

    def test_internal_linkage(self):
        self.check_global_variable_linkage('internal')

    def test_common_linkage(self):
        self.check_global_variable_linkage('common')

    def test_external_linkage(self):
        self.check_global_variable_linkage('external', has_undef=False)

    def test_available_externally_linkage(self):
        self.check_global_variable_linkage('available_externally')

    def test_private_linkage(self):
        self.check_global_variable_linkage('private')

    def test_linkonce_linkage(self):
        self.check_global_variable_linkage('linkonce')

    def test_weak_linkage(self):
        self.check_global_variable_linkage('weak')

    def test_appending_linkage(self):
        self.check_global_variable_linkage('appending')

    def test_extern_weak_linkage(self):
        self.check_global_variable_linkage('extern_weak', has_undef=False)

    def test_linkonce_odr_linkage(self):
        self.check_global_variable_linkage('linkonce_odr')

    def test_weak_odr_linkage(self):
        self.check_global_variable_linkage('weak_odr')


@unittest.skipUnless(platform.machine().startswith('x86'), "only on x86")
class TestInlineAsm(BaseTest):
    def test_inlineasm(self):
        llvm.initialize_native_asmparser()
        m = self.module(asm=asm_inlineasm)
        tm = self.target_machine(jit=False)
        asm = tm.emit_assembly(m)
        self.assertIn('nop', asm)


class TestObjectFile(BaseTest):

    mod_asm = """
        ;ModuleID = <string>
        target triple = "{triple}"

        declare i32 @sum(i32 %.1, i32 %.2)

        define i32 @sum_twice(i32 %.1, i32 %.2) {{
            %.3 = call i32 @sum(i32 %.1, i32 %.2)
            %.4 = call i32 @sum(i32 %.3, i32 %.3)
            ret i32 %.4
        }}
    """

    def test_object_file(self):
        target_machine = self.target_machine(jit=False)
        mod = self.module()
        obj_bin = target_machine.emit_object(mod)
        obj = llvm.ObjectFileRef.from_data(obj_bin)
        # Check that we have a text section, and that she has a name and data
        has_text_and_data = False
        last_address = -1
        for s in obj.sections():
            if (
                s.is_text()
                and len(s.data()) > 0
                and s.address() is not None
                and last_address < s.address()
            ):
                has_text_and_data = True
                last_address = s.address()
                break
        self.assertTrue(has_text_and_data)

    def test_add_object_file(self):
        target_machine = self.target_machine(jit=False)
        mod = self.module()
        obj_bin = target_machine.emit_object(mod)
        obj = llvm.ObjectFileRef.from_data(obj_bin)

        jit = llvm.create_mcjit_compiler(self.module(self.mod_asm),
                                         target_machine)

        jit.add_object_file(obj)

        sum_twice = CFUNCTYPE(c_int, c_int, c_int)(
            jit.get_function_address("sum_twice"))

        self.assertEqual(sum_twice(2, 3), 10)

    def test_add_object_file_from_filesystem(self):
        target_machine = self.target_machine(jit=False)
        mod = self.module()
        obj_bin = target_machine.emit_object(mod)
        temp_desc, temp_path = mkstemp()

        try:
            try:
                f = os.fdopen(temp_desc, "wb")
                f.write(obj_bin)
                f.flush()
            finally:
                f.close()

            jit = llvm.create_mcjit_compiler(self.module(self.mod_asm),
                                             target_machine)

            jit.add_object_file(temp_path)
        finally:
            os.unlink(temp_path)

        sum_twice = CFUNCTYPE(c_int, c_int, c_int)(
            jit.get_function_address("sum_twice"))

        self.assertEqual(sum_twice(2, 3), 10)

    def test_get_section_content(self):
        # See Issue #632 - section contents were getting truncated at null
        # bytes.
        elf = bytes.fromhex(issue_632_elf)
        obj = llvm.ObjectFileRef.from_data(elf)
        for s in obj.sections():
            if s.is_text():
                self.assertEqual(len(s.data()), 31)
                self.assertEqual(s.data().hex(), issue_632_text)


class TestLLVMLockCallbacks(BaseTest):
    def test_lock_callbacks(self):
        events = []

        def acq():
            events.append('acq')

        def rel():
            events.append('rel')

        # register callback
        llvm.ffi.register_lock_callback(acq, rel)

        # Check: events are initially empty
        self.assertFalse(events)
        # Call LLVM functions
        llvm.create_new_module_pass_manager()
        # Check: there must be at least one acq and one rel
        self.assertIn("acq", events)
        self.assertIn("rel", events)

        # unregister callback
        llvm.ffi.unregister_lock_callback(acq, rel)

        # Check: removing non-existent callbacks will trigger a ValueError
        with self.assertRaises(ValueError):
            llvm.ffi.unregister_lock_callback(acq, rel)


class TestPipelineTuningOptions(BaseTest):

    def pto(self):
        return llvm.create_pipeline_tuning_options()

    def test_close(self):
        pto = self.pto()
        pto.close()

    def test_speed_level(self):
        pto = self.pto()
        self.assertIsInstance(pto.speed_level, int)
        for i in range(4):
            pto.speed_level = i
            self.assertEqual(pto.speed_level, i)

    def test_size_level(self):
        pto = self.pto()
        self.assertIsInstance(pto.size_level, int)
        for i in range(3):
            pto.size_level = i
            self.assertEqual(pto.size_level, i)

    def test_inlining_threshold(self):
        pto = self.pto()
        self.assertIsInstance(pto.inlining_threshold, int)
        for i in (25, 80, 350):
            pto.inlining_threshold = i

    def test_loop_interleaving(self):
        pto = self.pto()
        self.assertIsInstance(pto.loop_interleaving, bool)
        for b in (True, False):
            pto.loop_interleaving = b
            self.assertEqual(pto.loop_interleaving, b)

    def test_loop_vectorization(self):
        pto = self.pto()
        self.assertIsInstance(pto.loop_vectorization, bool)
        for b in (True, False):
            pto.loop_vectorization = b
            self.assertEqual(pto.loop_vectorization, b)

    def test_slp_vectorization(self):
        pto = self.pto()
        self.assertIsInstance(pto.slp_vectorization, bool)
        for b in (True, False):
            pto.slp_vectorization = b
            self.assertEqual(pto.slp_vectorization, b)

    def test_loop_unrolling(self):
        pto = self.pto()
        self.assertIsInstance(pto.loop_unrolling, bool)
        for b in (True, False):
            pto.loop_unrolling = b
            self.assertEqual(pto.loop_unrolling, b)

    def test_speed_level_constraints(self):
        pto = self.pto()
        with self.assertRaises(ValueError):
            pto.speed_level = 4
        with self.assertRaises(ValueError):
            pto.speed_level = -1

    def test_size_level_constraints(self):
        pto = self.pto()
        with self.assertRaises(ValueError):
            pto.size_level = 3
        with self.assertRaises(ValueError):
            pto.speed_level = -1
        with self.assertRaises(ValueError):
            pto.speed_level = 3
            pto.size_level = 2


class NewPassManagerMixin(object):

    def pb(self, speed_level=0, size_level=0):
        tm = self.target_machine(jit=False)
        pto = llvm.create_pipeline_tuning_options(speed_level, size_level)
        pb = llvm.create_pass_builder(tm, pto)
        return pb


class TestPassBuilder(BaseTest, NewPassManagerMixin):

    def test_close(self):
        pb = self.pb()
        pb.close()

    def test_pto(self):
        tm = self.target_machine(jit=False)
        pto = llvm.create_pipeline_tuning_options(3, 0)
        pto.inlining_threshold = 2
        pto.loop_interleaving = True
        pto.loop_vectorization = True
        pto.slp_vectorization = True
        pto.loop_unrolling = False
        pb = llvm.create_pass_builder(tm, pto)
        pb.close()

    def test_get_module_pass_manager(self):
        pb = self.pb()
        mpm = pb.getModulePassManager()
        mpm.run(self.module(), pb)
        pb.close()

    def test_get_function_pass_manager(self):
        pb = self.pb()
        fpm = pb.getFunctionPassManager()
        fpm.run(self.module().get_function("sum"), pb)
        pb.close()

    def test_time_passes(self):
        """Test pass timing reports for O3 and O0 optimization levels"""
        def run_with_timing(speed_level):
            mod = self.module()
            pb = self.pb(speed_level=speed_level, size_level=0)
            pb.start_pass_timing()
            mpm = pb.getModulePassManager()
            mpm.run(mod, pb)
            report = pb.finish_pass_timing()
            pb.close()
            return report

        report_O3 = run_with_timing(3)
        report_O0 = run_with_timing(0)

        self.assertIsInstance(report_O3, str)
        self.assertIsInstance(report_O0, str)
        self.assertEqual(report_O3.count("Pass execution timing report"), 1)
        self.assertEqual(report_O0.count("Pass execution timing report"), 1)

    def test_empty_report(self):
        mod = self.module()
        pb = self.pb()
        mpm = pb.getModulePassManager()
        mpm.run(mod, pb)
        pb.start_pass_timing()
        report = pb.finish_pass_timing()
        pb.close()
        self.assertFalse(report)

    def test_multiple_timers_error(self):
        mod = self.module()
        pb = self.pb()
        pb.start_pass_timing()
        mpm = pb.getModulePassManager()
        mpm.run(mod, pb)
        pb.finish_pass_timing()
        with self.assertRaisesRegex(RuntimeError, "only be done once"):
            pb.start_pass_timing()
        pb.close()

    def test_empty_report_error(self):
        mod = self.module()
        pb = self.pb()
        mpm = pb.getModulePassManager()
        mpm.run(mod, pb)
        with self.assertRaisesRegex(RuntimeError, "not enabled"):
            pb.finish_pass_timing()
        pb.close()


class TestNewModulePassManager(BaseTest, NewPassManagerMixin):
    def pm(self):
        return llvm.create_new_module_pass_manager()

    def run_o_n(self, level):
        mod = self.module()
        orig_asm = str(mod)
        pb = self.pb(speed_level=level, size_level=0)
        mpm = pb.getModulePassManager()
        mpm.run(mod, pb)
        optimized_asm = str(mod)
        return orig_asm, optimized_asm

    def test_close(self):
        mpm = self.pm()
        mpm.close()

    def test_run_o3(self):
        orig_asm, optimized_asm = self.run_o_n(3)
        self.assertIn("%.4", orig_asm)
        self.assertNotIn("%.4", optimized_asm)

    def test_run_o0(self):
        orig_asm, optimized_asm = self.run_o_n(0)
        self.assertIn("%.4", orig_asm)
        self.assertIn("%.4", optimized_asm)

    def test_instcombine(self):
        pb = self.pb()
        mpm = self.pm()
        mpm.add_instruction_combine_pass()
        mod = self.module(asm_sum4)
        orig_asm = str(mod)
        mpm.run(mod, pb)
        optimized_asm = str(mod)
        self.assertIn("%.3", orig_asm)
        self.assertNotIn("%.3", optimized_asm)

    def test_optnone(self):
        pb = self.pb(speed_level=3, size_level=0)
        orig_asm = str(asm_alloca_optnone.replace("optnone ", ""))
        mod = llvm.parse_assembly(orig_asm)
        mpm = pb.getModulePassManager()
        mpm.run(mod, pb)
        optimized_asm = str(mod)
        self.assertIn("alloca", orig_asm)
        self.assertNotIn("alloca", optimized_asm)

        # Module shouldn't be optimized if the function has `optnone` attached
        orig_asm_optnone = str(asm_alloca_optnone)
        mpm = pb.getModulePassManager()
        mod = llvm.parse_assembly(orig_asm_optnone)
        mpm.run(mod, pb)
        optimized_asm_optnone = str(mod)
        self.assertIn("alloca", orig_asm_optnone)
        self.assertIn("alloca", optimized_asm_optnone)

    def test_add_passes(self):
        mpm = self.pm()
        mpm.add_argument_promotion_pass()
        mpm.add_post_order_function_attributes_pass()
        mpm.add_verifier()
        mpm.add_constant_merge_pass()
        mpm.add_dead_arg_elimination_pass()
        mpm.add_dot_call_graph_printer_pass()
        mpm.add_always_inliner_pass()
        mpm.add_rpo_function_attrs_pass()
        mpm.add_global_dead_code_eliminate_pass()
        mpm.add_global_opt_pass()
        mpm.add_ipsccp_pass()
        mpm.add_internalize_pass()
        mpm.add_loop_extract_pass()
        mpm.add_merge_functions_pass()
        mpm.add_partial_inliner_pass()
        mpm.add_strip_symbols_pass()
        mpm.add_strip_dead_debug_info_pass()
        mpm.add_strip_dead_prototype_pass()
        mpm.add_strip_debug_declare_pass()
        mpm.add_strip_non_debug_symbols_pass()
        mpm.add_aa_eval_pass()
        mpm.add_simplify_cfg_pass()
        mpm.add_loop_unroll_pass()
        mpm.add_instruction_combine_pass()
        mpm.add_jump_threading_pass()
        mpm.add_cfg_printer_pass()
        mpm.add_cfg_only_printer_pass()
        mpm.add_dom_printer_pass()
        mpm.add_dom_only_printer_pass()
        mpm.add_post_dom_printer_pass()
        mpm.add_post_dom_only_printer_pass()
        mpm.add_dom_viewer_pass()
        mpm.add_dom_only_printer_pass()
        mpm.add_post_dom_viewer_pass()
        mpm.add_post_dom_only_viewer_pass()
        mpm.add_lint_pass()
        mpm.add_aggressive_dce_pass()
        mpm.add_break_critical_edges_pass()
        mpm.add_dead_store_elimination_pass()
        mpm.add_dead_code_elimination_pass()
        mpm.add_aggressive_instcombine_pass()
        mpm.add_lcssa_pass()
        mpm.add_new_gvn_pass()
        mpm.add_loop_simplify_pass()
        mpm.add_loop_unroll_and_jam_pass()
        mpm.add_sccp_pass()
        mpm.add_lower_atomic_pass()
        mpm.add_lower_invoke_pass()
        mpm.add_lower_switch_pass()
        mpm.add_mem_copy_opt_pass()
        mpm.add_unify_function_exit_nodes_pass()
        mpm.add_reassociate_pass()
        mpm.add_register_to_memory_pass()
        mpm.add_sroa_pass()
        mpm.add_sinking_pass()
        mpm.add_tail_call_elimination_pass()
        mpm.add_instruction_namer_pass()
        mpm.add_loop_deletion_pass()
        mpm.add_loop_strength_reduce_pass()
        mpm.add_loop_rotate_pass()
        mpm.add_refprune_pass()


class TestNewFunctionPassManager(BaseTest, NewPassManagerMixin):
    def pm(self):
        return llvm.create_new_function_pass_manager()

    def test_close(self):
        fpm = self.pm()
        fpm.close()

    def run_o_n(self, level):
        mod = self.module()
        fun = mod.get_function("sum")
        orig_asm = str(fun)
        pb = self.pb(speed_level=level, size_level=0)
        fpm = pb.getFunctionPassManager()
        fpm.run(fun, pb)
        optimized_asm = str(fun)
        return orig_asm, optimized_asm

    def test_run_o3(self):
        orig_asm, optimized_asm = self.run_o_n(3)
        self.assertIn("%.4", orig_asm)
        self.assertNotIn("%.4", optimized_asm)

    def test_run_o0(self):
        orig_asm, optimized_asm = self.run_o_n(0)
        self.assertIn("%.4", orig_asm)
        self.assertIn("%.4", optimized_asm)

    def test_optnone(self):
        pb = self.pb(speed_level=3, size_level=0)
        orig_asm = str(asm_alloca_optnone.replace("optnone ", ""))
        fun = llvm.parse_assembly(orig_asm).get_function("foo")
        fpm = pb.getFunctionPassManager()
        fpm.run(fun, pb)
        optimized_asm = str(fun)
        self.assertIn("alloca", orig_asm)
        self.assertNotIn("alloca", optimized_asm)

        # Function shouldn't be optimized if the function has `optnone` attached
        orig_asm_optnone = str(asm_alloca_optnone)
        fun = llvm.parse_assembly(orig_asm_optnone).get_function("foo")
        fpm = pb.getFunctionPassManager()
        fpm.run(fun, pb)
        optimized_asm_optnone = str(fun)
        self.assertIn("alloca", orig_asm_optnone)
        self.assertIn("alloca", optimized_asm_optnone)

    def test_instcombine(self):
        pb = self.pb()
        fpm = self.pm()
        fun = self.module(asm_sum4).get_function("sum")
        fpm.add_instruction_combine_pass()
        orig_asm = str(fun)
        fpm.run(fun, pb)
        optimized_asm = str(fun)
        self.assertIn("%.3", orig_asm)
        self.assertNotIn("%.3", optimized_asm)

    # This should not crash
    def test_declarations(self):
        pb = self.pb(3)
        fpm = pb.getFunctionPassManager()
        for fun in llvm.parse_assembly(asm_declaration).functions:
            fpm.run(fun, pb)

    def test_add_passes(self):
        fpm = self.pm()
        fpm.add_aa_eval_pass()
        fpm.add_simplify_cfg_pass()
        fpm.add_loop_unroll_pass()
        fpm.add_instruction_combine_pass()
        fpm.add_jump_threading_pass()
        fpm.add_cfg_printer_pass()
        fpm.add_cfg_only_printer_pass()
        fpm.add_dom_printer_pass()
        fpm.add_dom_only_printer_pass()
        fpm.add_post_dom_printer_pass()
        fpm.add_post_dom_only_printer_pass()
        fpm.add_dom_viewer_pass()
        fpm.add_dom_only_printer_pass()
        fpm.add_post_dom_viewer_pass()
        fpm.add_post_dom_only_viewer_pass()
        fpm.add_lint_pass()
        fpm.add_aggressive_dce_pass()
        fpm.add_break_critical_edges_pass()
        fpm.add_dead_store_elimination_pass()
        fpm.add_dead_code_elimination_pass()
        fpm.add_aggressive_instcombine_pass()
        fpm.add_lcssa_pass()
        fpm.add_new_gvn_pass()
        fpm.add_loop_simplify_pass()
        fpm.add_loop_unroll_and_jam_pass()
        fpm.add_sccp_pass()
        fpm.add_lower_atomic_pass()
        fpm.add_lower_invoke_pass()
        fpm.add_lower_switch_pass()
        fpm.add_mem_copy_opt_pass()
        fpm.add_unify_function_exit_nodes_pass()
        fpm.add_reassociate_pass()
        fpm.add_register_to_memory_pass()
        fpm.add_sroa_pass()
        fpm.add_sinking_pass()
        fpm.add_tail_call_elimination_pass()
        fpm.add_instruction_namer_pass()
        fpm.add_loop_deletion_pass()
        fpm.add_loop_strength_reduce_pass()
        fpm.add_loop_rotate_pass()
        fpm.add_refprune_pass()


@unittest.skipUnless(os.environ.get('LLVMLITE_DIST_TEST'),
                     "Distribution-specific test")
@needs_lief
class TestBuild(TestCase):
    # These tests are for use by the Numba project maintainers to check that
    # package builds for which they are responsible are producing artifacts in
    # the expected way. If you are a package maintainer and these tests are
    # running, they shouldn't be by default. The only way they will run is if
    # the environment variable LLVMLITE_DIST_TEST is set. The things they are
    # checking are based on how the Numba project maintainers want to ship the
    # packages, this may be entirely different to how other maintainers wish to
    # ship. Basically, don't enable these tests unless you are sure they are
    # suitable for your use case.
    #
    # The llvmlite DSO is the foundation of Numba's JIT compiler stack and is
    # also used by other similar projects. It has to link against LLVM as that
    # is what provides the tooling to do e.g. IR generation and JIT compilation.
    # There are various options surrounding how to build LLVM and then how to
    # link it into llvmlite. There have been many occurences of surprising
    # linkages, symbol collisions and various other issues.
    #
    # The following tests are designed to try and test out some of the more
    # common combinations of package type and linkage.
    #
    # NOTE: For Numba project maintainers on packaging formats and expected
    # linkage. The following dictionaries capture the state of linkage as of
    # llvmlite release 0.44. This is not an indication that it is correct, just
    # that this is what is present in practice and clearly "works" to a large
    # degree by virtue of having fixed the few reported issues. If you need to
    # modify these dictionaries based on new information, that's fine, just make
    # sure that it is an understood action opposed to just capturing what
    # happened!

    wheel_expected = {"linux": {"x86_64": set(["pthread",
                                               "z",
                                               "dl",
                                               "m",
                                               "gcc_s",
                                               "c",
                                               "rt",
                                               "stdc++",
                                               "ld-linux-x86-64",]),
                                "aarch64":  set(["pthread",
                                                 "z",
                                                 "dl",
                                                 "m",
                                                 "gcc_s",
                                                 "c",
                                                 "rt",
                                                 "stdc++",]),
                                }, # end linux
                      # NOTE: on windows, this includes a "capture what is
                      # present and known to work and make sure it doesn"t
                      # change" approach.
                      "windows": {"amd64": set(["advapi32",
                                                "kernel32",
                                                "ntdll",
                                                "msvcp140",
                                                "vcruntime140",
                                                "vcruntime140_1",
                                                "api-ms-win-crt-convert-l1-1-0",
                                                "api-ms-win-crt-environment-l1-1-0", # noqa: E501
                                                "api-ms-win-crt-heap-l1-1-0",
                                                "api-ms-win-crt-locale-l1-1-0",
                                                "api-ms-win-crt-math-l1-1-0",
                                                "api-ms-win-crt-runtime-l1-1-0",
                                                "api-ms-win-crt-stdio-l1-1-0",
                                                "api-ms-win-crt-string-l1-1-0",
                                                "api-ms-win-crt-time-l1-1-0",
                                                "api-ms-win-crt-utility-l1-1-0",
                                                "shell32",  # this is delayed
                                                "ole32",]), # also delayed
                                  }, # end windows
                      "darwin": {"x86_64": set(["llvmlite",
                                                "system",
                                                "z",
                                                "corefoundation",
                                                "c++",]),
                                 "arm64": set(["llvmlite",
                                               "system",
                                               "z",
                                               "c++",]),
                                 },# end darwin
                      } # end wheel_expected

    conda_expected = {"linux": {"x86_64": set(["pthread",
                                               "z",
                                               "zstd",
                                               "dl",
                                               "m",
                                               "gcc_s",
                                               "c",
                                               # "stdc++", conda has static c++
                                               "ld-linux-x86-64",]),
                                "aarch64":  set(["pthread",
                                                 "z",
                                                 "zstd",
                                                 "dl",
                                                 "m",
                                                 "gcc_s",
                                                 "c",
                                                 # "stdc++", conda has static c++ # noqa: E501
                                                 "ld-linux-aarch64",]),
                                }, # end linux
                      # NOTE: on windows, this includes a "capture what is
                      # present and known to work and make sure it doesn"t
                      # change" approach.
                      "windows": {"amd64": set(["z",
                                                "zstd",
                                                "advapi32",
                                                "kernel32",
                                                "ntdll",
                                                "msvcp140",
                                                "vcruntime140",
                                                "vcruntime140_1",
                                                "api-ms-win-crt-convert-l1-1-0",
                                                "api-ms-win-crt-environment-l1-1-0", # noqa: E501
                                                "api-ms-win-crt-heap-l1-1-0",
                                                "api-ms-win-crt-locale-l1-1-0",
                                                "api-ms-win-crt-math-l1-1-0",
                                                "api-ms-win-crt-runtime-l1-1-0",
                                                "api-ms-win-crt-stdio-l1-1-0",
                                                "api-ms-win-crt-string-l1-1-0",
                                                "api-ms-win-crt-time-l1-1-0",
                                                "api-ms-win-crt-utility-l1-1-0",
                                                "shell32",  # this is delayed
                                                "ole32",]), # also delayed
                                  }, # end windows
                      "darwin": {"x86_64": set(["llvmlite",
                                                "system",
                                                "z",
                                                "zstd",
                                                "corefoundation",
                                                "c++",]),
                                 "arm64": set(["llvmlite",
                                               "system",
                                               "z",
                                               "zstd",
                                               "c++",]),
                                 },# end darwin
                      } # end wheel_expected

    def check_linkage(self, info, package_type):
        machine = platform.machine().lower()
        os_name = platform.system().lower()

        if package_type == "wheel":
            expected = self.wheel_expected[os_name][machine]
        elif package_type == "conda":
            expected = self.conda_expected[os_name][machine]
        else:
            raise ValueError(f"Unexpected package type: {package_type}")

        got = set(info["canonicalised_linked_libraries"])

        try:
            self.assertEqual(expected, got)
        except AssertionError as e:
            msg = ("Unexpected linkage encountered for libllvmlite:\n"
                   f"Expected: {sorted(expected)}\n"
                   f"     Got: {sorted(got)}\n\n"
                   f"Difference: {set.symmetric_difference(expected, got)}\n"
                   f"Only in Expected: {set.difference(expected, got)}\n"
                   f"Only in Got: {set.difference(got, expected)}\n")
            raise AssertionError(msg) from e

    @is_wheel_package
    def test_wheel_build(self):
        info = llvm.config.get_sysinfo()
        self.assertEqual(info['llvm_linkage_type'], "static")
        self.assertEqual(info['llvm_assertions_state'], "on")
        self.check_linkage(info, "wheel")

    @is_conda_package
    def test_conda_build(self):
        info = llvm.config.get_sysinfo()
        self.assertEqual(info['llvm_linkage_type'], "static")
        self.assertEqual(info['llvm_assertions_state'], "on")

        self.check_linkage(info, "conda")
        if platform.system().lower() == "linux":
            self.assertEqual(info['libstdcxx_linkage_type'], "static")


if __name__ == "__main__":
    unittest.main()
