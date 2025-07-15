"""
IR Construction Tests
"""

import copy
import itertools
import pickle
import re
import textwrap
import unittest

from . import TestCase
from llvmlite import ir
from llvmlite import binding as llvm

# FIXME: Remove me once typed pointers are no longer supported.
from llvmlite import opaque_pointers_enabled


int1 = ir.IntType(1)
int8 = ir.IntType(8)
int16 = ir.IntType(16)
int32 = ir.IntType(32)
int64 = ir.IntType(64)
hlf = ir.HalfType()
flt = ir.FloatType()
dbl = ir.DoubleType()


class TestBase(TestCase):
    """
    Utilities for IR tests.
    """

    def assertInText(self, pattern, text):
        """
        Assert *pattern* is in *text*, ignoring any whitespace differences
        (including newlines).
        """

        def escape(c):
            if not c.isalnum() and not c.isspace():
                return '\\' + c
            return c

        pattern = ''.join(map(escape, pattern))
        regex = re.sub(r'\s+', r'\\s*', pattern)
        self.assertRegex(text, regex)

    def assert_ir_line(self, line, mod):
        lines = [line.strip() for line in str(mod).splitlines()]
        self.assertIn(line, lines)

    def assert_valid_ir(self, mod):
        llvm.parse_assembly(str(mod))

    def assert_pickle_correctly(self, irobject):
        """Assert that the IR object pickles and unpickles correctly.
        The IR string is equal and that their type is equal
        """
        newobject = pickle.loads(pickle.dumps(irobject, protocol=-1))
        self.assertIs(irobject.__class__, newobject.__class__)
        self.assertEqual(str(irobject), str(newobject))
        return newobject

    def module(self):
        return ir.Module()

    def function(self, module=None, name='my_func'):
        module = module or self.module()
        fnty = ir.FunctionType(int32, (int32, int32, dbl,
                                       ir.PointerType(int32)))
        return ir.Function(module, fnty, name)

    def block(self, func=None, name=''):
        func = func or self.function()
        return func.append_basic_block(name)

    def descr(self, thing):
        buf = []
        thing.descr(buf)
        return "".join(buf)

    def _normalize_asm(self, asm):
        asm = textwrap.dedent(asm)
        # Normalize indent
        asm = asm.replace("\n    ", "\n  ")
        return asm

    def check_descr_regex(self, descr, asm):
        expected = self._normalize_asm(asm)
        self.assertRegex(descr, expected)

    def check_descr(self, descr, asm):
        expected = self._normalize_asm(asm)
        self.assertEqual(descr, expected)

    def check_block(self, block, asm):
        self.check_descr(self.descr(block), asm)

    def check_block_regex(self, block, asm):
        self.check_descr_regex(self.descr(block), asm)

    def check_module_body(self, module, asm):
        expected = self._normalize_asm(asm)
        actual = module._stringify_body()
        self.assertEqual(actual.strip(), expected.strip())

    def check_metadata(self, module, asm):
        """
        Check module metadata against *asm*.
        """
        expected = self._normalize_asm(asm)
        actual = module._stringify_metadata()
        self.assertEqual(actual.strip(), expected.strip())

    def check_func_body(self, func, asm):
        expected = self._normalize_asm(asm)
        actual = self.descr(func)
        actual = actual.partition('{')[2].rpartition('}')[0]
        self.assertEqual(actual.strip(), expected.strip())


class TestFunction(TestBase):

    # FIXME: Remove `else' once TP are no longer supported.
    proto = \
        """i32 @"my_func"(i32 %".1", i32 %".2", double %".3", ptr %".4")""" \
        if opaque_pointers_enabled else \
        """i32 @"my_func"(i32 %".1", i32 %".2", double %".3", i32* %".4")"""

    def test_declare(self):
        # A simple declaration
        func = self.function()
        asm = self.descr(func).strip()
        self.assertEqual(asm.strip(), "declare %s" % self.proto)

    def test_declare_attributes(self):
        # Now with function attributes
        func = self.function()
        func.attributes.add("optsize")
        func.attributes.add("alwaysinline")
        func.attributes.add("convergent")
        func.attributes.alignstack = 16
        tp_pers = ir.FunctionType(int8, (), var_arg=True)
        pers = ir.Function(self.module(), tp_pers, '__gxx_personality_v0')
        func.attributes.personality = pers
        asm = self.descr(func).strip()
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(asm,
                             ("declare %s alwaysinline convergent optsize "
                              "alignstack(16) "
                              "personality ptr @\"__gxx_personality_v0\"") %
                             self.proto)
        else:
            self.assertEqual(asm,
                             ("declare %s alwaysinline convergent optsize "
                              "alignstack(16) personality "
                              "i8 (...)* @\"__gxx_personality_v0\"") %
                             self.proto)
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_attributes(self):
        # Now with parameter attributes
        func = self.function()
        func.args[0].add_attribute("zeroext")
        func.args[1].attributes.dereferenceable = 5
        func.args[1].attributes.dereferenceable_or_null = 10
        func.args[3].attributes.align = 4
        func.args[3].add_attribute("nonnull")
        func.return_value.add_attribute("noalias")
        asm = self.descr(func).strip()
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(asm,
                             """declare noalias i32 @"my_func"(i32 zeroext %".1", i32 dereferenceable(5) dereferenceable_or_null(10) %".2", double %".3", ptr nonnull align 4 %".4")"""  # noqa E501
                             )
        else:
            self.assertEqual(asm,
                             """declare noalias i32 @"my_func"(i32 zeroext %".1", i32 dereferenceable(5) dereferenceable_or_null(10) %".2", double %".3", i32* nonnull align 4 %".4")"""  # noqa E501
                             )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_metadata(self):
        # Now with function metadata
        module = self.module()
        func = self.function(module)
        func.set_metadata('dbg', module.add_metadata([]))
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} !dbg !0'
                         )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_section(self):
        # Test function with section
        func = self.function()
        func.section = "a_section"
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} section "a_section"'
                         )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_section_meta(self):
        # Test function with section and metadata
        module = self.module()
        func = self.function(module)
        func.section = "a_section"
        func.set_metadata('dbg', module.add_metadata([]))
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} section "a_section" !dbg !0'
                         )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_attr_meta(self):
        # Test function with attributes and metadata
        module = self.module()
        func = self.function(module)
        func.attributes.add("alwaysinline")
        func.set_metadata('dbg', module.add_metadata([]))
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} alwaysinline !dbg !0'
                         )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_attr_section(self):
        # Test function with attributes and section
        func = self.function()
        func.attributes.add("optsize")
        func.section = "a_section"
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} optsize section "a_section"')
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_function_attr_section_meta(self):
        # Test function with attributes, section and metadata
        module = self.module()
        func = self.function(module)
        func.attributes.add("alwaysinline")
        func.section = "a_section"
        func.set_metadata('dbg', module.add_metadata([]))
        asm = self.descr(func).strip()
        self.assertEqual(asm,
                         f'declare {self.proto} alwaysinline section "a_section" !dbg !0'  # noqa E501
                         )
        # Check pickling
        self.assert_pickle_correctly(func)

    def test_define(self):
        # A simple definition
        func = self.function()
        func.attributes.add("alwaysinline")
        block = func.append_basic_block('my_block')
        builder = ir.IRBuilder(block)
        builder.ret_void()
        asm = self.descr(func)
        self.check_descr(asm, """\
            define {proto} alwaysinline
            {{
            my_block:
                ret void
            }}
            """.format(proto=self.proto))

    def test_declare_intrinsics(self):
        module = self.module()
        pint8 = int8.as_pointer()

        powi = module.declare_intrinsic('llvm.powi', [dbl])
        memset = module.declare_intrinsic('llvm.memset', [pint8, int32])
        memcpy = module.declare_intrinsic('llvm.memcpy', [pint8, pint8, int32])
        assume = module.declare_intrinsic('llvm.assume')
        self.check_descr(self.descr(powi).strip(), """\
            declare double @"llvm.powi.f64"(double %".1", i32 %".2")""")
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_descr(self.descr(memset).strip(), """\
                declare void @"llvm.memset.p0i8.i32"(ptr %".1", i8 %".2", i32 %".3", i1 %".4")""")  # noqa E501
            self.check_descr(self.descr(memcpy).strip(), """\
                declare void @"llvm.memcpy.p0i8.p0i8.i32"(ptr %".1", ptr %".2", i32 %".3", i1 %".4")""")  # noqa E501
        else:
            self.check_descr(self.descr(memset).strip(), """\
                declare void @"llvm.memset.p0i8.i32"(i8* %".1", i8 %".2", i32 %".3", i1 %".4")""")  # noqa E501
            self.check_descr(self.descr(memcpy).strip(), """\
                declare void @"llvm.memcpy.p0i8.p0i8.i32"(i8* %".1", i8* %".2", i32 %".3", i1 %".4")""")  # noqa E501
        self.check_descr(self.descr(assume).strip(), """\
            declare void @"llvm.assume"(i1 %".1")""")

    def test_redeclare_intrinsic(self):
        module = self.module()
        powi = module.declare_intrinsic('llvm.powi', [dbl])
        powi2 = module.declare_intrinsic('llvm.powi', [dbl])
        self.assertIs(powi, powi2)

    def test_pickling(self):
        fn = self.function()
        self.assert_pickle_correctly(fn)

    def test_alwaysinline_noinline_disallowed(self):
        module = self.module()
        func = self.function(module)
        func.attributes.add('alwaysinline')

        msg = "Can't have alwaysinline and noinline"
        with self.assertRaisesRegex(ValueError, msg):
            func.attributes.add('noinline')

    def test_noinline_alwaysinline_disallowed(self):
        module = self.module()
        func = self.function(module)
        func.attributes.add('noinline')

        msg = "Can't have alwaysinline and noinline"
        with self.assertRaisesRegex(ValueError, msg):
            func.attributes.add('alwaysinline')


class TestIR(TestBase):

    def test_unnamed_metadata(self):
        # An unnamed metadata node
        mod = self.module()
        mod.add_metadata([int32(123), int8(42)])
        self.assert_ir_line("!0 = !{ i32 123, i8 42 }", mod)
        self.assert_valid_ir(mod)

    def test_unnamed_metadata_2(self):
        # Several unnamed metadata nodes
        mod = self.module()
        # First node has a literal metadata string
        m0 = mod.add_metadata([int32(123), "kernel"])
        # Second node refers to the first one
        m1 = mod.add_metadata([int64(456), m0])
        # Third node is the same as the second one
        m2 = mod.add_metadata([int64(456), m0])
        self.assertIs(m2, m1)
        # Fourth node refers to the first three
        mod.add_metadata([m0, m1, m2])
        self.assert_ir_line('!0 = !{ i32 123, !"kernel" }', mod)
        self.assert_ir_line('!1 = !{ i64 456, !0 }', mod)
        self.assert_ir_line('!2 = !{ !0, !1, !1 }', mod)

    def test_unnamed_metadata_3(self):
        # Passing nested metadata as a sequence
        mod = self.module()
        mod.add_metadata([int32(123), [int32(456)], [int32(789)], [int32(456)]])
        self.assert_ir_line('!0 = !{ i32 456 }', mod)
        self.assert_ir_line('!1 = !{ i32 789 }', mod)
        self.assert_ir_line('!2 = !{ i32 123, !0, !1, !0 }', mod)

    def test_metadata_string(self):
        # Escaping contents of a metadata string
        mod = self.module()
        mod.add_metadata(["\"\\$"])
        self.assert_ir_line('!0 = !{ !"\\22\\5c$" }', mod)

    def test_named_metadata(self):
        # Add a named metadata node and add metadata values to it
        mod = self.module()
        m0 = mod.add_metadata([int32(123)])
        m1 = mod.add_metadata([int64(456)])
        nmd = mod.add_named_metadata("foo")
        nmd.add(m0)
        nmd.add(m1)
        nmd.add(m0)
        self.assert_ir_line("!foo = !{ !0, !1, !0 }", mod)
        self.assert_valid_ir(mod)
        # Check get_named_metadata()
        self.assertIs(nmd, mod.get_named_metadata("foo"))
        with self.assertRaises(KeyError):
            mod.get_named_metadata("bar")

    def test_named_metadata_2(self):
        # Add and set named metadata through a single add_named_metadata() call
        mod = self.module()
        m0 = mod.add_metadata([int32(123)])
        mod.add_named_metadata("foo", m0)
        mod.add_named_metadata("foo", [int64(456)])
        mod.add_named_metadata("foo", ["kernel"])
        mod.add_named_metadata("bar", [])
        self.assert_ir_line("!foo = !{ !0, !1, !2 }", mod)
        self.assert_ir_line("!0 = !{ i32 123 }", mod)
        self.assert_ir_line("!1 = !{ i64 456 }", mod)
        self.assert_ir_line('!2 = !{ !"kernel" }', mod)
        self.assert_ir_line("!bar = !{ !3 }", mod)
        self.assert_ir_line('!3 = !{  }', mod)
        self.assert_valid_ir(mod)

    def test_metadata_null(self):
        # A null metadata (typed) value
        mod = self.module()
        mod.add_metadata([int32.as_pointer()(None)])
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assert_ir_line("!0 = !{ ptr null }", mod)
        else:
            self.assert_ir_line("!0 = !{ i32* null }", mod)
        self.assert_valid_ir(mod)
        # A null metadata (untyped) value
        mod = self.module()
        mod.add_metadata([None, int32(123)])
        self.assert_ir_line("!0 = !{ null, i32 123 }", mod)
        self.assert_valid_ir(mod)

    def test_debug_info(self):
        # Add real world-looking debug information to a module
        # (with various value types)
        mod = self.module()
        di_file = mod.add_debug_info("DIFile", {
            "filename": "foo",
            "directory": "bar",
        })
        di_func_type = mod.add_debug_info("DISubroutineType", {
            # None as `null`
            "types": mod.add_metadata([None]),
        })
        di_compileunit = mod.add_debug_info("DICompileUnit", {
            "language": ir.DIToken("DW_LANG_Python"),
            "file": di_file,
            "producer": "ARTIQ",
            "runtimeVersion": 0,
            "isOptimized": True,
        }, is_distinct=True)
        mod.add_debug_info("DISubprogram", {
            "name": "my_func",
            "file": di_file,
            "line": 11,
            "type": di_func_type,
            "isLocal": False,
            "unit": di_compileunit,
        }, is_distinct=True)

        # Check output
        strmod = str(mod)
        self.assert_ir_line('!0 = !DIFile(directory: "bar", filename: "foo")',
                            strmod)
        self.assert_ir_line('!1 = !{ null }', strmod)
        self.assert_ir_line('!2 = !DISubroutineType(types: !1)', strmod)
        # self.assert_ir_line('!4 = !{ !3 }', strmod)
        self.assert_ir_line('!3 = distinct !DICompileUnit(file: !0, '
                            'isOptimized: true, language: DW_LANG_Python, '
                            'producer: "ARTIQ", runtimeVersion: 0)',
                            strmod)
        self.assert_ir_line('!4 = distinct !DISubprogram(file: !0, isLocal: '
                            'false, line: 11, name: "my_func", type: !2, unit: '
                            '!3)',
                            strmod)
        self.assert_valid_ir(mod)

    def test_debug_info_2(self):
        # Identical debug info nodes should be merged
        mod = self.module()
        di1 = mod.add_debug_info("DIFile",
                                 {"filename": "foo",
                                  "directory": "bar",
                                  })
        di2 = mod.add_debug_info("DIFile",
                                 {"filename": "foo",
                                  "directory": "bar",
                                  })
        di3 = mod.add_debug_info("DIFile",
                                 {"filename": "bar",
                                  "directory": "foo",
                                  })
        di4 = mod.add_debug_info("DIFile",
                                 {"filename": "foo",
                                  "directory": "bar",
                                  }, is_distinct=True)
        self.assertIs(di1, di2)
        self.assertEqual(len({di1, di2, di3, di4}), 3)
        # Check output
        strmod = str(mod)
        self.assert_ir_line('!0 = !DIFile(directory: "bar", filename: "foo")',
                            strmod)
        self.assert_ir_line('!1 = !DIFile(directory: "foo", filename: "bar")',
                            strmod)
        self.assert_ir_line('!2 = distinct !DIFile(directory: "bar", filename: '
                            '"foo")', strmod)
        self.assert_valid_ir(mod)

    def test_debug_info_gvar(self):
        # This test defines a module with a global variable named 'gvar'.
        # When the module is compiled and linked with a main function, gdb can
        # be used to interpret and print the the value of 'gvar'.
        mod = self.module()

        gvar = ir.GlobalVariable(mod, ir.FloatType(), 'gvar')
        gvar.initializer = ir.Constant(ir.FloatType(), 42)

        di_float = mod.add_debug_info("DIBasicType", {
            "name": "float",
            "size": 32,
            "encoding": ir.DIToken("DW_ATE_float")
        })
        di_gvar = mod.add_debug_info("DIGlobalVariableExpression", {
            "expr": mod.add_debug_info("DIExpression", {}),
            "var": mod.add_debug_info("DIGlobalVariable", {
                "name": gvar.name,
                "type": di_float,
                "isDefinition": True
            }, is_distinct=True)
        })
        gvar.set_metadata('dbg', di_gvar)

        # Check output
        strmod = str(mod)
        self.assert_ir_line('!0 = !DIBasicType(encoding: DW_ATE_float, '
                            'name: "float", size: 32)', strmod)
        self.assert_ir_line('!1 = !DIExpression()', strmod)
        self.assert_ir_line('!2 = distinct !DIGlobalVariable(isDefinition: '
                            'true, name: "gvar", type: !0)', strmod)
        self.assert_ir_line('!3 = !DIGlobalVariableExpression(expr: !1, '
                            'var: !2)', strmod)
        self.assert_ir_line('@"gvar" = global float 0x4045000000000000, '
                            '!dbg !3', strmod)

        # The remaining debug info is not part of the automated test, but
        # can be used to produce an object file that can be loaded into a
        # debugger to print the value of gvar. This can be done by printing the
        # module then compiling it with clang and inspecting with gdb:
        #
        #     clang test_debug_info_gvar.ll -c
        #     printf "file test_debug_info_gvar.o \n p gvar" | gdb
        #
        # Which should result in the output:
        #
        #     (gdb) $1 = 42

        dver = [ir.IntType(32)(2), 'Dwarf Version', ir.IntType(32)(4)]
        diver = [ir.IntType(32)(2), 'Debug Info Version', ir.IntType(32)(3)]
        dver = mod.add_metadata(dver)
        diver = mod.add_metadata(diver)
        flags = mod.add_named_metadata('llvm.module.flags')
        flags.add(dver)
        flags.add(diver)

        di_file = mod.add_debug_info("DIFile", {
            "filename": "foo",
            "directory": "bar",
        })
        di_cu = mod.add_debug_info("DICompileUnit", {
            "language": ir.DIToken("DW_LANG_Python"),
            "file": di_file,
            'emissionKind': ir.DIToken('FullDebug'),
            "globals": mod.add_metadata([di_gvar])
        }, is_distinct=True)
        mod.add_named_metadata('llvm.dbg.cu', di_cu)

    def test_debug_info_unicode_string(self):
        mod = self.module()
        mod.add_debug_info("DILocalVariable", {"name": "a∆"})
        # Check output
        strmod = str(mod)
        # The unicode character is utf8 encoded with \XX format, where XX is hex
        name = ''.join(map(lambda x: f"\\{x:02x}", "∆".encode()))
        self.assert_ir_line(f'!0 = !DILocalVariable(name: "a{name}")', strmod)

    def test_inline_assembly(self):
        mod = self.module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'foo')
        builder = ir.IRBuilder(foo.append_basic_block(''))
        asmty = ir.FunctionType(int32, [int32])
        asm = ir.InlineAsm(asmty, "mov $1, $2", "=r,r", side_effect=True)
        builder.call(asm, [int32(123)])
        builder.ret_void()
        pat = 'call i32 asm sideeffect "mov $1, $2", "=r,r" ( i32 123 )'
        self.assertInText(pat, str(mod))
        self.assert_valid_ir(mod)

    def test_builder_asm(self):
        mod = self.module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'foo')
        builder = ir.IRBuilder(foo.append_basic_block(''))
        asmty = ir.FunctionType(int32, [int32])
        builder.asm(asmty, "mov $1, $2", "=r,r", [int32(123)], side_effect=True)
        builder.ret_void()
        pat = 'call i32 asm sideeffect "mov $1, $2", "=r,r" ( i32 123 )'
        self.assertInText(pat, str(mod))
        self.assert_valid_ir(mod)

    def test_builder_load_reg(self):
        mod = self.module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'foo')
        builder = ir.IRBuilder(foo.append_basic_block(''))
        builder.load_reg(ir.IntType(64), "rax")
        builder.ret_void()
        pat = 'call i64 asm "", "={rax}"'
        self.assertInText(pat, str(mod))
        self.assert_valid_ir(mod)

    def test_builder_store_reg(self):
        mod = self.module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'foo')
        builder = ir.IRBuilder(foo.append_basic_block(''))
        builder.store_reg(int64(123), ir.IntType(64), "rax")
        builder.ret_void()
        pat = 'call void asm sideeffect "", "{rax}" ( i64 123 )'
        self.assertInText(pat, str(mod))
        self.assert_valid_ir(mod)


class TestGlobalValues(TestBase):

    def test_globals_access(self):
        mod = self.module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'foo')
        ir.Function(mod, ir.FunctionType(ir.VoidType(), []), 'bar')
        globdouble = ir.GlobalVariable(mod, ir.DoubleType(), 'globdouble')
        self.assertEqual(mod.get_global('foo'), foo)
        self.assertEqual(mod.get_global('globdouble'), globdouble)
        with self.assertRaises(KeyError):
            mod.get_global('kkk')
        # Globals should have a useful repr()
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(repr(globdouble),
                             "<ir.GlobalVariable 'globdouble' of type 'ptr'>")
        else:
            self.assertEqual(
                repr(globdouble),
                "<ir.GlobalVariable 'globdouble' of type 'double*'>")

    def test_functions_global_values_access(self):
        """
        Accessing functions and global values through Module.functions
        and Module.global_values.
        """
        mod = self.module()
        fty = ir.FunctionType(ir.VoidType(), [])
        foo = ir.Function(mod, fty, 'foo')
        bar = ir.Function(mod, fty, 'bar')
        globdouble = ir.GlobalVariable(mod, ir.DoubleType(), 'globdouble')
        self.assertEqual(set(mod.functions), set((foo, bar)))
        self.assertEqual(set(mod.global_values), set((foo, bar, globdouble)))

    def test_global_variables_ir(self):
        """
        IR serialization of global variables.
        """
        mod = self.module()
        # the following have side effects and write to self.module()
        a = ir.GlobalVariable(mod, int8, 'a')  # noqa F841
        b = ir.GlobalVariable(mod, int8, 'b', addrspace=42)  # noqa F841
        # Initialized global variable doesn't default to "external"
        c = ir.GlobalVariable(mod, int32, 'c')
        c.initializer = int32(123)
        d = ir.GlobalVariable(mod, int32, 'd')
        d.global_constant = True
        # Non-external linkage implies default "undef" initializer
        e = ir.GlobalVariable(mod, int32, 'e')
        e.linkage = "internal"
        f = ir.GlobalVariable(mod, int32, 'f', addrspace=456)
        f.unnamed_addr = True
        g = ir.GlobalVariable(mod, int32, 'g')
        g.linkage = "internal"
        g.initializer = int32(123)
        g.align = 16
        h = ir.GlobalVariable(mod, int32, 'h')
        h.linkage = "internal"
        h.initializer = int32(123)
        h.section = "h_section"
        i = ir.GlobalVariable(mod, int32, 'i')
        i.linkage = "internal"
        i.initializer = int32(456)
        i.align = 8
        i.section = "i_section"
        self.check_module_body(mod, """\
            @"a" = external global i8
            @"b" = external addrspace(42) global i8
            @"c" = global i32 123
            @"d" = external constant i32
            @"e" = internal global i32 undef
            @"f" = external unnamed_addr addrspace(456) global i32
            @"g" = internal global i32 123, align 16
            @"h" = internal global i32 123, section "h_section"
            @"i" = internal global i32 456, section "i_section", align 8
            """)

    def test_pickle(self):
        mod = self.module()
        self.assert_pickle_correctly(mod)


class TestBlock(TestBase):

    def test_attributes(self):
        func = self.function()
        block = ir.Block(parent=func, name='start')
        self.assertIs(block.parent, func)
        self.assertFalse(block.is_terminated)

    def test_descr(self):
        block = self.block(name='my_block')
        self.assertEqual(self.descr(block), "my_block:\n")
        block.instructions.extend(['a', 'b'])
        self.assertEqual(self.descr(block), "my_block:\n  a\n  b\n")

    def test_replace(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        c = builder.add(a, b, 'c')
        d = builder.sub(a, b, 'd')
        builder.mul(d, b, 'e')
        f = ir.Instruction(block, a.type, 'sdiv', (c, b), 'f')
        self.check_block(block, """\
            my_block:
                %"c" = add i32 %".1", %".2"
                %"d" = sub i32 %".1", %".2"
                %"e" = mul i32 %"d", %".2"
            """)
        block.replace(d, f)
        self.check_block(block, """\
            my_block:
                %"c" = add i32 %".1", %".2"
                %"f" = sdiv i32 %"c", %".2"
                %"e" = mul i32 %"f", %".2"
            """)

    def test_repr(self):
        """
        Blocks should have a useful repr()
        """
        func = self.function()
        block = ir.Block(parent=func, name='start')
        self.assertEqual(repr(block), "<ir.Block 'start' of type 'label'>")


class TestBuildInstructions(TestBase):
    """
    Test IR generation of LLVM instructions through the IRBuilder class.
    """

    maxDiff = 4000

    def test_simple(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        inst = builder.add(a, b, 'res')
        self.check_block(block, """\
            my_block:
                %"res" = add i32 %".1", %".2"
            """)
        # Instructions should have a useful repr()
        self.assertEqual(repr(inst),
                         "<ir.Instruction 'res' of type 'i32', opname 'add', "
                         "operands (<ir.Argument '.1' of type i32>, "
                         "<ir.Argument '.2' of type i32>)>")

    def test_binops(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b, ff = builder.function.args[:3]
        builder.add(a, b, 'c')
        builder.fadd(a, b, 'd')
        builder.sub(a, b, 'e')
        builder.fsub(a, b, 'f')
        builder.mul(a, b, 'g')
        builder.fmul(a, b, 'h')
        builder.udiv(a, b, 'i')
        builder.sdiv(a, b, 'j')
        builder.fdiv(a, b, 'k')
        builder.urem(a, b, 'l')
        builder.srem(a, b, 'm')
        builder.frem(a, b, 'n')
        builder.or_(a, b, 'o')
        builder.and_(a, b, 'p')
        builder.xor(a, b, 'q')
        builder.shl(a, b, 'r')
        builder.ashr(a, b, 's')
        builder.lshr(a, b, 't')
        with self.assertRaises(ValueError) as cm:
            builder.add(a, ff)
        self.assertEqual(str(cm.exception),
                         "Operands must be the same type, got (i32, double)")
        self.assertFalse(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"c" = add i32 %".1", %".2"
                %"d" = fadd i32 %".1", %".2"
                %"e" = sub i32 %".1", %".2"
                %"f" = fsub i32 %".1", %".2"
                %"g" = mul i32 %".1", %".2"
                %"h" = fmul i32 %".1", %".2"
                %"i" = udiv i32 %".1", %".2"
                %"j" = sdiv i32 %".1", %".2"
                %"k" = fdiv i32 %".1", %".2"
                %"l" = urem i32 %".1", %".2"
                %"m" = srem i32 %".1", %".2"
                %"n" = frem i32 %".1", %".2"
                %"o" = or i32 %".1", %".2"
                %"p" = and i32 %".1", %".2"
                %"q" = xor i32 %".1", %".2"
                %"r" = shl i32 %".1", %".2"
                %"s" = ashr i32 %".1", %".2"
                %"t" = lshr i32 %".1", %".2"
            """)

    def test_binop_flags(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        # As tuple
        builder.add(a, b, 'c', flags=('nuw',))
        # and as list
        builder.sub(a, b, 'd', flags=['nuw', 'nsw'])
        self.check_block(block, """\
            my_block:
                %"c" = add nuw i32 %".1", %".2"
                %"d" = sub nuw nsw i32 %".1", %".2"
            """)

    def test_binop_fastmath_flags(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        # As tuple
        builder.fadd(a, b, 'c', flags=('fast',))
        # and as list
        builder.fsub(a, b, 'd', flags=['ninf', 'nsz'])
        self.check_block(block, """\
            my_block:
                %"c" = fadd fast i32 %".1", %".2"
                %"d" = fsub ninf nsz i32 %".1", %".2"
            """)

    def test_binops_with_overflow(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        builder.sadd_with_overflow(a, b, 'c')
        builder.smul_with_overflow(a, b, 'd')
        builder.ssub_with_overflow(a, b, 'e')
        builder.uadd_with_overflow(a, b, 'f')
        builder.umul_with_overflow(a, b, 'g')
        builder.usub_with_overflow(a, b, 'h')
        self.check_block(block, """\
my_block:
    %"c" = call {i32, i1} @"llvm.sadd.with.overflow.i32"(i32 %".1", i32 %".2")
    %"d" = call {i32, i1} @"llvm.smul.with.overflow.i32"(i32 %".1", i32 %".2")
    %"e" = call {i32, i1} @"llvm.ssub.with.overflow.i32"(i32 %".1", i32 %".2")
    %"f" = call {i32, i1} @"llvm.uadd.with.overflow.i32"(i32 %".1", i32 %".2")
    %"g" = call {i32, i1} @"llvm.umul.with.overflow.i32"(i32 %".1", i32 %".2")
    %"h" = call {i32, i1} @"llvm.usub.with.overflow.i32"(i32 %".1", i32 %".2")
            """)

    def test_unary_ops(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b, c = builder.function.args[:3]
        builder.neg(a, 'd')
        builder.not_(b, 'e')
        builder.fneg(c, 'f')
        self.assertFalse(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"d" = sub i32 0, %".1"
                %"e" = xor i32 %".2", -1
                %"f" = fneg double %".3"
            """)

    def test_replace_operand(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        undef1 = ir.Constant(ir.IntType(32), ir.Undefined)
        undef2 = ir.Constant(ir.IntType(32), ir.Undefined)
        c = builder.add(undef1, undef2, 'c')
        self.check_block(block, """\
            my_block:
                %"c" = add i32 undef, undef
            """)
        c.replace_usage(undef1, a)
        c.replace_usage(undef2, b)
        self.check_block(block, """\
            my_block:
                %"c" = add i32 %".1", %".2"
            """)

    def test_integer_comparisons(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        builder.icmp_unsigned('==', a, b, 'c')
        builder.icmp_unsigned('!=', a, b, 'd')
        builder.icmp_unsigned('<', a, b, 'e')
        builder.icmp_unsigned('<=', a, b, 'f')
        builder.icmp_unsigned('>', a, b, 'g')
        builder.icmp_unsigned('>=', a, b, 'h')
        builder.icmp_signed('==', a, b, 'i')
        builder.icmp_signed('!=', a, b, 'j')
        builder.icmp_signed('<', a, b, 'k')
        builder.icmp_signed('<=', a, b, 'l')
        builder.icmp_signed('>', a, b, 'm')
        builder.icmp_signed('>=', a, b, 'n')
        with self.assertRaises(ValueError):
            builder.icmp_signed('uno', a, b, 'zz')
        with self.assertRaises(ValueError):
            builder.icmp_signed('foo', a, b, 'zz')
        self.assertFalse(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"c" = icmp eq i32 %".1", %".2"
                %"d" = icmp ne i32 %".1", %".2"
                %"e" = icmp ult i32 %".1", %".2"
                %"f" = icmp ule i32 %".1", %".2"
                %"g" = icmp ugt i32 %".1", %".2"
                %"h" = icmp uge i32 %".1", %".2"
                %"i" = icmp eq i32 %".1", %".2"
                %"j" = icmp ne i32 %".1", %".2"
                %"k" = icmp slt i32 %".1", %".2"
                %"l" = icmp sle i32 %".1", %".2"
                %"m" = icmp sgt i32 %".1", %".2"
                %"n" = icmp sge i32 %".1", %".2"
            """)

    def test_float_comparisons(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        builder.fcmp_ordered('==', a, b, 'c')
        builder.fcmp_ordered('!=', a, b, 'd')
        builder.fcmp_ordered('<', a, b, 'e')
        builder.fcmp_ordered('<=', a, b, 'f')
        builder.fcmp_ordered('>', a, b, 'g')
        builder.fcmp_ordered('>=', a, b, 'h')
        builder.fcmp_unordered('==', a, b, 'i')
        builder.fcmp_unordered('!=', a, b, 'j')
        builder.fcmp_unordered('<', a, b, 'k')
        builder.fcmp_unordered('<=', a, b, 'l')
        builder.fcmp_unordered('>', a, b, 'm')
        builder.fcmp_unordered('>=', a, b, 'n')
        # fcmp_ordered and fcmp_unordered are the same for these cases
        builder.fcmp_ordered('ord', a, b, 'u')
        builder.fcmp_ordered('uno', a, b, 'v')
        builder.fcmp_unordered('ord', a, b, 'w')
        builder.fcmp_unordered('uno', a, b, 'x')
        builder.fcmp_unordered('olt', a, b, 'y',
                               flags=['nnan', 'ninf', 'nsz', 'arcp', 'fast'])
        self.assertFalse(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"c" = fcmp oeq i32 %".1", %".2"
                %"d" = fcmp one i32 %".1", %".2"
                %"e" = fcmp olt i32 %".1", %".2"
                %"f" = fcmp ole i32 %".1", %".2"
                %"g" = fcmp ogt i32 %".1", %".2"
                %"h" = fcmp oge i32 %".1", %".2"
                %"i" = fcmp ueq i32 %".1", %".2"
                %"j" = fcmp une i32 %".1", %".2"
                %"k" = fcmp ult i32 %".1", %".2"
                %"l" = fcmp ule i32 %".1", %".2"
                %"m" = fcmp ugt i32 %".1", %".2"
                %"n" = fcmp uge i32 %".1", %".2"
                %"u" = fcmp ord i32 %".1", %".2"
                %"v" = fcmp uno i32 %".1", %".2"
                %"w" = fcmp ord i32 %".1", %".2"
                %"x" = fcmp uno i32 %".1", %".2"
                %"y" = fcmp nnan ninf nsz arcp fast olt i32 %".1", %".2"
            """)

    def test_misc_ops(self):
        block = self.block(name='my_block')
        t = ir.Constant(int1, True)
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        builder.select(t, a, b, 'c', flags=('arcp', 'nnan'))
        self.assertFalse(block.is_terminated)
        builder.unreachable()
        self.assertTrue(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"c" = select arcp nnan i1 true, i32 %".1", i32 %".2"
                unreachable
            """)

    def test_phi(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        bb2 = builder.function.append_basic_block('b2')
        bb3 = builder.function.append_basic_block('b3')
        phi = builder.phi(int32, 'my_phi', flags=('fast',))
        phi.add_incoming(a, bb2)
        phi.add_incoming(b, bb3)
        self.assertFalse(block.is_terminated)
        self.check_block(block, """\
            my_block:
                %"my_phi" = phi fast i32 [%".1", %"b2"], [%".2", %"b3"]
            """)

    def test_mem_ops(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b, z = builder.function.args[:3]
        c = builder.alloca(int32, name='c')
        d = builder.alloca(int32, size=42, name='d')  # noqa F841
        e = builder.alloca(dbl, size=a, name='e')
        e.align = 8
        self.assertEqual(e.type, ir.PointerType(dbl))
        ee = builder.store(z, e)
        self.assertEqual(ee.type, ir.VoidType())
        f = builder.store(b, c)
        self.assertEqual(f.type, ir.VoidType())
        g = builder.load(c, 'g')
        self.assertEqual(g.type, int32)
        # With alignment
        h = builder.store(b, c, align=1)
        self.assertEqual(h.type, ir.VoidType())
        i = builder.load(c, 'i', align=1)
        self.assertEqual(i.type, int32)
        # Atomics
        j = builder.store_atomic(b, c, ordering="seq_cst", align=4)
        self.assertEqual(j.type, ir.VoidType())
        k = builder.load_atomic(c, ordering="seq_cst", align=4, name='k')
        self.assertEqual(k.type, int32)
        if opaque_pointers_enabled:
            ptr = ir.Constant(ir.PointerType(), None)
        else:
            ptr = ir.Constant(ir.PointerType(int32), None)
        builder.store(ir.Constant(int32, 5), ptr)
        # Not pointer types
        with self.assertRaises(TypeError):
            builder.store(b, a)
        with self.assertRaises(TypeError):
            builder.load(b)
        # Mismatching pointer type
        with self.assertRaises(TypeError) as cm:
            builder.store(b, e)

        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(cm.exception),
                             "cannot store i32 to ptr: mismatching types")
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32
                    %"d" = alloca i32, i32 42
                    %"e" = alloca double, i32 %".1", align 8
                    store double %".3", ptr %"e"
                    store i32 %".2", ptr %"c"
                    %"g" = load i32, ptr %"c"
                    store i32 %".2", ptr %"c", align 1
                    %"i" = load i32, ptr %"c", align 1
                    store atomic i32 %".2", ptr %"c" seq_cst, align 4
                    %"k" = load atomic i32, ptr %"c" seq_cst, align 4
                    store i32 5, ptr null
                """)
        else:
            self.assertEqual(str(cm.exception),
                             "cannot store i32 to double*: mismatching types")
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32
                    %"d" = alloca i32, i32 42
                    %"e" = alloca double, i32 %".1", align 8
                    store double %".3", double* %"e"
                    store i32 %".2", i32* %"c"
                    %"g" = load i32, i32* %"c"
                    store i32 %".2", i32* %"c", align 1
                    %"i" = load i32, i32* %"c", align 1
                    store atomic i32 %".2", i32* %"c" seq_cst, align 4
                    %"k" = load atomic i32, i32* %"c" seq_cst, align 4
                    store i32 5, i32* null
                """)

    def test_gep(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        c = builder.alloca(ir.PointerType(int32), name='c')
        d = builder.gep(c, [ir.Constant(int32, 5), a], name='d')
        self.assertEqual(d.type, ir.PointerType(int32))
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca ptr
                    %"d" = getelementptr ptr, ptr %"c", i32 5, i32 %".1"
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32*
                    %"d" = getelementptr i32*, i32** %"c", i32 5, i32 %".1"
                """)
        # XXX test with more complex types

    def test_gep_castinstr(self):
        # similar to:
        # numba::runtime::nrtdynmod.py_define_nrt_meminfo_data()
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        int8ptr = int8.as_pointer()
        ls = ir.LiteralStructType([int64, int8ptr, int8ptr, int8ptr, int64])
        d = builder.bitcast(a, ls.as_pointer(), name='d')
        e = builder.gep(d, [ir.Constant(int32, x) for x in [0, 3]], name='e')
        self.assertEqual(e.type, ir.PointerType(int8ptr))
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"d" = bitcast i32 %".1" to ptr
                    %"e" = getelementptr {i64, ptr, ptr, ptr, i64}, ptr %"d", i32 0, i32 3
                """)  # noqa E501
        else:
            self.check_block(block, """\
                my_block:
                    %"d" = bitcast i32 %".1" to {i64, i8*, i8*, i8*, i64}*
                    %"e" = getelementptr {i64, i8*, i8*, i8*, i64}, {i64, i8*, i8*, i8*, i64}* %"d", i32 0, i32 3
                """)  # noqa E501

    def test_gep_castinstr_addrspace(self):
        # similar to:
        # numba::runtime::nrtdynmod.py_define_nrt_meminfo_data()
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        addrspace = 4
        int8ptr = int8.as_pointer()
        ls = ir.LiteralStructType([int64, int8ptr, int8ptr, int8ptr, int64])
        d = builder.bitcast(a, ls.as_pointer(addrspace=addrspace), name='d')
        e = builder.gep(d, [ir.Constant(int32, x) for x in [0, 3]], name='e')
        self.assertEqual(e.type.addrspace, addrspace)
        self.assertEqual(e.type, ir.PointerType(int8ptr, addrspace=addrspace))
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"d" = bitcast i32 %".1" to ptr addrspace(4)
                    %"e" = getelementptr {i64, ptr, ptr, ptr, i64}, ptr addrspace(4) %"d", i32 0, i32 3
                """)  # noqa E501
        else:
            self.check_block(block, """\
                my_block:
                    %"d" = bitcast i32 %".1" to {i64, i8*, i8*, i8*, i64} addrspace(4)*
                    %"e" = getelementptr {i64, i8*, i8*, i8*, i64}, {i64, i8*, i8*, i8*, i64} addrspace(4)* %"d", i32 0, i32 3
                """)  # noqa E501

    def test_gep_addrspace(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        addrspace = 4
        c = builder.alloca(ir.PointerType(int32, addrspace=addrspace), name='c')
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c.type), 'ptr')
        else:
            self.assertEqual(str(c.type), 'i32 addrspace(4)**')
        self.assertEqual(c.type.pointee.addrspace, addrspace)
        d = builder.gep(c, [ir.Constant(int32, 5), a], name='d')
        self.assertEqual(d.type.addrspace, addrspace)
        e = builder.gep(d, [ir.Constant(int32, 10)], name='e')
        self.assertEqual(e.type.addrspace, addrspace)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca ptr addrspace(4)
                    %"d" = getelementptr ptr addrspace(4), ptr %"c", i32 5, i32 %".1"
                    %"e" = getelementptr i32, ptr addrspace(4) %"d", i32 10
                """)  # noqa E501
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32 addrspace(4)*
                    %"d" = getelementptr i32 addrspace(4)*, i32 addrspace(4)** %"c", i32 5, i32 %".1"
                    %"e" = getelementptr i32, i32 addrspace(4)* %"d", i32 10
                """)  # noqa E501

    def test_extract_insert_value(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        tp_inner = ir.LiteralStructType([int32, int1])
        tp_outer = ir.LiteralStructType([int8, tp_inner])
        c_inner = ir.Constant(tp_inner, (ir.Constant(int32, 4),
                                         ir.Constant(int1, True)))
        # Flat structure
        c = builder.extract_value(c_inner, 0, name='c')  # noqa F841
        d = builder.insert_value(c_inner, a, 0, name='d')  # noqa F841
        e = builder.insert_value(d, ir.Constant(int1, False), 1, name='e')  # noqa F841 E501
        self.assertEqual(d.type, tp_inner)
        self.assertEqual(e.type, tp_inner)
        # Nested structure
        p_outer = builder.alloca(tp_outer, name='ptr')
        j = builder.load(p_outer, name='j')
        k = builder.extract_value(j, 0, name='k')
        l = builder.extract_value(j, 1, name='l')
        m = builder.extract_value(j, (1, 0), name='m')
        n = builder.extract_value(j, (1, 1), name='n')
        o = builder.insert_value(j, l, 1, name='o')
        p = builder.insert_value(j, a, (1, 0), name='p')
        self.assertEqual(k.type, int8)
        self.assertEqual(l.type, tp_inner)
        self.assertEqual(m.type, int32)
        self.assertEqual(n.type, int1)
        self.assertEqual(o.type, tp_outer)
        self.assertEqual(p.type, tp_outer)

        with self.assertRaises(TypeError):
            # Not an aggregate
            builder.extract_value(p_outer, 0)
        with self.assertRaises(TypeError):
            # Indexing too deep
            builder.extract_value(c_inner, (0, 0))
        with self.assertRaises(TypeError):
            # Index out of structure bounds
            builder.extract_value(c_inner, 5)
        with self.assertRaises(TypeError):
            # Not an aggregate
            builder.insert_value(a, b, 0)
        with self.assertRaises(TypeError):
            # Replacement value has the wrong type
            builder.insert_value(c_inner, a, 1)

        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = extractvalue {i32, i1} {i32 4, i1 true}, 0
                    %"d" = insertvalue {i32, i1} {i32 4, i1 true}, i32 %".1", 0
                    %"e" = insertvalue {i32, i1} %"d", i1 false, 1
                    %"ptr" = alloca {i8, {i32, i1}}
                    %"j" = load {i8, {i32, i1}}, ptr %"ptr"
                    %"k" = extractvalue {i8, {i32, i1}} %"j", 0
                    %"l" = extractvalue {i8, {i32, i1}} %"j", 1
                    %"m" = extractvalue {i8, {i32, i1}} %"j", 1, 0
                    %"n" = extractvalue {i8, {i32, i1}} %"j", 1, 1
                    %"o" = insertvalue {i8, {i32, i1}} %"j", {i32, i1} %"l", 1
                    %"p" = insertvalue {i8, {i32, i1}} %"j", i32 %".1", 1, 0
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = extractvalue {i32, i1} {i32 4, i1 true}, 0
                    %"d" = insertvalue {i32, i1} {i32 4, i1 true}, i32 %".1", 0
                    %"e" = insertvalue {i32, i1} %"d", i1 false, 1
                    %"ptr" = alloca {i8, {i32, i1}}
                    %"j" = load {i8, {i32, i1}}, {i8, {i32, i1}}* %"ptr"
                    %"k" = extractvalue {i8, {i32, i1}} %"j", 0
                    %"l" = extractvalue {i8, {i32, i1}} %"j", 1
                    %"m" = extractvalue {i8, {i32, i1}} %"j", 1, 0
                    %"n" = extractvalue {i8, {i32, i1}} %"j", 1, 1
                    %"o" = insertvalue {i8, {i32, i1}} %"j", {i32, i1} %"l", 1
                    %"p" = insertvalue {i8, {i32, i1}} %"j", i32 %".1", 1, 0
                """)

    def test_cast_ops(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b, fa, ptr = builder.function.args[:4]
        c = builder.trunc(a, int8, name='c')
        d = builder.zext(c, int32, name='d')  # noqa F841
        e = builder.sext(c, int32, name='e')  # noqa F841
        fb = builder.fptrunc(fa, flt, 'fb')
        fc = builder.fpext(fb, dbl, 'fc')  # noqa F841
        g = builder.fptoui(fa, int32, 'g')
        h = builder.fptosi(fa, int8, 'h')
        fd = builder.uitofp(g, flt, 'fd')  # noqa F841
        fe = builder.sitofp(h, dbl, 'fe')  # noqa F841
        i = builder.ptrtoint(ptr, int32, 'i')
        j = builder.inttoptr(i, ir.PointerType(int8), 'j')  # noqa F841
        k = builder.bitcast(a, flt, "k")  # noqa F841
        self.assertFalse(block.is_terminated)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = trunc i32 %".1" to i8
                    %"d" = zext i8 %"c" to i32
                    %"e" = sext i8 %"c" to i32
                    %"fb" = fptrunc double %".3" to float
                    %"fc" = fpext float %"fb" to double
                    %"g" = fptoui double %".3" to i32
                    %"h" = fptosi double %".3" to i8
                    %"fd" = uitofp i32 %"g" to float
                    %"fe" = sitofp i8 %"h" to double
                    %"i" = ptrtoint ptr %".4" to i32
                    %"j" = inttoptr i32 %"i" to ptr
                    %"k" = bitcast i32 %".1" to float
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = trunc i32 %".1" to i8
                    %"d" = zext i8 %"c" to i32
                    %"e" = sext i8 %"c" to i32
                    %"fb" = fptrunc double %".3" to float
                    %"fc" = fpext float %"fb" to double
                    %"g" = fptoui double %".3" to i32
                    %"h" = fptosi double %".3" to i8
                    %"fd" = uitofp i32 %"g" to float
                    %"fe" = sitofp i8 %"h" to double
                    %"i" = ptrtoint i32* %".4" to i32
                    %"j" = inttoptr i32 %"i" to i8*
                    %"k" = bitcast i32 %".1" to float
                """)

    def test_atomicrmw(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        c = builder.alloca(int32, name='c')
        d = builder.atomic_rmw('add', c, a, 'monotonic', 'd')
        self.assertEqual(d.type, int32)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32
                    %"d" = atomicrmw add ptr %"c", i32 %".1" monotonic
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32
                    %"d" = atomicrmw add i32* %"c", i32 %".1" monotonic
                """)

    def test_branch(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        bb_target = builder.function.append_basic_block(name='target')
        builder.branch(bb_target)
        self.assertTrue(block.is_terminated)
        self.check_block(block, """\
            my_block:
                br label %"target"
            """)

    def test_cbranch(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        bb_true = builder.function.append_basic_block(name='b_true')
        bb_false = builder.function.append_basic_block(name='b_false')
        builder.cbranch(ir.Constant(int1, False), bb_true, bb_false)
        self.assertTrue(block.is_terminated)
        self.check_block(block, """\
            my_block:
                br i1 false, label %"b_true", label %"b_false"
            """)

    def test_cbranch_weights(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        bb_true = builder.function.append_basic_block(name='b_true')
        bb_false = builder.function.append_basic_block(name='b_false')
        br = builder.cbranch(ir.Constant(int1, False), bb_true, bb_false)
        br.set_weights([5, 42])
        self.assertTrue(block.is_terminated)
        self.check_block(block, """\
            my_block:
                br i1 false, label %"b_true", label %"b_false", !prof !0
            """)
        self.check_metadata(builder.module, """\
            !0 = !{ !"branch_weights", i32 5, i32 42 }
            """)

    def test_branch_indirect(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        bb_1 = builder.function.append_basic_block(name='b_1')
        bb_2 = builder.function.append_basic_block(name='b_2')
        indirectbr = builder.branch_indirect(
            ir.BlockAddress(builder.function, bb_1))
        indirectbr.add_destination(bb_1)
        indirectbr.add_destination(bb_2)
        self.assertTrue(block.is_terminated)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    indirectbr ptr blockaddress(@"my_func", %"b_1"), [label %"b_1", label %"b_2"]
                """)  # noqa E501
        else:
            self.check_block(block, """\
                my_block:
                    indirectbr i8* blockaddress(@"my_func", %"b_1"), [label %"b_1", label %"b_2"]
                """)  # noqa E501

    def test_returns(self):
        def check(block, expected_ir):
            self.assertTrue(block.is_terminated)
            self.check_block(block, expected_ir)

        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        builder.ret_void()
        check(block, """\
            my_block:
                ret void
            """)

        block = self.block(name='other_block')
        builder = ir.IRBuilder(block)
        builder.ret(int32(5))
        check(block, """\
            other_block:
                ret i32 5
            """)

        # With metadata
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        inst = builder.ret_void()
        inst.set_metadata("dbg", block.module.add_metadata(()))
        check(block, """\
            my_block:
                ret void, !dbg !0
            """)

        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        inst = builder.ret(int32(6))
        inst.set_metadata("dbg", block.module.add_metadata(()))
        check(block, """\
            my_block:
                ret i32 6, !dbg !0
            """)

    def test_switch(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        bb_onzero = builder.function.append_basic_block(name='onzero')
        bb_onone = builder.function.append_basic_block(name='onone')
        bb_ontwo = builder.function.append_basic_block(name='ontwo')
        bb_else = builder.function.append_basic_block(name='otherwise')
        sw = builder.switch(a, bb_else)
        sw.add_case(ir.Constant(int32, 0), bb_onzero)
        sw.add_case(ir.Constant(int32, 1), bb_onone)
        # A plain Python value gets converted into the right IR constant
        sw.add_case(2, bb_ontwo)
        self.assertTrue(block.is_terminated)
        self.check_block(block, """\
            my_block:
                switch i32 %".1", label %"otherwise" [i32 0, label %"onzero" i32 1, label %"onone" i32 2, label %"ontwo"]
            """)  # noqa E501

    def test_call(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        tp_f = ir.FunctionType(flt, (int32, int32))
        tp_g = ir.FunctionType(dbl, (int32,), var_arg=True)
        tp_h = ir.FunctionType(hlf, (int32, int32))
        f = ir.Function(builder.function.module, tp_f, 'f')
        g = ir.Function(builder.function.module, tp_g, 'g')
        h = ir.Function(builder.function.module, tp_h, 'h')
        builder.call(f, (a, b), 'res_f')
        builder.call(g, (b, a), 'res_g')
        builder.call(h, (a, b), 'res_h')
        builder.call(f, (a, b), 'res_f_fast', cconv='fastcc')
        res_f_readonly = builder.call(f, (a, b), 'res_f_readonly')
        res_f_readonly.attributes.add('readonly')
        builder.call(f, (a, b), 'res_fast', fastmath='fast')
        builder.call(f, (a, b), 'res_nnan_ninf', fastmath=('nnan', 'ninf'))
        builder.call(f, (a, b), 'res_noinline', attrs='noinline')
        builder.call(f, (a, b), 'res_alwaysinline', attrs='alwaysinline')
        builder.call(f, (a, b), 'res_noinline_ro', attrs=('noinline',
                                                          'readonly'))
        builder.call(f, (a, b), 'res_convergent', attrs='convergent')
        self.check_block(block, """\
        my_block:
            %"res_f" = call float @"f"(i32 %".1", i32 %".2")
            %"res_g" = call double (i32, ...) @"g"(i32 %".2", i32 %".1")
            %"res_h" = call half @"h"(i32 %".1", i32 %".2")
            %"res_f_fast" = call fastcc float @"f"(i32 %".1", i32 %".2")
            %"res_f_readonly" = call float @"f"(i32 %".1", i32 %".2") readonly
            %"res_fast" = call fast float @"f"(i32 %".1", i32 %".2")
            %"res_nnan_ninf" = call ninf nnan float @"f"(i32 %".1", i32 %".2")
            %"res_noinline" = call float @"f"(i32 %".1", i32 %".2") noinline
            %"res_alwaysinline" = call float @"f"(i32 %".1", i32 %".2") alwaysinline
            %"res_noinline_ro" = call float @"f"(i32 %".1", i32 %".2") noinline readonly
            %"res_convergent" = call float @"f"(i32 %".1", i32 %".2") convergent
        """) # noqa E501

    def test_call_metadata(self):
        """
        Function calls with metadata arguments.
        """
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        dbg_declare_ty = ir.FunctionType(ir.VoidType(), [ir.MetaDataType()] * 3)
        dbg_declare = ir.Function(
            builder.module,
            dbg_declare_ty,
            'llvm.dbg.declare')
        a = builder.alloca(int32, name="a")
        b = builder.module.add_metadata(())
        builder.call(dbg_declare, (a, b, b))
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"a" = alloca i32
                    call void @"llvm.dbg.declare"(metadata ptr %"a", metadata !0, metadata !0)
                """)  # noqa E501
        else:
            self.check_block(block, """\
                my_block:
                    %"a" = alloca i32
                    call void @"llvm.dbg.declare"(metadata i32* %"a", metadata !0, metadata !0)
                """)  # noqa E501

    def test_call_attributes(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        fun_ty = ir.FunctionType(
            ir.VoidType(), (int32.as_pointer(), int32, int32.as_pointer()))
        fun = ir.Function(builder.function.module, fun_ty, 'fun')
        fun.args[0].add_attribute('sret')
        retval = builder.alloca(int32, name='retval')
        other = builder.alloca(int32, name='other')
        builder.call(
            fun,
            (retval, ir.Constant(int32, 42), other),
            arg_attrs={
                0: ('sret', 'noalias'),
                2: 'noalias'
            }
        )
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block_regex(block, """\
            my_block:
                %"retval" = alloca i32
                %"other" = alloca i32
                call void @"fun"\\(ptr noalias sret(\\(i32\\))? %"retval", i32 42, ptr noalias %"other"\\)
            """)  # noqa E501
        else:
            self.check_block_regex(block, """\
            my_block:
                %"retval" = alloca i32
                %"other" = alloca i32
                call void @"fun"\\(i32\\* noalias sret(\\(i32\\))? %"retval", i32 42, i32\\* noalias %"other"\\)
            """)  # noqa E501

    def test_call_tail(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        fun_ty = ir.FunctionType(ir.VoidType(), ())
        fun = ir.Function(builder.function.module, fun_ty, 'my_fun')

        builder.call(fun, ())
        builder.call(fun, (), tail=False)
        builder.call(fun, (), tail=True)
        builder.call(fun, (), tail='tail')
        builder.call(fun, (), tail='notail')
        builder.call(fun, (), tail='musttail')
        builder.call(fun, (), tail=[])  # This is a falsy value
        builder.call(fun, (), tail='not a marker')  # This is a truthy value

        self.check_block(block, """\
        my_block:
            call void @"my_fun"()
            call void @"my_fun"()
            tail call void @"my_fun"()
            tail call void @"my_fun"()
            notail call void @"my_fun"()
            musttail call void @"my_fun"()
            call void @"my_fun"()
            tail call void @"my_fun"()
        """)  # noqa E501

    def test_invalid_call_attributes(self):
        block = self.block()
        builder = ir.IRBuilder(block)
        fun_ty = ir.FunctionType(ir.VoidType(), ())
        fun = ir.Function(builder.function.module, fun_ty, 'fun')
        with self.assertRaises(ValueError):
            # The function has no arguments, so this should fail.
            builder.call(fun, (), arg_attrs={0: 'sret'})

    def test_invoke(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        tp_f = ir.FunctionType(flt, (int32, int32))
        f = ir.Function(builder.function.module, tp_f, 'f')
        bb_normal = builder.function.append_basic_block(name='normal')
        bb_unwind = builder.function.append_basic_block(name='unwind')
        builder.invoke(f, (a, b), bb_normal, bb_unwind, 'res_f')
        self.check_block(block, """\
            my_block:
                %"res_f" = invoke float @"f"(i32 %".1", i32 %".2")
                    to label %"normal" unwind label %"unwind"
            """)

    def test_invoke_attributes(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        fun_ty = ir.FunctionType(
            ir.VoidType(), (int32.as_pointer(), int32, int32.as_pointer()))
        fun = ir.Function(builder.function.module, fun_ty, 'fun')
        fun.calling_convention = "fastcc"
        fun.args[0].add_attribute('sret')
        retval = builder.alloca(int32, name='retval')
        other = builder.alloca(int32, name='other')
        bb_normal = builder.function.append_basic_block(name='normal')
        bb_unwind = builder.function.append_basic_block(name='unwind')
        builder.invoke(
            fun,
            (retval, ir.Constant(int32, 42), other),
            bb_normal,
            bb_unwind,
            cconv='fastcc',
            fastmath='fast',
            attrs='noinline',
            arg_attrs={
                0: ('sret', 'noalias'),
                2: 'noalias'
            }
        )
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block_regex(block, """\
            my_block:
                %"retval" = alloca i32
                %"other" = alloca i32
                invoke fast fastcc void @"fun"\\(ptr noalias sret(\\(i32\\))? %"retval", i32 42, ptr noalias %"other"\\) noinline
                    to label %"normal" unwind label %"unwind"
            """)  # noqa E501
        else:
            self.check_block_regex(block, """\
            my_block:
                %"retval" = alloca i32
                %"other" = alloca i32
                invoke fast fastcc void @"fun"\\(i32\\* noalias sret(\\(i32\\))? %"retval", i32 42, i32\\* noalias %"other"\\) noinline
                    to label %"normal" unwind label %"unwind"
            """)  # noqa E501

    def test_landingpad(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        lp = builder.landingpad(ir.LiteralStructType([int32,
                                                      int8.as_pointer()]), 'lp')
        int_typeinfo = ir.GlobalVariable(builder.function.module,
                                         int8.as_pointer(), "_ZTIi")
        int_typeinfo.global_constant = True
        lp.add_clause(ir.CatchClause(int_typeinfo))
        lp.add_clause(ir.FilterClause(ir.Constant(ir.ArrayType(
            int_typeinfo.type, 1), [int_typeinfo])))
        builder.resume(lp)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"lp" = landingpad {i32, ptr}
                        catch ptr @"_ZTIi"
                        filter [1 x ptr] [ptr @"_ZTIi"]
                    resume {i32, ptr} %"lp"
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"lp" = landingpad {i32, i8*}
                        catch i8** @"_ZTIi"
                        filter [1 x i8**] [i8** @"_ZTIi"]
                    resume {i32, i8*} %"lp"
                """)

    def test_assume(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        c = builder.icmp_signed('>', a, b, name='c')
        builder.assume(c)
        self.check_block(block, """\
            my_block:
                %"c" = icmp sgt i32 %".1", %".2"
                call void @"llvm.assume"(i1 %"c")
            """)

    def test_vector_ops(self):
        block = self.block(name='insert_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        a.name = 'a'
        b.name = 'b'

        vecty = ir.VectorType(a.type, 2)
        vec = ir.Constant(vecty, ir.Undefined)
        idxty = ir.IntType(32)
        vec = builder.insert_element(vec, a, idxty(0), name='vec1')
        vec = builder.insert_element(vec, b, idxty(1), name='vec2')

        self.check_block(block, """\
insert_block:
    %"vec1" = insertelement <2 x i32> <i32 undef, i32 undef>, i32 %"a", i32 0
    %"vec2" = insertelement <2 x i32> %"vec1", i32 %"b", i32 1
            """)

        block = builder.append_basic_block("shuffle_block")
        builder.branch(block)
        builder.position_at_end(block)

        mask = ir.Constant(vecty, [1, 0])
        builder.shuffle_vector(vec, vec, mask, name='shuf')

        self.check_block(block, """\
            shuffle_block:
                %"shuf" = shufflevector <2 x i32> %"vec2", <2 x i32> %"vec2", <2 x i32> <i32 1, i32 0>
            """)  # noqa E501

        block = builder.append_basic_block("add_block")
        builder.branch(block)
        builder.position_at_end(block)

        builder.add(vec, vec, name='sum')

        self.check_block(block, """\
            add_block:
                %"sum" = add <2 x i32> %"vec2", %"vec2"
            """)

        block = builder.append_basic_block("extract_block")
        builder.branch(block)
        builder.position_at_end(block)

        c = builder.extract_element(vec, idxty(0), name='ex1')
        d = builder.extract_element(vec, idxty(1), name='ex2')

        self.check_block(block, """\
            extract_block:
              %"ex1" = extractelement <2 x i32> %"vec2", i32 0
              %"ex2" = extractelement <2 x i32> %"vec2", i32 1
            """)

        builder.ret(builder.add(c, d))
        self.assert_valid_ir(builder.module)

    def test_bitreverse(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int64, 5)
        c = builder.bitreverse(a, name='c')
        builder.ret(c)
        self.check_block(block, """\
            my_block:
                %"c" = call i64 @"llvm.bitreverse.i64"(i64 5)
                ret i64 %"c"
            """)

    def test_bitreverse_wrongtype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5)

        with self.assertRaises(TypeError) as raises:
            builder.bitreverse(a, name='c')
        self.assertIn(
            "expected an integer type, got float",
            str(raises.exception))

    def test_fence(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        with self.assertRaises(ValueError) as raises:
            builder.fence("monotonic", None)
        self.assertIn(
            "Invalid fence ordering \"monotonic\"!",
            str(raises.exception))
        with self.assertRaises(ValueError) as raises:
            builder.fence(None, "monotonic")
        self.assertIn(
            "Invalid fence ordering \"None\"!",
            str(raises.exception))
        builder.fence("acquire", None)
        builder.fence("release", "singlethread")
        builder.fence("acq_rel", "singlethread")
        builder.fence("seq_cst")
        builder.ret_void()
        self.check_block(block, """\
            my_block:
                fence acquire
                fence syncscope("singlethread") release
                fence syncscope("singlethread") acq_rel
                fence seq_cst
                ret void
            """)

    def test_comment(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        with self.assertRaises(AssertionError):
            builder.comment("so\nmany lines")
        builder.comment("my comment")
        builder.ret_void()
        self.check_block(block, """\
            my_block:
                ; my comment
                ret void
            """)

    def test_bswap(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int32, 5)
        c = builder.bswap(a, name='c')
        builder.ret(c)
        self.check_block(block, """\
            my_block:
                %"c" = call i32 @"llvm.bswap.i32"(i32 5)
                ret i32 %"c"
            """)

    def test_ctpop(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int16, 5)
        c = builder.ctpop(a, name='c')
        builder.ret(c)
        self.check_block(block, """\
            my_block:
                %"c" = call i16 @"llvm.ctpop.i16"(i16 5)
                ret i16 %"c"
            """)

    def test_ctlz(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int16, 5)
        b = ir.Constant(int1, 1)
        c = builder.ctlz(a, b, name='c')
        builder.ret(c)
        self.check_block(block, """\
            my_block:
                %"c" = call i16 @"llvm.ctlz.i16"(i16 5, i1 1)
                ret i16 %"c"
            """)

    def test_convert_to_fp16_f32(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5.0)
        b = builder.convert_to_fp16(a, name='b')
        builder.ret(b)
        self.check_block(block, """\
            my_block:
                %"b" = call i16 @"llvm.convert.to.fp16.f32"(float 0x4014000000000000)
                ret i16 %"b"
            """)  # noqa E501

    def test_convert_to_fp16_f32_wrongtype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int16, 5)

        with self.assertRaises(TypeError) as raises:
            builder.convert_to_fp16(a, name='b')
        self.assertIn(
            "expected a float type, got i16",
            str(raises.exception))

    def test_convert_from_fp16_f32(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int16, 5)
        b = builder.convert_from_fp16(a, name='b', to=flt)
        builder.ret(b)
        self.check_block(block, """\
            my_block:
                %"b" = call float @"llvm.convert.from.fp16.f32"(i16 5)
                ret float %"b"
            """)

    def test_convert_from_fp16_f32_notype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5.5)

        with self.assertRaises(TypeError) as raises:
            builder.convert_from_fp16(a, name='b')
        self.assertIn(
            "expected a float return type",
            str(raises.exception))

    def test_convert_from_fp16_f32_wrongtype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5.5)

        with self.assertRaises(TypeError) as raises:
            builder.convert_from_fp16(a, name='b', to=flt)
        self.assertIn(
            "expected an i16 type, got float",
            str(raises.exception))

    def test_convert_from_fp16_f32_wrongtype2(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5.5)

        with self.assertRaises(TypeError) as raises:
            builder.convert_from_fp16(a, name='b', to=int16)
        self.assertIn(
            "expected a float type, got i16",
            str(raises.exception))

    def test_cttz(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int64, 5)
        b = ir.Constant(int1, 1)
        c = builder.cttz(a, b, name='c')
        builder.ret(c)
        self.check_block(block, """\
            my_block:
                %"c" = call i64 @"llvm.cttz.i64"(i64 5, i1 1)
                ret i64 %"c"
            """)

    def test_cttz_wrongflag(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int64, 5)
        b = ir.Constant(int32, 3)

        with self.assertRaises(TypeError) as raises:
            builder.cttz(a, b, name='c')
        self.assertIn(
            "expected an i1 type, got i32",
            str(raises.exception))

    def test_cttz_wrongtype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5)
        b = ir.Constant(int1, 1)

        with self.assertRaises(TypeError) as raises:
            builder.cttz(a, b, name='c')
        self.assertIn(
            "expected an integer type, got float",
            str(raises.exception))

    def test_fma(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5)
        b = ir.Constant(flt, 1)
        c = ir.Constant(flt, 2)
        fma = builder.fma(a, b, c, name='fma')
        builder.ret(fma)
        self.check_block(block, """\
            my_block:
                %"fma" = call float @"llvm.fma.f32"(float 0x4014000000000000, float 0x3ff0000000000000, float 0x4000000000000000)
                ret float %"fma"
            """)  # noqa E501

    def test_fma_wrongtype(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(int32, 5)
        b = ir.Constant(int32, 1)
        c = ir.Constant(int32, 2)

        with self.assertRaises(TypeError) as raises:
            builder.fma(a, b, c, name='fma')
        self.assertIn(
            "expected an floating point type, got i32",
            str(raises.exception))

    def test_fma_mixedtypes(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a = ir.Constant(flt, 5)
        b = ir.Constant(dbl, 1)
        c = ir.Constant(flt, 2)

        with self.assertRaises(TypeError) as raises:
            builder.fma(a, b, c, name='fma')
        self.assertIn(
            "expected types to be the same, got float, double, float",
            str(raises.exception))

    def test_arg_attributes(self):
        def gen_code(attr_name):
            fnty = ir.FunctionType(ir.IntType(32), [ir.IntType(32).as_pointer(),
                                                    ir.IntType(32)])
            module = ir.Module()

            func = ir.Function(module, fnty, name="sum")

            bb_entry = func.append_basic_block()
            bb_loop = func.append_basic_block()
            bb_exit = func.append_basic_block()

            builder = ir.IRBuilder()
            builder.position_at_end(bb_entry)

            builder.branch(bb_loop)
            builder.position_at_end(bb_loop)

            index = builder.phi(ir.IntType(32))
            index.add_incoming(ir.Constant(index.type, 0), bb_entry)
            accum = builder.phi(ir.IntType(32))
            accum.add_incoming(ir.Constant(accum.type, 0), bb_entry)

            func.args[0].add_attribute(attr_name)
            ptr = builder.gep(func.args[0], [index])
            value = builder.load(ptr)

            added = builder.add(accum, value)
            accum.add_incoming(added, bb_loop)

            indexp1 = builder.add(index, ir.Constant(index.type, 1))
            index.add_incoming(indexp1, bb_loop)

            cond = builder.icmp_unsigned('<', indexp1, func.args[1])
            builder.cbranch(cond, bb_loop, bb_exit)

            builder.position_at_end(bb_exit)
            builder.ret(added)

            return str(module)

        for attr_name in (
            'byref',
            'byval',
            'elementtype',
            'immarg',
            'inalloca',
            'inreg',
            'nest',
            'noalias',
            'nocapture',
            'nofree',
            'nonnull',
            'noundef',
            'preallocated',
            'returned',
            'signext',
            'swiftasync',
            'swifterror',
            'swiftself',
            'zeroext',
        ):
            # If this parses, we emitted the right byval attribute format
            llvm.parse_assembly(gen_code(attr_name))
        # sret doesn't fit this pattern and is tested in test_call_attributes


class TestBuilderMisc(TestBase):
    """
    Test various other features of the IRBuilder class.
    """

    def test_attributes(self):
        block = self.block(name='start')
        builder = ir.IRBuilder(block)
        self.assertIs(builder.function, block.parent)
        self.assertIsInstance(builder.function, ir.Function)
        self.assertIs(builder.module, block.parent.module)
        self.assertIsInstance(builder.module, ir.Module)

    def test_goto_block(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        a, b = builder.function.args[:2]
        builder.add(a, b, 'c')
        bb_new = builder.append_basic_block(name='foo')
        with builder.goto_block(bb_new):
            builder.fadd(a, b, 'd')
            with builder.goto_entry_block():
                builder.sub(a, b, 'e')
            builder.fsub(a, b, 'f')
            builder.branch(bb_new)
        builder.mul(a, b, 'g')
        with builder.goto_block(bb_new):
            builder.fmul(a, b, 'h')
        self.check_block(block, """\
            my_block:
                %"c" = add i32 %".1", %".2"
                %"e" = sub i32 %".1", %".2"
                %"g" = mul i32 %".1", %".2"
            """)
        self.check_block(bb_new, """\
            foo:
                %"d" = fadd i32 %".1", %".2"
                %"f" = fsub i32 %".1", %".2"
                %"h" = fmul i32 %".1", %".2"
                br label %"foo"
            """)

    def test_if_then(self):
        block = self.block(name='one')
        builder = ir.IRBuilder(block)
        z = ir.Constant(int1, 0)
        a = builder.add(z, z, 'a')
        with builder.if_then(a) as bbend:
            builder.add(z, z, 'b')
            # Block will be terminated implicitly
        self.assertIs(builder.block, bbend)
        c = builder.add(z, z, 'c')
        with builder.if_then(c):
            builder.add(z, z, 'd')
            builder.branch(block)
            # No implicit termination
        self.check_func_body(builder.function, """\
            one:
                %"a" = add i1 0, 0
                br i1 %"a", label %"one.if", label %"one.endif"
            one.if:
                %"b" = add i1 0, 0
                br label %"one.endif"
            one.endif:
                %"c" = add i1 0, 0
                br i1 %"c", label %"one.endif.if", label %"one.endif.endif"
            one.endif.if:
                %"d" = add i1 0, 0
                br label %"one"
            one.endif.endif:
            """)

    def test_if_then_nested(self):
        # Implicit termination in a nested if/then
        block = self.block(name='one')
        builder = ir.IRBuilder(block)
        z = ir.Constant(int1, 0)
        a = builder.add(z, z, 'a')
        with builder.if_then(a):
            b = builder.add(z, z, 'b')
            with builder.if_then(b):
                builder.add(z, z, 'c')
        builder.ret_void()
        self.check_func_body(builder.function, """\
            one:
                %"a" = add i1 0, 0
                br i1 %"a", label %"one.if", label %"one.endif"
            one.if:
                %"b" = add i1 0, 0
                br i1 %"b", label %"one.if.if", label %"one.if.endif"
            one.endif:
                ret void
            one.if.if:
                %"c" = add i1 0, 0
                br label %"one.if.endif"
            one.if.endif:
                br label %"one.endif"
            """)

    def test_if_then_long_label(self):
        full_label = 'Long' * 20
        block = self.block(name=full_label)
        builder = ir.IRBuilder(block)
        z = ir.Constant(int1, 0)
        a = builder.add(z, z, 'a')
        with builder.if_then(a):
            b = builder.add(z, z, 'b')
            with builder.if_then(b):
                builder.add(z, z, 'c')
        builder.ret_void()
        self.check_func_body(builder.function, """\
            {full_label}:
                %"a" = add i1 0, 0
                br i1 %"a", label %"{label}.if", label %"{label}.endif"
            {label}.if:
                %"b" = add i1 0, 0
                br i1 %"b", label %"{label}.if.if", label %"{label}.if.endif"
            {label}.endif:
                ret void
            {label}.if.if:
                %"c" = add i1 0, 0
                br label %"{label}.if.endif"
            {label}.if.endif:
                br label %"{label}.endif"
            """.format(full_label=full_label, label=full_label[:25] + '..'))

    def test_if_then_likely(self):
        def check(likely):
            block = self.block(name='one')
            builder = ir.IRBuilder(block)
            z = ir.Constant(int1, 0)
            with builder.if_then(z, likely=likely):
                pass
            self.check_block(block, """\
                one:
                    br i1 0, label %"one.if", label %"one.endif", !prof !0
                """)
            return builder
        builder = check(True)
        self.check_metadata(builder.module, """\
            !0 = !{ !"branch_weights", i32 99, i32 1 }
            """)
        builder = check(False)
        self.check_metadata(builder.module, """\
            !0 = !{ !"branch_weights", i32 1, i32 99 }
            """)

    def test_if_else(self):
        block = self.block(name='one')
        builder = ir.IRBuilder(block)
        z = ir.Constant(int1, 0)
        a = builder.add(z, z, 'a')
        with builder.if_else(a) as (then, otherwise):
            with then:
                builder.add(z, z, 'b')
            with otherwise:
                builder.add(z, z, 'c')
            # Each block will be terminated implicitly
        with builder.if_else(a) as (then, otherwise):
            with then:
                builder.branch(block)
            with otherwise:
                builder.ret_void()
            # No implicit termination
        self.check_func_body(builder.function, """\
            one:
                %"a" = add i1 0, 0
                br i1 %"a", label %"one.if", label %"one.else"
            one.if:
                %"b" = add i1 0, 0
                br label %"one.endif"
            one.else:
                %"c" = add i1 0, 0
                br label %"one.endif"
            one.endif:
                br i1 %"a", label %"one.endif.if", label %"one.endif.else"
            one.endif.if:
                br label %"one"
            one.endif.else:
                ret void
            one.endif.endif:
            """)

    def test_if_else_likely(self):
        def check(likely):
            block = self.block(name='one')
            builder = ir.IRBuilder(block)
            z = ir.Constant(int1, 0)
            with builder.if_else(z, likely=likely) as (then, otherwise):
                with then:
                    builder.branch(block)
                with otherwise:
                    builder.ret_void()
            self.check_func_body(builder.function, """\
                one:
                    br i1 0, label %"one.if", label %"one.else", !prof !0
                one.if:
                    br label %"one"
                one.else:
                    ret void
                one.endif:
                """)
            return builder
        builder = check(True)
        self.check_metadata(builder.module, """\
            !0 = !{ !"branch_weights", i32 99, i32 1 }
            """)
        builder = check(False)
        self.check_metadata(builder.module, """\
            !0 = !{ !"branch_weights", i32 1, i32 99 }
            """)

    def test_positioning(self):
        """
        Test IRBuilder.position_{before,after,at_start,at_end}.
        """
        func = self.function()
        builder = ir.IRBuilder()
        z = ir.Constant(int32, 0)
        bb_one = func.append_basic_block(name='one')
        bb_two = func.append_basic_block(name='two')
        bb_three = func.append_basic_block(name='three')
        # .at_start(empty block)
        builder.position_at_start(bb_one)
        builder.add(z, z, 'a')
        # .at_end(empty block)
        builder.position_at_end(bb_two)
        builder.add(z, z, 'm')
        builder.add(z, z, 'n')
        # .at_start(block)
        builder.position_at_start(bb_two)
        o = builder.add(z, z, 'o')
        builder.add(z, z, 'p')
        # .at_end(block)
        builder.position_at_end(bb_one)
        b = builder.add(z, z, 'b')
        # .after(instr)
        builder.position_after(o)
        builder.add(z, z, 'q')
        # .before(instr)
        builder.position_before(b)
        builder.add(z, z, 'c')
        self.check_block(bb_one, """\
            one:
                %"a" = add i32 0, 0
                %"c" = add i32 0, 0
                %"b" = add i32 0, 0
            """)
        self.check_block(bb_two, """\
            two:
                %"o" = add i32 0, 0
                %"q" = add i32 0, 0
                %"p" = add i32 0, 0
                %"m" = add i32 0, 0
                %"n" = add i32 0, 0
            """)
        self.check_block(bb_three, """\
            three:
            """)

    def test_instruction_removal(self):
        func = self.function()
        builder = ir.IRBuilder()
        blk = func.append_basic_block(name='entry')
        builder.position_at_end(blk)
        k = ir.Constant(int32, 1234)
        a = builder.add(k, k, 'a')
        retvoid = builder.ret_void()
        self.assertTrue(blk.is_terminated)
        builder.remove(retvoid)
        self.assertFalse(blk.is_terminated)
        b = builder.mul(a, a, 'b')
        c = builder.add(b, b, 'c')
        builder.remove(c)
        builder.ret_void()
        self.assertTrue(blk.is_terminated)
        self.check_block(blk, """\
            entry:
                %"a" = add i32 1234, 1234
                %"b" = mul i32 %"a", %"a"
                ret void
        """)

    def test_metadata(self):
        block = self.block(name='my_block')
        builder = ir.IRBuilder(block)
        builder.debug_metadata = builder.module.add_metadata([])
        builder.alloca(ir.PointerType(int32), name='c')
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca ptr, !dbg !0
                """)
        else:
            self.check_block(block, """\
                my_block:
                    %"c" = alloca i32*, !dbg !0
                """)


class TestTypes(TestBase):

    def has_logical_equality(self, ty):
        while isinstance(ty, ir.PointerType):
            ty = ty.pointee
        return not isinstance(ty, ir.LabelType)

    def assorted_types(self):
        """
        A bunch of mutually unequal types
        """
        # Avoid polluting the namespace
        context = ir.Context()
        types = [
            ir.LabelType(), ir.VoidType(),
            ir.FunctionType(int1, (int8, int8)), ir.FunctionType(int1, (int8,)),
            ir.FunctionType(int1, (int8,), var_arg=True),
            ir.FunctionType(int8, (int8,)),
            int1, int8, int32, flt, dbl,
            ir.ArrayType(flt, 5), ir.ArrayType(dbl, 5), ir.ArrayType(dbl, 4),
            ir.LiteralStructType((int1, int8)), ir.LiteralStructType((int8,
                                                                      int1)),
            context.get_identified_type("MyType1"),
            context.get_identified_type("MyType2"),
        ]
        types += [ir.PointerType(tp) for tp in types
                  if not isinstance(tp, (ir.VoidType, ir.LabelType))]

        return types

    def test_pickling(self):
        types = self.assorted_types()
        for ty in types:
            newty = self.assert_pickle_correctly(ty)
            if self.has_logical_equality(ty):
                self.assertEqual(newty, ty)

    def test_comparisons(self):
        types = self.assorted_types()
        for a, b in itertools.product(types, types):
            if a is not b:
                self.assertFalse(a == b, (a, b))
                self.assertTrue(a != b, (a, b))
        # We assume copy.copy() works fine here...
        for tp in types:
            other = copy.copy(tp)
            if self.has_logical_equality(tp):
                self.assertTrue(tp == other, (tp, other))
                self.assertFalse(tp != other, (tp, other))
            else:
                self.assertFalse(tp == other, (tp, other))
                self.assertTrue(tp != other, (tp, other))

    def test_str(self):
        """
        Test the string representation of types.
        """
        self.assertEqual(str(int1), 'i1')
        self.assertEqual(str(ir.IntType(29)), 'i29')
        self.assertEqual(str(flt), 'float')
        self.assertEqual(str(dbl), 'double')
        self.assertEqual(str(ir.VoidType()), 'void')
        self.assertEqual(str(ir.FunctionType(int1, ())), 'i1 ()')
        self.assertEqual(str(ir.FunctionType(int1, (flt,))), 'i1 (float)')
        self.assertEqual(str(ir.FunctionType(int1, (flt, dbl))),
                         'i1 (float, double)')
        self.assertEqual(str(ir.FunctionType(int1, (), var_arg=True)),
                         'i1 (...)')
        self.assertEqual(str(ir.FunctionType(int1, (flt,), var_arg=True)),
                         'i1 (float, ...)')
        self.assertEqual(str(ir.FunctionType(int1, (flt, dbl), var_arg=True)),
                         'i1 (float, double, ...)')
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(ir.PointerType(int32)), 'ptr')
            self.assertEqual(str(ir.PointerType(ir.PointerType(int32))), 'ptr')
        else:
            self.assertEqual(str(ir.PointerType(int32)), 'i32*')
            self.assertEqual(str(ir.PointerType(ir.PointerType(int32))),
                             'i32**')
        self.assertEqual(str(ir.ArrayType(int1, 5)), '[5 x i1]')
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(ir.ArrayType(ir.PointerType(int1), 5)),
                             '[5 x ptr]')
            self.assertEqual(str(ir.PointerType(ir.ArrayType(int1, 5))), 'ptr')
        else:
            self.assertEqual(str(ir.ArrayType(ir.PointerType(int1), 5)),
                             '[5 x i1*]')
            self.assertEqual(str(ir.PointerType(ir.ArrayType(int1, 5))),
                             '[5 x i1]*')
        self.assertEqual(str(ir.LiteralStructType((int1,))), '{i1}')
        self.assertEqual(str(ir.LiteralStructType((int1, flt))), '{i1, float}')
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(ir.LiteralStructType((
                ir.PointerType(int1), ir.LiteralStructType((int32, int8))))),
                '{ptr, {i32, i8}}')
        else:
            self.assertEqual(str(ir.LiteralStructType((
                ir.PointerType(int1), ir.LiteralStructType((int32, int8))))),
                '{i1*, {i32, i8}}')
        self.assertEqual(str(ir.LiteralStructType((int1,), packed=True)),
                         '<{i1}>')
        self.assertEqual(str(ir.LiteralStructType((int1, flt), packed=True)),
                         '<{i1, float}>')

        # Avoid polluting the namespace
        context = ir.Context()
        mytype = context.get_identified_type("MyType")
        self.assertEqual(str(mytype), "%\"MyType\"")
        mytype1 = context.get_identified_type("MyType\\")
        self.assertEqual(str(mytype1), "%\"MyType\\5c\"")
        mytype2 = context.get_identified_type("MyType\"")
        self.assertEqual(str(mytype2), "%\"MyType\\22\"")

    def test_hash(self):
        for typ in filter(self.has_logical_equality, self.assorted_types()):
            self.assertEqual(hash(typ), hash(copy.copy(typ)))

    def test_gep(self):
        def check_constant(tp, i, expected):
            actual = tp.gep(ir.Constant(int32, i))
            self.assertEqual(actual, expected)

        def check_index_type(tp):
            index = ir.Constant(dbl, 1.0)
            with self.assertRaises(TypeError):
                tp.gep(index)

        tp = ir.PointerType(dbl)
        for i in range(5):
            check_constant(tp, i, dbl)
        check_index_type(tp)

        tp = ir.ArrayType(int1, 3)
        for i in range(3):
            check_constant(tp, i, int1)
        check_index_type(tp)

        tp = ir.LiteralStructType((dbl, ir.LiteralStructType((int1, int8))))
        check_constant(tp, 0, dbl)
        check_constant(tp, 1, ir.LiteralStructType((int1, int8)))
        with self.assertRaises(IndexError):
            tp.gep(ir.Constant(int32, 2))
        check_index_type(tp)

        context = ir.Context()
        tp = ir.IdentifiedStructType(context, "MyType")
        tp.set_body(dbl, ir.LiteralStructType((int1, int8)))
        check_constant(tp, 0, dbl)
        check_constant(tp, 1, ir.LiteralStructType((int1, int8)))
        with self.assertRaises(IndexError):
            tp.gep(ir.Constant(int32, 2))
        check_index_type(tp)

    def test_abi_size(self):
        td = llvm.create_target_data("e-m:e-i64:64-f80:128-n8:16:32:64-S128")

        def check(tp, expected):
            self.assertEqual(tp.get_abi_size(td), expected)
        check(int8, 1)
        check(int32, 4)
        check(int64, 8)
        check(ir.ArrayType(int8, 5), 5)
        check(ir.ArrayType(int32, 5), 20)
        check(ir.LiteralStructType((dbl, flt, flt)), 16)

    def test_abi_alignment(self):
        td = llvm.create_target_data("e-m:e-i64:64-f80:128-n8:16:32:64-S128")

        def check(tp, expected):
            self.assertIn(tp.get_abi_alignment(td), expected)
        check(int8, (1, 2, 4))
        check(int32, (4,))
        check(int64, (8,))
        check(ir.ArrayType(int8, 5), (1, 2, 4))
        check(ir.ArrayType(int32, 5), (4,))
        check(ir.LiteralStructType((dbl, flt, flt)), (8,))

    def test_identified_struct(self):
        context = ir.Context()
        mytype = context.get_identified_type("MyType")
        module = ir.Module(context=context)
        self.assertTrue(mytype.is_opaque)
        self.assert_valid_ir(module)
        oldstr = str(module)
        mytype.set_body(ir.IntType(32), ir.IntType(64), ir.FloatType())
        self.assertFalse(mytype.is_opaque)
        self.assert_valid_ir(module)
        self.assertNotEqual(oldstr, str(module))

    def test_identified_struct_packed(self):
        td = llvm.create_target_data("e-m:e-i64:64-f80:128-n8:16:32:64-S128")
        context = ir.Context()
        mytype = context.get_identified_type("MyType", True)
        module = ir.Module(context=context)
        self.assertTrue(mytype.is_opaque)
        self.assert_valid_ir(module)
        oldstr = str(module)
        mytype.set_body(ir.IntType(16), ir.IntType(64), ir.FloatType())
        self.assertEqual(mytype.get_element_offset(td, 1, context), 2)
        self.assertFalse(mytype.is_opaque)
        self.assert_valid_ir(module)
        self.assertNotEqual(oldstr, str(module))

    def test_target_data_non_default_context(self):
        context = ir.Context()
        mytype = context.get_identified_type("MyType")
        mytype.elements = [ir.IntType(32)]
        td = llvm.create_target_data("e-m:e-i64:64-f80:128-n8:16:32:64-S128")
        self.assertEqual(mytype.get_abi_size(td, context=context), 4)

    def test_vector(self):
        vecty = ir.VectorType(ir.IntType(32), 8)
        self.assertEqual(str(vecty), "<8 x i32>")


def c32(i):
    return ir.Constant(int32, i)


class TestConstant(TestBase):

    def test_integers(self):
        c = ir.Constant(int32, 42)
        self.assertEqual(str(c), 'i32 42')
        c = ir.Constant(int1, 1)
        self.assertEqual(str(c), 'i1 1')
        c = ir.Constant(int1, 0)
        self.assertEqual(str(c), 'i1 0')
        c = ir.Constant(int1, True)
        self.assertEqual(str(c), 'i1 true')
        c = ir.Constant(int1, False)
        self.assertEqual(str(c), 'i1 false')
        c = ir.Constant(int1, ir.Undefined)
        self.assertEqual(str(c), 'i1 undef')
        c = ir.Constant(int1, None)
        self.assertEqual(str(c), 'i1 0')

    def test_reals(self):
        # XXX Test NaNs and infs
        c = ir.Constant(flt, 1.5)
        self.assertEqual(str(c), 'float 0x3ff8000000000000')
        c = ir.Constant(flt, -1.5)
        self.assertEqual(str(c), 'float 0xbff8000000000000')
        c = ir.Constant(dbl, 1.5)
        self.assertEqual(str(c), 'double 0x3ff8000000000000')
        c = ir.Constant(dbl, -1.5)
        self.assertEqual(str(c), 'double 0xbff8000000000000')
        c = ir.Constant(dbl, ir.Undefined)
        self.assertEqual(str(c), 'double undef')
        c = ir.Constant(dbl, None)
        self.assertEqual(str(c), 'double 0.0')

    def test_arrays(self):
        c = ir.Constant(ir.ArrayType(int32, 3), (c32(5), c32(6), c32(4)))
        self.assertEqual(str(c), '[3 x i32] [i32 5, i32 6, i32 4]')
        c = ir.Constant(ir.ArrayType(int32, 2), (c32(5), c32(ir.Undefined)))
        self.assertEqual(str(c), '[2 x i32] [i32 5, i32 undef]')

        c = ir.Constant.literal_array((c32(5), c32(6), c32(ir.Undefined)))
        self.assertEqual(str(c), '[3 x i32] [i32 5, i32 6, i32 undef]')
        with self.assertRaises(TypeError) as raises:
            ir.Constant.literal_array((c32(5), ir.Constant(flt, 1.5)))
        self.assertEqual(str(raises.exception),
                         "all elements must have the same type")

        c = ir.Constant(ir.ArrayType(int32, 2), ir.Undefined)
        self.assertEqual(str(c), '[2 x i32] undef')
        c = ir.Constant(ir.ArrayType(int32, 2), None)
        self.assertEqual(str(c), '[2 x i32] zeroinitializer')
        # Raw array syntax
        c = ir.Constant(ir.ArrayType(int8, 11), bytearray(b"foobar_123\x80"))
        self.assertEqual(str(c), r'[11 x i8] c"foobar_123\80"')
        c = ir.Constant(ir.ArrayType(int8, 4), bytearray(b"\x00\x01\x04\xff"))
        self.assertEqual(str(c), r'[4 x i8] c"\00\01\04\ff"')
        # Recursive instantiation of inner constants
        c = ir.Constant(ir.ArrayType(int32, 3), (5, ir.Undefined, 6))
        self.assertEqual(str(c), '[3 x i32] [i32 5, i32 undef, i32 6]')
        # Invalid number of args
        with self.assertRaises(ValueError):
            ir.Constant(ir.ArrayType(int32, 3), (5, 6))

    def test_vector(self):
        vecty = ir.VectorType(ir.IntType(32), 8)
        vals = [1, 2, 4, 3, 8, 6, 9, 7]
        vec = ir.Constant(vecty, vals)
        vec_repr = "<8 x i32> <{}>".format(
            ', '.join(map('i32 {}'.format, vals)))
        self.assertEqual(str(vec), vec_repr)

    def test_non_nullable_int(self):
        constant = ir.Constant(ir.IntType(32), None).constant
        self.assertEqual(constant, 0)

    def test_structs(self):
        st1 = ir.LiteralStructType((flt, int1))
        st2 = ir.LiteralStructType((int32, st1))
        c = ir.Constant(st1, (ir.Constant(ir.FloatType(), 1.5),
                              ir.Constant(int1, True)))
        self.assertEqual(str(c),
                         '{float, i1} {float 0x3ff8000000000000, i1 true}')
        c = ir.Constant.literal_struct((ir.Constant(ir.FloatType(), 1.5),
                                        ir.Constant(int1, True)))
        self.assertEqual(c.type, st1)
        self.assertEqual(str(c),
                         '{float, i1} {float 0x3ff8000000000000, i1 true}')
        c = ir.Constant.literal_struct((ir.Constant(ir.FloatType(), 1.5),
                                        ir.Constant(int1, ir.Undefined)))
        self.assertEqual(c.type, st1)
        self.assertEqual(str(c),
                         '{float, i1} {float 0x3ff8000000000000, i1 undef}')
        c = ir.Constant(st1, ir.Undefined)
        self.assertEqual(str(c), '{float, i1} undef')
        c = ir.Constant(st1, None)
        self.assertEqual(str(c), '{float, i1} zeroinitializer')
        # Recursive instantiation of inner constants
        c1 = ir.Constant(st1, (1.5, True))
        self.assertEqual(str(c1),
                         '{float, i1} {float 0x3ff8000000000000, i1 true}')
        c2 = ir.Constant(st2, (42, c1))
        self.assertEqual(str(c2), ('{i32, {float, i1}} {i32 42, {float, i1} '
                                   '{float 0x3ff8000000000000, i1 true}}'))
        c3 = ir.Constant(st2, (42, (1.5, True)))
        self.assertEqual(str(c3), str(c2))
        # Invalid number of args
        with self.assertRaises(ValueError):
            ir.Constant(st2, (4, 5, 6))

    def test_undefined_literal_struct_pickling(self):
        i8 = ir.IntType(8)
        st = ir.Constant(ir.LiteralStructType([i8, i8]), ir.Undefined)
        self.assert_pickle_correctly(st)

    def test_type_instantiaton(self):
        """
        Instantiating a type should create a constant.
        """
        c = int8(42)
        self.assertIsInstance(c, ir.Constant)
        self.assertEqual(str(c), 'i8 42')
        c = int1(True)
        self.assertIsInstance(c, ir.Constant)
        self.assertEqual(str(c), 'i1 true')
        # Arrays
        at = ir.ArrayType(int32, 3)
        c = at([c32(4), c32(5), c32(6)])
        self.assertEqual(str(c), '[3 x i32] [i32 4, i32 5, i32 6]')
        c = at([4, 5, 6])
        self.assertEqual(str(c), '[3 x i32] [i32 4, i32 5, i32 6]')
        c = at(None)
        self.assertEqual(str(c), '[3 x i32] zeroinitializer')
        with self.assertRaises(ValueError):
            at([4, 5, 6, 7])
        # Structs
        st1 = ir.LiteralStructType((flt, int1))
        st2 = ir.LiteralStructType((int32, st1))
        c = st1((1.5, True))
        self.assertEqual(str(c), ('{float, i1} {float 0x3ff8000000000000, i1 '
                                  'true}'))
        c = st2((42, (1.5, True)))
        self.assertEqual(str(c), ('{i32, {float, i1}} {i32 42, {float, i1} '
                                  '{float 0x3ff8000000000000, i1 true}}'))

    def test_repr(self):
        """
        Constants should have a useful repr().
        """
        c = int32(42)
        self.assertEqual(repr(c), "<ir.Constant type='i32' value=42>")

    def test_encoding_problem(self):
        c = ir.Constant(ir.ArrayType(ir.IntType(8), 256),
                        bytearray(range(256)))
        m = self.module()
        gv = ir.GlobalVariable(m, c.type, "myconstant")
        gv.global_constant = True
        gv.initializer = c
        # With utf-8, the following will cause:
        # UnicodeDecodeError: 'utf-8' codec can't decode byte 0xe0 in position
        # 136: invalid continuation byte
        parsed = llvm.parse_assembly(str(m))
        # Make sure the encoding does not modify the IR
        reparsed = llvm.parse_assembly(str(parsed))
        self.assertEqual(str(parsed), str(reparsed))

    def test_gep(self):
        m = self.module()
        tp = ir.LiteralStructType((flt, int1))
        gv = ir.GlobalVariable(m, tp, "myconstant")
        c = gv.gep([ir.Constant(int32, x) for x in (0, 1)])
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c),
                             'getelementptr ({float, i1}, ptr @"myconstant", i32 0, i32 1)')  # noqa E501
        else:
            self.assertEqual(str(c),
                             'getelementptr ({float, i1}, {float, i1}* @"myconstant", i32 0, i32 1)')  # noqa E501
        self.assertEqual(c.type, ir.PointerType(int1))

        const = ir.Constant(tp, None)
        with self.assertRaises(TypeError):
            const.gep([ir.Constant(int32, 0)])

        const_ptr = ir.Constant(tp.as_pointer(), None)
        c2 = const_ptr.gep([ir.Constant(int32, 0)])
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c2),
                             'getelementptr ({float, i1}, ptr null, i32 0)')  # noqa E501
        else:
            self.assertEqual(str(c2),
                             'getelementptr ({float, i1}, {float, i1}* null, i32 0)')  # noqa E501
        self.assertEqual(c.type, ir.PointerType(int1))

    def test_gep_addrspace_globalvar(self):
        m = self.module()
        tp = ir.LiteralStructType((flt, int1))
        addrspace = 4

        gv = ir.GlobalVariable(m, tp, "myconstant", addrspace=addrspace)
        self.assertEqual(gv.addrspace, addrspace)
        c = gv.gep([ir.Constant(int32, x) for x in (0, 1)])
        self.assertEqual(c.type.addrspace, addrspace)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c),
                             ('getelementptr ({float, i1}, ptr '
                              'addrspace(4) @"myconstant", i32 0, i32 1)'))
        else:
            self.assertEqual(str(c),
                             ('getelementptr ({float, i1}, {float, i1} '
                              'addrspace(4)* @"myconstant", i32 0, i32 1)'))
        self.assertEqual(c.type, ir.PointerType(int1, addrspace=addrspace))

    def test_trunc(self):
        c = ir.Constant(int64, 1).trunc(int32)
        self.assertEqual(str(c), 'trunc (i64 1 to i32)')

    def test_zext(self):
        c = ir.Constant(int32, 1).zext(int64)
        self.assertEqual(str(c), 'zext (i32 1 to i64)')

    def test_sext(self):
        c = ir.Constant(int32, -1).sext(int64)
        self.assertEqual(str(c), 'sext (i32 -1 to i64)')

    def test_fptrunc(self):
        c = ir.Constant(flt, 1).fptrunc(hlf)
        self.assertEqual(str(c), 'fptrunc (float 0x3ff0000000000000 to half)')

    def test_fpext(self):
        c = ir.Constant(flt, 1).fpext(dbl)
        self.assertEqual(str(c), 'fpext (float 0x3ff0000000000000 to double)')

    def test_bitcast(self):
        m = self.module()
        gv = ir.GlobalVariable(m, int32, "myconstant")
        c = gv.bitcast(int64.as_pointer())
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c), 'bitcast (ptr @"myconstant" to ptr)')
        else:
            self.assertEqual(str(c), 'bitcast (i32* @"myconstant" to i64*)')

    def test_fptoui(self):
        c = ir.Constant(flt, 1).fptoui(int32)
        self.assertEqual(str(c), 'fptoui (float 0x3ff0000000000000 to i32)')

    def test_uitofp(self):
        c = ir.Constant(int32, 1).uitofp(flt)
        self.assertEqual(str(c), 'uitofp (i32 1 to float)')

    def test_fptosi(self):
        c = ir.Constant(flt, 1).fptosi(int32)
        self.assertEqual(str(c), 'fptosi (float 0x3ff0000000000000 to i32)')

    def test_sitofp(self):
        c = ir.Constant(int32, 1).sitofp(flt)
        self.assertEqual(str(c), 'sitofp (i32 1 to float)')

    def test_ptrtoint_1(self):
        ptr = ir.Constant(int64.as_pointer(), None)
        one = ir.Constant(int32, 1)
        c = ptr.ptrtoint(int32)

        self.assertRaises(TypeError, one.ptrtoint, int64)
        self.assertRaises(TypeError, ptr.ptrtoint, flt)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c), 'ptrtoint (ptr null to i32)')
        else:
            self.assertEqual(str(c), 'ptrtoint (i64* null to i32)')

    def test_ptrtoint_2(self):
        m = self.module()
        gv = ir.GlobalVariable(m, int32, "myconstant")
        c = gv.ptrtoint(int64)
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c), 'ptrtoint (ptr @"myconstant" to i64)')

            self.assertRaisesRegex(
                TypeError,
                r"can only ptrtoint\(\) to integer type, not 'ptr'",
                gv.ptrtoint,
                int64.as_pointer())
        else:
            self.assertEqual(str(c), 'ptrtoint (i32* @"myconstant" to i64)')

            self.assertRaisesRegex(
                TypeError,
                r"can only ptrtoint\(\) to integer type, not 'i64\*'",
                gv.ptrtoint,
                int64.as_pointer())

        c2 = ir.Constant(int32, 0)
        self.assertRaisesRegex(
            TypeError,
            r"can only call ptrtoint\(\) on pointer type, not 'i32'",
            c2.ptrtoint,
            int64)

    def test_inttoptr(self):
        one = ir.Constant(int32, 1)
        pi = ir.Constant(flt, 3.14)
        c = one.inttoptr(int64.as_pointer())

        self.assertRaises(TypeError, one.inttoptr, int64)
        self.assertRaises(TypeError, pi.inttoptr, int64.as_pointer())
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertEqual(str(c), 'inttoptr (i32 1 to ptr)')
        else:
            self.assertEqual(str(c), 'inttoptr (i32 1 to i64*)')

    def test_neg(self):
        one = ir.Constant(int32, 1)
        self.assertEqual(str(one.neg()), 'sub (i32 0, i32 1)')

    def test_not(self):
        one = ir.Constant(int32, 1)
        self.assertEqual(str(one.not_()), 'xor (i32 1, i32 -1)')

    def test_fneg(self):
        one = ir.Constant(flt, 1)
        self.assertEqual(str(one.fneg()), 'fneg (float 0x3ff0000000000000)')

    def test_int_binops(self):
        one = ir.Constant(int32, 1)
        two = ir.Constant(int32, 2)

        oracle = {one.shl:  'shl',  one.lshr: 'lshr', one.ashr: 'ashr',
                  one.add:  'add',  one.sub:  'sub',  one.mul:  'mul',
                  one.udiv: 'udiv', one.sdiv: 'sdiv', one.urem: 'urem',
                  one.srem: 'srem', one.or_:  'or',   one.and_: 'and',
                  one.xor:  'xor'}
        for fn, irop in oracle.items():
            self.assertEqual(str(fn(two)), irop + ' (i32 1, i32 2)')

        # unsigned integer compare
        oracle = {'==': 'eq', '!=': 'ne', '>':
                  'ugt', '>=': 'uge', '<': 'ult', '<=': 'ule'}
        for cop, cond in oracle.items():
            actual = str(one.icmp_unsigned(cop, two))
            expected = 'icmp ' + cond + ' (i32 1, i32 2)'
            self.assertEqual(actual, expected)

        # signed integer compare
        oracle = {'==': 'eq', '!=': 'ne',
                  '>': 'sgt', '>=': 'sge', '<': 'slt', '<=': 'sle'}
        for cop, cond in oracle.items():
            actual = str(one.icmp_signed(cop, two))
            expected = 'icmp ' + cond + ' (i32 1, i32 2)'
            self.assertEqual(actual, expected)

    def test_flt_binops(self):
        one = ir.Constant(flt, 1)
        two = ir.Constant(flt, 2)

        oracle = {one.fadd: 'fadd', one.fsub: 'fsub', one.fmul:  'fmul',
                  one.fdiv: 'fdiv', one.frem: 'frem'}
        for fn, irop in oracle.items():
            actual = str(fn(two))
            expected = irop + (' (float 0x3ff0000000000000,'
                               ' float 0x4000000000000000)')
            self.assertEqual(actual, expected)

        # ordered float compare
        oracle = {'==': 'oeq', '!=': 'one', '>': 'ogt', '>=': 'oge',
                  '<': 'olt', '<=': 'ole'}
        for cop, cond in oracle.items():
            actual = str(one.fcmp_ordered(cop, two))
            expected = 'fcmp ' + cond + (' (float 0x3ff0000000000000,'
                                         ' float 0x4000000000000000)')
            self.assertEqual(actual, expected)

        # unordered float compare
        oracle = {'==': 'ueq', '!=': 'une', '>': 'ugt', '>=': 'uge',
                  '<': 'ult', '<=': 'ule'}
        for cop, cond in oracle.items():
            actual = str(one.fcmp_unordered(cop, two))
            expected = 'fcmp ' + cond + (' (float 0x3ff0000000000000,'
                                         ' float 0x4000000000000000)')
            self.assertEqual(actual, expected)


class TestTransforms(TestBase):
    def test_call_transform(self):
        mod = ir.Module()
        foo = ir.Function(mod, ir.FunctionType(ir.VoidType(), ()), "foo")
        bar = ir.Function(mod, ir.FunctionType(ir.VoidType(), ()), "bar")
        builder = ir.IRBuilder()
        builder.position_at_end(foo.append_basic_block())
        call = builder.call(foo, ())
        self.assertEqual(call.callee, foo)
        modified = ir.replace_all_calls(mod, foo, bar)
        self.assertIn(call, modified)
        self.assertNotEqual(call.callee, foo)
        self.assertEqual(call.callee, bar)


class TestSingleton(TestBase):
    def test_undefined(self):
        self.assertIs(ir.Undefined, ir.values._Undefined())
        self.assertIs(ir.Undefined, copy.copy(ir.Undefined))
        self.assertIs(ir.Undefined, copy.deepcopy(ir.Undefined))
        self.assert_pickle_correctly(ir.Undefined)


if __name__ == '__main__':
    unittest.main()
