from collections import namedtuple
import inspect
import re
import numpy as np
import math
from textwrap import dedent
import unittest
import warnings

from numba.tests.support import (TestCase, override_config,
                                 ignore_internal_warnings)
from numba import jit, njit
from numba.core import types
from numba.core.datamodel import default_manager
from numba.core.errors import NumbaDebugInfoWarning
import llvmlite.binding as llvm

#NOTE: These tests are potentially sensitive to changes in SSA or lowering
# behaviour and may need updating should changes be made to the corresponding
# algorithms.


class TestDebugInfo(TestCase):
    """
    These tests only checks the compiled assembly for debuginfo.
    """

    def _getasm(self, fn, sig):
        fn.compile(sig)
        return fn.inspect_asm(sig)

    def _check(self, fn, sig, expect):
        asm = self._getasm(fn, sig=sig)
        m = re.search(r"\.section.+debug", asm, re.I)
        got = m is not None
        self.assertEqual(expect, got, msg='debug info not found in:\n%s' % asm)

    def test_no_debuginfo_in_asm(self):
        @jit(nopython=True, debug=False)
        def foo(x):
            return x

        self._check(foo, sig=(types.int32,), expect=False)

    def test_debuginfo_in_asm(self):
        @jit(nopython=True, debug=True)
        def foo(x):
            return x

        self._check(foo, sig=(types.int32,), expect=True)

    def test_environment_override(self):
        with override_config('DEBUGINFO_DEFAULT', 1):
            # Using default value
            @jit(nopython=True)
            def foo(x):
                return x
            self._check(foo, sig=(types.int32,), expect=True)

            # User override default
            @jit(nopython=True, debug=False)
            def bar(x):
                return x
            self._check(bar, sig=(types.int32,), expect=False)

    def test_llvm_inliner_flag_conflict(self):
        # bar will be marked as 'alwaysinline', but when DEBUGINFO_DEFAULT is
        # set functions are marked as 'noinline' this results in a conflict.
        # baz will be marked as 'noinline' as a result of DEBUGINFO_DEFAULT

        @njit(forceinline=True)
        def bar(x):
            return math.sin(x)

        @njit(forceinline=False)
        def baz(x):
            return math.cos(x)

        @njit
        def foo(x):
            a = bar(x)
            b = baz(x)
            return a, b

        # check it compiles
        with override_config('DEBUGINFO_DEFAULT', 1):
            result = foo(np.pi)

        self.assertPreciseEqual(result, foo.py_func(np.pi))

        # check the LLVM IR has bar marked as 'alwaysinline' and baz as noinline
        full_ir = foo.inspect_llvm(foo.signatures[0])
        module = llvm.parse_assembly(full_ir)
        name = foo.overloads[foo.signatures[0]].fndesc.mangled_name
        funcs = [x for x in module.functions if x.name == name]
        self.assertEqual(len(funcs), 1)
        func = funcs[0]

        # find the function calls and save the associated statements
        f_names = []
        for blk in func.blocks:
            for stmt in blk.instructions:
                if stmt.opcode == 'call':
                    # stmt.function.name  This is the function being called
                    f_names.append(str(stmt).strip())

        # Need to check there's two specific things in the calls in the IR
        # 1. a call to the llvm.sin.f64 intrinsic, this is from the inlined bar
        # 2. a call to the baz function, this is from the noinline baz
        found_sin = False
        found_baz = False
        baz_name = baz.overloads[baz.signatures[0]].fndesc.mangled_name
        for x in f_names:
            if not found_sin and re.match('.*llvm.sin.f64.*', x):
                found_sin = True
            if not found_baz and re.match(f'.*{baz_name}.*', x):
                found_baz = True

        self.assertTrue(found_sin)
        self.assertTrue(found_baz)


class TestDebugInfoEmission(TestCase):
    """ Tests that debug info is emitted correctly.
    """

    _NUMBA_OPT_0_ENV = {'NUMBA_OPT': '0'}

    def _get_llvmir(self, fn, sig):
        with override_config('OPT', 0):
            fn.compile(sig)
            return fn.inspect_llvm(sig)

    def _get_metadata(self, fn, sig):
        ll = self._get_llvmir(fn, sig).splitlines()
        meta_re = re.compile(r'![0-9]+ =.*')
        metadata = []
        for line in ll:
            if meta_re.match(line):
                metadata.append(line)
        return metadata

    def _get_metadata_map(self, metadata):
        """Gets the map of DI label to md, e.g.
        '!33' -> '!{!"branch_weights", i32 1, i32 99}'
        """
        metadata_definition_map = dict()
        meta_definition_split = re.compile(r'(![0-9]+) = (.*)')
        for line in metadata:
            matched = meta_definition_split.match(line)
            if matched:
                dbg_val, info = matched.groups()
                metadata_definition_map[dbg_val] = info
        return metadata_definition_map

    def _get_lines_from_debuginfo(self, metadata):
        # Get the lines contained in the debug info
        md_def_map = self._get_metadata_map(metadata)

        lines = set()
        for md in md_def_map.values():
            m = re.match(r"!DILocation\(line: (\d+),", md)
            if m:
                ln = int(m.group(1))
                lines.add(ln)
        return lines

    def test_DW_LANG(self):

        @njit(debug=True)
        def foo():
            pass

        metadata = self._get_metadata(foo, sig=())
        DICompileUnit = metadata[0]
        self.assertEqual('!0', DICompileUnit[:2])
        self.assertIn('!DICompileUnit(language: DW_LANG_C_plus_plus',
                      DICompileUnit)
        self.assertIn('producer: "clang (Numba)"', DICompileUnit)

    def test_DILocation(self):
        """ Tests that DILocation information is reasonable.
        """
        @njit(debug=True, error_model='numpy')
        def foo(a):
            b = a + 1.23
            c = b * 2.34
            d = b / c
            print(d)
            return d

        # the above produces LLVM like:
        # define function() {
        # entry:
        #   alloca
        #   store 0 to alloca
        #   <arithmetic for doing the operations on b, c, d>
        #   setup for print
        #   branch
        # other_labels:
        # ... <elided>
        # }
        #
        # The following checks that:
        # * the alloca and store have no !dbg
        # * the arithmetic occurs in the order defined and with !dbg
        # * that the !dbg entries are monotonically increasing in value with
        #   source line number

        sig = (types.float64,)
        metadata = self._get_metadata(foo, sig=sig)
        full_ir = self._get_llvmir(foo, sig=sig)

        module = llvm.parse_assembly(full_ir)

        name = foo.overloads[foo.signatures[0]].fndesc.mangled_name
        funcs = [x for x in module.functions if x.name == name]
        self.assertEqual(len(funcs), 1)
        func = funcs[0]
        blocks = [x for x in func.blocks]
        self.assertGreater(len(blocks), 1)
        block = blocks[0]

        # Find non-call/non-memory instr and check the sequence is as expected
        instrs = [x for x in block.instructions if x.opcode not in
                  ['call', 'load', 'store']]
        op_expect = {'fadd', 'fmul', 'fdiv'}
        started = False
        for x in instrs:
            if x.opcode in op_expect:
                op_expect.remove(x.opcode)
                if not started:
                    started = True
            elif op_expect and started:
                self.fail("Math opcodes are not contiguous")
        self.assertFalse(op_expect, "Math opcodes were not found")

        # Parse out metadata from end of each line, check it monotonically
        # ascends with LLVM source line. Also store all the dbg references,
        # these will be checked later.
        line2dbg = set()
        re_dbg_ref = re.compile(r'.*!dbg (![0-9]+).*$')
        found = -1
        for instr in instrs:
            inst_as_str = str(instr)
            matched = re_dbg_ref.match(inst_as_str)
            if not matched:
                # if there's no match, ensure it is one of alloca or store,
                # it's important that the zero init/alloca instructions have
                # no dbg data
                accepted = ('alloca ', 'store ')
                self.assertTrue(any([x in inst_as_str for x in accepted]))
                continue
            groups = matched.groups()
            self.assertEqual(len(groups), 1)
            dbg_val = groups[0]
            int_dbg_val = int(dbg_val[1:])
            if found >= 0:
                self.assertTrue(int_dbg_val >= found)
            found = int_dbg_val
            # some lines will alias dbg info, this is fine, it's only used to
            # make sure that the line numbers are correct WRT python
            line2dbg.add(dbg_val)

        pysrc, pysrc_line_start = inspect.getsourcelines(foo)

        # build a map of dbg reference to DI* information
        metadata_definition_map = self._get_metadata_map(metadata)

        # Pull out metadata entries referred to by the llvm line end !dbg
        # check they match the python source, the +2 is for the @njit decorator
        # and the function definition line.
        offsets = [0,  # b = a + 1
                   1,  # a * 2.34
                   2,  # d = b / c
                   3,  # print(d)
                   ]
        pyln_range = [pysrc_line_start + 2 + x for x in offsets]

        # do the check
        for (k, line_no) in zip(sorted(line2dbg, key=lambda x: int(x[1:])),
                                pyln_range):
            dilocation_info = metadata_definition_map[k]
            self.assertIn(f'line: {line_no}', dilocation_info)

        # Check that variable "a" is declared as on the same line as function
        # definition.
        expr = r'.*!DILocalVariable\(name: "a",.*line: ([0-9]+),.*'
        match_local_var_a = re.compile(expr)
        for entry in metadata_definition_map.values():
            matched = match_local_var_a.match(entry)
            if matched:
                groups = matched.groups()
                self.assertEqual(len(groups), 1)
                dbg_line = int(groups[0])
                # +1 for the decorator.
                # Recall that Numba's DWARF refers to the "def" line, but
                # `inspect` uses the decorator as the first line.
                defline = pysrc_line_start + 1
                self.assertEqual(dbg_line, defline)
                break
        else:
            self.fail('Assertion on DILocalVariable not made')

    @TestCase.run_test_in_subprocess(envvars=_NUMBA_OPT_0_ENV)
    def test_DILocation_entry_blk(self):
        # Needs a subprocess as jitting literally anything at any point in the
        # lifetime of the process ends up with a codegen at opt 3. This is not
        # amenable to this test!
        # This test relies on the CFG not being simplified as it checks the jump
        # from the entry block to the first basic block. Force OPT as 0, if set
        # via the env var the targetmachine and various pass managers all end up
        # at OPT 0 and the IR is minimally transformed prior to lowering to ELF.
        #
        # This tests that the unconditional jump emitted at the tail of
        # the entry block has no debug metadata associated with it. In practice,
        # if debug metadata is associated with it, it manifests as the
        # prologue_end being associated with the end_sequence or similar (due to
        # the way code gen works for the entry block).

        @njit(debug=True)
        def foo(a):
            return a + 1
        foo(123)

        full_ir = foo.inspect_llvm(foo.signatures[0])
        # The above produces LLVM like:
        #
        # define function() {
        # entry:
        #   alloca
        #   store 0 to alloca
        #   unconditional jump to body:
        #
        # body:
        # ... <elided>
        # }

        module = llvm.parse_assembly(full_ir)
        name = foo.overloads[foo.signatures[0]].fndesc.mangled_name
        funcs = [x for x in module.functions if x.name == name]
        self.assertEqual(len(funcs), 1)
        func = funcs[0]
        blocks = [x for x in func.blocks]
        self.assertEqual(len(blocks), 2)
        entry_block, body_block = blocks

        # Assert that the tail of the entry block is an unconditional jump to
        # the body block and that the jump has no associated debug info.
        entry_instr = [x for x in entry_block.instructions]
        ujmp = entry_instr[-1]
        self.assertEqual(ujmp.opcode, 'br')
        ujmp_operands = [x for x in ujmp.operands]
        self.assertEqual(len(ujmp_operands), 1)
        target_data = ujmp_operands[0]
        target = str(target_data).split(':')[0].strip()
        # check the unconditional jump target is to the body block
        self.assertEqual(target, body_block.name)
        # check the uncondition jump instr itself has no metadata
        self.assertTrue(str(ujmp).endswith(target))

    @TestCase.run_test_in_subprocess(envvars=_NUMBA_OPT_0_ENV)
    def test_DILocation_decref(self):
        """ This tests that decref's generated from `ir.Del`s as variables go
        out of scope do not have debuginfo associated with them (the location of
        `ir.Del` is an implementation detail).
        """

        @njit(debug=True)
        def sink(*x):
            pass

        # This function has many decrefs!
        @njit(debug=True)
        def foo(a):
            x = (a, a)
            if a[0] == 0:
                sink(x)
                return 12
            z = x[0][0]
            return z

        sig = (types.float64[::1],)
        full_ir = self._get_llvmir(foo, sig=sig)

        # make sure decref lines end with `meminfo.<number>)` without !dbg info.
        count = 0
        for line in full_ir.splitlines():
            line_stripped = line.strip()
            if line_stripped.startswith('call void @NRT_decref'):
                self.assertRegex(line, r'.*meminfo\.[0-9]+\)$')
                count += 1
        self.assertGreater(count, 0) # make sure there were some decrefs!

    def test_DILocation_undefined(self):
        """ Tests that DILocation information for undefined vars is associated
        with the line of the function definition (so it ends up in the prologue)
        """
        @njit(debug=True)
        def foo(n):
            if n:
                if n > 0:
                    c = 0
                return c
            else:
                # variable c is not defined in this branch
                c += 1
                return c

        sig = (types.intp,)
        metadata = self._get_metadata(foo, sig=sig)
        pysrc, pysrc_line_start = inspect.getsourcelines(foo)
        # Looks for versions of variable "c" and captures the line number
        expr = r'.*!DILocalVariable\(name: "c\$?[0-9]?",.*line: ([0-9]+),.*'
        matcher = re.compile(expr)
        associated_lines = set()
        for md in metadata:
            match = matcher.match(md)
            if match:
                groups = match.groups()
                self.assertEqual(len(groups), 1)
                associated_lines.add(int(groups[0]))
        # 3 versions of 'c': `c = 0`, `return c`, `c+=1`
        self.assertEqual(len(associated_lines), 3)
        self.assertIn(pysrc_line_start, associated_lines)

    def test_DILocation_versioned_variables(self):
        """ Tests that DILocation information for versions of variables matches
        up to their definition site."""
        # Note: there's still something wrong in the DI/SSA naming, the ret c is
        # associated with the logically first definition.

        @njit(debug=True)
        def foo(n):
            if n:
                c = 5
            else:
                c = 1
            # prevents inline of return on py310
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return c

        sig = (types.intp,)
        metadata = self._get_metadata(foo, sig=sig)
        pysrc, pysrc_line_start = inspect.getsourcelines(foo)

        # Looks for SSA versioned names i.e. <basename>$<version id> of the
        # variable 'c' and captures the line
        expr = r'.*!DILocalVariable\(name: "c\$[0-9]?",.*line: ([0-9]+),.*'
        matcher = re.compile(expr)
        associated_lines = set()
        for md in metadata:
            match = matcher.match(md)
            if match:
                groups = match.groups()
                self.assertEqual(len(groups), 1)
                associated_lines.add(int(groups[0]))
        self.assertEqual(len(associated_lines), 2) # 2 SSA versioned names 'c'

        # Now find the `c = ` lines in the python source
        py_lines = set()
        for ix, pyln in enumerate(pysrc):
            if 'c = ' in pyln:
                py_lines.add(ix + pysrc_line_start)
        self.assertEqual(len(py_lines), 2) # 2 assignments to c

        # check that the DILocation from the DI for `c` matches the python src
        self.assertEqual(associated_lines, py_lines)

    def test_numeric_scalars(self):
        """ Tests that dwarf info is correctly emitted for numeric scalars."""

        DI = namedtuple('DI', 'name bits encoding')

        type_infos = {np.float32: DI("float32", 32, "DW_ATE_float"),
                      np.float64: DI("float64", 64, "DW_ATE_float"),
                      np.int8: DI("int8", 8, "DW_ATE_signed"),
                      np.int16: DI("int16", 16, "DW_ATE_signed"),
                      np.int32: DI("int32", 32, "DW_ATE_signed"),
                      np.int64: DI("int64", 64, "DW_ATE_signed"),
                      np.uint8: DI("uint8", 8, "DW_ATE_unsigned"),
                      np.uint16: DI("uint16", 16, "DW_ATE_unsigned"),
                      np.uint32: DI("uint32", 32, "DW_ATE_unsigned"),
                      np.uint64: DI("uint64", 64, "DW_ATE_unsigned"),
                      np.complex64: DI("complex64", 64,
                                       "DW_TAG_structure_type"),
                      np.complex128: DI("complex128", 128,
                                        "DW_TAG_structure_type"),}

        for ty, dwarf_info in type_infos.items():

            @njit(debug=True)
            def foo():
                a = ty(10)
                return a

            metadata = self._get_metadata(foo, sig=())
            metadata_definition_map = self._get_metadata_map(metadata)

            for k, v in metadata_definition_map.items():
                if 'DILocalVariable(name: "a"' in v:
                    lvar = metadata_definition_map[k]
                    break
            else:
                assert 0, "missing DILocalVariable 'a'"

            type_marker = re.match('.*type: (![0-9]+).*', lvar).groups()[0]
            type_decl = metadata_definition_map[type_marker]

            if 'DW_ATE' in dwarf_info.encoding:
                expected = (f'!DIBasicType(name: "{dwarf_info.name}", '
                            f'size: {dwarf_info.bits}, '
                            f'encoding: {dwarf_info.encoding})')
                self.assertEqual(type_decl, expected)
            else: # numerical complex type
                # Don't match the whole string, just the known parts
                raw_flt = 'float' if dwarf_info.bits == 64 else 'double'
                expected = (f'distinct !DICompositeType('
                            f'tag: {dwarf_info.encoding}, '
                            f'name: "{dwarf_info.name} '
                            f'({{{raw_flt}, {raw_flt}}})", '
                            f'size: {dwarf_info.bits}')
                self.assertIn(expected, type_decl)

    def test_arrays(self):

        @njit(debug=True)
        def foo():
            a = np.ones((2, 3), dtype=np.float64)
            return a

        metadata = self._get_metadata(foo, sig=())
        metadata_definition_map = self._get_metadata_map(metadata)

        for k, v in metadata_definition_map.items():
            if 'DILocalVariable(name: "a"' in v:
                lvar = metadata_definition_map[k]
                break
        else:
            assert 0, "missing DILocalVariable 'a'"

        type_marker = re.match('.*type: (![0-9]+).*', lvar).groups()[0]
        type_decl = metadata_definition_map[type_marker]

        # check type
        self.assertIn("!DICompositeType(tag: DW_TAG_structure_type", type_decl)
        # check name encoding
        self.assertIn(f'name: "{str(types.float64[:, ::1])}', type_decl)

        # pop out the "elements" of the composite type
        match_elements = re.compile(r'.*elements: (![0-9]+),.*')
        elem_matches = match_elements.match(type_decl).groups()
        self.assertEqual(len(elem_matches), 1)
        elem_match = elem_matches[0]
        # The match should be something like, it's the elements from an array
        # data model.
        # !{!35, !36, !37, !39, !40, !43, !45}'
        struct_markers = metadata_definition_map[elem_match]
        struct_pattern = '!{' + '(![0-9]+), ' * 6 + '(![0-9]+)}'
        match_struct = re.compile(struct_pattern)
        struct_member_matches = match_struct.match(struct_markers).groups()
        self.assertIsNotNone(struct_member_matches is not None)
        data_model = default_manager.lookup(types.float64[:, ::1])
        self.assertEqual(len(struct_member_matches), len(data_model._fields))

        ptr_size = types.intp.bitwidth
        ptr_re = (r'!DIDerivedType\(tag: DW_TAG_pointer_type, '
                  rf'baseType: ![0-9]+, size: {ptr_size}\)')
        int_re = (rf'!DIBasicType\(name: "int{ptr_size}", size: {ptr_size}, '
                  r'encoding: DW_ATE_signed\)')
        utuple_re = (r'!DICompositeType\(tag: DW_TAG_array_type, '
                     rf'name: "UniTuple\(int{ptr_size} x 2\) '
                     rf'\(\[2 x i{ptr_size}\]\)", baseType: ![0-9]+, '
                     rf'size: {2 * ptr_size}, elements: ![0-9]+, '
                     rf'identifier: "\[2 x i{ptr_size}\]"\)')
        expected = {'meminfo': ptr_re,
                    'parent': ptr_re,
                    'nitems': int_re,
                    'itemsize': int_re,
                    'data': ptr_re,
                    'shape': utuple_re,
                    'strides': utuple_re}

        # look for `baseType: <>` for the type
        base_type_pattern = r'!DIDerivedType\(.*, baseType: (![0-9]+),.*'
        base_type_matcher = re.compile(base_type_pattern)

        for ix, field in enumerate(data_model._fields):
            derived_type = metadata_definition_map[struct_member_matches[ix]]
            self.assertIn("DIDerivedType", derived_type)
            self.assertIn(f'name: "{field}"', derived_type)
            base_type_match = base_type_matcher.match(derived_type)
            base_type_matches = base_type_match.groups()
            self.assertEqual(len(base_type_matches), 1)
            base_type_marker = base_type_matches[0]
            data_type = metadata_definition_map[base_type_marker]
            self.assertRegex(data_type, expected[field])

    def test_debug_optnone(self):
        def get_debug_lines(fn):
            metadata = self._get_metadata(fn, fn.signatures[0])
            lines = self._get_lines_from_debuginfo(metadata)
            return lines

        def get_func_attrs(fn):
            cres = fn.overloads[fn.signatures[0]]
            lib = cres.library
            fn = lib._final_module.get_function(cres.fndesc.mangled_name)
            attrs = set(b' '.join(fn.attributes).split())
            return attrs

        def foo():
            n = 10
            c = 0
            for i in range(n):
                c += i
            return c

        foo_debug = njit(debug=True)(foo)
        foo_debug_optnone = njit(debug=True, _dbg_optnone=True)(foo)
        foo_debug_optnone_inline = njit(debug=True, _dbg_optnone=True,
                                        forceinline=True)(foo)

        firstline = foo.__code__.co_firstlineno

        expected_info = {}
        expected_info[foo_debug] = dict(
            # just the dummy line-0 and the line of the return statement
            lines={0, firstline + 5},
            must_have_attrs=set(),
            must_not_have_attrs=set([b"optnone"]),
        )
        expected_info[foo_debug_optnone] = dict(
            # all the lines should be included
            lines=set(range(firstline + 1, firstline + 6)),
            must_have_attrs=set([b"optnone"]),
            must_not_have_attrs=set(),
        )
        expected_info[foo_debug_optnone_inline] = dict(
            # optnone=True is overridden by forceinline, so this looks like the
            # foo_debug version
            lines={0, firstline + 5},
            must_have_attrs=set([b"alwaysinline"]),
            must_not_have_attrs=set([b"optnone"]),
        )

        expected_ret = foo()

        for udt, expected in expected_info.items():
            with self.subTest(udt.targetoptions):
                got = udt()
                self.assertEqual(got, expected_ret)

                # Compare the line locations in the debug info.
                self.assertEqual(get_debug_lines(udt), expected["lines"])

                # Check for attributes on the LLVM function
                attrs = get_func_attrs(udt)
                must_have = expected["must_have_attrs"]
                self.assertEqual(attrs & must_have, must_have)
                must_not_have = expected["must_not_have_attrs"]
                self.assertFalse(attrs & must_not_have)

    def test_omitted_arg(self):
        # See issue 7726
        @njit(debug=True)
        def foo(missing=None):
            pass

        # check that it will actually compile (verifies DI emission is ok)
        with override_config('DEBUGINFO_DEFAULT', 1):
            foo()

        metadata = self._get_metadata(foo, sig=(types.Omitted(None),))
        metadata_definition_map = self._get_metadata_map(metadata)

        # Find DISubroutineType
        tmp_disubr = []
        for md in metadata:
            if "DISubroutineType" in md:
                tmp_disubr.append(md)
        self.assertEqual(len(tmp_disubr), 1)
        disubr = tmp_disubr.pop()

        disubr_matched = re.match(r'.*!DISubroutineType\(types: ([!0-9]+)\)$',
                                  disubr)
        self.assertIsNotNone(disubr_matched)
        disubr_groups = disubr_matched.groups()
        self.assertEqual(len(disubr_groups), 1)
        disubr_meta = disubr_groups[0]

        # Find the types in the DISubroutineType arg list
        disubr_types = metadata_definition_map[disubr_meta]
        disubr_types_matched = re.match(r'!{(.*)}', disubr_types)
        self.assertIsNotNone(disubr_matched)
        disubr_types_groups = disubr_types_matched.groups()
        self.assertEqual(len(disubr_types_groups), 1)

        # fetch out and assert the last argument type, should be void *
        md_fn_arg = [x.strip() for x in disubr_types_groups[0].split(',')][-1]
        arg_ty = metadata_definition_map[md_fn_arg]
        expected_arg_ty = (r'^.*!DICompositeType\(tag: DW_TAG_structure_type, '
                           r'name: "Anonymous struct \({}\)", elements: '
                           r'(![0-9]+), identifier: "{}"\)')
        self.assertRegex(arg_ty, expected_arg_ty)
        md_base_ty = re.match(expected_arg_ty, arg_ty).groups()[0]
        base_ty = metadata_definition_map[md_base_ty]
        # expect ir.LiteralStructType([])
        self.assertEqual(base_ty, ('!{}'))

    def test_missing_source(self):
        strsrc = """
        def foo():
            return 1
        """
        l = dict()
        exec(dedent(strsrc), {}, l)
        foo = njit(debug=True)(l['foo'])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaDebugInfoWarning)
            ignore_internal_warnings()
            foo()

        self.assertEqual(len(w), 1)
        found = w[0]
        self.assertEqual(found.category, NumbaDebugInfoWarning)
        msg = str(found.message)
        # make sure the warning contains the right message
        self.assertIn('Could not find source for function', msg)
        # and refers to the offending function
        self.assertIn(str(foo.py_func), msg)

    def test_irregularly_indented_source(self):

        @njit(debug=True)
        def foo():
# NOTE: THIS COMMENT MUST START AT COLUMN 0 FOR THIS SAMPLE CODE TO BE VALID # noqa: E115, E501
            return 1

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaDebugInfoWarning)
            ignore_internal_warnings()
            foo()

        # No warnings
        self.assertEqual(len(w), 0)

        metadata = self._get_metadata(foo, foo.signatures[0])
        lines = self._get_lines_from_debuginfo(metadata)
        # Only one line
        self.assertEqual(len(lines), 1)


if __name__ == '__main__':
    unittest.main()
