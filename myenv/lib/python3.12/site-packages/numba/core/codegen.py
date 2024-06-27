import warnings
import functools
import locale
import weakref
import ctypes
import html
import textwrap

import llvmlite.binding as ll
import llvmlite.ir as llvmir

from abc import abstractmethod, ABCMeta
from numba.core import utils, config, cgutils
from numba.core.llvm_bindings import create_pass_manager_builder
from numba.core.runtime.nrtopt import remove_redundant_nrt_refct
from numba.core.runtime import rtsys
from numba.core.compiler_lock import require_global_compiler_lock
from numba.core.errors import NumbaInvalidConfigWarning
from numba.misc.inspection import disassemble_elf_to_cfg
from numba.misc.llvm_pass_timings import PassTimingsCollection


_x86arch = frozenset(['x86', 'i386', 'i486', 'i586', 'i686', 'i786',
                      'i886', 'i986'])


def _is_x86(triple):
    arch = triple.split('-')[0]
    return arch in _x86arch


def _parse_refprune_flags():
    """Parse refprune flags from the `config`.

    Invalid values are ignored an warn via a `NumbaInvalidConfigWarning`
    category.

    Returns
    -------
    flags : llvmlite.binding.RefPruneSubpasses
    """
    flags = config.LLVM_REFPRUNE_FLAGS.split(',')
    if not flags:
        return 0
    val = 0
    for item in flags:
        item = item.strip()
        try:
            val |= getattr(ll.RefPruneSubpasses, item.upper())
        except AttributeError:
            warnings.warn(f"invalid refprune flags {item!r}",
                          NumbaInvalidConfigWarning)
    return val


def dump(header, body, lang):
    if config.HIGHLIGHT_DUMPS:
        try:
            import pygments
        except ImportError:
            msg = "Please install pygments to see highlighted dumps"
            raise ValueError(msg)
        else:
            from pygments import highlight
            from pygments.lexers import GasLexer as gas_lexer
            from pygments.lexers import LlvmLexer as llvm_lexer
            from pygments.formatters import Terminal256Formatter
            from numba.misc.dump_style import by_colorscheme

            lexer_map = {'llvm': llvm_lexer, 'asm': gas_lexer}
            lexer = lexer_map[lang]
            def printer(arg):
                print(highlight(arg, lexer(),
                      Terminal256Formatter(style=by_colorscheme())))
    else:
        printer = print
    print('=' * 80)
    print(header.center(80, '-'))
    printer(body)
    print('=' * 80)


class _CFG(object):
    """
    Wraps the CFG graph for different display method.

    Instance of the class can be stringified (``__repr__`` is defined) to get
    the graph in DOT format.  The ``.display()`` method plots the graph in
    PDF.  If in IPython notebook, the returned image can be inlined.
    """
    def __init__(self, cres, name, py_func, **kwargs):
        self.cres = cres
        self.name = name
        self.py_func = py_func
        fn = cres.get_function(name)
        self.dot = ll.get_function_cfg(fn)
        self.kwargs = kwargs

    def pretty_printer(self, filename=None, view=None, render_format=None,
                       highlight=True,
                       interleave=False, strip_ir=False, show_key=True,
                       fontsize=10):
        """
        "Pretty" prints the DOT graph of the CFG.
        For explanation of the parameters see the docstring for
        numba.core.dispatcher::inspect_cfg.
        """
        import graphviz as gv
        import re
        import json
        import inspect
        from llvmlite import binding as ll
        from numba.typed import List
        from types import SimpleNamespace
        from collections import defaultdict

        _default = False
        _highlight = SimpleNamespace(incref=_default,
                                    decref=_default,
                                    returns=_default,
                                    raises=_default,
                                    meminfo=_default,
                                    branches=_default,
                                    llvm_intrin_calls=_default,
                                    function_calls=_default,)
        _interleave = SimpleNamespace(python=_default, lineinfo=_default)

        def parse_config(_config, kwarg):
            """ Parses the kwarg into a consistent format for use in configuring
            the Digraph rendering. _config is the configuration instance to
            update, kwarg is the kwarg on which to base the updates.
            """
            if isinstance(kwarg, bool):
                for attr in _config.__dict__:
                    setattr(_config, attr, kwarg)
            elif isinstance(kwarg, dict):
                for k, v in kwarg.items():
                    if k not in _config.__dict__:
                        raise ValueError("Unexpected key in kwarg: %s" % k)
                    if isinstance(v, bool):
                        setattr(_config, k, v)
                    else:
                        msg = "Unexpected value for key: %s, got:%s"
                        raise ValueError(msg % (k, v))
            elif isinstance(kwarg, set):
                for item in kwarg:
                    if item not in _config.__dict__:
                        raise ValueError("Unexpected key in kwarg: %s" % item)
                    else:
                        setattr(_config, item, True)
            else:
                msg = "Unhandled configuration type for kwarg %s"
                raise ValueError(msg % type(kwarg))

        parse_config(_highlight, highlight)
        parse_config(_interleave, interleave)

        # This is the colour scheme. The graphviz HTML label renderer only takes
        # names for colours: https://www.graphviz.org/doc/info/shapes.html#html
        cs = defaultdict(lambda: 'white') # default bg colour is white
        cs['marker'] = 'orange'
        cs['python'] = 'yellow'
        cs['truebr'] = 'green'
        cs['falsebr'] = 'red'
        cs['incref'] = 'cyan'
        cs['decref'] = 'turquoise'
        cs['raise'] = 'lightpink'
        cs['meminfo'] = 'lightseagreen'
        cs['return'] = 'purple'
        cs['llvm_intrin_calls'] = 'rosybrown'
        cs['function_calls'] = 'tomato'

        # Get the raw dot format information from LLVM and the LLVM IR
        fn = self.cres.get_function(self.name)
        #raw_dot = ll.get_function_cfg(fn).replace('\\l...', '')
        llvm_str = self.cres.get_llvm_str()

        def get_metadata(llvm_str):
            """ Gets the metadata entries from the LLVM IR, these look something
            like '!123 = INFORMATION'. Returns a map of metadata key to metadata
            value, i.e. from the example {'!123': INFORMATION}"""
            md = {}
            metadata_entry = re.compile(r'(^[!][0-9]+)(\s+=\s+.*)')
            for x in llvm_str.splitlines():
                match = metadata_entry.match(x)
                if match is not None:
                    g = match.groups()
                    if g is not None:
                        assert len(g) == 2
                        md[g[0]] = g[1]
            return md

        md = get_metadata(llvm_str)

        # setup digraph with initial properties
        def init_digraph(name, fname, fontsize):
            # name and fname are arbitrary graph and file names, they appear in
            # some rendering formats, the fontsize determines the output
            # fontsize.

            # truncate massive mangled names as file names as it causes OSError
            # when trying to render to pdf
            cmax = 200
            if len(fname) > cmax:
                wstr = (f'CFG output filename "{fname}" exceeds maximum '
                        f'supported length, it will be truncated.')
                warnings.warn(wstr, NumbaInvalidConfigWarning)
                fname = fname[:cmax]
            f = gv.Digraph(name, filename=fname)
            f.attr(rankdir='TB')
            f.attr('node', shape='none', fontsize='%s' % str(fontsize))
            return f

        f = init_digraph(self.name, self.name, fontsize)

        # A lot of regex is needed to parse the raw dot output. This output
        # contains a mix of LLVM IR in the labels, and also DOT markup.

        # DOT syntax, matches a "port" (where the tail of an edge starts)
        port_match = re.compile('.*{(.*)}.*')
        # DOT syntax, matches the "port" value from a found "port_match"
        port_jmp_match = re.compile('.*<(.*)>(.*)')
        # LLVM syntax, matches a LLVM debug marker
        metadata_marker = re.compile(r'.*!dbg\s+(![0-9]+).*')
        # LLVM syntax, matches a location entry
        location_expr = (r'.*!DILocation\(line:\s+([0-9]+),'
                         r'\s+column:\s+([0-9]),.*')
        location_entry = re.compile(location_expr)
        # LLVM syntax, matches LLVMs internal debug value calls
        dbg_value = re.compile(r'.*call void @llvm.dbg.value.*')
        # LLVM syntax, matches tokens for highlighting
        nrt_incref = re.compile(r"@NRT_incref\b")
        nrt_decref = re.compile(r"@NRT_decref\b")
        nrt_meminfo = re.compile("@NRT_MemInfo")
        ll_intrin_calls = re.compile(r".*call.*@llvm\..*")
        ll_function_call = re.compile(r".*call.*@.*")
        ll_raise = re.compile(r"store .*\!numba_exception_output.*")
        ll_return = re.compile("ret i32 [^1],?.*")

        # wrapper function for line wrapping LLVM lines
        def wrap(s):
            return textwrap.wrap(s, width=120, subsequent_indent='... ')

        # function to fix (sometimes escaped for DOT!) LLVM IR etc that needs to
        # be HTML escaped
        def clean(s):
            # Grab first 300 chars only, 1. this should be enough to identify
            # the token and it keeps names short. 2. graphviz/dot has a maximum
            # buffer size near 585?!, with additional transforms it's hard to
            # know if this would be exceeded. 3. hash of the token string is
            # written into the rendering to permit exact identification against
            # e.g. LLVM IR dump if necessary.
            n = 300
            if len(s) > n:
                hs = str(hash(s))
                s = '{}...<hash={}>'.format(s[:n], hs)
            s = html.escape(s) # deals with  &, < and >
            s = s.replace('\\{', "&#123;")
            s = s.replace('\\}', "&#125;")
            s = s.replace('\\', "&#92;")
            s = s.replace('%', "&#37;")
            s = s.replace('!', "&#33;")
            return s

        # These hold the node and edge ids from the raw dot information. They
        # are used later to wire up a new DiGraph that has the same structure
        # as the raw dot but with new nodes.
        node_ids = {}
        edge_ids = {}

        # Python source lines, used if python source interleave is requested
        if _interleave.python:
            src_code, firstlineno = inspect.getsourcelines(self.py_func)

        # This is the dot info from LLVM, it's in DOT form and has continuation
        # lines, strip them and then re-parse into `dot_json` form for use in
        # producing a formatted output.
        raw_dot = ll.get_function_cfg(fn).replace('\\l...', '')
        json_bytes = gv.Source(raw_dot).pipe(format='dot_json')
        jzon = json.loads(json_bytes.decode('utf-8'))

        idc = 0
        # Walk the "objects" (nodes) in the DOT output
        for obj in jzon['objects']:
            # These are used to keep tabs on the current line and column numbers
            # as per the markers. They are tracked so as to make sure a marker
            # is only emitted if there's a change in the marker.
            cur_line, cur_col = -1, -1
            label = obj['label']
            name = obj['name']
            gvid = obj['_gvid']
            node_ids[gvid] = name
            # Label is DOT format, it needs the head and tail removing and then
            # splitting for walking.
            label = label[1:-1]
            lines = label.split('\\l')

            # Holds the new lines
            new_lines = []

            # Aim is to produce an HTML table a bit like this:
            #
            # |------------|
            # | HEADER     | <-- this is the block header
            # |------------|
            # | LLVM SRC   | <--
            # | Marker?    | < this is the label/block body
            # | Python src?| <--
            # |------------|
            # | T   |  F   |  <-- this is the "ports", also determines col_span
            # --------------
            #

            # This is HTML syntax, its the column span. If there's a switch or a
            # branch at the bottom of the node this is rendered as multiple
            # columns in a table. First job is to go and render that and work
            # out how many columns are needed as that dictates how many columns
            # the rest of the source lines must span. In DOT syntax the places
            # that edges join nodes are referred to as "ports". Syntax in DOT
            # is like `node:port`.
            col_span = 1

            # First see if there is a port entry for this node
            port_line = ''
            matched = port_match.match(lines[-1])
            sliced_lines = lines
            if matched is not None:
                # There is a port
                ports = matched.groups()[0]
                ports_tokens = ports.split('|')
                col_span = len(ports_tokens)
                # Generate HTML table data cells, one for each port. If the
                # ports correspond to a branch then they can optionally
                # highlighted based on T/F.
                tdfmt = ('<td BGCOLOR="{}" BORDER="1" ALIGN="center" '
                         'PORT="{}">{}</td>')
                tbl_data = []
                if _highlight.branches:
                    colors = {'T': cs['truebr'], 'F': cs['falsebr']}
                else:
                    colors = {}
                for tok in ports_tokens:
                    target, value = port_jmp_match.match(tok).groups()
                    color = colors.get(value, 'white')
                    tbl_data.append(tdfmt.format(color, target, value))
                port_line = ''.join(tbl_data)
                # Drop the last line from the rest of the parse as it's the port
                # and just been dealt with.
                sliced_lines = lines[:-1]

            # loop peel the block header, it needs a HTML border
            fmtheader = ('<tr><td BGCOLOR="{}" BORDER="1" ALIGN="left" '
                         'COLSPAN="{}">{}</td></tr>')
            new_lines.append(fmtheader.format(cs['default'], col_span,
                                              clean(sliced_lines[0].strip())))

            # process rest of block creating the table row at a time.
            fmt = ('<tr><td BGCOLOR="{}" BORDER="0" ALIGN="left" '
                   'COLSPAN="{}">{}</td></tr>')

            def metadata_interleave(l, new_lines):
                """
                Search line `l` for metadata associated with python or line info
                and inject it into `new_lines` if requested.
                """
                matched = metadata_marker.match(l)
                if matched is not None:
                    # there's a metadata marker
                    g = matched.groups()
                    if g is not None:
                        assert len(g) == 1, g
                        marker = g[0]
                        debug_data = md.get(marker, None)
                        if debug_data is not None:
                            # and the metadata marker has a corresponding piece
                            # of metadata
                            ld = location_entry.match(debug_data)
                            if ld is not None:
                                # and the metadata is line info... proceed
                                assert len(ld.groups()) == 2, ld
                                line, col = ld.groups()
                                # only emit a new marker if the line number in
                                # the metadata is "new".
                                if line != cur_line or col != cur_col:
                                    if _interleave.lineinfo:
                                        mfmt = 'Marker %s, Line %s, column %s'
                                        mark_line = mfmt % (marker, line, col)
                                        ln = fmt.format(cs['marker'], col_span,
                                                        clean(mark_line))
                                        new_lines.append(ln)
                                    if _interleave.python:
                                        # TODO:
                                        # +1 for decorator, this probably needs
                                        # the same thing doing as for the
                                        # error messages where the decorator
                                        # is scanned for, its not always +1!
                                        lidx = int(line) - (firstlineno + 1)
                                        source_line = src_code[lidx + 1]
                                        ln = fmt.format(cs['python'], col_span,
                                                        clean(source_line))
                                        new_lines.append(ln)
                                    return line, col

            for l in sliced_lines[1:]:

                # Drop LLVM debug call entries
                if dbg_value.match(l):
                    continue

                # if requested generate interleaving of markers or python from
                # metadata
                if _interleave.lineinfo or _interleave.python:
                    updated_lineinfo = metadata_interleave(l, new_lines)
                    if updated_lineinfo is not None:
                        cur_line, cur_col = updated_lineinfo

                # Highlight other LLVM features if requested, HTML BGCOLOR
                # property is set by this.
                if _highlight.incref and nrt_incref.search(l):
                    colour = cs['incref']
                elif _highlight.decref and nrt_decref.search(l):
                    colour = cs['decref']
                elif _highlight.meminfo and nrt_meminfo.search(l):
                    colour = cs['meminfo']
                elif _highlight.raises and ll_raise.search(l):
                    # search for raise as its more specific than exit
                    colour = cs['raise']
                elif _highlight.returns and ll_return.search(l):
                    colour = cs['return']
                elif _highlight.llvm_intrin_calls and ll_intrin_calls.search(l):
                    colour = cs['llvm_intrin_calls']
                elif _highlight.function_calls and ll_function_call.search(l):
                    colour = cs['function_calls']
                else:
                    colour = cs['default']

                # Use the default coloring as a flag to force printing if a
                # special token print was requested AND LLVM ir stripping is
                # required
                if colour is not cs['default'] or not strip_ir:
                    for x in wrap(clean(l)):
                        new_lines.append(fmt.format(colour, col_span, x))

            # add in the port line at the end of the block if it was present
            # (this was built right at the top of the parse)
            if port_line:
                new_lines.append('<tr>{}</tr>'.format(port_line))

            # If there was data, create a table, else don't!
            dat = ''.join(new_lines)
            if dat:
                tab = (('<table id="%s" BORDER="1" CELLBORDER="0" '
                       'CELLPADDING="0" CELLSPACING="0">%s</table>') % (idc,
                                                                        dat))
                label = '<{}>'.format(tab)
            else:
                label = ''

            # finally, add a replacement node for the original with a new marked
            # up label.
            f.node(name, label=label)

        # Parse the edge data
        if 'edges' in jzon: # might be a single block, no edges
            for edge in jzon['edges']:
                gvid = edge['_gvid']
                tp = edge.get('tailport', None)
                edge_ids[gvid] = (edge['head'], edge['tail'], tp)

        # Write in the edge wiring with respect to the new nodes:ports.
        for gvid, edge in edge_ids.items():
            tail = node_ids[edge[1]]
            head = node_ids[edge[0]]
            port = edge[2]
            if port is not None:
                tail += ':%s' % port
            f.edge(tail, head)

        # Add a key to the graph if requested.
        if show_key:
            key_tab = []
            for k, v in cs.items():
                key_tab.append(('<tr><td BGCOLOR="{}" BORDER="0" ALIGN="center"'
                                '>{}</td></tr>').format(v, k))
            # The first < and last > are DOT syntax, rest is DOT HTML.
            f.node("Key", label=('<<table BORDER="1" CELLBORDER="1" '
                    'CELLPADDING="2" CELLSPACING="1"><tr><td BORDER="0">'
                    'Key:</td></tr>{}</table>>').format(''.join(key_tab)))

        # Render if required
        if filename is not None or view is not None:
            f.render(filename=filename, view=view, format=render_format)

        # Else pipe out a SVG
        return f.pipe(format='svg')

    def display(self, filename=None, format='pdf', view=False):
        """
        Plot the CFG.  In IPython notebook, the return image object can be
        inlined.

        The *filename* option can be set to a specific path for the rendered
        output to write to.  If *view* option is True, the plot is opened by
        the system default application for the image format (PDF). *format* can
        be any valid format string accepted by graphviz, default is 'pdf'.
        """
        rawbyt = self.pretty_printer(filename=filename, view=view,
                                     render_format=format, **self.kwargs)
        return rawbyt.decode('utf-8')

    def _repr_svg_(self):
        return self.pretty_printer(**self.kwargs).decode('utf-8')

    def __repr__(self):
        return self.dot


class CodeLibrary(metaclass=ABCMeta):
    """
    An interface for bundling LLVM code together and compiling it.
    It is tied to a *codegen* instance (e.g. JITCPUCodegen) that will
    determine how the LLVM code is transformed and linked together.
    """

    _finalized = False
    _object_caching_enabled = False
    _disable_inspection = False

    def __init__(self, codegen: "CPUCodegen", name: str):
        self._codegen = codegen
        self._name = name
        ptc_name = f"{self.__class__.__name__}({self._name!r})"
        self._recorded_timings = PassTimingsCollection(ptc_name)
        # Track names of the dynamic globals
        self._dynamic_globals = []

    @property
    def has_dynamic_globals(self):
        self._ensure_finalized()
        return len(self._dynamic_globals) > 0

    @property
    def recorded_timings(self):
        return self._recorded_timings

    @property
    def codegen(self):
        """
        The codegen object owning this library.
        """
        return self._codegen

    @property
    def name(self):
        return self._name

    def __repr__(self):
        return "<Library %r at 0x%x>" % (self.name, id(self))

    def _raise_if_finalized(self):
        if self._finalized:
            raise RuntimeError("operation impossible on finalized object %r"
                               % (self,))

    def _ensure_finalized(self):
        if not self._finalized:
            self.finalize()

    def create_ir_module(self, name):
        """
        Create an LLVM IR module for use by this library.
        """
        self._raise_if_finalized()
        ir_module = self._codegen._create_empty_module(name)
        return ir_module

    @abstractmethod
    def add_linking_library(self, library):
        """
        Add a library for linking into this library, without losing
        the original library.
        """

    @abstractmethod
    def add_ir_module(self, ir_module):
        """
        Add an LLVM IR module's contents to this library.
        """

    @abstractmethod
    def finalize(self):
        """
        Finalize the library.  After this call, nothing can be added anymore.
        Finalization involves various stages of code optimization and
        linking.
        """

    @abstractmethod
    def get_function(self, name):
        """
        Return the function named ``name``.
        """

    @abstractmethod
    def get_llvm_str(self):
        """
        Get the human-readable form of the LLVM module.
        """

    @abstractmethod
    def get_asm_str(self):
        """
        Get the human-readable assembly.
        """

    #
    # Object cache hooks and serialization
    #

    def enable_object_caching(self):
        self._object_caching_enabled = True
        self._compiled_object = None
        self._compiled = False

    def _get_compiled_object(self):
        if not self._object_caching_enabled:
            raise ValueError("object caching not enabled in %s" % (self,))
        if self._compiled_object is None:
            raise RuntimeError("no compiled object yet for %s" % (self,))
        return self._compiled_object

    def _set_compiled_object(self, value):
        if not self._object_caching_enabled:
            raise ValueError("object caching not enabled in %s" % (self,))
        if self._compiled:
            raise ValueError("library already compiled: %s" % (self,))
        self._compiled_object = value
        self._disable_inspection = True


class CPUCodeLibrary(CodeLibrary):

    def __init__(self, codegen, name):
        super().__init__(codegen, name)
        self._linking_libraries = []   # maintain insertion order
        self._final_module = ll.parse_assembly(
            str(self._codegen._create_empty_module(self.name)))
        self._final_module.name = cgutils.normalize_ir_text(self.name)
        self._shared_module = None

    def _optimize_functions(self, ll_module):
        """
        Internal: run function-level optimizations inside *ll_module*.
        """
        # Enforce data layout to enable layout-specific optimizations
        ll_module.data_layout = self._codegen._data_layout
        with self._codegen._function_pass_manager(ll_module) as fpm:
            # Run function-level optimizations to reduce memory usage and improve
            # module-level optimization.
            for func in ll_module.functions:
                k = f"Function passes on {func.name!r}"
                with self._recorded_timings.record(k):
                    fpm.initialize()
                    fpm.run(func)
                    fpm.finalize()

    def _optimize_final_module(self):
        """
        Internal: optimize this library's final module.
        """
        cheap_name = "Module passes (cheap optimization for refprune)"
        with self._recorded_timings.record(cheap_name):
            # A cheaper optimisation pass is run first to try and get as many
            # refops into the same function as possible via inlining
            self._codegen._mpm_cheap.run(self._final_module)
        # Refop pruning is then run on the heavily inlined function
        if not config.LLVM_REFPRUNE_PASS:
            self._final_module = remove_redundant_nrt_refct(self._final_module)
        full_name = "Module passes (full optimization)"
        with self._recorded_timings.record(full_name):
            # The full optimisation suite is then run on the refop pruned IR
            self._codegen._mpm_full.run(self._final_module)

    def _get_module_for_linking(self):
        """
        Internal: get a LLVM module suitable for linking multiple times
        into another library.  Exported functions are made "linkonce_odr"
        to allow for multiple definitions, inlining, and removal of
        unused exports.

        See discussion in https://github.com/numba/numba/pull/890
        """
        self._ensure_finalized()
        if self._shared_module is not None:
            return self._shared_module
        mod = self._final_module
        to_fix = []
        nfuncs = 0
        for fn in mod.functions:
            nfuncs += 1
            if not fn.is_declaration and fn.linkage == ll.Linkage.external:
                to_fix.append(fn.name)
        if nfuncs == 0:
            # This is an issue which can occur if loading a module
            # from an object file and trying to link with it, so detect it
            # here to make debugging easier.
            raise RuntimeError("library unfit for linking: "
                               "no available functions in %s"
                               % (self,))
        if to_fix:
            mod = mod.clone()
            for name in to_fix:
                # NOTE: this will mark the symbol WEAK if serialized
                # to an ELF file
                mod.get_function(name).linkage = 'linkonce_odr'
        self._shared_module = mod
        return mod

    def add_linking_library(self, library):
        library._ensure_finalized()
        self._linking_libraries.append(library)

    def add_ir_module(self, ir_module):
        self._raise_if_finalized()
        assert isinstance(ir_module, llvmir.Module)
        ir = cgutils.normalize_ir_text(str(ir_module))
        ll_module = ll.parse_assembly(ir)
        ll_module.name = ir_module.name
        ll_module.verify()
        self.add_llvm_module(ll_module)

    def add_llvm_module(self, ll_module):
        self._optimize_functions(ll_module)
        # TODO: we shouldn't need to recreate the LLVM module object
        if not config.LLVM_REFPRUNE_PASS:
            ll_module = remove_redundant_nrt_refct(ll_module)
        self._final_module.link_in(ll_module)

    def finalize(self):
        require_global_compiler_lock()

        # Report any LLVM-related problems to the user
        self._codegen._check_llvm_bugs()

        self._raise_if_finalized()

        if config.DUMP_FUNC_OPT:
            dump("FUNCTION OPTIMIZED DUMP %s" % self.name,
                 self.get_llvm_str(), 'llvm')

        # Link libraries for shared code
        seen = set()
        for library in self._linking_libraries:
            if library not in seen:
                seen.add(library)
                self._final_module.link_in(
                    library._get_module_for_linking(), preserve=True,
                )

        # Optimize the module after all dependences are linked in above,
        # to allow for inlining.
        self._optimize_final_module()

        self._final_module.verify()
        self._finalize_final_module()

    def _finalize_dynamic_globals(self):
        # Scan for dynamic globals
        for gv in self._final_module.global_variables:
            if gv.name.startswith('numba.dynamic.globals'):
                self._dynamic_globals.append(gv.name)

    def _verify_declare_only_symbols(self):
        # Verify that no declare-only function compiled by numba.
        for fn in self._final_module.functions:
            # We will only check for symbol name starting with '_ZN5numba'
            if fn.is_declaration and fn.name.startswith('_ZN5numba'):
                msg = 'Symbol {} not linked properly'
                raise AssertionError(msg.format(fn.name))

    def _finalize_final_module(self):
        """
        Make the underlying LLVM module ready to use.
        """
        self._finalize_dynamic_globals()
        self._verify_declare_only_symbols()

        # Remember this on the module, for the object cache hooks
        self._final_module.__library = weakref.proxy(self)

        # It seems add_module() must be done only here and not before
        # linking in other modules, otherwise get_pointer_to_function()
        # could fail.
        cleanup = self._codegen._add_module(self._final_module)
        if cleanup:
            weakref.finalize(self, cleanup)
        self._finalize_specific()

        self._finalized = True

        if config.DUMP_OPTIMIZED:
            dump("OPTIMIZED DUMP %s" % self.name, self.get_llvm_str(), 'llvm')

        if config.DUMP_ASSEMBLY:
            dump("ASSEMBLY %s" % self.name, self.get_asm_str(), 'asm')

    def get_defined_functions(self):
        """
        Get all functions defined in the library.  The library must have
        been finalized.
        """
        mod = self._final_module
        for fn in mod.functions:
            if not fn.is_declaration:
                yield fn

    def get_function(self, name):
        return self._final_module.get_function(name)

    def _sentry_cache_disable_inspection(self):
        if self._disable_inspection:
            warnings.warn('Inspection disabled for cached code. '
                          'Invalid result is returned.')

    def get_llvm_str(self):
        self._sentry_cache_disable_inspection()
        return str(self._final_module)

    def get_asm_str(self):
        self._sentry_cache_disable_inspection()
        return str(self._codegen._tm.emit_assembly(self._final_module))

    def get_function_cfg(self, name, py_func=None, **kwargs):
        """
        Get control-flow graph of the LLVM function
        """
        self._sentry_cache_disable_inspection()
        return _CFG(self, name, py_func, **kwargs)

    def get_disasm_cfg(self, mangled_name):
        """
        Get the CFG of the disassembly of the ELF object at symbol mangled_name.

        Requires python package: r2pipe
        Requires radare2 binary on $PATH.
        Notebook rendering requires python package: graphviz
        Optionally requires a compiler toolchain (via pycc) to link the ELF to
        get better disassembly results.
        """
        elf = self._get_compiled_object()
        return disassemble_elf_to_cfg(elf, mangled_name)

    @classmethod
    def _dump_elf(cls, buf):
        """
        Dump the symbol table of an ELF file.
        Needs pyelftools (https://github.com/eliben/pyelftools)
        """
        from elftools.elf.elffile import ELFFile
        from elftools.elf import descriptions
        from io import BytesIO
        f = ELFFile(BytesIO(buf))
        print("ELF file:")
        for sec in f.iter_sections():
            if sec['sh_type'] == 'SHT_SYMTAB':
                symbols = sorted(sec.iter_symbols(), key=lambda sym: sym.name)
                print("    symbols:")
                for sym in symbols:
                    if not sym.name:
                        continue
                    print("    - %r: size=%d, value=0x%x, type=%s, bind=%s"
                          % (sym.name.decode(),
                             sym['st_size'],
                             sym['st_value'],
                             descriptions.describe_symbol_type(sym['st_info']['type']),
                             descriptions.describe_symbol_bind(sym['st_info']['bind']),
                             ))
        print()

    @classmethod
    def _object_compiled_hook(cls, ll_module, buf):
        """
        `ll_module` was compiled into object code `buf`.
        """
        try:
            self = ll_module.__library
        except AttributeError:
            return
        if self._object_caching_enabled:
            self._compiled = True
            self._compiled_object = buf

    @classmethod
    def _object_getbuffer_hook(cls, ll_module):
        """
        Return a cached object code for `ll_module`.
        """
        try:
            self = ll_module.__library
        except AttributeError:
            return
        if self._object_caching_enabled and self._compiled_object:
            buf = self._compiled_object
            self._compiled_object = None
            return buf

    def serialize_using_bitcode(self):
        """
        Serialize this library using its bitcode as the cached representation.
        """
        self._ensure_finalized()
        return (self.name, 'bitcode', self._final_module.as_bitcode())

    def serialize_using_object_code(self):
        """
        Serialize this library using its object code as the cached
        representation.  We also include its bitcode for further inlining
        with other libraries.
        """
        self._ensure_finalized()
        data = (self._get_compiled_object(),
                self._get_module_for_linking().as_bitcode())
        return (self.name, 'object', data)

    @classmethod
    def _unserialize(cls, codegen, state):
        name, kind, data = state
        self = codegen.create_library(name)
        assert isinstance(self, cls)
        if kind == 'bitcode':
            # No need to re-run optimizations, just make the module ready
            self._final_module = ll.parse_bitcode(data)
            self._finalize_final_module()
            return self
        elif kind == 'object':
            object_code, shared_bitcode = data
            self.enable_object_caching()
            self._set_compiled_object(object_code)
            self._shared_module = ll.parse_bitcode(shared_bitcode)
            self._finalize_final_module()
            # Load symbols from cache
            self._codegen._engine._load_defined_symbols(self._shared_module)
            return self
        else:
            raise ValueError("unsupported serialization kind %r" % (kind,))


class AOTCodeLibrary(CPUCodeLibrary):

    def emit_native_object(self):
        """
        Return this library as a native object (a bytestring) -- for example
        ELF under Linux.

        This function implicitly calls .finalize().
        """
        self._ensure_finalized()
        return self._codegen._tm.emit_object(self._final_module)

    def emit_bitcode(self):
        """
        Return this library as LLVM bitcode (a bytestring).

        This function implicitly calls .finalize().
        """
        self._ensure_finalized()
        return self._final_module.as_bitcode()

    def _finalize_specific(self):
        pass


class JITCodeLibrary(CPUCodeLibrary):

    def get_pointer_to_function(self, name):
        """
        Generate native code for function named *name* and return a pointer
        to the start of the function (as an integer).

        This function implicitly calls .finalize().

        Returns
        -------
        pointer : int
            - zero (null) if no symbol of *name* is defined by this code
              library.
            - non-zero if the symbol is defined.
        """
        self._ensure_finalized()
        ee = self._codegen._engine
        if not ee.is_symbol_defined(name):
            return 0
        else:
            return self._codegen._engine.get_function_address(name)

    def _finalize_specific(self):
        self._codegen._scan_and_fix_unresolved_refs(self._final_module)
        with self._recorded_timings.record("Finalize object"):
            self._codegen._engine.finalize_object()


class RuntimeLinker(object):
    """
    For tracking unresolved symbols generated at runtime due to recursion.
    """
    PREFIX = '.numba.unresolved$'

    def __init__(self):
        self._unresolved = utils.UniqueDict()
        self._defined = set()
        self._resolved = []

    def scan_unresolved_symbols(self, module, engine):
        """
        Scan and track all unresolved external symbols in the module and
        allocate memory for it.
        """
        prefix = self.PREFIX

        for gv in module.global_variables:
            if gv.name.startswith(prefix):
                sym = gv.name[len(prefix):]
                # Avoid remapping to existing GV
                if engine.is_symbol_defined(gv.name):
                    continue
                # Allocate a memory space for the pointer
                abortfn = rtsys.library.get_pointer_to_function("nrt_unresolved_abort")
                ptr = ctypes.c_void_p(abortfn)
                engine.add_global_mapping(gv, ctypes.addressof(ptr))
                self._unresolved[sym] = ptr

    def scan_defined_symbols(self, module):
        """
        Scan and track all defined symbols.
        """
        for fn in module.functions:
            if not fn.is_declaration:
                self._defined.add(fn.name)

    def resolve(self, engine):
        """
        Fix unresolved symbols if they are defined.
        """
        # An iterator to get all unresolved but available symbols
        pending = [name for name in self._unresolved if name in self._defined]
        # Resolve pending symbols
        for name in pending:
            # Get runtime address
            fnptr = engine.get_function_address(name)
            # Fix all usage
            ptr = self._unresolved[name]
            ptr.value = fnptr
            self._resolved.append((name, ptr))   # keep ptr alive
            # Delete resolved
            del self._unresolved[name]

def _proxy(old):
    @functools.wraps(old)
    def wrapper(self, *args, **kwargs):
        return old(self._ee, *args, **kwargs)
    return wrapper


class JitEngine(object):
    """Wraps an ExecutionEngine to provide custom symbol tracking.
    Since the symbol tracking is incomplete  (doesn't consider
    loaded code object), we are not putting it in llvmlite.
    """
    def __init__(self, ee):
        self._ee = ee
        # Track symbol defined via codegen'd Module
        # but not any cached object.
        # NOTE: `llvm::ExecutionEngine` will catch duplicated symbols and
        # we are not going to protect against that.  A proper duplicated
        # symbol detection will need a more logic to check for the linkage
        # (e.g. like `weak` linkage symbol can override).   This
        # `_defined_symbols` set will be just enough to tell if a symbol
        # exists and will not cause the `EE` symbol lookup to `exit(1)`
        # when symbol-not-found.
        self._defined_symbols = set()

    def is_symbol_defined(self, name):
        """Is the symbol defined in this session?
        """
        return name in self._defined_symbols

    def _load_defined_symbols(self, mod):
        """Extract symbols from the module
        """
        for gsets in (mod.functions, mod.global_variables):
            self._defined_symbols |= {gv.name for gv in gsets
                                      if not gv.is_declaration}

    def add_module(self, module):
        """Override ExecutionEngine.add_module
        to keep info about defined symbols.
        """
        self._load_defined_symbols(module)
        return self._ee.add_module(module)

    def add_global_mapping(self, gv, addr):
        """Override ExecutionEngine.add_global_mapping
        to keep info about defined symbols.
        """
        self._defined_symbols.add(gv.name)
        return self._ee.add_global_mapping(gv, addr)

    #
    # The remaining methods are re-export of the ExecutionEngine APIs
    #
    set_object_cache = _proxy(ll.ExecutionEngine.set_object_cache)
    finalize_object = _proxy(ll.ExecutionEngine.finalize_object)
    get_function_address = _proxy(ll.ExecutionEngine.get_function_address)
    get_global_value_address = _proxy(
        ll.ExecutionEngine.get_global_value_address
        )


class Codegen(metaclass=ABCMeta):
    """
    Base Codegen class. It is expected that subclasses set the class attribute
    ``_library_class``, indicating the CodeLibrary class for the target.

    Subclasses should also initialize:

    ``self._data_layout``: the data layout for the target.
    ``self._target_data``: the binding layer ``TargetData`` for the target.
    """

    @abstractmethod
    def _create_empty_module(self, name):
        """
        Create a new empty module suitable for the target.
        """

    @abstractmethod
    def _add_module(self, module):
        """
        Add a module to the execution engine. Ownership of the module is
        transferred to the engine.
        """

    @property
    def target_data(self):
        """
        The LLVM "target data" object for this codegen instance.
        """
        return self._target_data

    def create_library(self, name, **kwargs):
        """
        Create a :class:`CodeLibrary` object for use with this codegen
        instance.
        """
        return self._library_class(self, name, **kwargs)

    def unserialize_library(self, serialized):
        return self._library_class._unserialize(self, serialized)


class CPUCodegen(Codegen):

    def __init__(self, module_name):
        initialize_llvm()

        self._data_layout = None
        self._llvm_module = ll.parse_assembly(
            str(self._create_empty_module(module_name)))
        self._llvm_module.name = "global_codegen_module"
        self._rtlinker = RuntimeLinker()
        self._init(self._llvm_module)

    def _init(self, llvm_module):
        assert list(llvm_module.global_variables) == [], "Module isn't empty"

        target = ll.Target.from_triple(ll.get_process_triple())
        tm_options = dict(opt=config.OPT)
        self._tm_features = self._customize_tm_features()
        self._customize_tm_options(tm_options)
        tm = target.create_target_machine(**tm_options)
        use_lmm = config.USE_LLVMLITE_MEMORY_MANAGER
        engine = ll.create_mcjit_compiler(llvm_module, tm, use_lmm=use_lmm)

        if config.ENABLE_PROFILING:
            engine.enable_jit_events()

        self._tm = tm
        self._engine = JitEngine(engine)
        self._target_data = engine.target_data
        self._data_layout = str(self._target_data)

        if config.OPT.is_opt_max:
            # If the OPT level is set to 'max' then the user is requesting that
            # compilation time is traded for potential performance gain. This
            # currently manifests as running the "cheap" pass at -O3
            # optimisation level with loop-vectorization enabled. There's no
            # guarantee that this will increase runtime performance, it may
            # detriment it, this is here to give the user an easily accessible
            # option to try.
            loopvect = True
            opt_level = 3
        else:
            # The default behaviour is to do an opt=0 pass to try and inline as
            # much as possible with the cheapest cost of doing so. This is so
            # that the ref-op pruner pass that runs after the cheap pass will
            # have the largest possible scope for working on pruning references.
            loopvect = False
            opt_level = 0

        self._mpm_cheap = self._module_pass_manager(loop_vectorize=loopvect,
                                                    slp_vectorize=False,
                                                    opt=opt_level,
                                                    cost="cheap")

        self._mpm_full = self._module_pass_manager()

        self._engine.set_object_cache(self._library_class._object_compiled_hook,
                                      self._library_class._object_getbuffer_hook)

    def _create_empty_module(self, name):
        ir_module = llvmir.Module(cgutils.normalize_ir_text(name))
        ir_module.triple = ll.get_process_triple()
        if self._data_layout:
            ir_module.data_layout = self._data_layout
        return ir_module

    def _module_pass_manager(self, **kwargs):
        pm = ll.create_module_pass_manager()
        pm.add_target_library_info(ll.get_process_triple())
        self._tm.add_analysis_passes(pm)
        cost = kwargs.pop("cost", None)
        with self._pass_manager_builder(**kwargs) as pmb:
            pmb.populate(pm)
        # If config.OPT==0 do not include these extra passes to help with
        # vectorization.
        if cost is not None and cost == "cheap" and config.OPT != 0:
            # This knocks loops into rotated form early to reduce the likelihood
            # of vectorization failing due to unknown PHI nodes.
            pm.add_loop_rotate_pass()
            # These passes are required to get SVML to vectorize tests
            pm.add_instruction_combining_pass()
            pm.add_jump_threading_pass()

        if config.LLVM_REFPRUNE_PASS:
            pm.add_refprune_pass(_parse_refprune_flags())
        return pm

    def _function_pass_manager(self, llvm_module, **kwargs):
        pm = ll.create_function_pass_manager(llvm_module)
        pm.add_target_library_info(llvm_module.triple)
        self._tm.add_analysis_passes(pm)
        with self._pass_manager_builder(**kwargs) as pmb:
            pmb.populate(pm)
        if config.LLVM_REFPRUNE_PASS:
            pm.add_refprune_pass(_parse_refprune_flags())
        return pm

    def _pass_manager_builder(self, **kwargs):
        """
        Create a PassManagerBuilder.

        Note: a PassManagerBuilder seems good only for one use, so you
        should call this method each time you want to populate a module
        or function pass manager.  Otherwise some optimizations will be
        missed...
        """
        opt_level = kwargs.pop('opt', config.OPT)
        loop_vectorize = kwargs.pop('loop_vectorize', config.LOOP_VECTORIZE)
        slp_vectorize = kwargs.pop('slp_vectorize', config.SLP_VECTORIZE)

        pmb = create_pass_manager_builder(opt=opt_level,
                                          loop_vectorize=loop_vectorize,
                                          slp_vectorize=slp_vectorize,
                                          **kwargs)

        return pmb

    def _check_llvm_bugs(self):
        """
        Guard against some well-known LLVM bug(s).
        """
        # Check the locale bug at https://github.com/numba/numba/issues/1569
        # Note we can't cache the result as locale settings can change
        # across a process's lifetime.  Also, for this same reason,
        # the check here is a mere heuristic (there may be a race condition
        # between now and actually compiling IR).
        ir = """
            define double @func()
            {
                ret double 1.23e+01
            }
            """
        mod = ll.parse_assembly(ir)
        ir_out = str(mod)
        if "12.3" in ir_out or "1.23" in ir_out:
            # Everything ok
            return
        if "1.0" in ir_out:
            loc = locale.getlocale()
            raise RuntimeError(
                "LLVM will produce incorrect floating-point code "
                "in the current locale %s.\nPlease read "
                "https://numba.readthedocs.io/en/stable/user/faq.html#llvm-locale-bug "
                "for more information."
                % (loc,))
        raise AssertionError("Unexpected IR:\n%s\n" % (ir_out,))

    def magic_tuple(self):
        """
        Return a tuple unambiguously describing the codegen behaviour.
        """
        return (self._llvm_module.triple, self._get_host_cpu_name(),
                self._tm_features)

    def _scan_and_fix_unresolved_refs(self, module):
        self._rtlinker.scan_unresolved_symbols(module, self._engine)
        self._rtlinker.scan_defined_symbols(module)
        self._rtlinker.resolve(self._engine)

    def insert_unresolved_ref(self, builder, fnty, name):
        voidptr = llvmir.IntType(8).as_pointer()
        ptrname = self._rtlinker.PREFIX + name
        llvm_mod = builder.module
        try:
            fnptr = llvm_mod.get_global(ptrname)
        except KeyError:
            # Not defined?
            fnptr = llvmir.GlobalVariable(llvm_mod, voidptr, name=ptrname)
            fnptr.linkage = 'external'
        return builder.bitcast(builder.load(fnptr), fnty.as_pointer())

    def _get_host_cpu_name(self):
        return (ll.get_host_cpu_name()
                if config.CPU_NAME is None
                else config.CPU_NAME)

    def _get_host_cpu_features(self):
        if config.CPU_FEATURES is not None:
            return config.CPU_FEATURES
        return get_host_cpu_features()


class AOTCPUCodegen(CPUCodegen):
    """
    A codegen implementation suitable for Ahead-Of-Time compilation
    (e.g. generation of object files).
    """

    _library_class = AOTCodeLibrary

    def __init__(self, module_name, cpu_name=None):
        # By default, use generic cpu model for the arch
        self._cpu_name = cpu_name or ''
        CPUCodegen.__init__(self, module_name)

    def _customize_tm_options(self, options):
        cpu_name = self._cpu_name
        if cpu_name == 'host':
            cpu_name = self._get_host_cpu_name()
        options['cpu'] = cpu_name
        options['reloc'] = 'pic'
        options['codemodel'] = 'default'
        options['features'] = self._tm_features

    def _customize_tm_features(self):
        # ISA features are selected according to the requested CPU model
        # in _customize_tm_options()
        return ''

    def _add_module(self, module):
        pass


class JITCPUCodegen(CPUCodegen):
    """
    A codegen implementation suitable for Just-In-Time compilation.
    """

    _library_class = JITCodeLibrary

    def _customize_tm_options(self, options):
        # As long as we don't want to ship the code to another machine,
        # we can specialize for this CPU.
        options['cpu'] = self._get_host_cpu_name()
        # LLVM 7 change: # https://reviews.llvm.org/D47211#inline-425406
        # JIT needs static relocation on x86*
        # native target is already initialized from base class __init__
        arch = ll.Target.from_default_triple().name
        if arch.startswith('x86'): # one of x86 or x86_64
            reloc_model = 'static'
        elif arch.startswith('ppc'):
            reloc_model = 'pic'
        else:
            reloc_model = 'default'
        options['reloc'] = reloc_model
        options['codemodel'] = 'jitdefault'

        # Set feature attributes (such as ISA extensions)
        # This overrides default feature selection by CPU model above
        options['features'] = self._tm_features

        # Deal with optional argument to ll.Target.create_target_machine
        sig = utils.pysignature(ll.Target.create_target_machine)
        if 'jit' in sig.parameters:
            # Mark that this is making a JIT engine
            options['jit'] = True

    def _customize_tm_features(self):
        # For JIT target, we will use LLVM to get the feature map
        return self._get_host_cpu_features()

    def _add_module(self, module):
        self._engine.add_module(module)
        # XXX: disabling remove module due to MCJIT engine leakage in
        #      removeModule.  The removeModule causes consistent access
        #      violation with certain test combinations.
        # # Early bind the engine method to avoid keeping a reference to self.
        # return functools.partial(self._engine.remove_module, module)

    def set_env(self, env_name, env):
        """Set the environment address.

        Update the GlobalVariable named *env_name* to the address of *env*.
        """
        gvaddr = self._engine.get_global_value_address(env_name)
        envptr = (ctypes.c_void_p * 1).from_address(gvaddr)
        envptr[0] = ctypes.c_void_p(id(env))


def initialize_llvm():
    """Safe to use multiple times.
    """
    ll.initialize()
    ll.initialize_native_target()
    ll.initialize_native_asmprinter()


def get_host_cpu_features():
    """Get host CPU features using LLVM.

    The features may be modified due to user setting.
    See numba.config.ENABLE_AVX.
    """
    try:
        features = ll.get_host_cpu_features()
    except RuntimeError:
        return ''
    else:
        if not config.ENABLE_AVX:
            # Disable all features with name starting with 'avx'
            for k in features:
                if k.startswith('avx'):
                    features[k] = False

        # Set feature attributes
        return features.flatten()
