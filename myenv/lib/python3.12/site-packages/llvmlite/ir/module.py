import collections

from llvmlite.ir import context, values, types, _utils


class Module(object):
    def __init__(self, name='', context=context.global_context):
        self.context = context
        self.name = name   # name is for debugging/informational
        self.data_layout = ""
        self.scope = _utils.NameScope()
        self.triple = 'unknown-unknown-unknown'
        self.globals = collections.OrderedDict()
        # Innamed metadata nodes.
        self.metadata = []
        # Named metadata nodes
        self.namedmetadata = {}
        # Cache for metadata node deduplication
        self._metadatacache = {}

    def _fix_metadata_operands(self, operands):
        fixed_ops = []
        for op in operands:
            if op is None:
                # A literal None creates a null metadata value
                op = types.MetaDataType()(None)
            elif isinstance(op, str):
                # A literal string creates a metadata string value
                op = values.MetaDataString(self, op)
            elif isinstance(op, (list, tuple)):
                # A sequence creates a metadata node reference
                op = self.add_metadata(op)
            fixed_ops.append(op)
        return fixed_ops

    def _fix_di_operands(self, operands):
        fixed_ops = []
        for name, op in operands:
            if isinstance(op, (list, tuple)):
                # A sequence creates a metadata node reference
                op = self.add_metadata(op)
            fixed_ops.append((name, op))
        return fixed_ops

    def add_metadata(self, operands):
        """
        Add an unnamed metadata to the module with the given *operands*
        (a sequence of values) or return a previous equivalent metadata.
        A MDValue instance is returned, it can then be associated to
        e.g. an instruction.
        """
        if not isinstance(operands, (list, tuple)):
            raise TypeError("expected a list or tuple of metadata values, "
                            "got %r" % (operands,))
        operands = self._fix_metadata_operands(operands)
        key = tuple(operands)
        if key not in self._metadatacache:
            n = len(self.metadata)
            md = values.MDValue(self, operands, name=str(n))
            self._metadatacache[key] = md
        else:
            md = self._metadatacache[key]
        return md

    def add_debug_info(self, kind, operands, is_distinct=False):
        """
        Add debug information metadata to the module with the given
        *operands* (a dict of values with string keys) or return
        a previous equivalent metadata.  *kind* is a string of the
        debug information kind (e.g. "DICompileUnit").

        A DIValue instance is returned, it can then be associated to e.g.
        an instruction.
        """
        operands = tuple(sorted(self._fix_di_operands(operands.items())))
        key = (kind, operands, is_distinct)
        if key not in self._metadatacache:
            n = len(self.metadata)
            di = values.DIValue(self, is_distinct, kind, operands, name=str(n))
            self._metadatacache[key] = di
        else:
            di = self._metadatacache[key]
        return di

    def add_named_metadata(self, name, element=None):
        """
        Add a named metadata node to the module, if it doesn't exist,
        or return the existing node.
        If *element* is given, it will append a new element to
        the named metadata node.  If *element* is a sequence of values
        (rather than a metadata value), a new unnamed node will first be
        created.

        Example::
            module.add_named_metadata("llvm.ident", ["llvmlite/1.0"])
        """
        if name in self.namedmetadata:
            nmd = self.namedmetadata[name]
        else:
            nmd = self.namedmetadata[name] = values.NamedMetaData(self)
        if element is not None:
            if not isinstance(element, values.Value):
                element = self.add_metadata(element)
            if not isinstance(element.type, types.MetaDataType):
                raise TypeError("wrong type for metadata element: got %r"
                                % (element,))
            nmd.add(element)
        return nmd

    def get_named_metadata(self, name):
        """
        Return the metadata node with the given *name*.  KeyError is raised
        if no such node exists (contrast with add_named_metadata()).
        """
        return self.namedmetadata[name]

    @property
    def functions(self):
        """
        A list of functions declared or defined in this module.
        """
        return [v for v in self.globals.values()
                if isinstance(v, values.Function)]

    @property
    def global_values(self):
        """
        An iterable of global values in this module.
        """
        return self.globals.values()

    def get_global(self, name):
        """
        Get a global value by name.
        """
        return self.globals[name]

    def add_global(self, globalvalue):
        """
        Add a new global value.
        """
        assert globalvalue.name not in self.globals
        self.globals[globalvalue.name] = globalvalue

    def get_unique_name(self, name=''):
        """
        Get a unique global name with the following *name* hint.
        """
        return self.scope.deduplicate(name)

    def declare_intrinsic(self, intrinsic, tys=(), fnty=None):
        def _error():
            raise NotImplementedError("unknown intrinsic %r with %d types"
                                      % (intrinsic, len(tys)))

        if intrinsic in {'llvm.cttz', 'llvm.ctlz', 'llvm.fma'}:
            suffixes = [tys[0].intrinsic_name]
        else:
            suffixes = [t.intrinsic_name for t in tys]
        name = '.'.join([intrinsic] + suffixes)
        if name in self.globals:
            return self.globals[name]

        if fnty is not None:
            # General case: function type is given
            pass
        # Compute function type if omitted for common cases
        elif len(tys) == 0 and intrinsic == 'llvm.assume':
            fnty = types.FunctionType(types.VoidType(), [types.IntType(1)])
        elif len(tys) == 1:
            if intrinsic == 'llvm.powi':
                fnty = types.FunctionType(tys[0], [tys[0], types.IntType(32)])
            elif intrinsic == 'llvm.pow':
                fnty = types.FunctionType(tys[0], tys * 2)
            elif intrinsic == 'llvm.convert.from.fp16':
                fnty = types.FunctionType(tys[0], [types.IntType(16)])
            elif intrinsic == 'llvm.convert.to.fp16':
                fnty = types.FunctionType(types.IntType(16), tys)
            else:
                fnty = types.FunctionType(tys[0], tys)
        elif len(tys) == 2:
            if intrinsic == 'llvm.memset':
                tys = [tys[0], types.IntType(8), tys[1],
                       types.IntType(1)]
                fnty = types.FunctionType(types.VoidType(), tys)
            elif intrinsic in {'llvm.cttz', 'llvm.ctlz'}:
                tys = [tys[0], types.IntType(1)]
                fnty = types.FunctionType(tys[0], tys)
            else:
                _error()
        elif len(tys) == 3:
            if intrinsic in ('llvm.memcpy', 'llvm.memmove'):
                tys = tys + [types.IntType(1)]
                fnty = types.FunctionType(types.VoidType(), tys)
            elif intrinsic == 'llvm.fma':
                tys = [tys[0]] * 3
                fnty = types.FunctionType(tys[0], tys)
            else:
                _error()
        else:
            _error()
        return values.Function(self, fnty, name=name)

    def get_identified_types(self):
        return self.context.identified_types

    def _get_body_lines(self):
        # Type declarations
        lines = [it.get_declaration()
                 for it in self.get_identified_types().values()]
        # Global values (including function definitions)
        lines += [str(v) for v in self.globals.values()]
        return lines

    def _get_metadata_lines(self):
        mdbuf = []
        for k, v in self.namedmetadata.items():
            mdbuf.append("!{name} = !{{ {operands} }}".format(
                name=k, operands=', '.join(i.get_reference()
                                           for i in v.operands)))
        for md in self.metadata:
            mdbuf.append(str(md))
        return mdbuf

    def _stringify_body(self):
        # For testing
        return "\n".join(self._get_body_lines())

    def _stringify_metadata(self):
        # For testing
        return "\n".join(self._get_metadata_lines())

    def __repr__(self):
        lines = []
        # Header
        lines += [
            '; ModuleID = "%s"' % (self.name,),
            'target triple = "%s"' % (self.triple,),
            'target datalayout = "%s"' % (self.data_layout,),
            '']
        # Body
        lines += self._get_body_lines()
        # Metadata
        lines += self._get_metadata_lines()

        return "\n".join(lines)
