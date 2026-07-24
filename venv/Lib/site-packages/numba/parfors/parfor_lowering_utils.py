from collections import namedtuple

from numba.core import types, ir
from numba.core.typing import signature


_CallableNode = namedtuple("BoundFunc", ["func", "sig"])


class ParforLoweringBuilder:
    """Helper class for building Numba-IR and lowering for Parfor.
    """
    def __init__(self, lowerer, scope, loc):
        self._lowerer = lowerer
        self._scope = scope
        self._loc = loc

    @property
    def _context(self):
        return self._lowerer.context

    @property
    def _typingctx(self):
        return self._context.typing_context

    @property
    def _typemap(self):
        return self._lowerer.fndesc.typemap

    @property
    def _calltypes(self):
        return self._lowerer.fndesc.calltypes

    def bind_global_function(self, fobj, ftype, args, kws=None):
        """Binds a global function to a variable.

        Parameters
        ----------
        fobj : object
            The function to be bound.
        ftype : types.Type
        args : Sequence[types.Type]
        kws : Mapping[str, types.Type]

        Returns
        -------
        callable: _CallableNode
        """
        if kws is None:
            kws = {}
        loc = self._loc
        varname = f"{fobj.__name__}_func"
        gvname = f"{fobj.__name__}"

        func_sig = self._typingctx.resolve_function_type(ftype, args, kws)
        func_var = self.assign(
            rhs=ir.Global(gvname, fobj, loc=loc), typ=ftype, name=varname
        )
        return _CallableNode(func=func_var, sig=func_sig)

    def make_const_variable(self, cval, typ, name="pf_const") -> ir.Var:
        """Makes a constant variable

        Parameters
        ----------
        cval : object
            The constant value
        typ : types.Type
            type of the value
        name : str
            variable name to store to

        Returns
        -------
        res : ir.Var
        """
        return self.assign(
            rhs=ir.Const(cval, loc=self._loc), typ=typ, name=name
        )

    def make_tuple_variable(self, varlist, name="pf_tuple") -> ir.Var:
        """Makes a tuple variable

        Parameters
        ----------
        varlist : Sequence[ir.Var]
            Variables containing the values to be stored.
        name : str
            variable name to store to

        Returns
        -------
        res : ir.Var
        """
        loc = self._loc
        vartys = [self._typemap[x.name] for x in varlist]
        tupty = types.Tuple.from_types(vartys)
        return self.assign(
            rhs=ir.Expr.build_tuple(varlist, loc), typ=tupty, name=name
        )

    def assign(self, rhs, typ, name="pf_assign") -> ir.Var:
        """Assign a value to a new variable

        Parameters
        ----------
        rhs : object
            The value
        typ : types.Type
            type of the value
        name : str
            variable name to store to

        Returns
        -------
        res : ir.Var
        """
        loc = self._loc
        var = self._scope.redefine(name, loc)
        self._typemap[var.name] = typ
        assign = ir.Assign(rhs, var, loc)
        self._lowerer.lower_inst(assign)
        return var

    def assign_inplace(self, rhs, typ, name) -> ir.Var:
        """Assign a value to a new variable or inplace if it already exist

        Parameters
        ----------
        rhs : object
            The value
        typ : types.Type
            type of the value
        name : str
            variable name to store to

        Returns
        -------
        res : ir.Var
        """
        loc = self._loc
        var = ir.Var(self._scope, name, loc)
        assign = ir.Assign(rhs, var, loc)
        self._typemap.setdefault(var.name, typ)
        self._lowerer.lower_inst(assign)
        return var

    def call(self, callable_node, args, kws=None) -> ir.Expr:
        """Call a bound callable

        Parameters
        ----------
        callable_node : _CallableNode
            The callee
        args : Sequence[ir.Var]
        kws : Mapping[str, ir.Var]

        Returns
        -------
        res : ir.Expr
            The expression node for the return value of the call
        """
        if kws is None:
            kws = {}
        call = ir.Expr.call(callable_node.func, args, kws, loc=self._loc)
        self._calltypes[call] = callable_node.sig
        return call

    def setitem(self, obj, index, val) -> ir.SetItem:
        """Makes a setitem call

        Parameters
        ----------
        obj : ir.Var
            the object being indexed
        index : ir.Var
            the index
        val : ir.Var
            the value to be stored

        Returns
        -------
        res : ir.SetItem
        """
        loc = self._loc
        tm = self._typemap
        setitem = ir.SetItem(obj, index, val, loc=loc)
        self._lowerer.fndesc.calltypes[setitem] = signature(
            types.none, tm[obj.name], tm[index.name], tm[val.name]
        )
        self._lowerer.lower_inst(setitem)
        return setitem

    def getitem(self, obj, index, typ) -> ir.Expr:
        """Makes a getitem call

        Parameters
        ----------
        obj : ir.Var
            the object being indexed
        index : ir.Var
            the index
        val : ir.Var
            the ty

        Returns
        -------
        res : ir.Expr
            the retrieved value
        """
        tm = self._typemap
        getitem = ir.Expr.getitem(obj, index, loc=self._loc)
        self._lowerer.fndesc.calltypes[getitem] = signature(
            typ, tm[obj.name], tm[index.name],
        )
        return getitem
