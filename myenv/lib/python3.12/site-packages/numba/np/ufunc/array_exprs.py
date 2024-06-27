import ast
from collections import defaultdict, OrderedDict
import contextlib
import sys
from types import SimpleNamespace

import numpy as np
import operator

from numba.core import types, targetconfig, ir, rewrites, compiler
from numba.core.typing import npydecl
from numba.np.ufunc.dufunc import DUFunc


def _is_ufunc(func):
    return isinstance(func, (np.ufunc, DUFunc))


@rewrites.register_rewrite('after-inference')
class RewriteArrayExprs(rewrites.Rewrite):
    '''The RewriteArrayExprs class is responsible for finding array
    expressions in Numba intermediate representation code, and
    rewriting those expressions to a single operation that will expand
    into something similar to a ufunc call.
    '''
    def __init__(self, state, *args, **kws):
        super(RewriteArrayExprs, self).__init__(state, *args, **kws)
        # Install a lowering hook if we are using this rewrite.
        special_ops = state.targetctx.special_ops
        if 'arrayexpr' not in special_ops:
            special_ops['arrayexpr'] = _lower_array_expr

    def match(self, func_ir, block, typemap, calltypes):
        """
        Using typing and a basic block, search the basic block for array
        expressions.
        Return True when one or more matches were found, False otherwise.
        """
        # We can trivially reject everything if there are no
        # calls in the type results.
        if len(calltypes) == 0:
            return False

        self.crnt_block = block
        self.typemap = typemap
        # { variable name: IR assignment (of a function call or operator) }
        self.array_assigns = OrderedDict()
        # { variable name: IR assignment (of a constant) }
        self.const_assigns = {}

        assignments = block.find_insts(ir.Assign)
        for instr in assignments:
            target_name = instr.target.name
            expr = instr.value
            # Does it assign an expression to an array variable?
            if (isinstance(expr, ir.Expr) and
                isinstance(typemap.get(target_name, None), types.Array)):
                self._match_array_expr(instr, expr, target_name)
            elif isinstance(expr, ir.Const):
                # Track constants since we might need them for an
                # array expression.
                self.const_assigns[target_name] = expr

        return len(self.array_assigns) > 0

    def _match_array_expr(self, instr, expr, target_name):
        """
        Find whether the given assignment (*instr*) of an expression (*expr*)
        to variable *target_name* is an array expression.
        """
        # We've matched a subexpression assignment to an
        # array variable.  Now see if the expression is an
        # array expression.
        expr_op = expr.op
        array_assigns = self.array_assigns

        if ((expr_op in ('unary', 'binop')) and (
                expr.fn in npydecl.supported_array_operators)):
            # It is an array operator that maps to a ufunc.
            # check that all args have internal types
            if all(self.typemap[var.name].is_internal
                   for var in expr.list_vars()):
                array_assigns[target_name] = instr

        elif ((expr_op == 'call') and (expr.func.name in self.typemap)):
            # It could be a match for a known ufunc call.
            func_type = self.typemap[expr.func.name]
            if isinstance(func_type, types.Function):
                func_key = func_type.typing_key
                if _is_ufunc(func_key):
                    # If so, check whether an explicit output is passed.
                    if not self._has_explicit_output(expr, func_key):
                        # If not, match it as a (sub)expression.
                        array_assigns[target_name] = instr

    def _has_explicit_output(self, expr, func):
        """
        Return whether the *expr* call to *func* (a ufunc) features an
        explicit output argument.
        """
        nargs = len(expr.args) + len(expr.kws)
        if expr.vararg is not None:
            # XXX *args unsupported here, assume there may be an explicit
            # output
            return True
        return nargs > func.nin

    def _get_array_operator(self, ir_expr):
        ir_op = ir_expr.op
        if ir_op in ('unary', 'binop'):
            return ir_expr.fn
        elif ir_op == 'call':
            return self.typemap[ir_expr.func.name].typing_key
        raise NotImplementedError(
            "Don't know how to find the operator for '{0}' expressions.".format(
                ir_op))

    def _get_operands(self, ir_expr):
        '''Given a Numba IR expression, return the operands to the expression
        in order they appear in the expression.
        '''
        ir_op = ir_expr.op
        if ir_op == 'binop':
            return ir_expr.lhs, ir_expr.rhs
        elif ir_op == 'unary':
            return ir_expr.list_vars()
        elif ir_op == 'call':
            return ir_expr.args
        raise NotImplementedError(
            "Don't know how to find the operands for '{0}' expressions.".format(
                ir_op))

    def _translate_expr(self, ir_expr):
        '''Translate the given expression from Numba IR to an array expression
        tree.
        '''
        ir_op = ir_expr.op
        if ir_op == 'arrayexpr':
            return ir_expr.expr
        operands_or_args = [self.const_assigns.get(op_var.name, op_var)
                            for op_var in self._get_operands(ir_expr)]
        return self._get_array_operator(ir_expr), operands_or_args

    def _handle_matches(self):
        '''Iterate over the matches, trying to find which instructions should
        be rewritten, deleted, or moved.
        '''
        replace_map = {}
        dead_vars = set()
        used_vars = defaultdict(int)
        for instr in self.array_assigns.values():
            expr = instr.value
            arr_inps = []
            arr_expr = self._get_array_operator(expr), arr_inps
            new_expr = ir.Expr(op='arrayexpr',
                               loc=expr.loc,
                               expr=arr_expr,
                               ty=self.typemap[instr.target.name])
            new_instr = ir.Assign(new_expr, instr.target, instr.loc)
            replace_map[instr] = new_instr
            self.array_assigns[instr.target.name] = new_instr
            for operand in self._get_operands(expr):
                operand_name = operand.name
                if operand.is_temp and operand_name in self.array_assigns:
                    child_assign = self.array_assigns[operand_name]
                    child_expr = child_assign.value
                    child_operands = child_expr.list_vars()
                    for operand in child_operands:
                        used_vars[operand.name] += 1
                    arr_inps.append(self._translate_expr(child_expr))
                    if child_assign.target.is_temp:
                        dead_vars.add(child_assign.target.name)
                        replace_map[child_assign] = None
                elif operand_name in self.const_assigns:
                    arr_inps.append(self.const_assigns[operand_name])
                else:
                    used_vars[operand.name] += 1
                    arr_inps.append(operand)
        return replace_map, dead_vars, used_vars

    def _get_final_replacement(self, replacement_map, instr):
        '''Find the final replacement instruction for a given initial
        instruction by chasing instructions in a map from instructions
        to replacement instructions.
        '''
        replacement = replacement_map[instr]
        while replacement in replacement_map:
            replacement = replacement_map[replacement]
        return replacement

    def apply(self):
        '''When we've found array expressions in a basic block, rewrite that
        block, returning a new, transformed block.
        '''
        # Part 1: Figure out what instructions should be rewritten
        # based on the matches found.
        replace_map, dead_vars, used_vars = self._handle_matches()
        # Part 2: Using the information above, rewrite the target
        # basic block.
        result = self.crnt_block.copy()
        result.clear()
        delete_map = {}
        for instr in self.crnt_block.body:
            if isinstance(instr, ir.Assign):
                if instr in replace_map:
                    replacement = self._get_final_replacement(
                        replace_map, instr)
                    if replacement:
                        result.append(replacement)
                        for var in replacement.value.list_vars():
                            var_name = var.name
                            if var_name in delete_map:
                                result.append(delete_map.pop(var_name))
                            if used_vars[var_name] > 0:
                                used_vars[var_name] -= 1

                else:
                    result.append(instr)
            elif isinstance(instr, ir.Del):
                instr_value = instr.value
                if used_vars[instr_value] > 0:
                    used_vars[instr_value] -= 1
                    delete_map[instr_value] = instr
                elif instr_value not in dead_vars:
                    result.append(instr)
            else:
                result.append(instr)
        if delete_map:
            for instr in delete_map.values():
                result.insert_before_terminator(instr)
        return result


_unaryops = {
    operator.pos: ast.UAdd,
    operator.neg: ast.USub,
    operator.invert: ast.Invert,
}

_binops = {
    operator.add: ast.Add,
    operator.sub: ast.Sub,
    operator.mul: ast.Mult,
    operator.truediv: ast.Div,
    operator.mod: ast.Mod,
    operator.or_: ast.BitOr,
    operator.rshift: ast.RShift,
    operator.xor: ast.BitXor,
    operator.lshift: ast.LShift,
    operator.and_: ast.BitAnd,
    operator.pow: ast.Pow,
    operator.floordiv: ast.FloorDiv,
}


_cmpops = {
    operator.eq: ast.Eq,
    operator.ne: ast.NotEq,
    operator.lt: ast.Lt,
    operator.le: ast.LtE,
    operator.gt: ast.Gt,
    operator.ge: ast.GtE,
}


def _arr_expr_to_ast(expr):
    '''Build a Python expression AST from an array expression built by
    RewriteArrayExprs.
    '''
    if isinstance(expr, tuple):
        op, arr_expr_args = expr
        ast_args = []
        env = {}
        for arg in arr_expr_args:
            ast_arg, child_env = _arr_expr_to_ast(arg)
            ast_args.append(ast_arg)
            env.update(child_env)
        if op in npydecl.supported_array_operators:
            if len(ast_args) == 2:
                if op in _binops:
                    return ast.BinOp(
                        ast_args[0], _binops[op](), ast_args[1]), env
                if op in _cmpops:
                    return ast.Compare(
                        ast_args[0], [_cmpops[op]()], [ast_args[1]]), env
            else:
                assert op in _unaryops
                return ast.UnaryOp(_unaryops[op](), ast_args[0]), env
        elif _is_ufunc(op):
            fn_name = "__ufunc_or_dufunc_{0}".format(
                hex(hash(op)).replace("-", "_"))
            fn_ast_name = ast.Name(fn_name, ast.Load())
            env[fn_name] = op # Stash the ufunc or DUFunc in the environment
            ast_call = ast.Call(fn_ast_name, ast_args, [])
            return ast_call, env
    elif isinstance(expr, ir.Var):
        return ast.Name(expr.name, ast.Load(),
                        lineno=expr.loc.line,
                        col_offset=expr.loc.col if expr.loc.col else 0), {}
    elif isinstance(expr, ir.Const):
        return ast.Constant(expr.value), {}
    raise NotImplementedError(
        "Don't know how to translate array expression '%r'" % (expr,))


@contextlib.contextmanager
def _legalize_parameter_names(var_list):
    """
    Legalize names in the variable list for use as a Python function's
    parameter names.
    """
    var_map = OrderedDict()
    for var in var_list:
        old_name = var.name
        new_name = var.scope.redefine(old_name, loc=var.loc).name
        new_name = new_name.replace("$", "_").replace(".", "_")
        # Caller should ensure the names are unique
        if new_name in var_map:
            raise AssertionError(f"{new_name!r} not unique")
        var_map[new_name] = var, old_name
        var.name = new_name
    param_names = list(var_map)
    try:
        yield param_names
    finally:
        # Make sure the old names are restored, to avoid confusing
        # other parts of Numba (see issue #1466)
        for var, old_name in var_map.values():
            var.name = old_name


class _EraseInvalidLineRanges(ast.NodeTransformer):
    def generic_visit(self, node: ast.AST) -> ast.AST:
        node = super().generic_visit(node)
        if hasattr(node, "lineno"):
            if getattr(node, "end_lineno", None) is not None:
                if node.lineno > node.end_lineno:
                    del node.lineno
                    del node.end_lineno
        return node


def _fix_invalid_lineno_ranges(astree: ast.AST):
    """Inplace fixes invalid lineno ranges.
    """
    # Make sure lineno and end_lineno are present
    ast.fix_missing_locations(astree)
    # Delete invalid lineno ranges
    _EraseInvalidLineRanges().visit(astree)
    # Make sure lineno and end_lineno are present
    ast.fix_missing_locations(astree)


def _lower_array_expr(lowerer, expr):
    '''Lower an array expression built by RewriteArrayExprs.
    '''
    expr_name = "__numba_array_expr_%s" % (hex(hash(expr)).replace("-", "_"))
    expr_filename = expr.loc.filename
    expr_var_list = expr.list_vars()
    # The expression may use a given variable several times, but we
    # should only create one parameter for it.
    expr_var_unique = sorted(set(expr_var_list), key=lambda var: var.name)

    # Arguments are the names external to the new closure
    expr_args = [var.name for var in expr_var_unique]

    # 1. Create an AST tree from the array expression.
    with _legalize_parameter_names(expr_var_unique) as expr_params:
        ast_args = [ast.arg(param_name, None)
                    for param_name in expr_params]
        # Parse a stub function to ensure the AST is populated with
        # reasonable defaults for the Python version.
        ast_module = ast.parse('def {0}(): return'.format(expr_name),
                               expr_filename, 'exec')
        assert hasattr(ast_module, 'body') and len(ast_module.body) == 1
        ast_fn = ast_module.body[0]
        ast_fn.args.args = ast_args
        ast_fn.body[0].value, namespace = _arr_expr_to_ast(expr.expr)
        _fix_invalid_lineno_ranges(ast_module)

    # 2. Compile the AST module and extract the Python function.
    code_obj = compile(ast_module, expr_filename, 'exec')
    exec(code_obj, namespace)
    impl = namespace[expr_name]

    # 3. Now compile a ufunc using the Python function as kernel.

    context = lowerer.context
    builder = lowerer.builder
    outer_sig = expr.ty(*(lowerer.typeof(name) for name in expr_args))
    inner_sig_args = []
    for argty in outer_sig.args:
        if isinstance(argty, types.Optional):
            argty = argty.type
        if isinstance(argty, types.Array):
            inner_sig_args.append(argty.dtype)
        else:
            inner_sig_args.append(argty)
    inner_sig = outer_sig.return_type.dtype(*inner_sig_args)

    flags = targetconfig.ConfigStack().top_or_none()
    flags = compiler.Flags() if flags is None else flags.copy() # make sure it's a clone or a fresh instance
    # Follow the Numpy error model.  Note this also allows e.g. vectorizing
    # division (issue #1223).
    flags.error_model = 'numpy'
    cres = context.compile_subroutine(builder, impl, inner_sig, flags=flags,
                                      caching=False)

    # Create kernel subclass calling our native function
    from numba.np import npyimpl

    class ExprKernel(npyimpl._Kernel):
        def generate(self, *args):
            arg_zip = zip(args, self.outer_sig.args, inner_sig.args)
            cast_args = [self.cast(val, inty, outty)
                         for val, inty, outty in arg_zip]
            result = self.context.call_internal(
                builder, cres.fndesc, inner_sig, cast_args)
            return self.cast(result, inner_sig.return_type,
                             self.outer_sig.return_type)

    # create a fake ufunc object which is enough to trick numpy_ufunc_kernel
    ufunc = SimpleNamespace(nin=len(expr_args), nout=1, __name__=expr_name)
    ufunc.nargs = ufunc.nin + ufunc.nout

    args = [lowerer.loadvar(name) for name in expr_args]
    return npyimpl.numpy_ufunc_kernel(
        context, builder, outer_sig, args, ufunc, ExprKernel)
