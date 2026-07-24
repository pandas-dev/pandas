"""Python parser that directly constructs a native AST (when compiled).

Use a Rust extension to generate a serialized AST, and deserialize the AST directly
to a mypy AST.

NOTE: This is work in progress. To use this, you need to manually build the
      ast_serialize Rust extension. See the README at https://github.com/mypyc/ast_serialize.

Expected benefits over mypy.fastparse:
 * No intermediate non-mypyc Python-level AST created, to improve performance
 * Parsing doesn't need GIL => use multithreading to construct serialized ASTs in parallel
 * Produce import dependencies without having to build an AST => helps parallel type checking
 * Support all Python syntax even if mypy is running on an older Python version
 * Generate an AST even if there are syntax errors
 * Potential to support incremental parsing (quickly process modified sections in a file)
 * Stripping function bodies in third-party code can happen earlier, for extra performance
"""

from __future__ import annotations

import os
import time
from typing import Final, cast

import ast_serialize
from librt.internal import (
    read_float as read_float_bare,
    read_int as read_int_bare,
    read_str as read_str_bare,
)

from mypy import message_registry, nodes, types
from mypy.cache import (
    DICT_STR_GEN,
    END_TAG,
    LIST_GEN,
    LIST_INT,
    LITERAL_FLOAT,
    LITERAL_NONE,
    LITERAL_STR,
    LOCATION,
    ReadBuffer,
    Tag,
    read_bool,
    read_int,
    read_str,
    read_str_opt,
    read_tag,
)
from mypy.nodes import (
    ARG_KINDS,
    ARG_POS,
    IMPORT_METADATA,
    IMPORTALL_METADATA,
    IMPORTFROM_METADATA,
    MISSING_FALLBACK,
    Argument,
    AssertStmt,
    AssignmentExpr,
    AssignmentStmt,
    AwaitExpr,
    Block,
    BreakStmt,
    BytesExpr,
    CallExpr,
    ClassDef,
    ComparisonExpr,
    ComplexExpr,
    ConditionalExpr,
    Context,
    ContinueStmt,
    Decorator,
    DelStmt,
    DictExpr,
    DictionaryComprehension,
    EllipsisExpr,
    Expression,
    ExpressionStmt,
    FileRawData,
    FloatExpr,
    ForStmt,
    FuncDef,
    GeneratorExpr,
    GlobalDecl,
    IfStmt,
    Import,
    ImportAll,
    ImportBase,
    ImportFrom,
    IndexExpr,
    IntExpr,
    LambdaExpr,
    ListComprehension,
    ListExpr,
    MatchStmt,
    MemberExpr,
    MypyFile,
    NameExpr,
    NonlocalDecl,
    OperatorAssignmentStmt,
    OpExpr,
    OverloadedFuncDef,
    OverloadPart,
    ParseError,
    PassStmt,
    RaiseStmt,
    RefExpr,
    ReturnStmt,
    SetComprehension,
    SetExpr,
    SliceExpr,
    StarExpr,
    Statement,
    StrExpr,
    SuperExpr,
    TemplateStrExpr,
    TempNode,
    TryStmt,
    TupleExpr,
    TypeAliasStmt,
    TypeParam,
    UnaryExpr,
    Var,
    WhileStmt,
    WithStmt,
    YieldExpr,
    YieldFromExpr,
)
from mypy.options import Options
from mypy.patterns import (
    AsPattern,
    ClassPattern,
    MappingPattern,
    OrPattern,
    Pattern,
    SequencePattern,
    SingletonPattern,
    StarredPattern,
    ValuePattern,
)
from mypy.reachability import infer_reachability_of_if_statement
from mypy.sharedparse import special_function_elide_names
from mypy.types import (
    AnyType,
    CallableArgument,
    CallableType,
    EllipsisType,
    Instance,
    ProperType,
    RawExpressionType,
    TupleType,
    Type,
    TypedDictType,
    TypeList,
    TypeOfAny,
    UnboundType,
    UnionType,
    UnpackType,
)
from mypy.util import unnamed_function

TypeIgnores = list[tuple[int, list[str]]]

# There is no way to create reasonable fallbacks at this stage,
# they must be patched later.
_dummy_fallback: Final = Instance(MISSING_FALLBACK, [], -1)


class State:
    def __init__(self, options: Options) -> None:
        self.options = options
        self.errors: list[ParseError] = []
        self.num_funcs = 0

    def add_error(
        self, message: str, line: int, column: int, *, blocker: bool = False, code: str
    ) -> None:
        """Report an error at a specific location."""
        self.errors.append(
            {"line": line, "column": column, "message": message, "blocker": blocker, "code": code}
        )


def native_parse(
    filename: str, options: Options, skip_function_bodies: bool = False
) -> tuple[MypyFile, list[ParseError], TypeIgnores]:
    """Parse a Python file using the native Rust-based parser.

    Return (MypyFile, errors, type_ignores).

    The returned tree is empty with actual serialized data stored in `raw_data`
    attribute. Use read_statements() and/or deserialize_imports() to de-serialize.

    The caller should set these additional attributes on the returned MypyFile:
      - ignored_lines: dict of type ignore comments (from the TypeIgnores return value)
      - is_stub: whether the file is a .pyi stub
    """
    # If the path is a directory, return empty AST (matching fastparse behavior)
    # This can happen for packages that only contain .pyc files without source
    if os.path.isdir(filename):
        node = MypyFile([], [])
        node.path = filename
        return node, [], []

    (
        b,
        errors,
        ignores,
        import_bytes,
        is_partial_package,
        uses_template_strings,
        source_hash,
        mypy_comments,
    ) = parse_to_binary_ast(filename, options, skip_function_bodies)
    node = MypyFile([], [])
    node.path = filename
    node.raw_data = FileRawData(
        b,
        import_bytes,
        errors,
        dict(ignores),
        is_partial_package,
        uses_template_strings,
        source_hash,
        mypy_comments,
    )
    return node, errors, ignores


def expect_end_tag(data: ReadBuffer) -> None:
    assert read_tag(data) == END_TAG


def expect_tag(data: ReadBuffer, tag: Tag) -> None:
    assert (actual := read_tag(data)) == tag, actual


def read_statements(state: State, data: ReadBuffer, n: int) -> list[Statement]:
    defs: list[Statement] = []
    old_num_funcs = state.num_funcs
    for _ in range(n):
        stmt = read_statement(state, data)
        defs.append(stmt)
    if state.num_funcs > old_num_funcs + 1:
        # There were at least two functions, so we may need to merge overloads.
        defs = fix_function_overloads(state, defs)
    return defs


def parse_to_binary_ast(
    filename: str, options: Options, skip_function_bodies: bool = False
) -> tuple[bytes, list[ParseError], TypeIgnores, bytes, bool, bool, str, list[tuple[int, str]]]:
    # This is a horrible hack to work around a mypyc bug where imported
    # module may be not ready in a thread sometimes.
    t0 = time.time()
    while ast_serialize is None:
        time.sleep(0.0001)  # type: ignore[unreachable]
        if time.time() - t0 > 10.0:
            raise ImportError("Cannot import ast_serialize")
    ast_bytes, errors, ignores, import_bytes, ast_data = ast_serialize.parse(
        filename,
        skip_function_bodies=skip_function_bodies,
        python_version=options.python_version,
        platform=options.platform,
        always_true=options.always_true,
        always_false=options.always_false,
        cache_version=3,
    )
    return (
        ast_bytes,
        errors,
        ignores,
        import_bytes,
        ast_data["is_partial_package"],
        ast_data["uses_template_strings"],
        ast_data["source_hash"],
        ast_data["mypy_comments"],
    )


def read_statement(state: State, data: ReadBuffer) -> Statement:
    # Branches ordered by frequency (based on mypy self-check)
    tag = read_tag(data)
    stmt: Statement
    if tag == nodes.ASSIGNMENT_STMT:
        lvalues = read_expression_list(state, data)
        rvalue = read_expression(state, data)
        has_type = read_bool(data)
        if has_type:
            type_annotation = read_type(state, data)
        else:
            type_annotation = None
        new_syntax = read_bool(data)
        a = AssignmentStmt(lvalues, rvalue, type=type_annotation, new_syntax=new_syntax)
        read_loc(data, a)
        # If rvalue is TempNode, copy location from AssignmentStmt
        if isinstance(rvalue, TempNode):
            set_line_column_range(rvalue, a)
        expect_end_tag(data)
        return a
    elif tag == nodes.EXPR_STMT:
        es = ExpressionStmt(read_expression(state, data))
        set_line_column_range(es, es.expr)
        expect_end_tag(data)
        return es
    elif tag == nodes.IF_STMT:
        expr = read_expression(state, data)
        body = read_block(state, data)

        num_elif = read_int(data)
        elif_exprs = []
        elif_bodies = []
        for i in range(num_elif):
            elif_exprs.append(read_expression(state, data))
            elif_bodies.append(read_block(state, data))

        has_else = read_bool(data)
        if has_else:
            else_body = read_block(state, data)
        else:
            else_body = None

        # Normalize elif into nested if/else statements
        # Build from the bottom up, starting with the final else body
        current_else = else_body

        for elif_expr, elif_body in reversed(list(zip(elif_exprs, elif_bodies))):
            elif_stmt = IfStmt([elif_expr], [elif_body], current_else)
            elif_stmt.line = elif_expr.line
            elif_stmt.column = elif_expr.column
            if current_else is not None:
                elif_stmt.end_line = current_else.end_line
                elif_stmt.end_column = current_else.end_column
            else:
                elif_stmt.end_line = elif_body.end_line
                elif_stmt.end_column = elif_body.end_column

            current_else = Block([elif_stmt])
            set_line_column_range(current_else, elif_stmt)

        if_stmt = IfStmt([expr], [body], current_else)
        read_loc(data, if_stmt)
        expect_end_tag(data)
        return if_stmt
    elif tag == nodes.RETURN_STMT:
        has_value = read_bool(data)
        if has_value:
            value = read_expression(state, data)
        else:
            value = None
        stmt = ReturnStmt(value)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.FUNC_DEF_STMT:
        return read_func_def(state, data)
    elif tag == nodes.IMPORT_FROM:
        relative = read_int(data)
        module_id = read_str(data)  # Empty string for "from . import x"
        n = read_int(data)
        names = []
        for _ in range(n):
            name = read_str(data)
            has_asname = read_bool(data)
            if has_asname:
                asname = read_str(data)
            else:
                asname = None
            names.append((name, asname))

        stmt = ImportFrom(module_id, relative, names)
        _read_and_set_import_metadata(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.FOR_STMT:
        index = read_expression(state, data)
        expr = read_expression(state, data)
        body = read_block(state, data)
        else_body = read_optional_block(state, data)
        is_async = read_bool(data)
        stmt = ForStmt(index, expr, body, else_body)
        stmt.is_async = is_async
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.ASSERT_STMT:
        test = read_expression(state, data)
        has_msg = read_bool(data)
        if has_msg:
            msg = read_expression(state, data)
        else:
            msg = None
        stmt = AssertStmt(test, msg)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.CLASS_DEF:
        return read_class_def(state, data)
    elif tag == nodes.DECORATOR:
        expect_tag(data, LIST_GEN)
        n_decorators = read_int_bare(data)
        decorators = [read_expression(state, data) for i in range(n_decorators)]
        line = read_int(data)
        column = read_int(data)
        fdef = read_statement(state, data)
        assert isinstance(fdef, FuncDef)
        fdef.is_decorated = True
        var = Var(fdef.name)
        var.line = fdef.line
        var.is_ready = False
        stmt = Decorator(fdef, decorators, var)
        stmt.line = line
        stmt.column = column
        stmt.end_line = fdef.end_line
        stmt.end_column = fdef.end_column
        # TODO: Adjust funcdef location to start after decorator?
        expect_end_tag(data)
        return stmt
    elif tag == nodes.IMPORT:
        n = read_int(data)
        ids = []
        for _ in range(n):
            name = read_str(data)
            has_asname = read_bool(data)
            if has_asname:
                asname = read_str(data)
            else:
                asname = None
            ids.append((name, asname))
        stmt = Import(ids)
        _read_and_set_import_metadata(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.RAISE_STMT:
        has_exc = read_bool(data)
        if has_exc:
            exc = read_expression(state, data)
        else:
            exc = None
        has_from = read_bool(data)
        if has_from:
            from_expr = read_expression(state, data)
        else:
            from_expr = None
        stmt = RaiseStmt(exc, from_expr)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.OPERATOR_ASSIGNMENT_STMT:
        op = read_str(data)
        lvalue = read_expression(state, data)
        rvalue = read_expression(state, data)
        stmt = OperatorAssignmentStmt(op, lvalue, rvalue)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.PASS_STMT:
        stmt = PassStmt()
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.CONTINUE_STMT:
        stmt = ContinueStmt()
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.WITH_STMT:
        n = read_int(data)
        expr_list = []
        target_list: list[Expression | None] = []
        for _ in range(n):
            context_expr = read_expression(state, data)
            expr_list.append(context_expr)
            has_target = read_bool(data)
            if has_target:
                target = read_expression(state, data)
                target_list.append(target)
            else:
                target_list.append(None)
        body = read_block(state, data)
        is_async = read_bool(data)
        stmt = WithStmt(expr_list, target_list, body)
        stmt.is_async = is_async
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.TRY_STMT:
        return read_try_stmt(state, data)
    elif tag == nodes.BREAK_STMT:
        stmt = BreakStmt()
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.WHILE_STMT:
        expr = read_expression(state, data)
        body = read_block(state, data)
        else_body = read_optional_block(state, data)
        stmt = WhileStmt(expr, body, else_body)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.DEL_STMT:
        expr = read_expression(state, data)
        stmt = DelStmt(expr)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.TYPE_ALIAS_STMT:
        return read_type_alias_stmt(state, data)
    elif tag == nodes.IMPORT_ALL:
        module_id = read_str(data)  # Empty string for "from . import *"
        relative = read_int(data)

        stmt = ImportAll(module_id, relative)
        _read_and_set_import_metadata(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.NONLOCAL_DECL:
        n = read_int(data)
        decl_names = []
        for _ in range(n):
            decl_names.append(read_str(data))
        stmt = NonlocalDecl(decl_names)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.GLOBAL_DECL:
        n = read_int(data)
        decl_names = []
        for _ in range(n):
            decl_names.append(read_str(data))
        stmt = GlobalDecl(decl_names)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    elif tag == nodes.MATCH_STMT:
        subject = read_expression(state, data)
        n_cases = read_int(data)
        patterns = []
        guards: list[Expression | None] = []
        bodies = []
        for _ in range(n_cases):
            pattern = read_pattern(state, data)
            patterns.append(pattern)
            has_guard = read_bool(data)
            if has_guard:
                guard = read_expression(state, data)
                guards.append(guard)
            else:
                guards.append(None)
            body = read_block(state, data)
            bodies.append(body)
        stmt = MatchStmt(subject, patterns, guards, bodies)
        read_loc(data, stmt)
        expect_end_tag(data)
        return stmt
    else:
        assert False, tag


def read_parameters(state: State, data: ReadBuffer) -> tuple[list[Argument], bool]:
    """Read function/lambda parameters.

    Return (parameters, has_annotations).
    """
    expect_tag(data, LIST_GEN)
    n_args = read_int_bare(data)
    arguments = []
    has_ann = False
    for _ in range(n_args):
        arg_name = read_str(data)
        arg_kind_int = read_int(data)
        arg_kind = ARG_KINDS[arg_kind_int]
        has_type = read_bool(data)
        if has_type:
            ann = read_type(state, data)
            has_ann = True
        else:
            ann = None
        has_default = read_bool(data)
        if has_default:
            default = read_expression(state, data)
        else:
            default = None
        pos_only = read_bool(data)

        if state.options.implicit_optional and ann is not None:
            optional = isinstance(default, NameExpr) and default.name == "None"
            if isinstance(ann, UnboundType):
                ann.optional = optional

        var = Var(arg_name, ann)
        var.is_inferred = False
        var.is_argument = True
        arg = Argument(var, ann, default, arg_kind, pos_only)
        read_loc(data, arg)
        set_line_column_range(var, arg)
        arguments.append(arg)

    return arguments, has_ann


def read_type_params(state: State, data: ReadBuffer) -> list[TypeParam]:
    """Read type parameters (PEP 695 generics)."""
    type_params: list[TypeParam] = []
    n = read_int_bare(data)
    for _ in range(n):
        kind = read_int(data)
        name = read_str(data)
        has_bound = read_bool(data)
        if has_bound:
            upper_bound = read_type(state, data)
        else:
            upper_bound = None

        expect_tag(data, LIST_GEN)
        n_values = read_int_bare(data)
        values = [read_type(state, data) for _ in range(n_values)]

        has_default = read_bool(data)
        if has_default:
            default = read_type(state, data)
        else:
            default = None

        type_params.append(TypeParam(name, kind, upper_bound, values, default))

    return type_params


def read_func_def(state: State, data: ReadBuffer) -> FuncDef:
    state.num_funcs += 1

    name = read_str(data)
    arguments, has_ann = read_parameters(state, data)

    if special_function_elide_names(name):
        for arg in arguments:
            arg.pos_only = True

    body = read_block(state, data)
    is_async = read_bool(data)

    # Type parameters (PEP 695)
    has_type_params = read_bool(data)
    if has_type_params:
        type_params = read_type_params(state, data)
    else:
        type_params = None

    has_return_type = read_bool(data)
    if has_return_type:
        return_type = read_type(state, data)
        has_ann = True
    else:
        return_type = None

    if has_ann:
        typ = CallableType(
            [
                arg.type_annotation if arg.type_annotation else AnyType(TypeOfAny.unannotated)
                for arg in arguments
            ],
            [arg.kind for arg in arguments],
            [None if arg.pos_only else arg.variable.name for arg in arguments],
            return_type if return_type else AnyType(TypeOfAny.unannotated),
            _dummy_fallback,
        )
    else:
        typ = None

    func_def = FuncDef(name, arguments, body, typ=typ, type_args=type_params)
    if is_async:
        func_def.is_coroutine = True
    read_loc(data, func_def)
    if typ:
        typ.line = func_def.line
        typ.column = func_def.column
        typ.definition = func_def
        # TODO: This seems wasteful, can we avoid it?
        func_def.unanalyzed_type = typ.copy_modified()
    expect_end_tag(data)
    return func_def


def read_class_def(state: State, data: ReadBuffer) -> ClassDef:
    name = read_str(data)
    body = read_block(state, data)
    base_type_exprs = read_expression_list(state, data)

    expect_tag(data, LIST_GEN)
    n_decorators = read_int_bare(data)
    decorators = [read_expression(state, data) for _ in range(n_decorators)]

    # Type parameters (PEP 695)
    has_type_params = read_bool(data)
    if has_type_params:
        type_params = read_type_params(state, data)
    else:
        type_params = None

    expect_tag(data, DICT_STR_GEN)
    n_keywords = read_int_bare(data)
    keywords = []
    for _ in range(n_keywords):
        key = read_str(data)
        value = read_expression(state, data)
        keywords.append((key, value))

    metaclass = dict(keywords).get("metaclass") if keywords else None

    class_def = ClassDef(
        name,
        body,
        base_type_exprs=base_type_exprs if base_type_exprs else None,
        metaclass=metaclass,
        # Note we keep metaclass in keywords as well, to match the old parser.
        keywords=keywords if keywords else None,
        type_args=type_params,
    )
    class_def.decorators = decorators
    read_loc(data, class_def)
    expect_end_tag(data)
    return class_def


def read_type_alias_stmt(state: State, data: ReadBuffer) -> TypeAliasStmt:
    """Read PEP 695 type alias statement."""
    name = read_expression(state, data)
    assert isinstance(name, NameExpr), f"Expected NameExpr for type alias name, got {type(name)}"

    n_type_params = read_int_bare(data)
    if n_type_params > 0:
        type_params = []
        for _ in range(n_type_params):
            kind = read_int(data)
            param_name = read_str(data)
            has_bound = read_bool(data)
            if has_bound:
                upper_bound = read_type(state, data)
            else:
                upper_bound = None

            # Read values (for constrained TypeVar)
            expect_tag(data, LIST_GEN)
            n_values = read_int_bare(data)
            values = [read_type(state, data) for _ in range(n_values)]

            has_default = read_bool(data)
            if has_default:
                default = read_type(state, data)
            else:
                default = None

            type_params.append(TypeParam(param_name, kind, upper_bound, values, default))
    else:
        type_params = []

    value_expr = read_expression(state, data)

    # Wrap the value expression in a LambdaExpr as expected by TypeAliasStmt
    # The LambdaExpr body is a Block with a single ReturnStmt
    return_stmt = ReturnStmt(value_expr)
    set_line_column_range(return_stmt, value_expr)

    block = Block([return_stmt])
    block.line = -1  # Synthetic block
    block.column = 0
    block.end_line = -1
    block.end_column = 0

    lambda_expr = LambdaExpr([], block)
    set_line_column_range(lambda_expr, value_expr)

    stmt = TypeAliasStmt(name, type_params, lambda_expr)
    read_loc(data, stmt)
    expect_end_tag(data)
    return stmt


def read_try_stmt(state: State, data: ReadBuffer) -> TryStmt:
    body = read_block(state, data)
    num_handlers = read_int(data)

    types_list: list[Expression | None] = []
    for _ in range(num_handlers):
        has_type = read_bool(data)
        if has_type:
            exc_type = read_expression(state, data)
            types_list.append(exc_type)
        else:
            types_list.append(None)

    vars_list: list[NameExpr | None] = []
    for _ in range(num_handlers):
        has_name = read_bool(data)
        if has_name:
            var_name = read_str(data)
            var_expr = NameExpr(var_name)
            read_loc(data, var_expr)
            vars_list.append(var_expr)
        else:
            vars_list.append(None)

    handlers = []
    for _ in range(num_handlers):
        handler_body = read_block(state, data)
        handlers.append(handler_body)

    has_else = read_bool(data)
    if has_else:
        else_body = read_block(state, data)
    else:
        else_body = None

    has_finally = read_bool(data)
    if has_finally:
        finally_body = read_block(state, data)
    else:
        finally_body = None

    # except* (Python 3.11+)
    is_star = read_bool(data)

    stmt = TryStmt(body, vars_list, types_list, handlers, else_body, finally_body)
    stmt.is_star = is_star
    read_loc(data, stmt)
    expect_end_tag(data)
    return stmt


def read_type(state: State, data: ReadBuffer) -> Type:
    tag = read_tag(data)
    if tag == types.UNBOUND_TYPE:
        name = read_str(data)
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        args = tuple(read_type(state, data) for i in range(n))
        empty_tuple_index = read_bool(data)
        t = read_tag(data)
        if t == LITERAL_NONE:
            original_str_expr = None
        elif t == LITERAL_STR:
            original_str_expr = read_str_bare(data)
        else:
            assert False, f"Unexpected tag for original_str_expr: {t}"
        t = read_tag(data)
        if t == LITERAL_NONE:
            original_str_fallback = None
        elif t == LITERAL_STR:
            original_str_fallback = read_str_bare(data)
        else:
            assert False, f"Unexpected tag for original_str_fallback: {t}"
        unbound = UnboundType(
            name,
            args,
            empty_tuple_index=empty_tuple_index,
            original_str_expr=original_str_expr,
            original_str_fallback=original_str_fallback,
        )
        read_loc(data, unbound)
        expect_end_tag(data)
        return unbound
    elif tag == types.UNION_TYPE:
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        items = [read_type(state, data) for i in range(n)]
        uses_pep604_syntax = read_bool(data)
        t = read_tag(data)
        if t == LITERAL_NONE:
            original_str_expr = None
        elif t == LITERAL_STR:
            original_str_expr = read_str_bare(data)
        else:
            assert False, f"Unexpected tag for original_str_expr: {t}"
        t = read_tag(data)
        if t == LITERAL_NONE:
            original_str_fallback = None
        elif t == LITERAL_STR:
            original_str_fallback = read_str_bare(data)
        else:
            assert False, f"Unexpected tag for original_str_fallback: {t}"
        union = UnionType(items, uses_pep604_syntax=uses_pep604_syntax)
        union.original_str_expr = original_str_expr
        union.original_str_fallback = original_str_fallback
        union.is_evaluated = read_bool(data)
        read_loc(data, union)
        expect_end_tag(data)
        return union
    elif tag == types.LIST_TYPE:
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        items = [read_type(state, data) for i in range(n)]
        type_list = TypeList(items)
        read_loc(data, type_list)
        expect_end_tag(data)
        return type_list
    elif tag == types.TUPLE_TYPE:
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        items = [read_type(state, data) for i in range(n)]
        implicit = read_bool(data)
        tuple_type = TupleType(items, _dummy_fallback, implicit=implicit)
        read_loc(data, tuple_type)
        expect_end_tag(data)
        return tuple_type
    elif tag == types.TYPED_DICT_TYPE:
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        keys = [read_str_opt(data) for i in range(n)]
        expect_tag(data, LIST_GEN)
        n = read_int_bare(data)
        values = [read_type(state, data) for i in range(n)]
        td_items = {}
        extra_items_from = []
        for key, val in zip(keys, values):
            if key is None:
                assert isinstance(val, ProperType)
                extra_items_from.append(val)
            else:
                td_items[key] = val
        typeddict_type = TypedDictType(td_items, set(), set(), _dummy_fallback)
        typeddict_type.extra_items_from = extra_items_from
        read_loc(data, typeddict_type)
        expect_end_tag(data)
        return typeddict_type
    elif tag == types.ELLIPSIS_TYPE:
        ellipsis_type = EllipsisType()
        read_loc(data, ellipsis_type)
        expect_end_tag(data)
        return ellipsis_type
    elif tag == types.RAW_EXPRESSION_TYPE:
        type_name = read_str(data)
        value: types.LiteralValue | str | None
        note: str | None = None
        if type_name == "builtins.bool":
            value = read_bool(data)
        elif type_name == "builtins.int":
            value = read_int(data)
        elif type_name == "builtins.str":
            value = read_str(data)
        elif type_name == "builtins.bytes":
            # Bytes literals are serialized as escaped strings
            value = read_str(data)
        elif type_name == "typing.Any":
            # Invalid type - read None value
            tag = read_tag(data)
            assert tag == LITERAL_NONE, f"Expected LITERAL_NONE for invalid type, got {tag}"
            value = None
            # Read optional note (cache_version >= 2)
            note = read_str_opt(data)
        else:
            assert False, f"Unsupported RawExpressionType: {type_name}"
        raw_type = RawExpressionType(value, type_name, note=note)
        read_loc(data, raw_type)
        expect_end_tag(data)
        return raw_type
    elif tag == types.UNPACK_TYPE:
        inner_type = read_type(state, data)
        from_star_syntax = read_bool(data)
        unpack = UnpackType(inner_type, from_star_syntax=from_star_syntax)
        read_loc(data, unpack)
        expect_end_tag(data)
        return unpack
    elif tag == types.CALL_TYPE:
        return read_call_type(state, data)
    else:
        assert False, tag


def stringify_type_name(typ: Type) -> str | None:
    if isinstance(typ, UnboundType):
        return typ.name
    return None


def extract_arg_name(typ: Type) -> str | None:
    if isinstance(typ, RawExpressionType) and typ.base_type_name == "builtins.str":
        return typ.literal_value  # type: ignore[return-value]
    elif isinstance(typ, UnboundType):
        if typ.name == "None":
            return None
        return typ.name
    return None  # Invalid, but let validation handle it


def read_call_type(state: State, data: ReadBuffer) -> Type:
    """Read Call in type context (Arg/DefaultArg/VarArg/KwArg constructor)."""
    callee_type = read_type(state, data)

    # Read positional arguments
    expect_tag(data, LIST_GEN)
    n_args = read_int_bare(data)
    args = [read_type(state, data) for _ in range(n_args)]

    # Read keyword arguments
    expect_tag(data, LIST_GEN)
    n_kwargs = read_int_bare(data)
    kwargs = []
    for _ in range(n_kwargs):
        tag_kw = read_tag(data)
        if tag_kw == LITERAL_NONE:
            kw_name = None
        elif tag_kw == LITERAL_STR:
            kw_name = read_str_bare(data)
        else:
            assert False, f"Unexpected tag for keyword name: {tag_kw}"
        kw_value = read_type(state, data)
        kwargs.append((kw_name, kw_value))

    constructor = stringify_type_name(callee_type)

    invalid = AnyType(TypeOfAny.from_error)
    read_loc(data, invalid)
    expect_end_tag(data)

    if not constructor:
        state.add_error(
            message_registry.ARG_CONSTRUCTOR_NAME_EXPECTED.value,
            invalid.line,
            invalid.column,
            blocker=True,
            code="misc",
        )
        return invalid

    # Extract type and name from arguments
    name: str | None = None
    name_set_from_positional = False
    default_type = AnyType(TypeOfAny.special_form)
    typ: Type = default_type
    typ_set_from_positional = False

    # Process positional arguments
    for i, arg in enumerate(args):
        if i == 0:
            typ = arg
            typ_set_from_positional = True
        elif i == 1:
            name = extract_arg_name(arg)
            name_set_from_positional = True
        else:
            state.add_error(
                message_registry.ARG_CONSTRUCTOR_TOO_MANY_ARGS.value,
                invalid.line,
                invalid.column,
                blocker=True,
                code="misc",
            )

    # Process keyword arguments
    for kw_name, kw_value in kwargs:
        if kw_name == "name":
            if name is not None and name_set_from_positional:
                state.add_error(
                    message_registry.MULTIPLE_VALUES_FOR_NAME_KWARG.format(constructor).value,
                    invalid.line,
                    invalid.column,
                    blocker=True,
                    code="misc",
                )
            name = extract_arg_name(kw_value)
        elif kw_name == "type":
            if typ is not default_type and typ_set_from_positional:
                state.add_error(
                    message_registry.MULTIPLE_VALUES_FOR_TYPE_KWARG.format(constructor).value,
                    invalid.line,
                    invalid.column,
                    blocker=True,
                    code="misc",
                )
            typ = kw_value
        else:
            state.add_error(
                message_registry.ARG_CONSTRUCTOR_UNEXPECTED_ARG.format(kw_name).value,
                invalid.line,
                invalid.column,
                blocker=True,
                code="misc",
            )

    call_arg = CallableArgument(typ, name, constructor)
    set_line_column_range(call_arg, invalid)
    return call_arg


def read_pattern(state: State, data: ReadBuffer) -> Pattern:
    tag = read_tag(data)
    if tag == nodes.AS_PATTERN:
        has_pattern = read_bool(data)
        if has_pattern:
            pattern = read_pattern(state, data)
        else:
            pattern = None
        has_name = read_bool(data)
        if has_name:
            name_str = read_str(data)
            name = NameExpr(name_str)
            read_loc(data, name)
        else:
            name = None
        as_pattern = AsPattern(pattern, name)
        read_loc(data, as_pattern)
        expect_end_tag(data)
        return as_pattern
    elif tag == nodes.OR_PATTERN:
        n = read_int(data)
        patterns = [read_pattern(state, data) for _ in range(n)]
        or_pattern = OrPattern(patterns)
        read_loc(data, or_pattern)
        expect_end_tag(data)
        return or_pattern
    elif tag == nodes.VALUE_PATTERN:
        expr = read_expression(state, data)
        value_pattern = ValuePattern(expr)
        read_loc(data, value_pattern)
        expect_end_tag(data)
        return value_pattern
    elif tag == nodes.SINGLETON_PATTERN:
        singleton_tag = read_tag(data)
        if singleton_tag == LITERAL_NONE:
            value = None
        else:
            # It's a boolean
            value = singleton_tag == 1  # TAG_LITERAL_TRUE
        singleton_pattern = SingletonPattern(value)
        read_loc(data, singleton_pattern)
        expect_end_tag(data)
        return singleton_pattern
    elif tag == nodes.SEQUENCE_PATTERN:
        n = read_int(data)
        patterns = [read_pattern(state, data) for _ in range(n)]
        sequence_pattern = SequencePattern(patterns)
        read_loc(data, sequence_pattern)
        expect_end_tag(data)
        return sequence_pattern
    elif tag == nodes.STARRED_PATTERN:
        has_name = read_bool(data)
        if has_name:
            name_str = read_str(data)
            name = NameExpr(name_str)
            read_loc(data, name)
        else:
            name = None
        starred_pattern = StarredPattern(name)
        read_loc(data, starred_pattern)
        expect_end_tag(data)
        return starred_pattern
    elif tag == nodes.MAPPING_PATTERN:
        n = read_int(data)
        keys = []
        values = []
        for _ in range(n):
            key = read_expression(state, data)
            value = read_pattern(state, data)
            keys.append(key)
            values.append(value)
        has_rest = read_bool(data)
        if has_rest:
            rest_str = read_str(data)
            rest = NameExpr(rest_str)
            read_loc(data, rest)
        else:
            rest = None
        mapping_pattern = MappingPattern(keys, values, rest)
        read_loc(data, mapping_pattern)
        expect_end_tag(data)
        return mapping_pattern
    elif tag == nodes.CLASS_PATTERN:
        class_ref = cast(RefExpr, read_expression(state, data))
        n_positional = read_int(data)
        positionals = [read_pattern(state, data) for _ in range(n_positional)]
        n_keywords = read_int(data)
        keyword_keys = []
        keyword_values = []
        for _ in range(n_keywords):
            key = read_str(data)
            value = read_pattern(state, data)
            keyword_keys.append(key)
            keyword_values.append(value)
        class_pattern = ClassPattern(class_ref, positionals, keyword_keys, keyword_values)
        read_loc(data, class_pattern)
        expect_end_tag(data)
        return class_pattern
    else:
        assert False, f"Unknown pattern tag: {tag}"


def read_block(state: State, data: ReadBuffer) -> Block:
    expect_tag(data, nodes.BLOCK)
    expect_tag(data, LIST_GEN)
    n = read_int_bare(data)
    is_unreachable = read_bool(data)
    if n == 0:
        # Empty block - read explicit location
        b = Block([], is_unreachable=is_unreachable)
        read_loc(data, b)
        expect_end_tag(data)
        return b
    else:
        # Non-empty block - read statements and set location from them
        a = read_statements(state, data, n)
        expect_end_tag(data)
        b = Block(a, is_unreachable=is_unreachable)
        b.line = a[0].line
        b.column = a[0].column
        b.end_line = a[-1].end_line
        b.end_column = a[-1].end_column
        return b


def read_optional_block(state: State, data: ReadBuffer) -> Block | None:
    expect_tag(data, nodes.BLOCK)
    expect_tag(data, LIST_GEN)
    n = read_int_bare(data)
    is_unreachable = read_bool(data)
    if n == 0:
        b = None
    else:
        a = [read_statement(state, data) for i in range(n)]
        b = Block(a, is_unreachable=is_unreachable)
        b.line = a[0].line
        b.column = a[0].column
        b.end_line = a[-1].end_line
        b.end_column = a[-1].end_column
    expect_end_tag(data)
    return b


bin_ops: Final = ["+", "-", "*", "@", "/", "%", "**", "<<", ">>", "|", "^", "&", "//"]
bool_ops: Final = ["and", "or"]
cmp_ops: Final = ["==", "!=", "<", "<=", ">", ">=", "is", "is not", "in", "not in"]
unary_ops: Final = ["~", "not", "+", "-"]


def read_expression(state: State, data: ReadBuffer) -> Expression:
    # Branches ordered by frequency (based on mypy self-check)
    tag = read_tag(data)
    expr: Expression
    if tag == nodes.NAME_EXPR:
        s = read_str(data)
        ne = NameExpr(s)
        read_loc(data, ne)
        expect_end_tag(data)
        return ne
    elif tag == nodes.MEMBER_EXPR:
        e = read_expression(state, data)
        attr = read_str(data)
        m = MemberExpr(e, attr)
        # Check if this is a super() call - if so, convert to SuperExpr
        if isinstance(e, CallExpr) and isinstance(e.callee, NameExpr) and e.callee.name == "super":
            result: Expression = SuperExpr(attr, e)
        else:
            result = m
        read_loc(data, result)
        expect_end_tag(data)
        return result
    elif tag == nodes.CALL_EXPR:
        callee = read_expression(state, data)
        args = read_expression_list(state, data)
        # Read argument kinds
        expect_tag(data, LIST_INT)
        n_kinds = read_int_bare(data)
        arg_kinds = [ARG_KINDS[read_int_bare(data)] for _ in range(n_kinds)]
        # Read argument names
        expect_tag(data, LIST_GEN)
        n_names = read_int_bare(data)
        arg_names: list[str | None] = []
        for _ in range(n_names):
            tag = read_tag(data)
            if tag == LITERAL_NONE:
                arg_names.append(None)
            elif tag == LITERAL_STR:
                arg_names.append(read_str_bare(data))
            else:
                assert False, f"Unexpected tag for arg_name: {tag}"
        ce = CallExpr(callee, args, arg_kinds, arg_names)
        read_loc(data, ce)
        expect_end_tag(data)
        return ce
    elif tag == nodes.STR_EXPR:
        se = StrExpr(read_str(data))
        read_loc(data, se)
        expect_end_tag(data)
        return se
    elif tag == nodes.COMPARISON_EXPR:
        left = read_expression(state, data)
        expect_tag(data, LIST_INT)
        n_ops = read_int_bare(data)
        ops = [cmp_ops[read_int_bare(data)] for _ in range(n_ops)]
        comparators = read_expression_list(state, data)
        assert len(ops) == len(comparators)
        expr = ComparisonExpr(ops, [left] + comparators)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.INT_EXPR:
        ie = IntExpr(read_int(data))
        read_loc(data, ie)
        expect_end_tag(data)
        return ie
    elif tag == nodes.INDEX_EXPR:
        base = read_expression(state, data)
        index = read_expression(state, data)
        expr = IndexExpr(base, index)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.LIST_EXPR:
        items = read_expression_list(state, data)
        expr = ListExpr(items)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.TUPLE_EXPR:
        items = read_expression_list(state, data)
        t = TupleExpr(items)
        read_loc(data, t)
        expect_end_tag(data)
        return t
    elif tag == nodes.BOOL_OP_EXPR:
        op = bool_ops[read_int(data)]
        values = read_expression_list(state, data)
        # Convert list of values to nested OpExpr nodes
        # E.g., [a, b, c] with "and" becomes OpExpr("and", a, OpExpr("and", b, c))
        # This matches the old parser behavior, on which we may implicitly rely.
        assert len(values) >= 2
        result = last = values[-1]
        for val in values[-2::-1]:
            result = OpExpr(op, val, result)
            result.line = val.line
            result.column = val.column
            result.end_line = last.end_line
            result.end_column = last.end_column
        read_loc(data, result)
        expect_end_tag(data)
        return result
    elif tag == nodes.UNARY_EXPR:
        op = unary_ops[read_int(data)]
        operand = read_expression(state, data)
        expr = UnaryExpr(op, operand)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.OP_EXPR:
        op = bin_ops[read_int(data)]
        left = read_expression(state, data)
        right = read_expression(state, data)
        o = OpExpr(op, left, right)
        # TODO: Store these explicitly?
        o.line = left.line
        o.column = left.column
        o.end_line = right.end_line
        o.end_column = right.end_column
        expect_end_tag(data)
        return o
    elif tag == nodes.FSTRING_EXPR:
        # F-strings are converted into nodes representing "".join([...]), to match
        # pre-existing behavior.
        nparts = read_int(data)
        fitems = []
        for _ in range(nparts):
            b = read_bool(data)
            if b:
                n = read_int(data)
                for i in range(n):
                    fitems.append(read_fstring_item(state, data))
            else:
                s = StrExpr(read_str(data))
                read_loc(data, s)
                fitems.append(s)
        expr = build_fstring_join(data, fitems)
        expect_end_tag(data)
        return expr
    elif tag == nodes.LIST_COMPREHENSION:
        generator = read_generator_expr(state, data)
        expr = ListComprehension(generator)
        read_loc(data, expr)
        set_line_column_range(generator, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.DICT_EXPR:
        expect_tag(data, LIST_GEN)
        n_keys = read_int_bare(data)
        keys: list[Expression | None] = []
        for _ in range(n_keys):
            has_key = read_bool(data)
            if has_key:
                keys.append(read_expression(state, data))
            else:
                keys.append(None)
        values = read_expression_list(state, data)
        items = list(zip(keys, values))
        expr = DictExpr(items)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.TEMP_NODE:
        temp = TempNode(AnyType(TypeOfAny.special_form), no_rhs=True)
        expect_end_tag(data)
        return temp
    elif tag == nodes.CONDITIONAL_EXPR:
        if_expr = read_expression(state, data)
        cond = read_expression(state, data)
        else_expr = read_expression(state, data)
        expr = ConditionalExpr(cond, if_expr, else_expr)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.SLICE_EXPR:
        has_begin = read_bool(data)
        begin_index = read_expression(state, data) if has_begin else None
        has_end = read_bool(data)
        end_index = read_expression(state, data) if has_end else None
        has_stride = read_bool(data)
        stride = read_expression(state, data) if has_stride else None
        expr = SliceExpr(begin_index, end_index, stride)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.GENERATOR_EXPR:
        expr = read_generator_expr(state, data)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.YIELD_EXPR:
        has_value = read_bool(data)
        if has_value:
            value = read_expression(state, data)
        else:
            value = None
        expr = YieldExpr(value)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.SET_EXPR:
        items = read_expression_list(state, data)
        expr = SetExpr(items)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.ELLIPSIS_EXPR:
        expr = EllipsisExpr()
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.TSTRING_EXPR:
        nparts = read_int(data)
        titems: list[Expression | tuple[Expression, str, str | None, Expression | None]] = []
        for _ in range(nparts):
            if read_bool(data):
                e = read_expression(state, data)
                s = read_str(data)
                if read_bool(data):
                    conv = read_str(data)
                else:
                    conv = None
                if read_bool(data):
                    # Parse format spec as a JoinedStr, this matches the old parser behavior.
                    format_spec = read_fstring_items(state, data)
                else:
                    format_spec = None
                titems.append((e, s, conv, format_spec))
            else:
                s = StrExpr(read_str(data))
                read_loc(data, s)
                titems.append(s)
        expr = TemplateStrExpr(titems)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.LAMBDA_EXPR:
        arguments, has_ann = read_parameters(state, data)
        body = read_block(state, data)

        if has_ann:
            typ = CallableType(
                [
                    arg.type_annotation if arg.type_annotation else AnyType(TypeOfAny.unannotated)
                    for arg in arguments
                ],
                [arg.kind for arg in arguments],
                [None if arg.pos_only else arg.variable.name for arg in arguments],
                AnyType(TypeOfAny.unannotated),
                _dummy_fallback,
            )
        else:
            typ = None

        expr = LambdaExpr(arguments, body)
        expr.type = typ
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.DICT_COMPREHENSION:
        key = read_expression(state, data)
        value = read_expression(state, data)
        n_generators = read_int(data)
        indices = [read_expression(state, data) for _ in range(n_generators)]
        sequences = [read_expression(state, data) for _ in range(n_generators)]
        condlists = [read_expression_list(state, data) for _ in range(n_generators)]
        is_async = [read_bool(data) for _ in range(n_generators)]
        expr = DictionaryComprehension(key, value, indices, sequences, condlists, is_async)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.SET_COMPREHENSION:
        generator = read_generator_expr(state, data)
        expr = SetComprehension(generator)
        read_loc(data, expr)
        set_line_column_range(generator, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.BYTES_EXPR:
        value = read_str(data)
        expr = BytesExpr(value)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.ASSIGNMENT_EXPR:
        target = read_expression(state, data)
        value = read_expression(state, data)
        assert isinstance(target, NameExpr), f"Expected NameExpr for target, got {type(target)}"
        expr = AssignmentExpr(target, value)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.FLOAT_EXPR:
        expect_tag(data, LITERAL_FLOAT)
        value = read_float_bare(data)
        fe = FloatExpr(value)
        read_loc(data, fe)
        expect_end_tag(data)
        return fe
    elif tag == nodes.STAR_EXPR:
        wrapped_expr = read_expression(state, data)
        expr = StarExpr(wrapped_expr)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.YIELD_FROM_EXPR:
        value = read_expression(state, data)
        expr = YieldFromExpr(value)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.AWAIT_EXPR:
        value = read_expression(state, data)
        expr = AwaitExpr(value)
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.COMPLEX_EXPR:
        expect_tag(data, LITERAL_FLOAT)
        real = read_float_bare(data)
        expect_tag(data, LITERAL_FLOAT)
        imag = read_float_bare(data)
        expr = ComplexExpr(complex(real, imag))
        read_loc(data, expr)
        expect_end_tag(data)
        return expr
    elif tag == nodes.BIG_INT_EXPR:
        strval = read_str(data)
        ie = IntExpr(int(strval, base=0))
        read_loc(data, ie)
        expect_end_tag(data)
        return ie
    else:
        assert False, tag


def read_fstring_items(state: State, data: ReadBuffer) -> Expression:
    n = read_int(data)
    items = [read_fstring_item(state, data) for _ in range(n)]
    return build_fstring_join(data, items)


def build_fstring_join(data: ReadBuffer, items: list[Expression]) -> Expression:
    items = collapse_consecutive_str_items(items)
    if len(items) == 1:
        expr = items[0]
        read_loc(data, expr)
        return expr
    args = ListExpr(items)
    str_expr = StrExpr("")
    member = MemberExpr(str_expr, "join")
    call = CallExpr(member, [args], [ARG_POS], [None])
    read_loc(data, call)
    set_line_column(args, call)
    set_line_column(str_expr, call)
    set_line_column(member, call)
    return call


def collapse_consecutive_str_items(items: list[Expression]) -> list[Expression]:
    if len(items) <= 1:
        return items
    last = items[0]
    new_items = [last]
    for item in items[1:]:
        if isinstance(last, StrExpr) and isinstance(item, StrExpr):
            last.value += item.value
            last.end_line = item.end_line
            last.end_column = item.end_column
        else:
            new_items.append(item)
            last = item
    return new_items


def read_fstring_item(state: State, data: ReadBuffer) -> Expression:
    t = read_tag(data)
    if t == LITERAL_STR:
        str_expr = StrExpr(read_str_bare(data))
        read_loc(data, str_expr)
        return str_expr
    elif t == nodes.FSTRING_INTERPOLATION:
        expr = read_expression(state, data)

        # Read conversion flag such as !r
        has_conv = read_bool(data)
        if has_conv:
            c = read_str(data)
            fmt = "{" + c + ":{}}"
        else:
            fmt = "{:{}}"

        # Read format spec such as <30 (which may have nested {...})
        has_spec = read_bool(data)
        if has_spec:
            spec = read_fstring_items(state, data)
        else:
            spec = StrExpr("")

        member = MemberExpr(StrExpr(fmt), "format")
        set_line_column(member, expr)
        call = CallExpr(member, [expr, spec], [ARG_POS, ARG_POS], [None, None])
        set_line_column(call, expr)
        expect_end_tag(data)
        return call
    else:
        raise ValueError(f"Unexpected tag {t}")


def set_line_column(target: Context, src: Context) -> None:
    target.line = src.line
    target.column = src.column


def set_line_column_range(target: Context, src: Context) -> None:
    target.line = src.line
    target.column = src.column
    target.end_line = src.end_line
    target.end_column = src.end_column


def read_expression_list(state: State, data: ReadBuffer) -> list[Expression]:
    expect_tag(data, LIST_GEN)
    n = read_int_bare(data)
    return [read_expression(state, data) for i in range(n)]


def read_generator_expr(state: State, data: ReadBuffer) -> GeneratorExpr:
    """Helper function to read comprehension data (shared by Generator, ListComp, SetComp)"""
    left_expr = read_expression(state, data)
    n_generators = read_int(data)
    indices = [read_expression(state, data) for _ in range(n_generators)]
    sequences = [read_expression(state, data) for _ in range(n_generators)]
    condlists = [read_expression_list(state, data) for _ in range(n_generators)]
    is_async = [read_bool(data) for _ in range(n_generators)]
    return GeneratorExpr(left_expr, indices, sequences, condlists, is_async)


def read_loc(data: ReadBuffer, node: Context) -> None:
    expect_tag(data, LOCATION)
    line = read_int_bare(data)
    node.line = line
    column = read_int_bare(data)
    node.column = column
    node.end_line = line + read_int_bare(data)
    node.end_column = column + read_int_bare(data)


def strip_contents_from_if_stmt(stmt: IfStmt) -> None:
    """Remove contents from IfStmt.

    Needed to still be able to check the conditions after the contents
    have been merged with the surrounding function overloads.
    """
    if len(stmt.body) == 1:
        stmt.body[0].body = []
    if stmt.else_body and len(stmt.else_body.body) == 1:
        if isinstance(stmt.else_body.body[0], IfStmt):
            strip_contents_from_if_stmt(stmt.else_body.body[0])
        else:
            stmt.else_body.body = []


def is_stripped_if_stmt(stmt: Statement) -> bool:
    """Check stmt to make sure it is a stripped IfStmt.

    See also: strip_contents_from_if_stmt
    """
    if not isinstance(stmt, IfStmt):
        return False

    if not (len(stmt.body) == 1 and len(stmt.body[0].body) == 0):
        # Body not empty
        return False

    if not stmt.else_body or len(stmt.else_body.body) == 0:
        # No or empty else_body
        return True

    # For elif, IfStmt are stored recursively in else_body
    return is_stripped_if_stmt(stmt.else_body.body[0])


def fail_merge_overload(state: State, node: IfStmt) -> None:
    """Report an error when overloads cannot be merged due to unknown condition."""
    state.add_error(
        message_registry.FAILED_TO_MERGE_OVERLOADS.value,
        node.line,
        node.column,
        blocker=False,
        code="misc",
    )


def check_ifstmt_for_overloads(
    stmt: IfStmt, current_overload_name: str | None = None
) -> str | None:
    """Check if IfStmt contains only overloads with the same name.
    Return overload_name if found, None otherwise.
    """
    # Check that block only contains a single Decorator, FuncDef, or OverloadedFuncDef.
    # Multiple overloads have already been merged as OverloadedFuncDef.
    if not (
        len(stmt.body[0].body) == 1
        and (
            isinstance(stmt.body[0].body[0], (Decorator, OverloadedFuncDef))
            or current_overload_name is not None
            and isinstance(stmt.body[0].body[0], FuncDef)
        )
        or len(stmt.body[0].body) > 1
        and isinstance(stmt.body[0].body[-1], OverloadedFuncDef)
        and all(is_stripped_if_stmt(if_stmt) for if_stmt in stmt.body[0].body[:-1])
    ):
        return None

    overload_name = cast(Decorator | FuncDef | OverloadedFuncDef, stmt.body[0].body[-1]).name
    if stmt.else_body is None or stmt.else_body.is_unreachable:
        return overload_name

    if len(stmt.else_body.body) == 1:
        # For elif: else_body contains an IfStmt itself -> do a recursive check.
        if (
            isinstance(stmt.else_body.body[0], (Decorator, FuncDef, OverloadedFuncDef))
            and stmt.else_body.body[0].name == overload_name
        ):
            return overload_name
        if (
            isinstance(stmt.else_body.body[0], IfStmt)
            and check_ifstmt_for_overloads(stmt.else_body.body[0], current_overload_name)
            == overload_name
        ):
            return overload_name

    return None


def get_executable_if_block_with_overloads(
    stmt: IfStmt, options: Options
) -> tuple[Block | None, IfStmt | None]:
    """Return block from IfStmt that will get executed.

    Return
        0 -> A block if sure that alternative blocks are unreachable.
        1 -> An IfStmt if the reachability of it can't be inferred,
             i.e. the truth value is unknown.
    """
    infer_reachability_of_if_statement(stmt, options)
    if stmt.else_body is None and stmt.body[0].is_unreachable is True:
        # always False condition with no else
        return None, None
    if (
        stmt.else_body is None
        or stmt.body[0].is_unreachable is False
        and stmt.else_body.is_unreachable is False
    ):
        # The truth value is unknown, thus not conclusive
        return None, stmt
    if stmt.else_body.is_unreachable:
        # else_body will be set unreachable if condition is always True
        return stmt.body[0], None
    if stmt.body[0].is_unreachable is True:
        # body will be set unreachable if condition is always False
        # else_body can contain an IfStmt itself (for elif) -> do a recursive check
        if isinstance(stmt.else_body.body[0], IfStmt):
            return get_executable_if_block_with_overloads(stmt.else_body.body[0], options)
        return stmt.else_body, None
    return None, stmt


def fix_function_overloads(state: State, stmts: list[Statement]) -> list[Statement]:
    """Merge consecutive function overloads into OverloadedFuncDef nodes.

    Also handles conditional overloads (overloads inside if statements).
    """
    ret: list[Statement] = []
    current_overload: list[OverloadPart] = []
    current_overload_name: str | None = None
    last_unconditional_func_def: str | None = None
    last_if_stmt: IfStmt | None = None
    last_if_overload: Decorator | FuncDef | OverloadedFuncDef | None = None
    last_if_stmt_overload_name: str | None = None
    last_if_unknown_truth_value: IfStmt | None = None
    skipped_if_stmts: list[IfStmt] = []
    for stmt in stmts:
        if_overload_name: str | None = None
        if_block_with_overload: Block | None = None
        if_unknown_truth_value: IfStmt | None = None
        if isinstance(stmt, IfStmt):
            # Check IfStmt block to determine if function overloads can be merged
            if_overload_name = check_ifstmt_for_overloads(stmt, current_overload_name)
            if if_overload_name is not None:
                if_block_with_overload, if_unknown_truth_value = (
                    get_executable_if_block_with_overloads(stmt, state.options)
                )

        if (
            current_overload_name is not None
            and isinstance(stmt, (Decorator, FuncDef))
            and stmt.name == current_overload_name
        ):
            if last_if_stmt is not None:
                skipped_if_stmts.append(last_if_stmt)
            if last_if_overload is not None:
                # Last stmt was an IfStmt with same overload name
                # Add overloads to current_overload
                if isinstance(last_if_overload, OverloadedFuncDef):
                    current_overload.extend(last_if_overload.items)
                else:
                    current_overload.append(last_if_overload)
                last_if_stmt, last_if_overload = None, None
            if last_if_unknown_truth_value:
                fail_merge_overload(state, last_if_unknown_truth_value)
                last_if_unknown_truth_value = None
            current_overload.append(stmt)
            if isinstance(stmt, FuncDef):
                # This is, strictly speaking, wrong: there might be a decorated
                # implementation. However, it only affects the error message we show:
                # ideally it's "already defined", but "implementation must come last"
                # is also reasonable.
                # TODO: can we get rid of this completely and just always emit
                # "implementation must come last" instead?
                last_unconditional_func_def = stmt.name
        elif (
            current_overload_name is not None
            and isinstance(stmt, IfStmt)
            and if_overload_name == current_overload_name
            and last_unconditional_func_def != current_overload_name
        ):
            # IfStmt only contains stmts relevant to current_overload.
            # Check if stmts are reachable and add them to current_overload,
            # otherwise skip IfStmt to allow subsequent overload
            # or function definitions.
            skipped_if_stmts.append(stmt)
            if if_block_with_overload is None:
                if if_unknown_truth_value is not None:
                    fail_merge_overload(state, if_unknown_truth_value)
                continue
            if last_if_overload is not None:
                # Last stmt was an IfStmt with same overload name
                # Add overloads to current_overload
                if isinstance(last_if_overload, OverloadedFuncDef):
                    current_overload.extend(last_if_overload.items)
                else:
                    current_overload.append(last_if_overload)
                last_if_stmt, last_if_overload = None, None
            if isinstance(if_block_with_overload.body[-1], OverloadedFuncDef):
                skipped_if_stmts.extend(cast(list[IfStmt], if_block_with_overload.body[:-1]))
                current_overload.extend(if_block_with_overload.body[-1].items)
            else:
                current_overload.append(cast(Decorator | FuncDef, if_block_with_overload.body[0]))
        else:
            if last_if_stmt is not None:
                ret.append(last_if_stmt)
                last_if_stmt_overload_name = current_overload_name
                last_if_stmt, last_if_overload = None, None
                last_if_unknown_truth_value = None

            if current_overload and current_overload_name == last_if_stmt_overload_name:
                # Remove last stmt (IfStmt) from ret if the overload names matched
                # Only happens if no executable block had been found in IfStmt
                popped = ret.pop()
                assert isinstance(popped, IfStmt)
                skipped_if_stmts.append(popped)
            if current_overload and skipped_if_stmts:
                # Add bare IfStmt (without overloads) to ret
                # Required for mypy to be able to still check conditions
                for if_stmt in skipped_if_stmts:
                    strip_contents_from_if_stmt(if_stmt)
                    ret.append(if_stmt)
                skipped_if_stmts = []
            if len(current_overload) == 1:
                ret.append(current_overload[0])
            elif len(current_overload) > 1:
                ret.append(OverloadedFuncDef(current_overload))

            # If we have multiple decorated functions named "_" next to each, we want to treat
            # them as a series of regular FuncDefs instead of one OverloadedFuncDef because
            # most of mypy/mypyc assumes that all the functions in an OverloadedFuncDef are
            # related, but multiple underscore functions next to each other aren't necessarily
            # related
            last_unconditional_func_def = None
            if isinstance(stmt, Decorator) and not unnamed_function(stmt.name):
                current_overload = [stmt]
                current_overload_name = stmt.name
            elif isinstance(stmt, IfStmt) and if_overload_name is not None:
                current_overload = []
                current_overload_name = if_overload_name
                last_if_stmt = stmt
                last_if_stmt_overload_name = None
                if if_block_with_overload is not None:
                    skipped_if_stmts.extend(cast(list[IfStmt], if_block_with_overload.body[:-1]))
                    last_if_overload = cast(
                        Decorator | FuncDef | OverloadedFuncDef, if_block_with_overload.body[-1]
                    )
                last_if_unknown_truth_value = if_unknown_truth_value
            else:
                current_overload = []
                current_overload_name = None
                ret.append(stmt)

    if current_overload and skipped_if_stmts:
        # Add bare IfStmt (without overloads) to ret
        # Required for mypy to be able to still check conditions
        for if_stmt in skipped_if_stmts:
            strip_contents_from_if_stmt(if_stmt)
            ret.append(if_stmt)
    if len(current_overload) == 1:
        ret.append(current_overload[0])
    elif len(current_overload) > 1:
        ret.append(OverloadedFuncDef(current_overload))
    elif last_if_overload is not None:
        ret.append(last_if_overload)
    elif last_if_stmt is not None:
        ret.append(last_if_stmt)
    return ret


def deserialize_imports(import_bytes: bytes) -> list[ImportBase]:
    """Deserialize import metadata from bytes into mypy AST nodes."""
    if not import_bytes:
        return []

    data = ReadBuffer(import_bytes)

    expect_tag(data, LIST_GEN)
    n_imports = read_int_bare(data)

    imports: list[ImportBase] = []

    for _ in range(n_imports):
        tag = read_tag(data)

        if tag == IMPORT_METADATA:
            name = read_str(data)
            relative = read_int(data)

            has_asname = read_bool(data)
            if has_asname:
                asname = read_str(data)
            else:
                asname = None

            # Note: relative imports are handled via ImportFrom, so relative should be 0 here
            stmt = Import([(name, asname)])
            _read_and_set_import_metadata(data, stmt)
            imports.append(stmt)

        elif tag == IMPORTFROM_METADATA:
            module = read_str(data)
            relative = read_int(data)

            expect_tag(data, LIST_GEN)
            n_names = read_int_bare(data)
            names: list[tuple[str, str | None]] = []

            for _ in range(n_names):
                name = read_str(data)
                has_asname = read_bool(data)
                if has_asname:
                    asname = read_str(data)
                else:
                    asname = None
                names.append((name, asname))

            stmt = ImportFrom(module, relative, names)
            _read_and_set_import_metadata(data, stmt)
            imports.append(stmt)

        elif tag == IMPORTALL_METADATA:
            module = read_str(data)
            relative = read_int(data)

            stmt = ImportAll(module, relative)
            _read_and_set_import_metadata(data, stmt)
            imports.append(stmt)

        else:
            raise ValueError(f"Unexpected tag in import metadata: {tag}")

    return imports


def _read_and_set_import_metadata(data: ReadBuffer, stmt: Import | ImportFrom | ImportAll) -> None:
    read_loc(data, stmt)

    # Metadata flags as a single integer bitfield
    flags = read_int(data)

    # Extract individual flags using bitwise operations
    # Bit 0: is_top_level
    # Bit 1: is_unreachable
    # Bit 2: is_mypy_only
    stmt.is_top_level = (flags & 0x01) != 0
    stmt.is_unreachable = (flags & 0x02) != 0
    stmt.is_mypy_only = (flags & 0x04) != 0
