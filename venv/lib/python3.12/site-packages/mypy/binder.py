from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from typing import NamedTuple, Optional, Union
from typing_extensions import TypeAlias as _TypeAlias

from mypy.erasetype import remove_instance_last_known_values
from mypy.literals import Key, extract_var_from_literal_hash, literal, literal_hash, subkeys
from mypy.nodes import Expression, IndexExpr, MemberExpr, NameExpr, RefExpr, TypeInfo, Var
from mypy.options import Options
from mypy.subtypes import is_same_type, is_subtype
from mypy.typeops import make_simplified_union
from mypy.types import (
    AnyType,
    Instance,
    NoneType,
    PartialType,
    ProperType,
    TupleType,
    Type,
    TypeOfAny,
    TypeType,
    TypeVarType,
    UnionType,
    UnpackType,
    find_unpack_in_list,
    get_proper_type,
)
from mypy.typevars import fill_typevars_with_any

BindableExpression: _TypeAlias = Union[IndexExpr, MemberExpr, NameExpr]


class CurrentType(NamedTuple):
    type: Type
    from_assignment: bool


class Frame:
    """A Frame represents a specific point in the execution of a program.

    It carries information about the current types of expressions at
    that point, arising either from assignments to those expressions
    or the result of isinstance checks and other type narrowing
    operations. It also records whether it is possible to reach that
    point at all.

    We add a new frame wherenever there is a new scope or control flow
    branching.

    This information is not copied into a new Frame when it is pushed
    onto the stack, so a given Frame only has information about types
    that were assigned in that frame.

    Expressions are stored in dicts using 'literal hashes' as keys (type
    "Key"). These are hashable values derived from expression AST nodes
    (only those that can be narrowed). literal_hash(expr) is used to
    calculate the hashes. Note that this isn't directly related to literal
    types -- the concept predates literal types.
    """

    def __init__(self, id: int, conditional_frame: bool = False) -> None:
        self.id = id
        self.types: dict[Key, CurrentType] = {}
        self.unreachable = False
        self.conditional_frame = conditional_frame
        self.suppress_unreachable_warnings = False

    def __repr__(self) -> str:
        return f"Frame({self.id}, {self.types}, {self.unreachable}, {self.conditional_frame})"


Assigns = defaultdict[Expression, list[tuple[Type, Optional[Type]]]]


class ConditionalTypeBinder:
    """Keep track of conditional types of variables.

    NB: Variables are tracked by literal hashes of expressions, so it is
    possible to confuse the binder when there is aliasing. Example:

        class A:
            a: int | str

        x = A()
        lst = [x]
        reveal_type(x.a)      # int | str
        x.a = 1
        reveal_type(x.a)      # int
        reveal_type(lst[0].a) # int | str
        lst[0].a = 'a'
        reveal_type(x.a)      # int
        reveal_type(lst[0].a) # str
    """

    # Stored assignments for situations with tuple/list lvalue and rvalue of union type.
    # This maps an expression to a list of bound types for every item in the union type.
    type_assignments: Assigns | None = None

    def __init__(self, options: Options) -> None:
        # Each frame gets an increasing, distinct id.
        self.next_id = 1

        # The stack of frames currently used.  These map
        # literal_hash(expr) -- literals like 'foo.bar' --
        # to types. The last element of this list is the
        # top-most, current frame. Each earlier element
        # records the state as of when that frame was last
        # on top of the stack.
        self.frames = [Frame(self._get_id())]

        # For frames higher in the stack, we record the set of
        # Frames that can escape there, either by falling off
        # the end of the frame or by a loop control construct
        # or raised exception. The last element of self.frames
        # has no corresponding element in this list.
        self.options_on_return: list[list[Frame]] = []

        # Maps literal_hash(expr) to get_declaration(expr)
        # for every expr stored in the binder
        self.declarations: dict[Key, Type | None] = {}
        # Set of other keys to invalidate if a key is changed, e.g. x -> {x.a, x[0]}
        # Whenever a new key (e.g. x.a.b) is added, we update this
        self.dependencies: dict[Key, set[Key]] = {}

        # Whether the last pop changed the newly top frame on exit
        self.last_pop_changed = False

        # These are used to track control flow in try statements and loops.
        self.try_frames: set[int] = set()
        self.break_frames: list[int] = []
        self.continue_frames: list[int] = []

        # If True, initial assignment to a simple variable (e.g. "x", but not "x.y")
        # is added to the binder. This allows more precise narrowing and more
        # flexible inference of variable types (--allow-redefinition-new).
        self.bind_all = options.allow_redefinition_new

    def _get_id(self) -> int:
        self.next_id += 1
        return self.next_id

    def _add_dependencies(self, key: Key, value: Key | None = None) -> None:
        if value is None:
            value = key
        else:
            self.dependencies.setdefault(key, set()).add(value)
        for elt in subkeys(key):
            self._add_dependencies(elt, value)

    def push_frame(self, conditional_frame: bool = False) -> Frame:
        """Push a new frame into the binder."""
        f = Frame(self._get_id(), conditional_frame)
        self.frames.append(f)
        self.options_on_return.append([])
        return f

    def _put(self, key: Key, type: Type, from_assignment: bool, index: int = -1) -> None:
        self.frames[index].types[key] = CurrentType(type, from_assignment)

    def _get(self, key: Key, index: int = -1) -> CurrentType | None:
        if index < 0:
            index += len(self.frames)
        for i in range(index, -1, -1):
            if key in self.frames[i].types:
                return self.frames[i].types[key]
        return None

    def put(self, expr: Expression, typ: Type, *, from_assignment: bool = True) -> None:
        """Directly set the narrowed type of expression (if it supports it).

        This is used for isinstance() etc. Assignments should go through assign_type().
        """
        if not isinstance(expr, (IndexExpr, MemberExpr, NameExpr)):
            return
        if not literal(expr):
            return
        key = literal_hash(expr)
        assert key is not None, "Internal error: binder tried to put non-literal"
        if key not in self.declarations:
            self.declarations[key] = get_declaration(expr)
            self._add_dependencies(key)
        self._put(key, typ, from_assignment)

    def unreachable(self) -> None:
        self.frames[-1].unreachable = True

    def suppress_unreachable_warnings(self) -> None:
        self.frames[-1].suppress_unreachable_warnings = True

    def get(self, expr: Expression) -> Type | None:
        key = literal_hash(expr)
        assert key is not None, "Internal error: binder tried to get non-literal"
        found = self._get(key)
        if found is None:
            return None
        return found.type

    def is_unreachable(self) -> bool:
        # TODO: Copy the value of unreachable into new frames to avoid
        # this traversal on every statement?
        return any(f.unreachable for f in self.frames)

    def is_unreachable_warning_suppressed(self) -> bool:
        return any(f.suppress_unreachable_warnings for f in self.frames)

    def cleanse(self, expr: Expression) -> None:
        """Remove all references to a Node from the binder."""
        key = literal_hash(expr)
        assert key is not None, "Internal error: binder tried cleanse non-literal"
        self._cleanse_key(key)

    def _cleanse_key(self, key: Key) -> None:
        """Remove all references to a key from the binder."""
        for frame in self.frames:
            if key in frame.types:
                del frame.types[key]

    def update_from_options(self, frames: list[Frame]) -> bool:
        """Update the frame to reflect that each key will be updated
        as in one of the frames.  Return whether any item changes.

        If a key is declared as AnyType, only update it if all the
        options are the same.
        """
        all_reachable = all(not f.unreachable for f in frames)
        frames = [f for f in frames if not f.unreachable]
        changed = False
        keys = {key for f in frames for key in f.types}

        for key in keys:
            current_value = self._get(key)
            resulting_values = [f.types.get(key, current_value) for f in frames]
            # Keys can be narrowed using two different semantics. The new semantics
            # is enabled for plain variables when bind_all is true, and it allows
            # variable types to be widened using subsequent assignments. This is
            # tricky to support for instance attributes (primarily due to deferrals),
            # so we don't use it for them.
            old_semantics = not self.bind_all or extract_var_from_literal_hash(key) is None
            if old_semantics and any(x is None for x in resulting_values):
                # We didn't know anything about key before
                # (current_value must be None), and we still don't
                # know anything about key in at least one possible frame.
                continue

            resulting_values = [x for x in resulting_values if x is not None]

            if all_reachable and all(
                x is not None and not x.from_assignment for x in resulting_values
            ):
                # Do not synthesize a new type if we encountered a conditional block
                # (if, while or match-case) without assignments.
                # See check-isinstance.test::testNoneCheckDoesNotMakeTypeVarOptional
                # This is a safe assumption: the fact that we checked something with `is`
                # or `isinstance` does not change the type of the value.
                continue

            current_type = resulting_values[0]
            assert current_type is not None
            type = current_type.type
            declaration_type = get_proper_type(self.declarations.get(key))
            if isinstance(declaration_type, AnyType):
                # At this point resulting values can't contain None, see continue above
                if not all(
                    t is not None and is_same_type(type, t.type) for t in resulting_values[1:]
                ):
                    type = AnyType(TypeOfAny.from_another_any, source_any=declaration_type)
            else:
                possible_types = []
                for t in resulting_values:
                    assert t is not None
                    possible_types.append(t.type)
                if len(possible_types) == 1:
                    # This is to avoid calling get_proper_type() unless needed, as this may
                    # interfere with our (hacky) TypeGuard support.
                    type = possible_types[0]
                else:
                    type = make_simplified_union(possible_types)
                    # Legacy guard for corner case when the original type is TypeVarType.
                    if isinstance(declaration_type, TypeVarType) and not is_subtype(
                        type, declaration_type
                    ):
                        type = declaration_type
                    # Try simplifying resulting type for unions involving variadic tuples.
                    # Technically, everything is still valid without this step, but if we do
                    # not do this, this may create long unions after exiting an if check like:
                    #     x: tuple[int, ...]
                    #     if len(x) < 10:
                    #         ...
                    # We want the type of x to be tuple[int, ...] after this block (if it is
                    # still equivalent to such type).
                    if isinstance(type, UnionType):
                        type = collapse_variadic_union(type)
                    if (
                        old_semantics
                        and isinstance(type, ProperType)
                        and isinstance(type, UnionType)
                    ):
                        # Simplify away any extra Any's that were added to the declared
                        # type when popping a frame.
                        simplified = UnionType.make_union(
                            [t for t in type.items if not isinstance(get_proper_type(t), AnyType)]
                        )
                        if simplified == self.declarations[key]:
                            type = simplified
            if current_value is None or not is_same_type(type, current_value.type):
                self._put(key, type, from_assignment=True)
                changed = True

        self.frames[-1].unreachable = not frames

        return changed

    def pop_frame(self, can_skip: bool, fall_through: int) -> Frame:
        """Pop a frame and return it.

        See frame_context() for documentation of fall_through.
        """

        if fall_through > 0:
            self.allow_jump(-fall_through)

        result = self.frames.pop()
        options = self.options_on_return.pop()

        if can_skip:
            options.insert(0, self.frames[-1])

        self.last_pop_changed = self.update_from_options(options)

        return result

    @contextmanager
    def accumulate_type_assignments(self) -> Iterator[Assigns]:
        """Push a new map to collect assigned types in multiassign from union.

        If this map is not None, actual binding is deferred until all items in
        the union are processed (a union of collected items is later bound
        manually by the caller).
        """
        old_assignments = None
        if self.type_assignments is not None:
            old_assignments = self.type_assignments
        self.type_assignments = defaultdict(list)
        yield self.type_assignments
        self.type_assignments = old_assignments

    def assign_type(self, expr: Expression, type: Type, declared_type: Type | None) -> None:
        """Narrow type of expression through an assignment.

        Do nothing if the expression doesn't support narrowing.

        When not narrowing though an assignment (isinstance() etc.), use put()
        directly. This omits some special-casing logic for assignments.
        """
        # We should erase last known value in binder, because if we are using it,
        # it means that the target is not final, and therefore can't hold a literal.
        type = remove_instance_last_known_values(type)

        if self.type_assignments is not None:
            # We are in a multiassign from union, defer the actual binding,
            # just collect the types.
            self.type_assignments[expr].append((type, declared_type))
            return
        if not isinstance(expr, (IndexExpr, MemberExpr, NameExpr)):
            return
        if not literal(expr):
            return
        self.invalidate_dependencies(expr)

        if declared_type is None:
            # Not sure why this happens.  It seems to mainly happen in
            # member initialization.
            return
        if not is_subtype(type, declared_type):
            # Pretty sure this is only happens when there's a type error.

            # Ideally this function wouldn't be called if the
            # expression has a type error, though -- do other kinds of
            # errors cause this function to get called at invalid
            # times?
            return

        p_declared = get_proper_type(declared_type)
        p_type = get_proper_type(type)
        if isinstance(p_type, AnyType):
            # Any type requires some special casing, for both historical reasons,
            # and to optimise user experience without sacrificing correctness too much.
            if isinstance(expr, RefExpr) and isinstance(expr.node, Var) and expr.node.is_inferred:
                # First case: a local/global variable without explicit annotation,
                # in this case we just assign Any (essentially following the SSA logic).
                self.put(expr, type)
            elif isinstance(p_declared, UnionType) and any(
                isinstance(get_proper_type(item), NoneType) for item in p_declared.items
            ):
                # Second case: explicit optional type, in this case we optimize for a common
                # pattern when an untyped value used as a fallback replacing None.
                new_items = [
                    type if isinstance(get_proper_type(item), NoneType) else item
                    for item in p_declared.items
                ]
                self.put(expr, UnionType(new_items))
            elif isinstance(p_declared, UnionType) and any(
                isinstance(get_proper_type(item), AnyType) for item in p_declared.items
            ):
                # Third case: a union already containing Any (most likely from an un-imported
                # name), in this case we allow assigning Any as well.
                self.put(expr, type)
            else:
                # In all other cases we don't narrow to Any to minimize false negatives.
                self.put(expr, declared_type)
        elif isinstance(p_declared, AnyType):
            # Mirroring the first case above, we don't narrow to a precise type if the variable
            # has an explicit `Any` type annotation.
            if isinstance(expr, RefExpr) and isinstance(expr.node, Var) and expr.node.is_inferred:
                self.put(expr, type)
            else:
                self.put(expr, declared_type)
        else:
            self.put(expr, type)

        for i in self.try_frames:
            # XXX This should probably not copy the entire frame, but
            # just copy this variable into a single stored frame.
            self.allow_jump(i)

    def invalidate_dependencies(self, expr: BindableExpression) -> None:
        """Invalidate knowledge of types that include expr, but not expr itself.

        For example, when expr is foo.bar, invalidate foo.bar.baz.

        It is overly conservative: it invalidates globally, including
        in code paths unreachable from here.
        """
        key = literal_hash(expr)
        assert key is not None
        for dep in self.dependencies.get(key, set()):
            self._cleanse_key(dep)

    def allow_jump(self, index: int) -> None:
        # self.frames and self.options_on_return have different lengths
        # so make sure the index is positive
        if index < 0:
            index += len(self.options_on_return)
        frame = Frame(self._get_id())
        for f in self.frames[index + 1 :]:
            frame.types.update(f.types)
            if f.unreachable:
                frame.unreachable = True
        self.options_on_return[index].append(frame)

    def handle_break(self) -> None:
        self.allow_jump(self.break_frames[-1])
        self.unreachable()

    def handle_continue(self) -> None:
        self.allow_jump(self.continue_frames[-1])
        self.unreachable()

    @contextmanager
    def frame_context(
        self,
        *,
        can_skip: bool,
        fall_through: int = 1,
        break_frame: int = 0,
        continue_frame: int = 0,
        conditional_frame: bool = False,
        try_frame: bool = False,
    ) -> Iterator[Frame]:
        """Return a context manager that pushes/pops frames on enter/exit.

        If can_skip is True, control flow is allowed to bypass the
        newly-created frame.

        If fall_through > 0, then it will allow control flow that
        falls off the end of the frame to escape to its ancestor
        `fall_through` levels higher. Otherwise control flow ends
        at the end of the frame.

        If break_frame > 0, then 'break' statements within this frame
        will jump out to the frame break_frame levels higher than the
        frame created by this call to frame_context. Similarly for
        continue_frame and 'continue' statements.

        If try_frame is true, then execution is allowed to jump at any
        point within the newly created frame (or its descendants) to
        its parent (i.e., to the frame that was on top before this
        call to frame_context).

        After the context manager exits, self.last_pop_changed indicates
        whether any types changed in the newly-topmost frame as a result
        of popping this frame.
        """
        assert len(self.frames) > 1

        if break_frame:
            self.break_frames.append(len(self.frames) - break_frame)
        if continue_frame:
            self.continue_frames.append(len(self.frames) - continue_frame)
        if try_frame:
            self.try_frames.add(len(self.frames) - 1)

        new_frame = self.push_frame(conditional_frame)
        if try_frame:
            # An exception may occur immediately
            self.allow_jump(-1)
        yield new_frame
        self.pop_frame(can_skip, fall_through)

        if break_frame:
            self.break_frames.pop()
        if continue_frame:
            self.continue_frames.pop()
        if try_frame:
            self.try_frames.remove(len(self.frames) - 1)

    @contextmanager
    def top_frame_context(self) -> Iterator[Frame]:
        """A variant of frame_context for use at the top level of
        a namespace (module, function, or class).
        """
        assert len(self.frames) == 1
        yield self.push_frame()
        self.pop_frame(True, 0)
        assert len(self.frames) == 1


def get_declaration(expr: BindableExpression) -> Type | None:
    """Get the declared or inferred type of a RefExpr expression.

    Return None if there is no type or the expression is not a RefExpr.
    This can return None if the type hasn't been inferred yet.
    """
    if isinstance(expr, RefExpr):
        if isinstance(expr.node, Var):
            type = expr.node.type
            if not isinstance(get_proper_type(type), PartialType):
                return type
        elif isinstance(expr.node, TypeInfo):
            return TypeType(fill_typevars_with_any(expr.node))
    return None


def collapse_variadic_union(typ: UnionType) -> Type:
    """Simplify a union involving variadic tuple if possible.

    This will collapse a type like e.g.
        tuple[X, Z] | tuple[X, Y, Z] | tuple[X, Y, Y, *tuple[Y, ...], Z]
    back to
        tuple[X, *tuple[Y, ...], Z]
    which is equivalent, but much simpler form of the same type.
    """
    tuple_items = []
    other_items = []
    for t in typ.items:
        p_t = get_proper_type(t)
        if isinstance(p_t, TupleType):
            tuple_items.append(p_t)
        else:
            other_items.append(t)
    if len(tuple_items) <= 1:
        # This type cannot be simplified further.
        return typ
    tuple_items = sorted(tuple_items, key=lambda t: len(t.items))
    first = tuple_items[0]
    last = tuple_items[-1]
    unpack_index = find_unpack_in_list(last.items)
    if unpack_index is None:
        return typ
    unpack = last.items[unpack_index]
    assert isinstance(unpack, UnpackType)
    unpacked = get_proper_type(unpack.type)
    if not isinstance(unpacked, Instance):
        return typ
    assert unpacked.type.fullname == "builtins.tuple"
    suffix = last.items[unpack_index + 1 :]

    # Check that first item matches the expected pattern and infer prefix.
    if len(first.items) < len(suffix):
        return typ
    if suffix and first.items[-len(suffix) :] != suffix:
        return typ
    if suffix:
        prefix = first.items[: -len(suffix)]
    else:
        prefix = first.items

    # Check that all middle types match the expected pattern as well.
    arg = unpacked.args[0]
    for i, it in enumerate(tuple_items[1:-1]):
        if it.items != prefix + [arg] * (i + 1) + suffix:
            return typ

    # Check the last item (the one with unpack), and choose an appropriate simplified type.
    if last.items != prefix + [arg] * (len(typ.items) - 1) + [unpack] + suffix:
        return typ
    if len(first.items) == 0:
        simplified: Type = unpacked.copy_modified()
    else:
        simplified = TupleType(prefix + [unpack] + suffix, fallback=last.partial_fallback)
    return UnionType.make_union([simplified] + other_items)
