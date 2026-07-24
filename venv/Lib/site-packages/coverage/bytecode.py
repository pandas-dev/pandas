# Licensed under the Apache License: http://www.apache.org/licenses/LICENSE-2.0
# For details: https://github.com/coveragepy/coveragepy/blob/main/NOTICE.txt

"""Bytecode analysis for coverage.py"""

from __future__ import annotations

import dis
from collections.abc import Iterable, Mapping
from types import CodeType

from coverage.types import TArc, TLineNo, TOffset


class ByteParser:
    """Parse bytecode to understand the structure of code."""

    def __init__(
        self,
        *,
        code: CodeType | None = None,
        text: str | None = None,
        filename: str | None = None,
    ) -> None:
        if code is None:
            assert text is not None
            code = compile(text, filename or "<string>", "exec", dont_inherit=True)
        self.code = code

    def _child_parsers(self) -> Iterable[ByteParser]:
        """Iterate over all the code objects nested within this one.

        The iteration includes `self` as its first value.

        We skip code objects named `__annotate__` since they are deferred
        annotations that usually are never run.  If there are errors in the
        annotations, they will be caught by type checkers or other tools that
        use annotations.

        """
        return (ByteParser(code=c) for c in self.code_objects() if c.co_name != "__annotate__")

    def code_objects(self) -> Iterable[CodeType]:
        """Iterate over all the code objects in `code`."""
        stack = [self.code]
        while stack:
            # We're going to return the code object on the stack, but first
            # push its children for later returning.
            code = stack.pop()
            for c in code.co_consts:
                if isinstance(c, CodeType):
                    stack.append(c)
            yield code

    def _line_numbers(self) -> Iterable[TLineNo]:
        """Yield the line numbers possible in this code object.

        Uses co_lines() to produce a sequence: l0, l1, ...
        """
        for _, _, line in self.code.co_lines():
            if line:
                yield line

    def find_statements(self) -> Iterable[TLineNo]:
        """Find the statements in `self.code`.

        Produce a sequence of line numbers that start statements.  Recurses
        into all code objects reachable from `self.code`.

        """
        for bp in self._child_parsers():
            # Get all of the lineno information from this code.
            yield from bp._line_numbers()


def bytes_to_lines(code: CodeType) -> dict[TOffset, TLineNo]:
    """Make a dict mapping byte code offsets to line numbers."""
    b2l = {}
    for bstart, bend, lineno in code.co_lines():
        if lineno is not None:
            for boffset in range(bstart, bend, 2):
                b2l[boffset] = lineno
    return b2l


def op_set(*op_names: str) -> set[int]:
    """Make a set of opcodes from instruction names.

    The names might not exist in this version of Python, skip those if not.
    """
    ops = {op for name in op_names if (op := dis.opmap.get(name))}
    assert ops, f"At least one opcode must exist: {op_names}"
    return ops


# Opcodes that are unconditional jumps elsewhere.
ALWAYS_JUMPS = op_set(
    "JUMP_BACKWARD",
    "JUMP_BACKWARD_NO_INTERRUPT",
    "JUMP_FORWARD",
)

# Opcodes that exit from a function.
RETURNS = op_set(
    "RETURN_VALUE",
    "RETURN_GENERATOR",
)


# CACHE doesn't exist in Python 3.10, but the branch resolver is only used
# on 3.14+, so a placeholder value is fine.
_CACHE = dis.opmap.get("CACHE", -1)
_EXTENDED_ARG = dis.opmap["EXTENDED_ARG"]

# All opcodes with a jump target.
JUMPS = set(dis.hasjrel) | set(dis.hasjabs)

# Opcodes that jump backwards.
BACKWARD_JUMPS = {op for op in JUMPS if "JUMP_BACKWARD" in dis.opname[op]}


class BranchArcResolver:
    """Resolve branch events to line arcs, one (source, dest) pair at a time.

    Branch events are one-shot (they are DISABLEd after firing), so each
    (source offset, destination offset) pair is resolved at most a couple of
    times per code object.  Resolving pairs on demand is much cheaper than
    precomputing trails for every branch in the code object, most of which
    never fire.  We walk the raw bytecode bytes so that we never need to
    disassemble whole code objects with `dis`.

    To resolve one pair:
    starting from the destination, follow the trail of instructions (through
    unconditional jumps) until we reach an instruction on a new source line
    (giving us the arc), a return (an arc to leaving the code object), or
    another branch possibility (no arc: that branch will produce its own
    events).

    """

    def __init__(
        self,
        code: CodeType,
        byte_to_line: Mapping[TOffset, TLineNo],
        multiline_map: Mapping[TLineNo, TLineNo],
    ) -> None:
        self.code = code
        # co_code re-copies the bytes on each access, so fetch it once.
        self.co_code = code.co_code
        self.byte_to_line = byte_to_line
        self.multiline_map = multiline_map

    def line_at(self, offset: TOffset) -> TLineNo | None:
        """The source line of the instruction at `offset`, de-multilined."""
        line = self.byte_to_line.get(offset)
        if line is not None:
            line = self.multiline_map.get(line, line)
        return line

    def resolve(self, source: TOffset, dest: TOffset) -> TArc | None:
        """Turn a branch event's (source, dest) offsets into an arc, or None."""
        from_line = self.line_at(source)
        if from_line is None:
            return None
        co_code = self.co_code
        max_offset = len(co_code)
        byte_to_line = self.byte_to_line
        multiline_map = self.multiline_map
        offset = dest
        ext_arg = 0
        seen: set[TOffset] = set()
        while 0 <= offset < max_offset and offset not in seen:
            seen.add(offset)
            op = co_code[offset]
            if op == _CACHE:
                offset += 2
                continue
            if op == _EXTENDED_ARG:
                ext_arg = (ext_arg | co_code[offset + 1]) << 8
                offset += 2
                continue
            line = byte_to_line.get(offset)
            if line is not None:
                line = multiline_map.get(line, line)
                if line and line != from_line:
                    return (from_line, line)
            if op in JUMPS:
                if op in ALWAYS_JUMPS:
                    arg = ext_arg | co_code[offset + 1]
                    # Jump distances are measured from the end of the
                    # instruction's inline CACHE entries, which appear in
                    # co_code as CACHE opcodes immediately following it.
                    next_offset = offset + 2
                    while next_offset < max_offset and co_code[next_offset] == _CACHE:
                        next_offset += 2
                    if op in BACKWARD_JUMPS:
                        offset = next_offset - 2 * arg
                    else:
                        offset = next_offset + 2 * arg
                    ext_arg = 0
                    continue
                # Another branch possibility: it will get its own events.
                return None
            if op in RETURNS:
                return (from_line, -self.code.co_firstlineno)
            ext_arg = 0
            offset += 2
        return None
