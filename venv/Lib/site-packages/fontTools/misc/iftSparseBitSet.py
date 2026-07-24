"""Sparse Bit Set encoding/decoding for IFT (Incremental Font Transfer).

Implements the sparse bit set format defined in the W3C IFT specification:
https://w3c.github.io/IFT/Overview.html#sparse-bit-set-decoding
"""

# Reference: https://github.com/googlefonts/fontations/blob/main/read-fonts/src/collections/int_set/sparse_bit_set.rs


import bisect
from collections import deque
from typing import Dict, Iterable, Optional, Set, Tuple

# Maximum tree height that fits within a 32-bit leaf node for each branching factor.
_BF_MAX_HEIGHT: Dict[int, int] = {2: 31, 4: 16, 8: 11, 32: 7}


class SparseBitSetDecodeError(Exception):
    pass


def decode(
    data: bytes, bias: int = 0, maxValue: int = 0xFFFFFFFF
) -> Tuple[Set[int], int]:
    """Decode a sparse bit set from binary data.

    Args:
        data: bytes-like object containing the sparse bit set encoding.
        bias: integer added to each decoded value.

    Returns:
        A tuple (values, bytesConsumed) where values is a set of integers
        and bytesConsumed is the number of bytes read from data.
    """
    if not data:
        raise SparseBitSetDecodeError("Empty data")

    branchFactor, height = _decodeHeader(data[0])

    maxHeight = _BF_MAX_HEIGHT[branchFactor]
    if height > maxHeight:
        raise SparseBitSetDecodeError(
            f"Height {height} exceeds max {maxHeight} for branch factor {branchFactor}"
        )

    return _decodeImpl(data, branchFactor, height, bias, maxValue)


def _decodeImpl(
    data: bytes, branchFactor: int, height: int, bias: int, maxValue: int
) -> Tuple[Set[int], int]:
    if height == 0:
        # 1 byte was used for the header.
        return (set(), 1)

    bitStream = _InputBitStream(data, branchFactor)
    result: Set[int] = set()
    # Queue entries are (startValue, depth), where startValue is the first
    # integer that could be covered by this node.
    queue: deque[Tuple[int, int]] = deque()
    queue.append((0, 1))

    while queue:
        start, depth = queue.popleft()
        bits = bitStream.next()
        if bits is None:
            raise SparseBitSetDecodeError("Unexpected end of data")

        # all bits were were zero which is a special command to completely fill
        # in all integers covered by this node.
        if bits == 0:
            exp = height - depth + 1
            nodeSize = branchFactor**exp
            fillStart = start + bias
            if fillStart > maxValue:
                continue
            fillEnd = min(maxValue, start + nodeSize - 1 + bias)
            if fillStart < 0:
                fillStart = 0
            if fillStart <= fillEnd:
                result.update(range(fillStart, fillEnd + 1))
            continue

        # Non-zero node: each set bit identifies a child/value.
        exp = height - depth
        nextNodeSize = branchFactor**exp
        while True:
            bitIndex = _trailingZeros(bits, 32)
            if bitIndex == 32:
                break
            if depth == height:
                val = start + bitIndex + bias
                if val > maxValue:
                    queue.clear()
                    break
                if val >= 0:
                    result.add(val)
            else:
                startDelta = bitIndex * nextNodeSize
                queue.append((start + startDelta, depth + 1))
            bits &= ~(1 << bitIndex)

    return (result, bitStream.bytesConsumed())


def encode(values: Iterable[int]) -> bytes:
    """Encode a set of integers as a sparse bit set.

    Tries all branching factors and returns the shortest encoding.

    Args:
        values: iterable of non-negative integers.

    Returns:
        bytes containing the sparse bit set encoding.
    """
    valuesSorted = sorted(set(values))
    if not valuesSorted:
        return _encodeHeader(2, 0)

    maxValue = valuesSorted[-1]
    valueSet = set(valuesSorted)

    best: Optional[bytes] = None
    for branchFactor in (2, 4, 8, 32):
        height = _treeHeight(branchFactor, maxValue)
        if height > _BF_MAX_HEIGHT[branchFactor]:
            continue
        encoded = _encodeWithBf(valueSet, branchFactor, height)
        if best is None or len(encoded) < len(best):
            best = encoded

    if best is None:
        raise ValueError(f"Cannot encode max value {maxValue}")
    return best


def _encodeHeader(branchFactor: int, height: int) -> bytes:
    branchFactorToId = {2: 0, 4: 1, 8: 2, 32: 3}
    return bytes([(height << 2) | branchFactorToId[branchFactor]])


def _decodeHeader(headerByte: int) -> Tuple[int, int]:
    id = headerByte & 0x03
    idToBranchFactor = {0: 2, 1: 4, 2: 8, 3: 32}
    height = (headerByte >> 2) & 0x1F
    return idToBranchFactor[id], height


def _treeHeight(branchFactor: int, maxValue: int) -> int:
    """Return the minimum tree height needed to represent maxValue."""
    height = 1
    capacity = branchFactor
    while capacity <= maxValue:
        capacity *= branchFactor
        height += 1
    return height


def _encodeWithBf(valueSet: Set[int], branchFactor: int, height: int) -> bytes:
    if height == 0:
        return _encodeHeader(branchFactor, 0)

    # Build layers bottom-up: layer 0 = leaves (individual values),
    # each higher layer groups bf children into one parent bitmask.
    layers: list[Dict[int, int]] = [{}]  # list of dicts: nodeIndex -> bitmask

    for v in valueSet:
        nodeIndex = v // branchFactor
        bitPos = v % branchFactor
        if nodeIndex not in layers[0]:
            layers[0][nodeIndex] = 0
        layers[0][nodeIndex] |= 1 << bitPos

    for _ in range(1, height):
        prevLayer = layers[-1]
        newLayer: Dict[int, int] = {}
        for nodeIndex, bitmask in prevLayer.items():
            parentIndex = nodeIndex // branchFactor
            bitPos = nodeIndex % branchFactor
            if parentIndex not in newLayer:
                newLayer[parentIndex] = 0
            newLayer[parentIndex] |= 1 << bitPos
        layers.append(newLayer)

    # For zero-node optimization: track count of values in sorted list.
    valuesSorted = sorted(valueSet)

    def rangeCount(lo: int, hi: int) -> int:
        return bisect.bisect_right(valuesSorted, hi) - bisect.bisect_left(
            valuesSorted, lo
        )

    # Emit nodes BFS order (root to leaves).
    # Queue entries: (nodeIndex, depthFromRoot, rangeStart, rangeEnd)
    stream = _OutputBitStream(branchFactor)
    subtreeSize = branchFactor**height
    queue: deque[Tuple[int, int, int, int]] = deque([(0, 0, 0, subtreeSize - 1)])

    while queue:
        nodeIndex, depth, rangeStart, rangeEnd = queue.popleft()
        layerIdx = height - 1 - depth  # layers[0]=leaves, layers[height-1]=root

        bitmask = (
            layers[layerIdx].get(nodeIndex, 0) if 0 <= layerIdx < len(layers) else 0
        )

        # Zero-node optimization: if entire range is filled on an INTERNAL node,
        # write 0 and skip children.  At leaf level we always write the explicit
        # bitmask so the encoding matches the reference.
        if (
            depth < height - 1
            and rangeCount(rangeStart, rangeEnd) == rangeEnd - rangeStart + 1
        ):
            stream.write(0)
            continue

        stream.write(bitmask)

        if bitmask != 0 and depth < height - 1:
            childSize = (rangeEnd - rangeStart + 1) // branchFactor
            bits = bitmask
            while bits:
                bitIndex = _trailingZeros(bits, 32)
                childIndex = nodeIndex * branchFactor + bitIndex
                childStart = rangeStart + bitIndex * childSize
                childEnd = childStart + childSize - 1
                queue.append((childIndex, depth + 1, childStart, childEnd))
                bits &= ~(1 << bitIndex)

    return _encodeHeader(branchFactor, height) + stream.toBytes()


def _trailingZeros(val: int, maxBits: int) -> int:
    if val == 0:
        return maxBits
    count = 0
    while (val & 1) == 0:
        val >>= 1
        count += 1
    return count


class _InputBitStream:
    """Reads bit nodes from a byte array, starting after the header byte."""

    def __init__(self, data: bytes, branchFactor: int):
        self.data = data
        self.branchFactor = branchFactor
        self.byteIndex = 1  # skip the header byte
        self.subIndex = 0  # bit offset within current byte (for bf=2,4)

    def next(self) -> Optional[int]:
        if self.branchFactor in (2, 4):
            if self.byteIndex >= len(self.data):
                return None
            mask = (1 << self.branchFactor) - 1
            val = (self.data[self.byteIndex] >> self.subIndex) & mask
            self.subIndex += self.branchFactor
            if self.subIndex >= 8:
                self.subIndex = 0
                self.byteIndex += 1
            return val
        elif self.branchFactor == 8:
            if self.byteIndex >= len(self.data):
                return None
            val = self.data[self.byteIndex]
            self.byteIndex += 1
            return val
        elif self.branchFactor == 32:
            if self.byteIndex + 3 >= len(self.data):
                return None
            b = self.data
            i = self.byteIndex
            val = b[i] | (b[i + 1] << 8) | (b[i + 2] << 16) | (b[i + 3] << 24)
            self.byteIndex += 4
            return val
        return None

    def bytesConsumed(self) -> int:
        return self.byteIndex + (1 if self.subIndex > 0 else 0)


class _OutputBitStream:
    """Writes bit nodes into a byte array."""

    def __init__(self, branchFactor: int):
        self.branchFactor = branchFactor
        self.data = bytearray()
        self.subIndex = 0  # bit offset within current byte (for bf=2,4)

    def write(self, value: int) -> None:
        if self.branchFactor in (2, 4):
            mask = (1 << self.branchFactor) - 1
            value &= mask
            if self.subIndex == 0:
                self.data.append(0)
            self.data[-1] |= value << self.subIndex
            self.subIndex += self.branchFactor
            if self.subIndex >= 8:
                self.subIndex = 0
        elif self.branchFactor == 8:
            self.data.append(value & 0xFF)
        elif self.branchFactor == 32:
            self.data.append(value & 0xFF)
            self.data.append((value >> 8) & 0xFF)
            self.data.append((value >> 16) & 0xFF)
            self.data.append((value >> 24) & 0xFF)

    def toBytes(self) -> bytes:
        return bytes(self.data)
