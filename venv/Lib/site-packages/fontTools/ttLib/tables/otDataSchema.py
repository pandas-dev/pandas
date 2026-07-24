"""fontTools.ttLib.tables.otDataSchema -- contains the FieldSpec dataclass, which is
used to define the schema of OpenType tables in otData.py.
"""

from typing import Any, NamedTuple


class FieldSpec(NamedTuple):
    """Represents a single field in an OpenType table schema.

    These are defined in `otData.py` and used to construct converters in
    `otConverters.py`.

    Attributes:
        type: The type of the field. This string is used by `otConverters.py` to
            select or construct the appropriate converter. It can be a simple
            type name (e.g. 'uint16', 'Tag', 'Offset'), a specialized type (e.g.
            'ValueRecord', 'DeltaValue'), or a more complex expression (e.g.
            'OffsetTo(AatStateTable)', 'AATLookup(UShort)') which is evaluated
            to create a converter factory.
        name: The name of the field as defined in the OpenType specification.
        repeat: If the field is an array, this defines how many times it
            repeats. It can be an integer, or the name of another field (usually
            a 'Count' field) that contains the number of repetitions. If it is a
            string, it is used to lookup the count value from the table.
        aux: Auxiliary data used by some converters. When `repeat` is not set
            and `aux` is a string, it is often compiled into a Python
            expression for evaluation. For repeating structs, it is often used
            as a constant or an expression.
        description: A human-readable description of the field's purpose.
    """

    type: str
    name: str
    repeat: Any = None
    aux: Any = None
    description: str = ""
