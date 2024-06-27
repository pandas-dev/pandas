import abc
import operator
import re
from datetime import date, datetime
from itertools import repeat
from typing import Any, Dict, List, Optional, Union

from pyparsing import (
    CaselessKeyword,
    Forward,
    OpAssoc,
    ParserElement,
    ParseResults,
    QuotedString,
    Suppress,
    Word,
    alphanums,
    exceptions,
    infix_notation,
    one_of,
    pyparsing_common,
)

try:
    # TODO import directly when depending on pyparsing>=3.1.0
    from pyparsing import DelimitedList
except ImportError:
    # delimited_list is deprecated in favor of DelimitedList in pyparsing 3.1.0
    from pyparsing import delimited_list as DelimitedList  # type: ignore[assignment]

from .exceptions import (
    InvalidInputException,
    InvalidStateException,
)


def _cast(type_: str, value: Any) -> Union[date, datetime, float, int, str]:
    # values are always cast from string to target type
    value = str(value)

    if type_ in ("bigint", "int", "smallint", "tinyint"):
        try:
            return int(value)  # no size is enforced
        except ValueError:
            raise ValueError(f'"{value}" is not an integer.')

    if type_ == "decimal":
        try:
            return float(value)
        except ValueError:
            raise ValueError(f"{value} is not a decimal.")

    if type_ in ("char", "string", "varchar"):
        return value  # no length is enforced

    if type_ == "date":
        try:
            return datetime.strptime(value, "%Y-%m-%d").date()
        except ValueError:
            raise ValueError(f"{value} is not a date.")

    if type_ == "timestamp":
        match = re.search(
            r"^(?P<timestamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})"
            r"(?P<nanos>\.\d{1,9})?$",
            value,
        )
        if match is None:
            raise ValueError(
                "Timestamp format must be yyyy-mm-dd hh:mm:ss[.fffffffff]"
                f" {value} is not a timestamp."
            )

        try:
            timestamp = datetime.strptime(match.group("timestamp"), "%Y-%m-%d %H:%M:%S")
        except ValueError:
            raise ValueError(
                "Timestamp format must be yyyy-mm-dd hh:mm:ss[.fffffffff]"
                f" {value} is not a timestamp."
            )

        # use nanosecond representation for timestamps
        posix_nanoseconds = int(timestamp.timestamp() * 1_000_000_000)

        nanos = match.group("nanos")
        if nanos is not None:
            # strip leading dot, reverse and left pad with zeros to nanoseconds
            nanos = "".join(reversed(nanos[1:])).zfill(9)
            for i, nanoseconds in enumerate(nanos):
                posix_nanoseconds += int(nanoseconds) * 10**i

        return posix_nanoseconds

    raise InvalidInputException("GetPartitions", f"Unknown type : '{type_}'")


def _escape_regex(pattern: str) -> str:
    """Taken from Python 3.7 to avoid escaping '%'."""
    _special_chars_map = {i: "\\" + chr(i) for i in b"()[]{}?*+-|^$\\.&~# \t\n\r\v\f"}
    return pattern.translate(_special_chars_map)


class _Expr(abc.ABC):
    @abc.abstractmethod
    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> Any:  # type: ignore[misc]
        raise NotImplementedError()


class _Ident(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: str = tokens[0]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> Any:
        try:
            return self._eval(part_keys, part_input)
        except ValueError as e:
            # existing partition values cannot be cast to current schema
            raise InvalidStateException("GetPartitions", str(e))

    def leval(self, part_keys: List[Dict[str, str]], literal: Any) -> Any:
        # evaluate literal by simulating partition input
        try:
            return self._eval(part_keys, part_input={"Values": repeat(literal)})
        except ValueError as e:
            # expression literal cannot be cast to current schema
            raise InvalidInputException("GetPartitions", str(e))

    def type_(self, part_keys: List[Dict[str, str]]) -> str:
        for key in part_keys:
            if self.ident == key["Name"]:
                return key["Type"]

        raise InvalidInputException("GetPartitions", f"Unknown column '{self.ident}'")

    def _eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> Any:
        for key, value in zip(part_keys, part_input["Values"]):
            if self.ident == key["Name"]:
                return _cast(key["Type"], value)

        # also raised for unpartitioned tables
        raise InvalidInputException("GetPartitions", f"Unknown column '{self.ident}'")


class _IsNull(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: _Ident = tokens[0]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return self.ident.eval(part_keys, part_input) is None


class _IsNotNull(_IsNull):
    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return not super().eval(part_keys, part_input)


class _BinOp(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: _Ident = tokens[0]
        self.bin_op: str = tokens[1]
        self.literal: Any = tokens[2]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        ident = self.ident.eval(part_keys, part_input)

        # simulate partition input for the lateral
        rhs = self.ident.leval(part_keys, self.literal)

        return {
            "<>": operator.ne,
            ">=": operator.ge,
            "<=": operator.le,
            ">": operator.gt,
            "<": operator.lt,
            "=": operator.eq,
        }[self.bin_op](ident, rhs)


class _Like(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: _Ident = tokens[0]
        self.literal: str = tokens[2]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        type_ = self.ident.type_(part_keys)
        if type_ in ("bigint", "int", "smallint", "tinyint"):
            raise InvalidInputException(
                "GetPartitions", "Integral data type doesn't support operation 'LIKE'"
            )

        if type_ in ("date", "decimal", "timestamp"):
            raise InvalidInputException(
                "GetPartitions",
                f"{type_[0].upper()}{type_[1:]} data type"
                " doesn't support operation 'LIKE'",
            )

        ident = self.ident.eval(part_keys, part_input)
        assert isinstance(ident, str)

        pattern = _cast("string", self.literal)

        # prepare SQL pattern for conversion to regex pattern
        pattern = _escape_regex(pattern)  # type: ignore

        # NOTE convert SQL wildcards to regex, no literal matches possible
        pattern = pattern.replace("_", ".").replace("%", ".*")

        # LIKE clauses always start at the beginning
        pattern = "^" + pattern + "$"

        return re.search(pattern, ident) is not None


class _NotLike(_Like):
    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return not super().eval(part_keys, part_input)


class _In(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: _Ident = tokens[0]
        self.values: List[Any] = tokens[2:]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        ident = self.ident.eval(part_keys, part_input)
        values = (self.ident.leval(part_keys, value) for value in self.values)

        return ident in values


class _NotIn(_In):
    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return not super().eval(part_keys, part_input)


class _Between(_Expr):
    def __init__(self, tokens: ParseResults):
        self.ident: _Ident = tokens[0]
        self.left: Any = tokens[2]
        self.right: Any = tokens[4]

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        ident = self.ident.eval(part_keys, part_input)
        left = self.ident.leval(part_keys, self.left)
        right = self.ident.leval(part_keys, self.right)

        return left <= ident <= right or left > ident > right


class _NotBetween(_Between):
    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return not super().eval(part_keys, part_input)


class _BoolAnd(_Expr):
    def __init__(self, tokens: ParseResults) -> None:
        self.operands: List[_Expr] = tokens[0][0::2]  # skip 'and' between tokens

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return all(operand.eval(part_keys, part_input) for operand in self.operands)


class _BoolOr(_Expr):
    def __init__(self, tokens: ParseResults) -> None:
        self.operands: List[_Expr] = tokens[0][0::2]  # skip 'or' between tokens

    def eval(self, part_keys: List[Dict[str, str]], part_input: Dict[str, Any]) -> bool:
        return any(operand.eval(part_keys, part_input) for operand in self.operands)


class _PartitionFilterExpressionCache:
    def __init__(self) -> None:
        # build grammar according to Glue.Client.get_partitions(Expression)
        lpar, rpar = map(Suppress, "()")

        # NOTE these are AWS Athena column name best practices
        ident = Forward().set_name("ident")
        ident <<= Word(alphanums + "._").set_parse_action(_Ident) | lpar + ident + rpar  # type: ignore

        number = Forward().set_name("number")
        number <<= pyparsing_common.number | lpar + number + rpar  # type: ignore

        string = Forward().set_name("string")
        string <<= QuotedString(quote_char="'", esc_quote="''") | lpar + string + rpar  # type: ignore

        literal = (number | string).set_name("literal")
        literal_list = DelimitedList(literal, min=1).set_name("list")

        bin_op = one_of("<> >= <= > < =").set_name("binary op")

        and_ = Forward()
        and_ <<= CaselessKeyword("and") | lpar + and_ + rpar

        or_ = Forward()
        or_ <<= CaselessKeyword("or") | lpar + or_ + rpar

        in_, between, like, not_, is_, null = map(
            CaselessKeyword, "in between like not is null".split()
        )
        not_ = Suppress(not_)  # type: ignore  # only needed for matching

        cond = (
            (ident + is_ + null).set_parse_action(_IsNull)
            | (ident + is_ + not_ + null).set_parse_action(_IsNotNull)
            | (ident + bin_op + literal).set_parse_action(_BinOp)
            | (ident + like + string).set_parse_action(_Like)
            | (ident + not_ + like + string).set_parse_action(_NotLike)
            | (ident + in_ + lpar + literal_list + rpar).set_parse_action(_In)
            | (ident + not_ + in_ + lpar + literal_list + rpar).set_parse_action(_NotIn)
            | (ident + between + literal + and_ + literal).set_parse_action(_Between)
            | (ident + not_ + between + literal + and_ + literal).set_parse_action(
                _NotBetween
            )
        ).set_name("cond")

        # conditions can be joined using 2-ary AND and/or OR
        expr = infix_notation(
            cond,
            [
                (and_, 2, OpAssoc.LEFT, _BoolAnd),
                (or_, 2, OpAssoc.LEFT, _BoolOr),
            ],
        )
        self._expr = expr.set_name("expr")

        self._cache: Dict[str, _Expr] = {}

    def get(self, expression: Optional[str]) -> Optional[_Expr]:
        if expression is None:
            return None

        if expression not in self._cache:
            ParserElement.enable_packrat()

            try:
                expr: ParseResults = self._expr.parse_string(expression, parse_all=True)
                self._cache[expression] = expr[0]
            except exceptions.ParseException:
                raise InvalidInputException(
                    "GetPartitions", f"Unsupported expression '{expression}'"
                )

        return self._cache[expression]


_PARTITION_FILTER_EXPRESSION_CACHE = _PartitionFilterExpressionCache()


class PartitionFilter:
    def __init__(self, expression: Optional[str], fake_table: Any):
        self.expression = expression
        self.fake_table = fake_table

    def __call__(self, fake_partition: Any) -> bool:
        expression = _PARTITION_FILTER_EXPRESSION_CACHE.get(self.expression)
        if expression is None:
            return True

        versions = list(self.fake_table.versions.values())
        return expression.eval(
            part_keys=versions[-1].get("PartitionKeys", []),
            part_input=fake_partition.partition_input,
        )
