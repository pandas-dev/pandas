import io as StringIO
import re
import string
from typing import Dict, Iterable, List, Match, Optional, TextIO, Tuple

from .metrics_core import Metric
from .samples import Sample
from .validation import (
    _is_valid_legacy_metric_name, _validate_labelname, _validate_metric_name,
)


def text_string_to_metric_families(text: str) -> Iterable[Metric]:
    """Parse Prometheus text format from a unicode string.

    See text_fd_to_metric_families.
    """
    yield from text_fd_to_metric_families(StringIO.StringIO(text))


ESCAPE_SEQUENCES = {
    '\\\\': '\\',
    '\\n': '\n',
    '\\"': '"',
}


def replace_escape_sequence(match: Match[str]) -> str:
    return ESCAPE_SEQUENCES[match.group(0)]


HELP_ESCAPING_RE = re.compile(r'\\[\\n]')
ESCAPING_RE = re.compile(r'\\[\\n"]')


def _replace_help_escaping(s: str) -> str:
    return HELP_ESCAPING_RE.sub(replace_escape_sequence, s)


def _replace_escaping(s: str) -> str:
    return ESCAPING_RE.sub(replace_escape_sequence, s)


def _is_character_escaped(s: str, charpos: int) -> bool:
    num_bslashes = 0
    while (charpos > num_bslashes
           and s[charpos - 1 - num_bslashes] == '\\'):
        num_bslashes += 1
    return num_bslashes % 2 == 1


def parse_labels(labels_string: str, openmetrics: bool = False) -> Dict[str, str]:
    labels: Dict[str, str] = {}

    # Copy original labels
    sub_labels = labels_string.strip()
    if openmetrics and sub_labels and sub_labels[0] == ',':
        raise ValueError("leading comma: " + labels_string)
    try:
        # Process one label at a time
        while sub_labels:
            # The label name is before the equal, or if there's no equal, that's the
            # metric name.
            
            term, sub_labels = _next_term(sub_labels, openmetrics)
            if not term:
                if openmetrics:
                    raise ValueError("empty term in line: " + labels_string)
                continue
            
            quoted_name = False
            operator_pos = _next_unquoted_char(term, '=')
            if operator_pos == -1:
                quoted_name = True
                label_name = "__name__"
            else:
                value_start = _next_unquoted_char(term, '=')
                label_name, quoted_name = _unquote_unescape(term[:value_start])
                term = term[value_start + 1:]
                
            if not quoted_name and not _is_valid_legacy_metric_name(label_name):
                raise ValueError("unquoted UTF-8 metric name")
                
            # Check for missing quotes 
            term = term.strip()
            if not term or term[0] != '"':
                raise ValueError

            # The first quote is guaranteed to be after the equal.
            # Find the last unescaped quote.
            i = 1
            while i < len(term):
                i = term.index('"', i)
                if not _is_character_escaped(term[:i], i):
                    break
                i += 1

            # The label value is between the first and last quote
            quote_end = i + 1
            if quote_end != len(term):
                raise ValueError("unexpected text after quote: " + labels_string)
            label_value, _ = _unquote_unescape(term[:quote_end])
            if label_name == '__name__':
                _validate_metric_name(label_name)
            else:
                _validate_labelname(label_name)
            if label_name in labels:
                raise ValueError("invalid line, duplicate label name: " + labels_string)
            labels[label_name] = label_value
        return labels
    except ValueError:
        raise ValueError("Invalid labels: " + labels_string)
    

def _next_term(text: str, openmetrics: bool) -> Tuple[str, str]:
    """Extract the next comma-separated label term from the text.
    
    Returns the stripped term and the stripped remainder of the string, 
    including the comma.
    
    Raises ValueError if the term is empty and we're in openmetrics mode.
    """
    
    # There may be a leading comma, which is fine here.
    if text[0] == ',':
        text = text[1:]
        if not text:
            return "", ""
        if text[0] == ',':
            raise ValueError("multiple commas")
    splitpos = _next_unquoted_char(text, ',}')
    if splitpos == -1:
        splitpos = len(text)
    term = text[:splitpos]
    if not term and openmetrics:
        raise ValueError("empty term:", term)
    
    sublabels = text[splitpos:]
    return term.strip(), sublabels.strip()


def _next_unquoted_char(text: str, chs: str, startidx: int = 0) -> int:
    """Return position of next unquoted character in tuple, or -1 if not found.
    
    It is always assumed that the first character being checked is not already
    inside quotes.
    """
    i = startidx
    in_quotes = False
    if chs is None:
        chs = string.whitespace
    while i < len(text):
        if text[i] == '"' and not _is_character_escaped(text, i):
            in_quotes = not in_quotes
        if not in_quotes:
            if text[i] in chs:
                return i
        i += 1
    return -1


def _last_unquoted_char(text: str, chs: str) -> int:
    """Return position of last unquoted character in list, or -1 if not found."""
    i = len(text) - 1
    in_quotes = False
    if chs is None:
        chs = string.whitespace
    while i > 0:
        if text[i] == '"' and not _is_character_escaped(text, i):
            in_quotes = not in_quotes
            
        if not in_quotes:
            if text[i] in chs:
                return i
        i -= 1
    return -1


def _split_quoted(text, separator, maxsplit=0):
    """Splits on split_ch similarly to strings.split, skipping separators if
    they are inside quotes.
    """

    tokens = ['']
    x = 0
    while x < len(text):
        split_pos = _next_unquoted_char(text, separator, x)
        if split_pos == -1:
            tokens[-1] = text[x:]
            x = len(text)
            continue
        if maxsplit > 0 and len(tokens) > maxsplit:
            tokens[-1] = text[x:]
            break
        tokens[-1] = text[x:split_pos]
        x = split_pos + 1
        tokens.append('')
    return tokens


def _unquote_unescape(text):
    """Returns the string, and true if it was quoted."""
    if not text:
        return text, False
    quoted = False
    text = text.strip()
    if text[0] == '"':
        if len(text) == 1 or text[-1] != '"':
            raise ValueError("missing close quote")
        text = text[1:-1]
        quoted = True
    if "\\" in text:
        text = _replace_escaping(text)
    return text, quoted


# If we have multiple values only consider the first
def _parse_value_and_timestamp(s: str) -> Tuple[float, Optional[float]]:
    s = s.lstrip()
    separator = " "
    if separator not in s:
        separator = "\t"
    values = [value.strip() for value in s.split(separator) if value.strip()]
    if not values:
        return float(s), None
    value = _parse_value(values[0])
    timestamp = (_parse_value(values[-1]) / 1000) if len(values) > 1 else None
    return value, timestamp


def _parse_value(value):
    value = ''.join(value)
    if value != value.strip() or '_' in value:
        raise ValueError(f"Invalid value: {value!r}")
    try:
        return int(value)
    except ValueError:
        return float(value)
    

def _parse_sample(text):
    separator = " # "
    # Detect the labels in the text
    label_start = _next_unquoted_char(text, '{')
    if label_start == -1 or separator in text[:label_start]:
        # We don't have labels, but there could be an exemplar.
        name_end = _next_unquoted_char(text, ' \t')
        name = text[:name_end].strip()
        if not _is_valid_legacy_metric_name(name):
            raise ValueError("invalid metric name:" + text)
        # Parse the remaining text after the name
        remaining_text = text[name_end + 1:]
        value, timestamp = _parse_value_and_timestamp(remaining_text)
        return Sample(name, {}, value, timestamp)
    name = text[:label_start].strip()
    label_end = _next_unquoted_char(text, '}')
    labels = parse_labels(text[label_start + 1:label_end], False)
    if not name:
        # Name might be in the labels
        if '__name__' not in labels:
            raise ValueError
        name = labels['__name__']
        del labels['__name__']
    elif '__name__' in labels:
        raise ValueError("metric name specified more than once")
    # Parsing labels succeeded, continue parsing the remaining text
    remaining_text = text[label_end + 1:]
    value, timestamp = _parse_value_and_timestamp(remaining_text)
    return Sample(name, labels, value, timestamp)


def text_fd_to_metric_families(fd: TextIO) -> Iterable[Metric]:
    """Parse Prometheus text format from a file descriptor.

    This is a laxer parser than the main Go parser,
    so successful parsing does not imply that the parsed
    text meets the specification.

    Yields Metric's.
    """
    name = ''
    documentation = ''
    typ = 'untyped'
    samples: List[Sample] = []
    allowed_names = []

    def build_metric(name: str, documentation: str, typ: str, samples: List[Sample]) -> Metric:
        # Munge counters into OpenMetrics representation
        # used internally.
        if typ == 'counter':
            if name.endswith('_total'):
                name = name[:-6]
            else:
                new_samples = []
                for s in samples:
                    new_samples.append(Sample(s[0] + '_total', *s[1:]))
                    samples = new_samples
        metric = Metric(name, documentation, typ)
        metric.samples = samples
        return metric

    for line in fd:
        line = line.strip()

        if line.startswith('#'):
            parts = _split_quoted(line, None, 3)
            if len(parts) < 2:
                continue
            candidate_name, quoted = '', False
            if len(parts) > 2:
                # Ignore comment tokens
                if parts[1] != 'TYPE' and parts[1] != 'HELP':
                    continue
                candidate_name, quoted = _unquote_unescape(parts[2])
                if not quoted and not _is_valid_legacy_metric_name(candidate_name):
                    raise ValueError
            if parts[1] == 'HELP':
                if candidate_name != name:
                    if name != '':
                        yield build_metric(name, documentation, typ, samples)
                    # New metric
                    name = candidate_name
                    typ = 'untyped'
                    samples = []
                    allowed_names = [candidate_name]
                if len(parts) == 4:
                    documentation = _replace_help_escaping(parts[3])
                else:
                    documentation = ''
            elif parts[1] == 'TYPE':
                if len(parts) < 4:
                    raise ValueError
                if candidate_name != name:
                    if name != '':
                        yield build_metric(name, documentation, typ, samples)
                    # New metric
                    name = candidate_name
                    documentation = ''
                    samples = []
                typ = parts[3]
                allowed_names = {
                    'counter': [''],
                    'gauge': [''],
                    'summary': ['_count', '_sum', ''],
                    'histogram': ['_count', '_sum', '_bucket'],
                }.get(typ, [''])
                allowed_names = [name + n for n in allowed_names]
        elif line == '':
            # Ignore blank lines
            pass
        else:
            sample = _parse_sample(line)
            if sample.name not in allowed_names:
                if name != '':
                    yield build_metric(name, documentation, typ, samples)
                # New metric, yield immediately as untyped singleton
                name = ''
                documentation = ''
                typ = 'untyped'
                samples = []
                allowed_names = []
                yield build_metric(sample[0], documentation, typ, [sample])
            else:
                samples.append(sample)

    if name != '':
        yield build_metric(name, documentation, typ, samples)
