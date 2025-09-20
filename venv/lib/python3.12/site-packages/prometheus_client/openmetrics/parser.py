#!/usr/bin/env python


import io as StringIO
import math
import re

from ..metrics_core import Metric
from ..parser import (
    _last_unquoted_char, _next_unquoted_char, _parse_value, _split_quoted,
    _unquote_unescape, parse_labels,
)
from ..samples import BucketSpan, Exemplar, NativeHistogram, Sample, Timestamp
from ..utils import floatToGoString
from ..validation import _is_valid_legacy_metric_name, _validate_metric_name


def text_string_to_metric_families(text):
    """Parse Openmetrics text format from a unicode string.

    See text_fd_to_metric_families.
    """
    yield from text_fd_to_metric_families(StringIO.StringIO(text))


_CANONICAL_NUMBERS = {float("inf")}


def _isUncanonicalNumber(s):
    f = float(s)
    if f not in _CANONICAL_NUMBERS:
        return False  # Only the canonical numbers are required to be canonical.
    return s != floatToGoString(f)


ESCAPE_SEQUENCES = {
    '\\\\': '\\',
    '\\n': '\n',
    '\\"': '"',
}


def _replace_escape_sequence(match):
    return ESCAPE_SEQUENCES[match.group(0)]


ESCAPING_RE = re.compile(r'\\[\\n"]')


def _replace_escaping(s):
    return ESCAPING_RE.sub(_replace_escape_sequence, s)


def _unescape_help(text):
    result = []
    slash = False

    for char in text:
        if slash:
            if char == '\\':
                result.append('\\')
            elif char == '"':
                result.append('"')
            elif char == 'n':
                result.append('\n')
            else:
                result.append('\\' + char)
            slash = False
        else:
            if char == '\\':
                slash = True
            else:
                result.append(char)

    if slash:
        result.append('\\')

    return ''.join(result)


def _parse_timestamp(timestamp):
    timestamp = ''.join(timestamp)
    if not timestamp:
        return None
    if timestamp != timestamp.strip() or '_' in timestamp:
        raise ValueError(f"Invalid timestamp: {timestamp!r}")
    try:
        # Simple int.
        return Timestamp(int(timestamp), 0)
    except ValueError:
        try:
            # aaaa.bbbb. Nanosecond resolution supported.
            parts = timestamp.split('.', 1)
            return Timestamp(int(parts[0]), int(parts[1][:9].ljust(9, "0")))
        except ValueError:
            # Float.
            ts = float(timestamp)
            if math.isnan(ts) or math.isinf(ts):
                raise ValueError(f"Invalid timestamp: {timestamp!r}")
            return ts


def _is_character_escaped(s, charpos):
    num_bslashes = 0
    while (charpos > num_bslashes
           and s[charpos - 1 - num_bslashes] == '\\'):
        num_bslashes += 1
    return num_bslashes % 2 == 1


def _parse_sample(text):
    separator = " # "
    # Detect the labels in the text
    label_start = _next_unquoted_char(text, '{')
    if label_start == -1 or separator in text[:label_start]:
        # We don't have labels, but there could be an exemplar.
        name_end = _next_unquoted_char(text, ' ')
        name = text[:name_end]
        if not _is_valid_legacy_metric_name(name):
            raise ValueError("invalid metric name:" + text)
        # Parse the remaining text after the name
        remaining_text = text[name_end + 1:]
        value, timestamp, exemplar = _parse_remaining_text(remaining_text)
        return Sample(name, {}, value, timestamp, exemplar)
    name = text[:label_start]
    label_end = _next_unquoted_char(text, '}')
    labels = parse_labels(text[label_start + 1:label_end], True)
    if not name:
        # Name might be in the labels
        if '__name__' not in labels:
            raise ValueError
        name = labels['__name__']
        del labels['__name__']
    elif '__name__' in labels:
        raise ValueError("metric name specified more than once")
    # Parsing labels succeeded, continue parsing the remaining text
    remaining_text = text[label_end + 2:]
    value, timestamp, exemplar = _parse_remaining_text(remaining_text)
    return Sample(name, labels, value, timestamp, exemplar)


def _parse_remaining_text(text):
    split_text = text.split(" ", 1)
    val = _parse_value(split_text[0])
    if len(split_text) == 1:
        # We don't have timestamp or exemplar
        return val, None, None  

    timestamp = []
    exemplar_value = []
    exemplar_timestamp = []
    exemplar_labels = None

    state = 'timestamp'
    text = split_text[1]

    it = iter(text)
    in_quotes = False
    for char in it:
        if char == '"':
            in_quotes = not in_quotes
        if in_quotes:
            continue
        if state == 'timestamp':
            if char == '#' and not timestamp:
                state = 'exemplarspace'
            elif char == ' ':
                state = 'exemplarhash'
            else:
                timestamp.append(char)
        elif state == 'exemplarhash':
            if char == '#':
                state = 'exemplarspace'
            else:
                raise ValueError("Invalid line: " + text)
        elif state == 'exemplarspace':
            if char == ' ':
                state = 'exemplarstartoflabels'
            else:
                raise ValueError("Invalid line: " + text)
        elif state == 'exemplarstartoflabels':
            if char == '{':
                label_start = _next_unquoted_char(text, '{')
                label_end = _last_unquoted_char(text, '}')
                exemplar_labels = parse_labels(text[label_start + 1:label_end], True)
                state = 'exemplarparsedlabels'
            else:
                raise ValueError("Invalid line: " + text)
        elif state == 'exemplarparsedlabels':
            if char == '}':
                state = 'exemplarvaluespace'
        elif state == 'exemplarvaluespace':
            if char == ' ':
                state = 'exemplarvalue'
            else:
                raise ValueError("Invalid line: " + text)
        elif state == 'exemplarvalue':
            if char == ' ' and not exemplar_value:
                raise ValueError("Invalid line: " + text)
            elif char == ' ':
                state = 'exemplartimestamp'
            else:
                exemplar_value.append(char)
        elif state == 'exemplartimestamp':
            exemplar_timestamp.append(char)

    # Trailing space after value.
    if state == 'timestamp' and not timestamp:
        raise ValueError("Invalid line: " + text)

    # Trailing space after value.
    if state == 'exemplartimestamp' and not exemplar_timestamp:
        raise ValueError("Invalid line: " + text)

    # Incomplete exemplar.
    if state in ['exemplarhash', 'exemplarspace', 'exemplarstartoflabels', 'exemplarparsedlabels']:
        raise ValueError("Invalid line: " + text)

    ts = _parse_timestamp(timestamp)
    exemplar = None
    if exemplar_labels is not None:
        exemplar_length = sum(len(k) + len(v) for k, v in exemplar_labels.items())
        if exemplar_length > 128:
            raise ValueError("Exemplar labels are too long: " + text)
        exemplar = Exemplar(
            exemplar_labels,
            _parse_value(exemplar_value),
            _parse_timestamp(exemplar_timestamp),
        )

    return val, ts, exemplar


def _parse_nh_sample(text, suffixes):
    """Determines if the line has a native histogram sample, and parses it if so."""
    labels_start = _next_unquoted_char(text, '{')
    labels_end = -1

    # Finding a native histogram sample requires careful parsing of
    # possibly-quoted text, which can appear in metric names, label names, and
    # values.
    # 
    # First, we need to determine if there are metric labels. Find the space
    # between the metric definition and the rest of the line. Look for unquoted
    # space or {.
    i = 0
    has_metric_labels = False
    i = _next_unquoted_char(text, ' {')
    if i == -1:
        return

    # If the first unquoted char was a {, then that is the metric labels (which
    # could contain a UTF-8 metric name).
    if text[i] == '{':
        has_metric_labels = True
        # Consume the labels -- jump ahead to the close bracket.
        labels_end = i = _next_unquoted_char(text, '}', i)
        if labels_end == -1:
            raise ValueError
    
    # If there is no subsequent unquoted {, then it's definitely not a nh.
    nh_value_start = _next_unquoted_char(text, '{', i + 1)
    if nh_value_start == -1:
        return
    
    # Edge case: if there is an unquoted # between the metric definition and the {,
    # then this is actually an exemplar
    exemplar = _next_unquoted_char(text, '#', i + 1)
    if exemplar != -1 and exemplar < nh_value_start:
        return
    
    nh_value_end = _next_unquoted_char(text, '}', nh_value_start)
    if nh_value_end == -1:
        raise ValueError
    
    if has_metric_labels:
        labelstext = text[labels_start + 1:labels_end]
        labels = parse_labels(labelstext, True)
        name_end = labels_start
        name = text[:name_end]
        if name.endswith(suffixes):
            raise ValueError("the sample name of a native histogram with labels should have no suffixes", name)
        if not name:
            # Name might be in the labels
            if '__name__' not in labels:
                raise ValueError
            name = labels['__name__']
            del labels['__name__']
            # Edge case: the only "label" is the name definition.
            if not labels:
                labels = None
             
        nh_value = text[nh_value_start:]
        nat_hist_value = _parse_nh_struct(nh_value)
        return Sample(name, labels, None, None, None, nat_hist_value)
    # check if it's a native histogram
    else:
        nh_value = text[nh_value_start:]
        name_end = nh_value_start - 1
        name = text[:name_end]
        if name.endswith(suffixes):
            raise ValueError("the sample name of a native histogram should have no suffixes", name)
        # Not possible for UTF-8 name here, that would have been caught as having a labelset.
        nat_hist_value = _parse_nh_struct(nh_value)
        return Sample(name, None, None, None, None, nat_hist_value)      


def _parse_nh_struct(text):
    pattern = r'(\w+):\s*([^,}]+)'
    re_spans = re.compile(r'(positive_spans|negative_spans):\[(\d+:\d+(,\d+:\d+)*)\]')
    re_deltas = re.compile(r'(positive_deltas|negative_deltas):\[(-?\d+(?:,-?\d+)*)\]')

    items = dict(re.findall(pattern, text))
    span_matches = re_spans.findall(text)
    deltas = dict(re_deltas.findall(text))

    count_value = int(items['count'])
    sum_value = int(items['sum'])
    schema = int(items['schema'])
    zero_threshold = float(items['zero_threshold'])
    zero_count = int(items['zero_count'])

    pos_spans = _compose_spans(span_matches, 'positive_spans')
    neg_spans = _compose_spans(span_matches, 'negative_spans')
    pos_deltas = _compose_deltas(deltas, 'positive_deltas')
    neg_deltas = _compose_deltas(deltas, 'negative_deltas')
      
    return NativeHistogram(
        count_value=count_value,
        sum_value=sum_value,
        schema=schema,
        zero_threshold=zero_threshold,
        zero_count=zero_count,
        pos_spans=pos_spans,
        neg_spans=neg_spans,
        pos_deltas=pos_deltas,
        neg_deltas=neg_deltas
    )
  

def _compose_spans(span_matches, spans_name):
    """Takes a list of span matches (expected to be a list of tuples) and a string 
    (the expected span list name) and processes the list so that the values extracted 
    from the span matches can be used to compose a tuple of BucketSpan objects"""
    spans = {}
    for match in span_matches:
        # Extract the key from the match (first element of the tuple).
        key = match[0]
        # Extract the value from the match (second element of the tuple).
        # Split the value string by commas to get individual pairs, 
        # split each pair by ':' to get start and end, and convert them to integers.
        value = [tuple(map(int, pair.split(':'))) for pair in match[1].split(',')]
        # Store the processed value in the spans dictionary with the key.
        spans[key] = value
    if spans_name not in spans:
        return None
    out_spans = []
    # Iterate over each start and end tuple in the list of tuples for the specified spans_name.
    for start, end in spans[spans_name]:
        # Compose a BucketSpan object with the start and end values 
        # and append it to the out_spans list.
        out_spans.append(BucketSpan(start, end))
        # Convert to tuple
    out_spans_tuple = tuple(out_spans)
    return out_spans_tuple


def _compose_deltas(deltas, deltas_name):
    """Takes a list of deltas matches (a dictionary) and a string (the expected delta list name),
    and processes its elements to compose a tuple of integers representing the deltas"""
    if deltas_name not in deltas:
        return None
    out_deltas = deltas.get(deltas_name)
    if out_deltas is not None and out_deltas.strip():
        elems = out_deltas.split(',')
    # Convert each element in the list elems to an integer 
    # after stripping whitespace and create a tuple from these integers.
    out_deltas_tuple = tuple(int(x.strip()) for x in elems)
    return out_deltas_tuple
        

def _group_for_sample(sample, name, typ):
    if typ == 'info':
        # We can't distinguish between groups for info metrics.
        return {}
    if typ == 'summary' and sample.name == name:
        d = sample.labels.copy()
        del d['quantile']
        return d
    if typ == 'stateset':
        d = sample.labels.copy()
        del d[name]
        return d
    if typ in ['histogram', 'gaugehistogram'] and sample.name == name + '_bucket':
        d = sample.labels.copy()
        del d['le']
        return d
    return sample.labels


def _check_histogram(samples, name):
    group = None
    timestamp = None

    def do_checks():
        if bucket != float('+Inf'):
            raise ValueError("+Inf bucket missing: " + name)
        if count is not None and value != count:
            raise ValueError("Count does not match +Inf value: " + name)
        if has_sum and count is None:
            raise ValueError("_count must be present if _sum is present: " + name)
        if has_gsum and count is None:
            raise ValueError("_gcount must be present if _gsum is present: " + name)
        if not (has_sum or has_gsum) and count is not None:
            raise ValueError("_sum/_gsum must be present if _count is present: " + name)
        if has_negative_buckets and has_sum:
            raise ValueError("Cannot have _sum with negative buckets: " + name)
        if not has_negative_buckets and has_negative_gsum:
            raise ValueError("Cannot have negative _gsum with non-negative buckets: " + name)

    for s in samples:
        suffix = s.name[len(name):]
        g = _group_for_sample(s, name, 'histogram')
        if len(suffix) == 0:
            continue
        if g != group or s.timestamp != timestamp:
            if group is not None:
                do_checks()
            count = None
            bucket = None
            has_negative_buckets = False
            has_sum = False
            has_gsum = False
            has_negative_gsum = False
            value = 0
        group = g
        timestamp = s.timestamp

        if suffix == '_bucket':
            b = float(s.labels['le'])
            if b < 0:
                has_negative_buckets = True
            if bucket is not None and b <= bucket:
                raise ValueError("Buckets out of order: " + name)
            if s.value < value:
                raise ValueError("Bucket values out of order: " + name)
            bucket = b
            value = s.value
        elif suffix in ['_count', '_gcount']:
            count = s.value
        elif suffix in ['_sum']:
            has_sum = True
        elif suffix in ['_gsum']:
            has_gsum = True
            if s.value < 0:
                has_negative_gsum = True

    if group is not None:
        do_checks()


def text_fd_to_metric_families(fd):
    """Parse Prometheus text format from a file descriptor.

    This is a laxer parser than the main Go parser,
    so successful parsing does not imply that the parsed
    text meets the specification.

    Yields Metric's.
    """
    name = None
    allowed_names = []
    eof = False

    seen_names = set()
    type_suffixes = {
        'counter': ['_total', '_created'],
        'summary': ['', '_count', '_sum', '_created'],
        'histogram': ['_count', '_sum', '_bucket', '_created'],
        'gaugehistogram': ['_gcount', '_gsum', '_bucket'],
        'info': ['_info'],
    }

    def build_metric(name, documentation, typ, unit, samples):
        if typ is None:
            typ = 'unknown'
        for suffix in set(type_suffixes.get(typ, []) + [""]):
            if name + suffix in seen_names:
                raise ValueError("Clashing name: " + name + suffix)
            seen_names.add(name + suffix)
        if documentation is None:
            documentation = ''
        if unit is None:
            unit = ''
        if unit and not name.endswith("_" + unit):
            raise ValueError("Unit does not match metric name: " + name)
        if unit and typ in ['info', 'stateset']:
            raise ValueError("Units not allowed for this metric type: " + name)
        if typ in ['histogram', 'gaugehistogram']:
            _check_histogram(samples, name)
        _validate_metric_name(name)
        metric = Metric(name, documentation, typ, unit)
        # TODO: check labelvalues are valid utf8
        metric.samples = samples
        return metric

    is_nh = False
    typ = None
    for line in fd:
        if line[-1] == '\n':
            line = line[:-1]

        if eof:
            raise ValueError("Received line after # EOF: " + line)

        if not line:
            raise ValueError("Received blank line")

        if line == '# EOF':
            eof = True
        elif line.startswith('#'):
            parts = _split_quoted(line, ' ', 3)
            if len(parts) < 4:
                raise ValueError("Invalid line: " + line)
            candidate_name, quoted = _unquote_unescape(parts[2])
            if not quoted and not _is_valid_legacy_metric_name(candidate_name):
                raise ValueError
            if candidate_name == name and samples:
                raise ValueError("Received metadata after samples: " + line)
            if candidate_name != name:
                if name is not None:
                    yield build_metric(name, documentation, typ, unit, samples)
                # New metric
                name = candidate_name
                unit = None
                typ = None
                documentation = None
                group = None
                seen_groups = set()
                group_timestamp = None
                group_timestamp_samples = set()
                samples = []
                allowed_names = [candidate_name]
            
            if parts[1] == 'HELP':
                if documentation is not None:
                    raise ValueError("More than one HELP for metric: " + line)
                documentation = _unescape_help(parts[3])
            elif parts[1] == 'TYPE':
                if typ is not None:
                    raise ValueError("More than one TYPE for metric: " + line)
                typ = parts[3]
                if typ == 'untyped':
                    raise ValueError("Invalid TYPE for metric: " + line)
                allowed_names = [name + n for n in type_suffixes.get(typ, [''])]
            elif parts[1] == 'UNIT':
                if unit is not None:
                    raise ValueError("More than one UNIT for metric: " + line)
                unit = parts[3]
            else:
                raise ValueError("Invalid line: " + line)
        else:
            if typ == 'histogram':
                # set to true to account for native histograms naming exceptions/sanitizing differences
                is_nh = True
                sample = _parse_nh_sample(line, tuple(type_suffixes['histogram']))
                # It's not a native histogram
                if sample is None:
                    is_nh = False
                    sample = _parse_sample(line)              
            else:
                is_nh = False
                sample = _parse_sample(line)
            if sample.name not in allowed_names and not is_nh:
                if name is not None:
                    yield build_metric(name, documentation, typ, unit, samples)
                # Start an unknown metric.
                candidate_name, quoted = _unquote_unescape(sample.name)
                if not quoted and not _is_valid_legacy_metric_name(candidate_name):
                    raise ValueError
                name = candidate_name
                documentation = None
                unit = None
                typ = 'unknown'
                samples = []
                group = None
                group_timestamp = None
                group_timestamp_samples = set()
                seen_groups = set()
                allowed_names = [sample.name]

            if typ == 'stateset' and name not in sample.labels:
                raise ValueError("Stateset missing label: " + line)
            if (name + '_bucket' == sample.name
                    and (sample.labels.get('le', "NaN") == "NaN"
                         or _isUncanonicalNumber(sample.labels['le']))):
                raise ValueError("Invalid le label: " + line)
            if (name + '_bucket' == sample.name
                    and (not isinstance(sample.value, int) and not sample.value.is_integer())):
                raise ValueError("Bucket value must be an integer: " + line)
            if ((name + '_count' == sample.name or name + '_gcount' == sample.name)
                    and (not isinstance(sample.value, int) and not sample.value.is_integer())):
                raise ValueError("Count value must be an integer: " + line)
            if (typ == 'summary' and name == sample.name
                    and (not (0 <= float(sample.labels.get('quantile', -1)) <= 1)
                         or _isUncanonicalNumber(sample.labels['quantile']))):
                raise ValueError("Invalid quantile label: " + line)

            if not is_nh:
                g = tuple(sorted(_group_for_sample(sample, name, typ).items()))
                if group is not None and g != group and g in seen_groups:
                    raise ValueError("Invalid metric grouping: " + line)
                if group is not None and g == group:
                    if (sample.timestamp is None) != (group_timestamp is None):
                        raise ValueError("Mix of timestamp presence within a group: " + line)
                    if group_timestamp is not None and group_timestamp > sample.timestamp and typ != 'info':
                        raise ValueError("Timestamps went backwards within a group: " + line)
                else:
                    group_timestamp_samples = set()

                series_id = (sample.name, tuple(sorted(sample.labels.items())))
                if sample.timestamp != group_timestamp or series_id not in group_timestamp_samples:
                    # Not a duplicate due to timestamp truncation.
                    samples.append(sample)
                group_timestamp_samples.add(series_id)

                group = g
                group_timestamp = sample.timestamp
                seen_groups.add(g)
            else:
                samples.append(sample)

            if typ == 'stateset' and sample.value not in [0, 1]:
                raise ValueError("Stateset samples can only have values zero and one: " + line)
            if typ == 'info' and sample.value != 1:
                raise ValueError("Info samples can only have value one: " + line)
            if typ == 'summary' and name == sample.name and sample.value < 0:
                raise ValueError("Quantile values cannot be negative: " + line)
            if sample.name[len(name):] in ['_total', '_sum', '_count', '_bucket', '_gcount', '_gsum'] and math.isnan(
                    sample.value):
                raise ValueError("Counter-like samples cannot be NaN: " + line)
            if sample.name[len(name):] in ['_total', '_sum', '_count', '_bucket', '_gcount'] and sample.value < 0:
                raise ValueError("Counter-like samples cannot be negative: " + line)
            if sample.exemplar and not (
                    (typ in ['histogram', 'gaugehistogram'] and sample.name.endswith('_bucket'))
                    or (typ in ['counter'] and sample.name.endswith('_total'))):
                raise ValueError("Invalid line only histogram/gaugehistogram buckets and counters can have exemplars: " + line)
    
    if name is not None:
        yield build_metric(name, documentation, typ, unit, samples)

    if not eof:
        raise ValueError("Missing # EOF at end")
