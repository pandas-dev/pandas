#!/usr/bin/env python

from io import StringIO
from sys import maxunicode
from typing import Callable

from ..utils import floatToGoString, parse_version
from ..validation import (
    _is_valid_legacy_labelname, _is_valid_legacy_metric_name,
)

CONTENT_TYPE_LATEST = 'application/openmetrics-text; version=1.0.0; charset=utf-8'
"""Content type of the latest OpenMetrics 1.0 text format"""
CONTENT_TYPE_LATEST_2_0 = 'application/openmetrics-text; version=2.0.0; charset=utf-8'
"""Content type of the OpenMetrics 2.0 text format"""
ESCAPING_HEADER_TAG = 'escaping'


ALLOWUTF8 = 'allow-utf-8'
UNDERSCORES = 'underscores'
DOTS = 'dots'
VALUES = 'values'


def _is_valid_exemplar_metric(metric, sample):
    if metric.type == 'counter' and sample.name.endswith('_total'):
        return True
    if metric.type in ('gaugehistogram') and sample.name.endswith('_bucket'):
        return True
    if metric.type in ('histogram') and sample.name.endswith('_bucket') or sample.name == metric.name:
        return True
    return False


def _compose_exemplar_string(metric, sample, exemplar):
    """Constructs an exemplar string."""
    if not _is_valid_exemplar_metric(metric, sample):
        raise ValueError(f"Metric {metric.name} has exemplars, but is not a histogram bucket or counter")
    labels = '{{{0}}}'.format(','.join(
        ['{}="{}"'.format(
            k, v.replace('\\', r'\\').replace('\n', r'\n').replace('"', r'\"'))
            for k, v in sorted(exemplar.labels.items())]))
    if exemplar.timestamp is not None:
        exemplarstr = ' # {} {} {}'.format(
            labels,
            floatToGoString(exemplar.value),
            exemplar.timestamp,
        )
    else:
        exemplarstr = ' # {} {}'.format(
            labels,
            floatToGoString(exemplar.value),
        )

    return exemplarstr


def generate_latest(registry, escaping=UNDERSCORES, version="1.0.0"):
    '''Returns the metrics from the registry in latest text format as a string.'''
    output = []
    for metric in registry.collect():
        try:
            mname = metric.name
            output.append('# HELP {} {}\n'.format(
                escape_metric_name(mname, escaping), _escape(metric.documentation, ALLOWUTF8, _is_legacy_labelname_rune)))
            output.append(f'# TYPE {escape_metric_name(mname, escaping)} {metric.type}\n')
            if metric.unit:
                output.append(f'# UNIT {escape_metric_name(mname, escaping)} {metric.unit}\n')
            for s in metric.samples:
                if escaping == ALLOWUTF8 and not _is_valid_legacy_metric_name(s.name):
                    labelstr = escape_metric_name(s.name, escaping)
                    if s.labels:
                        labelstr += ','
                else:
                    labelstr = ''

                if s.labels:
                    items = sorted(s.labels.items())
                    # Label values always support UTF-8
                    labelstr += ','.join(
                        ['{}="{}"'.format(
                            escape_label_name(k, escaping), _escape(v, ALLOWUTF8, _is_legacy_labelname_rune))
                            for k, v in items])
                if labelstr:
                    labelstr = "{" + labelstr + "}"
                if s.exemplar:
                    exemplarstr = _compose_exemplar_string(metric, s, s.exemplar)
                else:
                    exemplarstr = ''
                timestamp = ''
                if s.timestamp is not None:
                    timestamp = f' {s.timestamp}'
                
                # Skip native histogram samples entirely if version < 2.0.0
                if s.native_histogram and parse_version(version) < (2, 0, 0):
                    continue
                
                native_histogram = ''
                negative_spans = ''
                negative_deltas = ''
                positive_spans = ''
                positive_deltas = ''
                     
                if s.native_histogram:
                    # Initialize basic nh template
                    nh_sample_template = '{{count:{},sum:{},schema:{},zero_threshold:{},zero_count:{}'

                    args = [
                        s.native_histogram.count_value,
                        s.native_histogram.sum_value,
                        s.native_histogram.schema,
                        s.native_histogram.zero_threshold,
                        s.native_histogram.zero_count,
                    ]
                  
                    # If there are neg spans, append them and the neg deltas to the template and args
                    if s.native_histogram.neg_spans:
                        negative_spans = ','.join([f'{ns[0]}:{ns[1]}' for ns in s.native_histogram.neg_spans])
                        negative_deltas = ','.join(str(nd) for nd in s.native_histogram.neg_deltas)
                        nh_sample_template += ',negative_spans:[{}]'
                        args.append(negative_spans)
                        nh_sample_template += ',negative_deltas:[{}]'
                        args.append(negative_deltas)

                    # If there are pos spans, append them and the pos spans to the template and args
                    if s.native_histogram.pos_spans:
                        positive_spans = ','.join([f'{ps[0]}:{ps[1]}' for ps in s.native_histogram.pos_spans])
                        positive_deltas = ','.join(f'{pd}' for pd in s.native_histogram.pos_deltas)
                        nh_sample_template += ',positive_spans:[{}]'
                        args.append(positive_spans)
                        nh_sample_template += ',positive_deltas:[{}]'
                        args.append(positive_deltas)                       
                                                  
                    # Add closing brace
                    nh_sample_template += '}}'

                    # Format the template with the args
                    native_histogram = nh_sample_template.format(*args)
                    
                    if s.native_histogram.nh_exemplars:
                        for nh_ex in s.native_histogram.nh_exemplars:
                            nh_exemplarstr = _compose_exemplar_string(metric, s, nh_ex)
                            exemplarstr += nh_exemplarstr

                value = ''
                if s.native_histogram:
                    value = native_histogram
                elif s.value is not None:
                    value = floatToGoString(s.value)
                if (escaping != ALLOWUTF8) or _is_valid_legacy_metric_name(s.name):
                    output.append('{}{} {}{}{}\n'.format(
                        _escape(s.name, escaping, _is_legacy_labelname_rune),
                        labelstr,
                        value,
                        timestamp,
                        exemplarstr
                    ))
                else:
                    output.append('{} {}{}{}\n'.format(
                        labelstr,
                        value,
                        timestamp,
                        exemplarstr
                    ))
        except Exception as exception:
            exception.args = (exception.args or ('',)) + (metric,)
            raise

    output.append('# EOF\n')
    return ''.join(output).encode('utf-8')


def escape_metric_name(s: str, escaping: str = UNDERSCORES) -> str:
    """Escapes the metric name and puts it in quotes iff the name does not
    conform to the legacy Prometheus character set.
    """
    if len(s) == 0:
        return s
    if escaping == ALLOWUTF8:
        if not _is_valid_legacy_metric_name(s):
            return '"{}"'.format(_escape(s, escaping, _is_legacy_metric_rune))
        return _escape(s, escaping, _is_legacy_metric_rune)
    elif escaping == UNDERSCORES:
        if _is_valid_legacy_metric_name(s):
            return s
        return _escape(s, escaping, _is_legacy_metric_rune)
    elif escaping == DOTS:
        return _escape(s, escaping, _is_legacy_metric_rune)
    elif escaping == VALUES:
        if _is_valid_legacy_metric_name(s):
            return s
        return _escape(s, escaping, _is_legacy_metric_rune)
    return s


def escape_label_name(s: str, escaping: str = UNDERSCORES) -> str:
    """Escapes the label name and puts it in quotes iff the name does not
    conform to the legacy Prometheus character set.
    """
    if len(s) == 0:
        return s
    if escaping == ALLOWUTF8:
        if not _is_valid_legacy_labelname(s):
            return '"{}"'.format(_escape(s, escaping, _is_legacy_labelname_rune))
        return _escape(s, escaping, _is_legacy_labelname_rune)
    elif escaping == UNDERSCORES:
        if _is_valid_legacy_labelname(s):
            return s
        return _escape(s, escaping, _is_legacy_labelname_rune)
    elif escaping == DOTS:
        return _escape(s, escaping, _is_legacy_labelname_rune)
    elif escaping == VALUES:
        if _is_valid_legacy_labelname(s):
            return s
        return _escape(s, escaping, _is_legacy_labelname_rune)
    return s


def _escape(s: str, escaping: str, valid_rune_fn: Callable[[str, int], bool]) -> str:
    """Performs backslash escaping on backslash, newline, and double-quote characters.

    valid_rune_fn takes the input character and its index in the containing string."""
    if escaping == ALLOWUTF8:
        return s.replace('\\', r'\\').replace('\n', r'\n').replace('"', r'\"')
    elif escaping == UNDERSCORES:
        escaped = StringIO()
        for i, b in enumerate(s):
            if valid_rune_fn(b, i):
                escaped.write(b)
            else:
                escaped.write('_')
        return escaped.getvalue()
    elif escaping == DOTS:
        escaped = StringIO()
        for i, b in enumerate(s):
            if b == '_':
                escaped.write('__')
            elif b == '.':
                escaped.write('_dot_')
            elif valid_rune_fn(b, i):
                escaped.write(b)
            else:
                escaped.write('__')
        return escaped.getvalue()
    elif escaping == VALUES:
        escaped = StringIO()
        escaped.write("U__")
        for i, b in enumerate(s):
            if b == '_':
                escaped.write("__")
            elif valid_rune_fn(b, i):
                escaped.write(b)
            elif not _is_valid_utf8(b):
                escaped.write("_FFFD_")
            else:
                escaped.write('_')
                escaped.write(format(ord(b), 'x'))
                escaped.write('_')
        return escaped.getvalue()
    return s


def _is_legacy_metric_rune(b: str, i: int) -> bool:
    return _is_legacy_labelname_rune(b, i) or b == ':'


def _is_legacy_labelname_rune(b: str, i: int) -> bool:
    if len(b) != 1:
        raise ValueError("Input 'b' must be a single character.")
    return (
        ('a' <= b <= 'z')
        or ('A' <= b <= 'Z')
        or (b == '_')
        or ('0' <= b <= '9' and i > 0)
    )


_SURROGATE_MIN = 0xD800
_SURROGATE_MAX = 0xDFFF


def _is_valid_utf8(s: str) -> bool:
    if 0 <= ord(s) < _SURROGATE_MIN:
        return True
    if _SURROGATE_MAX < ord(s) <= maxunicode:
        return True
    return False
