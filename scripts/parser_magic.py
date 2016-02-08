from pandas.util.testing import set_trace
import pandas.util.testing as tm
import pandas.compat as compat

from pandas import *
import ast
import inspect
import sys


def merge(a, b):
    f, args, _ = parse_stmt(inspect.currentframe().f_back)
    return DataFrame({args[0]: a,
                      args[1]: b})


def parse_stmt(frame):
    info = inspect.getframeinfo(frame)
    call = info[-2][0]
    mod = ast.parse(call)
    body = mod.body[0]
    if isinstance(body, (ast.Assign, ast.Expr)):
        call = body.value
    elif isinstance(body, ast.Call):
        call = body
    return _parse_call(call)


def _parse_call(call):
    func = _maybe_format_attribute(call.func)

    str_args = []
    for arg in call.args:
        if isinstance(arg, ast.Name):
            str_args.append(arg.id)
        elif isinstance(arg, ast.Call):
            formatted = _format_call(arg)
            str_args.append(formatted)

    return func, str_args, {}


def _format_call(call):
    func, args, kwds = _parse_call(call)
    content = ''
    if args:
        content += ', '.join(args)
    if kwds:
        fmt_kwds = ['%s=%s' % item for item in compat.iteritems(kwds)]
        joined_kwds = ', '.join(fmt_kwds)
        if args:
            content = content + ', ' + joined_kwds
        else:
            content += joined_kwds
    return '%s(%s)' % (func, content)


def _maybe_format_attribute(name):
    if isinstance(name, ast.Attribute):
        return _format_attribute(name)
    return name.id


def _format_attribute(attr):
    obj = attr.value
    if isinstance(attr.value, ast.Attribute):
        obj = _format_attribute(attr.value)
    else:
        obj = obj.id
    return '.'.join((obj, attr.attr))

a = tm.makeTimeSeries()
b = tm.makeTimeSeries()
df = merge(a, b)
