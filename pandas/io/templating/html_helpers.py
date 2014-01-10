from __future__ import print_function

import sys, os, re, json, codecs
import uuid

from pandas.io.templating import *

__all__ = ["zebra", "hlrow", "hlcol"]

from collections import namedtuple


def zebra(s, color1, color2):
    style = [dict(selector="td.%s:nth-child(2n)" % DATA_CLASS,
                  props=[("background-color", color1)]),
             dict(selector="td.%s:nth-child(2n+1)" % DATA_CLASS,
                  props=[("background-color", color2)])]
    s.style.extend(style)


def hlcell(s, r, c, color="#aaa", with_headings=False):
    selector = "td.%s%d.%s%d" % (ROW_CLASS, r, COLUMN_CLASS, c)
    if not with_headings:
        selector += ".%s" % DATA_CLASS
    style = [dict(selector=selector,
                  props=[("background-color", color)])]
    s.style.extend(style)


def hlcol(s, n, color="#aaa", with_headings=False):
    selector = "td.%s%d" % (COLUMN_CLASS, n)
    if not with_headings:
        selector += ".%s" % DATA_CLASS
    style = [dict(selector=selector,
                  props=[("background-color", color)])]
    s.style.extend(style)


def hlrow(s, n, color="#ccc", with_headings=False):
    selector = "td.%s%d" % (ROW_CLASS, n)
    if not with_headings:
        selector += ".%s" % DATA_CLASS
    style = [dict(selector=selector,
                  props=[("background-color", color)])]
    s.style.extend(style)


def round_corners(s, radius):
    props_bl = [
        ("-moz-border-radius-bottomleft", "%dpx" % radius ),
        ("-webkit-border-bottom-left-radius", "%dpx" % radius ),
        ("border-bottom-left-radius", "%dpx" % radius )
    ]
    props_br = [
        ("-moz-border-radius-bottomright", "%dpx" % radius ),
        ("-webkit-border-bottom-right-radius", "%dpx" % radius ),
        ("border-bottom-right-radius", "%dpx" % radius )
    ]
    props_tl = [
        ("-moz-border-radius-topleft", "%dpx" % radius ),
        ("-webkit-border-top-left-radius", "%dpx" % radius ),
        ("border-top-left-radius", "%dpx" % radius )
    ]
    props_tr = [
        ("-moz-border-radius-topright", "%dpx" % radius ),
        ("-webkit-border-top-right-radius", "%dpx" % radius ),
        ("border-top-right-radius", "%dpx" % radius )
    ]

    style = [
        dict(selector="",
             props=[("border-collapse", "separate")]),
        dict(selector="td",
             props=[("border-width", "0px")]),
        dict(selector="th",
             props=[("border-width", "0px")]),
        dict(selector="td",
             props=[("border-left-width", "1px")]),
        dict(selector="td",
             props=[("border-top-width", "1px")]),

        dict(selector="tbody tr:last-child th",
             props=[("border-bottom-width", "1px")]),
        dict(selector="tr:last-child td",
             props=[("border-bottom-width", "1px")]),

        dict(selector="tr td:last-child",
             props=[("border-right-width", "1px")]),
        dict(selector="tr th:last-child",
             props=[("border-right-width", "1px")]),
        dict(selector="th",
             props=[("border-left-width", "1px")]),
        dict(selector="th",
             props=[("border-top-width", "1px")]),
        dict(selector="th td:last-child",
             props=[("border-right-width", "1px")]),


        dict(selector="tr:last-child th:first-child",
             props=props_bl),
        dict(selector="tr:last-child td:last-child",
             props=props_br),
        dict(selector="tr:first-child th.%s0" % COLUMN_CLASS,
             props=props_tl),
        dict(selector="tr:first-child th.%s0:first-child" % ROW_CLASS,
             props=props_tl),
        dict(selector="tr:first-child th:last-child",
             props=props_tr),
    ]

    s.style.extend(style)

#
# def rank_heatmap(s, df, row=None, col=None):
#     def color_class(cls, color):
#         return [dict(selector="td.%s" % cls,
#                      props=[("background-color", color)])]
#
#     def rank_col(n, ranking, u):
#         data = {i: {n: ["%s-%s" % (u, ranking[i])]}
#                 for i in range(len(ranking))}
#         return {"data": data}
#
#
#     # u = "U" + str(uuid.uuid1()).replace("-", "_")
#     # df = mkdf(9, 5, data_gen_f=lambda r, c: np.random.random())
#
#     ranking = df.iloc[:, 1].argsort().tolist()
#     cell_context = rank_col(1, ranking, u)
#
#     from .colors import red_scale
#
#     style = [color_class("%s-%s" % (u, intensity),
#                          red_scale[intensity])
#              for intensity in range(9)]
#
#     s.style.extend(style)
#     # s.cell_context.extend(s) # TODO
