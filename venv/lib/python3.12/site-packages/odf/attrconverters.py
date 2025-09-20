# -*- coding: utf-8 -*-
# Copyright (C) 2006-2013 Søren Roug, European Environment Agency
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contributor(s):
#

import sys, os.path
sys.path.append(os.path.dirname(__file__))
from odf.namespaces import *
import re, types

pattern_color =  re.compile(r'#[0-9a-fA-F]{6}')
pattern_vector3D = re.compile(r'\([ ]*-?([0-9]+(\.[0-9]*)?|\.[0-9]+)([ ]+-?([0-9]+(\.[0-9]*)?|\.[0-9]+)){2}[ ]*\)')

def make_NCName(arg):
    for c in (':',' '):
        arg = arg.replace(c,"_%x_" % ord(c))
    return arg

def cnv_angle(attribute, arg, element):
        if sys.version_info[0]==2:
            return unicode(arg)
        else:
            return str(arg)

def cnv_anyURI(attribute, arg, element):
    return str(arg)

def cnv_boolean(attribute, arg, element):
    """ XML Schema Part 2: Datatypes Second Edition
        An instance of a datatype that is defined as boolean can have the
        following legal literals {true, false, 1, 0}
    """
    if str(arg).lower() in ("0","false","no"):
        return "false"
    if str(arg).lower() in ("1","true","yes"):
        return "true"
    raise ValueError( "'%s' not allowed as Boolean value for %s" % (str(arg), attribute[1]))

# Potentially accept color values
def cnv_color(attribute, arg, element):
    """ A RGB color in conformance with §5.9.11 of [XSL], that is a RGB color in notation “#rrggbb”, where
        rr, gg and bb are 8-bit hexadecimal digits.
    """
    return str(arg)

def cnv_configtype(attribute, arg, element):
    if str(arg) not in ("boolean", "short", "int", "long",
    "double", "string", "datetime", "base64Binary"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_data_source_has_labels(attribute, arg, element):
    if str(arg) not in ("none","row","column","both"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

# Understand different date formats
def cnv_date(attribute, arg, element):
    """ A dateOrDateTime value is either an [xmlschema-2] date value or an [xmlschema-2] dateTime
        value.
    """
    return str(arg)

def cnv_dateTime(attribute, arg, element):
    """ A dateOrDateTime value is either an [xmlschema-2] date value or an [xmlschema-2] dateTime
        value.
    """
    return str(arg)

def cnv_double(attribute, arg, element):
    return str(arg)

def cnv_draw_aspect(attribute, arg, element):
    if str(arg) not in ("content", "thumbnail", "icon", "print-view"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_duration(attribute, arg, element):
    return str(arg)

def cnv_family(attribute, arg, element):
    """ A style family """
    if str(arg) not in ("text", "paragraph", "section", "ruby", "table", "table-column", "table-row", "table-cell",
      "graphic", "presentation", "drawing-page", "chart"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def __save_prefix(attribute, arg, element):
    prefix = arg.split(':',1)[0]
    if prefix == arg:
        return arg
    namespace = element.get_knownns(prefix)
    if namespace is None:
        #raise ValueError( "'%s' is an unknown prefix" % str(prefix))
        return str(arg)
    p = element.get_nsprefix(namespace)
    return str(arg)

def cnv_formula(attribute, arg, element):
    """ A string containing a formula. Formulas do not have a predefined syntax, but the string should
        begin with a namespace prefix, followed by a “:” (COLON, U+003A) separator, followed by the text
        of the formula. The namespace bound to the prefix determines the syntax and semantics of the
        formula.
    """
    return __save_prefix(attribute, arg, element)

def cnv_ID(attribute, arg, element):
    return str(arg)

def cnv_IDREF(attribute, arg, element):
    return str(arg)

def cnv_integer(attribute, arg, element):
    return str(arg)

pattern_language = re.compile(r'[a-zA-Z]{1,8}(-[a-zA-Z0-9]{1,8})*')

def cnv_language(attribute, arg, element):
    global pattern_language
    if not pattern_language.match(arg):
        raise ValueError( "'%s' is not a valid language token" % arg)
    return arg

def cnv_legend_position(attribute, arg, element):
    if str(arg) not in ("start", "end", "top", "bottom", "top-start", "bottom-start", "top-end", "bottom-end"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

pattern_length = re.compile(r'-?([0-9]+(\.[0-9]*)?|\.[0-9]+)((cm)|(mm)|(in)|(pt)|(pc)|(px))')

def cnv_length(attribute, arg, element):
    """ A (positive or negative) physical length, consisting of magnitude and unit, in conformance with the
        Units of Measure defined in §5.9.13 of [XSL].
    """
    global pattern_length
    if not pattern_length.match(arg):
        raise ValueError( "'%s' is not a valid length" % arg)
    return arg

def cnv_lengthorpercent(attribute, arg, element):
    failed = False
    try: return cnv_length(attribute, arg, element)
    except: failed = True
    try: return cnv_percent(attribute, arg, element)
    except: failed = True
    if failed:
        raise ValueError( "'%s' is not a valid length or percent" % arg)
    return arg

def cnv_list_linkage_type(attribute, arg, element):
    if arg not in ('selection','selection-indices'):
        raise ValueError( "'%s' is not either 'selection' or 'selection-indices'" % arg)
    return str(arg)

def cnv_metavaluetype(attribute, arg, element):
    if str(arg) not in ("float", "date", "time", "boolean", "string"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_major_minor(attribute, arg, element):
    if arg not in ('major','minor'):
        raise ValueError( "'%s' is not either 'minor' or 'major'" % arg)
    return str(arg)

pattern_namespacedToken = re.compile(r'[0-9a-zA-Z_]+:[0-9a-zA-Z._\-]+')

def cnv_namespacedToken(attribute, arg, element):
    global pattern_namespacedToken

    if not pattern_namespacedToken.match(arg):
        raise ValueError( "'%s' is not a valid namespaced token" % arg)
    return __save_prefix(attribute, arg, element)

def cnv_NCName(attribute, arg, element):
    """ NCName is defined in http://www.w3.org/TR/REC-xml-names/#NT-NCName
        Essentially an XML name minus ':'
    """
    if (sys.version_info[0]==3 and isinstance(arg, str)) or (sys.version_info[0]==2 and type(arg) in types.StringTypes):
        return make_NCName(arg)
    else:
        return arg.getAttrNS(STYLENS, 'name')

# This function takes either an instance of a style (preferred)
# or a text string naming the style. If it is a text string, then it must
# already have been converted to an NCName
# The text-string argument is mainly for when we build a structure from XML
def cnv_StyleNameRef(attribute, arg, element):
    try:
        return arg.getAttrNS(STYLENS, 'name')
    except:
        return arg

# This function takes either an instance of a style (preferred)
# or a text string naming the style. If it is a text string, then it must
# already have been converted to an NCName
# The text-string argument is mainly for when we build a structure from XML
def cnv_DrawNameRef(attribute, arg, element):
    try:
        return arg.getAttrNS(DRAWNS, 'name')
    except:
        return arg

# Must accept list of Style objects
def cnv_NCNames(attribute, arg, element):
    return ' '.join(arg)

def cnv_nonNegativeInteger(attribute, arg, element):
    return str(arg)

pattern_percent = re.compile(r'-?([0-9]+(\.[0-9]*)?|\.[0-9]+)%')

def cnv_percent(attribute, arg, element):
    global pattern_percent
    if not pattern_percent.match(arg):
        raise ValueError( "'%s' is not a valid length" % arg)
    return arg

# Real one doesn't allow floating point values
pattern_points = re.compile(r'-?[0-9]+,-?[0-9]+([ ]+-?[0-9]+,-?[0-9]+)*')
#pattern_points = re.compile(r'-?[0-9.]+,-?[0-9.]+([ ]+-?[0-9.]+,-?[0-9.]+)*')
def cnv_points(attribute, arg, element):
    global pattern_points
    if (sys.version_info[0]==3 and isinstance(arg, str)) or (sys.version_info[0]==2 and type(arg) in types.StringTypes):
        if not pattern_points.match(arg):
            raise ValueError( "x,y are separated by a comma and the points are separated by white spaces")
        return arg
    else:
        try:
            strarg = ' '.join([ "%d,%d" % p for p in arg])
        except:
            raise ValueError( "Points must be string or [(0,0),(1,1)] - not %s" % arg)
        return strarg

def cnv_positiveInteger(attribute, arg, element):
    return str(arg)

def cnv_rowOrCol(attribute, arg, element):
    if str(arg) not in ("row","column"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_string(attribute, arg, element):
    if sys.version_info[0]==2:
        return unicode(arg)
    else:
        return str(arg)

def cnv_stroke_linecap(attribute, arg, element):
    if str(arg) not in ("butt", "square", "round"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_textnoteclass(attribute, arg, element):
    if str(arg) not in ("footnote", "endnote"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

# Understand different time formats
def cnv_time(attribute, arg, element):
    return str(arg)

def cnv_token(attribute, arg, element):
    return str(arg)

pattern_viewbox = re.compile(r'-?[0-9]+([ ]+-?[0-9]+){3}$')

def cnv_viewbox(attribute, arg, element):
    global pattern_viewbox
    if not pattern_viewbox.match(arg):
        raise ValueError( "viewBox must be four integers separated by whitespaces")
    return arg

def cnv_xlinkshow(attribute, arg, element):
    if str(arg) not in ("new", "replace", "embed"):
        raise ValueError( "'%s' not allowed" % str(arg))
    return str(arg)

def cnv_xlinktype(attribute, arg, element):
    if arg != "simple":
        raise ValueError( "Value of '%s' must be 'simple'" % attribute[1])
    return arg


attrconverters = {
	((ANIMNS,u'audio-level'), None): cnv_double,
	((ANIMNS,u'color-interpolation'), None): cnv_string,
	((ANIMNS,u'color-interpolation-direction'), None): cnv_string,
	((ANIMNS,u'command'), None): cnv_string,
	((ANIMNS,u'formula'), None): cnv_string,
	((ANIMNS,u'id'), None): cnv_ID,
	((ANIMNS,u'iterate-interval'), None): cnv_duration,
	((ANIMNS,u'iterate-type'), None): cnv_string,
	((ANIMNS,u'name'), None): cnv_string,
	((ANIMNS,u'sub-item'), None): cnv_string,
	((ANIMNS,u'value'), None): cnv_string,
#	((DBNS,u'type'), None): cnv_namespacedToken,
	((CHARTNS,u'angle-offset'), None): cnv_angle,
	((CHARTNS,u'automatic-content'), None): cnv_boolean,
	((CHARTNS,u'auto-position'), None): cnv_boolean,
	((CHARTNS,u'auto-size'), None): cnv_boolean,
	((CHARTNS,u'axis-label-position'), None): cnv_string, # Multi-value
	((CHARTNS,u'axis-position'), None): cnv_string, # Multi-value
	((CHARTNS,u'attached-axis'), None): cnv_string,
	((CHARTNS,u'class'), (CHARTNS,u'grid')): cnv_major_minor,
	((CHARTNS,u'class'), None): cnv_namespacedToken,
	((CHARTNS,u'column-mapping'), None): cnv_string,
	((CHARTNS,u'connect-bars'), None): cnv_boolean,
	((CHARTNS,u'data-label-number'), None): cnv_string,
	((CHARTNS,u'data-label-symbol'), None): cnv_boolean,
	((CHARTNS,u'data-label-text'), None): cnv_boolean,
	((CHARTNS,u'data-source-has-labels'), None): cnv_data_source_has_labels,
	((CHARTNS,u'deep'), None): cnv_boolean,
	((CHARTNS,u'dimension'), None): cnv_string,
	((CHARTNS,u'display-equation'), None): cnv_boolean,
	((CHARTNS,u'display-label'), None): cnv_boolean,
	((CHARTNS,u'display-r-square'), None): cnv_boolean,
	((CHARTNS,u'error-category'), None): cnv_string,
	((CHARTNS,u'error-lower-indicator'), None): cnv_boolean,
	((CHARTNS,u'error-lower-limit'), None): cnv_string,
	((CHARTNS,u'error-margin'), None): cnv_string,
	((CHARTNS,u'error-percentage'), None): cnv_string,
	((CHARTNS,u'error-lower-range'), None): cnv_string,
	((CHARTNS,u'error-upper-indicator'), None): cnv_boolean,
	((CHARTNS,u'error-upper-limit'), None): cnv_string,
	((CHARTNS,u'error-upper-range'), None): cnv_string,
	((CHARTNS,u'gap-width'), None): cnv_string,
	((CHARTNS,u'group-bars-per-axis'), None): cnv_boolean,
	((CHARTNS,u'hole-size'), None): cnv_percent,
	((CHARTNS,u'include-hidden-cells'), None): cnv_boolean,
	((CHARTNS,u'interpolation'), None): cnv_string,
	((CHARTNS,u'interval-major'), None): cnv_string,
	((CHARTNS,u'interval-minor-divisor'), None): cnv_string,
	((CHARTNS,u'japanese-candle-stick'), None): cnv_boolean,
	((CHARTNS,u'label-arrangement'), None): cnv_string,
	((CHARTNS,u'label-cell-address'), None): cnv_string,
	((CHARTNS,u'label-position'), None): cnv_string, # Multi-value
	((CHARTNS,u'label-position-negative'), None): cnv_string, # Multi-value
	((CHARTNS,u'legend-align'), None): cnv_string,
	((CHARTNS,u'legend-position'), None): cnv_legend_position,
	((CHARTNS,u'lines'), None): cnv_boolean,
	((CHARTNS,u'link-data-style-to-source'), None): cnv_boolean,
	((CHARTNS,u'logarithmic'), None): cnv_boolean,
	((CHARTNS,u'maximum'), None): cnv_string,
	((CHARTNS,u'mean-value'), None): cnv_boolean,
	((CHARTNS,u'minimum'), None): cnv_string,
	((CHARTNS,u'name'), None): cnv_string,
	((CHARTNS,u'origin'), None): cnv_string,
	((CHARTNS,u'overlap'), None): cnv_string,
	((CHARTNS,u'percentage'), None): cnv_boolean,
	((CHARTNS,u'pie-offset'), None): cnv_string,
	((CHARTNS,u'regression-type'), None): cnv_string,
	((CHARTNS,u'repeated'), None): cnv_nonNegativeInteger,
	((CHARTNS,u'reverse-direction'), None): cnv_boolean,
	((CHARTNS,u'right-angled-axes'), None): cnv_boolean,
	((CHARTNS,u'row-mapping'), None): cnv_string,
	((CHARTNS,u'scale-text'), None): cnv_boolean,
	((CHARTNS,u'series-source'), None): cnv_string,
	((CHARTNS,u'solid-type'), None): cnv_string,
	((CHARTNS,u'sort-by-x-values'), None): cnv_boolean,
	((CHARTNS,u'spline-order'), None): cnv_string,
	((CHARTNS,u'spline-resolution'), None): cnv_string,
	((CHARTNS,u'stacked'), None): cnv_boolean,
	((CHARTNS,u'style-name'), None): cnv_StyleNameRef,
	((CHARTNS,u'symbol-height'), None): cnv_string,
	((CHARTNS,u'symbol-name'), None): cnv_string,
	((CHARTNS,u'symbol-type'), None): cnv_string,
	((CHARTNS,u'symbol-width'), None): cnv_string,
	((CHARTNS,u'text-overlap'), None): cnv_boolean,
	((CHARTNS,u'three-dimensional'), None): cnv_boolean,
	((CHARTNS,u'tick-mark-position'), None): cnv_string, # Multi-value
	((CHARTNS,u'tick-marks-major-inner'), None): cnv_boolean,
	((CHARTNS,u'tick-marks-major-outer'), None): cnv_boolean,
	((CHARTNS,u'tick-marks-minor-inner'), None): cnv_boolean,
	((CHARTNS,u'tick-marks-minor-outer'), None): cnv_boolean,
	((CHARTNS,u'treat-empty-cells'), None): cnv_string, # Multi-value
	((CHARTNS,u'values-cell-range-address'), None): cnv_string,
	((CHARTNS,u'vertical'), None): cnv_boolean,
	((CHARTNS,u'visible'), None): cnv_boolean,
	((CONFIGNS,u'name'), None): cnv_formula,
	((CONFIGNS,u'type'), None): cnv_configtype,
	((DR3DNS,u'ambient-color'), None): cnv_string,
	((DR3DNS,u'back-scale'), None): cnv_string,
	((DR3DNS,u'backface-culling'), None): cnv_string,
	((DR3DNS,u'center'), None): cnv_string,
	((DR3DNS,u'close-back'), None): cnv_boolean,
	((DR3DNS,u'close-front'), None): cnv_boolean,
	((DR3DNS,u'depth'), None): cnv_length,
	((DR3DNS,u'diffuse-color'), None): cnv_string,
	((DR3DNS,u'direction'), None): cnv_string,
	((DR3DNS,u'distance'), None): cnv_length,
	((DR3DNS,u'edge-rounding'), None): cnv_string,
	((DR3DNS,u'edge-rounding-mode'), None): cnv_string,
	((DR3DNS,u'emissive-color'), None): cnv_string,
	((DR3DNS,u'enabled'), None): cnv_boolean,
	((DR3DNS,u'end-angle'), None): cnv_string,
	((DR3DNS,u'focal-length'), None): cnv_length,
	((DR3DNS,u'horizontal-segments'), None): cnv_string,
	((DR3DNS,u'lighting-mode'), None): cnv_boolean,
	((DR3DNS,u'max-edge'), None): cnv_string,
	((DR3DNS,u'min-edge'), None): cnv_string,
	((DR3DNS,u'normals-direction'), None): cnv_string,
	((DR3DNS,u'normals-kind'), None): cnv_string,
	((DR3DNS,u'projection'), None): cnv_string,
	((DR3DNS,u'shade-mode'), None): cnv_string,
	((DR3DNS,u'shadow'), None): cnv_string,
	((DR3DNS,u'shadow-slant'), None): cnv_nonNegativeInteger,
	((DR3DNS,u'shininess'), None): cnv_string,
	((DR3DNS,u'size'), None): cnv_string,
	((DR3DNS,u'specular'), None): cnv_boolean,
	((DR3DNS,u'specular-color'), None): cnv_string,
	((DR3DNS,u'texture-filter'), None): cnv_string,
	((DR3DNS,u'texture-generation-mode-x'), None): cnv_string,
	((DR3DNS,u'texture-generation-mode-y'), None): cnv_string,
	((DR3DNS,u'texture-kind'), None): cnv_string,
	((DR3DNS,u'texture-mode'), None): cnv_string,
	((DR3DNS,u'transform'), None): cnv_string,
	((DR3DNS,u'vertical-segments'), None): cnv_string,
	((DR3DNS,u'vpn'), None): cnv_string,
	((DR3DNS,u'vrp'), None): cnv_string,
	((DR3DNS,u'vup'), None): cnv_string,
	((DRAWNS,u'align'), None): cnv_string,
	((DRAWNS,u'angle'), None): cnv_integer,
	((DRAWNS,u'archive'), None): cnv_string,
	((DRAWNS,u'auto-grow-height'), None): cnv_boolean,
	((DRAWNS,u'auto-grow-width'), None): cnv_boolean,
	((DRAWNS,u'background-size'), None): cnv_string,
	((DRAWNS,u'blue'), None): cnv_string,
	((DRAWNS,u'border'), None): cnv_string,
	((DRAWNS,u'caption-angle'), None): cnv_string,
	((DRAWNS,u'caption-angle-type'), None): cnv_string,
	((DRAWNS,u'caption-escape'), None): cnv_string,
	((DRAWNS,u'caption-escape-direction'), None): cnv_string,
	((DRAWNS,u'caption-fit-line-length'), None): cnv_boolean,
	((DRAWNS,u'caption-gap'), None): cnv_string,
	((DRAWNS,u'caption-line-length'), None): cnv_length,
	((DRAWNS,u'caption-point-x'), None): cnv_string,
	((DRAWNS,u'caption-point-y'), None): cnv_string,
	((DRAWNS,u'caption-id'), None): cnv_IDREF,
	((DRAWNS,u'caption-type'), None): cnv_string,
	((DRAWNS,u'chain-next-name'), None): cnv_string,
	((DRAWNS,u'class-id'), None): cnv_string,
	((DRAWNS,u'class-names'), None): cnv_NCNames,
	((DRAWNS,u'code'), None): cnv_string,
	((DRAWNS,u'color'), None): cnv_string,
	((DRAWNS,u'color-inversion'), None): cnv_boolean,
	((DRAWNS,u'color-mode'), None): cnv_string,
	((DRAWNS,u'concave'), None): cnv_string,
	((DRAWNS,u'concentric-gradient-fill-allowed'), None): cnv_boolean,
	((DRAWNS,u'contrast'), None): cnv_string,
	((DRAWNS,u'control'), None): cnv_IDREF,
	((DRAWNS,u'copy-of'), None): cnv_string,
	((DRAWNS,u'corner-radius'), None): cnv_length,
	((DRAWNS,u'corners'), None): cnv_positiveInteger,
	((DRAWNS,u'cx'), None): cnv_string,
	((DRAWNS,u'cy'), None): cnv_string,
	((DRAWNS,u'data'), None): cnv_string,
	((DRAWNS,u'decimal-places'), None): cnv_string,
	((DRAWNS,u'display'), None): cnv_string,
	((DRAWNS,u'display-name'), None): cnv_string,
	((DRAWNS,u'distance'), None): cnv_lengthorpercent,
	((DRAWNS,u'dots1'), None): cnv_integer,
	((DRAWNS,u'dots1-length'), None): cnv_lengthorpercent,
	((DRAWNS,u'dots2'), None): cnv_integer,
	((DRAWNS,u'dots2-length'), None): cnv_lengthorpercent,
	((DRAWNS,u'draw-aspect'), None): cnv_draw_aspect,
	((DRAWNS,u'end-angle'), None): cnv_angle,
	((DRAWNS,u'end'), None): cnv_string,
	((DRAWNS,u'end-color'), None): cnv_string,
	((DRAWNS,u'end-glue-point'), None): cnv_nonNegativeInteger,
	((DRAWNS,u'end-guide'), None): cnv_length,
	((DRAWNS,u'end-intensity'), None): cnv_string,
	((DRAWNS,u'end-line-spacing-horizontal'), None): cnv_string,
	((DRAWNS,u'end-line-spacing-vertical'), None): cnv_string,
	((DRAWNS,u'end-shape'), None): cnv_IDREF,
	((DRAWNS,u'engine'), None): cnv_namespacedToken,
	((DRAWNS,u'enhanced-path'), None): cnv_string,
	((DRAWNS,u'escape-direction'), None): cnv_string,
	((DRAWNS,u'extrusion-allowed'), None): cnv_boolean,
	((DRAWNS,u'extrusion-brightness'), None): cnv_string,
	((DRAWNS,u'extrusion'), None): cnv_boolean,
	((DRAWNS,u'extrusion-color'), None): cnv_boolean,
	((DRAWNS,u'extrusion-depth'), None): cnv_double,
	((DRAWNS,u'extrusion-diffusion'), None): cnv_string,
	((DRAWNS,u'extrusion-first-light-direction'), None): cnv_string,
	((DRAWNS,u'extrusion-first-light-harsh'), None): cnv_boolean,
	((DRAWNS,u'extrusion-first-light-level'), None): cnv_string,
	((DRAWNS,u'extrusion-light-face'), None): cnv_boolean,
	((DRAWNS,u'extrusion-metal'), None): cnv_boolean,
	((DRAWNS,u'extrusion-number-of-line-segments'), None): cnv_integer,
	((DRAWNS,u'extrusion-origin'), None): cnv_double,
	((DRAWNS,u'extrusion-rotation-angle'), None): cnv_double,
	((DRAWNS,u'extrusion-rotation-center'), None): cnv_string,
	((DRAWNS,u'extrusion-second-light-direction'), None): cnv_string,
	((DRAWNS,u'extrusion-second-light-harsh'), None): cnv_boolean,
	((DRAWNS,u'extrusion-second-light-level'), None): cnv_string,
	((DRAWNS,u'extrusion-shininess'), None): cnv_string,
	((DRAWNS,u'extrusion-skew'), None): cnv_double,
	((DRAWNS,u'extrusion-specularity'), None): cnv_string,
	((DRAWNS,u'extrusion-viewpoint'), None): cnv_string,
	((DRAWNS,u'fill'), None): cnv_string,
	((DRAWNS,u'fill-color'), None): cnv_string,
	((DRAWNS,u'fill-gradient-name'), None): cnv_string,
	((DRAWNS,u'fill-hatch-name'), None): cnv_string,
	((DRAWNS,u'fill-hatch-solid'), None): cnv_boolean,
	((DRAWNS,u'fill-image-height'), None): cnv_lengthorpercent,
	((DRAWNS,u'fill-image-name'), None): cnv_DrawNameRef,
	((DRAWNS,u'fill-image-ref-point'), None): cnv_string,
	((DRAWNS,u'fill-image-ref-point-x'), None): cnv_string,
	((DRAWNS,u'fill-image-ref-point-y'), None): cnv_string,
	((DRAWNS,u'fill-image-width'), None): cnv_lengthorpercent,
	((DRAWNS,u'filter-name'), None): cnv_string,
	((DRAWNS,u'fit-to-contour'), None): cnv_boolean,
	((DRAWNS,u'fit-to-size'), None): cnv_string,  # ODF 1.2 says boolean
	((DRAWNS,u'formula'), None): cnv_string,
	((DRAWNS,u'frame-display-border'), None): cnv_boolean,
	((DRAWNS,u'frame-display-scrollbar'), None): cnv_boolean,
	((DRAWNS,u'frame-margin-horizontal'), None): cnv_string,
	((DRAWNS,u'frame-margin-vertical'), None): cnv_string,
	((DRAWNS,u'frame-name'), None): cnv_string,
	((DRAWNS,u'gamma'), None): cnv_string,
	((DRAWNS,u'glue-point-leaving-directions'), None): cnv_string,
	((DRAWNS,u'glue-point-type'), None): cnv_string,
	((DRAWNS,u'glue-points'), None): cnv_string,
	((DRAWNS,u'gradient-step-count'), None): cnv_string,
	((DRAWNS,u'green'), None): cnv_string,
	((DRAWNS,u'guide-distance'), None): cnv_string,
	((DRAWNS,u'guide-overhang'), None): cnv_length,
	((DRAWNS,u'handle-mirror-horizontal'), None): cnv_boolean,
	((DRAWNS,u'handle-mirror-vertical'), None): cnv_boolean,
	((DRAWNS,u'handle-polar'), None): cnv_string,
	((DRAWNS,u'handle-position'), None): cnv_string,
	((DRAWNS,u'handle-radius-range-maximum'), None): cnv_string,
	((DRAWNS,u'handle-radius-range-minimum'), None): cnv_string,
	((DRAWNS,u'handle-range-x-maximum'), None): cnv_string,
	((DRAWNS,u'handle-range-x-minimum'), None): cnv_string,
	((DRAWNS,u'handle-range-y-maximum'), None): cnv_string,
	((DRAWNS,u'handle-range-y-minimum'), None): cnv_string,
	((DRAWNS,u'handle-switched'), None): cnv_boolean,
#	((DRAWNS,u'id'), None): cnv_ID,
#	((DRAWNS,u'id'), None): cnv_nonNegativeInteger,   # ?? line 6581 in RNG
	((DRAWNS,u'id'), None): cnv_string,
	((DRAWNS,u'image-opacity'), None): cnv_string,
	((DRAWNS,u'kind'), None): cnv_string,
	((DRAWNS,u'layer'), None): cnv_string,
	((DRAWNS,u'line-distance'), None): cnv_string,
	((DRAWNS,u'line-skew'), None): cnv_string,
	((DRAWNS,u'luminance'), None): cnv_string,
	((DRAWNS,u'marker-end-center'), None): cnv_boolean,
	((DRAWNS,u'marker-end'), None): cnv_string,
	((DRAWNS,u'marker-end-width'), None): cnv_length,
	((DRAWNS,u'marker-start-center'), None): cnv_boolean,
	((DRAWNS,u'marker-start'), None): cnv_string,
	((DRAWNS,u'marker-start-width'), None): cnv_length,
	((DRAWNS,u'master-page-name'), None): cnv_StyleNameRef,
	((DRAWNS,u'may-script'), None): cnv_boolean,
	((DRAWNS,u'measure-align'), None): cnv_string,
	((DRAWNS,u'measure-vertical-align'), None): cnv_string,
	((DRAWNS,u'mime-type'), None): cnv_string,
	((DRAWNS,u'mirror-horizontal'), None): cnv_boolean,
	((DRAWNS,u'mirror-vertical'), None): cnv_boolean,
	((DRAWNS,u'modifiers'), None): cnv_string,
	((DRAWNS,u'name'), None): cnv_NCName,
#	((DRAWNS,u'name'), None): cnv_string,
	((DRAWNS,u'nav-order'), None): cnv_IDREF,
	((DRAWNS,u'nohref'), None): cnv_string,
	((DRAWNS,u'notify-on-update-of-ranges'), None): cnv_string,
	((DRAWNS,u'object'), None): cnv_string,
	((DRAWNS,u'ole-draw-aspect'), None): cnv_string,
	((DRAWNS,u'opacity'), None): cnv_string,
	((DRAWNS,u'opacity-name'), None): cnv_string,
	((DRAWNS,u'page-number'), None): cnv_positiveInteger,
	((DRAWNS,u'parallel'), None): cnv_boolean,
	((DRAWNS,u'path-stretchpoint-x'), None): cnv_double,
	((DRAWNS,u'path-stretchpoint-y'), None): cnv_double,
	((DRAWNS,u'placing'), None): cnv_string,
	((DRAWNS,u'points'), None): cnv_points,
	((DRAWNS,u'protected'), None): cnv_boolean,
	((DRAWNS,u'recreate-on-edit'), None): cnv_boolean,
	((DRAWNS,u'red'), None): cnv_string,
	((DRAWNS,u'rotation'), None): cnv_integer,
	((DRAWNS,u'secondary-fill-color'), None): cnv_string,
	((DRAWNS,u'shadow'), None): cnv_string,
	((DRAWNS,u'shadow-color'), None): cnv_string,
	((DRAWNS,u'shadow-offset-x'), None): cnv_length,
	((DRAWNS,u'shadow-offset-y'), None): cnv_length,
	((DRAWNS,u'shadow-opacity'), None): cnv_string,
	((DRAWNS,u'shape-id'), None): cnv_IDREF,
	((DRAWNS,u'sharpness'), None): cnv_string,
	((DRAWNS,u'show-unit'), None): cnv_boolean,
	((DRAWNS,u'start-angle'), None): cnv_angle,
	((DRAWNS,u'start'), None): cnv_string,
	((DRAWNS,u'start-color'), None): cnv_string,
	((DRAWNS,u'start-glue-point'), None): cnv_nonNegativeInteger,
	((DRAWNS,u'start-guide'), None): cnv_length,
	((DRAWNS,u'start-intensity'), None): cnv_string,
	((DRAWNS,u'start-line-spacing-horizontal'), None): cnv_string,
	((DRAWNS,u'start-line-spacing-vertical'), None): cnv_string,
	((DRAWNS,u'start-shape'), None): cnv_IDREF,
	((DRAWNS,u'stroke'), None): cnv_string,
	((DRAWNS,u'stroke-dash'), None): cnv_string,
	((DRAWNS,u'stroke-dash-names'), None): cnv_string,
	((DRAWNS,u'stroke-linejoin'), None): cnv_string,
	((DRAWNS,u'style'), None): cnv_string,
	((DRAWNS,u'style-name'), None): cnv_StyleNameRef,
	((DRAWNS,u'symbol-color'), None): cnv_string,
	((DRAWNS,u'text-areas'), None): cnv_string,
	((DRAWNS,u'text-path-allowed'), None): cnv_boolean,
	((DRAWNS,u'text-path'), None): cnv_boolean,
	((DRAWNS,u'text-path-mode'), None): cnv_string,
	((DRAWNS,u'text-path-same-letter-heights'), None): cnv_boolean,
	((DRAWNS,u'text-path-scale'), None): cnv_string,
	((DRAWNS,u'text-rotate-angle'), None): cnv_double,
	((DRAWNS,u'text-style-name'), None): cnv_StyleNameRef,
	((DRAWNS,u'textarea-horizontal-align'), None): cnv_string,
	((DRAWNS,u'textarea-vertical-align'), None): cnv_string,
	((DRAWNS,u'tile-repeat-offset'), None): cnv_string,
	((DRAWNS,u'transform'), None): cnv_string,
	((DRAWNS,u'type'), None): cnv_string,
	((DRAWNS,u'unit'), None): cnv_string,
	((DRAWNS,u'value'), None): cnv_string,
	((DRAWNS,u'visible-area-height'), None): cnv_string,
	((DRAWNS,u'visible-area-left'), None): cnv_string,
	((DRAWNS,u'visible-area-top'), None): cnv_string,
	((DRAWNS,u'visible-area-width'), None): cnv_string,
	((DRAWNS,u'wrap-influence-on-position'), None): cnv_string,
	((DRAWNS,u'z-index'), None): cnv_nonNegativeInteger,
	((FONS,u'background-color'), None): cnv_string,
	((FONS,u'border-bottom'), None): cnv_string,
	((FONS,u'border'), None): cnv_string,
	((FONS,u'border-left'), None): cnv_string,
	((FONS,u'border-right'), None): cnv_string,
	((FONS,u'border-top'), None): cnv_string,
	((FONS,u'break-after'), None): cnv_string,
	((FONS,u'break-before'), None): cnv_string,
	((FONS,u'clip'), None): cnv_string,
	((FONS,u'color'), None): cnv_string,
	((FONS,u'column-count'), None): cnv_positiveInteger,
	((FONS,u'column-gap'), None): cnv_length,
	((FONS,u'country'), None): cnv_token,
	((FONS,u'end-indent'), None): cnv_length,
	((FONS,u'font-family'), None): cnv_string,
	((FONS,u'font-size'), None): cnv_string,
	((FONS,u'font-style'), None): cnv_string,
	((FONS,u'font-variant'), None): cnv_string,
	((FONS,u'font-weight'), None): cnv_string,
	((FONS,u'height'), None): cnv_string,
	((FONS,u'hyphenate'), None): cnv_boolean,
	((FONS,u'hyphenation-keep'), None): cnv_string,
	((FONS,u'hyphenation-ladder-count'), None): cnv_string,
	((FONS,u'hyphenation-push-char-count'), None): cnv_string,
	((FONS,u'hyphenation-remain-char-count'), None): cnv_string,
	((FONS,u'keep-together'), None): cnv_string,
	((FONS,u'keep-with-next'), None): cnv_string,
	((FONS,u'language'), None): cnv_token,
	((FONS,u'letter-spacing'), None): cnv_string,
	((FONS,u'line-height'), None): cnv_string,
	((FONS,u'margin-bottom'), None): cnv_string,
	((FONS,u'margin'), None): cnv_string,
	((FONS,u'margin-left'), None): cnv_string,
	((FONS,u'margin-right'), None): cnv_string,
	((FONS,u'margin-top'), None): cnv_string,
	((FONS,u'max-height'), None): cnv_string,
	((FONS,u'max-width'), None): cnv_string,
	((FONS,u'min-height'), None): cnv_length,
	((FONS,u'min-width'), None): cnv_string,
	((FONS,u'orphans'), None): cnv_string,
	((FONS,u'padding-bottom'), None): cnv_string,
	((FONS,u'padding'), None): cnv_string,
	((FONS,u'padding-left'), None): cnv_string,
	((FONS,u'padding-right'), None): cnv_string,
	((FONS,u'padding-top'), None): cnv_string,
	((FONS,u'page-height'), None): cnv_length,
	((FONS,u'page-width'), None): cnv_length,
	((FONS,u'script'), None): cnv_token,
	((FONS,u'space-after'), None): cnv_length,
	((FONS,u'space-before'), None): cnv_length,
	((FONS,u'start-indent'), None): cnv_length,
	((FONS,u'text-align'), None): cnv_string,
	((FONS,u'text-align-last'), None): cnv_string,
	((FONS,u'text-indent'), None): cnv_string,
	((FONS,u'text-shadow'), None): cnv_string,
	((FONS,u'text-transform'), None): cnv_string,
	((FONS,u'widows'), None): cnv_string,
	((FONS,u'width'), None): cnv_string,
	((FONS,u'wrap-option'), None): cnv_string,
	((FORMNS,u'allow-deletes'), None): cnv_boolean,
	((FORMNS,u'allow-inserts'), None): cnv_boolean,
	((FORMNS,u'allow-updates'), None): cnv_boolean,
	((FORMNS,u'apply-design-mode'), None): cnv_boolean,
	((FORMNS,u'apply-filter'), None): cnv_boolean,
	((FORMNS,u'auto-complete'), None): cnv_boolean,
	((FORMNS,u'automatic-focus'), None): cnv_boolean,
	((FORMNS,u'bound-column'), None): cnv_string,
	((FORMNS,u'button-type'), None): cnv_string,
	((FORMNS,u'command'), None): cnv_string,
	((FORMNS,u'command-type'), None): cnv_string,
	((FORMNS,u'control-implementation'), None): cnv_namespacedToken,
	((FORMNS,u'convert-empty-to-null'), None): cnv_boolean,
	((FORMNS,u'current-selected'), None): cnv_boolean,
	((FORMNS,u'current-state'), None): cnv_string,
#	((FORMNS,u'current-value'), None): cnv_date,
#	((FORMNS,u'current-value'), None): cnv_double,
	((FORMNS,u'current-value'), None): cnv_string,
#	((FORMNS,u'current-value'), None): cnv_time,
	((FORMNS,u'data-field'), None): cnv_string,
	((FORMNS,u'datasource'), None): cnv_string,
	((FORMNS,u'default-button'), None): cnv_boolean,
	((FORMNS,u'delay-for-repeat'), None): cnv_duration,
	((FORMNS,u'detail-fields'), None): cnv_string,
	((FORMNS,u'disabled'), None): cnv_boolean,
	((FORMNS,u'dropdown'), None): cnv_boolean,
	((FORMNS,u'echo-char'), None): cnv_string,
	((FORMNS,u'enctype'), None): cnv_string,
	((FORMNS,u'escape-processing'), None): cnv_boolean,
	((FORMNS,u'filter'), None): cnv_string,
	((FORMNS,u'focus-on-click'), None): cnv_boolean,
	((FORMNS,u'for'), None): cnv_string,
	((FORMNS,u'id'), None): cnv_ID,
	((FORMNS,u'ignore-result'), None): cnv_boolean,
	((FORMNS,u'image-align'), None): cnv_string,
	((FORMNS,u'image-data'), None): cnv_anyURI,
	((FORMNS,u'image-position'), None): cnv_string,
	((FORMNS,u'is-tristate'), None): cnv_boolean,
	((FORMNS,u'label'), None): cnv_string,
	((FORMNS,u'linked-cell'), None): cnv_string,
	((FORMNS,u'list-linkage-type'), None): cnv_list_linkage_type,
	((FORMNS,u'list-source'), None): cnv_string,
	((FORMNS,u'list-source-type'), None): cnv_string,
	((FORMNS,u'master-fields'), None): cnv_string,
	((FORMNS,u'max-length'), None): cnv_nonNegativeInteger,
#	((FORMNS,u'max-value'), None): cnv_date,
#	((FORMNS,u'max-value'), None): cnv_double,
	((FORMNS,u'max-value'), None): cnv_string,
#	((FORMNS,u'max-value'), None): cnv_time,
	((FORMNS,u'method'), None): cnv_string,
#	((FORMNS,u'min-value'), None): cnv_date,
#	((FORMNS,u'min-value'), None): cnv_double,
	((FORMNS,u'min-value'), None): cnv_string,
#	((FORMNS,u'min-value'), None): cnv_time,
	((FORMNS,u'multi-line'), None): cnv_boolean,
	((FORMNS,u'multiple'), None): cnv_boolean,
	((FORMNS,u'name'), None): cnv_string,
	((FORMNS,u'navigation-mode'), None): cnv_string,
	((FORMNS,u'order'), None): cnv_string,
	((FORMNS,u'orientation'), None): cnv_string,
	((FORMNS,u'page-step-size'), None): cnv_positiveInteger,
	((FORMNS,u'printable'), None): cnv_boolean,
	((FORMNS,u'property-name'), None): cnv_string,
	((FORMNS,u'readonly'), None): cnv_boolean,
	((FORMNS,u'repeat'), None): cnv_boolean,
	((FORMNS,u'selected'), None): cnv_boolean,
	((FORMNS,u'size'), None): cnv_nonNegativeInteger,
	((FORMNS,u'source-cell-range'), None): cnv_string,
	((FORMNS,u'spin-button'), None): cnv_boolean,
	((FORMNS,u'state'), None): cnv_string,
	((FORMNS,u'step-size'), None): cnv_positiveInteger,
	((FORMNS,u'tab-cycle'), None): cnv_string,
	((FORMNS,u'tab-index'), None): cnv_nonNegativeInteger,
	((FORMNS,u'tab-stop'), None): cnv_boolean,
	((FORMNS,u'text-style-name'), None): cnv_StyleNameRef,
	((FORMNS,u'title'), None): cnv_string,
	((FORMNS,u'toggle'), None): cnv_boolean,
	((FORMNS,u'validation'), None): cnv_boolean,
#	((FORMNS,u'value'), None): cnv_date,
#	((FORMNS,u'value'), None): cnv_double,
	((FORMNS,u'value'), None): cnv_string,
#	((FORMNS,u'value'), None): cnv_time,
	((FORMNS,u'visual-effect'), None): cnv_string,
	((FORMNS,u'xforms-list-source'), None): cnv_string,
	((FORMNS,u'xforms-submission'), None): cnv_string,
        ((GRDDLNS,u'transformation'), None): cnv_string,
	((LOEXTNS,u'contextual-spacing'), None): cnv_boolean,
        ((LOEXTNS,u'scale-to-X'), None): cnv_string,
        ((LOEXTNS,u'scale-to-Y'), None): cnv_string,
	((MANIFESTNS,u'algorithm-name'), None): cnv_string,
	((MANIFESTNS,u'checksum'), None): cnv_string,
	((MANIFESTNS,u'checksum-type'), None): cnv_string,
	((MANIFESTNS,u'full-path'), None): cnv_string,
	((MANIFESTNS,u'initialisation-vector'), None): cnv_string,
	((MANIFESTNS,u'iteration-count'), None): cnv_nonNegativeInteger,
	((MANIFESTNS,u'key-derivation-name'), None): cnv_string,
	((MANIFESTNS,u'media-type'), None): cnv_string,
	((MANIFESTNS,u'preferred-view-mode'), None): cnv_string,
	((MANIFESTNS,u'salt'), None): cnv_string,
	((MANIFESTNS,u'size'), None): cnv_nonNegativeInteger,
	((MANIFESTNS,u'version'), None): cnv_string,
	((METANS,u'cell-count'), None): cnv_nonNegativeInteger,
	((METANS,u'character-count'), None): cnv_nonNegativeInteger,
	((METANS,u'date'), None): cnv_dateTime,
	((METANS,u'delay'), None): cnv_duration,
	((METANS,u'draw-count'), None): cnv_nonNegativeInteger,
	((METANS,u'frame-count'), None): cnv_nonNegativeInteger,
	((METANS,u'image-count'), None): cnv_nonNegativeInteger,
	((METANS,u'name'), None): cnv_string,
	((METANS,u'non-whitespace-character-count'), None): cnv_nonNegativeInteger,
	((METANS,u'object-count'), None): cnv_nonNegativeInteger,
	((METANS,u'ole-object-count'), None): cnv_nonNegativeInteger,
	((METANS,u'page-count'), None): cnv_nonNegativeInteger,
	((METANS,u'paragraph-count'), None): cnv_nonNegativeInteger,
	((METANS,u'row-count'), None): cnv_nonNegativeInteger,
	((METANS,u'sentence-count'), None): cnv_nonNegativeInteger,
	((METANS,u'syllable-count'), None): cnv_nonNegativeInteger,
	((METANS,u'table-count'), None): cnv_nonNegativeInteger,
	((METANS,u'value-type'), None): cnv_metavaluetype,
	((METANS,u'word-count'), None): cnv_nonNegativeInteger,
	((NUMBERNS,u'automatic-order'), None): cnv_boolean,
	((NUMBERNS,u'calendar'), None): cnv_string,
	((NUMBERNS,u'country'), None): cnv_token,
	((NUMBERNS,u'decimal-places'), None): cnv_integer,
	((NUMBERNS,u'decimal-replacement'), None): cnv_string,
	((NUMBERNS,u'denominator-value'), None): cnv_integer,
	((NUMBERNS,u'display-factor'), None): cnv_double,
	((NUMBERNS,u'format-source'), None): cnv_string,
	((NUMBERNS,u'grouping'), None): cnv_boolean,
	((NUMBERNS,u'language'), None): cnv_token,
	((NUMBERNS,u'min-denominator-digits'), None): cnv_integer,
	((NUMBERNS,u'min-exponent-digits'), None): cnv_integer,
	((NUMBERNS,u'min-integer-digits'), None): cnv_integer,
	((NUMBERNS,u'min-numerator-digits'), None): cnv_integer,
	((NUMBERNS,u'position'), None): cnv_integer,
	((NUMBERNS,u'possessive-form'), None): cnv_boolean,
	((NUMBERNS,u'rfc-language-tag'), None): cnv_language,
	((NUMBERNS,u'script'), None): cnv_token,
	((NUMBERNS,u'style'), None): cnv_string,
	((NUMBERNS,u'textual'), None): cnv_boolean,
	((NUMBERNS,u'title'), None): cnv_string,
	((NUMBERNS,u'transliteration-country'), None): cnv_token,
	((NUMBERNS,u'transliteration-format'), None): cnv_string,
	((NUMBERNS,u'transliteration-language'), None): cnv_token,
	((NUMBERNS,u'transliteration-style'), None): cnv_string,
	((NUMBERNS,u'truncate-on-overflow'), None): cnv_boolean,
	((OFFICENS,u'automatic-update'), None): cnv_boolean,
	((OFFICENS,u'boolean-value'), None): cnv_boolean,
	((OFFICENS,u'conversion-mode'), None): cnv_string,
	((OFFICENS,u'currency'), None): cnv_string,
	((OFFICENS,u'date-value'), None): cnv_dateTime,
	((OFFICENS,u'dde-application'), None): cnv_string,
	((OFFICENS,u'dde-item'), None): cnv_string,
	((OFFICENS,u'dde-topic'), None): cnv_string,
	((OFFICENS,u'display'), None): cnv_boolean,
	((OFFICENS,u'mimetype'), None): cnv_string,
	((OFFICENS,u'name'), None): cnv_string,
	((OFFICENS,u'process-content'), None): cnv_boolean,
	((OFFICENS,u'server-map'), None): cnv_boolean,
	((OFFICENS,u'string-value'), None): cnv_string,
	((OFFICENS,u'target-frame'), None): cnv_string,
	((OFFICENS,u'target-frame-name'), None): cnv_string,
	((OFFICENS,u'time-value'), None): cnv_duration,
	((OFFICENS,u'title'), None): cnv_string,
	((OFFICENS,u'value'), None): cnv_double,
	((OFFICENS,u'value-type'), None): cnv_string,
	((OFFICENS,u'version'), None): cnv_string,
	((PRESENTATIONNS,u'action'), None): cnv_string,
	((PRESENTATIONNS,u'animations'), None): cnv_string,
	((PRESENTATIONNS,u'background-objects-visible'), None): cnv_boolean,
	((PRESENTATIONNS,u'background-visible'), None): cnv_boolean,
	((PRESENTATIONNS,u'class'), None): cnv_string,
	((PRESENTATIONNS,u'class-names'), None): cnv_NCNames,
	((PRESENTATIONNS,u'delay'), None): cnv_duration,
	((PRESENTATIONNS,u'direction'), None): cnv_string,
	((PRESENTATIONNS,u'display-date-time'), None): cnv_boolean,
	((PRESENTATIONNS,u'display-footer'), None): cnv_boolean,
	((PRESENTATIONNS,u'display-header'), None): cnv_boolean,
	((PRESENTATIONNS,u'display-page-number'), None): cnv_boolean,
	((PRESENTATIONNS,u'duration'), None): cnv_string,
	((PRESENTATIONNS,u'effect'), None): cnv_string,
	((PRESENTATIONNS,u'endless'), None): cnv_boolean,
	((PRESENTATIONNS,u'force-manual'), None): cnv_boolean,
	((PRESENTATIONNS,u'full-screen'), None): cnv_boolean,
	((PRESENTATIONNS,u'group-id'), None): cnv_string,
	((PRESENTATIONNS,u'master-element'), None): cnv_IDREF,
	((PRESENTATIONNS,u'mouse-as-pen'), None): cnv_boolean,
	((PRESENTATIONNS,u'mouse-visible'), None): cnv_boolean,
	((PRESENTATIONNS,u'name'), None): cnv_string,
	((PRESENTATIONNS,u'node-type'), None): cnv_string,
	((PRESENTATIONNS,u'object'), None): cnv_string,
	((PRESENTATIONNS,u'pages'), None): cnv_string,
	((PRESENTATIONNS,u'path-id'), None): cnv_string,
	((PRESENTATIONNS,u'pause'), None): cnv_duration,
	((PRESENTATIONNS,u'placeholder'), None): cnv_boolean,
	((PRESENTATIONNS,u'play-full'), None): cnv_boolean,
	((PRESENTATIONNS,u'presentation-page-layout-name'), None): cnv_StyleNameRef,
	((PRESENTATIONNS,u'preset-class'), None): cnv_string,
	((PRESENTATIONNS,u'preset-id'), None): cnv_string,
	((PRESENTATIONNS,u'preset-sub-type'), None): cnv_string,
	((PRESENTATIONNS,u'show'), None): cnv_string,
	((PRESENTATIONNS,u'show-end-of-presentation-slide'), None): cnv_boolean,
	((PRESENTATIONNS,u'show-logo'), None): cnv_boolean,
	((PRESENTATIONNS,u'source'), None): cnv_string,
	((PRESENTATIONNS,u'speed'), None): cnv_string,
	((PRESENTATIONNS,u'start-page'), None): cnv_string,
	((PRESENTATIONNS,u'start-scale'), None): cnv_string,
	((PRESENTATIONNS,u'start-with-navigator'), None): cnv_boolean,
	((PRESENTATIONNS,u'stay-on-top'), None): cnv_boolean,
	((PRESENTATIONNS,u'style-name'), None): cnv_StyleNameRef,
	((PRESENTATIONNS,u'transition-on-click'), None): cnv_string,
	((PRESENTATIONNS,u'transition-speed'), None): cnv_string,
	((PRESENTATIONNS,u'transition-style'), None): cnv_string,
	((PRESENTATIONNS,u'transition-type'), None): cnv_string,
	((PRESENTATIONNS,u'use-date-time-name'), None): cnv_string,
	((PRESENTATIONNS,u'use-footer-name'), None): cnv_string,
	((PRESENTATIONNS,u'use-header-name'), None): cnv_string,
	((PRESENTATIONNS,u'user-transformed'), None): cnv_boolean,
	((PRESENTATIONNS,u'verb'), None): cnv_nonNegativeInteger,
	((PRESENTATIONNS,u'visibility'), None): cnv_string,
	((SCRIPTNS,u'event-name'), None): cnv_formula,
	((SCRIPTNS,u'language'), None): cnv_formula,
	((SCRIPTNS,u'macro-name'), None): cnv_string,
	((SMILNS,u'accelerate'), None): cnv_double,
	((SMILNS,u'accumulate'), None): cnv_string,
	((SMILNS,u'additive'), None): cnv_string,
	((SMILNS,u'attributeName'), None): cnv_string,
	((SMILNS,u'autoReverse'), None): cnv_boolean,
	((SMILNS,u'begin'), None): cnv_string,
	((SMILNS,u'by'), None): cnv_string,
	((SMILNS,u'calcMode'), None): cnv_string,
	((SMILNS,u'decelerate'), None): cnv_double,
	((SMILNS,u'direction'), None): cnv_string,
	((SMILNS,u'dur'), None): cnv_string,
	((SMILNS,u'end'), None): cnv_string,
	((SMILNS,u'endsync'), None): cnv_string,
	((SMILNS,u'fadeColor'), None): cnv_string,
	((SMILNS,u'fill'), None): cnv_string,
	((SMILNS,u'fillDefault'), None): cnv_string,
	((SMILNS,u'from'), None): cnv_string,
	((SMILNS,u'keySplines'), None): cnv_string,
	((SMILNS,u'keyTimes'), None): cnv_string,
	((SMILNS,u'mode'), None): cnv_string,
	((SMILNS,u'repeatCount'), None): cnv_nonNegativeInteger,
	((SMILNS,u'repeatDur'), None): cnv_string,
	((SMILNS,u'restart'), None): cnv_string,
	((SMILNS,u'restartDefault'), None): cnv_string,
	((SMILNS,u'subtype'), None): cnv_string,
	((SMILNS,u'targetElement'), None): cnv_IDREF,
	((SMILNS,u'to'), None): cnv_string,
	((SMILNS,u'type'), None): cnv_string,
	((SMILNS,u'values'), None): cnv_string,
	((STYLENS,u'adjustment'), None): cnv_string,
	((STYLENS,u'apply-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'auto-text-indent'), None): cnv_boolean,
	((STYLENS,u'auto-update'), None): cnv_boolean,
	((STYLENS,u'background-transparency'), None): cnv_string,
	((STYLENS,u'base-cell-address'), None): cnv_string,
	((STYLENS,u'border-line-width-bottom'), None): cnv_string,
	((STYLENS,u'border-line-width'), None): cnv_string,
	((STYLENS,u'border-line-width-left'), None): cnv_string,
	((STYLENS,u'border-line-width-right'), None): cnv_string,
	((STYLENS,u'border-line-width-top'), None): cnv_string,
	((STYLENS,u'cell-protect'), None): cnv_string,
	((STYLENS,u'char'), None): cnv_string,
	((STYLENS,u'class'), None): cnv_string,
	((STYLENS,u'color'), None): cnv_string,
	((STYLENS,u'column-width'), None): cnv_string,
	((STYLENS,u'condition'), None): cnv_string,
	((STYLENS,u'country-asian'), None): cnv_string,
	((STYLENS,u'country-complex'), None): cnv_string,
	((STYLENS,u'data-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'decimal-places'), None): cnv_string,
	((STYLENS,u'default-outline-level'), None): cnv_positiveInteger,
	((STYLENS,u'diagonal-bl-tr'), None): cnv_string,
	((STYLENS,u'diagonal-bl-tr-widths'), None): cnv_string,
	((STYLENS,u'diagonal-tl-br'), None): cnv_string,
	((STYLENS,u'diagonal-tl-br-widths'), None): cnv_string,
	((STYLENS,u'direction'), None): cnv_string,
	((STYLENS,u'display'), None): cnv_boolean,
	((STYLENS,u'display-name'), None): cnv_string,
	((STYLENS,u'distance-after-sep'), None): cnv_length,
	((STYLENS,u'distance-before-sep'), None): cnv_length,
	((STYLENS,u'distance'), None): cnv_length,
	((STYLENS,u'dynamic-spacing'), None): cnv_boolean,
	((STYLENS,u'editable'), None): cnv_boolean,
	((STYLENS,u'family'), None): cnv_family,
	((STYLENS,u'filter-name'), None): cnv_string,
	((STYLENS,u'first-page-number'), None): cnv_string,
	((STYLENS,u'flow-with-text'), None): cnv_boolean,
	((STYLENS,u'font-adornments'), None): cnv_string,
	((STYLENS,u'font-charset'), None): cnv_string,
	((STYLENS,u'font-charset-asian'), None): cnv_string,
	((STYLENS,u'font-charset-complex'), None): cnv_string,
	((STYLENS,u'font-family-asian'), None): cnv_string,
	((STYLENS,u'font-family-complex'), None): cnv_string,
	((STYLENS,u'font-family-generic-asian'), None): cnv_string,
	((STYLENS,u'font-family-generic'), None): cnv_string,
	((STYLENS,u'font-family-generic-complex'), None): cnv_string,
	((STYLENS,u'font-independent-line-spacing'), None): cnv_boolean,
	((STYLENS,u'font-name-asian'), None): cnv_string,
	((STYLENS,u'font-name'), None): cnv_string,
	((STYLENS,u'font-name-complex'), None): cnv_string,
	((STYLENS,u'font-pitch-asian'), None): cnv_string,
	((STYLENS,u'font-pitch'), None): cnv_string,
	((STYLENS,u'font-pitch-complex'), None): cnv_string,
	((STYLENS,u'font-relief'), None): cnv_string,
	((STYLENS,u'font-size-asian'), None): cnv_string,
	((STYLENS,u'font-size-complex'), None): cnv_string,
	((STYLENS,u'font-size-rel-asian'), None): cnv_length,
	((STYLENS,u'font-size-rel'), None): cnv_length,
	((STYLENS,u'font-size-rel-complex'), None): cnv_length,
	((STYLENS,u'font-style-asian'), None): cnv_string,
	((STYLENS,u'font-style-complex'), None): cnv_string,
	((STYLENS,u'font-style-name-asian'), None): cnv_string,
	((STYLENS,u'font-style-name'), None): cnv_string,
	((STYLENS,u'font-style-name-complex'), None): cnv_string,
	((STYLENS,u'font-weight-asian'), None): cnv_string,
	((STYLENS,u'font-weight-complex'), None): cnv_string,
	((STYLENS,u'footnote-max-height'), None): cnv_length,
	((STYLENS,u'glyph-orientation-vertical'), None): cnv_string,
	((STYLENS,u'height'), None): cnv_string,
	((STYLENS,u'horizontal-pos'), None): cnv_string,
	((STYLENS,u'horizontal-rel'), None): cnv_string,
	((STYLENS,u'join-border'), None): cnv_boolean,
	((STYLENS,u'justify-single-word'), None): cnv_boolean,
	((STYLENS,u'language-asian'), None): cnv_string,
	((STYLENS,u'language-complex'), None): cnv_string,
	((STYLENS,u'layout-grid-base-height'), None): cnv_length,
	((STYLENS,u'layout-grid-base-width'), None): cnv_length,
	((STYLENS,u'layout-grid-color'), None): cnv_string,
	((STYLENS,u'layout-grid-display'), None): cnv_boolean,
	((STYLENS,u'layout-grid-lines'), None): cnv_string,
	((STYLENS,u'layout-grid-mode'), None): cnv_string,
	((STYLENS,u'layout-grid-print'), None): cnv_boolean,
	((STYLENS,u'layout-grid-ruby-below'), None): cnv_boolean,
	((STYLENS,u'layout-grid-ruby-height'), None): cnv_length,
	((STYLENS,u'layout-grid-snap-to'), None): cnv_boolean,
	((STYLENS,u'layout-grid-standard-mode'), None): cnv_boolean,
	((STYLENS,u'leader-char'), None): cnv_string,
	((STYLENS,u'leader-color'), None): cnv_string,
	((STYLENS,u'leader-style'), None): cnv_string,
	((STYLENS,u'leader-text'), None): cnv_string,
	((STYLENS,u'leader-text-style'), None): cnv_StyleNameRef,
	((STYLENS,u'leader-type'), None): cnv_string,
	((STYLENS,u'leader-width'), None): cnv_string,
	((STYLENS,u'legend-expansion-aspect-ratio'), None): cnv_double,
	((STYLENS,u'legend-expansion'), None): cnv_string,
	((STYLENS,u'length'), None): cnv_positiveInteger,
	((STYLENS,u'letter-kerning'), None): cnv_boolean,
	((STYLENS,u'line-break'), None): cnv_string,
	((STYLENS,u'line-height-at-least'), None): cnv_string,
	((STYLENS,u'line-spacing'), None): cnv_length,
	((STYLENS,u'line-style'), None): cnv_string,
	((STYLENS,u'lines'), None): cnv_positiveInteger,
	((STYLENS,u'list-level'), None): cnv_positiveInteger,
	((STYLENS,u'list-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'master-page-name'), None): cnv_StyleNameRef,
	((STYLENS,u'may-break-between-rows'), None): cnv_boolean,
	((STYLENS,u'min-row-height'), None): cnv_string,
	((STYLENS,u'mirror'), None): cnv_string,
	((STYLENS,u'name'), None): cnv_NCName,
 	((STYLENS,u'name'), (STYLENS,u'font-face')): cnv_string,
	((STYLENS,u'next-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'num-format'), None): cnv_string,
	((STYLENS,u'num-letter-sync'), None): cnv_boolean,
	((STYLENS,u'num-prefix'), None): cnv_string,
	((STYLENS,u'num-suffix'), None): cnv_string,
	((STYLENS,u'number-wrapped-paragraphs'), None): cnv_string,
	((STYLENS,u'overflow-behavior'), None): cnv_string,
	((STYLENS,u'page-layout-name'), None): cnv_StyleNameRef,
	((STYLENS,u'page-number'), None): cnv_string,
	((STYLENS,u'page-usage'), None): cnv_string,
	((STYLENS,u'paper-tray-name'), None): cnv_string,
	((STYLENS,u'parent-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'percentage-data-style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'position'), (STYLENS,u'tab-stop')): cnv_length,
	((STYLENS,u'position'), None): cnv_string,
	((STYLENS,u'print'), None): cnv_string,
	((STYLENS,u'print-content'), None): cnv_boolean,
	((STYLENS,u'print-orientation'), None): cnv_string,
	((STYLENS,u'print-page-order'), None): cnv_string,
	((STYLENS,u'protect'), (STYLENS,u'section-properties')): cnv_boolean,
	((STYLENS,u'protect'), (STYLENS,u'graphic-properties')): cnv_string,
#	((STYLENS,u'protect'), None): cnv_boolean,
	((STYLENS,u'punctuation-wrap'), None): cnv_string,
	((STYLENS,u'register-true'), None): cnv_boolean,
	((STYLENS,u'register-truth-ref-style-name'), None): cnv_string,
	((STYLENS,u'rel-column-width'), None): cnv_string,
	((STYLENS,u'rel-height'), None): cnv_string,
	((STYLENS,u'rel-width'), None): cnv_string,
	((STYLENS,u'repeat'), None): cnv_string,
	((STYLENS,u'repeat-content'), None): cnv_boolean,
	((STYLENS,u'rfc-language-tag'), None): cnv_language,
	((STYLENS,u'rfc-language-tag-asian'), None): cnv_language,
	((STYLENS,u'rfc-language-tag-complex'), None): cnv_language,
	((STYLENS,u'rotation-align'), None): cnv_string,
	((STYLENS,u'rotation-angle'), None): cnv_string,
	((STYLENS,u'row-height'), None): cnv_string,
	((STYLENS,u'ruby-align'), None): cnv_string,
	((STYLENS,u'ruby-position'), None): cnv_string,
	((STYLENS,u'run-through'), None): cnv_string,
	((STYLENS,u'scale-to'), None): cnv_string,
	((STYLENS,u'scale-to-pages'), None): cnv_string,
	((STYLENS,u'script-asian'), None): cnv_string,
	((STYLENS,u'script-complex'), None): cnv_string,
	((STYLENS,u'script-type'), None): cnv_string,
	((STYLENS,u'shadow'), None): cnv_string,
	((STYLENS,u'shrink-to-fit'), None): cnv_boolean,
	((STYLENS,u'snap-to-layout-grid'), None): cnv_boolean,
	((STYLENS,u'style'), None): cnv_string,
	((STYLENS,u'style-name'), None): cnv_StyleNameRef,
	((STYLENS,u'tab-stop-distance'), None): cnv_string,
	((STYLENS,u'table-centering'), None): cnv_string,
	((STYLENS,u'text-align-source'), None): cnv_string,
	((STYLENS,u'text-autospace'), None): cnv_string,
	((STYLENS,u'text-blinking'), None): cnv_boolean,
	((STYLENS,u'text-combine'), None): cnv_string,
	((STYLENS,u'text-combine-end-char'), None): cnv_string,
	((STYLENS,u'text-combine-start-char'), None): cnv_string,
	((STYLENS,u'text-emphasize'), None): cnv_string,
	((STYLENS,u'text-line-through-color'), None): cnv_string,
	((STYLENS,u'text-line-through-mode'), None): cnv_string,
	((STYLENS,u'text-line-through-style'), None): cnv_string,
	((STYLENS,u'text-line-through-text'), None): cnv_string,
	((STYLENS,u'text-line-through-text-style'), None): cnv_string,
	((STYLENS,u'text-line-through-type'), None): cnv_string,
	((STYLENS,u'text-line-through-width'), None): cnv_string,
	((STYLENS,u'text-outline'), None): cnv_boolean,
	((STYLENS,u'text-overline-color'), None): cnv_string,
	((STYLENS,u'text-overline-mode'), None): cnv_string,
	((STYLENS,u'text-overline-style'), None): cnv_string,
	((STYLENS,u'text-overline-type'), None): cnv_string,
	((STYLENS,u'text-overline-width'), None): cnv_string,
	((STYLENS,u'text-position'), None): cnv_string,
	((STYLENS,u'text-rotation-angle'), None): cnv_string,
	((STYLENS,u'text-rotation-scale'), None): cnv_string,
	((STYLENS,u'text-scale'), None): cnv_string,
	((STYLENS,u'text-underline-color'), None): cnv_string,
	((STYLENS,u'text-underline-mode'), None): cnv_string,
	((STYLENS,u'text-underline-style'), None): cnv_string,
	((STYLENS,u'text-underline-type'), None): cnv_string,
	((STYLENS,u'text-underline-width'), None): cnv_string,
	((STYLENS,u'type'), None): cnv_string,
	((STYLENS,u'use-optimal-column-width'), None): cnv_boolean,
	((STYLENS,u'use-optimal-row-height'), None): cnv_boolean,
	((STYLENS,u'use-window-font-color'), None): cnv_boolean,
	((STYLENS,u'vertical-align'), None): cnv_string,
	((STYLENS,u'vertical-pos'), None): cnv_string,
	((STYLENS,u'vertical-rel'), None): cnv_string,
	((STYLENS,u'volatile'), None): cnv_boolean,
	((STYLENS,u'width'), None): cnv_string,
	((STYLENS,u'wrap'), None): cnv_string,
	((STYLENS,u'wrap-contour'), None): cnv_boolean,
	((STYLENS,u'wrap-contour-mode'), None): cnv_string,
	((STYLENS,u'wrap-dynamic-threshold'), None): cnv_length,
	((STYLENS,u'writing-mode-automatic'), None): cnv_boolean,
	((STYLENS,u'writing-mode'), None): cnv_string,
	((SVGNS,u'accent-height'), None): cnv_integer,
	((SVGNS,u'alphabetic'), None): cnv_integer,
	((SVGNS,u'ascent'), None): cnv_integer,
	((SVGNS,u'bbox'), None): cnv_string,
	((SVGNS,u'cap-height'), None): cnv_integer,
	((SVGNS,u'cx'), None): cnv_string,
	((SVGNS,u'cy'), None): cnv_string,
	((SVGNS,u'd'), None): cnv_string,
	((SVGNS,u'descent'), None): cnv_integer,
	((SVGNS,u'fill-rule'), None): cnv_string,
	((SVGNS,u'font-family'), None): cnv_string,
	((SVGNS,u'font-size'), None): cnv_string,
	((SVGNS,u'font-stretch'), None): cnv_string,
	((SVGNS,u'font-style'), None): cnv_string,
	((SVGNS,u'font-variant'), None): cnv_string,
	((SVGNS,u'font-weight'), None): cnv_string,
	((SVGNS,u'fx'), None): cnv_string,
	((SVGNS,u'fy'), None): cnv_string,
	((SVGNS,u'gradientTransform'), None): cnv_string,
	((SVGNS,u'gradientUnits'), None): cnv_string,
	((SVGNS,u'hanging'), None): cnv_integer,
	((SVGNS,u'height'), None): cnv_length,
	((SVGNS,u'ideographic'), None): cnv_integer,
	((SVGNS,u'mathematical'), None): cnv_integer,
	((SVGNS,u'name'), None): cnv_string,
	((SVGNS,u'offset'), None): cnv_string,
	((SVGNS,u'origin'), None): cnv_string,
	((SVGNS,u'overline-position'), None): cnv_integer,
	((SVGNS,u'overline-thickness'), None): cnv_integer,
	((SVGNS,u'panose-1'), None): cnv_string,
	((SVGNS,u'path'), None): cnv_string,
	((SVGNS,u'r'), None): cnv_length,
	((SVGNS,u'rx'), None): cnv_length,
	((SVGNS,u'ry'), None): cnv_length,
	((SVGNS,u'slope'), None): cnv_integer,
	((SVGNS,u'spreadMethod'), None): cnv_string,
	((SVGNS,u'stemh'), None): cnv_integer,
	((SVGNS,u'stemv'), None): cnv_integer,
	((SVGNS,u'stop-color'), None): cnv_string,
	((SVGNS,u'stop-opacity'), None): cnv_double,
	((SVGNS,u'strikethrough-position'), None): cnv_integer,
	((SVGNS,u'strikethrough-thickness'), None): cnv_integer,
	((SVGNS,u'string'), None): cnv_string,
	((SVGNS,u'stroke-color'), None): cnv_string,
	((SVGNS,u'stroke-linecap'), None): cnv_stroke_linecap,
	((SVGNS,u'stroke-opacity'), None): cnv_string,
	((SVGNS,u'stroke-width'), None): cnv_length,
	((SVGNS,u'type'), None): cnv_string,
	((SVGNS,u'underline-position'), None): cnv_integer,
	((SVGNS,u'underline-thickness'), None): cnv_integer,
	((SVGNS,u'unicode-range'), None): cnv_string,
	((SVGNS,u'units-per-em'), None): cnv_integer,
	((SVGNS,u'v-alphabetic'), None): cnv_integer,
	((SVGNS,u'v-hanging'), None): cnv_integer,
	((SVGNS,u'v-ideographic'), None): cnv_integer,
	((SVGNS,u'v-mathematical'), None): cnv_integer,
	((SVGNS,u'viewBox'), None): cnv_viewbox,
	((SVGNS,u'width'), None): cnv_length,
	((SVGNS,u'widths'), None): cnv_string,
	((SVGNS,u'x'), None): cnv_length,
	((SVGNS,u'x-height'), None): cnv_integer,
	((SVGNS,u'x1'), None): cnv_lengthorpercent,
	((SVGNS,u'x2'), None): cnv_lengthorpercent,
	((SVGNS,u'y'), None): cnv_length,
	((SVGNS,u'y1'), None): cnv_lengthorpercent,
	((SVGNS,u'y2'), None): cnv_lengthorpercent,
	((TABLENS,u'acceptance-state'), None): cnv_string,
	((TABLENS,u'add-empty-lines'), None): cnv_boolean,
	((TABLENS,u'algorithm'), None): cnv_formula,
	((TABLENS,u'align'), None): cnv_string,
	((TABLENS,u'allow-empty-cell'), None): cnv_boolean,
	((TABLENS,u'application-data'), None): cnv_string,
	((TABLENS,u'automatic-find-labels'), None): cnv_boolean,
	((TABLENS,u'base-cell-address'), None): cnv_string,
	((TABLENS,u'bind-styles-to-content'), None): cnv_boolean,
	((TABLENS,u'border-color'), None): cnv_string,
	((TABLENS,u'border-model'), None): cnv_string,
	((TABLENS,u'buttons'), None): cnv_string,
	((TABLENS,u'buttons'), None): cnv_string,
	((TABLENS,u'case-sensitive'), None): cnv_boolean,
	((TABLENS,u'case-sensitive'), None): cnv_string,
	((TABLENS,u'cell-address'), None): cnv_string,
	((TABLENS,u'cell-range-address'), None): cnv_string,
	((TABLENS,u'cell-range-address'), None): cnv_string,
	((TABLENS,u'cell-range'), None): cnv_string,
	((TABLENS,u'column'), None): cnv_integer,
	((TABLENS,u'comment'), None): cnv_string,
	((TABLENS,u'condition'), None): cnv_formula,
	((TABLENS,u'condition-source'), None): cnv_string,
	((TABLENS,u'condition-source-range-address'), None): cnv_string,
	((TABLENS,u'contains-error'), None): cnv_boolean,
	((TABLENS,u'contains-header'), None): cnv_boolean,
	((TABLENS,u'content-validation-name'), None): cnv_string,
	((TABLENS,u'copy-back'), None): cnv_boolean,
	((TABLENS,u'copy-formulas'), None): cnv_boolean,
	((TABLENS,u'copy-styles'), None): cnv_boolean,
	((TABLENS,u'count'), None): cnv_positiveInteger,
	((TABLENS,u'country'), None): cnv_token,
	((TABLENS,u'data-cell-range-address'), None): cnv_string,
	((TABLENS,u'data-field'), None): cnv_string,
	((TABLENS,u'data-type'), None): cnv_string,
	((TABLENS,u'database-name'), None): cnv_string,
	((TABLENS,u'database-table-name'), None): cnv_string,
	((TABLENS,u'date-end'), None): cnv_string,
	((TABLENS,u'date-start'), None): cnv_string,
	((TABLENS,u'date-value'), None): cnv_date,
	((TABLENS,u'default-cell-style-name'), None): cnv_StyleNameRef,
	((TABLENS,u'direction'), None): cnv_string,
	((TABLENS,u'display-border'), None): cnv_boolean,
	((TABLENS,u'display'), None): cnv_boolean,
	((TABLENS,u'display-duplicates'), None): cnv_boolean,
	((TABLENS,u'display-filter-buttons'), None): cnv_boolean,
	((TABLENS,u'display-list'), None): cnv_string,
	((TABLENS,u'display-member-mode'), None): cnv_string,
	((TABLENS,u'drill-down-on-double-click'), None): cnv_boolean,
	((TABLENS,u'embedded-number-behavior'), None): cnv_string,
	((TABLENS,u'enabled'), None): cnv_boolean,
	((TABLENS,u'end-cell-address'), None): cnv_string,
	((TABLENS,u'end'), None): cnv_string,
	((TABLENS,u'end-column'), None): cnv_integer,
	((TABLENS,u'end-position'), None): cnv_integer,
	((TABLENS,u'end-row'), None): cnv_integer,
	((TABLENS,u'end-table'), None): cnv_integer,
	((TABLENS,u'end-x'), None): cnv_length,
	((TABLENS,u'end-y'), None): cnv_length,
	((TABLENS,u'execute'), None): cnv_boolean,
	((TABLENS,u'expression'), None): cnv_formula,
	((TABLENS,u'field-name'), None): cnv_string,
	((TABLENS,u'field-number'), None): cnv_nonNegativeInteger,
	((TABLENS,u'field-number'), None): cnv_string,
	((TABLENS,u'filter-name'), None): cnv_string,
	((TABLENS,u'filter-options'), None): cnv_string,
	((TABLENS,u'first-row-end-column'), None): cnv_rowOrCol,
	((TABLENS,u'first-row-start-column'), None): cnv_rowOrCol,
	((TABLENS,u'formula'), None): cnv_formula,
	((TABLENS,u'function'), None): cnv_string,
	((TABLENS,u'function'), None): cnv_string,
	((TABLENS,u'grand-total'), None): cnv_string,
	((TABLENS,u'group-by-field-number'), None): cnv_nonNegativeInteger,
	((TABLENS,u'grouped-by'), None): cnv_string,
	((TABLENS,u'has-persistent-data'), None): cnv_boolean,
	((TABLENS,u'id'), None): cnv_string,
	((TABLENS,u'identify-categories'), None): cnv_boolean,
	((TABLENS,u'ignore-empty-rows'), None): cnv_boolean,
	((TABLENS,u'index'), None): cnv_nonNegativeInteger,
	((TABLENS,u'is-active'), None): cnv_boolean,
	((TABLENS,u'is-data-layout-field'), None): cnv_string,
	((TABLENS,u'is-selection'), None): cnv_boolean,
	((TABLENS,u'is-sub-table'), None): cnv_boolean,
	((TABLENS,u'label-cell-range-address'), None): cnv_string,
	((TABLENS,u'language'), None): cnv_token,
	((TABLENS,u'language'), None): cnv_token,
	((TABLENS,u'last-column-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'last-row-end-column'), None): cnv_rowOrCol,
	((TABLENS,u'last-row-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'last-row-start-column'), None): cnv_rowOrCol,
	((TABLENS,u'layout-mode'), None): cnv_string,
	((TABLENS,u'link-to-source-data'), None): cnv_boolean,
	((TABLENS,u'marked-invalid'), None): cnv_boolean,
	((TABLENS,u'matrix-covered'), None): cnv_boolean,
	((TABLENS,u'maximum-difference'), None): cnv_double,
	((TABLENS,u'member-count'), None): cnv_nonNegativeInteger,
	((TABLENS,u'member-name'), None): cnv_string,
	((TABLENS,u'member-type'), None): cnv_string,
	((TABLENS,u'message-type'), None): cnv_string,
	((TABLENS,u'mode'), None): cnv_string,
	((TABLENS,u'multi-deletion-spanned'), None): cnv_integer,
	((TABLENS,u'name'), None): cnv_string,
	((TABLENS,u'name'), None): cnv_string,
	((TABLENS,u'null-year'), None): cnv_positiveInteger,
	((TABLENS,u'number-columns-repeated'), None): cnv_positiveInteger,
	((TABLENS,u'number-columns-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'number-matrix-columns-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'number-matrix-rows-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'number-rows-repeated'), None): cnv_positiveInteger,
	((TABLENS,u'number-rows-spanned'), None): cnv_positiveInteger,
	((TABLENS,u'object-name'), None): cnv_string,
	((TABLENS,u'on-update-keep-size'), None): cnv_boolean,
	((TABLENS,u'on-update-keep-styles'), None): cnv_boolean,
	((TABLENS,u'operator'), None): cnv_string,
	((TABLENS,u'operator'), None): cnv_string,
	((TABLENS,u'order'), None): cnv_string,
	((TABLENS,u'orientation'), None): cnv_string,
	((TABLENS,u'orientation'), None): cnv_string,
	((TABLENS,u'page-breaks-on-group-change'), None): cnv_boolean,
	((TABLENS,u'paragraph-style-name'), None): cnv_StyleNameRef,
	((TABLENS,u'parse-sql-statement'), None): cnv_boolean,
	((TABLENS,u'password'), None): cnv_string,
	((TABLENS,u'position'), None): cnv_integer,
	((TABLENS,u'precision-as-shown'), None): cnv_boolean,
	((TABLENS,u'print'), None): cnv_boolean,
	((TABLENS,u'print-ranges'), None): cnv_string,
	((TABLENS,u'protect'), None): cnv_boolean,
	((TABLENS,u'protected'), None): cnv_boolean,
	((TABLENS,u'protection-key'), None): cnv_string,
	((TABLENS,u'protection-key-digest-algorithm'), None): cnv_anyURI,
	((TABLENS,u'query-name'), None): cnv_string,
	((TABLENS,u'range-usable-as'), None): cnv_string,
	((TABLENS,u'rfc-language-tag'), None): cnv_language,
	((TABLENS,u'refresh-delay'), None): cnv_boolean,
	((TABLENS,u'refresh-delay'), None): cnv_duration,
	((TABLENS,u'rejecting-change-id'), None): cnv_string,
	((TABLENS,u'row'), None): cnv_integer,
	((TABLENS,u'scenario-ranges'), None): cnv_string,
	((TABLENS,u'script'), None): cnv_string,
	((TABLENS,u'search-criteria-must-apply-to-whole-cell'), None): cnv_boolean,
	((TABLENS,u'selected-page'), None): cnv_string,
	((TABLENS,u'show-details'), None): cnv_boolean,
	((TABLENS,u'show-empty'), None): cnv_boolean,
	((TABLENS,u'show-empty'), None): cnv_string,
	((TABLENS,u'show-filter-button'), None): cnv_boolean,
	((TABLENS,u'sort-mode'), None): cnv_string,
	((TABLENS,u'source-cell-range-addresses'), None): cnv_string,
	((TABLENS,u'source-cell-range-addresses'), None): cnv_string,
	((TABLENS,u'source-field-name'), None): cnv_string,
	((TABLENS,u'source-field-name'), None): cnv_string,
	((TABLENS,u'source-name'), None): cnv_string,
	((TABLENS,u'sql-statement'), None): cnv_string,
	((TABLENS,u'start'), None): cnv_string,
	((TABLENS,u'start-column'), None): cnv_integer,
	((TABLENS,u'start-position'), None): cnv_integer,
	((TABLENS,u'start-row'), None): cnv_integer,
	((TABLENS,u'start-table'), None): cnv_integer,
	((TABLENS,u'status'), None): cnv_string,
	((TABLENS,u'step'), None): cnv_double,
	((TABLENS,u'steps'), None): cnv_positiveInteger,
	((TABLENS,u'structure-protected'), None): cnv_boolean,
	((TABLENS,u'style-name'), None): cnv_StyleNameRef,
	((TABLENS,u'table-background'), None): cnv_boolean,
	((TABLENS,u'table'), None): cnv_integer,
	((TABLENS,u'table-name'), None): cnv_string,
	((TABLENS,u'target-cell-address'), None): cnv_string,
	((TABLENS,u'target-cell-address'), None): cnv_string,
	((TABLENS,u'target-range-address'), None): cnv_string,
	((TABLENS,u'target-range-address'), None): cnv_string,
	((TABLENS,u'template-name'), None): cnv_string,
	((TABLENS,u'title'), None): cnv_string,
	((TABLENS,u'track-changes'), None): cnv_boolean,
	((TABLENS,u'type'), None): cnv_string,
	((TABLENS,u'use-banding-columns-styles'), None): cnv_boolean,
	((TABLENS,u'use-banding-rows-styles'), None): cnv_boolean,
	((TABLENS,u'use-first-column-styles'), None): cnv_boolean,
	((TABLENS,u'use-first-row-styles'), None): cnv_boolean,
	((TABLENS,u'use-labels'), None): cnv_string,
	((TABLENS,u'use-last-column-styles'), None): cnv_boolean,
	((TABLENS,u'use-last-row-styles'), None): cnv_boolean,
	((TABLENS,u'use-regular-expressions'), None): cnv_boolean,
	((TABLENS,u'use-wildcards'), None): cnv_boolean,
	((TABLENS,u'used-hierarchy'), None): cnv_integer,
	((TABLENS,u'user-name'), None): cnv_string,
	((TABLENS,u'value'), None): cnv_string,
	((TABLENS,u'value'), None): cnv_string,
	((TABLENS,u'value-type'), None): cnv_string,
	((TABLENS,u'visibility'), None): cnv_string,
	((TEXTNS,u'active'), None): cnv_boolean,
	((TEXTNS,u'address'), None): cnv_string,
	((TEXTNS,u'alphabetical-separators'), None): cnv_boolean,
	((TEXTNS,u'anchor-page-number'), None): cnv_positiveInteger,
	((TEXTNS,u'anchor-type'), None): cnv_string,
	((TEXTNS,u'animation'), None): cnv_string,
	((TEXTNS,u'animation-delay'), None): cnv_string,
	((TEXTNS,u'animation-direction'), None): cnv_string,
	((TEXTNS,u'animation-repeat'), None): cnv_string,
	((TEXTNS,u'animation-start-inside'), None): cnv_boolean,
	((TEXTNS,u'animation-steps'), None): cnv_length,
	((TEXTNS,u'animation-stop-inside'), None): cnv_boolean,
	((TEXTNS,u'annote'), None): cnv_string,
	((TEXTNS,u'author'), None): cnv_string,
	((TEXTNS,u'bibliography-data-field'), None): cnv_string,
	((TEXTNS,u'bibliography-type'), None): cnv_string,
	((TEXTNS,u'booktitle'), None): cnv_string,
	((TEXTNS,u'bullet-char'), None): cnv_string,
	((TEXTNS,u'bullet-relative-size'), None): cnv_string,
	((TEXTNS,u'c'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'capitalize-entries'), None): cnv_boolean,
	((TEXTNS,u'caption-sequence-format'), None): cnv_string,
	((TEXTNS,u'caption-sequence-name'), None): cnv_string,
	((TEXTNS,u'change-id'), None): cnv_IDREF,
	((TEXTNS,u'chapter'), None): cnv_string,
	((TEXTNS,u'citation-body-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'citation-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'class-names'), None): cnv_NCNames,
	((TEXTNS,u'column-name'), None): cnv_string,
	((TEXTNS,u'combine-entries'), None): cnv_boolean,
	((TEXTNS,u'combine-entries-with-dash'), None): cnv_boolean,
	((TEXTNS,u'combine-entries-with-pp'), None): cnv_boolean,
	((TEXTNS,u'comma-separated'), None): cnv_boolean,
	((TEXTNS,u'cond-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'condition'), None): cnv_formula,
	((TEXTNS,u'connection-name'), None): cnv_string,
	((TEXTNS,u'consecutive-numbering'), None): cnv_boolean,
	((TEXTNS,u'continue-list'), None): cnv_IDREF,
	((TEXTNS,u'continue-numbering'), None): cnv_boolean,
	((TEXTNS,u'copy-outline-levels'), None): cnv_boolean,
	((TEXTNS,u'count-empty-lines'), None): cnv_boolean,
	((TEXTNS,u'count-in-text-boxes'), None): cnv_boolean,
	((TEXTNS,u'current-value'), None): cnv_boolean,
	((TEXTNS,u'custom1'), None): cnv_string,
	((TEXTNS,u'custom2'), None): cnv_string,
	((TEXTNS,u'custom3'), None): cnv_string,
	((TEXTNS,u'custom4'), None): cnv_string,
	((TEXTNS,u'custom5'), None): cnv_string,
	((TEXTNS,u'database-name'), None): cnv_string,
	((TEXTNS,u'date-adjust'), None): cnv_duration,
	((TEXTNS,u'date-value'), None): cnv_date,
#	((TEXTNS,u'date-value'), None): cnv_dateTime,
	((TEXTNS,u'default-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'description'), None): cnv_string,
	((TEXTNS,u'display'), None): cnv_string,
	((TEXTNS,u'display-levels'), None): cnv_positiveInteger,
	((TEXTNS,u'display-outline-level'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'dont-balance-text-columns'), None): cnv_boolean,
	((TEXTNS,u'duration'), None): cnv_duration,
	((TEXTNS,u'edition'), None): cnv_string,
	((TEXTNS,u'editor'), None): cnv_string,
	((TEXTNS,u'filter-name'), None): cnv_string,
	((TEXTNS,u'fixed'), None): cnv_boolean,
	((TEXTNS,u'footnotes-position'), None): cnv_string,
	((TEXTNS,u'formula'), None): cnv_formula,
	((TEXTNS,u'global'), None): cnv_boolean,
	((TEXTNS,u'howpublished'), None): cnv_string,
	((TEXTNS,u'id'), None): cnv_ID,
#	((TEXTNS,u'id'), None): cnv_string,
	((TEXTNS,u'identifier'), None): cnv_string,
	((TEXTNS,u'ignore-case'), None): cnv_boolean,
	((TEXTNS,u'increment'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'index-name'), None): cnv_string,
	((TEXTNS,u'index-scope'), None): cnv_string,
	((TEXTNS,u'institution'), None): cnv_string,
	((TEXTNS,u'is-hidden'), None): cnv_boolean,
	((TEXTNS,u'is-list-header'), None): cnv_boolean,
	((TEXTNS,u'isbn'), None): cnv_string,
	((TEXTNS,u'issn'), None): cnv_string,
	((TEXTNS,u'issn'), None): cnv_string,
	((TEXTNS,u'journal'), None): cnv_string,
	((TEXTNS,u'key'), None): cnv_string,
	((TEXTNS,u'key1'), None): cnv_string,
	((TEXTNS,u'key1-phonetic'), None): cnv_string,
	((TEXTNS,u'key2'), None): cnv_string,
	((TEXTNS,u'key2-phonetic'), None): cnv_string,
	((TEXTNS,u'kind'), None): cnv_string,
	((TEXTNS,u'label'), None): cnv_string,
	((TEXTNS,u'label-followed-by'), None): cnv_string,
	((TEXTNS,u'level'), None): cnv_positiveInteger,
	((TEXTNS,u'line-break'), None): cnv_boolean,
	((TEXTNS,u'line-number'), None): cnv_string,
	((TEXTNS,u'list-id'), None): cnv_NCName,
	((TEXTNS,u'list-level-position-and-space-mode'), None): cnv_string,
	((TEXTNS,u'list-tab-stop-position'), None): cnv_length,
	((TEXTNS,u'main-entry'), None): cnv_boolean,
	((TEXTNS,u'main-entry-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'master-page-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'min-label-distance'), None): cnv_string,
	((TEXTNS,u'min-label-width'), None): cnv_string,
	((TEXTNS,u'month'), None): cnv_string,
	((TEXTNS,u'name'), None): cnv_string,
	((TEXTNS,u'note-class'), None): cnv_textnoteclass,
	((TEXTNS,u'note'), None): cnv_string,
	((TEXTNS,u'number'), None): cnv_string,
	((TEXTNS,u'number-lines'), None): cnv_boolean,
	((TEXTNS,u'number-position'), None): cnv_string,
	((TEXTNS,u'numbered-entries'), None): cnv_boolean,
	((TEXTNS,u'offset'), None): cnv_string,
	((TEXTNS,u'organizations'), None): cnv_string,
	((TEXTNS,u'outline-level'), None): cnv_string,
	((TEXTNS,u'page-adjust'), None): cnv_integer,
	((TEXTNS,u'pages'), None): cnv_string,
	((TEXTNS,u'placeholder-type'), None): cnv_string,
	((TEXTNS,u'prefix'), None): cnv_string,
	((TEXTNS,u'protected'), None): cnv_boolean,
	((TEXTNS,u'protection-key'), None): cnv_string,
	((TEXTNS,u'protection-key-digest-algorithm'), None): cnv_anyURI,
	((TEXTNS,u'publisher'), None): cnv_string,
	((TEXTNS,u'ref-name'), None): cnv_string,
	((TEXTNS,u'reference-format'), None): cnv_string,
	((TEXTNS,u'relative-tab-stop-position'), None): cnv_boolean,
	((TEXTNS,u'report-type'), None): cnv_string,
	((TEXTNS,u'restart-numbering'), None): cnv_boolean,
	((TEXTNS,u'restart-on-page'), None): cnv_boolean,
	((TEXTNS,u'row-number'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'school'), None): cnv_string,
	((TEXTNS,u'section-name'), None): cnv_string,
	((TEXTNS,u'select-page'), None): cnv_string,
	((TEXTNS,u'separation-character'), None): cnv_string,
	((TEXTNS,u'series'), None): cnv_string,
	((TEXTNS,u'sort-algorithm'), None): cnv_string,
	((TEXTNS,u'sort-ascending'), None): cnv_boolean,
	((TEXTNS,u'sort-by-position'), None): cnv_boolean,
	((TEXTNS,u'space-before'), None): cnv_string,
	((TEXTNS,u'start-numbering-at'), None): cnv_string,
	((TEXTNS,u'start-value'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'start-value'), None): cnv_positiveInteger,
	((TEXTNS,u'string-value'), None): cnv_string,
	((TEXTNS,u'string-value-if-false'), None): cnv_string,
	((TEXTNS,u'string-value-if-true'), None): cnv_string,
	((TEXTNS,u'string-value-phonetic'), None): cnv_string,
	((TEXTNS,u'style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'style-override'), None): cnv_StyleNameRef,
	((TEXTNS,u'suffix'), None): cnv_string,
	((TEXTNS,u'tab-ref'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'table-name'), None): cnv_string,
	((TEXTNS,u'table-type'), None): cnv_string,
	((TEXTNS,u'time-adjust'), None): cnv_duration,
	((TEXTNS,u'time-value'), None): cnv_dateTime,
	((TEXTNS,u'time-value'), None): cnv_time,
	((TEXTNS,u'title'), None): cnv_string,
	((TEXTNS,u'track-changes'), None): cnv_boolean,
	((TEXTNS,u'url'), None): cnv_string,
	((TEXTNS,u'use-caption'), None): cnv_boolean,
	((TEXTNS,u'use-chart-objects'), None): cnv_boolean,
	((TEXTNS,u'use-draw-objects'), None): cnv_boolean,
	((TEXTNS,u'use-floating-frames'), None): cnv_boolean,
	((TEXTNS,u'use-graphics'), None): cnv_boolean,
	((TEXTNS,u'use-index-marks'), None): cnv_boolean,
	((TEXTNS,u'use-index-source-styles'), None): cnv_boolean,
	((TEXTNS,u'use-keys-as-entries'), None): cnv_boolean,
	((TEXTNS,u'use-math-objects'), None): cnv_boolean,
	((TEXTNS,u'use-objects'), None): cnv_boolean,
	((TEXTNS,u'use-other-objects'), None): cnv_boolean,
	((TEXTNS,u'use-outline-level'), None): cnv_boolean,
	((TEXTNS,u'use-soft-page-breaks'), None): cnv_boolean,
	((TEXTNS,u'use-spreadsheet-objects'), None): cnv_boolean,
	((TEXTNS,u'use-tables'), None): cnv_boolean,
	((TEXTNS,u'value'), None): cnv_nonNegativeInteger,
	((TEXTNS,u'visited-style-name'), None): cnv_StyleNameRef,
	((TEXTNS,u'volume'), None): cnv_string,
	((TEXTNS,u'year'), None): cnv_string,
	((XFORMSNS,u'bind'), None): cnv_string,
	((XHTMLNS,u'about'), None): cnv_anyURI,
	((XHTMLNS,u'content'), None): cnv_string,
	((XHTMLNS,u'datatype'), None): cnv_anyURI,
	((XHTMLNS,u'property'), None): cnv_anyURI,
	((XLINKNS,u'actuate'), None): cnv_string,
	((XLINKNS,u'href'), None): cnv_anyURI,
	((XLINKNS,u'show'), None): cnv_xlinkshow,
	((XLINKNS,u'title'), None): cnv_string,
	((XLINKNS,u'type'), None): cnv_xlinktype,
	((XMLNS,u'id'), None): cnv_NCName,
}

class AttrConverters:
    def convert(self, attribute, value, element):
        """ Based on the element, figures out how to check/convert the attribute value
            All values are converted to string
        """
        conversion = attrconverters.get((attribute, element.qname), None)
        if conversion is not None:
            return conversion(attribute, value, element)
        else:
            conversion = attrconverters.get((attribute, None), None)
            if conversion is not None:
                return conversion(attribute, value, element)
        if sys.version_info[0]==2:
            return unicode(value)
        else:
            return str(value)
