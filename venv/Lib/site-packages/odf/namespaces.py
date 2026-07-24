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
__version__ = "1.4.1"

TOOLSVERSION = u"ODFPY/" + __version__

ANIMNS         = u"urn:oasis:names:tc:opendocument:xmlns:animation:1.0"
CHARTNS        = u"urn:oasis:names:tc:opendocument:xmlns:chart:1.0"
CHARTOOONS     = u"http://openoffice.org/2010/chart"
CONFIGNS       = u"urn:oasis:names:tc:opendocument:xmlns:config:1.0"
CSS3TNS        = u"http://www.w3.org/TR/css3-text/"
#DBNS           = u"http://openoffice.org/2004/database"
DBNS           = u"urn:oasis:names:tc:opendocument:xmlns:database:1.0"
DCNS           = u"http://purl.org/dc/elements/1.1/"
DOMNS          = u"http://www.w3.org/2001/xml-events"
DR3DNS         = u"urn:oasis:names:tc:opendocument:xmlns:dr3d:1.0"
DRAWNS         = u"urn:oasis:names:tc:opendocument:xmlns:drawing:1.0"
FIELDNS        = u"urn:openoffice:names:experimental:ooo-ms-interop:xmlns:field:1.0"
FONS           = u"urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0"
FORMNS         = u"urn:oasis:names:tc:opendocument:xmlns:form:1.0"
FORMXNS        = u"urn:openoffice:names:experimental:ooxml-odf-interop:xmlns:form:1.0"
GRDDLNS        = u"http://www.w3.org/2003/g/data-view#"
KOFFICENS      = u"http://www.koffice.org/2005/"
LOEXTNS        = u"urn:org:documentfoundation:names:experimental:office:xmlns:loext:1.0"
MANIFESTNS     = u"urn:oasis:names:tc:opendocument:xmlns:manifest:1.0"
MATHNS         = u"http://www.w3.org/1998/Math/MathML"
METANS         = u"urn:oasis:names:tc:opendocument:xmlns:meta:1.0"
NUMBERNS       = u"urn:oasis:names:tc:opendocument:xmlns:datastyle:1.0"
OFFICENS       = u"urn:oasis:names:tc:opendocument:xmlns:office:1.0"
OFNS           = u"urn:oasis:names:tc:opendocument:xmlns:of:1.2"
OOOCNS         = u"http://openoffice.org/2004/calc"
OOONS          = u"http://openoffice.org/2004/office"
OOOWNS         = u"http://openoffice.org/2004/writer"
PRESENTATIONNS = u"urn:oasis:names:tc:opendocument:xmlns:presentation:1.0"
RDFANS         = u"http://docs.oasis-open.org/opendocument/meta/rdfa#"
RPTNS          = u"http://openoffice.org/2005/report"
SCRIPTNS       = u"urn:oasis:names:tc:opendocument:xmlns:script:1.0"
SMILNS         = u"urn:oasis:names:tc:opendocument:xmlns:smil-compatible:1.0"
STYLENS        = u"urn:oasis:names:tc:opendocument:xmlns:style:1.0"
SVGNS          = u"urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"
TABLENS        = u"urn:oasis:names:tc:opendocument:xmlns:table:1.0"
TABLEOOONS     = u"http://openoffice.org/2009/table"
TEXTNS         = u"urn:oasis:names:tc:opendocument:xmlns:text:1.0"
XFORMSNS       = u"http://www.w3.org/2002/xforms"
XHTMLNS        = u"http://www.w3.org/1999/xhtml"
XLINKNS        = u"http://www.w3.org/1999/xlink"
XMLNS          = u"http://www.w3.org/XML/1998/namespace"
XSDNS          = u"http://www.w3.org/2001/XMLSchema"
XSINS          = u"http://www.w3.org/2001/XMLSchema-instance"

nsdict = {
   ANIMNS: u'anim',
   CHARTNS: u'chart',
   CHARTOOONS: u'chartooo',
   CONFIGNS: u'config',
   CSS3TNS: u'css3t',
   DBNS: u'db',
   DCNS: u'dc',
   DOMNS: u'dom',
   DR3DNS: u'dr3d',
   DRAWNS: u'draw',
   FIELDNS: u'field',
   FONS: u'fo',
   FORMNS: u'form',
   FORMXNS: u'formx',
   GRDDLNS: u'grddl',
   KOFFICENS: u'koffice',
   LOEXTNS: u'loext',
   MANIFESTNS: u'manifest',
   MATHNS: u'math',
   METANS: u'meta',
   NUMBERNS: u'number',
   OFFICENS: u'office',
   OFNS: u'of',
   OOONS: u'ooo',
   OOOWNS: u'ooow',
   OOOCNS: u'oooc',
   PRESENTATIONNS: u'presentation',
   RDFANS: u'rdfa',
   RPTNS:  u'rpt',
   SCRIPTNS: u'script',
   SMILNS: u'smil',
   STYLENS: u'style',
   SVGNS: u'svg',
   TABLENS: u'table',
   TABLEOOONS: u'tableooo',
   TEXTNS: u'text',
   XFORMSNS: u'xforms',
   XLINKNS: u'xlink',
   XHTMLNS: u'xhtml',
   XMLNS: u'xml',
   XSDNS: u'xsd',
   XSINS: u'xsi',
}
