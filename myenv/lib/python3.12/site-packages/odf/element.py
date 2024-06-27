#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2007-2010 Søren Roug, European Environment Agency
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

# Note: This script has copied a lot of text from xml.dom.minidom.
# Whatever license applies to that file also applies to this file.
#
import sys, os.path
sys.path.append(os.path.dirname(__file__))
import re
import xml.dom
from xml.dom.minicompat import *
from odf.namespaces import nsdict
import odf.grammar as grammar
from odf.attrconverters import AttrConverters

if sys.version_info[0] == 3:
    unicode=str # unicode function does not exist
    unichr=chr  # unichr does not exist

_xml11_illegal_ranges = (
    (0x0, 0x0,),
    (0xd800, 0xdfff,),
    (0xfffe, 0xffff,),
)

_xml10_illegal_ranges = _xml11_illegal_ranges + (
    (0x01, 0x08,),
    (0x0b, 0x0c,),
    (0x0e, 0x1f,),
)

_xml_discouraged_ranges = (
    (0x7f, 0x84,),
    (0x86, 0x9f,),
)

if sys.maxunicode >= 0x10000:
    # modern or "wide" python build
    _xml_discouraged_ranges = _xml_discouraged_ranges + (
        (0x1fffe, 0x1ffff,),
        (0x2fffe, 0x2ffff,),
        (0x3fffe, 0x3ffff,),
        (0x4fffe, 0x4ffff,),
        (0x5fffe, 0x5ffff,),
        (0x6fffe, 0x6ffff,),
        (0x7fffe, 0x7ffff,),
        (0x8fffe, 0x8ffff,),
        (0x9fffe, 0x9ffff,),
        (0xafffe, 0xaffff,),
        (0xbfffe, 0xbffff,),
        (0xcfffe, 0xcffff,),
        (0xdfffe, 0xdffff,),
        (0xefffe, 0xeffff,),
        (0xffffe, 0xfffff,),
        (0x10fffe, 0x10ffff,),
    )
# else "narrow" python build - only possible with old versions

def _range_seq_to_re(range_seq):
    # range pairs are specified as closed intervals
    return re.compile(u"[{}]".format(
        u"".join(
            u"{}-{}".format(re.escape(unichr(lo)), re.escape(unichr(hi)))
            for lo, hi in range_seq
        )
    ), flags=re.UNICODE)

_xml_filtered_chars_re = _range_seq_to_re(_xml10_illegal_ranges + _xml_discouraged_ranges)

def _handle_unrepresentable(data):
    return _xml_filtered_chars_re.sub(u"\ufffd", data)

# The following code is pasted form xml.sax.saxutils
# Tt makes it possible to run the code without the xml sax package installed
# To make it possible to have <rubbish> in your text elements, it is necessary to escape the texts
def _escape(data, entities={}):
    """ Escape &, <, and > in a string of data.

        You can escape other strings of data by passing a dictionary as
        the optional entities parameter.  The keys and values must all be
        strings; each key will be replaced with its corresponding value.
    """
    data = data.replace("&", "&amp;")
    data = data.replace("<", "&lt;")
    data = data.replace(">", "&gt;")
    for chars, entity in entities.items():
        data = data.replace(chars, entity)
    return data

def _sanitize(data, entities={}):
    return _escape(_handle_unrepresentable(data), entities=entities)

def _quoteattr(data, entities={}):
    """ Escape and quote an attribute value.

        Escape &, <, and > in a string of data, then quote it for use as
        an attribute value.  The \" character will be escaped as well, if
        necessary.

        You can escape other strings of data by passing a dictionary as
        the optional entities parameter.  The keys and values must all be
        strings; each key will be replaced with its corresponding value.
    """
    entities['\n']='&#10;'
    entities['\r']='&#12;'
    data = _sanitize(data, entities)
    if '"' in data:
        if "'" in data:
            data = '"%s"' % data.replace('"', "&quot;")
        else:
            data = "'%s'" % data
    else:
        data = '"%s"' % data
    return data

def _nssplit(qualifiedName):
    """ Split a qualified name into namespace part and local part.  """
    fields = qualifiedName.split(':', 1)
    if len(fields) == 2:
        return fields
    else:
        return (None, fields[0])

def _nsassign(namespace):
    return nsdict.setdefault(namespace,"ns" + str(len(nsdict)))


# Exceptions
class IllegalChild(Exception):
    """ Complains if you add an element to a parent where it is not allowed """
class IllegalText(Exception):
    """ Complains if you add text or cdata to an element where it is not allowed """

class Node(xml.dom.Node):
    """ super class for more specific nodes """
    parentNode = None
    nextSibling = None
    previousSibling = None

    def hasChildNodes(self):
        """ Tells whether this element has any children; text nodes,
            subelements, whatever.
        """
        if self.childNodes:
            return True
        else:
            return False

    def _get_childNodes(self):
        return self.childNodes

    def _get_firstChild(self):
        if self.childNodes:
            return self.childNodes[0]

    def _get_lastChild(self):
        if self.childNodes:
            return self.childNodes[-1]

    def insertBefore(self, newChild, refChild):
        """ Inserts the node newChild before the existing child node refChild.
            If refChild is null, insert newChild at the end of the list of children.
        """
        if newChild.nodeType not in self._child_node_types:
            raise IllegalChild( "%s cannot be child of %s" % (newChild.tagName, self.tagName))
        if newChild.parentNode is not None:
            newChild.parentNode.removeChild(newChild)
        if refChild is None:
            self.appendChild(newChild)
        else:
            try:
                index = self.childNodes.index(refChild)
            except ValueError:
                raise xml.dom.NotFoundErr()
            self.childNodes.insert(index, newChild)
            newChild.nextSibling = refChild
            refChild.previousSibling = newChild
            if index:
                node = self.childNodes[index-1]
                node.nextSibling = newChild
                newChild.previousSibling = node
            else:
                newChild.previousSibling = None
            newChild.parentNode = self
        return newChild

    def appendChild(self, newChild):
        """ Adds the node newChild to the end of the list of children of this node.
            If the newChild is already in the tree, it is first removed.
        """
        if newChild.nodeType == self.DOCUMENT_FRAGMENT_NODE:
            for c in tuple(newChild.childNodes):
                self.appendChild(c)
            ### The DOM does not clearly specify what to return in this case
            return newChild
        if newChild.nodeType not in self._child_node_types:
            raise IllegalChild( "<%s> is not allowed in %s" % ( newChild.tagName, self.tagName))
        if newChild.parentNode is not None:
            newChild.parentNode.removeChild(newChild)
        _append_child(self, newChild)
        newChild.nextSibling = None
        return newChild

    def removeChild(self, oldChild):
        """ Removes the child node indicated by oldChild from the list of children, and returns it.
        """
        #FIXME: update ownerDocument.element_dict or find other solution
        try:
            self.childNodes.remove(oldChild)
        except ValueError:
            raise xml.dom.NotFoundErr()
        if oldChild.nextSibling is not None:
            oldChild.nextSibling.previousSibling = oldChild.previousSibling
        if oldChild.previousSibling is not None:
            oldChild.previousSibling.nextSibling = oldChild.nextSibling
        oldChild.nextSibling = oldChild.previousSibling = None
        if self.ownerDocument:
            self.ownerDocument.remove_from_caches(oldChild)
        oldChild.parentNode = None
        return oldChild

    def __str__(self):
        val = []
        for c in self.childNodes:
            val.append(str(c))
        return ''.join(val)

    def __unicode__(self):
        val = []
        for c in self.childNodes:
            val.append(unicode(c))
        return u''.join(val)

defproperty(Node, "firstChild", doc="First child node, or None.")
defproperty(Node, "lastChild",  doc="Last child node, or None.")

def _append_child(self, node):
    # fast path with less checks; usable by DOM builders if careful
    childNodes = self.childNodes
    if childNodes:
        last = childNodes[-1]
        node.__dict__["previousSibling"] = last
        last.__dict__["nextSibling"] = node
    childNodes.append(node)
    node.__dict__["parentNode"] = self

class Childless:
    """ Mixin that makes childless-ness easy to implement and avoids
        the complexity of the Node methods that deal with children.
    """

    attributes = None
    childNodes = EmptyNodeList()
    firstChild = None
    lastChild = None

    def _get_firstChild(self):
        return None

    def _get_lastChild(self):
        return None

    def appendChild(self, node):
        """ Raises an error """
        raise xml.dom.HierarchyRequestErr(
            self.tagName + " nodes cannot have children")

    def hasChildNodes(self):
        return False

    def insertBefore(self, newChild, refChild):
        """ Raises an error """
        raise xml.dom.HierarchyRequestErr(
            self.tagName + " nodes do not have children")

    def removeChild(self, oldChild):
        """ Raises an error """
        raise xml.dom.NotFoundErr(
            self.tagName + " nodes do not have children")

    def replaceChild(self, newChild, oldChild):
        """ Raises an error """
        raise xml.dom.HierarchyRequestErr(
            self.tagName + " nodes do not have children")

class Text(Childless, Node):
    nodeType = Node.TEXT_NODE
    tagName = "Text"

    def __init__(self, data):
        self.data = data

    def __str__(self):
        return self.data

    def __unicode__(self):
        return self.data

    def toXml(self,level,f):
        """ Write XML in UTF-8 """
        if self.data:
            f.write(_sanitize(unicode(self.data)))
    
class CDATASection(Text, Childless):
    nodeType = Node.CDATA_SECTION_NODE

    def toXml(self,level,f):
        """ Generate XML output of the node. If the text contains "]]>", then
            escape it by going out of CDATA mode (]]>), then write the string
            and then go into CDATA mode again. (<![CDATA[)
        """
        if self.data:
            f.write('<![CDATA[%s]]>' % self.data.replace(']]>',']]>]]><![CDATA['))

class Element(Node):
    """ Creates a arbitrary element and is intended to be subclassed not used on its own.
        This element is the base of every element it defines a class which resembles
        a xml-element. The main advantage of this kind of implementation is that you don't
        have to create a toXML method for every different object. Every element
        consists of an attribute, optional subelements, optional text and optional cdata.
    """

    nodeType = Node.ELEMENT_NODE
    namespaces = {}  # Due to shallow copy this is a static variable

    _child_node_types = (Node.ELEMENT_NODE,
                         Node.PROCESSING_INSTRUCTION_NODE,
                         Node.COMMENT_NODE,
                         Node.TEXT_NODE,
                         Node.CDATA_SECTION_NODE,
                         Node.ENTITY_REFERENCE_NODE)

    def __init__(self, attributes=None, text=None, cdata=None, qname=None, qattributes=None, check_grammar=True, **args):
        if qname is not None:
            self.qname = qname
        assert(hasattr(self, 'qname'))
        self.ownerDocument = None
        self.childNodes=[]
        self.allowed_children = grammar.allowed_children.get(self.qname)
        prefix = self.get_nsprefix(self.qname[0])
        self.tagName = prefix + ":" + self.qname[1]
        if text is not None:
            self.addText(text)
        if cdata is not None:
            self.addCDATA(cdata)

        allowed_attrs = self.allowed_attributes()
        if allowed_attrs is not None:
            allowed_args = [ a[1].lower().replace('-','') for a in allowed_attrs]
        self.attributes={}
        # Load the attributes from the 'attributes' argument
        if attributes:
            for attr, value in attributes.items():
                self.setAttribute(attr, value)
        # Load the qualified attributes
        if qattributes:
            for attr, value in qattributes.items():
                self.setAttrNS(attr[0], attr[1], value)
        if allowed_attrs is not None:
            # Load the attributes from the 'args' argument
            for arg in args.keys():
                self.setAttribute(arg, args[arg])
        else:
            for arg in args.keys():  # If any attribute is allowed
                self.attributes[arg]=args[arg]
        if not check_grammar:
            return
        # Test that all mandatory attributes have been added.
        required = grammar.required_attributes.get(self.qname)
        if required:
            for r in required:
                if self.getAttrNS(r[0],r[1]) is None:
                    raise AttributeError( "Required attribute missing: %s in <%s>" % (r[1].lower().replace('-',''), self.tagName))

    def get_knownns(self, prefix):
        """ Odfpy maintains a list of known namespaces. In some cases a prefix is used, and
            we need to know which namespace it resolves to.
        """
        global nsdict
        for ns,p in nsdict.items():
            if p == prefix: return ns
        return None

    def get_nsprefix(self, namespace):
        """ Odfpy maintains a list of known namespaces. In some cases we have a namespace URL,
            and needs to look up or assign the prefix for it.
        """
        if namespace is None: namespace = ""
        prefix = _nsassign(namespace)
        if not namespace in self.namespaces:
            self.namespaces[namespace] = prefix
        return prefix

    def allowed_attributes(self):
        return grammar.allowed_attributes.get(self.qname)

    def _setOwnerDoc(self, element):
        element.ownerDocument = self.ownerDocument
        for child in element.childNodes:
            self._setOwnerDoc(child)

    def addElement(self, element, check_grammar=True):
        """ adds an element to an Element

            Element.addElement(Element)
        """
        if check_grammar and self.allowed_children is not None:
            if element.qname not in self.allowed_children:
                raise IllegalChild( "<%s> is not allowed in <%s>" % ( element.tagName, self.tagName))
        self.appendChild(element)
        self._setOwnerDoc(element)
        if self.ownerDocument:
            self.ownerDocument.rebuild_caches(element)

    def addText(self, text, check_grammar=True):
        """ Adds text to an element
            Setting check_grammar=False turns off grammar checking
        """
        if check_grammar and self.qname not in grammar.allows_text:
            raise IllegalText( "The <%s> element does not allow text" % self.tagName)
        else:
            if text != '':
                self.appendChild(Text(text))

    def addCDATA(self, cdata, check_grammar=True):
        """ Adds CDATA to an element
            Setting check_grammar=False turns off grammar checking
        """
        if check_grammar and self.qname not in grammar.allows_text:
            raise IllegalText( "The <%s> element does not allow text" % self.tagName)
        else:
            self.appendChild(CDATASection(cdata))

    def removeAttribute(self, attr, check_grammar=True):
        """ Removes an attribute by name. """
        allowed_attrs = self.allowed_attributes()
        if allowed_attrs is None:
            if type(attr) == type(()):
                prefix, localname = attr
                self.removeAttrNS(prefix, localname)
            else:
                raise AttributeError( "Unable to add simple attribute - use (namespace, localpart)")
        else:
            # Construct a list of allowed arguments
            allowed_args = [ a[1].lower().replace('-','') for a in allowed_attrs]
            if check_grammar and attr not in allowed_args:
                raise AttributeError( "Attribute %s is not allowed in <%s>" % ( attr, self.tagName))
            i = allowed_args.index(attr)
            self.removeAttrNS(allowed_attrs[i][0], allowed_attrs[i][1])

    def setAttribute(self, attr, value, check_grammar=True):
        """ Add an attribute to the element
            This is sort of a convenience method. All attributes in ODF have
            namespaces. The library knows what attributes are legal and then allows
            the user to provide the attribute as a keyword argument and the
            library will add the correct namespace.
            Must overwrite, If attribute already exists.
        """
        if attr == 'parent' and value is not None:
            value.addElement(self)
        else:
            allowed_attrs = self.allowed_attributes()
            if allowed_attrs is None:
                if type(attr) == type(()):
                    prefix, localname = attr
                    self.setAttrNS(prefix, localname, value)
                else:
                    raise AttributeError( "Unable to add simple attribute - use (namespace, localpart)")
            else:
                # Construct a list of allowed arguments
                allowed_args = [ a[1].lower().replace('-','') for a in allowed_attrs]
                if check_grammar and attr not in allowed_args:
                    raise AttributeError( "Attribute %s is not allowed in <%s>" % ( attr, self.tagName))
                i = allowed_args.index(attr)
                self.setAttrNS(allowed_attrs[i][0], allowed_attrs[i][1], value)

    def setAttrNS(self, namespace, localpart, value):
        """ Add an attribute to the element
            In case you need to add an attribute the library doesn't know about
            then you must provide the full qualified name
            It will not check that the attribute is legal according to the schema.
            Must overwrite, If attribute already exists.
        """
        allowed_attrs = self.allowed_attributes()
        prefix = self.get_nsprefix(namespace)
#       if allowed_attrs and (namespace, localpart) not in allowed_attrs:
#           raise AttributeError( "Attribute %s:%s is not allowed in element <%s>" % ( prefix, localpart, self.tagName))
        c = AttrConverters()
        self.attributes[(namespace, localpart)] = c.convert((namespace, localpart), value, self)

    def getAttrNS(self, namespace, localpart):
        """
        gets an attribute, given a namespace and a key
        @param namespace a unicode string or a bytes: the namespace
        @param localpart a unicode string or a bytes:
        the key to get the attribute
        @return an attribute as a unicode string or a bytes: if both paramters
        are byte strings, it will be a bytes; if both attributes are
        unicode strings, it will be a unicode string
        """
        prefix = self.get_nsprefix(namespace)
        result = self.attributes.get((namespace, localpart))

        assert(
            (type(namespace), type(namespace), type(namespace) == \
                 type(b""), type(b""), type(b"")) or
            (type(namespace), type(namespace), type(namespace) == \
                 type(u""), type(u""), type(u""))
            )

        return result

    def removeAttrNS(self, namespace, localpart):
        del self.attributes[(namespace, localpart)]

    def getAttribute(self, attr):
        """ Get an attribute value. The method knows which namespace the attribute is in
        """
        allowed_attrs = self.allowed_attributes()
        if allowed_attrs is None:
            if type(attr) == type(()):
                prefix, localname = attr
                return self.getAttrNS(prefix, localname)
            else:
                raise AttributeError( "Unable to get simple attribute - use (namespace, localpart)")
        else:
            # Construct a list of allowed arguments
            allowed_args = [ a[1].lower().replace('-','') for a in allowed_attrs]
            i = allowed_args.index(attr)
            return self.getAttrNS(allowed_attrs[i][0], allowed_attrs[i][1])

    def write_open_tag(self, level, f):
        f.write(('<'+self.tagName))
        if level == 0:
            for namespace, prefix in self.namespaces.items():
                f.write(u' xmlns:' + prefix + u'="'+ _sanitize(str(namespace))+'"')
        for qname in self.attributes.keys():
            prefix = self.get_nsprefix(qname[0])
            f.write(u' '+_sanitize(str(prefix+u':'+qname[1]))+u'='+_quoteattr(unicode(self.attributes[qname])))
        f.write(u'>')

    def write_close_tag(self, level, f):
        f.write('</'+self.tagName+'>')

    def toXml(self, level, f):
        """
        Generate an XML stream out of the tree structure
        @param level integer: level in the XML tree; zero at root of the tree
        @param f an open writable file able to accept unicode strings
        """
        f.write(u'<'+self.tagName)
        if level == 0:
            for namespace, prefix in self.namespaces.items():
                f.write(u' xmlns:' + prefix + u'="'+ _sanitize(str(namespace))+u'"')
        for qname in self.attributes.keys():
            prefix = self.get_nsprefix(qname[0])
            f.write(u' '+_sanitize(unicode(prefix+':'+qname[1]))+u'='+_quoteattr(unicode(self.attributes[qname])))
        if self.childNodes:
            f.write(u'>')
            for element in self.childNodes:
                element.toXml(level+1,f)
            f.write(u'</'+self.tagName+'>')
        else:
            f.write(u'/>')

    def _getElementsByObj(self, obj, accumulator):
        if self.qname == obj.qname:
            accumulator.append(self)
        for e in self.childNodes:
            if e.nodeType == Node.ELEMENT_NODE:
                accumulator = e._getElementsByObj(obj, accumulator)
        return accumulator

    def getElementsByType(self, element):
        """ Gets elements based on the type, which is function from text.py, draw.py etc. """
        obj = element(check_grammar=False)
        return self._getElementsByObj(obj,[])

    def isInstanceOf(self, element):
        """ This is a check to see if the object is an instance of a type """
        obj = element(check_grammar=False)
        return self.qname == obj.qname
