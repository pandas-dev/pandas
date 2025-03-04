===
XML
===

.. _io.read_xml:

Reading XML
'''''''''''

.. versionadded:: 1.3.0

The top-level :func:`~pandas.io.xml.read_xml` function can accept an XML
string/file/URL and will parse nodes and attributes into a pandas ``DataFrame``.

.. note::

   Since there is no standard XML structure where design types can vary in
   many ways, ``read_xml`` works best with flatter, shallow versions. If
   an XML document is deeply nested, use the ``stylesheet`` feature to
   transform XML into a flatter version.

Let's look at a few examples.

Read an XML string:

.. ipython:: python

    from io import StringIO
   xml = """<?xml version="1.0" encoding="UTF-8"?>
   <bookstore>
     <book category="cooking">
       <title lang="en">Everyday Italian</title>
       <author>Giada De Laurentiis</author>
       <year>2005</year>
       <price>30.00</price>
     </book>
     <book category="children">
       <title lang="en">Harry Potter</title>
       <author>J K. Rowling</author>
       <year>2005</year>
       <price>29.99</price>
     </book>
     <book category="web">
       <title lang="en">Learning XML</title>
       <author>Erik T. Ray</author>
       <year>2003</year>
       <price>39.95</price>
     </book>
   </bookstore>"""

   df = pd.read_xml(StringIO(xml))
   df

Read a URL with no options:

.. ipython:: python

   df = pd.read_xml("https://www.w3schools.com/xml/books.xml")
   df

Read in the content of the "books.xml" file and pass it to ``read_xml``
as a string:

.. ipython:: python

   from io import StringIO

   file_path = "books.xml"
   with open(file_path, "w") as f:
       f.write(xml)

   with open(file_path, "r") as f:
       df = pd.read_xml(StringIO(f.read()))
   df

Read in the content of the "books.xml" as instance of ``StringIO`` or
``BytesIO`` and pass it to ``read_xml``:

.. ipython:: python

   from io import StringIO

   with open(file_path, "r") as f:
       sio = StringIO(f.read())

   df = pd.read_xml(sio)
   df

.. ipython:: python

   with open(file_path, "rb") as f:
       bio = BytesIO(f.read())

   df = pd.read_xml(bio)
   df

Even read XML from AWS S3 buckets such as NIH NCBI PMC Article Datasets providing
Biomedical and Life Science Journals:

.. code-block:: python

   >>> df = pd.read_xml(
   ...    "s3://pmc-oa-opendata/oa_comm/xml/all/PMC1236943.xml",
   ...    xpath=".//journal-meta",
   ...)
   >>> df
         journal-id  journal-title  issn  publisher
   0 Cardiovasc Ultrasound Cardiovascular Ultrasound 1476-7120 NaN

With `lxml`_ as default ``parser``, you access the full-featured XML library
that extends Python's ElementTree API. One powerful tool is ability to query
nodes selectively or conditionally with more expressive XPath:

.. _lxml: https://lxml.de

.. ipython:: python

   df = pd.read_xml(file_path, xpath="//book[year=2005]")
   df

Specify only elements or only attributes to parse:

.. ipython:: python

   df = pd.read_xml(file_path, elems_only=True)
   df

.. ipython:: python

   df = pd.read_xml(file_path, attrs_only=True)
   df

.. ipython:: python
   :suppress:

   import os
   os.remove("books.xml")

XML documents can have namespaces with prefixes and default namespaces without
prefixes both of which are denoted with a special attribute ``xmlns``. In order
to parse by node under a namespace context, ``xpath`` must reference a prefix.

For example, below XML contains a namespace with prefix, ``doc``, and URI at
``https://example.com``. In order to parse ``doc:row`` nodes,
``namespaces`` must be used.

.. ipython:: python

   from io import StringIO

   xml = """<?xml version='1.0' encoding='utf-8'?>
   <doc:data xmlns:doc="https://example.com">
     <doc:row>
       <doc:shape>square</doc:shape>
       <doc:degrees>360</doc:degrees>
       <doc:sides>4.0</doc:sides>
     </doc:row>
     <doc:row>
       <doc:shape>circle</doc:shape>
       <doc:degrees>360</doc:degrees>
       <doc:sides/>
     </doc:row>
     <doc:row>
       <doc:shape>triangle</doc:shape>
       <doc:degrees>180</doc:degrees>
       <doc:sides>3.0</doc:sides>
     </doc:row>
   </doc:data>"""

   df = pd.read_xml(StringIO(xml),
                    xpath="//doc:row",
                    namespaces={"doc": "https://example.com"})
   df

Similarly, an XML document can have a default namespace without prefix. Failing
to assign a temporary prefix will return no nodes and raise a ``ValueError``.
But assigning *any* temporary name to correct URI allows parsing by nodes.

.. ipython:: python

   from io import StringIO

   xml = """<?xml version='1.0' encoding='utf-8'?>
   <data xmlns="https://example.com">
    <row>
      <shape>square</shape>
      <degrees>360</degrees>
      <sides>4.0</sides>
    </row>
    <row>
      <shape>circle</shape>
      <degrees>360</degrees>
      <sides/>
    </row>
    <row>
      <shape>triangle</shape>
      <degrees>180</degrees>
      <sides>3.0</sides>
    </row>
   </data>"""

   df = pd.read_xml(StringIO(xml),
                    xpath="//pandas:row",
                    namespaces={"pandas": "https://example.com"})
   df

However, if XPath does not reference node names such as default, ``/*``, then
``namespaces`` is not required.

.. note::

   Since ``xpath`` identifies the parent of content to be parsed, only immediate
   descendants which include child nodes or current attributes are parsed.
   Therefore, ``read_xml`` will not parse the text of grandchildren or other
   descendants and will not parse attributes of any descendant. To retrieve
   lower level content, adjust xpath to lower level. For example,

   .. ipython:: python
        :okwarning:

      from io import StringIO

      xml = """
      <data>
        <row>
          <shape sides="4">square</shape>
          <degrees>360</degrees>
        </row>
        <row>
          <shape sides="0">circle</shape>
          <degrees>360</degrees>
        </row>
        <row>
          <shape sides="3">triangle</shape>
          <degrees>180</degrees>
        </row>
      </data>"""

      df = pd.read_xml(StringIO(xml), xpath="./row")
      df

   shows the attribute ``sides`` on ``shape`` element was not parsed as
   expected since this attribute resides on the child of ``row`` element
   and not ``row`` element itself. In other words, ``sides`` attribute is a
   grandchild level descendant of ``row`` element. However, the ``xpath``
   targets ``row`` element which covers only its children and attributes.

With `lxml`_ as parser, you can flatten nested XML documents with an XSLT
script which also can be string/file/URL types. As background, `XSLT`_ is
a special-purpose language written in a special XML file that can transform
original XML documents into other XML, HTML, even text (CSV, JSON, etc.)
using an XSLT processor.

.. _lxml: https://lxml.de
.. _XSLT: https://www.w3.org/TR/xslt/

For example, consider this somewhat nested structure of Chicago "L" Rides
where station and rides elements encapsulate data in their own sections.
With below XSLT, ``lxml`` can transform original nested document into a flatter
output (as shown below for demonstration) for easier parse into ``DataFrame``:

.. ipython:: python

   from io import StringIO

   xml = """<?xml version='1.0' encoding='utf-8'?>
    <response>
     <row>
       <station id="40850" name="Library"/>
       <month>2020-09-01T00:00:00</month>
       <rides>
         <avg_weekday_rides>864.2</avg_weekday_rides>
         <avg_saturday_rides>534</avg_saturday_rides>
         <avg_sunday_holiday_rides>417.2</avg_sunday_holiday_rides>
       </rides>
     </row>
     <row>
       <station id="41700" name="Washington/Wabash"/>
       <month>2020-09-01T00:00:00</month>
       <rides>
         <avg_weekday_rides>2707.4</avg_weekday_rides>
         <avg_saturday_rides>1909.8</avg_saturday_rides>
         <avg_sunday_holiday_rides>1438.6</avg_sunday_holiday_rides>
       </rides>
     </row>
     <row>
       <station id="40380" name="Clark/Lake"/>
       <month>2020-09-01T00:00:00</month>
       <rides>
         <avg_weekday_rides>2949.6</avg_weekday_rides>
         <avg_saturday_rides>1657</avg_saturday_rides>
         <avg_sunday_holiday_rides>1453.8</avg_sunday_holiday_rides>
       </rides>
     </row>
    </response>"""

   xsl = """<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
      <xsl:output method="xml" omit-xml-declaration="no" indent="yes"/>
      <xsl:strip-space elements="*"/>
      <xsl:template match="/response">
         <xsl:copy>
           <xsl:apply-templates select="row"/>
         </xsl:copy>
      </xsl:template>
      <xsl:template match="row">
         <xsl:copy>
           <station_id><xsl:value-of select="station/@id"/></station_id>
           <station_name><xsl:value-of select="station/@name"/></station_name>
           <xsl:copy-of select="month|rides/*"/>
         </xsl:copy>
      </xsl:template>
    </xsl:stylesheet>"""

   output = """<?xml version='1.0' encoding='utf-8'?>
    <response>
      <row>
         <station_id>40850</station_id>
         <station_name>Library</station_name>
         <month>2020-09-01T00:00:00</month>
         <avg_weekday_rides>864.2</avg_weekday_rides>
         <avg_saturday_rides>534</avg_saturday_rides>
         <avg_sunday_holiday_rides>417.2</avg_sunday_holiday_rides>
      </row>
      <row>
         <station_id>41700</station_id>
         <station_name>Washington/Wabash</station_name>
         <month>2020-09-01T00:00:00</month>
         <avg_weekday_rides>2707.4</avg_weekday_rides>
         <avg_saturday_rides>1909.8</avg_saturday_rides>
         <avg_sunday_holiday_rides>1438.6</avg_sunday_holiday_rides>
      </row>
      <row>
         <station_id>40380</station_id>
         <station_name>Clark/Lake</station_name>
         <month>2020-09-01T00:00:00</month>
         <avg_weekday_rides>2949.6</avg_weekday_rides>
         <avg_saturday_rides>1657</avg_saturday_rides>
         <avg_sunday_holiday_rides>1453.8</avg_sunday_holiday_rides>
      </row>
    </response>"""

   df = pd.read_xml(StringIO(xml), stylesheet=StringIO(xsl))
   df

For very large XML files that can range in hundreds of megabytes to gigabytes, :func:`pandas.read_xml`
supports parsing such sizeable files using `lxml's iterparse`_ and `etree's iterparse`_
which are memory-efficient methods to iterate through an XML tree and extract specific elements and attributes.
without holding entire tree in memory.

.. versionadded:: 1.5.0

.. _`lxml's iterparse`: https://lxml.de/3.2/parsing.html#iterparse-and-iterwalk
.. _`etree's iterparse`: https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.iterparse

To use this feature, you must pass a physical XML file path into ``read_xml`` and use the ``iterparse`` argument.
Files should not be compressed or point to online sources but stored on local disk. Also, ``iterparse`` should be
a dictionary where the key is the repeating nodes in document (which become the rows) and the value is a list of
any element or attribute that is a descendant (i.e., child, grandchild) of repeating node. Since XPath is not
used in this method, descendants do not need to share same relationship with one another. Below shows example
of reading in Wikipedia's very large (12 GB+) latest article data dump.

.. code-block:: ipython

    In [1]: df = pd.read_xml(
    ...         "/path/to/downloaded/enwikisource-latest-pages-articles.xml",
    ...         iterparse = {"page": ["title", "ns", "id"]}
    ...     )
    ...     df
    Out[2]:
                                                         title   ns        id
    0                                       Gettysburg Address    0     21450
    1                                                Main Page    0     42950
    2                            Declaration by United Nations    0      8435
    3             Constitution of the United States of America    0      8435
    4                     Declaration of Independence (Israel)    0     17858
    ...                                                    ...  ...       ...
    3578760               Page:Black cat 1897 07 v2 n10.pdf/17  104    219649
    3578761               Page:Black cat 1897 07 v2 n10.pdf/43  104    219649
    3578762               Page:Black cat 1897 07 v2 n10.pdf/44  104    219649
    3578763      The History of Tom Jones, a Foundling/Book IX    0  12084291
    3578764  Page:Shakespeare of Stratford (1926) Yale.djvu/91  104     21450

    [3578765 rows x 3 columns]

.. _io.xml:

Writing XML
'''''''''''

.. versionadded:: 1.3.0

``DataFrame`` objects have an instance method ``to_xml`` which renders the
contents of the ``DataFrame`` as an XML document.

.. note::

   This method does not support special properties of XML including DTD,
   CData, XSD schemas, processing instructions, comments, and others.
   Only namespaces at the root level is supported. However, ``stylesheet``
   allows design changes after initial output.

Let's look at a few examples.

Write an XML without options:

.. ipython:: python

   geom_df = pd.DataFrame(
       {
           "shape": ["square", "circle", "triangle"],
           "degrees": [360, 360, 180],
           "sides": [4, np.nan, 3],
       }
   )

   print(geom_df.to_xml())


Write an XML with new root and row name:

.. ipython:: python

   print(geom_df.to_xml(root_name="geometry", row_name="objects"))

Write an attribute-centric XML:

.. ipython:: python

   print(geom_df.to_xml(attr_cols=geom_df.columns.tolist()))

Write a mix of elements and attributes:

.. ipython:: python

   print(
       geom_df.to_xml(
           index=False,
           attr_cols=['shape'],
           elem_cols=['degrees', 'sides'])
   )

Any ``DataFrames`` with hierarchical columns will be flattened for XML element names
with levels delimited by underscores:

.. ipython:: python

   ext_geom_df = pd.DataFrame(
       {
           "type": ["polygon", "other", "polygon"],
           "shape": ["square", "circle", "triangle"],
           "degrees": [360, 360, 180],
           "sides": [4, np.nan, 3],
       }
   )

   pvt_df = ext_geom_df.pivot_table(index='shape',
                                    columns='type',
                                    values=['degrees', 'sides'],
                                    aggfunc='sum')
   pvt_df

   print(pvt_df.to_xml())

Write an XML with default namespace:

.. ipython:: python

   print(geom_df.to_xml(namespaces={"": "https://example.com"}))

Write an XML with namespace prefix:

.. ipython:: python

   print(
       geom_df.to_xml(namespaces={"doc": "https://example.com"},
                      prefix="doc")
   )

Write an XML without declaration or pretty print:

.. ipython:: python

   print(
       geom_df.to_xml(xml_declaration=False,
                      pretty_print=False)
   )

Write an XML and transform with stylesheet:

.. ipython:: python

   from io import StringIO

   xsl = """<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
      <xsl:output method="xml" omit-xml-declaration="no" indent="yes"/>
      <xsl:strip-space elements="*"/>
      <xsl:template match="/data">
        <geometry>
          <xsl:apply-templates select="row"/>
        </geometry>
      </xsl:template>
      <xsl:template match="row">
        <object index="{index}">
          <xsl:if test="shape!='circle'">
              <xsl:attribute name="type">polygon</xsl:attribute>
          </xsl:if>
          <xsl:copy-of select="shape"/>
          <property>
              <xsl:copy-of select="degrees|sides"/>
          </property>
        </object>
      </xsl:template>
    </xsl:stylesheet>"""

   print(geom_df.to_xml(stylesheet=StringIO(xsl)))


XML Final Notes
'''''''''''''''

* All XML documents adhere to `W3C specifications`_. Both ``etree`` and ``lxml``
  parsers will fail to parse any markup document that is not well-formed or
  follows XML syntax rules. Do be aware HTML is not an XML document unless it
  follows XHTML specs. However, other popular markup types including KML, XAML,
  RSS, MusicML, MathML are compliant `XML schemas`_.

* For above reason, if your application builds XML prior to pandas operations,
  use appropriate DOM libraries like ``etree`` and ``lxml`` to build the necessary
  document and not by string concatenation or regex adjustments. Always remember
  XML is a *special* text file with markup rules.

* With very large XML files (several hundred MBs to GBs), XPath and XSLT
  can become memory-intensive operations. Be sure to have enough available
  RAM for reading and writing to large XML files (roughly about 5 times the
  size of text).

* Because XSLT is a programming language, use it with caution since such scripts
  can pose a security risk in your environment and can run large or infinite
  recursive operations. Always test scripts on small fragments before full run.

* The `etree`_ parser supports all functionality of both ``read_xml`` and
  ``to_xml`` except for complex XPath and any XSLT. Though limited in features,
  ``etree`` is still a reliable and capable parser and tree builder. Its
  performance may trail ``lxml`` to a certain degree for larger files but
  relatively unnoticeable on small to medium size files.

.. _`W3C specifications`: https://www.w3.org/TR/xml/
.. _`XML schemas`: https://en.wikipedia.org/wiki/List_of_types_of_XML_schemas
.. _`etree`: https://docs.python.org/3/library/xml.etree.elementtree.html
