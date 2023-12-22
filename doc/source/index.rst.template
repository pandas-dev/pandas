:notoc:

.. pandas documentation master file, created by

.. module:: pandas

********************
pandas documentation
********************

**Date**: |today| **Version**: |version|

**Download documentation**: `Zipped HTML <pandas.zip>`__

**Previous versions**: Documentation of previous pandas versions is available at
`pandas.pydata.org <https://pandas.pydata.org/>`__.

**Useful links**:
`Binary Installers <https://pypi.org/project/pandas>`__ |
`Source Repository <https://github.com/pandas-dev/pandas>`__ |
`Issues & Ideas <https://github.com/pandas-dev/pandas/issues>`__ |
`Q&A Support <https://stackoverflow.com/questions/tagged/pandas>`__ |
`Mailing List <https://groups.google.com/g/pydata>`__

:mod:`pandas` is an open source, BSD-licensed library providing high-performance,
easy-to-use data structures and data analysis tools for the `Python <https://www.python.org/>`__
programming language.

{% if not single_doc or single_doc == "index.rst" -%}
.. grid:: 1 2 2 2
    :gutter: 4
    :padding: 2 2 0 0
    :class-container: sd-text-center

    .. grid-item-card:: Getting started
        :img-top: _static/index_getting_started.svg
        :class-card: intro-card
        :shadow: md

        New to *pandas*? Check out the getting started guides. They contain an
        introduction to *pandas'* main concepts and links to additional tutorials.

        +++

        .. button-ref:: getting_started
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the getting started guides

    .. grid-item-card::  User guide
        :img-top: _static/index_user_guide.svg
        :class-card: intro-card
        :shadow: md

        The user guide provides in-depth information on the
        key concepts of pandas with useful background information and explanation.

        +++

        .. button-ref:: user_guide
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the user guide

    .. grid-item-card::  API reference
        :img-top: _static/index_api.svg
        :class-card: intro-card
        :shadow: md

        The reference guide contains a detailed description of
        the pandas API. The reference describes how the methods work and which parameters can
        be used. It assumes that you have an understanding of the key concepts.

        +++

        .. button-ref:: api
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the reference guide

    .. grid-item-card::  Developer guide
        :img-top: _static/index_contribute.svg
        :class-card: intro-card
        :shadow: md

        Saw a typo in the documentation? Want to improve
        existing functionalities? The contributing guidelines will guide
        you through the process of improving pandas.

        +++

        .. button-ref:: development
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the development guide

{% endif %}
{% if single_doc and single_doc.endswith('.rst') -%}
.. toctree::
    :maxdepth: 3
    :titlesonly:

    {{ single_doc[:-4] }}
{% elif single_doc and single_doc.count('.') <= 1 %}
.. autosummary::
    :toctree: reference/api/

    {{ single_doc }}
{% elif single_doc %}
.. autosummary::
    :toctree: reference/api/
    :template: autosummary/accessor_method.rst

    {{ single_doc }}
{% else -%}
.. toctree::
    :maxdepth: 3
    :hidden:
    :titlesonly:
{% endif %}
{% if not single_doc %}
    getting_started/index
    user_guide/index
    {% endif -%}
    {% if include_api -%}
    reference/index
    {% endif -%}
    {% if not single_doc -%}
    development/index
    whatsnew/index
{% endif %}
