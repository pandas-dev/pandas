{{ header }}

.. _api.style:

=====
Style
=====
.. currentmodule:: pandas.io.formats.style

``Styler`` objects are returned by :attr:`pandas.DataFrame.style`.

Styler constructor
------------------
.. autosummary::
   :toctree: api/

   Styler
   Styler.from_custom_template

Styler properties
-----------------
.. autosummary::
   :toctree: api/

   Styler.env
   Styler.template_html
   Styler.template_html_style
   Styler.template_html_table
   Styler.template_latex
   Styler.template_string
   Styler.loader

Style application
-----------------
.. autosummary::
   :toctree: api/

   Styler.apply
   Styler.map
   Styler.apply_index
   Styler.map_index
   Styler.format
   Styler.format_index
   Styler.format_index_names
   Styler.relabel_index
   Styler.hide
   Styler.concat
   Styler.set_td_classes
   Styler.set_table_styles
   Styler.set_table_attributes
   Styler.set_tooltips
   Styler.set_caption
   Styler.set_sticky
   Styler.set_properties
   Styler.set_uuid
   Styler.clear
   Styler.pipe

Builtin styles
--------------
.. autosummary::
   :toctree: api/

   Styler.highlight_null
   Styler.highlight_max
   Styler.highlight_min
   Styler.highlight_between
   Styler.highlight_quantile
   Styler.background_gradient
   Styler.text_gradient
   Styler.bar

Style export and import
-----------------------
.. autosummary::
   :toctree: api/

   Styler.to_html
   Styler.to_latex
   Styler.to_excel
   Styler.to_string
   Styler.export
   Styler.use
