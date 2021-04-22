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
   Styler.loader

Style application
-----------------
.. autosummary::
   :toctree: api/

   Styler.apply
   Styler.applymap
   Styler.where
   Styler.format
   Styler.set_td_classes
   Styler.set_table_styles
   Styler.set_table_attributes
   Styler.set_tooltips
   Styler.set_caption
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
   Styler.background_gradient
   Styler.bar

Style export and import
-----------------------
.. autosummary::
   :toctree: api/

   Styler.render
   Styler.export
   Styler.use
   Styler.to_excel
