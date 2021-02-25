{# Update the template_structure.html documentation too #}
{% if no_styles %}
{% extends "html_excl.tpl" %}
{% else %}
{% extends "html_all.tpl" %}
{% endif %}
