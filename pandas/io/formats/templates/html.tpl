{# Update the template_structure.html documentation too #}
{% if exclude_styles %}
{% extends "html_basic.tpl" %}
{% else %}
{% extends "html_styles.tpl" %}
{% endif %}
