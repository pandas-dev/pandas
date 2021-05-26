{# Update the template_structure.html documentation too #}
{% if doctype_html %}
<!DOCTYPE html>
<html>
<head>
<meta charset="{{encoding}}">
{% if not exclude_styles %}{% include "html_style.tpl" %}{% endif %}
</head>
<body>
{% include "html_table.tpl" %}
</body>
</html>
{% elif not doctype_html %}
{% if not exclude_styles %}{% include "html_style.tpl" %}{% endif %}
{% include "html_table.tpl" %}
{% endif %}
