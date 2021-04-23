{# Update the template_structure.html documentation too #}
{% if doctype_html %}
<!DOCTYPE html>
<html>
<head>
<meta charset="{{encoding}}">
{% if not exclude_styles %}{% include "html_style.tpl" %}{% endif %}
</head>
<body>
{% if not exclude_styles %}{% include "html_body_inc.tpl" %}{% else %}{% include "html_body_exc.tpl" %}{% endif %}
</body>
</html>
{% elif not doctype_html %}
{% if exclude_styles %}
{% include "html_body_exc.tpl" %}
{% else %}
{% include "html_style.tpl" %}
{% include "html_body_inc.tpl" %}
{% endif %}
{% endif %}
