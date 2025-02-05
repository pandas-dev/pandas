#table(
  columns: {{ head[0] | length }},
{% for r in head %}
  {% for c in r %}[{% if c["is_visible"] %}{{ c["display_value"] }}{% endif %}],{% if not loop.last %} {% endif%}{% endfor %}

{% endfor %}

{% for r in body %}
  {% for c in r %}[{% if c["is_visible"] %}{{ c["display_value"] }}{% endif %}],{% if not loop.last %} {% endif%}{% endfor %}

{% endfor %}
)
