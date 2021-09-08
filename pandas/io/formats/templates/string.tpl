{% for r in head %}
{% for c in r %}
{% if c.is_visible != False %}{{c.display_value}}{% if not loop.last %} {% endif %}{% endif %}
{% endfor %}

{% endfor %}
{% for r in body %}
{% for c in r %}
{% if c.is_visible != False %}{{c.display_value}}{% if not loop.last %} {% endif %}{% endif %}
{% endfor %}

{% endfor %}
