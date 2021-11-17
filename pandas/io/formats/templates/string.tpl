{% for r in head %}{% for c in r %}
{{ c["display_value"] }}{% if not loop.last %}{{ delimiter }}{% endif %}
{% endfor %}{% endfor %}
{% if hrules %}
{% for val in col_max_char %}{{ "-" * (val if val < max_colwidth else max_colwidth + delimiter|length) }}{% endfor %}
{% endif %}
{% for r in body %}{{ parse_row(r, delimiter, line_width) }}{%  endfor %}
