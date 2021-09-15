{% for r in head %}{{ parse_row(r, delimiter, line_width) }}{% endfor %}
{% for r in body %}{{ parse_row(r, delimiter, line_width) }}{%  endfor %}
