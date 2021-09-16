{% for r in head %}{{ parse_row(r, delimiter, line_width) }}{% endfor %}
{% for val in col_max_char %}{{ "-" * (val if val < max_colwidth else max_colwidth + delimiter|length) }}{% endfor %}

{% for r in body %}{{ parse_row(r, delimiter, line_width) }}{%  endfor %}

{% for val in col_max_char %}{{ val }}{% endfor %}
