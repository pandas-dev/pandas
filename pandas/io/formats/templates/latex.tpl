\begin{table}
{%- set position = parse(table_styles, 'position') %}
{%- if position is not none %}
[{{position}}]
{% endif %}
{% set float = parse(table_styles, 'float') %}
{% if float is not none%}
\{{float}}
{% endif %}
{% if caption %}
\caption{% raw %}{{% endraw %}{{caption}}{% raw %}}{% endraw %}

{% endif %}
{% set label = parse(table_styles, 'label') %}
{% if label is not none %}
\label{% raw %}{{% endraw %}{{label}}{% raw %}}{% endraw %}

{% endif %}
\begin{tabular}{% raw %}{{% endraw %}{% for c in head[0] %}{% if c.is_visible != False %}l{% endif %}{% endfor %}{% raw %}}{% endraw %}

{% set toprule = parse(table_styles, 'toprule') %}
{% if toprule is not none %}
\{{toprule}}
{% endif %}
{% for item in head[0] %}{%- if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
{% set midrule = parse(table_styles, 'midrule') %}
{% if midrule is not none %}
\{{midrule}}
{% endif %}
{% for row in body %}
{% for item in row %}{% if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
{% endfor %}
{% set bottomrule = parse(table_styles, 'bottomrule') %}
{% if bottomrule is not none %}
\{{bottomrule}}
{% endif %}
\end{tabular}
\end{table}
