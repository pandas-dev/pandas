\begin{table}{% if parse(table_styles, 'position') is not none %}[{{parse(table_styles, 'position')}}]{% endif %}
{%- if parse(table_styles, 'float') is not none%}

\{{parse(table_styles, 'float')}}
{%- endif %}

\begin{tabular}{% raw %}{{% endraw %}{% for item in head[0] %}l{% endfor %}{% raw %}}{% endraw %}

{% for item in head[0] %}{%- if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
\hline
{% for row in body %}
{% for item in row %}{% if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
{% endfor %}
\end{tabular}
{%- if caption %}

\caption{% raw %}{{% endraw %}{{caption}}{% raw %}}{% endraw %}
{%- endif %}
{%- if parse(table_styles, 'label') is not none %}

\label{% raw %}{{% endraw %}{{parse(table_styles, 'label')}}{% raw %}}{% endraw %}
{%- endif %}

\end{table}
