\begin{table}
\centering
\begin{tabular}{% raw %}{{% endraw %}{% for item in head[0] %}l{% endfor %}{% raw %}}{% endraw %}

{% for item in head[0] %}{% if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
\hline
{% for row in body %}
{% for item in row %}{% if not loop.first %} & {% endif %}{{item.display_value}}{% endfor %} \\
{% endfor %}
\end{tabular}
\end{table}
