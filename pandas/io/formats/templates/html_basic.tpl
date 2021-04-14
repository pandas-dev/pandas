{# this template excludes all style, ids and classes #}
{% block before_table %}{% endblock before_table %}
{% block table %}
<table>
{% block caption %}
{% if caption %}
  <caption>{{caption}}</caption>
{% endif %}
{% endblock caption %}
{% block thead %}
  <thead>
{% block before_head_rows %}{% endblock %}
{% for r in head %}
{% block head_tr scoped %}
    <tr>
{% for c in r %}
{% if c.is_visible != False %}
      <{{c.type}} {{c.attributes}}>{{c.value}}</{{c.type}}>
{% endif %}
{% endfor %}
    </tr>
{% endblock head_tr %}
{% endfor %}
{% block after_head_rows %}{% endblock %}
  </thead>
{% endblock thead %}
{% block tbody %}
  <tbody>
{% block before_rows %}{% endblock before_rows %}
{% for r in body %}
{% block tr scoped %}
    <tr>
{% for c in r %}
{% if c.is_visible != False %}
      <{{c.type}} {{c.attributes}}>{{c.display_value}}</{{c.type}}>
{% endif %}
{% endfor %}
    </tr>
{% endblock tr %}
{% endfor %}
{% block after_rows %}{% endblock after_rows %}
  </tbody>
{% endblock tbody %}
</table>
{% endblock table %}
{% block after_table %}{% endblock after_table %}
