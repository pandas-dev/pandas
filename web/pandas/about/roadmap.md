# Roadmap

This page provides an overview of the major themes in pandas'
development. Each of these items requires a relatively large amount of
effort to implement. These may be achieved more quickly with dedicated
funding or interest from contributors.

An item being on the roadmap does not mean that it will *necessarily*
happen, even with unlimited funding. During the implementation period we
may discover issues preventing the adoption of the feature.

Additionally, an item *not* being on the roadmap does not exclude it
from inclusion in pandas. The roadmap is intended for larger,
fundamental changes to the project that are likely to take months or
years of developer time. Smaller-scoped items will continue to be
tracked on our [issue tracker](https://github.com/pandas-dev/pandas/issues).

The roadmap is defined as a set of major enhancement proposals named PDEPs.
For more information about PDEPs, and how to submit one, please refer to
[PEDP-1]({{ base_url }}pdeps/0001-purpose-and-guidelines.html).

## PDEPs

{% for pdep_type in ["Under discussion", "Accepted", "Implemented", "Rejected"] %}

<h3 id="pdeps-{{pdep_type}}">{{ pdep_type.replace("_", " ").capitalize() }}</h3>

<ul>
{% for pdep in pdeps[pdep_type] %}
    <li><a href="{% if not pdep.url.startswith("http") %}{{ base_url }}{% endif %}{{ pdep.url }}">{{ pdep.title }}</a></li>
{% else %}
    <li>There are currently no PDEPs with this status</li>
{% endfor %}
</ul>

{% endfor %}
