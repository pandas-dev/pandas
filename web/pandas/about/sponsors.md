# Sponsors

## NumFOCUS

![](https://numfocus.org/wp-content/uploads/2018/01/optNumFocus_LRG.png)

_pandas_ is a Sponsored Project of [NumFOCUS](https://numfocus.org/), a 501(c)(3) nonprofit charity in the United States.
NumFOCUS provides _pandas_ with fiscal, legal, and administrative support to help ensure the
health and sustainability of the project. Visit numfocus.org for more information.

Donations to _pandas_ are managed by NumFOCUS. For donors in the United States, your gift is tax-deductible
to the extent provided by law. As with any donation, you should consult with your tax adviser about your particular tax situation.

## Become a sponsor

As a free and open source project, _pandas_ relies on the support of the community of users for its development.
If you work for an organization that uses and benefits from _pandas_, please consider supporting pandas. There
are different ways, such as employing people to work on pandas, funding the project, or becoming a
[NumFOCUS sponsor](https://numfocus.org/sponsors) to support the broader ecosystem. Please contact us at
[admin@numfocus.org](mailto:admin@numfocus.org) to discuss.

## Institutional partners

Institutional partners are companies and universities that support the project by employing contributors.
Current institutional partners include:

<ul>
    {% for company in sponsors.active if company.kind == "partner" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## Sponsors

Sponsors are organizations that provide funding for pandas. Current sponsors include:

<ul>
    {% for company in sponsors.active if company.kind == "regular" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## In-kind sponsors

In-kind sponsors are organizations that support pandas development with goods or services.
Current in-kind sponsors include:

<ul>
    {% for company in sponsors.inkind %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## Past institutional partners

<ul>
    {% for company in sponsors.past if company.kind == "partner" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a></li>
    {% endfor %}
</ul>
