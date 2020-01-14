# Sponsors

## NumFOCUS

![](https://numfocus.org/wp-content/uploads/2018/01/optNumFocus_LRG.png)

_pandas_ is a Sponsored Project of [NumFOCUS](https://numfocus.org/), a 501(c)(3) nonprofit charity in the United States.
NumFOCUS provides _pandas_ with fiscal, legal, and administrative support to help ensure the
health and sustainability of the project. Visit numfocus.org for more information.

Donations to _pandas_ are managed by NumFOCUS. For donors in the United States, your gift is tax-deductible
to the extent provided by law. As with any donation, you should consult with your tax adviser about your particular tax situation.

## Tidelift

_pandas_ is part of the [Tidelift subscription](https://tidelift.com/subscription/pkg/pypi-pandas?utm_source=pypi-pandas&utm_medium=referral&utm_campaign=readme).
You can support pandas by becoming a Tidelift subscriber.

## Institutional partners

Institutional Partners are companies and universities that support the project by employing contributors.
Current Institutional Partners include:

<ul>
    {% for company in partners.active if company.employs %}
        <li><a href="{{ company.url }}">{{ company.name }}</a> ({{ company.employs }})</li>
    {% endfor %}
</ul>

## In-kind sponsors

- [OVH](https://us.ovhcloud.com/): Hosting
- [Indeed](https://opensource.indeedeng.io/): Logo and website design

## Past institutional partners

<ul>
    {% for company in partners.past %}
        <li><a href="{{ company.url }}">{{ company.name }}</a></li>
    {% endfor %}
</ul>
