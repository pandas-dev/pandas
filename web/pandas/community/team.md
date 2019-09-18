# Team

## Contributors

_pandas_ is made with love by more than [1,500 volunteer contributors](https://github.com/pandas-dev/pandas/graphs/contributors).

If you want to support pandas development, you can find information in the [donations page](../donate.html).

## Maintainers

<div class="row maintainers">
    {% for row in maintainers.people | batch(6, "") %}
        <div class="card-group maintainers">
            {% for person in row %}
                {% if person %}
                    <div class="card">
                        <img class="card-img-top" alt="" src="{{ person.avatar_url }}"/>
                        <div class="card-body">
                            <h6 class="card-title">
                                {% if person.blog %}
                                    <a href="{{ person.blog }}">
                                        {{ person.name or person.login }}
                                    </a>
                                {% else %}
                                    {{ person.name or person.login }}
                                {% endif %}
                            </h6>
                            <p class="card-text small"><a href="{{ person.html_url }}">{{ person.login }}</a></p>
                        </div>
                    </div>
                {% else %}
                    <div class="card border-0"></div>
                {% endif %}
            {% endfor %}
        </div>
    {% endfor %}
</div>

## BDFL

Wes McKinney is the Benevolent Dictator for Life (BDFL).

## Governance

The project governance is available in the [project governance documents](https://github.com/pandas-dev/pandas-governance).

## NumFOCUS

![](https://numfocus.org/wp-content/uploads/2018/01/optNumFocus_LRG.png)

_pandas_ is a Sponsored Project of [NumFOCUS](https://numfocus.org/), a 501(c)(3) nonprofit charity in the United States.
NumFOCUS provides _pandas_ with fiscal, legal, and administrative support to help ensure the
health and sustainability of the project. Visit numfocus.org for more information.

Donations to _pandas_ are managed by NumFOCUS. For donors in the United States, your gift is tax-deductible
to the extent provided by law. As with any donation, you should consult with your tax adviser about your particular tax situation.

## Code of conduct committee

<ul>
    {% for person in maintainers.coc %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>

## NumFOCUS committee

<ul>
    {% for person in maintainers.numfocus %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>

## Institutional partners

<ul>
    {% for company in partners.active if company.employs %}
        <li><a href="{{ company.url }}">{{ company.name }}</a> ({{ company.employs }})</li>
    {% endfor %}
</ul>

In-kind sponsors

- [Indeed](https://opensource.indeedeng.io/): Logo and website design
- Can we find a donor for the hosting (website, benchmarks,...?)

## Emeritus maintainers

<ul>
    {% for person in maintainers.emeritus %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>

## Past institutional partners

<ul>
    {% for company in partners.past %}
        <li><a href="{{ company.url }}">{{ company.name }}</a></li>
    {% endfor %}
</ul>
