# Team

## Contributors

_pandas_ is made with love by more than [2,000 volunteer contributors](https://github.com/pandas-dev/pandas/graphs/contributors).

If you want to support pandas development, you can find information in the [donations page](../donate.html).

## Active maintainers

<div class="card-group maintainers">
    {% for person in maintainers.active_with_github_info %}
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
    {% endfor %}
</div>

## Diversity and Inclusion

> _pandas_ expressly welcomes and encourages contributions from anyone who faces under-representation, discrimination in the technology industry
> or anyone willing to increase the diversity of our team.
> We have identified visible gaps and obstacles in sustaining diversity and inclusion in the open-source communities and we are proactive in increasing
> the diversity of our team.
> We have a [code of conduct](../community/coc.html) to ensure a friendly and welcoming environment.
> Please send an email to [pandas-code-of-conduct-committee](mailto:pandas-coc@googlegroups.com), if you think we can do a
> better job at achieving this goal.

## Governance

Wes McKinney is the Benevolent Dictator for Life (BDFL).

The project governance is available in the [project governance page](governance.html).

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

## Inactive maintainers

<ul>
    {% for person in maintainers.inactive_with_github_info %}
        <li>
            <a href="{{ person.blog or person.html_url }}">
                {{ person.name or person.login }}
            </a>
        </li>
    {% endfor %}
</ul>
