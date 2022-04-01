Title: Write up of the NumFOCUS grant to improve pandas benchmarks and diversity
Date: 2022-04-01

# Write up of the NumFOCUS grant to improve pandas benchmarks and diversity

*By Lucy Jiménez and Dorothy Kabarozi B.*

We want to share our experience working on **Improvements to the**
**ASV benchmarking framework and diversity efforts** sponsored by
[NumFOCUS](https://numfocus.org/) to the [pandas](https://pandas.pydata.org/)
project.

This grant focused on two aspects: the first one is to improve the
[asv library](https://asv.readthedocs.io/en/stable/), a tool used by
benchmarking Python packages and used by pandas; this project was
unmaintained, and the codebase was quite old; additionally, it didn't
adhere to modern standards, had Python 2 compatibility code that could
be removed, and also the CI could be improved. The second aspect is
encouraging more underrepresented groups to contribute to open source
projects. This grant was held over 10 weeks, working around 20 hours a
week. It was developed by Dorothy Kabarozi B. from Uganda and Lucy
Jiménez from Colombia, under the mentoring of Marc Garcia.

## Why were we part of the grant?

Even when we come from different backgrounds, Dorothy from systems
engineering and Lucy from computational chemistry, we have always been
interested in participating and contributing to open source software
projects. For that reason, we have been running the PyLadies meetups in
our communities ([PyLadies Kampala](https://twitter.com/pyladieskla),
[PyLadies Colombia](https://twitter.com/pyladies_co)) and have always
been on the lookout for any opportunities that lead us to contribute.

It all happened through Marc Garcia; he had put out a call ​through a post
on social media to mentor ladies from diverse backgrounds. Dorothy got to
be part of the pandas mentorship group. At the same time, Lucy was
co-organizer of the SciPy Latam conference, and it is from here she met
Marc, who was the speaker at that conference, and through this mutual
connection, we were able to learn about this benchmarks grant.

In brief, by attending conferences, meetups, and social media, you can
make connections and links that will lead you to these opportunities.

## Learning from the source code

At the beginning of the grant, we started from the basics. We noticed that
we could improve our skills in managing Git and GitHub. For example, we had
some troubles with the git workflow, so we had to read and practice more
about it. One of the valuable resources was the explanation from Marc about
[how to make an open source contribution](https://tubedu.org/w/kjnHEg72j76StmSFmjzbnE),
which we invite you to take a look at it.

We learned a lot from the source code and gained immense knowledge about
best practices and code quality through this grant. We have been working
on: updating the code to improve the style to follow the PEP-8 guidelines,
removing Python 2 compatibility code and six dependencies, and finding
unused code and removing it. We also learned about GitHub actions, and we
started building the CI on GitHub actions for the asv package; for that we
have been working on add linting with Flake8, testing with pytest, building
docs, and running CI on different python versions.

Additionally, we were able to identify bugs in the source code, review
pull request from other contributors, and create new issues, something we
thought only maintainers could do but not contributors. Finally, not only
is reviewing the code itself a learning experience, but also the structure
and folder hierarchy in the project started to be more transparent.

## Our experience

For this grant, we had a fantastic Mentor, Marc Garcia. He was always
willing to share his knowledge, explain unclear concepts and share helpful
feedback. Whenever we would implement that feedback, it felt easier to work
on more issues faster. We felt the growth from the time we started on this
project, and we will carry it along as we contribute to more open source
projects; this all goes back to Marc for his fantastic mentorship. It is
also important to note that we received feedback from other contributors,
stakeholders, and core devs during this process, which gave us a broader
look at the work in open source projects.

We also built a strong teamwork partnership. We helped each other a lot as
we had numerous one-on-one calls to understand the tasks better. We always
looked for ways to support each other from the technical side and encouraged
each other when needed. For us, it was professional and human growth.

## Running an open source software sprint

The knowledge and experience acquired in this process allowed us to
organize two virtual sprints. The events were carried out in the company
of local PyLadies communities; the first one was on February 26th with
[PyLadies Kampala](https://twitter.com/pyladieskla) and on March 21
with [PyLadies Colombia](https://bit.ly/sprint-asv).

While organizing these events, we learned how to organize and conduct a
virtual sprint. Some participants in the sprint ultimately had no idea
about open source, and it was great explaining open source concepts and
taking them through the Git workflow. Finally, they were able to make their
first contribution. We learned how to follow up on contributors, helping
them along the way until their PRs were merged and by reviewing their
contributions on GitHub.

The most outstanding achievement was mentoring new contributors and
sharing the knowledge acquired from this grant with others participants
in our respective communities. Most new contributors after the experience
have gone ahead to apply for outreach and the upcoming
[Google Summer of Code](https://summerofcode.withgoogle.com/)
to apply the skills they learned from these sprints.

## Conclusion

In conclusion, we learned a lot from this experience from the code part,
the workflow on the open source projects, how to be resilient in difficult
moments, and encouraging more women and people from our local communities
to contribute to open source projects.

Finally, if you want to be part of an open source project, we invite you
to check out GitHub repos for different projects you are interested in and
search for the easy issues to work on and get started. Also, you can contact
the maintainers of the projects with specific questions, search for the
open source communities in your country or contact us for more help.

## Acknowledgments

Many thanks to [NumFOCUS](https://numfocus.org/) for giving us this support
through [Small Development Grants](https://numfocus.org/programs/small-development-grants)
and Marc for the excellent mentoring he generously gave us throughout these
weeks.

We are looking forward to contributing more and impacting our communities
and the open source community!

___
If you want to know more, please don't hesitate to connect with us through
these channels:

*Lucy Jiménez*
* [Twitter](https://twitter.com/JimenezLucyJ)
* [LinkedIn](https://www.linkedin.com/in/lucy-j/)

*Dorothy Kabarozi*
* [Twitter](https://twitter.com/kizdorothy)
* [LinkedIn](https://www.linkedin.com/in/dorothy-kabarozi/)
