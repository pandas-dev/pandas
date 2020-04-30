Title: pandas 1.0
Date: 2020-01-29

# pandas 1.0

Today pandas celebrates its 1.0.0 release. In many ways this is just a normal release with a host of new features, performance improvements, and bug fixes, which are documented in our [release notes](https://pandas.pydata.org/pandas-docs/version/1.0.0/whatsnew/v1.0.0.html). But it’s also something a bit more — a milestone for the project beyond just the commits. We wanted to take some time to reflect on where we've been and where we're going.

## Reflections

The world of scientific Python has changed a lot since pandas was started.  In 2011, [the ecosystem was fragmented](https://wesmckinney.com/blog/a-roadmap-for-rich-scientific-data-structures-in-python/): a standard *rich* data structure for statistics and data science had yet to emerge. This echos a similar story for NumPy, which consolidated array efforts that were [previously fragmented](https://numpy.org/old_array_packages.html).

Over the subsequent years, pandas emerged as a *de facto* standard. It’s used by data scientists and analysts and as a data structure for other libraries to build on top of. StackOverflow [cited pandas](https://stackoverflow.blog/2017/09/14/python-growing-quickly/) as one of the reasons for Python being the fastest growing major programming language.

![Growth of pandas](https://149351115.v2.pressablecdn.com/wp-content/uploads/2017/09/related_tags_over_time-1-1000x1000.png)

Today, the ecosystem is in another phase of exploration.
Several new DataFrame implementations are cropping up to fill needs not met by pandas.
We're [working with those projects](https://datapythonista.me/blog/dataframe-summit-at-euroscipy.html) to establish shared standards and semantics for rich data structures.

## Community and Project Health

This release cycle is the first to involve any kind of grant funding for pandas. [Pandas received funding](https://chanzuckerberg.com/eoss/proposals/) as part of the CZI’s [*Essential Open Source Software for Science*](https://medium.com/@cziscience/the-invisible-foundations-of-biomedicine-4ab7f8d4f5dd) [program](https://medium.com/@cziscience/the-invisible-foundations-of-biomedicine-4ab7f8d4f5dd). The pandas project relies overwhelmingly on volunteer contributors. These volunteer contributions are shepherded and augmented by some maintainers who are given time from their employers — our [institutional partners](https://github.com/pandas-dev/pandas-governance/blob/master/people.md#institutional-partners). The largest work item in our grant award was library maintenance, which specifically includes working with community members to address our large backlog of open issues and pull requests.

While a “1.0.0” version might seem arbitrary or anti-climactic (given that pandas as a codebase is nearly 12 years old), we see it as a symbolic milestone celebrating the growth of our core developer team and depth of our contributor base.  Few open source projects are ever truly “done” and pandas is no different. We recognize the essential role that pandas now occupies, and we intend to continue to evolve the project and adapt to the needs of the world’s data wranglers.

## Going Forward

Our [roadmap](https://pandas.pydata.org/pandas-docs/version/1.0.0/development/roadmap.html) contains an up-to-date listing of where we see the project heading over the next few years.
Needless to say, there's still plenty to do.

Check out the [release notes](https://pandas.pydata.org/pandas-docs/version/1.0.0/whatsnew/v1.0.0.html) and visit the [installation page](https://pandas.pydata.org/pandas-docs/version/1.0.0/getting_started/install.html) for instructions on updating to pandas 1.0.
