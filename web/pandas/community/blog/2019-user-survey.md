Title: 2019 pandas user survey
Date: 2019-08-22

<style type="text/css">
table td {
    background: none;
}

table tr.even td {
    background: none;
}

table {
	text-shadow: none;
}

</style>

# 2019 pandas user survey

Pandas recently conducted a user survey to help guide future development.
Thanks to everyone who participated! This post presents the high-level results.

This analysis and the raw data can be found [on GitHub](https://github.com/pandas-dev/pandas-user-surveys) and run on Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pandas-dev/pandas-user-surveys/master?filepath=2019.ipynb)


We had about 1250 repsonses over the 15 days we ran the survey in the summer of 2019.

## About the Respondents

There was a fair amount of representation across pandas experience and frequeny of use, though the majority of respondents are on the more experienced side.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_4_0.png)




![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_5_0.png)


We included a few questions that were also asked in the [Python Developers Survey](https://www.jetbrains.com/research/python-developers-survey-2018/) so we could compare Pandas' population to Python's.

90% of our respondents use Python as a primary language (compared with 84% from the PSF survey).





    Yes    90.67%
    No      9.33%
    Name: Is Python your main language?, dtype: object



Windows users are well represented (see [Steve Dower's talk](https://www.youtube.com/watch?v=uoI57uMdDD4) on this topic).





    Linux      61.57%
    Windows    60.21%
    MacOS      42.75%
    Name: What Operating Systems do you use?, dtype: object



For environment isolation, [conda](https://conda.io/en/latest/) was the most popular.




![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_13_0.png)


Most repondents are Python 3 only.





    3        92.39%
    2 & 3     6.80%
    2         0.81%
    Name: Python 2 or 3?, dtype: object



## Pandas APIs

It can be hard for open source projects to know what features are actually used. We asked a few questions to get an idea.

CSV and Excel are (for better or worse) the most popular formats.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_18_0.png)


In preperation for a possible refactor of pandas internals, we wanted to get a sense for
how common wide (100s of columns or more) DataFrames are.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_20_0.png)


Pandas is slowly growing new exentension types. Categoricals are the most popular,
and the nullable integer type is already almost as popular as datetime with timezone.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_22_0.png)


More and better examples seem to be a high-priority development item.
Pandas recently received a NumFOCUS grant to improve our documentation,
which we're using to write tutorial-style documentation, which should help
meet this need.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_24_0.png)


We also asked about specific, commonly-requested features.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_26_0.png)


Of these, the clear standout is "scaling" to large datasets. A couple observations:

1. Perhaps pandas' documentation should do a better job of promoting libraries that provide scalable dataframes (like [Dask](https://dask.org), [vaex](https://dask.org), and [modin](https://modin.readthedocs.io/en/latest/))
2. Memory efficiency (perhaps from a native string data type, fewer internal copies, etc.) is a valuable goal.

After that, the next-most critical improvement is integer missing values. Those were actually added in [Pandas 0.24](https://pandas.pydata.org/pandas-docs/stable/whatsnew/v0.24.0.html#optional-integer-na-support), but they're not the default, and there's still some incompatibilites with the rest of pandas API.

Pandas is a less conservative library than, say, NumPy. We're approaching 1.0, but on the way we've made many deprecations and some outright API breaking changes. Fortunately, most people are OK with the tradeoff.





    Yes    94.89%
    No      5.11%
    Name: Is Pandas stable enough for you?, dtype: object



There's a perception (which is shared by many of the pandas maintainers) that the pandas API is too large. To measure that, we asked whether users thought that pandas' API was too large, too small, or just right.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_31_0.png)


Finally, we asked for an overall satisfaction with the library, from 1 (not very unsatisfied) to 5 (very satisfied).



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_33_0.png)


Most people are very satisfied. The average response is 4.39. I look forward to tracking this number over time.

If you're analyzing the raw data, be sure to share the results with us [@pandas_dev](https://twitter.com/pandas_dev).
