.. _tutorials:

*********
Tutorials
*********

This is a guide to many pandas tutorials, geared mainly for new users.

Internal Guides
===============

pandas' own :ref:`10 Minutes to pandas<10min>`.

More complex recipes are in the :ref:`Cookbook<cookbook>`.

A handy pandas `cheat sheet <http://pandas.pydata.org/Pandas_Cheat_Sheet.pdf>`_.

Community Guides
================

pandas Cookbook
---------------

The goal of this 2015 cookbook (by `Julia Evans <http://jvns.ca>`_) is to
give you some concrete examples for getting started with pandas. These
are examples with real-world data, and all the bugs and weirdness that
entails.

Here are links to the v0.2 release. For an up-to-date table of contents, see the `pandas-cookbook GitHub
repository <http://github.com/jvns/pandas-cookbook>`_. To run the examples in this tutorial, you'll need to
clone the GitHub repository and get IPython Notebook running.
See `How to use this cookbook <https://github.com/jvns/pandas-cookbook#how-to-use-this-cookbook>`_.

*  `A quick tour of the IPython Notebook: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/A%20quick%20tour%20of%20IPython%20Notebook.ipynb>`_
   Shows off IPython's awesome tab completion and magic functions.
*  `Chapter 1: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%201%20-%20Reading%20from%20a%20CSV.ipynb>`_
   Reading your data into pandas is pretty much the easiest thing. Even
   when the encoding is wrong!
*  `Chapter 2: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%202%20-%20Selecting%20data%20%26%20finding%20the%20most%20common%20complaint%20type.ipynb>`_
   It's not totally obvious how to select data from a pandas dataframe.
   Here we explain the basics (how to take slices and get columns)
*  `Chapter 3: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%203%20-%20Which%20borough%20has%20the%20most%20noise%20complaints%20%28or%2C%20more%20selecting%20data%29.ipynb>`_
   Here we get into serious slicing and dicing and learn how to filter
   dataframes in complicated ways, really fast.
*  `Chapter 4: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%204%20-%20Find%20out%20on%20which%20weekday%20people%20bike%20the%20most%20with%20groupby%20and%20aggregate.ipynb>`_
   Groupby/aggregate is seriously my favorite thing about pandas
   and I use it all the time. You should probably read this.
*  `Chapter 5:  <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%205%20-%20Combining%20dataframes%20and%20scraping%20Canadian%20weather%20data.ipynb>`_
   Here you get to find out if it's cold in Montreal in the winter
   (spoiler: yes). Web scraping with pandas is fun! Here we combine dataframes.
*  `Chapter 6:  <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%206%20-%20String%20Operations-%20Which%20month%20was%20the%20snowiest.ipynb>`_
   Strings with pandas are great. It has all these vectorized string
   operations and they're the best. We will turn a bunch of strings
   containing "Snow" into vectors of numbers in a trice.
*  `Chapter 7: <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%207%20-%20Cleaning%20up%20messy%20data.ipynb>`_
   Cleaning up messy data is never a joy, but with pandas it's easier.
*  `Chapter 8:  <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%208%20-%20How%20to%20deal%20with%20timestamps.ipynb>`_
   Parsing Unix timestamps is confusing at first but it turns out
   to be really easy.
*  `Chapter 9:  <http://nbviewer.jupyter.org/github/jvns/pandas-cookbook/blob/v0.2/cookbook/Chapter%209%20-%20Loading%20data%20from%20SQL%20databases.ipynb>`_
   Reading data from SQL databases.

.. _tutorial-qa-video-series:

pandas Q&A video series
-----------------------

This is a series of 31 video tutorials suitable for pandas beginners. Each
video teaches a related set of pandas functions using a real-world dataset.
All of the code from the series is available as a
`Jupyter Notebook <http://nbviewer.jupyter.org/github/justmarkham/pandas-videos/blob/master/pandas.ipynb>`_,
and all of the datasets are available on
`GitHub <https://github.com/justmarkham/pandas-videos>`_.

* 1 - `What is pandas? <https://www.youtube.com/watch?v=yzIMircGU5I&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=1>`_
* 2 - `How do I read a tabular data file into pandas? <https://www.youtube.com/watch?v=5_QXMwezPJE&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=2>`_
* 3 - `How do I select a pandas Series from a DataFrame? <https://www.youtube.com/watch?v=zxqjeyKP2Tk&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=3>`_
* 4 - `Why do some pandas commands end with parentheses (and others don't)? <https://www.youtube.com/watch?v=hSrDViyKWVk&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=4>`_
* 5 - `How do I rename columns in a pandas DataFrame? <https://www.youtube.com/watch?v=0uBirYFhizE&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=5>`_
* 6 - `How do I remove columns from a pandas DataFrame? <https://www.youtube.com/watch?v=gnUKkS964WQ&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=6>`_
* 7 - `How do I sort a pandas DataFrame or a Series? <https://www.youtube.com/watch?v=zY4doF6xSxY&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=7>`_
* 8 - `How do I filter rows of a pandas DataFrame by column value? <https://www.youtube.com/watch?v=2AFGPdNn4FM&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=8>`_
* 9 - `How do I apply multiple filter criteria to a pandas DataFrame? <https://www.youtube.com/watch?v=YPItfQ87qjM&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=9>`_
* 10 - `Your pandas questions answered! <https://www.youtube.com/watch?v=B-r9VuK80dk&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=10>`_
* 11 - `How do I use the "axis" parameter in pandas? <https://www.youtube.com/watch?v=PtO3t6ynH-8&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=11>`_
* 12 - `How do I use string methods in pandas? <https://www.youtube.com/watch?v=bofaC0IckHo&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=12>`_
* 13 - `How do I change the data type of a pandas Series? <https://www.youtube.com/watch?v=V0AWyzVMf54&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=13>`_
* 14 - `When should I use a "groupby" in pandas? <https://www.youtube.com/watch?v=qy0fDqoMJx8&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=14>`_
* 15 - `How do I explore a pandas Series? <https://www.youtube.com/watch?v=QTVTq8SPzxM&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=15>`_
* 16 - `How do I handle missing values in pandas? <https://www.youtube.com/watch?v=fCMrO_VzeL8&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=16>`_
* 17 - `What do I need to know about the pandas index? (Part 1) <https://www.youtube.com/watch?v=OYZNk7Z9s6I&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=17>`_
* 18 - `What do I need to know about the pandas index? (Part 2) <https://www.youtube.com/watch?v=15q-is8P_H4&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=18>`_
* 19 - `How do I select multiple rows and columns from a pandas DataFrame? <https://www.youtube.com/watch?v=xvpNA7bC8cs&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=19>`_
* 20 - `When should I use the "inplace" parameter in pandas? <https://www.youtube.com/watch?v=XaCSdr7pPmY&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=20>`_
* 21 - `How do I make my pandas DataFrame smaller and faster? <https://www.youtube.com/watch?v=wDYDYGyN_cw&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=21>`_
* 22 - `How do I use pandas with scikit-learn to create Kaggle submissions? <https://www.youtube.com/watch?v=ylRlGCtAtiE&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=22>`_
* 23 - `More of your pandas questions answered! <https://www.youtube.com/watch?v=oH3wYKvwpJ8&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=23>`_
* 24 - `How do I create dummy variables in pandas? <https://www.youtube.com/watch?v=0s_1IsROgDc&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=24>`_
* 25 - `How do I work with dates and times in pandas? <https://www.youtube.com/watch?v=yCgJGsg0Xa4&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=25>`_
* 26 - `How do I find and remove duplicate rows in pandas? <https://www.youtube.com/watch?v=ht5buXUMqkQ&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=26>`_
* 27 - `How do I avoid a SettingWithCopyWarning in pandas? <https://www.youtube.com/watch?v=4R4WsDJ-KVc&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=27>`_
* 28 - `How do I change display options in pandas? <https://www.youtube.com/watch?v=yiO43TQ4xvc&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=28>`_
* 29 - `How do I create a pandas DataFrame from another object? <https://www.youtube.com/watch?v=-Ov1N1_FbP8&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=29>`_
* 30 - `How do I apply a function to a pandas Series or DataFrame? <https://www.youtube.com/watch?v=P_q0tkYqvSk&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=30>`_
* 31 - `How do I use the MultiIndex in pandas? <https://www.youtube.com/watch?v=tcRGa2soc-c&list=PL5-da3qGB5ICCsgW1MxlZ0Hq8LL5U3u9y&index=31>`_

Lessons for new pandas users
----------------------------

For more resources, please visit the main `repository <https://bitbucket.org/hrojas/learn-pandas>`__.

* `01 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/01%20-%20Lesson.ipynb>`_
    * Importing libraries
    * Creating data sets
    * Creating data frames
    * Reading from CSV
    * Exporting to CSV
    * Finding maximums
    * Plotting data

* `02 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/02%20-%20Lesson.ipynb>`_
    * Reading from TXT
    * Exporting to TXT
    * Selecting top/bottom records
    * Descriptive statistics
    * Grouping/sorting data

* `03 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/03%20-%20Lesson.ipynb>`_
    * Creating functions
    * Reading from EXCEL
    * Exporting to EXCEL
    * Outliers
    * Lambda functions
    * Slice and dice data

* `04 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/04%20-%20Lesson.ipynb>`_
    * Adding/deleting columns
    * Index operations

* `05 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/05%20-%20Lesson.ipynb>`_
    * Stack/Unstack/Transpose functions

* `06 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/06%20-%20Lesson.ipynb>`_
    * GroupBy function

* `07 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/07%20-%20Lesson.ipynb>`_
    * Ways to calculate outliers

* `08 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/08%20-%20Lesson.ipynb>`_
    * Read from Microsoft SQL databases

* `09 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/09%20-%20Lesson.ipynb>`_
    * Export to CSV/EXCEL/TXT

* `10 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/10%20-%20Lesson.ipynb>`_
    * Converting between different kinds of formats

* `11 - Lesson: <http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/11%20-%20Lesson.ipynb>`_
    * Combining data from various sources


Practical data analysis with Python
-----------------------------------

This `guide <http://wavedatalab.github.io/datawithpython>`_ is a comprehensive introduction to the data analysis process using the Python data ecosystem and an interesting open dataset.
There are four sections covering selected topics as follows:

* `Munging Data <http://wavedatalab.github.io/datawithpython/munge.html>`_

* `Aggregating Data <http://wavedatalab.github.io/datawithpython/aggregate.html>`_

* `Visualizing Data <http://wavedatalab.github.io/datawithpython/visualize.html>`_

* `Time Series <http://wavedatalab.github.io/datawithpython/timeseries.html>`_

.. _tutorial-exercises-new-users:

Exercises for new users
-----------------------
Practice your skills with real data sets and exercises.
For more resources, please visit the main `repository <https://github.com/guipsamora/pandas_exercises>`__.

* `01 - Getting & Knowing Your Data <https://github.com/guipsamora/pandas_exercises/tree/master/01_Getting_%26_Knowing_Your_Data>`_

* `02 - Filtering & Sorting <https://github.com/guipsamora/pandas_exercises/tree/master/02_Filtering_%26_Sorting>`_

* `03 - Grouping <https://github.com/guipsamora/pandas_exercises/tree/master/03_Grouping>`_

* `04 - Apply <https://github.com/guipsamora/pandas_exercises/tree/master/04_Apply>`_

* `05 - Merge <https://github.com/guipsamora/pandas_exercises/tree/master/05_Merge>`_

* `06 - Stats <https://github.com/guipsamora/pandas_exercises/tree/master/06_Stats>`_

* `07 - Visualization <https://github.com/guipsamora/pandas_exercises/tree/master/07_Visualization>`_

* `08 - Creating Series and DataFrames <https://github.com/guipsamora/pandas_exercises/tree/master/08_Creating_Series_and_DataFrames/Pokemon>`_

* `09 - Time Series <https://github.com/guipsamora/pandas_exercises/tree/master/09_Time_Series>`_

* `10 - Deleting <https://github.com/guipsamora/pandas_exercises/tree/master/10_Deleting>`_

.. _tutorial-modern:

Modern pandas
-------------

Tutorial series written in 2016 by 
`Tom Augspurger <https://github.com/TomAugspurger>`_.
The source may be found in the GitHub repository
`TomAugspurger/effective-pandas <https://github.com/TomAugspurger/effective-pandas>`_.

* `Modern Pandas <http://tomaugspurger.github.io/modern-1-intro.html>`_
* `Method Chaining <http://tomaugspurger.github.io/method-chaining.html>`_
* `Indexes <http://tomaugspurger.github.io/modern-3-indexes.html>`_
* `Performance <http://tomaugspurger.github.io/modern-4-performance.html>`_
* `Tidy Data <http://tomaugspurger.github.io/modern-5-tidy.html>`_
* `Visualization <http://tomaugspurger.github.io/modern-6-visualization.html>`_
* `Timeseries <http://tomaugspurger.github.io/modern-7-timeseries.html>`_

Excel charts with pandas, vincent and xlsxwriter
------------------------------------------------

*  `Using Pandas and XlsxWriter to create Excel charts <https://pandas-xlsxwriter-charts.readthedocs.io/>`_

Video Tutorials
---------------

* `Pandas From The Ground Up <https://www.youtube.com/watch?v=5JnMutdy6Fw>`_
  (2015) (2:24)
  `GitHub repo <https://github.com/brandon-rhodes/pycon-pandas-tutorial>`__
* `Introduction Into Pandas <https://www.youtube.com/watch?v=-NR-ynQg0YM>`_
  (2016) (1:28)
  `GitHub repo <https://github.com/chendaniely/2016-pydata-carolinas-pandas>`__
* `Pandas: .head() to .tail() <https://www.youtube.com/watch?v=7vuO9QXDN50>`_
  (2016) (1:26)
  `GitHub repo <https://github.com/TomAugspurger/pydata-chi-h2t>`__


Various Tutorials
-----------------

* `Wes McKinney's (pandas BDFL) blog <http://blog.wesmckinney.com/>`_
* `Statistical analysis made easy in Python with SciPy and pandas DataFrames, by Randal Olson <http://www.randalolson.com/2012/08/06/statistical-analysis-made-easy-in-python/>`_
* `Statistical Data Analysis in Python, tutorial videos, by Christopher Fonnesbeck from SciPy 2013 <http://conference.scipy.org/scipy2013/tutorial_detail.php?id=109>`_
* `Financial analysis in Python, by Thomas Wiecki <http://nbviewer.ipython.org/github/twiecki/financial-analysis-python-tutorial/blob/master/1.%20Pandas%20Basics.ipynb>`_
* `Intro to pandas data structures, by Greg Reda <http://www.gregreda.com/2013/10/26/intro-to-pandas-data-structures/>`_
* `Pandas and Python: Top 10, by Manish Amde <http://manishamde.github.io/blog/2013/03/07/pandas-and-python-top-10/>`_
* `Pandas DataFrames Tutorial, by Karlijn Willems <http://www.datacamp.com/community/tutorials/pandas-tutorial-dataframe-python>`_
* `A concise tutorial with real life examples <https://tutswiki.com/pandas-cookbook/chapter1>`_
