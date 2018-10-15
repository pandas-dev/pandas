.. _cookbook:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import pandas as pd
   import numpy as np
   from pandas.compat import StringIO

   import random
   import os
   import itertools
   import functools
   import datetime

   np.random.seed(123456)

   pd.options.display.max_rows=15

   import matplotlib
   # matplotlib.style.use('default')

   np.set_printoptions(precision=4, suppress=True)


********
Cookbook
********

This is a repository for *short and sweet* examples and links for useful pandas recipes.
We encourage users to add to this documentation.

Adding interesting links and/or inline examples to this section is a great *First Pull Request*.

Simplified, condensed, new-user friendly, in-line examples have been inserted where possible to
augment the Stack-Overflow and GitHub links.  Many of the links contain expanded information,
above what the in-line examples offer.

Pandas (pd) and Numpy (np) are the only two abbreviated imported modules. The rest are kept
explicitly imported for newer users.

These examples are written for Python 3.  Minor tweaks might be necessary for earlier python
versions.

Idioms
------

.. _cookbook.idioms:

These are some neat pandas ``idioms``

`if-then/if-then-else on one column, and assignment to another one or more columns:
<http://stackoverflow.com/questions/17128302/python-pandas-idiom-for-if-then-else>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

if-then...
**********

An if-then on one column

.. ipython:: python

   df.loc[df.AAA >= 5,'BBB'] = -1; df

An if-then with assignment to 2 columns:

.. ipython:: python

   df.loc[df.AAA >= 5,['BBB','CCC']] = 555; df

Add another line with different logic, to do the -else

.. ipython:: python

   df.loc[df.AAA < 5,['BBB','CCC']] = 2000; df

Or use pandas where after you've set up a mask

.. ipython:: python

   df_mask = pd.DataFrame({'AAA' : [True] * 4, 'BBB' : [False] * 4,'CCC' : [True,False] * 2})
   df.where(df_mask,-1000)

`if-then-else using numpy's where()
<http://stackoverflow.com/questions/19913659/pandas-conditional-creation-of-a-series-dataframe-column>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

   df['logic'] = np.where(df['AAA'] > 5,'high','low'); df

Splitting
*********

`Split a frame with a boolean criterion
<http://stackoverflow.com/questions/14957116/how-to-split-a-dataframe-according-to-a-boolean-criterion>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

   dflow = df[df.AAA <= 5]; dflow
   dfhigh = df[df.AAA > 5]; dfhigh

Building Criteria
*****************

`Select with multi-column criteria
<http://stackoverflow.com/questions/15315452/selecting-with-complex-criteria-from-pandas-dataframe>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

...and (without assignment returns a Series)

.. ipython:: python

   newseries = df.loc[(df['BBB'] < 25) & (df['CCC'] >= -40), 'AAA']; newseries

...or (without assignment returns a Series)

.. ipython:: python

   newseries = df.loc[(df['BBB'] > 25) | (df['CCC'] >= -40), 'AAA']; newseries

...or (with assignment modifies the DataFrame.)

.. ipython:: python

   df.loc[(df['BBB'] > 25) | (df['CCC'] >= 75), 'AAA'] = 0.1; df

`Select rows with data closest to certain value using argsort
<http://stackoverflow.com/questions/17758023/return-rows-in-a-dataframe-closest-to-a-user-defined-number>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

   aValue = 43.0
   df.loc[(df.CCC-aValue).abs().argsort()]

`Dynamically reduce a list of criteria using a binary operators
<http://stackoverflow.com/questions/21058254/pandas-boolean-operation-in-a-python-list/21058331>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

   Crit1 = df.AAA <= 5.5
   Crit2 = df.BBB == 10.0
   Crit3 = df.CCC > -40.0

One could hard code:

.. ipython:: python

   AllCrit = Crit1 & Crit2 & Crit3

...Or it can be done with a list of dynamically built criteria

.. ipython:: python

   CritList = [Crit1,Crit2,Crit3]
   AllCrit = functools.reduce(lambda x,y: x & y, CritList)

   df[AllCrit]

.. _cookbook.selection:

Selection
---------

DataFrames
**********

The :ref:`indexing <indexing>` docs.

`Using both row labels and value conditionals
<http://stackoverflow.com/questions/14725068/pandas-using-row-labels-in-boolean-indexing>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}); df

   df[(df.AAA <= 6) & (df.index.isin([0,2,4]))]

`Use loc for label-oriented slicing and iloc positional slicing
<https://github.com/pandas-dev/pandas/issues/2904>`__

.. ipython:: python

   data = {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40],'CCC' : [100,50,-30,-50]}
   df = pd.DataFrame(data=data,index=['foo','bar','boo','kar']); df

There are 2 explicit slicing methods, with a third general case

1. Positional-oriented (Python slicing style : exclusive of end)
2. Label-oriented (Non-Python slicing style : inclusive of end)
3. General (Either slicing style : depends on if the slice contains labels or positions)

.. ipython:: python
   df.iloc[0:3] #Positional

   df.loc['bar':'kar'] #Label

   # Generic
   df.iloc[0:3]
   df.loc['bar':'kar']

Ambiguity arises when an index consists of integers with a non-zero start or non-unit increment.

.. ipython:: python

   df2 = pd.DataFrame(data=data,index=[1,2,3,4]); #Note index starts at 1.

   df2.iloc[1:3] #Position-oriented

   df2.loc[1:3] #Label-oriented

`Using inverse operator (~) to take the complement of a mask
<http://stackoverflow.com/questions/14986510/picking-out-elements-based-on-complement-of-indices-in-python-pandas>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [4,5,6,7], 'BBB' : [10,20,30,40], 'CCC' : [100,50,-30,-50]}); df

   df[~((df.AAA <= 6) & (df.index.isin([0,2,4])))]

Panels
******

`Extend a panel frame by transposing, adding a new dimension, and transposing back to the original dimensions
<http://stackoverflow.com/questions/15364050/extending-a-pandas-panel-frame-along-the-minor-axis>`__

.. ipython:: python

   rng = pd.date_range('1/1/2013',periods=100,freq='D')
   data = np.random.randn(100, 4)
   cols = ['A','B','C','D']
   df1, df2, df3 = pd.DataFrame(data, rng, cols), pd.DataFrame(data, rng, cols), pd.DataFrame(data, rng, cols)

   pf = pd.Panel({'df1':df1,'df2':df2,'df3':df3});pf

   pf.loc[:,:,'F'] = pd.DataFrame(data, rng, cols);pf

`Mask a panel by using np.where and then reconstructing the panel with the new masked values
<http://stackoverflow.com/questions/14650341/boolean-mask-in-pandas-panel>`__

New Columns
***********

`Efficiently and dynamically creating new columns using applymap
<http://stackoverflow.com/questions/16575868/efficiently-creating-additional-columns-in-a-pandas-dataframe-using-map>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [1,2,1,3], 'BBB' : [1,1,2,2], 'CCC' : [2,1,3,1]}); df

   source_cols = df.columns # or some subset would work too.
   new_cols = [str(x) + "_cat" for x in source_cols]
   categories = {1 : 'Alpha', 2 : 'Beta', 3 : 'Charlie' }

   df[new_cols] = df[source_cols].applymap(categories.get);df

`Keep other columns when using min() with groupby
<http://stackoverflow.com/questions/23394476/keep-other-columns-when-using-min-with-groupby>`__

.. ipython:: python

   df = pd.DataFrame(
        {'AAA' : [1,1,1,2,2,2,3,3], 'BBB' : [2,1,3,4,5,1,2,3]}); df

Method 1 : idxmin() to get the index of the minimums

.. ipython:: python

   df.loc[df.groupby("AAA")["BBB"].idxmin()]

Method 2 : sort then take first of each

.. ipython:: python

   df.sort_values(by="BBB").groupby("AAA", as_index=False).first()

Notice the same results, with the exception of the index.

.. _cookbook.multi_index:

MultiIndexing
-------------

The :ref:`multindexing <advanced.hierarchical>` docs.

`Creating a MultiIndex from a labeled frame
<http://stackoverflow.com/questions/14916358/reshaping-dataframes-in-pandas-based-on-column-labels>`__

.. ipython:: python

   df = pd.DataFrame({'row' : [0,1,2],
                      'One_X' : [1.1,1.1,1.1],
                      'One_Y' : [1.2,1.2,1.2],
                      'Two_X' : [1.11,1.11,1.11],
                      'Two_Y' : [1.22,1.22,1.22]}); df

   # As Labelled Index
   df = df.set_index('row');df
   # With Hierarchical Columns
   df.columns = pd.MultiIndex.from_tuples([tuple(c.split('_')) for c in df.columns]);df
   # Now stack & Reset
   df = df.stack(0).reset_index(1);df
   # And fix the labels (Notice the label 'level_1' got added automatically)
   df.columns = ['Sample','All_X','All_Y'];df

Arithmetic
**********

`Performing arithmetic with a MultiIndex that needs broadcasting
<http://stackoverflow.com/questions/19501510/divide-entire-pandas-multiindex-dataframe-by-dataframe-variable/19502176#19502176>`__

.. ipython:: python

   cols = pd.MultiIndex.from_tuples([ (x,y) for x in ['A','B','C'] for y in ['O','I']])
   df = pd.DataFrame(np.random.randn(2,6),index=['n','m'],columns=cols); df
   df = df.div(df['C'],level=1); df

Slicing
*******

`Slicing a MultiIndex with xs
<http://stackoverflow.com/questions/12590131/how-to-slice-multindex-columns-in-pandas-dataframes>`__

.. ipython:: python

   coords = [('AA','one'),('AA','six'),('BB','one'),('BB','two'),('BB','six')]
   index = pd.MultiIndex.from_tuples(coords)
   df = pd.DataFrame([11,22,33,44,55],index,['MyData']); df

To take the cross section of the 1st level and 1st axis the index:

.. ipython:: python

   df.xs('BB',level=0,axis=0)  #Note : level and axis are optional, and default to zero

...and now the 2nd level of the 1st axis.

.. ipython:: python

   df.xs('six',level=1,axis=0)

`Slicing a MultiIndex with xs, method #2
<http://stackoverflow.com/questions/14964493/multiindex-based-indexing-in-pandas>`__

.. ipython:: python

   index = list(itertools.product(['Ada','Quinn','Violet'],['Comp','Math','Sci']))
   headr = list(itertools.product(['Exams','Labs'],['I','II']))

   indx = pd.MultiIndex.from_tuples(index,names=['Student','Course'])
   cols = pd.MultiIndex.from_tuples(headr) #Notice these are un-named

   data = [[70+x+y+(x*y)%3 for x in range(4)] for y in range(9)]

   df = pd.DataFrame(data,indx,cols); df

   All = slice(None)

   df.loc['Violet']
   df.loc[(All,'Math'),All]
   df.loc[(slice('Ada','Quinn'),'Math'),All]
   df.loc[(All,'Math'),('Exams')]
   df.loc[(All,'Math'),(All,'II')]

`Setting portions of a MultiIndex with xs
<http://stackoverflow.com/questions/19319432/pandas-selecting-a-lower-level-in-a-dataframe-to-do-a-ffill>`__

Sorting
*******

`Sort by specific column or an ordered list of columns, with a MultiIndex
<http://stackoverflow.com/questions/14733871/mutli-index-sorting-in-pandas>`__

.. ipython:: python

   df.sort_values(by=('Labs', 'II'), ascending=False)

`Partial Selection, the need for sortedness;
<https://github.com/pandas-dev/pandas/issues/2995>`__

Levels
******

`Prepending a level to a multiindex
<http://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex>`__

`Flatten Hierarchical columns
<http://stackoverflow.com/questions/14507794/python-pandas-how-to-flatten-a-hierarchical-index-in-columns>`__

.. _cookbook.missing_data:

Missing Data
------------

The :ref:`missing data<missing_data>` docs.

Fill forward a reversed timeseries

.. ipython:: python

   df = pd.DataFrame(np.random.randn(6,1), index=pd.date_range('2013-08-01', periods=6, freq='B'), columns=list('A'))
   df.loc[df.index[3], 'A'] = np.nan
   df
   df.reindex(df.index[::-1]).ffill()

`cumsum reset at NaN values
<http://stackoverflow.com/questions/18196811/cumsum-reset-at-nan>`__

Replace
*******

`Using replace with backrefs
<http://stackoverflow.com/questions/16818871/extracting-value-and-creating-new-column-out-of-it>`__

.. _cookbook.grouping:

Grouping
--------

The :ref:`grouping <groupby>` docs.

`Basic grouping with apply
<http://stackoverflow.com/questions/15322632/python-pandas-df-groupy-agg-column-reference-in-agg>`__

Unlike agg, apply's callable is passed a sub-DataFrame which gives you access to all the columns

.. ipython:: python

   df = pd.DataFrame({'animal': 'cat dog cat fish dog cat cat'.split(),
                      'size': list('SSMMMLL'),
                      'weight': [8, 10, 11, 1, 20, 12, 12],
                      'adult' : [False] * 5 + [True] * 2}); df

   #List the size of the animals with the highest weight.
   df.groupby('animal').apply(lambda subf: subf['size'][subf['weight'].idxmax()])

`Using get_group
<http://stackoverflow.com/questions/14734533/how-to-access-pandas-groupby-dataframe-by-key>`__

.. ipython:: python

   gb = df.groupby(['animal'])

   gb.get_group('cat')

`Apply to different items in a group
<http://stackoverflow.com/questions/15262134/apply-different-functions-to-different-items-in-group-object-python-pandas>`__

.. ipython:: python

   def GrowUp(x):
      avg_weight =  sum(x[x['size'] == 'S'].weight * 1.5)
      avg_weight += sum(x[x['size'] == 'M'].weight * 1.25)
      avg_weight += sum(x[x['size'] == 'L'].weight)
      avg_weight /= len(x)
      return pd.Series(['L',avg_weight,True], index=['size', 'weight', 'adult'])

   expected_df = gb.apply(GrowUp)

   expected_df

`Expanding Apply
<http://stackoverflow.com/questions/14542145/reductions-down-a-column-in-pandas>`__

.. ipython:: python

   S = pd.Series([i / 100.0 for i in range(1,11)])

   def CumRet(x,y):
      return x * (1 + y)

   def Red(x):
      return functools.reduce(CumRet,x,1.0)

   S.expanding().apply(Red, raw=True)


`Replacing some values with mean of the rest of a group
<http://stackoverflow.com/questions/14760757/replacing-values-with-groupby-means>`__

.. ipython:: python

   df = pd.DataFrame({'A' : [1, 1, 2, 2], 'B' : [1, -1, 1, 2]})
   gb = df.groupby('A')

   def replace(g):
       mask = g < 0
       return g.where(mask, g[~mask].mean())

   gb.transform(replace)

`Sort groups by aggregated data
<http://stackoverflow.com/questions/14941366/pandas-sort-by-group-aggregate-and-column>`__

.. ipython:: python

   df = pd.DataFrame({'code': ['foo', 'bar', 'baz'] * 2,
                      'data': [0.16, -0.21, 0.33, 0.45, -0.59, 0.62],
                      'flag': [False, True] * 3})

   code_groups = df.groupby('code')

   agg_n_sort_order = code_groups[['data']].transform(sum).sort_values(by='data')

   sorted_df = df.loc[agg_n_sort_order.index]

   sorted_df

`Create multiple aggregated columns
<http://stackoverflow.com/questions/14897100/create-multiple-columns-in-pandas-aggregation-function>`__

.. ipython:: python

   rng = pd.date_range(start="2014-10-07",periods=10,freq='2min')
   ts = pd.Series(data = list(range(10)), index = rng)

   def MyCust(x):
      if len(x) > 2:
         return x[1] * 1.234
      return pd.NaT

   mhc = {'Mean' : np.mean, 'Max' : np.max, 'Custom' : MyCust}
   ts.resample("5min").apply(mhc)
   ts

`Create a value counts column and reassign back to the DataFrame
<http://stackoverflow.com/questions/17709270/i-want-to-create-a-column-of-value-counts-in-my-pandas-dataframe>`__

.. ipython:: python

   df = pd.DataFrame({'Color': 'Red Red Red Blue'.split(),
                      'Value': [100, 150, 50, 50]}); df
   df['Counts'] = df.groupby(['Color']).transform(len)
   df

`Shift groups of the values in a column based on the index
<http://stackoverflow.com/q/23198053/190597>`__

.. ipython:: python

   df = pd.DataFrame(
      {u'line_race': [10, 10, 8, 10, 10, 8],
       u'beyer': [99, 102, 103, 103, 88, 100]},
       index=[u'Last Gunfighter', u'Last Gunfighter', u'Last Gunfighter',
              u'Paynter', u'Paynter', u'Paynter']); df
   df['beyer_shifted'] = df.groupby(level=0)['beyer'].shift(1)
   df

`Select row with maximum value from each group
<http://stackoverflow.com/q/26701849/190597>`__

.. ipython:: python

   df = pd.DataFrame({'host':['other','other','that','this','this'],
                      'service':['mail','web','mail','mail','web'],
                      'no':[1, 2, 1, 2, 1]}).set_index(['host', 'service'])
   mask = df.groupby(level=0).agg('idxmax')
   df_count = df.loc[mask['no']].reset_index()
   df_count

`Grouping like Python's itertools.groupby
<http://stackoverflow.com/q/29142487/846892>`__

.. ipython:: python

   df = pd.DataFrame([0, 1, 0, 1, 1, 1, 0, 1, 1], columns=['A'])
   df.A.groupby((df.A != df.A.shift()).cumsum()).groups
   df.A.groupby((df.A != df.A.shift()).cumsum()).cumsum()

Expanding Data
**************

`Alignment and to-date
<http://stackoverflow.com/questions/15489011/python-time-series-alignment-and-to-date-functions>`__

`Rolling Computation window based on values instead of counts
<http://stackoverflow.com/questions/14300768/pandas-rolling-computation-with-window-based-on-values-instead-of-counts>`__

`Rolling Mean by Time Interval
<http://stackoverflow.com/questions/15771472/pandas-rolling-mean-by-time-interval>`__

Splitting
*********

`Splitting a frame
<http://stackoverflow.com/questions/13353233/best-way-to-split-a-dataframe-given-an-edge/15449992#15449992>`__

Create a list of dataframes, split using a delineation based on logic included in rows.

.. ipython:: python

   df = pd.DataFrame(data={'Case' : ['A','A','A','B','A','A','B','A','A'],
                           'Data' : np.random.randn(9)})

   dfs = list(zip(*df.groupby((1*(df['Case']=='B')).cumsum().rolling(window=3,min_periods=1).median())))[-1]

   dfs[0]
   dfs[1]
   dfs[2]

.. _cookbook.pivot:

Pivot
*****
The :ref:`Pivot <reshaping.pivot>` docs.

`Partial sums and subtotals
<http://stackoverflow.com/questions/15570099/pandas-pivot-tables-row-subtotals/15574875#15574875>`__

.. ipython:: python

   df = pd.DataFrame(data={'Province' : ['ON','QC','BC','AL','AL','MN','ON'],
                            'City' : ['Toronto','Montreal','Vancouver','Calgary','Edmonton','Winnipeg','Windsor'],
                            'Sales' : [13,6,16,8,4,3,1]})
   table = pd.pivot_table(df,values=['Sales'],index=['Province'],columns=['City'],aggfunc=np.sum,margins=True)
   table.stack('City')

`Frequency table like plyr in R
<http://stackoverflow.com/questions/15589354/frequency-tables-in-pandas-like-plyr-in-r>`__

.. ipython:: python

   grades = [48,99,75,80,42,80,72,68,36,78]
   df = pd.DataFrame( {'ID': ["x%d" % r for r in range(10)],
                       'Gender' : ['F', 'M', 'F', 'M', 'F', 'M', 'F', 'M', 'M', 'M'],
                       'ExamYear': ['2007','2007','2007','2008','2008','2008','2008','2009','2009','2009'],
                       'Class': ['algebra', 'stats', 'bio', 'algebra', 'algebra', 'stats', 'stats', 'algebra', 'bio', 'bio'],
                       'Participated': ['yes','yes','yes','yes','no','yes','yes','yes','yes','yes'],
                       'Passed': ['yes' if x > 50 else 'no' for x in grades],
                       'Employed': [True,True,True,False,False,False,False,True,True,False],
                       'Grade': grades})

   df.groupby('ExamYear').agg({'Participated': lambda x: x.value_counts()['yes'],
                       'Passed': lambda x: sum(x == 'yes'),
                       'Employed' : lambda x : sum(x),
                       'Grade' : lambda x : sum(x) / len(x)})

`Plot pandas DataFrame with year over year data
<http://stackoverflow.com/questions/30379789/plot-pandas-data-frame-with-year-over-year-data>`__

To create year and month cross tabulation:

.. ipython:: python

   df = pd.DataFrame({'value': np.random.randn(36)},
                     index=pd.date_range('2011-01-01', freq='M', periods=36))

   pd.pivot_table(df, index=df.index.month, columns=df.index.year,
                  values='value', aggfunc='sum')

Apply
*****

`Rolling Apply to Organize - Turning embedded lists into a MultiIndex frame
<http://stackoverflow.com/questions/17349981/converting-pandas-dataframe-with-categorical-values-into-binary-values>`__

.. ipython:: python

   df = pd.DataFrame(data={'A' : [[2,4,8,16],[100,200],[10,20,30]], 'B' : [['a','b','c'],['jj','kk'],['ccc']]},index=['I','II','III'])

   def SeriesFromSubList(aList):
      return pd.Series(aList)

   df_orgz = pd.concat(dict([ (ind,row.apply(SeriesFromSubList)) for ind,row in df.iterrows() ]))

`Rolling Apply with a DataFrame returning a Series
<http://stackoverflow.com/questions/19121854/using-rolling-apply-on-a-dataframe-object>`__

Rolling Apply to multiple columns where function calculates a Series before a Scalar from the Series is returned

.. ipython:: python

   df = pd.DataFrame(data=np.random.randn(2000,2)/10000,
                     index=pd.date_range('2001-01-01',periods=2000),
                     columns=['A','B']); df

   def gm(aDF,Const):
      v = ((((aDF.A+aDF.B)+1).cumprod())-1)*Const
      return (aDF.index[0],v.iloc[-1])

   S = pd.Series(dict([ gm(df.iloc[i:min(i+51,len(df)-1)],5) for i in range(len(df)-50) ])); S

`Rolling apply with a DataFrame returning a Scalar
<http://stackoverflow.com/questions/21040766/python-pandas-rolling-apply-two-column-input-into-function/21045831#21045831>`__

Rolling Apply to multiple columns where function returns a Scalar (Volume Weighted Average Price)

.. ipython:: python

   rng = pd.date_range(start = '2014-01-01',periods = 100)
   df = pd.DataFrame({'Open' : np.random.randn(len(rng)),
                      'Close' : np.random.randn(len(rng)),
                      'Volume' : np.random.randint(100,2000,len(rng))}, index=rng); df

   def vwap(bars): return ((bars.Close*bars.Volume).sum()/bars.Volume.sum())
   window = 5
   s = pd.concat([ (pd.Series(vwap(df.iloc[i:i+window]), index=[df.index[i+window]])) for i in range(len(df)-window) ]);
   s.round(2)

Timeseries
----------

`Between times
<http://stackoverflow.com/questions/14539992/pandas-drop-rows-outside-of-time-range>`__

`Using indexer between time
<http://stackoverflow.com/questions/17559885/pandas-dataframe-mask-based-on-index>`__

`Constructing a datetime range that excludes weekends and includes only certain times
<http://stackoverflow.com/questions/24010830/pandas-generate-sequential-timestamp-with-jump/24014440#24014440?>`__

`Vectorized Lookup
<http://stackoverflow.com/questions/13893227/vectorized-look-up-of-values-in-pandas-dataframe>`__

`Aggregation and plotting time series
<http://nipunbatra.github.io/2015/06/timeseries/>`__

Turn a matrix with hours in columns and days in rows into a continuous row sequence in the form of a time series.
`How to rearrange a Python pandas DataFrame?
<http://stackoverflow.com/questions/15432659/how-to-rearrange-a-python-pandas-dataframe>`__

`Dealing with duplicates when reindexing a timeseries to a specified frequency
<http://stackoverflow.com/questions/22244383/pandas-df-refill-adding-two-columns-of-different-shape>`__

Calculate the first day of the month for each entry in a DatetimeIndex

.. ipython:: python

   dates = pd.date_range('2000-01-01', periods=5)
   dates.to_period(freq='M').to_timestamp()

.. _cookbook.resample:

Resampling
**********

The :ref:`Resample <timeseries.resampling>` docs.

`Using Grouper instead of TimeGrouper for time grouping of values
<https://stackoverflow.com/questions/15297053/how-can-i-divide-single-values-of-a-dataframe-by-monthly-averages>`__

`Time grouping with some missing values
<https://stackoverflow.com/questions/33637312/pandas-grouper-by-frequency-with-completeness-requirement>`__

`Valid frequency arguments to Grouper
<http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__

`Grouping using a MultiIndex
<https://stackoverflow.com/questions/41483763/pandas-timegrouper-on-multiindex>`__

`Using TimeGrouper and another grouping to create subgroups, then apply a custom function
<https://github.com/pandas-dev/pandas/issues/3791>`__

`Resampling with custom periods
<http://stackoverflow.com/questions/15408156/resampling-with-custom-periods>`__

`Resample intraday frame without adding new days
<http://stackoverflow.com/questions/14898574/resample-intrday-pandas-dataframe-without-add-new-days>`__

`Resample minute data
<http://stackoverflow.com/questions/14861023/resampling-minute-data>`__

`Resample with groupby <http://stackoverflow.com/q/18677271/564538>`__

.. _cookbook.merge:

Merge
-----

The :ref:`Concat <merging.concatenation>` docs. The :ref:`Join <merging.join>` docs.

`Append two dataframes with overlapping index (emulate R rbind)
<http://stackoverflow.com/questions/14988480/pandas-version-of-rbind>`__

.. ipython:: python

   rng = pd.date_range('2000-01-01', periods=6)
   df1 = pd.DataFrame(np.random.randn(6, 3), index=rng, columns=['A', 'B', 'C'])
   df2 = df1.copy()

Depending on df construction, ``ignore_index`` may be needed

.. ipython:: python

   df = df1.append(df2,ignore_index=True); df

`Self Join of a DataFrame
<https://github.com/pandas-dev/pandas/issues/2996>`__

.. ipython:: python

   df = pd.DataFrame(data={'Area' : ['A'] * 5 + ['C'] * 2,
                           'Bins' : [110] * 2 + [160] * 3 + [40] * 2,
                           'Test_0' : [0, 1, 0, 1, 2, 0, 1],
                           'Data' : np.random.randn(7)});df

   df['Test_1'] = df['Test_0'] - 1

   pd.merge(df, df, left_on=['Bins', 'Area','Test_0'], right_on=['Bins', 'Area','Test_1'],suffixes=('_L','_R'))

`How to set the index and join
<http://stackoverflow.com/questions/14341805/pandas-merge-pd-merge-how-to-set-the-index-and-join>`__

`KDB like asof join
<http://stackoverflow.com/questions/12322289/kdb-like-asof-join-for-timeseries-data-in-pandas/12336039#12336039>`__

`Join with a criteria based on the values
<http://stackoverflow.com/questions/15581829/how-to-perform-an-inner-or-outer-join-of-dataframes-with-pandas-on-non-simplisti>`__

`Using searchsorted to merge based on values inside a range
<http://stackoverflow.com/questions/25125626/pandas-merge-with-logic/2512764>`__

.. _cookbook.plotting:

Plotting
--------

The :ref:`Plotting <visualization>` docs.

`Make Matplotlib look like R
<http://stackoverflow.com/questions/14349055/making-matplotlib-graphs-look-like-r-by-default>`__

`Setting x-axis major and minor labels
<http://stackoverflow.com/questions/12945971/pandas-timeseries-plot-setting-x-axis-major-and-minor-ticks-and-labels>`__

`Plotting multiple charts in an ipython notebook
<http://stackoverflow.com/questions/16392921/make-more-than-one-chart-in-same-ipython-notebook-cell>`__

`Creating a multi-line plot
<http://stackoverflow.com/questions/16568964/make-a-multiline-plot-from-csv-file-in-matplotlib>`__

`Plotting a heatmap
<http://stackoverflow.com/questions/17050202/plot-timeseries-of-histograms-in-python>`__

`Annotate a time-series plot
<http://stackoverflow.com/questions/11067368/annotate-time-series-plot-in-matplotlib>`__

`Annotate a time-series plot #2
<http://stackoverflow.com/questions/17891493/annotating-points-from-a-pandas-dataframe-in-matplotlib-plot>`__

`Generate Embedded plots in excel files using Pandas, Vincent and xlsxwriter
<https://pandas-xlsxwriter-charts.readthedocs.io/>`__

`Boxplot for each quartile of a stratifying variable
<http://stackoverflow.com/questions/23232989/boxplot-stratified-by-column-in-python-pandas>`__

.. ipython:: python

   df = pd.DataFrame(
        {u'stratifying_var': np.random.uniform(0, 100, 20),
         u'price': np.random.normal(100, 5, 20)})

   df[u'quartiles'] = pd.qcut(
       df[u'stratifying_var'],
       4,
       labels=[u'0-25%', u'25-50%', u'50-75%', u'75-100%'])

   @savefig quartile_boxplot.png
   df.boxplot(column=u'price', by=u'quartiles')

Data In/Out
-----------

`Performance comparison of SQL vs HDF5
<http://stackoverflow.com/questions/16628329/hdf5-and-sqlite-concurrency-compression-i-o-performance>`__

.. _cookbook.csv:

CSV
***

The :ref:`CSV <io.read_csv_table>` docs

`read_csv in action <http://wesmckinney.com/blog/update-on-upcoming-pandas-v0-10-new-file-parser-other-performance-wins/>`__

`appending to a csv
<http://stackoverflow.com/questions/17134942/pandas-dataframe-output-end-of-csv>`__

`Reading a csv chunk-by-chunk
<http://stackoverflow.com/questions/11622652/large-persistent-dataframe-in-pandas/12193309#12193309>`__

`Reading only certain rows of a csv chunk-by-chunk
<http://stackoverflow.com/questions/19674212/pandas-data-frame-select-rows-and-clear-memory>`__

`Reading the first few lines of a frame
<http://stackoverflow.com/questions/15008970/way-to-read-first-few-lines-for-pandas-dataframe>`__

Reading a file that is compressed but not by ``gzip/bz2`` (the native compressed formats which ``read_csv`` understands).
This example shows a ``WinZipped`` file, but is a general application of opening the file within a context manager and
using that handle to read.
`See here
<http://stackoverflow.com/questions/17789907/pandas-convert-winzipped-csv-file-to-data-frame>`__

`Inferring dtypes from a file
<http://stackoverflow.com/questions/15555005/get-inferred-dataframe-types-iteratively-using-chunksize>`__

`Dealing with bad lines
<http://github.com/pandas-dev/pandas/issues/2886>`__

`Dealing with bad lines II
<http://nipunbatra.github.io/2013/06/reading-unclean-data-csv-using-pandas/>`__

`Reading CSV with Unix timestamps and converting to local timezone
<http://nipunbatra.github.io/2013/06/pandas-reading-csv-with-unix-timestamps-and-converting-to-local-timezone/>`__

`Write a multi-row index CSV without writing duplicates
<http://stackoverflow.com/questions/17349574/pandas-write-multiindex-rows-with-to-csv>`__

.. _cookbook.csv.multiple_files:

Reading multiple files to create a single DataFrame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The best way to combine multiple files into a single DataFrame is to read the individual frames one by one, put all
of the individual frames into a list, and then combine the frames in the list using :func:`pd.concat`:

.. ipython:: python

    for i in range(3):
        data = pd.DataFrame(np.random.randn(10, 4))
        data.to_csv('file_{}.csv'.format(i))

    files = ['file_0.csv', 'file_1.csv', 'file_2.csv']
    result = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

You can use the same approach to read all files matching a pattern.  Here is an example using ``glob``:

.. ipython:: python

    import glob
    files = glob.glob('file_*.csv')
    result = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

Finally, this strategy will work with the other ``pd.read_*(...)`` functions described in the :ref:`io docs<io>`.

.. ipython:: python
    :suppress:

    for i in range(3):
        os.remove('file_{}.csv'.format(i))

Parsing date components in multi-columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parsing date components in multi-columns is faster with a format

.. code-block:: python

    In [30]: i = pd.date_range('20000101',periods=10000)

    In [31]: df = pd.DataFrame(dict(year = i.year, month = i.month, day = i.day))

    In [32]: df.head()
    Out[32]:
       day  month  year
    0    1      1  2000
    1    2      1  2000
    2    3      1  2000
    3    4      1  2000
    4    5      1  2000

    In [33]: %timeit pd.to_datetime(df.year*10000+df.month*100+df.day,format='%Y%m%d')
    100 loops, best of 3: 7.08 ms per loop

    # simulate combinging into a string, then parsing
    In [34]: ds = df.apply(lambda x: "%04d%02d%02d" % (x['year'],x['month'],x['day']),axis=1)

    In [35]: ds.head()
    Out[35]:
    0    20000101
    1    20000102
    2    20000103
    3    20000104
    4    20000105
    dtype: object

    In [36]: %timeit pd.to_datetime(ds)
    1 loops, best of 3: 488 ms per loop

Skip row between header and data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ipython:: python

    data = """;;;;
     ;;;;
     ;;;;
     ;;;;
     ;;;;
     ;;;;
    ;;;;
     ;;;;
     ;;;;
    ;;;;
    date;Param1;Param2;Param4;Param5
        ;m²;°C;m²;m
    ;;;;
    01.01.1990 00:00;1;1;2;3
    01.01.1990 01:00;5;3;4;5
    01.01.1990 02:00;9;5;6;7
    01.01.1990 03:00;13;7;8;9
    01.01.1990 04:00;17;9;10;11
    01.01.1990 05:00;21;11;12;13
    """

Option 1: pass rows explicitly to skip rows
"""""""""""""""""""""""""""""""""""""""""""

.. ipython:: python

    pd.read_csv(StringIO(data), sep=';', skiprows=[11,12],
            index_col=0, parse_dates=True, header=10)

Option 2: read column names and then data
"""""""""""""""""""""""""""""""""""""""""

.. ipython:: python

    pd.read_csv(StringIO(data), sep=';', header=10, nrows=10).columns
    columns = pd.read_csv(StringIO(data), sep=';', header=10, nrows=10).columns
    pd.read_csv(StringIO(data), sep=';', index_col=0,
                header=12, parse_dates=True, names=columns)


.. _cookbook.sql:

SQL
***

The :ref:`SQL <io.sql>` docs

`Reading from databases with SQL
<http://stackoverflow.com/questions/10065051/python-pandas-and-databases-like-mysql>`__

.. _cookbook.excel:

Excel
*****

The :ref:`Excel <io.excel>` docs

`Reading from a filelike handle
<http://stackoverflow.com/questions/15588713/sheets-of-excel-workbook-from-a-url-into-a-pandas-dataframe>`__

`Modifying formatting in XlsxWriter output
<http://pbpython.com/improve-pandas-excel-output.html>`__

.. _cookbook.html:

HTML
****

`Reading HTML tables from a server that cannot handle the default request
header <http://stackoverflow.com/a/18939272/564538>`__

.. _cookbook.hdf:

HDFStore
********

The :ref:`HDFStores <io.hdf5>` docs

`Simple Queries with a Timestamp Index
<http://stackoverflow.com/questions/13926089/selecting-columns-from-pandas-hdfstore-table>`__

`Managing heterogeneous data using a linked multiple table hierarchy
<http://github.com/pandas-dev/pandas/issues/3032>`__

`Merging on-disk tables with millions of rows
<http://stackoverflow.com/questions/14614512/merging-two-tables-with-millions-of-rows-in-python/14617925#14617925>`__

`Avoiding inconsistencies when writing to a store from multiple processes/threads
<http://stackoverflow.com/a/29014295/2858145>`__

De-duplicating a large store by chunks, essentially a recursive reduction operation. Shows a function for taking in data from
csv file and creating a store by chunks, with date parsing as well.
`See here
<http://stackoverflow.com/questions/16110252/need-to-compare-very-large-files-around-1-5gb-in-python/16110391#16110391>`__

`Creating a store chunk-by-chunk from a csv file
<http://stackoverflow.com/questions/20428355/appending-column-to-frame-of-hdf-file-in-pandas/20428786#20428786>`__

`Appending to a store, while creating a unique index
<http://stackoverflow.com/questions/16997048/how-does-one-append-large-amounts-of-data-to-a-pandas-hdfstore-and-get-a-natural/16999397#16999397>`__

`Large Data work flows
<http://stackoverflow.com/questions/14262433/large-data-work-flows-using-pandas>`__

`Reading in a sequence of files, then providing a global unique index to a store while appending
<http://stackoverflow.com/questions/16997048/how-does-one-append-large-amounts-of-data-to-a-pandas-hdfstore-and-get-a-natural>`__

`Groupby on a HDFStore with low group density
<http://stackoverflow.com/questions/15798209/pandas-group-by-query-on-large-data-in-hdfstore>`__

`Groupby on a HDFStore with high group density
<http://stackoverflow.com/questions/25459982/trouble-with-grouby-on-millions-of-keys-on-a-chunked-file-in-python-pandas/25471765#25471765>`__

`Hierarchical queries on a HDFStore
<http://stackoverflow.com/questions/22777284/improve-query-performance-from-a-large-hdfstore-table-with-pandas/22820780#22820780>`__

`Counting with a HDFStore
<http://stackoverflow.com/questions/20497897/converting-dict-of-dicts-into-pandas-dataframe-memory-issues>`__

`Troubleshoot HDFStore exceptions
<http://stackoverflow.com/questions/15488809/how-to-trouble-shoot-hdfstore-exception-cannot-find-the-correct-atom-type>`__

`Setting min_itemsize with strings
<http://stackoverflow.com/questions/15988871/hdfstore-appendstring-dataframe-fails-when-string-column-contents-are-longer>`__

`Using ptrepack to create a completely-sorted-index on a store
<http://stackoverflow.com/questions/17893370/ptrepack-sortby-needs-full-index>`__

Storing Attributes to a group node

.. ipython:: python

   df = pd.DataFrame(np.random.randn(8,3))
   store = pd.HDFStore('test.h5')
   store.put('df',df)

   # you can store an arbitrary Python object via pickle
   store.get_storer('df').attrs.my_attribute = dict(A = 10)
   store.get_storer('df').attrs.my_attribute

.. ipython:: python
   :suppress:

   store.close()
   os.remove('test.h5')

.. _cookbook.binary:

Binary Files
************

pandas readily accepts NumPy record arrays, if you need to read in a binary
file consisting of an array of C structs. For example, given this C program
in a file called ``main.c`` compiled with ``gcc main.c -std=gnu99`` on a
64-bit machine,

.. code-block:: c

   #include <stdio.h>
   #include <stdint.h>

   typedef struct _Data
   {
       int32_t count;
       double avg;
       float scale;
   } Data;

   int main(int argc, const char *argv[])
   {
       size_t n = 10;
       Data d[n];

       for (int i = 0; i < n; ++i)
       {
           d[i].count = i;
           d[i].avg = i + 1.0;
           d[i].scale = (float) i + 2.0f;
       }

       FILE *file = fopen("binary.dat", "wb");
       fwrite(&d, sizeof(Data), n, file);
       fclose(file);

       return 0;
   }

the following Python code will read the binary file ``'binary.dat'`` into a
pandas ``DataFrame``, where each element of the struct corresponds to a column
in the frame:

.. code-block:: python

   names = 'count', 'avg', 'scale'

   # note that the offsets are larger than the size of the type because of
   # struct padding
   offsets = 0, 8, 16
   formats = 'i4', 'f8', 'f4'
   dt = np.dtype({'names': names, 'offsets': offsets, 'formats': formats},
                 align=True)
   df = pd.DataFrame(np.fromfile('binary.dat', dt))

.. note::

   The offsets of the structure elements may be different depending on the
   architecture of the machine on which the file was created. Using a raw
   binary file format like this for general data storage is not recommended, as
   it is not cross platform. We recommended either HDF5 or msgpack, both of
   which are supported by pandas' IO facilities.

Computation
-----------

`Numerical integration (sample-based) of a time series
<http://nbviewer.ipython.org/5720498>`__

Correlation
***********

The `method` argument within `DataFrame.corr` can accept a callable in addition to the named correlation types.  Here we compute the `distance correlation <https://en.wikipedia.org/wiki/Distance_correlation>`__ matrix for a `DataFrame` object.

.. ipython:: python

    def distcorr(x, y):
        n = len(x)
        a = np.zeros(shape=(n, n))
        b = np.zeros(shape=(n, n))

        for i in range(n):
            for j in range(i + 1, n):
                a[i, j] = abs(x[i] - x[j])
                b[i, j] = abs(y[i] - y[j])

        a += a.T
        b += b.T

        a_bar = np.vstack([np.nanmean(a, axis=0)] * n)
        b_bar = np.vstack([np.nanmean(b, axis=0)] * n)

        A = a - a_bar - a_bar.T + np.full(shape=(n, n), fill_value=a_bar.mean())
        B = b - b_bar - b_bar.T + np.full(shape=(n, n), fill_value=b_bar.mean())

        cov_ab = np.sqrt(np.nansum(A * B)) / n
        std_a = np.sqrt(np.sqrt(np.nansum(A**2)) / n)
        std_b = np.sqrt(np.sqrt(np.nansum(B**2)) / n)

        return cov_ab / std_a / std_b

    df = pd.DataFrame(np.random.normal(size=(100, 3)))

    df.corr(method=distcorr)

Timedeltas
----------

The :ref:`Timedeltas <timedeltas.timedeltas>` docs.

`Using timedeltas
<http://github.com/pandas-dev/pandas/pull/2899>`__

.. ipython:: python

   s  = pd.Series(pd.date_range('2012-1-1', periods=3, freq='D'))

   s - s.max()

   s.max() - s

   s - datetime.datetime(2011,1,1,3,5)

   s + datetime.timedelta(minutes=5)

   datetime.datetime(2011,1,1,3,5) - s

   datetime.timedelta(minutes=5) + s

`Adding and subtracting deltas and dates
<http://stackoverflow.com/questions/16385785/add-days-to-dates-in-dataframe>`__

.. ipython:: python

   deltas = pd.Series([ datetime.timedelta(days=i) for i in range(3) ])

   df = pd.DataFrame(dict(A = s, B = deltas)); df

   df['New Dates'] = df['A'] + df['B'];

   df['Delta'] = df['A'] - df['New Dates']; df

   df.dtypes

`Another example
<http://stackoverflow.com/questions/15683588/iterating-through-a-pandas-dataframe>`__

Values can be set to NaT using np.nan, similar to datetime

.. ipython:: python

   y = s - s.shift(); y

   y[1] = np.nan; y

Aliasing Axis Names
-------------------

To globally provide aliases for axis names, one can define these 2 functions:

.. ipython:: python

   def set_axis_alias(cls, axis, alias):
      if axis not in cls._AXIS_NUMBERS:
         raise Exception("invalid axis [%s] for alias [%s]" % (axis, alias))
      cls._AXIS_ALIASES[alias] = axis

.. ipython:: python

   def clear_axis_alias(cls, axis, alias):
      if axis not in cls._AXIS_NUMBERS:
         raise Exception("invalid axis [%s] for alias [%s]" % (axis, alias))
      cls._AXIS_ALIASES.pop(alias,None)

.. ipython:: python

   set_axis_alias(pd.DataFrame,'columns', 'myaxis2')
   df2 = pd.DataFrame(np.random.randn(3,2),columns=['c1','c2'],index=['i1','i2','i3'])
   df2.sum(axis='myaxis2')
   clear_axis_alias(pd.DataFrame,'columns', 'myaxis2')

Creating Example Data
---------------------

To create a dataframe from every combination of some given values, like R's ``expand.grid()``
function, we can create a dict where the keys are column names and the values are lists
of the data values:

.. ipython:: python


   def expand_grid(data_dict):
      rows = itertools.product(*data_dict.values())
      return pd.DataFrame.from_records(rows, columns=data_dict.keys())

   df = expand_grid(
      {'height': [60, 70],
       'weight': [100, 140, 180],
       'sex': ['Male', 'Female']})
   df
