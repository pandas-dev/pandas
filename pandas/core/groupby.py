import numpy as np

from cStringIO import StringIO

from pandas.core.frame import DataFrame
from pandas.core.matrix import DataMatrix
from pandas.core.series import Series
from pandas.lib.tseries import isnull

def groupby_withnull(index, mapper):
    index = np.asarray(index)
    mapped_index = np.array([mapper(x) for x in index])

    # A little hack here
    if issubclass(mapped_index.dtype.type, basestring):
        mapped_index = mapped_index.astype(object)

    result = GroupDict()

    mask = isnull(mapped_index)
    nullkeys = index[mask]

    if nullkeys is not None and len(nullkeys) > 0:
        result[np.NaN] = nullkeys

    notmask = -mask
    index = index[notmask]
    mapped_index = mapped_index[notmask]

    for idx, key in zip(index, mapped_index):
        result.setdefault(key, []).append(idx)

    return result

class GroupDict(dict):
    def __repr__(self):
        stringDict = dict((str(x), x) for x in self)
        sortedKeys = sorted(stringDict)

        maxLen = max(len(x) for x in stringDict)

        output = StringIO()
        output.write(str(self.__class__))

        for k in sortedKeys:
            key = stringDict[k]
            size = len(self[key])
            output.write('\n %s -> %d values' % (str(k).ljust(maxLen), size))

        return output.getvalue()

def groupby(obj, grouper):
    """
    Intercepts creation and dispatches to the appropriate class based
    on type.
    """
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataMatrix):
        klass = DataMatrixGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else: # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, grouper)

class GroupBy(object):
    """
    Class for grouping and aggregating relational data.

    Supported classes
    -----------------
    Series / TimeSeries
    DataFrame / DataMatrix (and derivatives thereof)
    """
    _groups = None

    def __init__(self, obj, grouper):
        self.obj = obj
        if hasattr(grouper, 'get'):
            grouper = grouper.get
        self.grouper = grouper

    @property
    def groups(self):
        if self._groups is None:
            self._groups = groupby_withnull(self.obj.index, self.grouper)

        return self._groups

    def getGroup(self, groupList):
        return self.obj.reindex(groupList)

    def __iter__(self):
        """
        Groupby iterator

        Yields
        ------
        Sequence of (groupName, subsetted object) for each group
        """
        try:
            groupNames = sorted(self.groups)
        except Exception: # pragma: no cover
            groupNames = self.groups.keys()

        for groupName in groupNames:
            groupList = self.groups[groupName]
            yield groupName, self.getGroup(groupList)

    def aggregate(self, func):
        raise NotImplementedError

    def agg(self, func):
        return self.aggregate(func)

    def transform(self, func):
        raise NotImplementedError

    def __getitem__(self, key):
        return self.getGroup(self.groups[key])

class SeriesGroupBy(GroupBy):

    def aggregate(self, applyfunc):
        """
        See doc for DataFrame.groupby, group series using mapper (dict or key
        function, apply given function to group, return result as series).

        Main difference here is that applyfunc must return a value, so that the
        result is a sensible series.

        Parameters
        ----------
        mapper : function
            Called on each element of the Series index to determine the groups
        applyfunc : function
            Function to use to aggregate each group

        Returns
        -------
        Series or DataFrame
        """
        from pandas.lib.tseries import groupbyfunc

        if hasattr(applyfunc,'__iter__'):

            if not isinstance(applyfunc, dict):
                applyfunc = dict((func.__name__, func) for func in applyfunc)

            results = {}

            for name, func in applyfunc.iteritems():
                result = self.aggregate(func)
                results[name] = result

            retVal = DataFrame(results)
        else:
            try:
                result = groupbyfunc(self.obj.index, self.obj,
                                     self.grouper, applyfunc)
            except Exception:
                result = {}
                theUnion = set([])
                for groupName, groupList in self.groups.iteritems():
                    groupList = list(groupList)
                    theUnion.update(groupList)
                    grp = self.getGroup(groupList)
                    grp.groupName = groupName
                    output = applyfunc(grp)

                    if isinstance(output, Series):
                        raise Exception('Given applyfunc did not return a '
                                        'value from the subseries as expected!')

                    result[groupName] = output

                for missingIdx in (self.obj.index - list(theUnion)):
                    result[missingIdx] = np.nan

            retVal = Series(result)

        return retVal

    def transform(self, applyfunc):
        """
        For given Series, group index by given mapper function or dict, take
        the sub-Series (reindex) for this group and call apply(applyfunc)
        on this sub-Series. Return a Series of the results for each
        key.

        Parameters
        ----------
        mapper : function
            on being called on each element of the Series
            index, determines the groups.

        applyfunc : function to apply to each group

        Note
        ----
        This function does not aggregate like groupby/tgroupby,
        the results of the given function on the subSeries should be another
        Series.

        Example
        -------
        series.fgroupby(lambda x: mapping[x],
                        lambda x: (x - mean(x)) / std(x))

        Returns
        -------
        Series standardized by each unique value of mapping
        """

        result = {}
        for val, subseries in iter(self):
            subseries.groupName = val
            result[val] = applyfunc(subseries)

        if not isinstance(result.values()[0], Series):
            raise Exception('Given applyfunc did not return a Series from '
                            ' the subseries as expected!')

        allSeries = {}
        for val, subseries in result.iteritems():
            allSeries.update(subseries.iteritems())

        return Series(allSeries)

class DataFrameGroupBy(GroupBy):
    def __init__(self, obj, grouper):
        if isinstance(grouper, basestring) and grouper in obj:
            grouper = obj[grouper].get

        GroupBy.__init__(self, obj, grouper)

    def aggregate(self, applyfunc):
        """
        For given DataFrame, group index by given mapper function or dict, take
        the sub-DataFrame (reindex) for this group and call apply(applyfunc)
        on this sub-DataFrame. Return a DataFrame of the results for each
        key.

        Parameters
        ----------
        mapper : function, dict-like, or string
            Mapping or mapping function. If string given, must be a column
            name in the frame
        applyfunc : function
            Function to use for aggregating groups

        N.B.: applyfunc must produce one value from a Series, otherwise
        an error will occur.

        Optional: provide set mapping as dictionary
        """
        result = {}
        theUnion = set([])
        for val, group in self.groups.iteritems():
            groupList = list(group)
            theUnion.update(groupList)
            result[val] = self.getGroup(groupList).apply(applyfunc)
            assert(isinstance(result[val], Series))
        for missingIdx in (set(self.obj.index) - theUnion):
            result[missingIdx] = np.repeat(np.nan, len(self.obj.cols()))

        return DataFrame(data=result).T

    def transform(self, applyfunc):
        """
        For given DataFrame, group index by given mapper function or dict, take
        the sub-DataFrame (reindex) for this group and call apply(applyfunc)
        on this sub-DataFrame. Return a DataFrame of the results for each
        key.

        Note: this function does not aggregate like groupby/tgroupby,
        the results of the given function on the subDataFrame should be another
        DataFrame.

        Parameters
        ----------
        mapper : function, dict-like, or string
            Mapping or mapping function. If string given, must be a column
            name in the frame
        applyfunc : function
            Function to apply to each subframe

        Note
        ----
        Each subframe is endowed the attribute 'groupName' in case
        you need to know which group you are working on.

        Example
        --------
        df.fgroupby(lambda x: mapping[x],
                    lambda x: (x - x.apply(mean)) / x.apply(std))

            DataFrame standardized by each unique value of mapping
        """
        result = {}
        for val, group in self.groups.iteritems():
            subframe = self.obj.reindex(list(group))
            subframe.groupName = val
            result[val] = applyfunc(subframe).T

        allSeries = {}
        for val, frame in result.iteritems():
            allSeries.update(frame._series)

        return DataFrame(data=allSeries).T


class DataMatrixGroupBy(DataFrameGroupBy):
    def aggregate(self, applyfunc):
        """
        Group series using mapper function (dict or key function, apply
        given function to group, return result as series)

        Main difference here is that applyfunc must return a value, so that the
        result is a sensible series.

        Parameters
        ----------
        mapper : function, dict-like, or string
            Mapping or mapping function. If string given, must be a column
            name in the frame
        applyfunc : function
            Function to use to aggregate each group
        """
        result = {}
        theUnion = set([])
        for val, group in self.groups.iteritems():
            groupList = list(group)
            theUnion.update(groupList)
            result[val] = self.getGroup(groupList).apply(applyfunc)
            assert(isinstance(result[val], Series))
        for missingIdx in (set(self.obj.index) - theUnion):
            result[missingIdx] = np.repeat(np.nan, len(self.obj.cols()))

        return DataMatrix(data=result).T

    def transform(self, applyfunc):
        """
        For given DataMatrix, group index by given mapper function or dict, take
        the sub-DataMatrix (reindex) for this group and call apply(applyfunc)
        on this sub-DataMatrix. Return a DataMatrix of the results for each
        key.

        Note: this function does not aggregate like groupby/tgroupby,
        the results of the given function on the subDataMatrix should be another
        DataMatrix.

        Parameters
        ----------
        mapper : function
            Function to apply to each index value
        applyfunc : function
            Function to apply to each subframe

        Note
        ----
        Each subframe is endowed the attribute 'groupName' in case
        you need to know which group you are working on.

        Example
        --------
        df.fgroupby(lambda x: mapping[x],
                    lambda x: (x - x.apply(mean)) / x.apply(std))

        DataMatrix standardized by each unique value of mapping
        """
        result = {}
        for val, group in self.groups.iteritems():
            subframe = self.obj.reindex(list(group))
            subframe.groupName = val
            result[val] = applyfunc(subframe).T

        allSeries = {}
        for val, frame in result.iteritems():
            allSeries.update(frame._series)

        return DataMatrix(data = allSeries).T
