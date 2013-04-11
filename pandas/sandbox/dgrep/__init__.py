
# TODO, add axis argument
def dgrep(self,pred,cols=None,C=0,B=0,A=0,split=False,keys=True):
    """Select rows by regex match or predicate function, against *data*.

    This is an unindexed operation, and is substantially slower then
    index-based selection for large datasets.

    cols : string or sequence of str
             name or sequence of column named if running against a DataFrame,
             ignored otherwise.
    pred : string regex or f(val) -> bool or value to test equality against.

             if the predicate function expects *args or multiple unnamed
             arguments, the row values for the specified columns will be passed
             in to the the predicate function as a list, one call per row.

    A/B,C : int, grep-like argument, context lines (A)fter/(B)efore or  (C)entered (C)ontext
    split: bool , False returns a  slice of the current object, if context lines overlap
                between matches, they will only appear once. a True value will return
                a list of frames or (matched_index_label, self_sliced) pairs (default),
                depending on the the value of `keys`. Similar to the groupby API.
    keys: bool, if split==True, keys=False will make the function return
             a list of frames, rather then a list of (label, dataframe) pairs.

    Usage examples:

    from pandas.util.testing import makeCustomDataframe as mkdf

    df=mkdf(30,4,r_idx_nlevels=3)
    df.index=range(30)
    df.iloc[5,0] = "supercool"
    df.iloc[6,0] = "supercool"
    df.iloc[29,0] = "supercool"
    df.iloc[15,1] = "supercool"
    df.iloc[17,2] = "supercool"
    # accepts colname and regex string
    df.dgrep(".cool$","C_l0_g0")

    df.dgrep(".cool$",["C_l0_g0","C_l0_g1"])
    # specifying C=2 (or A/B=) does a grep context , providing
    # context lines around the hit
    # NB overlapping context lines do not cause line duplication (*)
    df.dgrep(".cool$",["C_l0_g0"],C=2)
    # also accepts lambda
    # NB, last match is at end, so only previous line of context displayed
    df.dgrep(lambda x: bool(re.search(".cool$",x)),["C_l0_g0"],C=3)
    # split=True returns a series of (index_label_matched, dataframe)
    # pairs, similar to groupby
    # NB some lines appear in more then one group in this case (*)
    df.dgrep(".cool$",["C_l0_g0"],split=True,C=3)

    # works on series too
    df.C_l0_g0.dgrep(".cool$",C=3)

    # can also get the values "applied" onto the function
    df.dgrep(lambda c1,c2: "cool" in c1 or "cool" in c2,df.columns[:2])

    # which also works with *args
    df.dgrep(lambda *args: "supercool" in args,df.columns[:3])
    """
    from pandas import DataFrame
    from pandas.core.common import _is_sequence
    import inspect

    if _is_sequence(cols):
        cols = list(cols)   # convert index to list, from slice such as df.columns[:3]
    if not isinstance(cols,(list,tuple)):
        cols = [cols]

    combine=False
    if callable(pred):
        fargs=inspect.getargspec(pred)
        if fargs.varargs:
           combine=True

        elif len(fargs.args) > 1:
           if len(fargs.args) !=  len(cols):
               raise ValueError("predicate function argcount doesn't match num. of cols")
           combine=True

    elif isinstance(pred,basestring):
        _pat = pred
        def f1(x):
            import re
            return bool(re.search(_pat,unicode(x)))
        pred=f1
    else:
        def f2(x):
            return x == pred
        pred=f2

    indicies =  set()
    if isinstance(self,DataFrame):
        if  combine:
            # print( [ self.irow(i).ix[cols] for i in range(len(self)) ])
            indicies.update([ i for i,x in enumerate(self.index) if pred(*self.irow(i).ix[cols])])
        else:
            for name in cols:
                s = self[name]
                indicies.update([ i for i,x in enumerate(s)
                                if pred(x)])
    else:
            indicies.update([ i for i,x in enumerate(self)
                                if pred(x)])

    return self.neighbours(self.index[list(sorted(indicies))],C=C,B=B,A=A,split=split,keys=keys)

def neighbours(self,labels,C=None,B=None,A=None,split=False,keys=True):
    """Takes a list of labels and return one ore more frame/series with the indicated
    rows + surrounding rows as determined by the (A)fter/(B)efore or
    (C)entered (C)ontext.

    see the `dgrep` docstring for more details about the identical arguments.

    """
    if C:
        B = C//2
        A = C-B-1

    indicies = map(self.index.get_loc,labels)
    if split:
        #list of (hit_label,sliced frame)
        def g(x):
            return (x,range(max(0,x-B),min(x+A+1,len(self.index))))

        indicies_grps = map(g,indicies)
        results = []
        for i,indicies in indicies_grps:
            if keys:
                results.append((self.index[i],self.iloc[indicies]))
            else:
                results.append(self.iloc[indicies])
        return results
    else:
        indicies=reduce(lambda acc,x: acc+range(max(0,x-B),min(x+A+1,len(self.index))),
                    indicies,[])
        # there's just one, and return just the sliced frame, not the hit label
        return self.iloc[sorted(set(indicies))]
