''' Provides composition class, MetaDataFrame, which is an ordinary python object that stores a Dataframe and 
attempts to promote attributes and methods to the instance level (eg self.x instead of self.df.x).  This object
can be subclassed and ensures persistence of custom attributes.  The goal of this MetaDataFrame is to provide a 
subclassing api beyond monkey patching (which currently fails in persisting attributes upon most method returns 
and upon derialization.'''

from types import MethodType
import copy
import functools
import cPickle
import collections

from pandas.core.indexing import _NDFrameIndexer

from pandas import DataFrame, Series

## for testing
from numpy.random import randn

#----------------------------------------------------------------------
# Store attributes/methods of dataframe for later inspection with __setattr__
# Note: This is preferred to a storing individual instances of self._df with custom 
#       attr as if user tried self.a and self._df.a existed, it would call this...
_dfattrs=[x for x in dir(DataFrame) if '__' not in x]

#----------------------------------------------------------------------
# Loading (perhaps change name?) ... Doesn't work correctly as instance methods

def mload(inname):
    ''' Load MetaDataFrame from file'''
    if isinstance(inname, basestring):
        inname=open(inname, 'r')
    return cPickle.load(inname)

def mloads(string):
    ''' Load a MetaDataFrame from string stored in memory.'''
    return cPickle.loads(string)        
        

class MetaDataFrame(object):
    ''' Base composition for subclassing dataframe.'''

    def __init__(self, *dfargs, **dfkwargs):
        ''' Stores a dataframe under reserved attribute name, self._df'''      
        self._df=DataFrame(*dfargs, **dfkwargs)
                
    ### Save methods    
    def save(self, outname):
        ''' Takes in str or opened file and saves. cPickle.dump wrapper.'''
        if isinstance(outname, basestring):
            outname=open(outname, 'w')
        cPickle.dump(self, outname)

    def dumps(self):
        ''' Output TimeSpectra into a pickled string in memory.'''
        return cPickle.dumps(self)

    def deepcopy(self):
        ''' Make a deepcopy of self, including the dataframe.'''
        return copy.deepcopy(self)   

    def as_dataframe(self):
        ''' Convience method to return a raw dataframe, self._df'''
        return self._df    

    #----------------------------------------------------------------------
    # Overwrite Dataframe methods and operators

    def __getitem__(self, key):
        ''' Item lookup.  If output is an interable, _transfer is called.  
        Sometimes __getitem__ returns a float (indexing a series) at which 
        point we just want to return that.'''

        dfout=self._df.__getitem__(key)

        try:
            iter(dfout)  #Test if iterable without forcing user to have collections package.
        except TypeError:
            return dfout
        else:
            return self._transfer(self._df.__getitem__(key) )               

    def __setitem__(self, key, value):
        self._df.__setitem__(key, value)    

    ### These tell python to ignore __getattr__ when pickling; hence, treat this like a normal class    
    def __getstate__(self): return self.__dict__
    def __setstate__(self, d): self.__dict__.update(d)    

    def __getattr__(self, attr, *fcnargs, **fcnkwargs):
        ''' Tells python how to handle all attributes that are not found.  Basic attributes 
        are directly referenced to self._df; however, instance methods (like df.corr() ) are
        handled specially using a special private parsing method, _dfgetattr().'''

        ### Return basic attribute
        
        try:
            refout=getattr(self._df, attr)
        except AttributeError:
            raise AttributeError('Could not find attribute "%s" in %s or its underlying DataFrame'%(attr, self.__class__.__name__))           
           
        if not isinstance(refout, MethodType):
            return refout

        ### Handle instance methods using _dfgetattr().
        ### see http://stackoverflow.com/questions/3434938/python-allowing-methods-not-specifically-defined-to-be-called-ala-getattr
        else:         
            return functools.partial(self._dfgetattr, attr, *fcnargs, **fcnkwargs)
            ### This is a reference to the fuction (aka a wrapper) not the function itself
            
    def __setattr__(self, name, value):
        ''' When user sets an attribute, this tries to intercept any name conflicts.  For example, if user attempts to set
        self.columns=50, this will actually try self._df.columns=50, which throws an error.  The behavior is acheived by
        using dir() on the data frame created upon initialization, filtering __x__ type methods.   Not guaranteed to work 100%
        of the time due to implicit possible issues with dir() and inspection in Python.  Best practice is for users to avoid name
        conflicts when possible.'''
        
        super(MetaDataFrame, self).__setattr__(name, value)        
        if name in _dfattrs:
            setattr(self._df, name, value)
        else:
            self.__dict__[name]=value


    def _transfer(self, dfnew):
        ''' Copies all attribtues into a new object except has to store current dataframe
        in memory as this can't be copied correctly using copy.deepcopy.  Probably a quicker way...

        dfnew is used if one wants to pass a new dataframe in.  This is used primarily in calls from __getattr__.'''
        ### Store old value of df and remove current df to copy operation will take
        olddf=self._df.copy() #Removed deep=True because series return could not implement it
        self._df=None

        ### Create new object and apply new df 
        newobj=copy.deepcopy(self)  #This looks like None, but is it type (MetaDataFrame, just __union__ prints None
        newobj._df=dfnew

        ### Restore old value of df and return new object
        self._df=olddf
        return newobj


    def _dfgetattr(self, attr, *fcnargs, **fcnkwargs):
        ''' Called by __getattr__ as a wrapper, this private method is used to ensure that any
        DataFrame method that returns a new DataFrame will actually return a TimeSpectra object
        instead.  It does so by typechecking the return of attr().

        **kwargs: use_base - If true, program attempts to call attribute on the baseline.  Baseline ought
        to be maintained as a series, and Series/Dataframe API's must be same.

        *fcnargs and **fcnkwargs are passed to the dataframe method.

        Note: tried to ad an as_new keyword to do this operation in place, but doing self=dfout instead of return dfout
        didn't work.  Could try to add this at the __getattr__ level; however, may not be worth it.'''

        out=getattr(self._df, attr)(*fcnargs, **fcnkwargs)

        ### If operation returns a dataframe, return new TimeSpectra
        if isinstance(out, DataFrame):
            dfout=self._transfer(out)
            return dfout

        ### Otherwise return whatever the method return would be
        else:
            return out

    def __repr__(self):
        return self._df.__repr__()

    ### Operator overloading ####
    ### In place operations need to overwrite self._df
    def __add__(self, x):
        return self._transfer(self._df.__add__(x))

    def __sub__(self, x):
        return self._transfer(self._df.__sub__(x))

    def __mul__(self, x):
        return self._transfer(self._df.__mul__(x))

    def __div__(self, x):
        return self._transfer(self._df.__div__(x))

    def __truediv__(self, x):
        return self._transfer(self._df.__truediv__(x))

    ### From what I can tell, __pos__(), __abs__() builtin to df, just __neg__()    
    def __neg__(self):  
        return self._transfer(self._df.__neg__() )

    ### Object comparison operators
    def __lt__(self, x):
        return self._transfer(self._df.__lt__(x))

    def __le__(self, x):
        return self._transfer(self._df.__le__(x))

    def __eq__(self, x):
        return self._transfer(self._df.__eq__(x))

    def __ne__(self, x):
        return self._transfer(self._df.__ne__(x))

    def __ge__(self, x):
        return self._transfer(self._df.__ge__(x))

    def __gt__(self, x):
        return self._transfer(self._df.__gt__(x))     

    def __len__(self):
        return self._df.__len__()

    def __nonzero__(self):
        return self._df.__nonzero__()

    def __contains__(self, x):
        return self._df.__contains__(x)

    def __iter__(self):
        return self._df.__iter__()


    ## Fancy indexing
    _ix=None     
        
    @property	  	
    def ix(self, *args, **kwargs):      	
        ''' Pandas Indexing.  Note, this has been modified to ensure that series returns (eg ix[3])
        still maintain attributes.  To remove this behavior, replace the following:
        
        self._ix = _MetaIndexer(self, _NDFrameIndexer(self) ) --> self._ix=_NDFrameIndexer(self)
        
        The above works because slicing preserved attributes because the _NDFrameIndexer is a python object 
        subclass.'''
        if self._ix is None:
            self._ix=_MetaIndexer(self)
        return self._ix        
            
class _MetaIndexer(_NDFrameIndexer):
    ''' Intercepts the slicing of ix so Series returns can be handled properly.  In addition,
        it makes sure that the new index is assigned properly.
        
        Notes:
        -----
          Under the hood pandas called upon _NDFrameIndexer methods, so this merely overwrites the
          ___getitem__() method and leaves all the rest intact'''
    
    def __getitem__(self, key):
        out=super(_MetaIndexer, self).__getitem__(key)   

        ### Series returns transformed to MetaDataFrame
        if isinstance(out, Series):
            df=DataFrame(out)
            return self.obj._transfer(out)

        ### Make sure the new object's index property is syched to its ._df index.
        else:
            return out  
    


class SubFoo(MetaDataFrame):
    ''' Shows an example of how to subclass MetaDataFrame with custom attributes, a and b.'''

    def __init__(self, a, b, *dfargs, **dfkwargs):
        self.a = a
        self.b = b    

        super(SubFoo, self).__init__(*dfargs, **dfkwargs)

    def __repr__(self):
        return "Hi I'm SubFoo. I'm not really a DataFrame, but I quack like one."

    @property
    def data(self):
        ''' Return underyling dataframe attribute self._df'''
        return self._data


#### TESTING ###
if __name__ == '__main__':

    ### Create a MetaDataFrame
    meta_df=MetaDataFrame(abs(randn(3,3)), index=['A','B','C'], columns=['c11','c22', 'c33'])    
    
    meta_df.to_csv('deletejunkme')

    ### Add some new attributes
    meta_df.a=50
    meta_df.b='Pamela'
    print 'See the original metadataframe\n'
    print meta_df
    print '\nI can operate on it (+ - / *) and call dataframe methods like rank()'
    
    meta_df.ix[0]

    ### Perform some intrinsic DF operations
    new=meta_df*50.0
    new=new.rank()
    print '\nSee modified dataframe:\n'
    print new

    ### Verify attribute persistence
    print '\nAttributes a = %s and b = %s will persist when new metadataframes are returned.'%(new.a, new.b)

    ### Demonstrate subclassing by invoking SubFoo class
    print '\nI can subclass a dataframe an overwrite its __repr__() or more carefully __bytes__()/__unicode__() method(s)\n'
    subclass=SubFoo(50, 200, abs(randn(3,3)), index=['A','B','C'], columns=['c11','c22', 'c33'])    
    print subclass
    ### Access underlying dataframe
    print '\nMy underlying dataframe is stored in the "data" attribute.\n'
    print subclass.data

    ### Pickle
    print '\nSave me by using x.save() / x.dumps() and load using mload(x) / mloads(x).'
#    df.save('outpath')
#    f=open('outpath', 'r')
#    df2=load(f)    



