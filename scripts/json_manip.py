"""

Tasks
-------

Search and transform jsonable structures, specifically to make it 'easy' to make tabular/csv output for other consumers.

Example
~~~~~~~~~~~~~

    *give me a list of all the fields called 'id' in this stupid, gnarly
    thing*

    >>> Q('id',gnarly_data)
    ['id1','id2','id3']


Observations:
---------------------

1) 'simple data structures' exist and are common.  They are tedious
   to search.

2)  The DOM is another nested / treeish structure, and jQuery selector is
    a good tool for that.

3a) R, Numpy, Excel and other analysis tools want 'tabular' data.  These
    analyses are valuable and worth doing.

3b) Dot/Graphviz, NetworkX, and some other analyses *like* treeish/dicty
    things, and those analyses are also worth doing!

3c) Some analyses are best done using 'one-off' and custom code in C, Python,
    or another 'real' programming language.

4)  Arbitrary transforms are tedious and error prone.  SQL is one solution,
    XSLT is another,

5)  the XPATH/XML/XSLT family is.... not universally loved :)  They are
    very complete, and the completeness can make simple cases... gross.

6)  For really complicated data structures, we can write one-off code.  Getting
    80% of the way is mostly okay.  There will always have to be programmers
    in the loop.

7)  Re-inventing SQL is probably a failure mode.  So is reinventing XPATH, XSLT
    and the like.  Be wary of mission creep!  Re-use when possible (e.g., can
    we put the thing into a DOM using

8)  If the interface is good, people can improve performance later.


Simplifying
---------------


1)  Assuming 'jsonable' structures

2)  keys are strings or stringlike.  Python allows any hashable to be a key.
    for now, we pretend that doesn't happen.

3)  assumes most dicts are 'well behaved'.  DAG, no cycles!

4)  assume that if people want really specialized transforms, they can do it
    themselves.

"""
from __future__ import print_function

from collections import namedtuple
import csv
import itertools
from itertools import product
from operator import attrgetter as aget, itemgetter as iget
import operator
import sys
from pandas.compat import map, u, callable, Counter
import pandas.compat as compat


##  note 'url' appears multiple places and not all extensions have same struct
ex1 = {
    'name': 'Gregg',
    'extensions': [
        {'id':'hello',
         'url':'url1'},
        {'id':'gbye',
         'url':'url2',
         'more': dict(url='url3')},
    ]
}

## much longer example
ex2 = {u('metadata'): {u('accessibilities'): [{u('name'): u('accessibility.tabfocus'),
    u('value'): 7},
   {u('name'): u('accessibility.mouse_focuses_formcontrol'), u('value'): False},
   {u('name'): u('accessibility.browsewithcaret'), u('value'): False},
   {u('name'): u('accessibility.win32.force_disabled'), u('value'): False},
   {u('name'): u('accessibility.typeaheadfind.startlinksonly'), u('value'): False},
   {u('name'): u('accessibility.usebrailledisplay'), u('value'): u('')},
   {u('name'): u('accessibility.typeaheadfind.timeout'), u('value'): 5000},
   {u('name'): u('accessibility.typeaheadfind.enabletimeout'), u('value'): True},
   {u('name'): u('accessibility.tabfocus_applies_to_xul'), u('value'): False},
   {u('name'): u('accessibility.typeaheadfind.flashBar'), u('value'): 1},
   {u('name'): u('accessibility.typeaheadfind.autostart'), u('value'): True},
   {u('name'): u('accessibility.blockautorefresh'), u('value'): False},
   {u('name'): u('accessibility.browsewithcaret_shortcut.enabled'),
    u('value'): True},
   {u('name'): u('accessibility.typeaheadfind.enablesound'), u('value'): True},
   {u('name'): u('accessibility.typeaheadfind.prefillwithselection'),
    u('value'): True},
   {u('name'): u('accessibility.typeaheadfind.soundURL'), u('value'): u('beep')},
   {u('name'): u('accessibility.typeaheadfind'), u('value'): False},
   {u('name'): u('accessibility.typeaheadfind.casesensitive'), u('value'): 0},
   {u('name'): u('accessibility.warn_on_browsewithcaret'), u('value'): True},
   {u('name'): u('accessibility.usetexttospeech'), u('value'): u('')},
   {u('name'): u('accessibility.accesskeycausesactivation'), u('value'): True},
   {u('name'): u('accessibility.typeaheadfind.linksonly'), u('value'): False},
   {u('name'): u('isInstantiated'), u('value'): True}],
  u('extensions'): [{u('id'): u('216ee7f7f4a5b8175374cd62150664efe2433a31'),
    u('isEnabled'): True},
   {u('id'): u('1aa53d3b720800c43c4ced5740a6e82bb0b3813e'), u('isEnabled'): False},
   {u('id'): u('01ecfac5a7bd8c9e27b7c5499e71c2d285084b37'), u('isEnabled'): True},
   {u('id'): u('1c01f5b22371b70b312ace94785f7b0b87c3dfb2'), u('isEnabled'): True},
   {u('id'): u('fb723781a2385055f7d024788b75e959ad8ea8c3'), u('isEnabled'): True}],
  u('fxVersion'): u('9.0'),
  u('location'): u('zh-CN'),
  u('operatingSystem'): u('WINNT Windows NT 5.1'),
  u('surveyAnswers'): u(''),
  u('task_guid'): u('d69fbd15-2517-45b5-8a17-bb7354122a75'),
  u('tpVersion'): u('1.2'),
  u('updateChannel'): u('beta')},
 u('survey_data'): {
  u('extensions'): [{u('appDisabled'): False,
    u('id'): u('testpilot?labs.mozilla.com'),
    u('isCompatible'): True,
    u('isEnabled'): True,
    u('isPlatformCompatible'): True,
    u('name'): u('Test Pilot')},
   {u('appDisabled'): True,
    u('id'): u('dict?www.youdao.com'),
    u('isCompatible'): False,
    u('isEnabled'): False,
    u('isPlatformCompatible'): True,
    u('name'): u('Youdao Word Capturer')},
   {u('appDisabled'): False,
    u('id'): u('jqs?sun.com'),
    u('isCompatible'): True,
    u('isEnabled'): True,
    u('isPlatformCompatible'): True,
    u('name'): u('Java Quick Starter')},
   {u('appDisabled'): False,
    u('id'): u('?20a82645-c095-46ed-80e3-08825760534b?'),
    u('isCompatible'): True,
    u('isEnabled'): True,
    u('isPlatformCompatible'): True,
    u('name'): u('Microsoft .NET Framework Assistant')},
   {u('appDisabled'): False,
    u('id'): u('?a0d7ccb3-214d-498b-b4aa-0e8fda9a7bf7?'),
    u('isCompatible'): True,
    u('isEnabled'): True,
    u('isPlatformCompatible'): True,
    u('name'): u('WOT')}],
  u('version_number'): 1}}

# class SurveyResult(object):

#     def __init__(self, record):
#         self.record = record
#         self.metadata, self.survey_data = self._flatten_results()

#     def _flatten_results(self):
#         survey_data = self.record['survey_data']
#         extensions = DataFrame(survey_data['extensions'])

def denorm(queries,iterable_of_things,default=None):
    """
    'repeat', or 'stutter' to 'tableize' for downstream.
    (I have no idea what a good word for this is!)

    Think ``kronecker`` products, or:

    ``SELECT single,multiple FROM table;``

    single   multiple
    -------  ---------
    id1      val1
    id1      val2


    Args:

        queries:  iterable of ``Q`` queries.
        iterable_of_things:  to be queried.

    Returns:

        list of 'stuttered' output, where if a query returns
        a 'single', it gets repeated appropriately.


    """

    def _denorm(queries,thing):
        fields = []
        results = []
        for q in queries:
            #print(q)
            r = Ql(q,thing)
            #print("-- result: ", r)
            if not r:
                r = [default]
            if isinstance(r[0], type({})):
                fields.append(sorted(r[0].keys()))  # dicty answers
            else:
                fields.append([q])  # stringy answer

            results.append(r)

        #print(results)
        #print(fields)
        flist =  list(flatten(*map(iter,fields)))

        prod = itertools.product(*results)
        for p in prod:
            U = dict()
            for (ii,thing) in enumerate(p):
                #print(ii,thing)
                if isinstance(thing, type({})):
                    U.update(thing)
                else:
                    U[fields[ii][0]] = thing

            yield U

    return list(flatten(*[_denorm(queries,thing) for thing in iterable_of_things]))


def default_iget(fields,default=None,):
    """ itemgetter with 'default' handling, that *always* returns lists

    API CHANGES from ``operator.itemgetter``

    Note: Sorry to break the iget api... (fields vs *fields)
    Note: *always* returns a list... unlike itemgetter,
        which can return tuples or 'singles'
    """
    myiget = operator.itemgetter(*fields)
    L = len(fields)
    def f(thing):
        try:
            ans = list(myiget(thing))
            if L < 2:
                ans = [ans,]
            return ans
        except KeyError:
            # slower!
            return [thing.get(x,default) for x in fields]

    f.__doc__ = "itemgetter with default %r for fields %r" %(default,fields)
    f.__name__ = "default_itemgetter"
    return f


def flatten(*stack):
    """
    helper function for flattening iterables of generators in a
    sensible way.
    """
    stack = list(stack)
    while stack:
        try: x = next(stack[0])
        except StopIteration:
            stack.pop(0)
            continue
        if hasattr(x,'next') and callable(getattr(x,'next')):
            stack.insert(0, x)

        #if isinstance(x, (GeneratorType,listerator)):
        else: yield x


def _Q(filter_, thing):
    """ underlying machinery for Q function recursion """
    T = type(thing)
    if isinstance({}, T):
        for k,v in compat.iteritems(thing):
            #print(k,v)
            if filter_ == k:
                if isinstance(v, type([])):
                    yield iter(v)
                else:
                    yield v

            if type(v)  in (type({}),type([])):
                yield Q(filter_,v)

    elif isinstance([], T):
        for k in thing:
            #print(k)
            yield Q(filter_,k)

    else:
        # no recursion.
        pass

def Q(filter_,thing):
    """
    type(filter):
    - list:  a flattened list of all searches (one list)
    - dict:  dict with vals each of which is that search

    Notes:

    [1] 'parent thing', with space, will do a descendent
    [2] this will come back 'flattened' jQuery style
    [3] returns a generator.  Use ``Ql`` if you want a list.

    """
    if isinstance(filter_, type([])):
        return flatten(*[_Q(x,thing) for x in filter_])
    elif isinstance(filter_, type({})):
        d = dict.fromkeys(list(filter_.keys()))
        #print(d)
        for k in d:
            #print(flatten(Q(k,thing)))
            d[k] = Q(k,thing)

        return d

    else:
        if " " in filter_:   # i.e. "antecendent post"
            parts = filter_.strip().split()
            r = None
            for p in parts:
                r = Ql(p,thing)
                thing = r

            return r

        else:  # simple.
            return flatten(_Q(filter_,thing))

def Ql(filter_,thing):
    """ same as Q, but returns a list, not a generator """
    res = Q(filter_,thing)

    if isinstance(filter_, type({})):
        for k in res:
            res[k] = list(res[k])
        return res

    else:
        return list(res)



def countit(fields,iter_of_iter,default=None):
    """
    note: robust to fields not being in i_of_i, using ``default``
    """
    C = Counter()  # needs hashables
    T = namedtuple("Thing",fields)
    get = default_iget(*fields,default=default)
    return Counter(
        (T(*get(thing)) for thing in iter_of_iter)
    )


## right now this works for one row...
def printout(queries,things,default=None, f=sys.stdout, **kwargs):
    """ will print header and objects

    **kwargs go to csv.DictWriter

    help(csv.DictWriter) for more.
    """

    results = denorm(queries,things,default=None)
    fields = set(itertools.chain(*(x.keys() for x in results)))

    W = csv.DictWriter(f=f,fieldnames=fields,**kwargs)
    #print("---prod---")
    #print(list(prod))
    W.writeheader()
    for r in results:
        W.writerow(r)


def test_run():
    print("\n>>> print(list(Q('url',ex1)))")
    print(list(Q('url',ex1)))
    assert  list(Q('url',ex1)) == ['url1','url2','url3']
    assert Ql('url',ex1) == ['url1','url2','url3']

    print("\n>>>  print(list(Q(['name','id'],ex1)))")
    print(list(Q(['name','id'],ex1)))
    assert Ql(['name','id'],ex1) == ['Gregg','hello','gbye']


    print("\n>>> print(Ql('more url',ex1))")
    print(Ql('more url',ex1))


    print("\n>>> list(Q('extensions',ex1))")
    print(list(Q('extensions',ex1)))

    print("\n>>> print(Ql('extensions',ex1))")
    print(Ql('extensions',ex1))

    print("\n>>> printout(['name','extensions'],[ex1,], extrasaction='ignore')")
    printout(['name','extensions'],[ex1,], extrasaction='ignore')

    print("\n\n")

    from pprint import pprint as pp

    print("-- note that the extension fields are also flattened!  (and N/A) -- ")
    pp(denorm(['location','fxVersion','notthere','survey_data extensions'],[ex2,], default="N/A")[:2])


if __name__ == "__main__":
    pass
