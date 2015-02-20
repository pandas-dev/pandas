#!/usr/bin/python3
'''
Useful tools for working with the python library pandas
Written in 2015 by Garrett Berg <garrett@cloudformdesign.com>

© Creative Commons 0
To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. THIS SOFTWARE IS DISTRIBUTED WITHOUT ANY WARRANTY.
<http://creativecommons.org/publicdomain/zero/1.0/>
'''


def _dataframe_dict(data, index=None, filler='', header=None):
    if isinstance(data, dict):
        try:
            if depth(data, isiter=True) < 2:
                return data
        except TypeError:
            return data
    if not isinstance(data, dict):
        header = resolve_header(header)
        if header is None:
            header = get_header(data[0])
        data = unpack(data, header)
    data = flatten(data)
    data = fill_keys(data, filler)
    return data


def dataframe_dict(data, index=None, filler='', header=None):
    '''General loader of dataframes from python objects. Can either be a
    dict of lists or a list of dicts.
    Header is detected automatically and will be multiindex if the dict
    is nested'''
    data = _dataframe_dict(data, index, filler, header)
    data = pd.DataFrame.from_dict(data)
    if index is not None:
        data.set_index(index, inplace=True)
        data.sort_index(inplace=True)
    return data


# Helper funcitons
def resolve_header(header):
    if header is None:
        return None
    if isinstance(header, dict):
        return get_header(header)
    else:
        return header


#!/usr/bin/python3
'''
Useful tools to inspect and work with dictionaries
Written in 2015 by Garrett Berg <garrett@cloudformdesign.com>

© Creative Commons 0
To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. THIS SOFTWARE IS DISTRIBUTED WITHOUT ANY WARRANTY.
<http://creativecommons.org/publicdomain/zero/1.0/>
'''

# builtins
import itertools
import collections

def _isiter(obj, exclude=(str, bytes, bytearray)):
    '''Returns True if object is an iterator.
    Returns False for str, bytes and bytearray objects
    by default'''
    return (False if isinstance(obj, exclude)
            else True if hasattr(obj, '__iter__')
            else False)

def consume(iterator, n=None):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(itertools.islice(iterator, n, n), None)


def throw(exception):
    '''Raises an exception. Can be used inside compressions'''
    raise exception

def depth(d, deep=0, isiter=False):
    '''Find the depth of a nested dictionary'''
    if not isinstance(d, dict) or not d:  # not a dict or an empty dict
        throw(TypeError) if isiter and \
            not _isiter(d) else None
        return deep
    return max(depth(v, deep + 1, isiter) for k, v in d.items())


# Dictionary methods
def get_header(item, extra_levels=None, filler=''):
    '''Returns the header of a nested dictionary
    The header is a list of tuples detailing the structure of the dictionary
    Very useful in pandas'''
    levels = extra_levels
    if levels is None:
        levels = depth(item)
    keys = []
    for key, value in item.items():
        if isinstance(value, dict):
            keys.extend((key,) + v for v in
                        get_header(value, levels - 1, filler))
        else:
            keys.append((key,))
    return keys


def getitem(dic, item):
    '''Dictionary item access with tuples'''
    for i in item:
        dic = dic[i]
    return dic


def setitem(dic, item, value):
    '''Dictionary item setting with tuples'''
    for i, k in enumerate(item):
        if i < len(item) - 1:
            if k not in dic:
                dic[k] = type(dic)()
            dic = dic[k]
        else:
            dic[k] = value
            break
    else:
        assert False
    return dic


def unpack(data, header=None):
    '''Unpacks a list of dictionaries into a dictionary of lists
    according to the header'''
    if header is None:
        header = get_header(data[0])
    out = type(data[0])()
    for key in header:
        setitem(out, key, [])
    for d in data:
        for h in header:
            getitem(out, h).append(getitem(d, h))
    return out


def flatten(data, start=()):
    '''Flattens a dictionary so that the keys are all tuples of keys'''
    flat = type(data)()
    for key, value in data.items():
        if isinstance(value, dict):
            flat.update(flatten(value, start=start + (key,)))
        else:
            flat[start + (key,)] = value
    return flat


def fill_keys(data, filler=None):
    '''Makes all dictionary keys tuples of the same length'''
    keys, values = zip(*data.items())
    # convert all keys to tuples
    keys = tuple(key if isinstance(key, tuple) else (key,) for key in keys)
    maxlen = max(map(len, keys))
    return type(data)((key + ((filler,) * (maxlen - len(key))), value) for
                       (key, value) in zip(keys, values))


def update(todict, fromdict, keys=None):
    '''Copy only keys from one dictionary to another

    keys=None is equivalent to todict.update(fromdict)
    '''
    todict.update(fromdict if keys is None else
                  {key: fromdict[key] for key in keys})


def remove(obj, keys, check=True):
    '''remove unwanted keys
    Arguments:
        obj -- object on which keys should be removed
        keys -- iterator of keys to remove
        check -- whether to check whether keys exist
        '''
    if check:
        consume(map(obj.pop, keys))
    else:
        consume(map(obj.pop, keys, itertools.repeat(None)))
