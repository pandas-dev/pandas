"""
Expressions
-----------

Offer fast expression evaluation thru numexpr

"""
import numpy as np

try:
    import numexpr as ne
    _NUMEXPR_INSTALLED = True
except ImportError:  # pragma: no cover
    _NUMEXPR_INSTALLED = False

_USE_NUMEXPR = _NUMEXPR_INSTALLED
_evaluate    = None
_where       = None

# the set of dtypes that we will allow pass to numexpr
_ALLOWED_DTYPES = dict(evaluate = set(['int64','int32','float64','float32','bool']),
                       where    = set(['int64','float64','bool']))

# the minimum prod shape that we will use numexpr
_MIN_ELEMENTS   = 10000

def set_use_numexpr(v = True):
    # set/unset to use numexpr
    global _USE_NUMEXPR
    if _NUMEXPR_INSTALLED:
        _USE_NUMEXPR = v

    # choose what we are going to do
    global _evaluate, _where
    if not _USE_NUMEXPR:
        _evaluate = _evaluate_standard
        _where    = _where_standard
    else:
        _evaluate = _evaluate_numexpr
        _where    = _where_numexpr

def set_numexpr_threads(n = None):
    # if we are using numexpr, set the threads to n
    # otherwise reset
    try:
        if _NUMEXPR_INSTALLED and _USE_NUMEXPR:
            if n is None:
                n = ne.detect_number_of_cores()
            ne.set_num_threads(n)
    except:
        pass


def _evaluate_standard(op, op_str, a, b, raise_on_error=True, **eval_kwargs):
    """ standard evaluation """
    return op(a,b)

def _can_use_numexpr(op, op_str, a, b, dtype_check):
    """ return a boolean if we WILL be using numexpr """
    if op_str is not None:
        
        # required min elements (otherwise we are adding overhead)
        if np.prod(a.shape) > _MIN_ELEMENTS:

            # check for dtype compatiblity
            dtypes = set()
            for o in [ a, b ]:
                if hasattr(o,'get_dtype_counts'):
                    s = o.get_dtype_counts()
                    if len(s) > 1:
                        return False
                    dtypes |= set(s.index)
                elif isinstance(o,np.ndarray):
                    dtypes |= set([o.dtype.name])

            # allowed are a superset
            if not len(dtypes) or _ALLOWED_DTYPES[dtype_check] >= dtypes:
                return True

    return False

def _evaluate_numexpr(op, op_str, a, b, raise_on_error = False, **eval_kwargs):
    result = None

    if _can_use_numexpr(op, op_str, a, b, 'evaluate'):
        try:
            a_value, b_value = a, b
            if hasattr(a_value,'values'):
                a_value = a_value.values
            if hasattr(b_value,'values'):
                b_value = b_value.values
            result = ne.evaluate('a_value %s b_value' % op_str, 
                                 local_dict={ 'a_value' : a_value, 
                                              'b_value' : b_value }, 
                                 casting='safe', **eval_kwargs)
        except (ValueError), detail:
            if 'unknown type object' in str(detail):
                pass
        except (Exception), detail:
            if raise_on_error:
                raise TypeError(str(detail))

    if result is None:
        result = _evaluate_standard(op,op_str,a,b,raise_on_error)

    return result

def _where_standard(cond, a, b, raise_on_error=True):           
    return np.where(cond, a, b)

def _where_numexpr(cond, a, b, raise_on_error = False):
    result = None

    if _can_use_numexpr(None, 'where', a, b, 'where'):

        try:
            cond_value, a_value, b_value = cond, a, b
            if hasattr(cond_value,'values'):
                cond_value = cond_value.values
            if hasattr(a_value,'values'):
                a_value = a_value.values
            if hasattr(b_value,'values'):
                b_value = b_value.values
            result = ne.evaluate('where(cond_value,a_value,b_value)',
                                 local_dict={ 'cond_value' : cond_value,
                                              'a_value' : a_value, 
                                              'b_value' : b_value }, 
                                 casting='safe')
        except (ValueError), detail:
            if 'unknown type object' in str(detail):
                pass
        except (Exception), detail:
            if raise_on_error:
                raise TypeError(str(detail))

    if result is None:
        result = _where_standard(cond,a,b,raise_on_error)

    return result


# turn myself on
set_use_numexpr(True)

def evaluate(op, op_str, a, b, raise_on_error=False, use_numexpr=True, **eval_kwargs):
    """ evaluate and return the expression of the op on a and b

        Parameters
        ----------

        op :    the actual operand
        op_str: the string version of the op
        a :     left operand
        b :     right operand
        raise_on_error : pass the error to the higher level if indicated (default is False),
                         otherwise evaluate the op with and return the results
        use_numexpr : whether to try to use numexpr (default True)
        """

    if use_numexpr:
        return _evaluate(op, op_str, a, b, raise_on_error=raise_on_error, **eval_kwargs)
    return _evaluate_standard(op, op_str, a, b, raise_on_error=raise_on_error)

def where(cond, a, b, raise_on_error=False, use_numexpr=True):
    """ evaluate the where condition cond on a and b

        Parameters
        ----------

        cond : a boolean array
        a :    return if cond is True
        b :    return if cond is False
        raise_on_error : pass the error to the higher level if indicated (default is False),
                         otherwise evaluate the op with and return the results
        use_numexpr : whether to try to use numexpr (default True)
        """

    if use_numexpr:
        return _where(cond, a, b, raise_on_error=raise_on_error)
    return _where_standard(cond, a, b, raise_on_error=raise_on_error)
