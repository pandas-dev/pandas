import collections
import inspect
from functools import partial

import numpy as np

from pandas._libs import groupby as libgroupby
from pandas.core.dtypes.missing import isnull


class CythonDispatcher(object):

    def __init__(self, groupby):
        self.groupby = groupby
        self.func_nm = None
        self.obj = None

    @property
    def func_metadata(self):
        """
        Stores the metadata required to dispatch each function.

        The format of the dict is as follows:

        attr_name : {
            'application': {'aggregate', 'transform'}
            'cython_nm': ...  # Name of the Cython function to call
            'extra_kwargs': {...}  # Extra kwargs to pass to Cython
            'type_blacklist': [...]  # Dtypes for which func should raise
            'result_type': ...  # dtype of result from Cython
            'conversion_in': ...  # dtype or callable for conversion pre-Cython
            'needs_values': ...  # Whether the obj values should pass to Cython
            'needs_mask': ...   # Whether a mask of NA values should be passed
            'conversion_out': ...  # dtype or callable for conv post-Cython
        }
        """
        return {
            'any': {
                'application': 'aggregate',
                'cython_nm' : 'group_any_all',
                'extra_kwargs': {'val_test': 'any'},
                'type_blacklist': [],
                'result_type': np.uint8,
                'conversion_in': self._any_all_convertor,
                'needs_values': True,
                'needs_mask': True,
                'conversion_out': np.bool
            },
            'all': {
                'application': 'aggregate',
                'cython_nm' : 'group_any_all',
                'extra_kwargs': {'val_test': 'all'},
                'type_blacklist': [],
                'result_type': np.uint8,
                'conversion_in': self._any_all_convertor,
                'needs_values': True,
                'needs_mask': True,
                'conversion_out': np.bool
            }
        }

    @property
    def application_type(self):
        return self.func_metadata[self.func_nm]['application']

    def _any_all_convertor(self, vals):
        """
        Converts objects to appropriate type for any/all calculations.
        """
        try:
            vals = vals.astype(np.bool)
        except ValueError:  # for objects
            vals = np.array([bool(x) for x in vals])

        return vals.view(np.uint8)

    def _validate_types(self):
        """
        Validate that the types of the `grp_by` object.

        Raises
        ------
        ``TypeError`` if the `grp_by` dtypes are not valid for `func_nm`.
        """
        if self.obj.values.dtype in self.func_metadata[
                self.func_nm]['type_blacklist']:
            raise TypeError("'{}' cannot be applied to a dtype of {}".format(
                self.func_nm, self.obj.values.dtype))

    def _get_result(self, **kwargs):
        """
        Fetch the result from the Cython layer.

        Parameters
        ----------
        kwargs
            Extra arguments to bind to the `func_nm` Cython signature.

        Resolve function name in case of templating use.
        """
        # Since this func is called in a loop, the below might be better
        # served outside of the loop and passed in?
        labels, _, ngroups = self.groupby.grouper.group_info

        if self.application_type == 'aggregate':
            res_sz = ngroups
        elif self.application_type == 'transform':
            res_sz = len(labels)

        res_type = self.func_metadata[self.func_nm].get('result_type',
                                                        self.obj.values.dtype)

        result = np.zeros(res_sz, dtype=res_type)
        base_func = getattr(libgroupby,
                            self.func_metadata[self.func_nm]['cython_nm'])
        func = partial(base_func, result, labels)

        if self.func_metadata[self.func_nm].get('needs_values'):
            conv_in = self.func_metadata[self.func_nm].get('conversion_in')
            vals = self.obj.values
            # Below conditional needs refactoring but essentially want
            # to differentiate callables from dtypes
            if callable(conv_in) and not inspect.isclass(conv_in):
                vals = conv_in(self.obj.values)
            elif conv_in:  # is a type to convert to
                vals = self.obj.values.astype(conv, copy=False)
            func = partial(func, vals)

        if self.func_metadata[self.func_nm].get('needs_values'):
            mask = isnull(self.obj.values).view(np.uint8)
            func = partial(func, mask)

        # Not backwards compatible (py>=3.5 only)
        cy_kwargs = {**kwargs, **self.func_metadata[self.func_nm].get(
            'extra_kwargs', {})}
        func(**cy_kwargs)

        conv_out = self.func_metadata[self.func_nm].get('conversion_out')
        # Just like before, this needs refactoring
        if callable(conv_out) and not inspect.isclass(conv_out):
            result = conv_out(result)
        elif conv_out:
            result = result.astype(conv_out, copy=False)

        return result

    def _wrap_output(self, output):
        """
        Bind and apply the appropriate wrap func from `self.groupby`.
        """
        if self.application_type == 'aggregate':
            return getattr(self.groupby, '_wrap_aggregated_output')(output)
        elif self.application_type == 'transform':
            return getattr(self.groupby, '_wrap_transformed_output')(output)

        raise ValueError("Unknown application type for {}".format(
            self.func_nm))

    def dispatch(self, func_nm, **kwargs):
        """
        Dispatch the `func_nm` appropriately to the Cython layer.

        Will resolve any type and conversion dependencies, as well as apply
        any post-Cython conversions required for the given `func_nm`.

        Parameters
        ----------
        func_nm : str
            Conceptual name of the function to be applied.
        kwargs
            Extra arguments to bind to the `func_nm` Cython signature.

        Returns
        -------
        ndarray
            Result of Cython operation with appropriate conversions applied.
        """
        self.func_nm = func_nm

        output = collections.OrderedDict()
        for name, obj in self.groupby._iterate_slices():
            self.obj = obj
            self._validate_types()
            output[name] = self._get_result(**kwargs)

        return self._wrap_output(output)
