"""Engine classes for :func:`~pandas.eval`
"""

import abc

from pandas import compat
from pandas.core import common as com
from pandas.computation.align import _align, _reconstruct_object
from pandas.computation.ops import UndefinedVariableError


class AbstractEngine(object):

    """Object serving as a base class for all engines."""

    __metaclass__ = abc.ABCMeta

    has_neg_frac = False

    def __init__(self, expr):
        self.expr = expr
        self.aligned_axes = None
        self.result_type = None

    def convert(self):
        """Convert an expression for evaluation.

        Defaults to return the expression as a string.
        """
        return com.pprint_thing(self.expr)

    def pre_evaluate(self):
        self.expr.check_name_clashes()

    def evaluate(self):
        """Run the engine on the expression

        This method performs alignment which is necessary no matter what engine
        is being used, thus its implementation is in the base class.

        Returns
        -------
        obj : object
            The result of the passed expression.
        """
        if not self._is_aligned:
            self.result_type, self.aligned_axes = _align(self.expr.terms)

        # make sure no names in resolvers and locals/globals clash
        self.pre_evaluate()
        res = self._evaluate()
        return _reconstruct_object(self.result_type, res, self.aligned_axes,
                                   self.expr.terms.return_type)

    @property
    def _is_aligned(self):
        return self.aligned_axes is not None and self.result_type is not None

    @abc.abstractmethod
    def _evaluate(self):
        """Return an evaluated expression.

        Parameters
        ----------
        env : Scope
            The local and global environment in which to evaluate an
            expression.

        Notes
        -----
        Must be implemented by subclasses.
        """
        pass


class NumExprEngine(AbstractEngine):

    """NumExpr engine class"""
    has_neg_frac = True

    def __init__(self, expr):
        super(NumExprEngine, self).__init__(expr)

    def convert(self):
        return str(super(NumExprEngine, self).convert())

    def _evaluate(self):
        import numexpr as ne

        # add the resolvers to locals
        self.expr.add_resolvers_to_locals()

        # convert the expression to a valid numexpr expression
        s = self.convert()

        try:
            return ne.evaluate(s, local_dict=self.expr.env.locals,
                               global_dict=self.expr.env.globals,
                               truediv=self.expr.truediv)
        except KeyError as e:
            # python 3 compat kludge
            try:
                msg = e.message
            except AttributeError:
                msg = compat.text_type(e)
            raise UndefinedVariableError(msg)


class PythonEngine(AbstractEngine):

    """Evaluate an expression in Python space.

    Mostly for testing purposes.
    """
    has_neg_frac = False

    def __init__(self, expr):
        super(PythonEngine, self).__init__(expr)

    def evaluate(self):
        self.pre_evaluate()
        return self.expr()

    def _evaluate(self):
        pass


_engines = {'numexpr': NumExprEngine, 'python': PythonEngine}
