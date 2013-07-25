import abc

from pandas.core import common as com
from pandas.computation.align import _align, _reconstruct_object


class AbstractEngine(object):
    """AbstractEngine object serving as a base class for all engines."""
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

        res = self._evaluate(self.expr.env)
        return _reconstruct_object(self.result_type, res, self.aligned_axes,
                                   self.expr.terms.return_type)

    @property
    def _is_aligned(self):
        return self.aligned_axes is not None and self.result_type is not None

    @abc.abstractmethod
    def _evaluate(self, env):
        """Return an evaluated expression.

        Parameters
        ----------
        env : Scope
            The local and global environment in which to evaluate an
            expression.

        Notes
        -----
        This method must be implemented by any class the subclasses this class.
        """
        pass


class NumExprEngine(AbstractEngine):
    """NumExpr engine class"""
    has_neg_frac = True

    def __init__(self, expr):
        super(NumExprEngine, self).__init__(expr)

    def _evaluate(self, env):
        import numexpr as ne

        try:
            s = self.convert()
            return ne.evaluate(s, local_dict=env.locals,
                               global_dict=env.globals,
                               truediv=self.expr.truediv)
        except KeyError as e:
            raise NameError('{0!r} is not defined'.format(e.message))


class PythonEngine(AbstractEngine):
    """Evaluate an expression in Python space.

    Mostly for testing purposes.
    """
    has_neg_frac = False

    def __init__(self, expr):
        super(PythonEngine, self).__init__(expr)

    def evaluate(self):
        return self.expr(self.expr.env)

    def _evaluate(self, env):
        pass


_engines = {'numexpr': NumExprEngine, 'python': PythonEngine}
