import abc

from pandas.computation.align import _align, _reconstruct_object


class AbstractEngine(object):
    """"""
    __metaclass__ = abc.ABCMeta

    has_neg_frac = False

    def __init__(self, expr):
        self.expr = expr
        self.aligned_axes = None
        self.result_type = None

    @abc.abstractmethod
    def convert(self):
        """Convert an expression for evaluation."""
        pass

    def evaluate(self):
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
        """Return an evaluated expression."""
        pass


class NumExprEngine(AbstractEngine):
    """NumExpr engine class"""
    has_neg_frac = True

    def __init__(self, expr):
        super(NumExprEngine, self).__init__(expr)

    def convert(self):
        """Return a string"""
        return '%s' % self.expr

    def _evaluate(self, env):
        import numexpr as ne

        try:
            return ne.evaluate(self.convert(), local_dict=env.locals,
                               global_dict=env.globals,
                               truediv=self.expr.truediv)
        except KeyError as e:
            raise NameError('{0!r} is not defined'.format(e.message))


class PythonEngine(AbstractEngine):
    """Use NumPy even if numexpr is installed"""
    has_neg_frac = False

    def __init__(self, expr):
        super(PythonEngine, self).__init__(expr)

    def convert(self):
        pass

    def evaluate(self):
        return self.expr(self.expr.env)

    def _evaluate(self, env):
        pass


_engines = {'numexpr': NumExprEngine, 'python': PythonEngine }
