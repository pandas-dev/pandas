# flake8: noqa

from pandas.computation.eval import eval


# deprecation, xref #13790
def Expr(*args, **kwargs):
    import warnings

    warnings.warn("pd.Expr is deprecated as it is not "
                  "applicable to user code",
                  FutureWarning, stacklevel=2)
    from pandas.computation.expr import Expr
    return Expr(*args, **kwargs)
