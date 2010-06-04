import numpy as np

from pandas import DataFrame, DataMatrix

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

r = robjects.r

def load_data(name, package=None):
    if package:
        importr(package)

    r.data(name)
    return _convert_robj(r[name])

def _convert_robj(robj):
    if isinstance(robj, robjects.DataFrame):
        return _from_DataFrame(robj)
    elif isinstance(robj, robjects.Matrix):
        return _from_Matrix(robj)

def _from_DataFrame(rdf):
    columns = list(rdf.colnames)

    data = {}
    for i, col in enumerate(columns):
        vec = rdf.rx2(i + 1)
        data[col] = list(vec)

    return DataFrame(data)

def _from_Matrix(mat):
    columns = list(mat.colnames)

    return DataMatrix(np.array(mat), columns=columns)

