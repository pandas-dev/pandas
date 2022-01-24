import pandas as pd
import numpy as np
from astropy import units as u

class QSeries(pd.Series):

    _metadata = ["unit"]

    def __init__(self, *args, unit=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.unit = unit

    @property
    def _constructor(self):
        return QSeries

    @property
    def _constructor_expanddim(self):
        return QDataFrame

    def to(self, unit):
        transformed_values = self.values * self.unit.to(unit)
        return QSeries(data=transformed_values, unit=unit)

class QDataFrame(pd.DataFrame):
    # normal properties
    _metadata = ["units"]

    def __init__(self, *args, units=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.units = units

    @property
    def units(self):
        return {column:self[column].unit 
                for column in self.columns 
                if hasattr(self[column], 'unit')}

    @units.setter
    def units(self, new_units):

        if new_units is None:
            return
        for column in self.columns:
            if column in new_units:
                self[column].unit = new_units[column]

    @property
    def _constructor(self):
        return QDataFrame

    @property
    def _constructor_sliced(self):
        return QSeries

    def __finalize__(self, other, method=None, **kwargs):
        """ 
        Code taken from: https://github.com/geopandas/geopandas/blob/master/geopandas/geodataframe.py
        """
        super().__finalize__(other, method=method, **kwargs)
        if hasattr(other, 'units'):
            self.units = other.units

        # merge operation: using metadata of the left object
        if method == "merge":
            for name in self._metadata:
                object.__setattr__(self, name, getattr(other.left, name, None))
        # concat operation: using metadata of the first object
        elif method == "concat":
            for name in self._metadata:
                object.__setattr__(self, name, getattr(other.objs[0], name, None))

        return self

    def _repr_html_(self) -> str:
        repr_html = super()._repr_html_()
        for col_name, col_unit in self.units.items():
            if col_unit is None:
                continue
            col_unit_latex = col_unit.to_string('latex')
            repr_html = repr_html.replace(f'<th>{col_name}</th>', f'<th>{col_name} [ {col_unit_latex} ]</th>')
        return repr_html 