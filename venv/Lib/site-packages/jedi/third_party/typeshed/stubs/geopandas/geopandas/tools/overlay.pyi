from ..geodataframe import GeoDataFrame

def overlay(
    df1: GeoDataFrame, df2: GeoDataFrame, how: str = "intersection", keep_geom_type: bool | None = None, make_valid: bool = True
) -> GeoDataFrame: ...
