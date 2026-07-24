from .reference import Reference

def SeriesFactory(
    values: Reference | str,
    xvalues: Reference | str | None = None,
    zvalues: Reference | str | None = None,
    title: object = None,
    title_from_data: bool = False,
): ...
