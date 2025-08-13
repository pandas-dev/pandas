def get_marker(item, name):
    """Get a marker from a pytest item.

    This is here in order to stay compatible with pytest < 3.6 and not produce
    any deprecation warnings with >= 3.6.
    """
    try:
        return item.get_closest_marker(name)
    except AttributeError:
        # pytest < 3.6
        return item.get_marker(name)
