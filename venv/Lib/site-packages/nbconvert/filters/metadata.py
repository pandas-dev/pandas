"""filters for metadata"""


def get_metadata(output, key, mimetype=None):
    """Resolve an output metadata key

    If mimetype given, resolve at mimetype level first,
    then fallback to top-level.
    Otherwise, just resolve at top-level.
    Returns None if no data found.
    """
    md = output.get("metadata") or {}
    if mimetype and mimetype in md:
        value = md[mimetype].get(key)
        if value is not None:
            return value
    return md.get(key)
