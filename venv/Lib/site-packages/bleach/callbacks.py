"""A set of basic callbacks for bleach.linkify."""


def nofollow(attrs, new=False):
    href_key = (None, "href")

    if href_key not in attrs:
        return attrs

    if attrs[href_key].startswith("mailto:"):
        return attrs

    rel_key = (None, "rel")
    rel_values = [val for val in attrs.get(rel_key, "").split(" ") if val]
    if "nofollow" not in [rel_val.lower() for rel_val in rel_values]:
        rel_values.append("nofollow")
    attrs[rel_key] = " ".join(rel_values)

    return attrs


def target_blank(attrs, new=False):
    href_key = (None, "href")

    if href_key not in attrs:
        return attrs

    if attrs[href_key].startswith("mailto:"):
        return attrs

    attrs[(None, "target")] = "_blank"
    return attrs
