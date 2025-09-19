# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

def doc_subst(snippets):
    """ Substitute format strings in class or function docstring """
    def decorator(cls):
        # Running python with -OO will discard docstrings (__doc__ is None).
        if cls.__doc__ is not None:
            # Strip the snippets to avoid trailing new lines and whitespace
            stripped_snippets = {
                key: snippet.strip() for (key, snippet) in snippets.items()
            }
            cls.__doc__ = cls.__doc__.format(**stripped_snippets)
        return cls
    return decorator
