"""Simple registry of singletons."""


class Globals:
    """Accessors to doit singletons.

    :cvar dep_manager: (doit.dependency.Dependency) The doit dependency manager,
                       holding all persistent task data.
    """
    dep_manager = None
