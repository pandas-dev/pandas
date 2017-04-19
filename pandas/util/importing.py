class _UnSubclassable(type):
    """
    Metaclass to raise an ImportError when subclassed
    """
    msg = ""

    def __init__(cls, name, bases, clsdict):
        if len(cls.mro()) > 2:
            raise ImportError(cls.msg)
        super(_UnSubclassable, cls).__init__(name, bases, clsdict)
