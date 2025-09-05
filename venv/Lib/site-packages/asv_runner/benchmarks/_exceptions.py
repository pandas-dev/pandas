class NotRequired(ImportError):
    """
    Exception raised when a requirement is not met.

    This exception inherits from `ImportError`. It's typically used when a particular
    package, module or other dependency that is not essential for the overall function
    of the program is not found or doesn't meet specific requirements.

    #### Attributes
    **message** (`str`)
    : A string that provides a more detailed explanation of the error.

    #### Example
    This exception might be used in a scenario where an optional feature of a program
    relies on a specific package that is not installed:


    ```{code-block} python
    try:
        import optional_package
    except ImportError:
        raise NotRequired("optional_package is not installed.")
    ```
    """

    def __init__(self, message):
        """
        Initialize a new instance of `NotRequired`.

        #### Parameters
        **message** (`str`)
        : A string that provides a more detailed explanation of the error.
        """
        self.message = message
        super().__init__(self.message)
