class MissingOperands(ValueError):
    def __init__(self, expr, missing_ops):
        self.expr = expr
        self.missing_ops = missing_ops

        message = f'Lazy expression "{expr}" with missing operands: {missing_ops}'
        super().__init__(message)
