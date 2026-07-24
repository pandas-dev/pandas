from _typeshed import Incomplete

class XRayMiddleware:
    app: Incomplete
    in_lambda_ctx: bool
    def __init__(self, app, recorder) -> None: ...
