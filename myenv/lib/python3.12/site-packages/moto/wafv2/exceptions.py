from moto.core.exceptions import JsonRESTError


class WAFv2ClientError(JsonRESTError):
    code = 400


class WAFV2DuplicateItemException(WAFv2ClientError):
    def __init__(self) -> None:
        super().__init__(
            "WafV2DuplicateItem",
            "AWS WAF could not perform the operation because some resource in your request is a duplicate of an existing one.",
        )


class WAFNonexistentItemException(WAFv2ClientError):
    def __init__(self) -> None:
        super().__init__(
            "WAFNonexistentItemException",
            "AWS WAF couldn’t perform the operation because your resource doesn’t exist.",
        )


class WAFOptimisticLockException(WAFv2ClientError):
    def __init__(self) -> None:
        super().__init__(
            "WAFOptimisticLockException",
            "AWS WAF couldn’t save your changes because someone changed the resource after you started to edit it. Reapply your changes.",
        )
