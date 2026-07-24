from ratelimit.decorators import RateLimitDecorator, sleep_and_retry
from ratelimit.exception import RateLimitException

limits = RateLimitDecorator
rate_limited = RateLimitDecorator

__all__ = ["RateLimitException", "limits", "rate_limited", "sleep_and_retry"]
