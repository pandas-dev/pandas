"""
Neptune is a bit of an odd duck.
It shares almost everything with RDS: the endpoint URL, and the features. Only the parameters to these features can be different.

Because the endpoint URL is the same (rds.amazonaws.com), every request is intercepted by the RDS service.
RDS then has to determine whether any incoming call was meant for RDS, or for neptune.
"""

from .models import neptune_backends  # noqa: F401
