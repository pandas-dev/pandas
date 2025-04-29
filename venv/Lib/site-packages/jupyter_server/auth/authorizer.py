"""An Authorizer for use in the Jupyter server.

The default authorizer (AllowAllAuthorizer)
allows all authenticated requests

.. versionadded:: 2.0
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from typing import TYPE_CHECKING

from traitlets import Instance
from traitlets.config import LoggingConfigurable

from .identity import IdentityProvider, User

if TYPE_CHECKING:
    from collections.abc import Awaitable

    from jupyter_server.base.handlers import JupyterHandler


class Authorizer(LoggingConfigurable):
    """Base class for authorizing access to resources
    in the Jupyter Server.

    All authorizers used in Jupyter Server
    should inherit from this base class and, at the very minimum,
    implement an ``is_authorized`` method with the
    same signature as in this base class.

    The ``is_authorized`` method is called by the ``@authorized`` decorator
    in JupyterHandler. If it returns True, the incoming request
    to the server is accepted; if it returns False, the server
    returns a 403 (Forbidden) error code.

    The authorization check will only be applied to requests
    that have already been authenticated.

    .. versionadded:: 2.0
    """

    identity_provider = Instance(IdentityProvider)

    def is_authorized(
        self, handler: JupyterHandler, user: User, action: str, resource: str
    ) -> Awaitable[bool] | bool:
        """A method to determine if ``user`` is authorized to perform ``action``
        (read, write, or execute) on the ``resource`` type.

        Parameters
        ----------
        user : jupyter_server.auth.User
            An object representing the authenticated user,
            as returned by :meth:`jupyter_server.auth.IdentityProvider.get_user`.

        action : str
            the category of action for the current request: read, write, or execute.

        resource : str
            the type of resource (i.e. contents, kernels, files, etc.) the user is requesting.

        Returns
        -------
        bool
            True if user authorized to make request; False, otherwise
        """
        raise NotImplementedError


class AllowAllAuthorizer(Authorizer):
    """A no-op implementation of the Authorizer

    This authorizer allows all authenticated requests.

    .. versionadded:: 2.0
    """

    def is_authorized(
        self, handler: JupyterHandler, user: User, action: str, resource: str
    ) -> bool:
        """This method always returns True.

        All authenticated users are allowed to do anything in the Jupyter Server.
        """
        return True
