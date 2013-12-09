from __future__ import print_function
# see LICENSES directory for copyright and license
import os
import sys
import logging

import httplib2

import apiclient.discovery as gapi
import gflags
import oauth2client.file as auth_file
import oauth2client.client as oauth
import oauth2client.tools as tools
OOB_CALLBACK_URN = oauth.OOB_CALLBACK_URN


class AuthenticationConfigError(ValueError):
    pass

FLOWS = {}
FLAGS = gflags.FLAGS
DEFAULT_SECRETS = os.path.join(
    os.path.dirname(__file__), 'client_secrets.json')
DEFAULT_SCOPE = 'https://www.googleapis.com/auth/analytics.readonly'
DEFAULT_TOKEN_FILE = os.path.join(os.path.dirname(__file__), 'analytics.dat')
MISSING_CLIENT_MSG = """
WARNING: Please configure OAuth 2.0

You need to populate the client_secrets.json file found at:

   %s

with information from the APIs Console <https://code.google.com/apis/console>.

"""
DOC_URL = ('https://developers.google.com/api-client-library/python/guide/'
           'aaa_client_secrets')

gflags.DEFINE_enum('logging_level', 'ERROR',
                   ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                   'Set the level of logging detail.')

# Name of file that will store the access and refresh tokens to access
# the API without having to login each time. Make sure this file is in
# a secure place.


def process_flags(flags=[]):
    """Uses the command-line flags to set the logging level.

    Args:
    argv: List of command line arguments passed to the python script.
    """

    # Let the gflags module process the command-line arguments.
    try:
        FLAGS(flags)
    except gflags.FlagsError as e:
        print('%s\nUsage: %s ARGS\n%s' % (e, str(flags), FLAGS))
        sys.exit(1)

    # Set the logging according to the command-line flag.
    logging.getLogger().setLevel(getattr(logging, FLAGS.logging_level))


def get_flow(secret, scope, redirect):
    """
    Retrieve an authentication flow object based on the given
    configuration in the secret file name, the authentication scope,
    and a redirect URN
    """
    key = (secret, scope, redirect)
    flow = FLOWS.get(key, None)
    if flow is None:
        msg = MISSING_CLIENT_MSG % secret
        if not os.path.exists(secret):
            raise AuthenticationConfigError(msg)
        flow = oauth.flow_from_clientsecrets(secret, scope,
                                             redirect_uri=redirect,
                                             message=msg)
        FLOWS[key] = flow
    return flow


def make_token_store(fpath=None):
    """create token storage from give file name"""
    if fpath is None:
        fpath = DEFAULT_TOKEN_FILE
    return auth_file.Storage(fpath)


def authenticate(flow, storage=None):
    """
    Try to retrieve a valid set of credentials from the token store if possible
    Otherwise use the given authentication flow to obtain new credentials
    and return an authenticated http object

    Parameters
    ----------
    flow : authentication workflow
    storage: token storage, default None
    """
    http = httplib2.Http()

    # Prepare credentials, and authorize HTTP object with them.
    credentials = storage.get()
    if credentials is None or credentials.invalid:
        credentials = tools.run(flow, storage)

    http = credentials.authorize(http)
    return http


def init_service(http):
    """
    Use the given http object to build the analytics service object
    """
    return gapi.build('analytics', 'v3', http=http)


def reset_default_token_store():
    import os
    os.remove(DEFAULT_TOKEN_FILE)
