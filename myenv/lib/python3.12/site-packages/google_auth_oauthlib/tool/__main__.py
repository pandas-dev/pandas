# Copyright (C) 2017 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Command-line tool for obtaining authorization and credentials from a user.

This tool uses the OAuth 2.0 Authorization Code grant as described in
`section 1.3.1 of RFC6749`_ and implemeted by
:class:`google_auth_oauthlib.flow.Flow`.

This tool is intended for assist developers in obtaining credentials
for testing applications where it may not be possible or easy to run a
complete OAuth 2.0 authorization flow, especially in the case of code
samples or embedded devices without input / display capabilities.

This is not intended for production use where a combination of
companion and on-device applications should complete the OAuth 2.0
authorization flow to get authorization from the users.

.. _section 1.3.1 of RFC6749: https://tools.ietf.org/html/rfc6749#section-1.3.1
"""

import json
import os
import os.path

import click

import google_auth_oauthlib.flow


APP_NAME = "google-oauthlib-tool"
DEFAULT_CREDENTIALS_FILENAME = "credentials.json"


@click.command()
@click.option(
    "--client-secrets",
    metavar="<client_secret_json_file>",
    required=True,
    help="Path to OAuth2 client secret JSON file.",
)
@click.option(
    "--scope",
    multiple=True,
    metavar="<oauth2 scope>",
    required=True,
    help="API scopes to authorize access for.",
)
@click.option(
    "--save",
    is_flag=True,
    metavar="<save_mode>",
    show_default=True,
    default=False,
    help="Save the credentials to file.",
)
@click.option(
    "--credentials",
    metavar="<oauth2_credentials>",
    show_default=True,
    default=os.path.join(click.get_app_dir(APP_NAME), DEFAULT_CREDENTIALS_FILENAME),
    help="Path to store OAuth2 credentials.",
)
def main(client_secrets, scope, save, credentials):
    """Command-line tool for obtaining authorization and credentials from a user.

    This tool uses the OAuth 2.0 Authorization Code grant as described
    in section 1.3.1 of RFC6749:
    https://tools.ietf.org/html/rfc6749#section-1.3.1

    This tool is intended for assist developers in obtaining credentials
    for testing applications or samples.

    This is not intended for production use where a combination of
    companion and on-device applications should complete the OAuth 2.0
    authorization flow to get authorization from the users.

    """

    flow = google_auth_oauthlib.flow.InstalledAppFlow.from_client_secrets_file(
        client_secrets, scopes=scope
    )

    creds = flow.run_local_server()

    creds_data = {
        "token": creds.token,
        "refresh_token": creds.refresh_token,
        "token_uri": creds.token_uri,
        "client_id": creds.client_id,
        "client_secret": creds.client_secret,
        "scopes": creds.scopes,
    }

    if save:
        del creds_data["token"]

        config_path = os.path.dirname(credentials)
        if config_path and not os.path.isdir(config_path):
            os.makedirs(config_path)

        with open(credentials, "w") as outfile:
            json.dump(creds_data, outfile)

        click.echo("credentials saved: %s" % credentials)

    else:
        click.echo(json.dumps(creds_data))


if __name__ == "__main__":
    # pylint doesn't realize that click has changed the function signature.
    main()  # pylint: disable=no-value-for-parameter
