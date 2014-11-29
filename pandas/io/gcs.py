from StringIO import StringIO
import logging

import pandas as pd

_GOOGLE_API_CLIENT_INSTALLED = True
_GOOGLE_API_CLIENT_VALID_VERSION = True
_GOOGLE_FLAGS_INSTALLED = True
_GOOGLE_FLAGS_VALID_VERSION = True
_HTTPLIB2_INSTALLED = True
_SETUPTOOLS_INSTALLED = True

try:
    import pkg_resources
    _SETUPTOOLS_INSTALLED = True
except ImportError:
    _SETUPTOOLS_INSTALLED = False

if _SETUPTOOLS_INSTALLED:
    try:
        from apiclient.discovery import build
        from apiclient.errors import HttpError
        from apiclient.http import MediaFileUpload, MediaIoBaseUpload

        from oauth2client.client import AccessTokenRefreshError
        from oauth2client.client import OAuth2WebServerFlow
        from oauth2client.client import SignedJwtAssertionCredentials
        from oauth2client.client import flow_from_clientsecrets
        from oauth2client.file import Storage
        from oauth2client.tools import run

        _GOOGLE_API_CLIENT_INSTALLED=True
        _GOOGLE_API_CLIENT_VERSION = pkg_resources.get_distribution('google-api-python-client').version

        #if LooseVersion(_GOOGLE_API_CLIENT_VERSION) >= '1.2.0':
            #_GOOGLE_API_CLIENT_VALID_VERSION = True

    except ImportError:
        _GOOGLE_API_CLIENT_INSTALLED = False


    try:
        import gflags as flags
        _GOOGLE_FLAGS_INSTALLED = True

        _GOOGLE_FLAGS_VERSION = pkg_resources.get_distribution('python-gflags').version

        #if LooseVersion(_GOOGLE_FLAGS_VERSION) >= '2.0':
            #_GOOGLE_FLAGS_VALID_VERSION = True

    except ImportError:
        _GOOGLE_FLAGS_INSTALLED = False

    try:
        import httplib2
        _HTTPLIB2_INSTALLED = True
    except ImportError:
        _HTTPLIB2_INSTALLED = False

def _test_imports():
    _GOOGLE_API_CLIENT_INSTALLED
    _GOOGLE_API_CLIENT_VALID_VERSION
    _GOOGLE_FLAGS_INSTALLED
    _GOOGLE_FLAGS_VALID_VERSION
    _HTTPLIB2_INSTALLED
    _SETUPTOOLS_INSTALLED

    #if compat.PY3:
        #raise NotImplementedError("Google's libraries do not support Python 3 yet")

    if not _SETUPTOOLS_INSTALLED:
        raise ImportError('Could not import pkg_resources (setuptools).')

    if not _GOOGLE_API_CLIENT_INSTALLED:
        raise ImportError('Could not import Google API Client.')

    if not _GOOGLE_FLAGS_INSTALLED:
        raise ImportError('Could not import Google Command Line Flags Module.')

    if not _GOOGLE_API_CLIENT_VALID_VERSION:
        raise ImportError("pandas requires google-api-python-client >= 1.2.0 for Google "
                          "BigQuery support, current version " + _GOOGLE_API_CLIENT_VERSION)

    if not _GOOGLE_FLAGS_VALID_VERSION:
        raise ImportError("pandas requires python-gflags >= 2.0.0 for Google "
                          "BigQuery support, current version " + _GOOGLE_FLAGS_VERSION)

    if not _HTTPLIB2_INSTALLED:
        raise ImportError("pandas requires httplib2 for Google BigQuery support")

CHUNK_SIZE = 2 * 1024 * 1024

logger = logging.getLogger('pandas.io.gcs')
logger.setLevel(logging.DEBUG)

class GCSConnector(object):

    def __init__(self, project_id, service_account=None, private_key=None):
        self.project_id = project_id

        self.service_account = service_account
        self.private_key = private_key

        self.credentials = self.get_credentials(service_account, private_key)
        self.service = self.get_service(self.credentials)

    def get_service(self, credentials):

        http = httplib2.Http()
        http = credentials.authorize(http)

        service = build('storage', 'v1', http=http)

        return service

    def get_credentials(self, service_account, private_key):

        scope = "https://www.googleapis.com/auth/devstorage.read_write"

        if service_account and private_key:
            flow = SignedJwtAssertionCredentials(service_account, private_key, scope)
        else:
            flow = OAuth2WebServerFlow(client_id='495642085510-k0tmvj2m941jhre2nbqka17vqpjfddtd.apps.googleusercontent.com',
                                       client_secret='kOc9wMptUtxkcIFbtZCcrEAc',
                                       scope=scope,
                                       redirect_uri='urn:ietf:wg:oauth:2.0:oob')

        storage = Storage("gcs_credentials.dat")
        credentials = storage.get()

        if credentials is None or credentials.invalid:
            credentials = run(flow, storage)

        return credentials

    def upload(self, dataframe, bucket, name, project_id,
               to_gcs_kwargs=None, to_csv_kwargs=None):

        resumable = to_gcs_kwargs.get("resumable", False)

        string_df = StringIO()
        dataframe.to_csv(string_df, **to_csv_kwargs)

        media = MediaIoBaseUpload(string_df, mimetype='text/csv', **to_gcs_kwargs)

        request = self.service.objects().insert(bucket=bucket,
                                                name=name,
                                                media_body=media)

        if resumable:
            response = None
            while response is None:
                progress, response = request.next_chunk()

        else:
            request.execute()

    def download(self, bucket, name, project_id):

        request = (self.service.objects()
                        .get_media(bucket="th-code", object="pandas.csv")
                        .execute())

        in_memory_df = StringIO(request)

        df = pd.read_csv(in_memory_df)

        return df

def from_gcs(bucket, name, project_id, service_account, private_key):
    """
    Read a DataFrame from Google Cloud Storage.

    Parameters
    ----------
    bucket : string
        Bucket in GCS where the object resides
    name : string
        Name of the object, or regex matching the name of the object
    project_id : string
        ProjectId in google
    service_account : string
        Service account email
    private_key : string
        Path to private key file

    Returns
    -------
    dataframe : Dataframe
    """

    _test_imports()

    g = GCSConnector(project_id=project_id,
                     service_account=service_account,
                     private_key=private_key)

    return g.download(bucket, name, project_id)


def to_gcs(dataframe, bucket, name, project_id,
           service_account=None, private_key=None, to_gcs_kwargs=None, to_csv_kwargs=None):
    """
    Write a DataFrame to Google Cloud Storage.

    Parameters
    ----------
    dataframe : DataFrame
        DataFrame to be written
    bucket : string
        Bucket in GCS where the dataframe will be written
    name : string
        Object name in GCS
    project_id : string
        ProjectId in google
    service_account : string
        Service account email
    private_key : string
        Path to private key file
    to_gcs_kwargs : dict
        Dictionary of keywords passed directly to insert objects
    to_csv_kwargs : dict
        Dictionary of keywords passed  directly `to_csv`
    """

    _test_imports()

    if to_gcs_kwargs is None:
        to_gcs_kwargs = {}

    if to_csv_kwargs is None:
        to_csv_kwargs = {}

    gcs_connection = GCSConnector(project_id=project_id,
                                  service_account=service_account,
                                  private_key=private_key)

    gcs_connection.upload(dataframe=dataframe,
                          bucket=bucket,
                          name=name,
                          project_id=project_id,
                          to_gcs_kwargs=to_gcs_kwargs,
                          to_csv_kwargs=to_csv_kwargs)
