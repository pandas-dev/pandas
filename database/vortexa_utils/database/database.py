# @Author: richard
# @Date:   2018-12-04T17:58:19+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T17:58:19+00:00
import json
import logging
import os
from dataclasses import dataclass, field
from typing import Dict

import boto3
import requests
from sqlalchemy import create_engine

logger = logging.getLogger(__name__)

secretsmanager = boto3.client("secretsmanager")

DEFAULT_CERT_URL = (
    "https://s3.amazonaws.com/rds-downloads/rds-combined-ca-bundle.pem"
)
DEFAULT_CERT_PATH = "/tmp/vortexa_utils_py/rds/ca-bundle.pem"

DEFAULT_CREDENTIAL = "rds/dev/default"
DEFAULT_CREDENTIAL_MAPPING = dict(
    host="host", username="user", port="port", password="password"
)


@dataclass
class DatabaseFactory(object):
    """DatabaseFactory Class.

    Class for createing a database engine factory.

    usage::

        factory = DatabaseFactory()
        engine = factory.engine()

    Parameters
    ----------
    secret_id : str
        `secret_id` of the database credential.
        (the default is 'rds/dev/default' wich points to the dev database host)
    cert_file : str
        The location to store the ssl certificate file
    cert_url : str
        The url to fetch the aws rds ssl certificates from
    credential_mapping : Dict[str, str]
        A mapping between the `psycopg` connection args and the credential keys
    """

    secret_id: str = DEFAULT_CREDENTIAL
    cert_file: str = DEFAULT_CERT_PATH
    cert_url: str = DEFAULT_CERT_URL
    credential_mapping: Dict[str, str] = field(
        default_factory=lambda: dict(DEFAULT_CREDENTIAL_MAPPING)
    )

    def __post_init__(self):
        logger.debug(f"Created {self.secret_id} factory object")

    def fetch_cert(self, force: bool = False):
        if not os.path.isfile(self.cert_file) or force:
            logger.info("getting cert")
            os.makedirs(os.path.dirname(self.cert_file), exist_ok=True)
            cert = requests.get(self.cert_url)
            with open(self.cert_file, "w") as f:
                f.write(cert.text)
        return self.cert_file

    def get_credential(self):
        secret = secretsmanager.get_secret_value(SecretId=self.secret_id)
        return json.loads(secret["SecretString"])

    def engine(self, dbname: str = None, echo: bool = False, **kwargs):
        # type (...) -> sqlalchemy.engine.Engine
        """`sqlalchemy.engine.Engine` instance factory.

        Parameters
        ----------
        dbname : str
            database name `dbname` to connect to.
            (the default is `None`, which will use the dbname in the secret
             credential).
        echo : bool
             `echo` (the default is False).

        Returns
        -------
        sqlalchemy.engine.Engine
            SQLalchemy connection engine

        Examples
        -------
        >>> factory = DatabaseFactory()
        >>> engine = factory.engine()

        """
        cert_filename = self.fetch_cert()
        credential = self.get_credential()
        connect_args = {
            v: credential[k] for k, v in self.credential_mapping.items()
        }

        dbname = dbname or os.environ.get("DBNAME") or credential["dbname"]
        host = connect_args.pop("host")
        port = connect_args.pop("port")

        connect_args.update(sslmode="verify-full", sslrootcert=cert_filename)
        engine = create_engine(
            f"postgresql://{host}:{port}/{dbname}",
            echo=echo,
            connect_args=connect_args,
            **kwargs,
        )
        return engine
