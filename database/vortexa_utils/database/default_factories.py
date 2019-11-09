from dataclasses import dataclass

from .database import DatabaseFactory


@dataclass
class DevFactory(DatabaseFactory):
    secret_id: str = "rds/dev/default"


@dataclass
class ProdFactory(DatabaseFactory):
    secret_id: str = "rds/prod/default"


@dataclass
class RedFactory(DatabaseFactory):
    cert_url: str = "https://s3.amazonaws.com/redshift-downloads/redshift-ca-bundle.crt"
    cert_file: str = "/tmp/vortexa_utils_py/rds/redshift-ca-bundle.pem"
    secret_id: str = "redshift/prod/default"
