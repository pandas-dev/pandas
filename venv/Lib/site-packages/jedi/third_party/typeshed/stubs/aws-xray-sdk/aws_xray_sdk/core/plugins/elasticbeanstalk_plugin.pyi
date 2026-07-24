from logging import Logger
from typing import Final

log: Logger
CONF_PATH: Final = "/var/elasticbeanstalk/xray/environment.conf"
SERVICE_NAME: Final = "elastic_beanstalk"
ORIGIN: Final = "AWS::ElasticBeanstalk::Environment"

def initialize() -> None: ...
