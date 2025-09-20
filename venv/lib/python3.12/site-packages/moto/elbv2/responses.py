from typing import Optional

from moto.core.exceptions import RESTError
from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.serialize import return_if_not_empty

from .exceptions import ListenerOrBalancerMissingError, TargetGroupNotFoundError
from .models import ELBv2Backend, elbv2_backends

SSL_POLICIES = [
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
            {"name": "AES128-GCM-SHA256", "priority": 13},
            {"name": "AES128-SHA256", "priority": 14},
            {"name": "AES128-SHA", "priority": 15},
            {"name": "AES256-GCM-SHA384", "priority": 16},
            {"name": "AES256-SHA256", "priority": 17},
            {"name": "AES256-SHA", "priority": 18},
        ],
        "name": "ELBSecurityPolicy-2016-08",
        "ssl_protocols": ["TLSv1", "TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 6},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 7},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 9},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 11},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-2-2021-06",
        "ssl_protocols": ["TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 6},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 7},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-2-Res-2021-06",
        "ssl_protocols": ["TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 6},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 7},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 9},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 11},
            {"name": "AES128-GCM-SHA256", "priority": 12},
            {"name": "AES128-SHA256", "priority": 13},
            {"name": "AES256-GCM-SHA384", "priority": 14},
            {"name": "AES256-SHA256", "priority": 15},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-2-Ext1-2021-06",
        "ssl_protocols": ["TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 6},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 7},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 8},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 9},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 12},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 13},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 14},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 15},
            {"name": "AES128-GCM-SHA256", "priority": 16},
            {"name": "AES128-SHA256", "priority": 17},
            {"name": "AES128-SHA", "priority": 18},
            {"name": "AES256-GCM-SHA384", "priority": 19},
            {"name": "AES256-SHA256", "priority": 20},
            {"name": "AES256-SHA", "priority": 21},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-2-Ext2-2021-06",
        "ssl_protocols": ["TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 6},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 7},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 8},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 9},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 12},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 13},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 14},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 15},
            {"name": "AES128-GCM-SHA256", "priority": 16},
            {"name": "AES128-SHA256", "priority": 17},
            {"name": "AES128-SHA", "priority": 18},
            {"name": "AES256-GCM-SHA384", "priority": 19},
            {"name": "AES256-SHA256", "priority": 20},
            {"name": "AES256-SHA", "priority": 21},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-1-2021-06",
        "ssl_protocols": ["TLSv1.1", "TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 4},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 5},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 6},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 7},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 8},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 9},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 12},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 13},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 14},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 15},
            {"name": "AES128-GCM-SHA256", "priority": 16},
            {"name": "AES128-SHA256", "priority": 17},
            {"name": "AES128-SHA", "priority": 18},
            {"name": "AES256-GCM-SHA384", "priority": 19},
            {"name": "AES256-SHA256", "priority": 20},
            {"name": "AES256-SHA", "priority": 21},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-0-2021-06",
        "ssl_protocols": ["TLSv1", "TLSv1.1", "TLSv1.2", "TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "TLS_AES_128_GCM_SHA256", "priority": 1},
            {"name": "TLS_AES_256_GCM_SHA384", "priority": 2},
            {"name": "TLS_CHACHA20_POLY1305_SHA256", "priority": 3},
        ],
        "name": "ELBSecurityPolicy-TLS13-1-3-2021-06",
        "ssl_protocols": ["TLSv1.3"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 5},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 8},
            {"name": "AES128-GCM-SHA256", "priority": 9},
            {"name": "AES128-SHA256", "priority": 10},
            {"name": "AES256-GCM-SHA384", "priority": 11},
            {"name": "AES256-SHA256", "priority": 12},
        ],
        "name": "ELBSecurityPolicy-TLS-1-2-2017-01",
        "ssl_protocols": ["TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
            {"name": "AES128-GCM-SHA256", "priority": 13},
            {"name": "AES128-SHA256", "priority": 14},
            {"name": "AES128-SHA", "priority": 15},
            {"name": "AES256-GCM-SHA384", "priority": 16},
            {"name": "AES256-SHA256", "priority": 17},
            {"name": "AES256-SHA", "priority": 18},
        ],
        "name": "ELBSecurityPolicy-TLS-1-1-2017-01",
        "ssl_protocols": ["TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
            {"name": "AES128-GCM-SHA256", "priority": 13},
            {"name": "AES128-SHA256", "priority": 14},
            {"name": "AES128-SHA", "priority": 15},
            {"name": "AES256-GCM-SHA384", "priority": 16},
            {"name": "AES256-SHA256", "priority": 17},
            {"name": "AES256-SHA", "priority": 18},
        ],
        "name": "ELBSecurityPolicy-TLS-1-2-Ext-2018-06",
        "ssl_protocols": ["TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
        ],
        "name": "ELBSecurityPolicy-FS-2018-06",
        "ssl_protocols": ["TLSv1", "TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
            {"name": "AES128-GCM-SHA256", "priority": 13},
            {"name": "AES128-SHA256", "priority": 14},
            {"name": "AES128-SHA", "priority": 15},
            {"name": "AES256-GCM-SHA384", "priority": 16},
            {"name": "AES256-SHA256", "priority": 17},
            {"name": "AES256-SHA", "priority": 18},
        ],
        "name": "ELBSecurityPolicy-2015-05",
        "ssl_protocols": ["TLSv1", "TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
            {"name": "AES128-GCM-SHA256", "priority": 13},
            {"name": "AES128-SHA256", "priority": 14},
            {"name": "AES128-SHA", "priority": 15},
            {"name": "AES256-GCM-SHA384", "priority": 16},
            {"name": "AES256-SHA256", "priority": 17},
            {"name": "AES256-SHA", "priority": 18},
            {"name": "DES-CBC3-SHA", "priority": 19},
        ],
        "name": "ELBSecurityPolicy-TLS-1-0-2015-04",
        "ssl_protocols": ["TLSv1", "TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 5},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 8},
        ],
        "name": "ELBSecurityPolicy-FS-1-2-Res-2019-08",
        "ssl_protocols": ["TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
        ],
        "name": "ELBSecurityPolicy-FS-1-1-2019-08",
        "ssl_protocols": ["TLSv1.1", "TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES128-SHA256", "priority": 3},
            {"name": "ECDHE-RSA-AES128-SHA256", "priority": 4},
            {"name": "ECDHE-ECDSA-AES128-SHA", "priority": 5},
            {"name": "ECDHE-RSA-AES128-SHA", "priority": 6},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 7},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 8},
            {"name": "ECDHE-ECDSA-AES256-SHA384", "priority": 9},
            {"name": "ECDHE-RSA-AES256-SHA384", "priority": 10},
            {"name": "ECDHE-RSA-AES256-SHA", "priority": 11},
            {"name": "ECDHE-ECDSA-AES256-SHA", "priority": 12},
        ],
        "name": "ELBSecurityPolicy-FS-1-2-2019-08",
        "ssl_protocols": ["TLSv1.2"],
    },
    {
        "ciphers": [
            {"name": "ECDHE-ECDSA-AES128-GCM-SHA256", "priority": 1},
            {"name": "ECDHE-RSA-AES128-GCM-SHA256", "priority": 2},
            {"name": "ECDHE-ECDSA-AES256-GCM-SHA384", "priority": 3},
            {"name": "ECDHE-RSA-AES256-GCM-SHA384", "priority": 4},
        ],
        "name": "ELBSecurityPolicy-FS-1-2-Res-2020-10",
        "ssl_protocols": ["TLSv1.2"],
    },
]


def transform_dict(data: dict[str, str]) -> list[dict[str, str]]:
    transformed = [{"Key": key, "Value": value} for key, value in data.items()]
    return transformed


def transform_certificates(data: list[str]) -> Optional[list[dict[str, str]]]:
    return [{"CertificateArn": cert} for cert in data] or None


class ELBV2Response(BaseResponse):
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "CreateListenerOutput.Listeners.Listener.Certificates": transform_certificates,
        "CreateTargetGroupOutput.TargetGroups.TargetGroup.LoadBalancerArns": return_if_not_empty,
        "DescribeListenerAttributesOutput.Attributes": transform_dict,
        "DescribeListenersOutput.Listeners.Listener.Certificates": transform_certificates,
        "DescribeLoadBalancerAttributesOutput.Attributes": transform_dict,
        "DescribeTagsOutput.TagDescriptions.TagDescription.Tags": transform_dict,
        "DescribeTargetGroupAttributesOutput.Attributes": transform_dict,
        "ModifyListenerAttributesOutput.Attributes": transform_dict,
        "ModifyListenerOutput.Listeners.Listener.Certificates": transform_certificates,
        "ModifyLoadBalancerAttributesOutput.Attributes": transform_dict,
        "ModifyTargetGroupAttributesOutput.Attributes": transform_dict,
    }

    def __init__(self) -> None:
        super().__init__(service_name="elbv2")
        self.automated_parameter_parsing = True

    @property
    def elbv2_backend(self) -> ELBv2Backend:
        return elbv2_backends[self.current_account][self.region]

    def create_load_balancer(self) -> ActionResult:
        params = self._get_params()
        load_balancer_name = params.get("Name")
        subnet_ids = params.get("Subnets", [])
        subnet_mappings = params.get("SubnetMappings", [])
        security_groups = params.get("SecurityGroups", [])
        scheme = params.get("Scheme")
        loadbalancer_type = params.get("Type")
        tags = params.get("Tags")
        load_balancer = self.elbv2_backend.create_load_balancer(
            name=load_balancer_name,  # type: ignore
            security_groups=security_groups,
            subnet_ids=subnet_ids,
            subnet_mappings=subnet_mappings,
            scheme=scheme,  # type: ignore
            loadbalancer_type=loadbalancer_type,
            tags=tags,
        )
        result = {"LoadBalancers": [load_balancer]}
        return ActionResult(result)

    def create_rule(self) -> ActionResult:
        params = self._get_params()
        rule = self.elbv2_backend.create_rule(
            listener_arn=params["ListenerArn"],
            conditions=params["Conditions"],
            priority=params["Priority"],
            actions=params["Actions"],
            tags=params.get("Tags"),
        )
        result = {"Rules": [rule]}
        return ActionResult(result)

    def create_target_group(self) -> ActionResult:
        params = self._get_params()
        name = params.get("Name")
        vpc_id = params.get("VpcId")
        protocol = params.get("Protocol")
        protocol_version = params.get("ProtocolVersion")
        port = params.get("Port")
        healthcheck_protocol = self._get_param("HealthCheckProtocol")
        healthcheck_port = self._get_param("HealthCheckPort")
        healthcheck_path = self._get_param("HealthCheckPath")
        healthcheck_interval_seconds = self._get_int_param("HealthCheckIntervalSeconds")
        healthcheck_timeout_seconds = self._get_int_param("HealthCheckTimeoutSeconds")
        healthcheck_enabled = self._get_param("HealthCheckEnabled")
        healthy_threshold_count = self._get_param("HealthyThresholdCount")
        unhealthy_threshold_count = self._get_param("UnhealthyThresholdCount")
        matcher = params.get("Matcher")
        target_type = params.get("TargetType", "instance")
        ip_address_type = params.get("IpAddressType")
        tags = params.get("Tags")

        target_group = self.elbv2_backend.create_target_group(
            name,  # type: ignore
            vpc_id=vpc_id,
            protocol=protocol,
            protocol_version=protocol_version,
            port=port,
            healthcheck_protocol=healthcheck_protocol,
            healthcheck_port=healthcheck_port,
            healthcheck_path=healthcheck_path,
            healthcheck_interval_seconds=healthcheck_interval_seconds,
            healthcheck_timeout_seconds=healthcheck_timeout_seconds,
            healthcheck_enabled=healthcheck_enabled,
            healthy_threshold_count=healthy_threshold_count,
            unhealthy_threshold_count=unhealthy_threshold_count,
            matcher=matcher,
            ip_address_type=ip_address_type,
            target_type=target_type,
            tags=tags,
        )
        result = {"TargetGroups": [target_group]}
        return ActionResult(result)

    def create_listener(self) -> ActionResult:
        params = self._get_params()
        load_balancer_arn = self._get_param("LoadBalancerArn")
        protocol = self._get_param("Protocol")
        port = self._get_param("Port")
        ssl_policy = self._get_param("SslPolicy", "ELBSecurityPolicy-2016-08")
        certificates = self._get_param("Certificates", [])
        if certificates:
            certificate = certificates[0].get("certificate_arn")
        else:
            certificate = None
        default_actions = params.get("DefaultActions", [])
        alpn_policy = params.get("AlpnPolicy", [])
        tags = params.get("Tags")

        listener = self.elbv2_backend.create_listener(
            load_balancer_arn=load_balancer_arn,
            protocol=protocol,
            port=port,
            ssl_policy=ssl_policy,
            certificate=certificate,
            actions=default_actions,
            alpn_policy=alpn_policy,
            tags=tags,
        )

        result = {"Listeners": [listener]}
        return ActionResult(result)

    def describe_load_balancers(self) -> ActionResult:
        arns = self._get_param("LoadBalancerArns", [])
        names = self._get_param("Names", [])
        all_load_balancers = list(
            self.elbv2_backend.describe_load_balancers(arns, names)
        )
        marker = self._get_param("Marker")
        all_names = [balancer.name for balancer in all_load_balancers]
        if marker:
            start = all_names.index(marker) + 1
        else:
            start = 0
        page_size = self._get_int_param(
            "PageSize", 50
        )  # the default is 400, but using 50 to make testing easier
        load_balancers_resp = all_load_balancers[start : start + page_size]
        next_marker = None
        if len(all_load_balancers) > start + page_size:
            next_marker = load_balancers_resp[-1].name
        result = {"LoadBalancers": load_balancers_resp, "NextMarker": next_marker}
        return ActionResult(result)

    def describe_rules(self) -> ActionResult:
        listener_arn = self._get_param("ListenerArn")
        rule_arns = self._get_param("RuleArns")
        all_rules = list(self.elbv2_backend.describe_rules(listener_arn, rule_arns))
        all_arns = [rule.arn for rule in all_rules]
        page_size = self._get_int_param("PageSize", 50)  # set 50 for temporary

        marker = self._get_param("Marker")
        if marker:
            start = all_arns.index(marker) + 1
        else:
            start = 0
        rules_resp = all_rules[start : start + page_size]
        next_marker = None

        if len(all_rules) > start + page_size:
            next_marker = rules_resp[-1].arn
        result = {"Rules": rules_resp, "NextMarker": next_marker}
        return ActionResult(result)

    def describe_target_groups(self) -> ActionResult:
        load_balancer_arn = self._get_param("LoadBalancerArn")
        target_group_arns = self._get_param("TargetGroupArns", [])
        names = self._get_param("Names", [])
        target_groups = self.elbv2_backend.describe_target_groups(
            load_balancer_arn, target_group_arns, names
        )
        result = {"TargetGroups": target_groups}
        return ActionResult(result)

    def describe_target_group_attributes(self) -> ActionResult:
        target_group_arn = self._get_param("TargetGroupArn")
        target_group = self.elbv2_backend.target_groups.get(target_group_arn)
        if not target_group:
            raise TargetGroupNotFoundError()
        result = {"Attributes": target_group.attributes}
        return ActionResult(result)

    def describe_listeners(self) -> ActionResult:
        load_balancer_arn = self._get_param("LoadBalancerArn")
        listener_arns = self._get_param("ListenerArns", [])
        if not load_balancer_arn and not listener_arns:
            raise ListenerOrBalancerMissingError()

        listeners = self.elbv2_backend.describe_listeners(
            load_balancer_arn, listener_arns
        )
        result = {"Listeners": listeners}
        return ActionResult(result)

    def delete_load_balancer(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        self.elbv2_backend.delete_load_balancer(arn)
        return EmptyResult()

    def delete_rule(self) -> ActionResult:
        arn = self._get_param("RuleArn")
        self.elbv2_backend.delete_rule(arn)
        return EmptyResult()

    def delete_target_group(self) -> ActionResult:
        arn = self._get_param("TargetGroupArn")
        self.elbv2_backend.delete_target_group(arn)
        return EmptyResult()

    def delete_listener(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        self.elbv2_backend.delete_listener(arn)
        return EmptyResult()

    def modify_rule(self) -> ActionResult:
        rule_arn = self._get_param("RuleArn")
        params = self._get_params()
        conditions = params.get("Conditions", [])
        actions = params.get("Actions", [])
        rule = self.elbv2_backend.modify_rule(
            rule_arn=rule_arn, conditions=conditions, actions=actions
        )
        result = {"Rules": [rule]}
        return ActionResult(result)

    def modify_target_group_attributes(self) -> ActionResult:
        target_group_arn = self._get_param("TargetGroupArn")
        attrs = self._get_param("Attributes", [])
        attributes = {attr["key"]: attr["value"] for attr in attrs}
        self.elbv2_backend.modify_target_group_attributes(target_group_arn, attributes)
        result = {"Attributes": attributes}
        return ActionResult(result)

    def register_targets(self) -> ActionResult:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_param("Targets", [])
        self.elbv2_backend.register_targets(target_group_arn, targets)
        return EmptyResult()

    def deregister_targets(self) -> ActionResult:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_param("Targets", [])
        self.elbv2_backend.deregister_targets(target_group_arn, targets)
        return EmptyResult()

    def describe_target_health(self) -> ActionResult:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_param("Targets", [])
        target_health_descriptions = self.elbv2_backend.describe_target_health(
            target_group_arn, targets
        )
        result = {"TargetHealthDescriptions": target_health_descriptions}
        return ActionResult(result)

    def set_rule_priorities(self) -> ActionResult:
        rule_priorities = self._get_param("RulePriorities", [])
        rules = self.elbv2_backend.set_rule_priorities(rule_priorities)
        result = {"Rules": rules}
        return ActionResult(result)

    def add_tags(self) -> ActionResult:
        resource_arns = self._get_param("ResourceArns", [])
        tags = self._get_param("Tags", [])
        self.elbv2_backend.add_tags(resource_arns, tags)  # type: ignore
        return EmptyResult()

    def remove_tags(self) -> ActionResult:
        resource_arns = self._get_param("ResourceArns", [])
        tag_keys = self._get_param("TagKeys", [])
        self.elbv2_backend.remove_tags(resource_arns, tag_keys)
        return EmptyResult()

    def describe_tags(self) -> ActionResult:
        resource_arns = self._get_param("ResourceArns", [])
        resource_tags = self.elbv2_backend.describe_tags(resource_arns)
        result = {
            "TagDescriptions": [
                {"ResourceArn": arn, "Tags": tags}
                for arn, tags in resource_tags.items()
            ]
        }
        return ActionResult(result)

    def describe_account_limits(self) -> ActionResult:
        # Supports paging but not worth implementing yet
        # marker = self._get_param('Marker')
        # page_size = self._get_int_param('PageSize')

        limit_data = {
            "application-load-balancers": 20,
            "target-groups": 3000,
            "targets-per-application-load-balancer": 30,
            "listeners-per-application-load-balancer": 50,
            "rules-per-application-load-balancer": 100,
            "network-load-balancers": 20,
            "targets-per-network-load-balancer": 200,
            "listeners-per-network-load-balancer": 50,
            "certificates-per-application-load-balancer": 25,
        }
        limits = [{"Name": k, "Max": v} for k, v in limit_data.items()]
        result = {"Limits": limits}
        return ActionResult(result)

    def describe_ssl_policies(self) -> ActionResult:
        names = self._get_param("Names", [])
        # Supports paging but not worth implementing yet
        # marker = self._get_param('Marker')
        # page_size = self._get_int_param('PageSize')

        policies = SSL_POLICIES
        if names:
            policies = filter(lambda policy: policy["name"] in names, policies)  # type: ignore

        result = {"SslPolicies": policies}
        return ActionResult(result)

    def set_ip_address_type(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        ip_type = self._get_param("IpAddressType")
        self.elbv2_backend.set_ip_address_type(arn, ip_type)
        result = {"IpAddressType": ip_type}
        return ActionResult(result)

    def set_security_groups(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        sec_groups = self._get_param("SecurityGroups", [])
        self.elbv2_backend.set_security_groups(arn, sec_groups)
        result = {"SecurityGroups": sec_groups}
        return ActionResult(result)

    def set_subnets(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        subnets = self._get_param("Subnets", [])
        subnet_mappings = self._get_param("SubnetMappings", [])
        subnet_zone_list = self.elbv2_backend.set_subnets(arn, subnets, subnet_mappings)
        result = {"AvailabilityZones": subnet_zone_list}
        return ActionResult(result)

    def modify_load_balancer_attributes(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        attrs = self._get_param("Attributes", [])
        attrs = {attr["Key"]: attr["Value"] for attr in attrs}
        all_attrs = self.elbv2_backend.modify_load_balancer_attributes(arn, attrs)
        result = {"Attributes": all_attrs}
        return ActionResult(result)

    def describe_load_balancer_attributes(self) -> ActionResult:
        arn = self._get_param("LoadBalancerArn")
        attrs = self.elbv2_backend.describe_load_balancer_attributes(arn)
        result = {"Attributes": attrs}
        return ActionResult(result)

    def modify_target_group(self) -> ActionResult:
        arn = self._get_param("TargetGroupArn")

        health_check_proto = self._get_param(
            "HealthCheckProtocol"
        )  # 'HTTP' | 'HTTPS' | 'TCP',
        health_check_port = self._get_param("HealthCheckPort")
        health_check_path = self._get_param("HealthCheckPath")
        health_check_interval = self._get_int_param("HealthCheckIntervalSeconds")
        health_check_timeout = self._get_int_param("HealthCheckTimeoutSeconds")
        health_check_enabled = self._get_param("HealthCheckEnabled")
        healthy_threshold_count = self._get_param("HealthyThresholdCount")
        unhealthy_threshold_count = self._get_param("UnhealthyThresholdCount")
        http_codes = self._get_param("Matcher.HttpCode")
        target_group = self.elbv2_backend.modify_target_group(
            arn,
            health_check_proto,
            health_check_port,
            health_check_path,
            health_check_interval,
            health_check_timeout,
            healthy_threshold_count,
            unhealthy_threshold_count,
            http_codes,
            health_check_enabled=health_check_enabled,
        )
        result = {"TargetGroups": [target_group]}
        return ActionResult(result)

    def modify_listener(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        port = self._get_param("Port")
        protocol = self._get_param("Protocol")
        ssl_policy = self._get_param("SslPolicy")
        certificates = self._get_param("Certificates", [])
        default_actions = self._get_param("DefaultActions", [])

        # Should really move SSL Policies to models
        if ssl_policy is not None and ssl_policy not in [
            item["name"] for item in SSL_POLICIES
        ]:
            raise RESTError("SSLPolicyNotFound", f"Policy {ssl_policy} not found")

        listener = self.elbv2_backend.modify_listener(
            arn, port, protocol, ssl_policy, certificates, default_actions
        )

        result = {"Listeners": [listener]}
        return ActionResult(result)

    def add_listener_certificates(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        certificates = self._get_param("Certificates", [])
        certificate_arns = self.elbv2_backend.add_listener_certificates(
            arn, certificates
        )
        result = {"Certificates": [{"CertificateArn": arn} for arn in certificate_arns]}
        return ActionResult(result)

    def describe_listener_certificates(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        certificates = self.elbv2_backend.describe_listener_certificates(arn)
        result = {"Certificates": [{"CertificateArn": arn} for arn in certificates]}
        return ActionResult(result)

    def remove_listener_certificates(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        certificates = self._get_param("Certificates", [])
        self.elbv2_backend.remove_listener_certificates(arn, certificates)
        return EmptyResult()

    def describe_listener_attributes(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        attrs = self.elbv2_backend.describe_listener_attributes(arn)
        result = {"Attributes": attrs}
        return ActionResult(result)

    def modify_listener_attributes(self) -> ActionResult:
        arn = self._get_param("ListenerArn")
        attrs = self._get_params()["Attributes"]
        updated_attrs = self.elbv2_backend.modify_listener_attributes(
            listener_arn=arn, attrs=attrs
        )
        result = {"Attributes": updated_attrs}
        return ActionResult(result)

    def describe_capacity_reservation(self) -> ActionResult:
        result = {"CapacityReservationState": [{"State": {"Code": "provisioned"}}]}
        return ActionResult(result)
