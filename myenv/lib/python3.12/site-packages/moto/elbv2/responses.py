from moto.core.exceptions import RESTError
from moto.core.responses import BaseResponse

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


class ELBV2Response(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="elbv2")

    @property
    def elbv2_backend(self) -> ELBv2Backend:
        return elbv2_backends[self.current_account][self.region]

    def create_load_balancer(self) -> str:
        params = self._get_params()
        load_balancer_name = params.get("Name")
        subnet_ids = self._get_multi_param("Subnets.member")
        subnet_mappings = params.get("SubnetMappings", [])
        security_groups = self._get_multi_param("SecurityGroups.member")
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
        template = self.response_template(CREATE_LOAD_BALANCER_TEMPLATE)
        return template.render(load_balancer=load_balancer)

    def create_rule(self) -> str:
        params = self._get_params()
        rules = self.elbv2_backend.create_rule(
            listener_arn=params["ListenerArn"],
            conditions=params["Conditions"],
            priority=params["Priority"],
            actions=params["Actions"],
            tags=params.get("Tags"),
        )

        template = self.response_template(CREATE_RULE_TEMPLATE)
        return template.render(rules=rules)

    def create_target_group(self) -> str:
        params = self._get_params()
        name = params.get("Name")
        vpc_id = params.get("VpcId")
        protocol = params.get("Protocol")
        protocol_version = params.get("ProtocolVersion")
        port = params.get("Port")
        healthcheck_protocol = self._get_param("HealthCheckProtocol")
        healthcheck_port = self._get_param("HealthCheckPort")
        healthcheck_path = self._get_param("HealthCheckPath")
        healthcheck_interval_seconds = self._get_param("HealthCheckIntervalSeconds")
        healthcheck_timeout_seconds = self._get_param("HealthCheckTimeoutSeconds")
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

        template = self.response_template(CREATE_TARGET_GROUP_TEMPLATE)
        return template.render(target_group=target_group)

    def create_listener(self) -> str:
        params = self._get_params()
        load_balancer_arn = self._get_param("LoadBalancerArn")
        protocol = self._get_param("Protocol")
        port = self._get_param("Port")
        ssl_policy = self._get_param("SslPolicy", "ELBSecurityPolicy-2016-08")
        certificates = self._get_list_prefix("Certificates.member")
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

        template = self.response_template(CREATE_LISTENER_TEMPLATE)
        return template.render(listener=listener)

    def describe_load_balancers(self) -> str:
        arns = self._get_multi_param("LoadBalancerArns.member")
        names = self._get_multi_param("Names.member")
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

        template = self.response_template(DESCRIBE_LOAD_BALANCERS_TEMPLATE)
        return template.render(load_balancers=load_balancers_resp, marker=next_marker)

    def describe_rules(self) -> str:
        listener_arn = self._get_param("ListenerArn")
        rule_arns = (
            self._get_multi_param("RuleArns.member")
            if any(
                k
                for k in list(self.querystring.keys())
                if k.startswith("RuleArns.member")
            )
            else None
        )
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
        template = self.response_template(DESCRIBE_RULES_TEMPLATE)
        return template.render(rules=rules_resp, marker=next_marker)

    def describe_target_groups(self) -> str:
        load_balancer_arn = self._get_param("LoadBalancerArn")
        target_group_arns = self._get_multi_param("TargetGroupArns.member")
        names = self._get_multi_param("Names.member")

        target_groups = self.elbv2_backend.describe_target_groups(
            load_balancer_arn, target_group_arns, names
        )
        template = self.response_template(DESCRIBE_TARGET_GROUPS_TEMPLATE)
        return template.render(target_groups=target_groups)

    def describe_target_group_attributes(self) -> str:
        target_group_arn = self._get_param("TargetGroupArn")
        target_group = self.elbv2_backend.target_groups.get(target_group_arn)
        if not target_group:
            raise TargetGroupNotFoundError()
        template = self.response_template(DESCRIBE_TARGET_GROUP_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=target_group.attributes)

    def describe_listeners(self) -> str:
        load_balancer_arn = self._get_param("LoadBalancerArn")
        listener_arns = self._get_multi_param("ListenerArns.member")
        if not load_balancer_arn and not listener_arns:
            raise ListenerOrBalancerMissingError()

        listeners = self.elbv2_backend.describe_listeners(
            load_balancer_arn, listener_arns
        )
        template = self.response_template(DESCRIBE_LISTENERS_TEMPLATE)
        return template.render(listeners=listeners)

    def delete_load_balancer(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        self.elbv2_backend.delete_load_balancer(arn)
        template = self.response_template(DELETE_LOAD_BALANCER_TEMPLATE)
        return template.render()

    def delete_rule(self) -> str:
        arn = self._get_param("RuleArn")
        self.elbv2_backend.delete_rule(arn)
        template = self.response_template(DELETE_RULE_TEMPLATE)
        return template.render()

    def delete_target_group(self) -> str:
        arn = self._get_param("TargetGroupArn")
        self.elbv2_backend.delete_target_group(arn)
        template = self.response_template(DELETE_TARGET_GROUP_TEMPLATE)
        return template.render()

    def delete_listener(self) -> str:
        arn = self._get_param("ListenerArn")
        self.elbv2_backend.delete_listener(arn)
        template = self.response_template(DELETE_LISTENER_TEMPLATE)
        return template.render()

    def modify_rule(self) -> str:
        rule_arn = self._get_param("RuleArn")
        params = self._get_params()
        conditions = params.get("Conditions", [])
        actions = params.get("Actions", [])
        rules = self.elbv2_backend.modify_rule(
            rule_arn=rule_arn, conditions=conditions, actions=actions
        )
        template = self.response_template(MODIFY_RULE_TEMPLATE)
        return template.render(rules=rules)

    def modify_target_group_attributes(self) -> str:
        target_group_arn = self._get_param("TargetGroupArn")
        attrs = self._get_list_prefix("Attributes.member")
        attributes = {attr["key"]: attr["value"] for attr in attrs}
        self.elbv2_backend.modify_target_group_attributes(target_group_arn, attributes)

        template = self.response_template(MODIFY_TARGET_GROUP_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=attributes)

    def register_targets(self) -> str:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_list_prefix("Targets.member")
        self.elbv2_backend.register_targets(target_group_arn, targets)

        template = self.response_template(REGISTER_TARGETS_TEMPLATE)
        return template.render()

    def deregister_targets(self) -> str:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_list_prefix("Targets.member")
        self.elbv2_backend.deregister_targets(target_group_arn, targets)

        template = self.response_template(DEREGISTER_TARGETS_TEMPLATE)
        return template.render()

    def describe_target_health(self) -> str:
        target_group_arn = self._get_param("TargetGroupArn")
        targets = self._get_list_prefix("Targets.member")
        target_health_descriptions = self.elbv2_backend.describe_target_health(
            target_group_arn, targets
        )

        template = self.response_template(DESCRIBE_TARGET_HEALTH_TEMPLATE)
        return template.render(target_health_descriptions=target_health_descriptions)

    def set_rule_priorities(self) -> str:
        rule_priorities = self._get_list_prefix("RulePriorities.member")
        for rule_priority in rule_priorities:
            rule_priority["priority"] = int(rule_priority["priority"])
        rules = self.elbv2_backend.set_rule_priorities(rule_priorities)
        template = self.response_template(SET_RULE_PRIORITIES_TEMPLATE)
        return template.render(rules=rules)

    def add_tags(self) -> str:
        resource_arns = self._get_multi_param("ResourceArns.member")
        tags = self._get_params().get("Tags")

        self.elbv2_backend.add_tags(resource_arns, tags)  # type: ignore

        template = self.response_template(ADD_TAGS_TEMPLATE)
        return template.render()

    def remove_tags(self) -> str:
        resource_arns = self._get_multi_param("ResourceArns.member")
        tag_keys = self._get_multi_param("TagKeys.member")

        self.elbv2_backend.remove_tags(resource_arns, tag_keys)

        template = self.response_template(REMOVE_TAGS_TEMPLATE)
        return template.render()

    def describe_tags(self) -> str:
        resource_arns = self._get_multi_param("ResourceArns.member")
        resource_tags = self.elbv2_backend.describe_tags(resource_arns)

        template = self.response_template(DESCRIBE_TAGS_TEMPLATE)
        return template.render(resource_tags=resource_tags)

    def describe_account_limits(self) -> str:
        # Supports paging but not worth implementing yet
        # marker = self._get_param('Marker')
        # page_size = self._get_int_param('PageSize')

        limits = {
            "application-load-balancers": 20,
            "target-groups": 3000,
            "targets-per-application-load-balancer": 30,
            "listeners-per-application-load-balancer": 50,
            "rules-per-application-load-balancer": 100,
            "network-load-balancers": 20,
            "targets-per-network-load-balancer": 200,
            "listeners-per-network-load-balancer": 50,
        }

        template = self.response_template(DESCRIBE_LIMITS_TEMPLATE)
        return template.render(limits=limits)

    def describe_ssl_policies(self) -> str:
        names = self._get_multi_param("Names.member.")
        # Supports paging but not worth implementing yet
        # marker = self._get_param('Marker')
        # page_size = self._get_int_param('PageSize')

        policies = SSL_POLICIES
        if names:
            policies = filter(lambda policy: policy["name"] in names, policies)  # type: ignore

        template = self.response_template(DESCRIBE_SSL_POLICIES_TEMPLATE)
        return template.render(policies=policies)

    def set_ip_address_type(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        ip_type = self._get_param("IpAddressType")

        self.elbv2_backend.set_ip_address_type(arn, ip_type)

        template = self.response_template(SET_IP_ADDRESS_TYPE_TEMPLATE)
        return template.render(ip_type=ip_type)

    def set_security_groups(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        sec_groups = self._get_multi_param("SecurityGroups.member.")

        self.elbv2_backend.set_security_groups(arn, sec_groups)

        template = self.response_template(SET_SECURITY_GROUPS_TEMPLATE)
        return template.render(sec_groups=sec_groups)

    def set_subnets(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        subnets = self._get_multi_param("Subnets.member.")
        subnet_mappings = self._get_params().get("SubnetMappings", [])

        subnet_zone_list = self.elbv2_backend.set_subnets(arn, subnets, subnet_mappings)

        template = self.response_template(SET_SUBNETS_TEMPLATE)
        return template.render(subnets=subnet_zone_list)

    def modify_load_balancer_attributes(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        attrs = self._get_map_prefix(
            "Attributes.member", key_end="Key", value_end="Value"
        )

        all_attrs = self.elbv2_backend.modify_load_balancer_attributes(arn, attrs)

        template = self.response_template(MODIFY_LOADBALANCER_ATTRS_TEMPLATE)
        return template.render(attrs=all_attrs)

    def describe_load_balancer_attributes(self) -> str:
        arn = self._get_param("LoadBalancerArn")
        attrs = self.elbv2_backend.describe_load_balancer_attributes(arn)

        template = self.response_template(DESCRIBE_LOADBALANCER_ATTRS_TEMPLATE)
        return template.render(attrs=attrs)

    def modify_target_group(self) -> str:
        arn = self._get_param("TargetGroupArn")

        health_check_proto = self._get_param(
            "HealthCheckProtocol"
        )  # 'HTTP' | 'HTTPS' | 'TCP',
        health_check_port = self._get_param("HealthCheckPort")
        health_check_path = self._get_param("HealthCheckPath")
        health_check_interval = self._get_param("HealthCheckIntervalSeconds")
        health_check_timeout = self._get_param("HealthCheckTimeoutSeconds")
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

        template = self.response_template(MODIFY_TARGET_GROUP_TEMPLATE)
        return template.render(target_group=target_group)

    def modify_listener(self) -> str:
        arn = self._get_param("ListenerArn")
        port = self._get_param("Port")
        protocol = self._get_param("Protocol")
        ssl_policy = self._get_param("SslPolicy")
        certificates = self._get_list_prefix("Certificates.member")
        default_actions = self._get_params().get("DefaultActions", [])

        # Should really move SSL Policies to models
        if ssl_policy is not None and ssl_policy not in [
            item["name"] for item in SSL_POLICIES
        ]:
            raise RESTError("SSLPolicyNotFound", f"Policy {ssl_policy} not found")

        listener = self.elbv2_backend.modify_listener(
            arn, port, protocol, ssl_policy, certificates, default_actions
        )

        template = self.response_template(MODIFY_LISTENER_TEMPLATE)
        return template.render(listener=listener)

    def add_listener_certificates(self) -> str:
        arn = self._get_param("ListenerArn")
        certificates = self._get_list_prefix("Certificates.member")
        certificate_arns = self.elbv2_backend.add_listener_certificates(
            arn, certificates
        )

        template = self.response_template(ADD_LISTENER_CERTIFICATES_TEMPLATE)
        return template.render(certificates=certificate_arns)

    def describe_listener_certificates(self) -> str:
        arn = self._get_param("ListenerArn")
        certificates = self.elbv2_backend.describe_listener_certificates(arn)

        template = self.response_template(DESCRIBE_LISTENER_CERTIFICATES_TEMPLATE)
        return template.render(certificates=certificates)

    def remove_listener_certificates(self) -> str:
        arn = self._get_param("ListenerArn")
        certificates = self._get_list_prefix("Certificates.member")
        self.elbv2_backend.remove_listener_certificates(arn, certificates)

        template = self.response_template(REMOVE_LISTENER_CERTIFICATES_TEMPLATE)
        return template.render()


ADD_TAGS_TEMPLATE = """<AddTagsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <AddTagsResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</AddTagsResponse>"""

REMOVE_TAGS_TEMPLATE = """<RemoveTagsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <RemoveTagsResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</RemoveTagsResponse>"""

DESCRIBE_TAGS_TEMPLATE = """<DescribeTagsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeTagsResult>
    <TagDescriptions>
      {% for resource_arn, tags in resource_tags.items() %}
      <member>
        <ResourceArn>{{ resource_arn }}</ResourceArn>
        <Tags>
          {% for key, value in tags.items() %}
          <member>
            <Value>{{ value }}</Value>
            <Key>{{ key }}</Key>
          </member>
          {% endfor %}
        </Tags>
      </member>
      {% endfor %}
    </TagDescriptions>
  </DescribeTagsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeTagsResponse>"""

CREATE_LOAD_BALANCER_TEMPLATE = """<CreateLoadBalancerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <CreateLoadBalancerResult>
    <LoadBalancers>
      <member>
        <LoadBalancerArn>{{ load_balancer.arn }}</LoadBalancerArn>
        <Scheme>{{ load_balancer.scheme }}</Scheme>
        <LoadBalancerName>{{ load_balancer.name }}</LoadBalancerName>
        <VpcId>{{ load_balancer.vpc_id }}</VpcId>
        <CanonicalHostedZoneId>Z2P70J7EXAMPLE</CanonicalHostedZoneId>
        <CreatedTime>{{ load_balancer.created_time }}</CreatedTime>
        <AvailabilityZones>
          {% for subnet in load_balancer.subnets %}
          <member>
            <SubnetId>{{ subnet.id }}</SubnetId>
            <ZoneName>{{ subnet.availability_zone }}</ZoneName>
          </member>
          {% endfor %}
        </AvailabilityZones>
        <SecurityGroups>
          {% for security_group in load_balancer.security_groups %}
          <member>{{ security_group }}</member>
          {% endfor %}
        </SecurityGroups>
        <DNSName>{{ load_balancer.dns_name }}</DNSName>
        <State>
          <Code>{{ load_balancer.state }}</Code>
        </State>
        <Type>{{ load_balancer.loadbalancer_type }}</Type>
      </member>
    </LoadBalancers>
  </CreateLoadBalancerResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateLoadBalancerResponse>"""

CREATE_RULE_TEMPLATE = """<CreateRuleResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <CreateRuleResult>
    <Rules>
      <member>
        <IsDefault>{{ "true" if rules.is_default else "false" }}</IsDefault>
        <Conditions>
          {% for condition in rules.conditions %}
          <member>
            <Field>{{ condition["Field"] }}</Field>
            {% if "Values" in condition %}
            <Values>
              {% for value in condition["Values"] %}
              <member>{{ value }}</member>
              {% endfor %}
            </Values>
            {% endif %}
            {% if "HttpHeaderConfig" in condition %}
            <HttpHeaderConfig>
              <HttpHeaderName>{{ condition["HttpHeaderConfig"]["HttpHeaderName"] }}</HttpHeaderName>
              <Values>
                {% for value in condition["HttpHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpHeaderConfig>
            {% endif %}
            {% if "HttpRequestMethodConfig" in condition %}
            <HttpRequestMethodConfig>
              <Values>
                {% for value in condition["HttpRequestMethodConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpRequestMethodConfig>
            {% endif %}
            {% if "QueryStringConfig" in condition %}
            <QueryStringConfig>
              <Values>
                {% for value in condition["QueryStringConfig"]["Values"] %}
                <member>
                    <Key>{{ value["Key"] }}</Key>
                    <Value>{{ value["Value"] }}</Value>
                </member>
                {% endfor %}
              </Values>
            </QueryStringConfig>
            {% endif %}
            {% if "SourceIpConfig" in condition %}
            <SourceIpConfig>
              <Values>
                {% for value in condition["SourceIpConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </SourceIpConfig>
            {% endif %}
            {% if "PathPatternConfig" in condition %}
            <PathPatternConfig>
              <Values>
                {% for value in condition["PathPatternConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </PathPatternConfig>
            {% endif %}
            {% if "HostHeaderConfig" in condition %}
            <HostHeaderConfig>
              <Values>
                {% for value in condition["HostHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HostHeaderConfig>
            {% endif %}
          </member>
          {% endfor %}
        </Conditions>
        <Priority>{{ rules.priority }}</Priority>
        <RuleArn>{{ rules.arn }}</RuleArn>
        <Actions>
          {% for action in rules.actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </Actions>
      </member>
    </Rules>
  </CreateRuleResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateRuleResponse>"""

CREATE_TARGET_GROUP_TEMPLATE = """<CreateTargetGroupResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <CreateTargetGroupResult>
    <TargetGroups>
      <member>
        <TargetGroupArn>{{ target_group.arn }}</TargetGroupArn>
        <TargetGroupName>{{ target_group.name }}</TargetGroupName>
        {% if target_group.protocol %}
        <Protocol>{{ target_group.protocol }}</Protocol>
        {% if target_group.protocol_version %}
        <ProtocolVersion>{{ target_group.protocol_version }}</ProtocolVersion>
        {% endif %}
        {% endif %}
        {% if target_group.port %}
        <Port>{{ target_group.port }}</Port>
        {% endif %}
        {% if target_group.vpc_id %}
        <VpcId>{{ target_group.vpc_id }}</VpcId>
        {% endif %}
        {% if target_group.healthcheck_enabled %}
        {% if target_group.healthcheck_port %}
        <HealthCheckPort>{{ target_group.healthcheck_port }}</HealthCheckPort>
        {% endif %}
        {% if target_group.healthcheck_protocol %}
        <HealthCheckProtocol>{{ target_group.healthcheck_protocol or "None" }}</HealthCheckProtocol>
        {% endif %}
        {% endif %}
        {% if target_group.healthcheck_path %}
        <HealthCheckPath>{{ target_group.healthcheck_path or '' }}</HealthCheckPath>
        {% endif %}
        <HealthCheckIntervalSeconds>{{ target_group.healthcheck_interval_seconds }}</HealthCheckIntervalSeconds>
        <HealthCheckTimeoutSeconds>{{ target_group.healthcheck_timeout_seconds }}</HealthCheckTimeoutSeconds>
        <HealthyThresholdCount>{{ target_group.healthy_threshold_count }}</HealthyThresholdCount>
        <UnhealthyThresholdCount>{{ target_group.unhealthy_threshold_count }}</UnhealthyThresholdCount>
        <HealthCheckEnabled>{{ target_group.healthcheck_enabled and 'true' or 'false' }}</HealthCheckEnabled>
        {% if target_group.matcher %}
        <Matcher>
          {% if target_group.matcher.get("HttpCode") %}<HttpCode>{{ target_group.matcher['HttpCode'] }}</HttpCode>{% endif %}
          {% if target_group.matcher.get("GrpcCode") %}<GrpcCode>{{ target_group.matcher['GrpcCode'] }}</GrpcCode>{% endif %}
        </Matcher>
        {% endif %}
        {% if target_group.target_type %}
        <TargetType>{{ target_group.target_type }}</TargetType>
        {% endif %}
        {% if target_group.ip_address_type %}
        <IpAddressType>{{ target_group.ip_address_type }}</IpAddressType>
        {% endif %}
      </member>
    </TargetGroups>
  </CreateTargetGroupResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateTargetGroupResponse>"""

CREATE_LISTENER_TEMPLATE = """<CreateListenerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <CreateListenerResult>
    <Listeners>
      <member>
        <LoadBalancerArn>{{ listener.load_balancer_arn }}</LoadBalancerArn>
        <Protocol>{{ listener.protocol }}</Protocol>
        {% if listener.certificates %}
        <Certificates>
          {% for cert in listener.certificates %}
          <member>
            <CertificateArn>{{ cert }}</CertificateArn>
          </member>
          {% endfor %}
        </Certificates>
        {% endif %}
        {% if listener.port %}
        <Port>{{ listener.port }}</Port>
        {% endif %}
        <SslPolicy>{{ listener.ssl_policy }}</SslPolicy>
        <ListenerArn>{{ listener.arn }}</ListenerArn>
        <DefaultActions>
          {% for action in listener.default_actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </DefaultActions>
        <AlpnPolicy>
          {% for policy in listener.alpn_policy %}
          <member>{{ policy }}</member>
          {% endfor %}
        </AlpnPolicy>
      </member>
    </Listeners>
  </CreateListenerResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateListenerResponse>"""

DELETE_LOAD_BALANCER_TEMPLATE = """<DeleteLoadBalancerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DeleteLoadBalancerResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeleteLoadBalancerResponse>"""

DELETE_RULE_TEMPLATE = """<DeleteRuleResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DeleteRuleResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeleteRuleResponse>"""

DELETE_TARGET_GROUP_TEMPLATE = """<DeleteTargetGroupResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DeleteTargetGroupResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeleteTargetGroupResponse>"""

DELETE_LISTENER_TEMPLATE = """<DeleteListenerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DeleteListenerResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeleteListenerResponse>"""

DESCRIBE_LOAD_BALANCERS_TEMPLATE = """<DescribeLoadBalancersResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeLoadBalancersResult>
    <LoadBalancers>
      {% for load_balancer in load_balancers %}
      <member>
        <LoadBalancerArn>{{ load_balancer.arn }}</LoadBalancerArn>
        <Scheme>{{ load_balancer.scheme }}</Scheme>
        <LoadBalancerName>{{ load_balancer.name }}</LoadBalancerName>
        <VpcId>{{ load_balancer.vpc_id }}</VpcId>
        <CanonicalHostedZoneId>Z2P70J7EXAMPLE</CanonicalHostedZoneId>
        <CreatedTime>{{ load_balancer.created_time }}</CreatedTime>
        <AvailabilityZones>
          {% for subnet in load_balancer.subnets %}
          <member>
            <SubnetId>{{ subnet.id }}</SubnetId>
            <ZoneName>{{ subnet.availability_zone }}</ZoneName>
          </member>
          {% endfor %}
        </AvailabilityZones>
        <SecurityGroups>
          {% for security_group in load_balancer.security_groups %}
          <member>{{ security_group }}</member>
          {% endfor %}
        </SecurityGroups>
        <DNSName>{{ load_balancer.dns_name }}</DNSName>
        <State>
          <Code>{{ load_balancer.state }}</Code>
        </State>
        <Type>{{ load_balancer.loadbalancer_type }}</Type>
        <IpAddressType>ipv4</IpAddressType>
      </member>
      {% endfor %}
    </LoadBalancers>
    {% if marker %}
    <NextMarker>{{ marker }}</NextMarker>
    {% endif %}
  </DescribeLoadBalancersResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeLoadBalancersResponse>"""

DESCRIBE_RULES_TEMPLATE = """<DescribeRulesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeRulesResult>
    <Rules>
      {% for rule in rules %}
      <member>
        <IsDefault>{{ "true" if rule.is_default else "false" }}</IsDefault>
        <Conditions>
          {% for condition in rule.conditions %}
          <member>
            <Field>{{ condition["Field"] }}</Field>
            {% if "HttpHeaderConfig" in condition %}
            <HttpHeaderConfig>
              <HttpHeaderName>{{ condition["HttpHeaderConfig"]["HttpHeaderName"] }}</HttpHeaderName>
              <Values>
                {% for value in condition["HttpHeaderConfig"]["Values"] %}
                  <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpHeaderConfig>
            {% endif %}
            {% if "HttpRequestMethodConfig" in condition %}
            <HttpRequestMethodConfig>
              <Values>
                {% for value in condition["HttpRequestMethodConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpRequestMethodConfig>
            {% endif %}
            {% if "QueryStringConfig" in condition %}
            <QueryStringConfig>
              <Values>
                {% for value in condition["QueryStringConfig"]["Values"] %}
                <member>
                    <Key>{{ value["Key"] }}</Key>
                    <Value>{{ value["Value"] }}</Value>
                </member>
                {% endfor %}
              </Values>
            </QueryStringConfig>
            {% endif %}
            {% if "SourceIpConfig" in condition %}
            <SourceIpConfig>
              <Values>
                {% for value in condition["SourceIpConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </SourceIpConfig>
            {% endif %}
            {% if "PathPatternConfig" in condition %}
            <PathPatternConfig>
              <Values>
                {% for value in condition["PathPatternConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </PathPatternConfig>
            {% endif %}
            {% if "HostHeaderConfig" in condition %}
            <HostHeaderConfig>
              <Values>
                {% for value in condition["HostHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HostHeaderConfig>
            {% endif %}
            {% if "Values" in condition %}
            <Values>
              {% for value in condition["Values"] %}
              <member>{{ value }}</member>
              {% endfor %}
            </Values>
            {% endif %}
          </member>
          {% endfor %}
        </Conditions>
        <Priority>{{ rule.priority }}</Priority>
        <RuleArn>{{ rule.arn }}</RuleArn>
        <Actions>
          {% for action in rule.actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </Actions>
      </member>
      {% endfor %}
    </Rules>
    {% if marker %}
    <NextMarker>{{ marker }}</NextMarker>
    {% endif %}
  </DescribeRulesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeRulesResponse>"""

DESCRIBE_TARGET_GROUPS_TEMPLATE = """<DescribeTargetGroupsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeTargetGroupsResult>
    <TargetGroups>
      {% for target_group in target_groups %}
      <member>
        <TargetGroupArn>{{ target_group.arn }}</TargetGroupArn>
        <TargetGroupName>{{ target_group.name }}</TargetGroupName>
        {% if target_group.protocol %}
        <Protocol>{{ target_group.protocol }}</Protocol>
        <ProtocolVersion>{{ target_group.protocol_version }}</ProtocolVersion>
        {% endif %}
        {% if target_group.port %}
        <Port>{{ target_group.port }}</Port>
        {% endif %}
        {% if target_group.vpc_id %}
        <VpcId>{{ target_group.vpc_id }}</VpcId>
        {% endif %}
        <IpAddressType>{{ target_group.ip_address_type }}</IpAddressType>
        <HealthCheckProtocol>{{ target_group.healthcheck_protocol }}</HealthCheckProtocol>
        {% if target_group.healthcheck_port %}<HealthCheckPort>{{ target_group.healthcheck_port }}</HealthCheckPort>{% endif %}
        <HealthCheckPath>{{ target_group.healthcheck_path or '' }}</HealthCheckPath>
        <HealthCheckIntervalSeconds>{{ target_group.healthcheck_interval_seconds }}</HealthCheckIntervalSeconds>
        <HealthCheckTimeoutSeconds>{{ target_group.healthcheck_timeout_seconds }}</HealthCheckTimeoutSeconds>
        <HealthCheckEnabled>{{ target_group.healthcheck_enabled and 'true' or 'false' }}</HealthCheckEnabled>
        <HealthyThresholdCount>{{ target_group.healthy_threshold_count }}</HealthyThresholdCount>
        <UnhealthyThresholdCount>{{ target_group.unhealthy_threshold_count }}</UnhealthyThresholdCount>
        {% if target_group.matcher %}
        <Matcher>
            {% if target_group.matcher.get("HttpCode") %}<HttpCode>{{ target_group.matcher['HttpCode'] }}</HttpCode>{% endif %}
            {% if target_group.matcher.get("GrpcCode") %}<GrpcCode>{{ target_group.matcher['GrpcCode'] }}</GrpcCode>{% endif %}
        </Matcher>
        {% endif %}
        {% if target_group.target_type %}
        <TargetType>{{ target_group.target_type }}</TargetType>
        {% endif %}
        <LoadBalancerArns>
          {% for load_balancer_arn in target_group.load_balancer_arns %}
          <member>{{ load_balancer_arn }}</member>
          {% endfor %}
        </LoadBalancerArns>
      </member>
      {% endfor %}
    </TargetGroups>
  </DescribeTargetGroupsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeTargetGroupsResponse>"""

DESCRIBE_TARGET_GROUP_ATTRIBUTES_TEMPLATE = """<DescribeTargetGroupAttributesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeTargetGroupAttributesResult>
    <Attributes>
      {% for key, value in attributes.items() %}
      <member>
        <Key>{{ key }}</Key>
        <Value>{{ value }}</Value>
      </member>
      {% endfor %}
    </Attributes>
  </DescribeTargetGroupAttributesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeTargetGroupAttributesResponse>"""

DESCRIBE_LISTENERS_TEMPLATE = """<DescribeListenersResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeListenersResult>
    <Listeners>
      {% for listener in listeners %}
      <member>
        <LoadBalancerArn>{{ listener.load_balancer_arn }}</LoadBalancerArn>
        <Protocol>{{ listener.protocol }}</Protocol>
        {% if listener.certificate %}
        <Certificates>
          <member>
            <CertificateArn>{{ listener.certificate }}</CertificateArn>
          </member>
        </Certificates>
        {% endif %}
        {% if listener.port %}<Port>{{ listener.port }}</Port>{% endif %}
        <SslPolicy>{{ listener.ssl_policy }}</SslPolicy>
        <ListenerArn>{{ listener.arn }}</ListenerArn>
        <DefaultActions>
          {% for action in listener.default_actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </DefaultActions>
        <AlpnPolicy>
          {% for policy in listener.alpn_policy %}
          <member>{{ policy }}</member>
          {% endfor %}
        </AlpnPolicy>
      </member>
      {% endfor %}
    </Listeners>
  </DescribeListenersResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeListenersResponse>"""

CONFIGURE_HEALTH_CHECK_TEMPLATE = """<ConfigureHealthCheckResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ConfigureHealthCheckResult>
    <HealthCheck>
      <Interval>{{ check.interval }}</Interval>
      <Target>{{ check.target }}</Target>
      <HealthyThreshold>{{ check.healthy_threshold }}</HealthyThreshold>
      <Timeout>{{ check.timeout }}</Timeout>
      <UnhealthyThreshold>{{ check.unhealthy_threshold }}</UnhealthyThreshold>
    </HealthCheck>
  </ConfigureHealthCheckResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ConfigureHealthCheckResponse>"""

MODIFY_RULE_TEMPLATE = """<ModifyRuleResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ModifyRuleResult>
    <Rules>
      <member>
        <IsDefault>{{ "true" if rules.is_default else "false" }}</IsDefault>
        <Conditions>
          {% for condition in rules.conditions %}
          <member>
            <Field>{{ condition["Field"] }}</Field>
            {% if "PathPatternConfig" in condition %}
            <PathPatternConfig>
              <Values>
                {% for value in condition["PathPatternConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </PathPatternConfig>
            {% endif %}
            {% if "HostHeaderConfig" in condition %}
            <HostHeaderConfig>
              <Values>
                {% for value in condition["HostHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HostHeaderConfig>
            {% endif %}
            {% if "HttpHeaderConfig" in condition %}
            <HttpHeaderConfig>
              <HttpHeaderName>{{ condition["HttpHeaderConfig"]["HttpHeaderName"] }}</HttpHeaderName>
              <Values>
                {% for value in condition["HttpHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpHeaderConfig>
            {% endif %}
            {% if "HttpRequestMethodConfig" in condition %}
            <HttpRequestMethodConfig>
              <Values>
                {% for value in condition["HttpRequestMethodConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpRequestMethodConfig>
            {% endif %}
            {% if "QueryStringConfig" in condition %}
            <QueryStringConfig>
              <Values>
                {% for value in condition["QueryStringConfig"]["Values"] %}
                <member>
                    <Key>{{ value["Key"] }}</Key>
                    <Value>{{ value["Value"] }}</Value>
                </member>
                {% endfor %}
              </Values>
            </QueryStringConfig>
            {% endif %}
            {% if "SourceIpConfig" in condition %}
            <SourceIpConfig>
              <Values>
                {% for value in condition["SourceIpConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </SourceIpConfig>
            {% endif %}
            {% if "Values" in condition %}
            <Values>
              {% for value in condition["Values"] %}
              <member>{{ value }}</member>
              {% endfor %}
            </Values>
            {% endif %}
          </member>
          {% endfor %}
        </Conditions>
        <Priority>{{ rules.priority }}</Priority>
        <RuleArn>{{ rules.arn }}</RuleArn>
        <Actions>
          {% for action in rules.actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </Actions>
      </member>
    </Rules>
  </ModifyRuleResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ModifyRuleResponse>"""

MODIFY_TARGET_GROUP_ATTRIBUTES_TEMPLATE = """<ModifyTargetGroupAttributesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ModifyTargetGroupAttributesResult>
    <Attributes>
      {% for key, value in attributes.items() %}
      <member>
        <Key>{{ key }}</Key>
        <Value>{{ value }}</Value>
      </member>
      {% endfor %}
    </Attributes>
  </ModifyTargetGroupAttributesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ModifyTargetGroupAttributesResponse>"""

REGISTER_TARGETS_TEMPLATE = """<RegisterTargetsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <RegisterTargetsResult>
  </RegisterTargetsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</RegisterTargetsResponse>"""

DEREGISTER_TARGETS_TEMPLATE = """<DeregisterTargetsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DeregisterTargetsResult>
  </DeregisterTargetsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeregisterTargetsResponse>"""

SET_LOAD_BALANCER_SSL_CERTIFICATE = """<SetLoadBalancerListenerSSLCertificateResponse xmlns="http://elasticloadbalan cing.amazonaws.com/doc/2015-12-01/">
 <SetLoadBalancerListenerSSLCertificateResult/>
<ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
</ResponseMetadata>
</SetLoadBalancerListenerSSLCertificateResponse>"""

DELETE_LOAD_BALANCER_LISTENERS = """<DeleteLoadBalancerListenersResponse xmlns="http://elasticloadbalan cing.amazonaws.com/doc/2015-12-01/">
 <DeleteLoadBalancerListenersResult/>
<ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
</ResponseMetadata>
</DeleteLoadBalancerListenersResponse>"""

DESCRIBE_ATTRIBUTES_TEMPLATE = """<DescribeLoadBalancerAttributesResponse  xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeLoadBalancerAttributesResult>
    <LoadBalancerAttributes>
      <AccessLog>
        <Enabled>{{ attributes.access_log.enabled }}</Enabled>
        {% if attributes.access_log.enabled %}
        <S3BucketName>{{ attributes.access_log.s3_bucket_name }}</S3BucketName>
        <S3BucketPrefix>{{ attributes.access_log.s3_bucket_prefix }}</S3BucketPrefix>
        <EmitInterval>{{ attributes.access_log.emit_interval }}</EmitInterval>
        {% endif %}
      </AccessLog>
      <ConnectionSettings>
        <IdleTimeout>{{ attributes.connecting_settings.idle_timeout }}</IdleTimeout>
      </ConnectionSettings>
      <CrossZoneLoadBalancing>
        <Enabled>{{ attributes.cross_zone_load_balancing.enabled }}</Enabled>
      </CrossZoneLoadBalancing>
      <ConnectionDraining>
        {% if attributes.connection_draining.enabled %}
        <Enabled>true</Enabled>
        <Timeout>{{ attributes.connection_draining.timeout }}</Timeout>
        {% else %}
        <Enabled>false</Enabled>
        {% endif %}
      </ConnectionDraining>
    </LoadBalancerAttributes>
  </DescribeLoadBalancerAttributesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeLoadBalancerAttributesResponse>
"""

CREATE_LOAD_BALANCER_POLICY_TEMPLATE = """<CreateLoadBalancerPolicyResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <CreateLoadBalancerPolicyResult/>
  <ResponseMetadata>
      <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateLoadBalancerPolicyResponse>
"""

SET_LOAD_BALANCER_POLICIES_OF_LISTENER_TEMPLATE = """<SetLoadBalancerPoliciesOfListenerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
    <SetLoadBalancerPoliciesOfListenerResult/>
    <ResponseMetadata>
        <RequestId>{{ request_id }}</RequestId>
    </ResponseMetadata>
</SetLoadBalancerPoliciesOfListenerResponse>
"""

SET_LOAD_BALANCER_POLICIES_FOR_BACKEND_SERVER_TEMPLATE = """<SetLoadBalancerPoliciesForBackendServerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
    <SetLoadBalancerPoliciesForBackendServerResult/>
    <ResponseMetadata>
        <RequestId>{{ request_id }}</RequestId>
    </ResponseMetadata>
</SetLoadBalancerPoliciesForBackendServerResponse>
"""

DESCRIBE_TARGET_HEALTH_TEMPLATE = """<DescribeTargetHealthResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeTargetHealthResult>
    <TargetHealthDescriptions>
      {% for target_health in target_health_descriptions %}
      <member>
        {% if target_health.health_port %}<HealthCheckPort>{{ target_health.health_port }}</HealthCheckPort>{% endif %}
        <TargetHealth>
          <State>{{ target_health.status }}</State>
          {% if target_health.reason %}
            <Reason>{{ target_health.reason }}</Reason>
          {% endif %}
          {% if target_health.description %}
            <Description>{{ target_health.description }}</Description>
          {% endif %}
        </TargetHealth>
        <Target>
          {% if target_health.port %}
          <Port>{{ target_health.port }}</Port>
          {% endif %}
          <Id>{{ target_health.instance_id }}</Id>
        </Target>
      </member>
      {% endfor %}
    </TargetHealthDescriptions>
  </DescribeTargetHealthResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeTargetHealthResponse>"""

SET_RULE_PRIORITIES_TEMPLATE = """<SetRulePrioritiesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <SetRulePrioritiesResult>
    <Rules>
      {% for rule in rules %}
      <member>
        <IsDefault>{{ "true" if rule.is_default else "false" }}</IsDefault>
        <Conditions>
          {% for condition in rule.conditions %}
          <member>
            <Field>{{ condition["Field"] }}</Field>
            {% if "Values" in condition %}
            <Values>
              {% for value in condition["Values"] %}
              <member>{{ value }}</member>
              {% endfor %}
            </Values>
            {% endif %}
            {% if "HttpHeaderConfig" in condition %}
            <HttpHeaderConfig>
              <HttpHeaderName>{{ condition["HttpHeaderConfig"]["HttpHeaderName"] }}</HttpHeaderName>
              <Values>
                {% for value in condition["HttpHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpHeaderConfig>
            {% endif %}
            {% if "HttpRequestMethodConfig" in condition %}
            <HttpRequestMethodConfig>
              <Values>
                {% for value in condition["HttpRequestMethodConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HttpRequestMethodConfig>
            {% endif %}
            {% if "QueryStringConfig" in condition %}
            <QueryStringConfig>
              <Values>
                {% for value in condition["QueryStringConfig"]["Values"] %}
                <member>
                    <Key>{{ value["Key"] }}</Key>
                    <Value>{{ value["Value"] }}</Value>
                </member>
                {% endfor %}
              </Values>
            </QueryStringConfig>
            {% endif %}
            {% if "SourceIpConfig" in condition %}
            <SourceIpConfig>
              <Values>
                {% for value in condition["SourceIpConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </SourceIpConfig>
            {% endif %}
            {% if "PathPatternConfig" in condition %}
            <PathPatternConfig>
              <Values>
                {% for value in condition["PathPatternConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </PathPatternConfig>
            {% endif %}
            {% if "HostHeaderConfig" in condition %}
            <HostHeaderConfig>
              <Values>
                {% for value in condition["HostHeaderConfig"]["Values"] %}
                <member>{{ value }}</member>
                {% endfor %}
              </Values>
            </HostHeaderConfig>
            {% endif %}
          </member>
          {% endfor %}
        </Conditions>
        <Priority>{{ rule.priority }}</Priority>
        <RuleArn>{{ rule.arn }}</RuleArn>
        <Actions>
          {% for action in rule.actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </Actions>
      </member>
      {% endfor %}
    </Rules>
  </SetRulePrioritiesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</SetRulePrioritiesResponse>"""

DESCRIBE_LIMITS_TEMPLATE = """<DescribeAccountLimitsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeAccountLimitsResult>
    <Limits>
      {% for key, value in limits.items() %}
      <member>
        <Name>{{ key }}</Name>
        <Max>{{ value }}</Max>
      </member>
      {% endfor %}
    </Limits>
  </DescribeAccountLimitsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeAccountLimitsResponse>"""

DESCRIBE_SSL_POLICIES_TEMPLATE = """<DescribeSSLPoliciesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeSSLPoliciesResult>
    <SslPolicies>
      {% for policy in policies %}
      <member>
        <Name>{{ policy['name'] }}</Name>
        <Ciphers>
          {% for cipher in policy['ciphers'] %}
          <member>
            <Name>{{ cipher['name'] }}</Name>
            <Priority>{{ cipher['priority'] }}</Priority>
          </member>
          {% endfor %}
        </Ciphers>
        <SslProtocols>
          {% for proto in policy['ssl_protocols'] %}
          <member>{{ proto }}</member>
          {% endfor %}
        </SslProtocols>
      </member>
      {% endfor %}
    </SslPolicies>
  </DescribeSSLPoliciesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeSSLPoliciesResponse>"""

SET_IP_ADDRESS_TYPE_TEMPLATE = """<SetIpAddressTypeResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <SetIpAddressTypeResult>
    <IpAddressType>{{ ip_type }}</IpAddressType>
  </SetIpAddressTypeResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</SetIpAddressTypeResponse>"""

SET_SECURITY_GROUPS_TEMPLATE = """<SetSecurityGroupsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <SetSecurityGroupsResult>
    <SecurityGroupIds>
      {% for group in sec_groups %}
      <member>{{ group }}</member>
      {% endfor %}
    </SecurityGroupIds>
  </SetSecurityGroupsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</SetSecurityGroupsResponse>"""

SET_SUBNETS_TEMPLATE = """<SetSubnetsResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <SetSubnetsResult>
    <AvailabilityZones>
      {% for zone_id, subnet_id in subnets.items() %}
      <member>
        <SubnetId>{{ subnet_id }}</SubnetId>
        <ZoneName>{{ zone_id }}</ZoneName>
      </member>
      {% endfor %}
    </AvailabilityZones>
  </SetSubnetsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</SetSubnetsResponse>"""

MODIFY_LOADBALANCER_ATTRS_TEMPLATE = """<ModifyLoadBalancerAttributesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ModifyLoadBalancerAttributesResult>
    <Attributes>
      {% for key, value in attrs.items() %}
      <member>
        {% if value == None %}<Value />{% else %}<Value>{{ value }}</Value>{% endif %}
        <Key>{{ key }}</Key>
      </member>
      {% endfor %}
    </Attributes>
  </ModifyLoadBalancerAttributesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ModifyLoadBalancerAttributesResponse>"""

DESCRIBE_LOADBALANCER_ATTRS_TEMPLATE = """<DescribeLoadBalancerAttributesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeLoadBalancerAttributesResult>
    <Attributes>
      {% for key, value in attrs.items() %}
      <member>
        {% if value == None %}<Value />{% else %}<Value>{{ value }}</Value>{% endif %}
        <Key>{{ key }}</Key>
      </member>
      {% endfor %}
    </Attributes>
  </DescribeLoadBalancerAttributesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeLoadBalancerAttributesResponse>"""


MODIFY_TARGET_GROUP_TEMPLATE = """<ModifyTargetGroupResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ModifyTargetGroupResult>
    <TargetGroups>
      <member>
        <TargetGroupArn>{{ target_group.arn }}</TargetGroupArn>
        <TargetGroupName>{{ target_group.name }}</TargetGroupName>
        <Protocol>{{ target_group.protocol }}</Protocol>
        {% if target_group.port %}<Port>{{ target_group.port }}</Port>{% endif %}
        <VpcId>{{ target_group.vpc_id }}</VpcId>
        <HealthCheckProtocol>{{ target_group.healthcheck_protocol }}</HealthCheckProtocol>
        {% if target_group.healthcheck_port %}<HealthCheckPort>{{ target_group.healthcheck_port }}</HealthCheckPort>{% endif %}
        <HealthCheckPath>{{ target_group.healthcheck_path }}</HealthCheckPath>
        <HealthCheckIntervalSeconds>{{ target_group.healthcheck_interval_seconds }}</HealthCheckIntervalSeconds>
        <HealthCheckTimeoutSeconds>{{ target_group.healthcheck_timeout_seconds }}</HealthCheckTimeoutSeconds>
        <HealthyThresholdCount>{{ target_group.healthy_threshold_count }}</HealthyThresholdCount>
        <UnhealthyThresholdCount>{{ target_group.unhealthy_threshold_count }}</UnhealthyThresholdCount>
        {% if target_group.protocol in ["HTTP", "HTTPS"] %}
        <Matcher>
          <HttpCode>{{ target_group.matcher['HttpCode'] }}</HttpCode>
        </Matcher>
        {% endif %}
        <LoadBalancerArns>
          {% for load_balancer_arn in target_group.load_balancer_arns %}
          <member>{{ load_balancer_arn }}</member>
          {% endfor %}
        </LoadBalancerArns>
      </member>
    </TargetGroups>
  </ModifyTargetGroupResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ModifyTargetGroupResponse>"""

MODIFY_LISTENER_TEMPLATE = """<ModifyListenerResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <ModifyListenerResult>
    <Listeners>
      <member>
        <LoadBalancerArn>{{ listener.load_balancer_arn }}</LoadBalancerArn>
        <Protocol>{{ listener.protocol }}</Protocol>
        {% if listener.certificates %}
        <Certificates>
          {% for cert in listener.certificates %}
          <member>
            <CertificateArn>{{ cert["certificate_arn"] }}</CertificateArn>
          </member>
          {% endfor %}
        </Certificates>
        {% endif %}
        <Port>{{ listener.port }}</Port>
        <SslPolicy>{{ listener.ssl_policy }}</SslPolicy>
        <ListenerArn>{{ listener.arn }}</ListenerArn>
        <DefaultActions>
          {% for action in listener.default_actions %}
          <member>
            {{ action.to_xml() }}
          </member>
          {% endfor %}
        </DefaultActions>
      </member>
    </Listeners>
  </ModifyListenerResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ModifyListenerResponse>"""

ADD_LISTENER_CERTIFICATES_TEMPLATE = """<AddListenerCertificatesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <AddListenerCertificatesResult>
    <Certificates>
      {% for cert in certificates %}
      <member>
        <CertificateArn>{{ cert }}</CertificateArn>
      </member>
      {% endfor %}
    </Certificates>
  </AddListenerCertificatesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</AddListenerCertificatesResponse>"""

DESCRIBE_LISTENER_CERTIFICATES_TEMPLATE = """<DescribeListenerCertificatesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <DescribeListenerCertificatesResult>
    <Certificates>
      {% for cert in certificates %}
      <member>
        <CertificateArn>{{ cert }}</CertificateArn>
      </member>
      {% endfor %}
    </Certificates>
  </DescribeListenerCertificatesResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DescribeListenerCertificatesResponse>"""

REMOVE_LISTENER_CERTIFICATES_TEMPLATE = """<RemoveListenerCertificatesResponse xmlns="http://elasticloadbalancing.amazonaws.com/doc/2015-12-01/">
  <RemoveListenerCertificatesResult />
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</RemoveListenerCertificatesResponse>"""
