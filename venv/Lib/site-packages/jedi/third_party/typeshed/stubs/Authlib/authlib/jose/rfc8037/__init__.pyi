from .jws_eddsa import register_jws_rfc8037 as register_jws_rfc8037
from .okp_key import OKPKey as OKPKey

__all__ = ["register_jws_rfc8037", "OKPKey"]
