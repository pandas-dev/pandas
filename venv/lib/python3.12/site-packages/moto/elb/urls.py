from typing import Any

from moto.core.request import normalize_request
from moto.elb.responses import ELBResponse
from moto.elbv2.responses import ELBV2Response


def api_version_elb_backend(*args: Any, **kwargs: Any) -> Any:
    """
    ELB and ELBV2 (Classic and Application load balancers) use the same
    hostname and url space. To differentiate them we must read the
    `Version` parameter out of the url-encoded request body. TODO: There
    has _got_ to be a better way to do this. Please help us think of
    one.
    """
    request = normalize_request(args[0])
    version = request.values.get("Version")
    if "2012-06-01" == version:
        return ELBResponse.dispatch(*args, **kwargs)
    elif "2015-12-01" == version:
        return ELBV2Response.dispatch(*args, **kwargs)
    else:
        raise Exception(f"Unknown ELB API version: {version}")


url_bases = [r"https?://elasticloadbalancing\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": api_version_elb_backend}
