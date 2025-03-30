from collections import namedtuple

Arn = namedtuple(
    "Arn", ["account", "region", "service", "resource_type", "resource_id"]
)


def parse_arn(arn: str) -> Arn:
    # http://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html
    # this method needs probably some more fine tuning,
    # when also other targets are supported
    _, _, service, region, account, resource = arn.split(":", 5)

    if ":" in resource and "/" in resource:
        if resource.index(":") < resource.index("/"):
            resource_type, resource_id = resource.split(":", 1)
        else:
            resource_type, resource_id = resource.split("/", 1)
    elif ":" in resource:
        resource_type, resource_id = resource.split(":", 1)
    elif "/" in resource:
        resource_type, resource_id = resource.split("/", 1)
    else:
        resource_type = None
        resource_id = resource

    return Arn(
        account=account,
        region=region,
        service=service,
        resource_type=resource_type,
        resource_id=resource_id,
    )
