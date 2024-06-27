from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_eval import (
    ResourceEval,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_eval_s3 import (
    ResourceEvalS3,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    Resource,
    ServiceResource,
)


def resource_eval_for(resource: Resource) -> ResourceEval:
    if isinstance(resource, ServiceResource):
        if resource.service_name == "s3":
            return ResourceEvalS3(resource=resource)
    raise ValueError(
        f"ItemReader's Resource fields must be states service resource, instead got '{resource.resource_arn}'."
    )
