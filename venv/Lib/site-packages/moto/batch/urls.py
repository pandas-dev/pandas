from .responses import BatchResponse

url_bases = [r"https?://batch\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/v1/createcomputeenvironment$": BatchResponse.dispatch,
    "{0}/v1/describecomputeenvironments$": BatchResponse.dispatch,
    "{0}/v1/deletecomputeenvironment": BatchResponse.dispatch,
    "{0}/v1/updatecomputeenvironment": BatchResponse.dispatch,
    "{0}/v1/createjobqueue": BatchResponse.dispatch,
    "{0}/v1/describejobqueues": BatchResponse.dispatch,
    "{0}/v1/updatejobqueue": BatchResponse.dispatch,
    "{0}/v1/deletejobqueue": BatchResponse.dispatch,
    "{0}/v1/registerjobdefinition": BatchResponse.dispatch,
    "{0}/v1/deregisterjobdefinition": BatchResponse.dispatch,
    "{0}/v1/describejobdefinitions": BatchResponse.dispatch,
    "{0}/v1/createschedulingpolicy": BatchResponse.dispatch,
    "{0}/v1/describeschedulingpolicies": BatchResponse.dispatch,
    "{0}/v1/listschedulingpolicies": BatchResponse.dispatch,
    "{0}/v1/deleteschedulingpolicy": BatchResponse.dispatch,
    "{0}/v1/updateschedulingpolicy": BatchResponse.dispatch,
    "{0}/v1/submitjob": BatchResponse.dispatch,
    "{0}/v1/describejobs": BatchResponse.dispatch,
    "{0}/v1/listjobs": BatchResponse.dispatch,
    "{0}/v1/terminatejob": BatchResponse.dispatch,
    "{0}/v1/canceljob": BatchResponse.dispatch,
    "{0}/v1/tags/(?P<arn_part_1>[^/]+)/(?P<arn_part_2>[^/]+)/?$": BatchResponse.dispatch,
    "{0}/v1/tags/(?P<arn>[^/]+)/?$": BatchResponse.dispatch,
}
