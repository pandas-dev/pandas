# Sourced from https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-lambda-function-code-cfnresponsemodule.html
# 01/Nov/2021

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

from __future__ import print_function

import json
from typing import Any

import urllib3

SUCCESS = "SUCCESS"
FAILED = "FAILED"

http = urllib3.PoolManager()


def send(
    event: Any,
    context: Any,
    responseStatus: Any,
    responseData: Any,
    physicalResourceId: Any = None,
    noEcho: bool = False,
    reason: Any = None,
) -> None:
    responseUrl = event["ResponseURL"]

    print(responseUrl)  # noqa: T201

    responseBody = {
        "Status": responseStatus,
        "Reason": reason
        or "See the details in CloudWatch Log Stream: {}".format(
            context.log_stream_name
        ),
        "PhysicalResourceId": physicalResourceId or context.log_stream_name,
        "StackId": event["StackId"],
        "RequestId": event["RequestId"],
        "LogicalResourceId": event["LogicalResourceId"],
        "NoEcho": noEcho,
        "Data": responseData,
    }

    json_responseBody = json.dumps(responseBody)

    print("Response body:")  # noqa: T201
    print(json_responseBody)  # noqa: T201

    headers = {"content-type": "", "content-length": str(len(json_responseBody))}

    try:
        response = http.request(  # type: ignore
            "PUT", responseUrl, headers=headers, body=json_responseBody
        )
        print("Status code:", response.status)  # noqa: T201

    except Exception as e:
        print("send(..) failed executing http.request(..):", e)  # noqa: T201
