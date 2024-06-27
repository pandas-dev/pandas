# This will only exist in responses >= 0.17
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlunparse

import responses

from .custom_responses_mock import CallbackResponse, not_implemented_callback
from .utils import get_equivalent_url_in_aws_domain


class CustomRegistry(responses.registries.FirstMatchRegistry):
    """
    Custom Registry that returns requests in an order that makes sense for Moto:
     - Implemented callbacks take precedence over non-implemented-callbacks
     - CallbackResponses are not discarded after first use - users can mock the same URL as often as they like
     - CallbackResponses are persisted in a dictionary, with the request-method as key
       This reduces the number of possible responses that we need to search
    """

    def __init__(self) -> None:
        self._registered: Dict[str, List[responses.BaseResponse]] = defaultdict(list)

    @property
    def registered(self) -> List[responses.BaseResponse]:
        res = []
        for resps in self._registered.values():
            res += resps
        return res

    def add(self, response: responses.BaseResponse) -> responses.BaseResponse:
        self._registered[response.method].append(response)
        return response

    def reset(self) -> None:
        self._registered.clear()

    def find(self, request: Any) -> Tuple[Optional[responses.BaseResponse], List[str]]:
        # We don't have to search through all possible methods - only the ones registered for this particular method
        all_possibles = (
            responses._default_mock._registry.registered
            + self._registered[request.method]
        )
        found = []
        match_failed_reasons = []

        # Handle non-standard AWS endpoint hostnames from ISO regions or custom S3 endpoints.
        parsed_url, url_was_modified = get_equivalent_url_in_aws_domain(request.url)
        if url_was_modified:
            url_with_standard_aws_domain = urlunparse(parsed_url)
            request_with_standard_aws_domain = request.copy()
            request_with_standard_aws_domain.prepare_url(
                url_with_standard_aws_domain, {}
            )
        else:
            request_with_standard_aws_domain = request

        for response in all_possibles:
            match_result, reason = response.matches(request_with_standard_aws_domain)
            if match_result:
                found.append(response)
            else:
                match_failed_reasons.append(reason)

        # Look for implemented callbacks first
        implemented_matches = [
            m
            for m in found
            if type(m) is not CallbackResponse or m.callback != not_implemented_callback
        ]
        if implemented_matches:
            # Prioritize responses.CallbackResponse
            # Usecases:
            #
            #  - Callbacks created by APIGateway to intercept URL requests to *.execute-api.amazonaws.com
            #
            for match in implemented_matches:
                if type(match) == responses.CallbackResponse:
                    return match, match_failed_reasons
            # Return moto.core.custom_responses_mock.CallbackResponse otherwise
            return implemented_matches[0], match_failed_reasons
        elif found:
            # We had matches, but all were of type not_implemented_callback
            return found[0], match_failed_reasons
        else:
            return None, match_failed_reasons
