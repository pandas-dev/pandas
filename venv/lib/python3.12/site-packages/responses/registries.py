import copy
from typing import TYPE_CHECKING
from typing import List
from typing import Optional
from typing import Tuple

if TYPE_CHECKING:  # pragma: no cover
    # import only for linter run
    from requests import PreparedRequest

    from responses import BaseResponse


class FirstMatchRegistry:
    def __init__(self) -> None:
        self._responses: List["BaseResponse"] = []

    @property
    def registered(self) -> List["BaseResponse"]:
        return self._responses

    def reset(self) -> None:
        self._responses = []

    def find(
        self, request: "PreparedRequest"
    ) -> Tuple[Optional["BaseResponse"], List[str]]:
        found = None
        found_match = None
        match_failed_reasons = []
        for i, response in enumerate(self.registered):
            match_result, reason = response.matches(request)
            if match_result:
                if found is None:
                    found = i
                    found_match = response
                else:
                    if self.registered[found].call_count > 0:
                        # that assumes that some responses were added between calls
                        self.registered.pop(found)
                        found_match = response
                        break
                    # Multiple matches found.  Remove & return the first response.
                    return self.registered.pop(found), match_failed_reasons
            else:
                match_failed_reasons.append(reason)
        return found_match, match_failed_reasons

    def add(self, response: "BaseResponse") -> "BaseResponse":
        if any(response is resp for resp in self.registered):
            # if user adds multiple responses that reference the same instance.
            # do a comparison by memory allocation address.
            # see https://github.com/getsentry/responses/issues/479
            response = copy.deepcopy(response)

        self.registered.append(response)
        return response

    def remove(self, response: "BaseResponse") -> List["BaseResponse"]:
        removed_responses = []
        while response in self.registered:
            self.registered.remove(response)
            removed_responses.append(response)
        return removed_responses

    def replace(self, response: "BaseResponse") -> "BaseResponse":
        try:
            index = self.registered.index(response)
        except ValueError:
            raise ValueError(f"Response is not registered for URL {response.url}")
        self.registered[index] = response
        return response


class OrderedRegistry(FirstMatchRegistry):
    """Registry where `Response` objects are dependent on the insertion order and invocation index.

    OrderedRegistry applies the rule of first in - first out. Responses should be invoked in
    the same order in which they were added to the registry. Otherwise, an error is returned.
    """

    def find(
        self, request: "PreparedRequest"
    ) -> Tuple[Optional["BaseResponse"], List[str]]:
        """Find the next registered `Response` and check if it matches the request.

        Search is performed by taking the first element of the registered responses list
        and removing this object (popping from the list).

        Parameters
        ----------
        request : PreparedRequest
            Request that was caught by the custom adapter.

        Returns
        -------
        Tuple[Optional["BaseResponse"], List[str]]
            Matched `Response` object and empty list in case of match.
            Otherwise, None and a list with reasons for not finding a match.

        """

        if not self.registered:
            return None, ["No more registered responses"]

        response = self.registered.pop(0)
        match_result, reason = response.matches(request)
        if not match_result:
            self.reset()
            self.add(response)
            reason = (
                "Next 'Response' in the order doesn't match "
                f"due to the following reason: {reason}."
            )
            return None, [reason]

        return response, []
