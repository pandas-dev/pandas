import json
import re
import string
from typing import Any, Dict, List, Optional
from urllib.parse import parse_qs, urlparse

from moto.moto_api._internal import mock_random as random


def networkid_from_managedblockchain_url(full_url: str) -> str:
    id_search = re.search(r"\/n-[A-Z0-9]{26}", full_url, re.IGNORECASE)
    return_id = None
    if id_search:
        return_id = id_search.group(0).replace("/", "")
    return return_id  # type: ignore[return-value]


def get_network_id() -> str:
    return "n-" + "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
    )


def memberid_from_managedblockchain_request(full_url: str, body: Dict[str, Any]) -> str:
    id_search = re.search(r"\/m-[A-Z0-9]{26}", full_url, re.IGNORECASE)
    return_id = None
    if id_search:
        return_id = id_search.group(0).replace("/", "")
    else:
        # >= botocore 1.19.41 can add the memberId as a query parameter, or in the body
        parsed_url = urlparse(full_url)
        qs = parse_qs(parsed_url.query)
        if "memberId" in qs:
            return_id = qs.get("memberId")[0]  # type: ignore
        elif body:
            body = json.loads(body)  # type: ignore
            return_id = body["MemberId"]
    return return_id  # type: ignore[return-value]


def get_member_id() -> str:
    return "m-" + "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
    )


def proposalid_from_managedblockchain_url(full_url: str) -> str:
    id_search = re.search(r"\/p-[A-Z0-9]{26}", full_url, re.IGNORECASE)
    return_id = None
    if id_search:
        return_id = id_search.group(0).replace("/", "")
    return return_id  # type: ignore[return-value]


def get_proposal_id() -> str:
    return "p-" + "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
    )


def invitationid_from_managedblockchain_url(full_url: str) -> str:
    id_search = re.search(r"\/in-[A-Z0-9]{26}", full_url, re.IGNORECASE)
    return_id = None
    if id_search:
        return_id = id_search.group(0).replace("/", "")
    return return_id  # type: ignore[return-value]


def get_invitation_id() -> str:
    return "in-" + "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
    )


def member_name_exist_in_network(
    members: Dict[str, Any], networkid: str, membername: str
) -> bool:
    for member in members.values():
        if member.network_id == networkid:
            if member.name == membername:
                return True
    return False


def number_of_members_in_network(
    members: Dict[str, Any], networkid: str, member_status: Optional[str] = None
) -> int:
    return len(
        [
            member
            for member in members.values()
            if member.network_id == networkid
            and (member_status is None or member.member_status == member_status)
        ]
    )


def admin_password_ok(password: str) -> bool:
    if not re.search("[a-z]", password):
        return False
    elif not re.search("[A-Z]", password):
        return False
    elif not re.search("[0-9]", password):
        return False
    elif re.search("['\"@\\/]", password):
        return False
    else:
        return True


def nodeid_from_managedblockchain_url(full_url: str) -> str:
    id_search = re.search(r"\/nd-[A-Z0-9]{26}", full_url, re.IGNORECASE)
    return_id = None
    if id_search:
        return_id = id_search.group(0).replace("/", "")
    return return_id  # type: ignore[return-value]


def get_node_id() -> str:
    return "nd-" + "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
    )


def number_of_nodes_in_member(
    nodes: Dict[str, Any], memberid: str, node_status: Optional[str] = None
) -> int:
    return len(
        [
            node
            for node in nodes.values()
            if node.member_id == memberid
            and (node_status is None or node.node_status == node_status)
        ]
    )


def nodes_in_member(nodes: Dict[str, Any], memberid: str) -> List[str]:
    return [nodid for nodid in nodes if nodes[nodid].member_id == memberid]
