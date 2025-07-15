import datetime
import re
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow

from .exceptions import (
    BadRequestException,
    InvalidRequestException,
    ResourceAlreadyExistsException,
    ResourceLimitExceededException,
    ResourceNotFoundException,
)
from .utils import (
    admin_password_ok,
    get_invitation_id,
    get_member_id,
    get_network_id,
    get_node_id,
    get_proposal_id,
    member_name_exist_in_network,
    nodes_in_member,
    number_of_members_in_network,
    number_of_nodes_in_member,
)

FRAMEWORKS = [
    "HYPERLEDGER_FABRIC",
]

FRAMEWORKVERSIONS = [
    "1.2",
]

EDITIONS: Dict[str, Any] = {
    "STARTER": {
        "MaxMembers": 5,
        "MaxNodesPerMember": 2,
        "AllowedNodeInstanceTypes": ["bc.t3.small", "bc.t3.medium"],
    },
    "STANDARD": {
        "MaxMembers": 14,
        "MaxNodesPerMember": 3,
        "AllowedNodeInstanceTypes": ["bc.t3", "bc.m5", "bc.c5"],
    },
}

VOTEVALUES = ["YES", "NO"]


class ManagedBlockchainNetwork(BaseModel):
    def __init__(
        self,
        network_id: str,
        name: str,
        framework: str,
        frameworkversion: str,
        frameworkconfiguration: Dict[str, Any],
        voting_policy: Dict[str, Any],
        member_configuration: Dict[str, Any],
        region: str,
        description: Optional[str] = None,
    ):
        self.creationdate = utcnow()
        self.id = network_id
        self.name = name
        self.description = description
        self.framework = framework
        self.frameworkversion = frameworkversion
        self.frameworkconfiguration = frameworkconfiguration
        self.voting_policy = voting_policy
        self.member_configuration = member_configuration
        self.region = region

    @property
    def network_name(self) -> str:
        return self.name

    @property
    def network_framework(self) -> str:
        return self.framework

    @property
    def network_framework_version(self) -> str:
        return self.frameworkversion

    @property
    def network_creationdate(self) -> str:
        return self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z")

    @property
    def network_description(self) -> Optional[str]:
        return self.description

    @property
    def network_edition(self) -> str:
        return self.frameworkconfiguration["Fabric"]["Edition"]

    @property
    def vote_pol_proposal_duration(self) -> float:
        return self.voting_policy["ApprovalThresholdPolicy"]["ProposalDurationInHours"]

    @property
    def vote_pol_threshold_percentage(self) -> float:
        return self.voting_policy["ApprovalThresholdPolicy"]["ThresholdPercentage"]

    @property
    def vote_pol_threshold_comparator(self) -> str:
        return self.voting_policy["ApprovalThresholdPolicy"]["ThresholdComparator"]

    def to_dict(self) -> Dict[str, Any]:
        # Format for list_networks
        d = {
            "Id": self.id,
            "Name": self.name,
            "Framework": self.framework,
            "FrameworkVersion": self.frameworkversion,
            "Status": "AVAILABLE",
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
        }
        if self.description is not None:
            d["Description"] = self.description
        return d

    def get_format(self) -> Dict[str, Any]:
        # Format for get_network
        frameworkattributes = {
            "Fabric": {
                "OrderingServiceEndpoint": f"orderer.{self.id.lower()}.managedblockchain.{self.region}.amazonaws.com:30001",
                "Edition": self.frameworkconfiguration["Fabric"]["Edition"],
            }
        }

        vpcendpointname = (
            f"com.amazonaws.{self.region}.managedblockchain.{self.id.lower()}"
        )

        d = {
            "Id": self.id,
            "Name": self.name,
            "Framework": self.framework,
            "FrameworkVersion": self.frameworkversion,
            "FrameworkAttributes": frameworkattributes,
            "VpcEndpointServiceName": vpcendpointname,
            "VotingPolicy": self.voting_policy,
            "Status": "AVAILABLE",
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
        }
        if self.description is not None:
            d["Description"] = self.description
        return d


class ManagedBlockchainProposal(BaseModel):
    def __init__(
        self,
        proposal_id: str,
        networkid: str,
        memberid: str,
        membername: str,
        numofmembers: int,
        actions: Dict[str, Any],
        network_expiration: float,
        network_threshold: float,
        network_threshold_comp: str,
        description: Optional[str] = None,
    ):
        # In general, passing all values instead of creating
        # an apparatus to look them up
        self.id = proposal_id
        self.networkid = networkid
        self.memberid = memberid
        self.membername = membername
        self.numofmembers = numofmembers
        self.actions = actions
        self.network_expiration = network_expiration
        self.network_threshold = network_threshold
        self.network_threshold_comp = network_threshold_comp
        self.description = description

        self.creationdate = utcnow()
        self.expirationdate = self.creationdate + datetime.timedelta(
            hours=network_expiration
        )
        self.yes_vote_count = 0
        self.no_vote_count = 0
        self.outstanding_vote_count = self.numofmembers
        self.status = "IN_PROGRESS"
        self.votes: Dict[str, Dict[str, str]] = {}

    @property
    def network_id(self) -> str:
        return self.networkid

    @property
    def proposal_status(self) -> str:
        return self.status

    @property
    def proposal_votes(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self.votes

    def proposal_actions(self, action_type: str) -> List[Dict[str, Any]]:
        if action_type.lower() == "invitations":
            if "Invitations" in self.actions:
                return self.actions["Invitations"]
        elif action_type.lower() == "removals":
            if "Removals" in self.actions:
                return self.actions["Removals"]
        return []

    def check_to_expire_proposal(self) -> None:
        if utcnow() > self.expirationdate:
            self.status = "EXPIRED"

    def to_dict(self) -> Dict[str, Any]:
        # Format for list_proposals
        return {
            "ProposalId": self.id,
            "ProposedByMemberId": self.memberid,
            "ProposedByMemberName": self.membername,
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "ExpirationDate": self.expirationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
        }

    def get_format(self) -> Dict[str, Any]:
        # Format for get_proposal
        d = {
            "ProposalId": self.id,
            "NetworkId": self.networkid,
            "Actions": self.actions,
            "ProposedByMemberId": self.memberid,
            "ProposedByMemberName": self.membername,
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "ExpirationDate": self.expirationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "YesVoteCount": self.yes_vote_count,
            "NoVoteCount": self.no_vote_count,
            "OutstandingVoteCount": self.outstanding_vote_count,
        }
        if self.description is not None:
            d["Description"] = self.description
        return d

    def set_vote(self, votermemberid: str, votermembername: str, vote: str) -> None:
        if vote.upper() == "YES":
            self.yes_vote_count += 1
        else:
            self.no_vote_count += 1
        self.outstanding_vote_count -= 1

        perct_yes = (self.yes_vote_count / self.numofmembers) * 100
        perct_no = (self.no_vote_count / self.numofmembers) * 100
        self.votes[votermemberid] = {
            "MemberId": votermemberid,
            "MemberName": votermembername,
            "Vote": vote.upper(),
        }

        if self.network_threshold_comp == "GREATER_THAN_OR_EQUAL_TO":
            if perct_yes >= self.network_threshold:
                self.status = "APPROVED"
            elif perct_no >= self.network_threshold:
                self.status = "REJECTED"
        else:
            if perct_yes > self.network_threshold:
                self.status = "APPROVED"
            elif perct_no > self.network_threshold:
                self.status = "REJECTED"

        # It is a tie - reject
        if (
            self.status == "IN_PROGRESS"
            and self.network_threshold_comp == "GREATER_THAN"
            and self.outstanding_vote_count == 0
            and perct_yes == perct_no
        ):
            self.status = "REJECTED"


class ManagedBlockchainInvitation(BaseModel):
    def __init__(
        self,
        invitation_id: str,
        networkid: str,
        networkname: str,
        networkframework: str,
        networkframeworkversion: str,
        networkcreationdate: str,
        region: str,
        networkdescription: Optional[str] = None,
    ):
        self.id = invitation_id
        self.networkid = networkid
        self.networkname = networkname
        self.networkdescription = networkdescription
        self.networkframework = networkframework
        self.networkframeworkversion = networkframeworkversion
        self.networkstatus = "AVAILABLE"
        self.networkcreationdate = networkcreationdate
        self.status = "PENDING"
        self.region = region

        self.creationdate = utcnow()
        self.expirationdate = self.creationdate + datetime.timedelta(days=7)

    @property
    def invitation_status(self) -> str:
        return self.status

    @property
    def invitation_networkid(self) -> str:
        return self.networkid

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "InvitationId": self.id,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "ExpirationDate": self.expirationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "Status": self.status,
            "NetworkSummary": {
                "Id": self.networkid,
                "Name": self.networkname,
                "Framework": self.networkframework,
                "FrameworkVersion": self.networkframeworkversion,
                "Status": self.networkstatus,
                "CreationDate": self.networkcreationdate,
            },
        }
        if self.networkdescription is not None:
            d["NetworkSummary"]["Description"] = self.networkdescription
        return d

    def accept_invitation(self) -> None:
        self.status = "ACCEPTED"

    def reject_invitation(self) -> None:
        self.status = "REJECTED"

    def set_network_status(self, network_status: str) -> None:
        self.networkstatus = network_status


class ManagedBlockchainMember(BaseModel):
    def __init__(
        self,
        member_id: str,
        networkid: str,
        member_configuration: Dict[str, Any],
        region: str,
    ):
        self.creationdate = utcnow()
        self.id = member_id
        self.networkid = networkid
        self.member_configuration = member_configuration
        self.status = "AVAILABLE"
        self.region = region
        self.description = None

    @property
    def network_id(self) -> str:
        return self.networkid

    @property
    def name(self) -> str:
        return self.member_configuration["Name"]

    @property
    def member_status(self) -> str:
        return self.status

    def to_dict(self) -> Dict[str, Any]:
        # Format for list_members
        d = {
            "Id": self.id,
            "Name": self.member_configuration["Name"],
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "IsOwned": True,
        }
        if "Description" in self.member_configuration:
            self.description = self.member_configuration["Description"]
        return d

    def get_format(self) -> Dict[str, Any]:
        # Format for get_member
        frameworkattributes = {
            "Fabric": {
                "AdminUsername": self.member_configuration["FrameworkConfiguration"][
                    "Fabric"
                ]["AdminUsername"],
                "CaEndpoint": f"ca.{self.id.lower()}.{self.networkid.lower()}.managedblockchain.{self.region}.amazonaws.com:30002",
            }
        }

        d = {
            "NetworkId": self.networkid,
            "Id": self.id,
            "Name": self.name,
            "FrameworkAttributes": frameworkattributes,
            "LogPublishingConfiguration": self.member_configuration[
                "LogPublishingConfiguration"
            ],
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
        }
        if "Description" in self.member_configuration:
            d["Description"] = self.description
        return d

    def delete(self) -> None:
        self.status = "DELETED"

    def update(self, logpublishingconfiguration: Dict[str, Any]) -> None:
        self.member_configuration["LogPublishingConfiguration"] = (
            logpublishingconfiguration
        )


class ManagedBlockchainNode(BaseModel):
    def __init__(
        self,
        node_id: str,
        networkid: str,
        memberid: str,
        availabilityzone: str,
        instancetype: str,
        logpublishingconfiguration: Dict[str, Any],
        region: str,
    ):
        self.creationdate = utcnow()
        self.id = node_id
        self.instancetype = instancetype
        self.networkid = networkid
        self.memberid = memberid
        self.logpublishingconfiguration = logpublishingconfiguration
        self.region = region
        self.status = "AVAILABLE"
        self.availabilityzone = availabilityzone

    @property
    def member_id(self) -> str:
        return self.memberid

    @property
    def node_status(self) -> str:
        return self.status

    def to_dict(self) -> Dict[str, Any]:
        # Format for list_nodes
        return {
            "Id": self.id,
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
            "AvailabilityZone": self.availabilityzone,
            "InstanceType": self.instancetype,
        }

    def get_format(self) -> Dict[str, Any]:
        # Format for get_node
        frameworkattributes = {
            "Fabric": {
                "PeerEndpoint": f"{self.id.lower()}.{self.networkid.lower()}.{self.memberid.lower()}.managedblockchain.{self.region}.amazonaws.com:30003",
                "PeerEventEndpoint": f"{self.id.lower()}.{self.networkid.lower()}.{self.memberid.lower()}.managedblockchain.{self.region}.amazonaws.com:30004",
            }
        }

        return {
            "NetworkId": self.networkid,
            "MemberId": self.memberid,
            "Id": self.id,
            "InstanceType": self.instancetype,
            "AvailabilityZone": self.availabilityzone,
            "FrameworkAttributes": frameworkattributes,
            "LogPublishingConfiguration": self.logpublishingconfiguration,
            "Status": self.status,
            "CreationDate": self.creationdate.strftime("%Y-%m-%dT%H:%M:%S.%f%z"),
        }

    def delete(self) -> None:
        self.status = "DELETED"

    def update(self, logpublishingconfiguration: Dict[str, Any]) -> None:
        self.logpublishingconfiguration = logpublishingconfiguration


class ManagedBlockchainBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.networks: Dict[str, ManagedBlockchainNetwork] = {}
        self.members: Dict[str, ManagedBlockchainMember] = {}
        self.proposals: Dict[str, ManagedBlockchainProposal] = {}
        self.invitations: Dict[str, ManagedBlockchainInvitation] = {}
        self.nodes: Dict[str, ManagedBlockchainNode] = {}

    def create_network(
        self,
        name: str,
        framework: str,
        frameworkversion: str,
        frameworkconfiguration: Dict[str, Any],
        voting_policy: Dict[str, Any],
        member_configuration: Dict[str, Any],
        description: Optional[str] = None,
    ) -> Dict[str, str]:
        # Check framework
        if framework not in FRAMEWORKS:
            raise BadRequestException("CreateNetwork", "Invalid request body")

        # Check framework version
        if frameworkversion not in FRAMEWORKVERSIONS:
            raise BadRequestException(
                "CreateNetwork",
                f"Invalid version {frameworkversion} requested for framework HYPERLEDGER_FABRIC",
            )

        # Check edition
        if frameworkconfiguration["Fabric"]["Edition"] not in EDITIONS:
            raise BadRequestException("CreateNetwork", "Invalid request body")

        # Generate network ID
        network_id = get_network_id()

        # Generate memberid ID and initial member
        member_id = get_member_id()
        self.members[member_id] = ManagedBlockchainMember(
            member_id=member_id,
            networkid=network_id,
            member_configuration=member_configuration,
            region=self.region_name,
        )

        self.networks[network_id] = ManagedBlockchainNetwork(
            network_id=network_id,
            name=name,
            framework=framework,
            frameworkversion=frameworkversion,
            frameworkconfiguration=frameworkconfiguration,
            voting_policy=voting_policy,
            member_configuration=member_configuration,
            region=self.region_name,
            description=description,
        )

        # Return the network and member ID
        return {"NetworkId": network_id, "MemberId": member_id}

    def list_networks(self) -> List[ManagedBlockchainNetwork]:
        return list(self.networks.values())

    def get_network(self, network_id: str) -> ManagedBlockchainNetwork:
        if network_id not in self.networks:
            raise ResourceNotFoundException(
                "GetNetwork", f"Network {network_id} not found."
            )
        return self.networks[network_id]

    def create_proposal(
        self,
        networkid: str,
        memberid: str,
        actions: Dict[str, Any],
        description: Optional[str] = None,
    ) -> Dict[str, str]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "CreateProposal", f"Network {networkid} not found."
            )

        # Check if member exists
        if memberid not in self.members:
            raise ResourceNotFoundException(
                "CreateProposal", f"Member {memberid} not found."
            )

        # CLI docs say that Invitations and Removals cannot both be passed - but it does
        # not throw an error and can be performed
        if "Invitations" in actions:
            for propinvitation in actions["Invitations"]:
                if re.match("[0-9]{12}", propinvitation["Principal"]) is None:
                    raise InvalidRequestException(
                        "CreateProposal",
                        "Account ID format specified in proposal is not valid.",
                    )

        if "Removals" in actions:
            for propmember in actions["Removals"]:
                if propmember["MemberId"] not in self.members:
                    raise InvalidRequestException(
                        "CreateProposal",
                        "Member ID format specified in proposal is not valid.",
                    )

        # Generate proposal ID
        proposal_id = get_proposal_id()

        self.proposals[proposal_id] = ManagedBlockchainProposal(
            proposal_id=proposal_id,
            networkid=networkid,
            memberid=memberid,
            membername=self.members[memberid].name,
            numofmembers=number_of_members_in_network(self.members, networkid),
            actions=actions,
            network_expiration=self.networks[networkid].vote_pol_proposal_duration,
            network_threshold=self.networks[networkid].vote_pol_threshold_percentage,
            network_threshold_comp=self.networks[
                networkid
            ].vote_pol_threshold_comparator,
            description=description,
        )

        # Return the proposal ID
        return {"ProposalId": proposal_id}

    def list_proposals(self, networkid: str) -> List[ManagedBlockchainProposal]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "ListProposals", f"Network {networkid} not found."
            )

        proposalsfornetwork = []
        for proposal_id in self.proposals:
            if self.proposals[proposal_id].network_id == networkid:
                # See if any are expired
                self.proposals[proposal_id].check_to_expire_proposal()
                proposalsfornetwork.append(self.proposals[proposal_id])
        return proposalsfornetwork

    def get_proposal(
        self, networkid: str, proposalid: str
    ) -> ManagedBlockchainProposal:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "GetProposal", f"Network {networkid} not found."
            )

        if proposalid not in self.proposals:
            raise ResourceNotFoundException(
                "GetProposal", f"Proposal {proposalid} not found."
            )

        # See if it needs to be set to expipred
        self.proposals[proposalid].check_to_expire_proposal()
        return self.proposals[proposalid]

    def vote_on_proposal(
        self, networkid: str, proposalid: str, votermemberid: str, vote: str
    ) -> None:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "VoteOnProposal", f"Network {networkid} not found."
            )

        if proposalid not in self.proposals:
            raise ResourceNotFoundException(
                "VoteOnProposal", f"Proposal {proposalid} not found."
            )

        if votermemberid not in self.members:
            raise ResourceNotFoundException(
                "VoteOnProposal", f"Member {votermemberid} not found."
            )

        if vote.upper() not in VOTEVALUES:
            raise BadRequestException("VoteOnProposal", "Invalid request body")

        # See if it needs to be set to expipred
        self.proposals[proposalid].check_to_expire_proposal()

        # Exception if EXPIRED
        if self.proposals[proposalid].proposal_status == "EXPIRED":
            raise InvalidRequestException(
                "VoteOnProposal",
                f"Proposal {proposalid} is expired and you cannot vote on it.",
            )

        # Check if IN_PROGRESS
        if self.proposals[proposalid].proposal_status != "IN_PROGRESS":
            raise InvalidRequestException(
                "VoteOnProposal",
                f"Proposal {proposalid} has status {self.proposals[proposalid].proposal_status} and you cannot vote on it.",
            )

        # Check to see if this member already voted
        if votermemberid in self.proposals[proposalid].proposal_votes:
            raise ResourceAlreadyExistsException(
                "VoteOnProposal",
                f"Member {votermemberid} has already voted on proposal {proposalid}.",
            )

        # Cast vote
        self.proposals[proposalid].set_vote(
            votermemberid, self.members[votermemberid].name, vote.upper()
        )

        if self.proposals[proposalid].proposal_status == "APPROVED":
            # Generate invitations
            for _ in self.proposals[proposalid].proposal_actions("Invitations"):
                invitation_id = get_invitation_id()
                self.invitations[invitation_id] = ManagedBlockchainInvitation(
                    invitation_id=invitation_id,
                    networkid=networkid,
                    networkname=self.networks[networkid].network_name,
                    networkframework=self.networks[networkid].network_framework,
                    networkframeworkversion=self.networks[
                        networkid
                    ].network_framework_version,
                    networkcreationdate=self.networks[networkid].network_creationdate,
                    region=self.region_name,
                    networkdescription=self.networks[networkid].network_description,
                )

            # Delete members
            for propmember in self.proposals[proposalid].proposal_actions("Removals"):
                self.delete_member(networkid, propmember["MemberId"])

    def list_proposal_votes(self, networkid: str, proposalid: str) -> List[str]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "ListProposalVotes", f"Network {networkid} not found."
            )

        if proposalid not in self.proposals:
            raise ResourceNotFoundException(
                "ListProposalVotes", f"Proposal {proposalid} not found."
            )

        # Output the vote summaries
        proposalvotesfornetwork = []
        for proposal in self.proposals.values():
            if proposal.network_id == networkid:
                for proposal_vote in proposal.proposal_votes.values():
                    proposalvotesfornetwork.append(proposal_vote)
        return proposalvotesfornetwork

    def list_invitations(self) -> List[ManagedBlockchainInvitation]:
        return list(self.invitations.values())

    def reject_invitation(self, invitationid: str) -> None:
        if invitationid not in self.invitations:
            raise ResourceNotFoundException(
                "RejectInvitation", f"InvitationId {invitationid} not found."
            )
        self.invitations[invitationid].reject_invitation()

    def create_member(
        self, invitationid: str, networkid: str, member_configuration: Dict[str, Any]
    ) -> Dict[str, str]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "CreateMember", f"Network {networkid} not found."
            )

        if invitationid not in self.invitations:
            raise InvalidRequestException(
                "CreateMember", f"Invitation {invitationid} not valid"
            )

        if self.invitations[invitationid].invitation_status != "PENDING":
            raise InvalidRequestException(
                "CreateMember", f"Invitation {invitationid} not valid"
            )

        if member_name_exist_in_network(
            self.members, networkid, member_configuration["Name"]
        ):
            raise InvalidRequestException(
                "CreateMember",
                f"Member name {member_configuration['Name']} already exists in network {networkid}.",
            )

        networkedition = self.networks[networkid].network_edition
        if (
            number_of_members_in_network(self.members, networkid)
            >= EDITIONS[networkedition]["MaxMembers"]
        ):
            raise ResourceLimitExceededException(
                "CreateMember",
                f"You cannot create a member in network {networkid}.{EDITIONS[networkedition]['MaxMembers']} is the maximum number of members allowed in a {networkedition} Edition network.",
            )

        memberadminpassword = member_configuration["FrameworkConfiguration"]["Fabric"][
            "AdminPassword"
        ]
        if admin_password_ok(memberadminpassword) is False:
            raise BadRequestException("CreateMember", "Invalid request body")

        member_id = get_member_id()
        self.members[member_id] = ManagedBlockchainMember(
            member_id=member_id,
            networkid=networkid,
            member_configuration=member_configuration,
            region=self.region_name,
        )

        # Accept the invitaiton
        self.invitations[invitationid].accept_invitation()

        # Return the member ID
        return {"MemberId": member_id}

    def list_members(self, networkid: str) -> List[ManagedBlockchainMember]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "ListMembers", f"Network {networkid} not found."
            )

        membersfornetwork = []
        for member in self.members.values():
            if member.network_id == networkid:
                membersfornetwork.append(member)
        return membersfornetwork

    def get_member(self, networkid: str, memberid: str) -> ManagedBlockchainMember:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "GetMember", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "GetMember", f"Member {memberid} not found."
            )

        # Cannot get a member than has been deleted (it does show up in the list)
        if self.members[memberid].member_status == "DELETED":
            raise ResourceNotFoundException(
                "GetMember", f"Member {memberid} not found."
            )

        return self.members[memberid]

    def delete_member(self, networkid: str, memberid: str) -> None:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "DeleteMember", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "DeleteMember", f"Member {memberid} not found."
            )

        self.members[memberid].delete()

        # Is this the last member in the network? (all set to DELETED)
        if number_of_members_in_network(
            self.members, networkid, member_status="DELETED"
        ) == len(self.members):
            # Set network status to DELETED for all invitations
            for invitation in self.invitations.values():
                if invitation.invitation_networkid == networkid:
                    invitation.set_network_status("DELETED")

            # Remove network
            del self.networks[networkid]

        # Remove any nodes associated
        for nodeid in nodes_in_member(self.nodes, memberid):
            del self.nodes[nodeid]

    def update_member(
        self, networkid: str, memberid: str, logpublishingconfiguration: Dict[str, Any]
    ) -> None:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "UpdateMember", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "UpdateMember", f"Member {memberid} not found."
            )

        self.members[memberid].update(logpublishingconfiguration)

    def create_node(
        self,
        networkid: str,
        memberid: str,
        availabilityzone: str,
        instancetype: str,
        logpublishingconfiguration: Dict[str, Any],
    ) -> Dict[str, str]:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "CreateNode", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "CreateNode", f"Member {memberid} not found."
            )

        networkedition = self.networks[networkid].network_edition
        if (
            number_of_nodes_in_member(self.nodes, memberid)
            >= EDITIONS[networkedition]["MaxNodesPerMember"]
        ):
            raise ResourceLimitExceededException(
                "CreateNode",
                f"Maximum number of nodes exceeded in member {memberid}. The maximum number of nodes you can have in a member in a {networkedition} Edition network is {EDITIONS[networkedition]['MaxNodesPerMember']}",
            )

        # See if the instance family is correct
        correctinstancefamily = False
        for chkinsttypepre in EDITIONS["STANDARD"]["AllowedNodeInstanceTypes"]:
            chkinsttypepreregex = chkinsttypepre + ".*"
            if re.match(chkinsttypepreregex, instancetype, re.IGNORECASE):
                correctinstancefamily = True
                break

        if correctinstancefamily is False:
            raise InvalidRequestException(
                "CreateNode",
                f"Requested instance {instancetype} isn't supported.",
            )

        # Check for specific types for starter
        if networkedition == "STARTER":
            if instancetype not in EDITIONS["STARTER"]["AllowedNodeInstanceTypes"]:
                raise InvalidRequestException(
                    "CreateNode",
                    f"Instance type {instancetype} is not supported with STARTER Edition networks.",
                )

        # Simple availability zone check
        chkregionpreregex = self.region_name + "[a-z]"
        if re.match(chkregionpreregex, availabilityzone, re.IGNORECASE) is None:
            raise InvalidRequestException(
                "CreateNode", "Availability Zone is not valid"
            )

        node_id = get_node_id()
        self.nodes[node_id] = ManagedBlockchainNode(
            node_id=node_id,
            networkid=networkid,
            memberid=memberid,
            availabilityzone=availabilityzone,
            instancetype=instancetype,
            logpublishingconfiguration=logpublishingconfiguration,
            region=self.region_name,
        )

        # Return the node ID
        return {"NodeId": node_id}

    def list_nodes(
        self, networkid: str, memberid: str, status: Optional[str] = None
    ) -> List[ManagedBlockchainNode]:
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "ListNodes", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "ListNodes", f"Member {memberid} not found."
            )

        # If member is deleted, cannot list nodes
        if self.members[memberid].member_status == "DELETED":
            raise ResourceNotFoundException(
                "ListNodes", f"Member {memberid} not found."
            )

        nodesformember = []
        for node in self.nodes.values():
            if node.member_id == memberid and (
                status is None or node.node_status == status
            ):
                nodesformember.append(node)
        return nodesformember

    def get_node(
        self, networkid: str, memberid: str, nodeid: str
    ) -> ManagedBlockchainNode:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "GetNode", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException("GetNode", f"Member {memberid} not found.")

        if nodeid not in self.nodes:
            raise ResourceNotFoundException("GetNode", f"Node {nodeid} not found.")

        # Cannot get a node than has been deleted (it does show up in the list)
        if self.nodes[nodeid].node_status == "DELETED":
            raise ResourceNotFoundException("GetNode", f"Node {nodeid} not found.")

        return self.nodes[nodeid]

    def delete_node(self, networkid: str, memberid: str, nodeid: str) -> None:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "DeleteNode", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "DeleteNode", f"Member {memberid} not found."
            )

        if nodeid not in self.nodes:
            raise ResourceNotFoundException("DeleteNode", f"Node {nodeid} not found.")

        self.nodes[nodeid].delete()

    def update_node(
        self,
        networkid: str,
        memberid: str,
        nodeid: str,
        logpublishingconfiguration: Dict[str, Any],
    ) -> None:
        # Check if network exists
        if networkid not in self.networks:
            raise ResourceNotFoundException(
                "UpdateNode", f"Network {networkid} not found."
            )

        if memberid not in self.members:
            raise ResourceNotFoundException(
                "UpdateNode", f"Member {memberid} not found."
            )

        if nodeid not in self.nodes:
            raise ResourceNotFoundException("UpdateNode", f"Node {nodeid} not found.")

        self.nodes[nodeid].update(logpublishingconfiguration)


managedblockchain_backends = BackendDict(ManagedBlockchainBackend, "managedblockchain")
