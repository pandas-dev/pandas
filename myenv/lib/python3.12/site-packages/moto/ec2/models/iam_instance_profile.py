from typing import Any, Dict, List, Optional, Tuple

from moto.core.common_models import CloudFormationModel
from moto.ec2.models.instances import Instance
from moto.iam.models import InstanceProfile

from ..exceptions import (
    IncorrectStateIamProfileAssociationError,
    InvalidAssociationIDIamProfileAssociationError,
)
from ..utils import (
    filter_iam_instance_profile_associations,
    filter_iam_instance_profiles,
    random_iam_instance_profile_association_id,
)


class IamInstanceProfileAssociation(CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        association_id: str,
        instance: Instance,
        iam_instance_profile: InstanceProfile,
    ):
        self.ec2_backend = ec2_backend
        self.id = association_id
        self.instance = instance
        self.iam_instance_profile = iam_instance_profile
        self.state = "associated"
        ec2_backend.modify_instance_attribute(
            instance.id,
            "iam_instance_profile",
            {"Arn": self.iam_instance_profile.arn, "Id": association_id},
        )


class IamInstanceProfileAssociationBackend:
    def __init__(self) -> None:
        self.iam_instance_profile_associations: Dict[
            str, IamInstanceProfileAssociation
        ] = {}

    def associate_iam_instance_profile(
        self,
        instance_id: str,
        iam_instance_profile_name: Optional[str] = None,
        iam_instance_profile_arn: Optional[str] = None,
    ) -> IamInstanceProfileAssociation:
        iam_association_id = random_iam_instance_profile_association_id()

        instance_profile = filter_iam_instance_profiles(
            self.account_id,  # type: ignore[attr-defined]
            partition=self.partition,  # type: ignore[attr-defined]
            iam_instance_profile_arn=iam_instance_profile_arn,
            iam_instance_profile_name=iam_instance_profile_name,
        )

        if instance_id in self.iam_instance_profile_associations.keys():
            raise IncorrectStateIamProfileAssociationError(instance_id)

        iam_instance_profile_association = IamInstanceProfileAssociation(
            self,
            iam_association_id,
            self.get_instance(instance_id) if instance_id else None,  # type: ignore[attr-defined]
            instance_profile,
        )
        # Regarding to AWS there can be only one association with ec2.
        self.iam_instance_profile_associations[instance_id] = (
            iam_instance_profile_association
        )
        return iam_instance_profile_association

    def describe_iam_instance_profile_associations(
        self,
        association_ids: List[str],
        filters: Any = None,
        max_results: int = 100,
        next_token: Optional[str] = None,
    ) -> Tuple[List[IamInstanceProfileAssociation], Optional[str]]:
        associations_list: List[IamInstanceProfileAssociation] = []
        if association_ids:
            for association in self.iam_instance_profile_associations.values():
                if association.id in association_ids:
                    associations_list.append(association)
        else:
            # That's mean that no association id were given. Showing all.
            associations_list.extend(self.iam_instance_profile_associations.values())

        associations_list = filter_iam_instance_profile_associations(
            associations_list, filters
        )

        starting_point = int(next_token or 0)
        ending_point = starting_point + int(max_results or 100)
        associations_page = associations_list[starting_point:ending_point]
        new_next_token = (
            str(ending_point) if ending_point < len(associations_list) else None
        )

        return associations_page, new_next_token

    def disassociate_iam_instance_profile(
        self, association_id: str
    ) -> IamInstanceProfileAssociation:
        iam_instance_profile_association = None
        for association_key in self.iam_instance_profile_associations.keys():
            if (
                self.iam_instance_profile_associations[association_key].id
                == association_id
            ):
                iam_instance_profile_association = (
                    self.iam_instance_profile_associations[association_key]
                )
                del self.iam_instance_profile_associations[association_key]
                # Deleting once and avoiding `RuntimeError: dictionary changed size during iteration`
                break

        if not iam_instance_profile_association:
            raise InvalidAssociationIDIamProfileAssociationError(association_id)

        return iam_instance_profile_association

    def replace_iam_instance_profile_association(
        self,
        association_id: str,
        iam_instance_profile_name: Optional[str] = None,
        iam_instance_profile_arn: Optional[str] = None,
    ) -> IamInstanceProfileAssociation:
        instance_profile = filter_iam_instance_profiles(
            self.account_id,  # type: ignore[attr-defined]
            partition=self.partition,  # type: ignore[attr-defined]
            iam_instance_profile_arn=iam_instance_profile_arn,
            iam_instance_profile_name=iam_instance_profile_name,
        )

        iam_instance_profile_association = None
        for association_key in self.iam_instance_profile_associations.keys():
            if (
                self.iam_instance_profile_associations[association_key].id
                == association_id
            ):
                self.iam_instance_profile_associations[
                    association_key
                ].iam_instance_profile = instance_profile
                iam_instance_profile_association = (
                    self.iam_instance_profile_associations[association_key]
                )
                break

        if not iam_instance_profile_association:
            raise InvalidAssociationIDIamProfileAssociationError(association_id)

        return iam_instance_profile_association
