from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class IamInstanceProfiles(EC2BaseResponse):
    def associate_iam_instance_profile(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        iam_instance_profile_name = self._get_param("IamInstanceProfile.Name")
        iam_instance_profile_arn = self._get_param("IamInstanceProfile.Arn")
        iam_association = self.ec2_backend.associate_iam_instance_profile(
            instance_id, iam_instance_profile_name, iam_instance_profile_arn
        )
        result = {"IamInstanceProfileAssociation": iam_association}
        return ActionResult(result)

    def describe_iam_instance_profile_associations(self) -> ActionResult:
        association_ids = self._get_multi_param("AssociationId")
        filters = self._filters_from_querystring()
        max_items = self._get_param("MaxItems")
        next_token = self._get_param("NextToken")
        (
            iam_associations,
            next_token,
        ) = self.ec2_backend.describe_iam_instance_profile_associations(
            association_ids, filters, max_items, next_token
        )
        result = {
            "IamInstanceProfileAssociations": iam_associations,
            "NextToken": next_token,
        }
        return ActionResult(result)

    def disassociate_iam_instance_profile(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        iam_association = self.ec2_backend.disassociate_iam_instance_profile(
            association_id
        )
        result = {"IamInstanceProfileAssociation": iam_association}
        return ActionResult(result)

    def replace_iam_instance_profile_association(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        iam_instance_profile_name = self._get_param("IamInstanceProfile.Name")
        iam_instance_profile_arn = self._get_param("IamInstanceProfile.Arn")
        iam_association = self.ec2_backend.replace_iam_instance_profile_association(
            association_id, iam_instance_profile_name, iam_instance_profile_arn
        )
        result = {"IamInstanceProfileAssociation": iam_association}
        return ActionResult(result)
