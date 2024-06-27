from ._base_response import EC2BaseResponse


class IamInstanceProfiles(EC2BaseResponse):
    def associate_iam_instance_profile(self) -> str:
        instance_id = self._get_param("InstanceId")
        iam_instance_profile_name = self._get_param("IamInstanceProfile.Name")
        iam_instance_profile_arn = self._get_param("IamInstanceProfile.Arn")
        iam_association = self.ec2_backend.associate_iam_instance_profile(
            instance_id, iam_instance_profile_name, iam_instance_profile_arn
        )
        template = self.response_template(IAM_INSTANCE_PROFILE_RESPONSE)
        return template.render(iam_association=iam_association, state="associating")

    def describe_iam_instance_profile_associations(self) -> str:
        association_ids = self._get_multi_param("AssociationId")
        filters = self._get_object_map("Filter")
        max_items = self._get_param("MaxItems")
        next_token = self._get_param("NextToken")
        (
            iam_associations,
            next_token,
        ) = self.ec2_backend.describe_iam_instance_profile_associations(
            association_ids, filters, max_items, next_token
        )
        template = self.response_template(DESCRIBE_IAM_INSTANCE_PROFILE_RESPONSE)
        return template.render(iam_associations=iam_associations, next_token=next_token)

    def disassociate_iam_instance_profile(self) -> str:
        association_id = self._get_param("AssociationId")
        iam_association = self.ec2_backend.disassociate_iam_instance_profile(
            association_id
        )
        template = self.response_template(IAM_INSTANCE_PROFILE_RESPONSE)
        return template.render(iam_association=iam_association, state="disassociating")

    def replace_iam_instance_profile_association(self) -> str:
        association_id = self._get_param("AssociationId")
        iam_instance_profile_name = self._get_param("IamInstanceProfile.Name")
        iam_instance_profile_arn = self._get_param("IamInstanceProfile.Arn")
        iam_association = self.ec2_backend.replace_iam_instance_profile_association(
            association_id, iam_instance_profile_name, iam_instance_profile_arn
        )
        template = self.response_template(IAM_INSTANCE_PROFILE_RESPONSE)
        return template.render(iam_association=iam_association, state="associating")


# https://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_AssociateIamInstanceProfile.html
IAM_INSTANCE_PROFILE_RESPONSE = """
<AssociateIamInstanceProfileResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>e10deeaf-7cda-48e7-950b-example</requestId>
    <iamInstanceProfileAssociation>
        <associationId>{{ iam_association.id }}</associationId>
        {% if iam_association.iam_instance_profile %}
        <iamInstanceProfile>
            <arn>{{ iam_association.iam_instance_profile.arn }}</arn>
            <id>{{ iam_association.iam_instance_profile.id }}</id>
        </iamInstanceProfile>
        {% endif %}
        <instanceId>{{ iam_association.instance.id }}</instanceId>
        <state>{{ state }}</state>
    </iamInstanceProfileAssociation>
</AssociateIamInstanceProfileResponse>
"""


# https://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_DescribeIamInstanceProfileAssociations.html
# Note: this API description page contains an error! Provided `iamInstanceProfileAssociations` doesn't work, you
# should use `iamInstanceProfileAssociationSet` instead.
DESCRIBE_IAM_INSTANCE_PROFILE_RESPONSE = """
<DescribeIamInstanceProfileAssociationsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>84c2d2a6-12dc-491f-a9ee-example</requestId>
    {% if next_token %}<nextToken>{{ next_token }}</nextToken>{% endif %}
    <iamInstanceProfileAssociationSet>
         {% for iam_association in iam_associations %}
            <item>
                <associationId>{{ iam_association.id }}</associationId>
                <iamInstanceProfile>
                    <arn>{{ iam_association.iam_instance_profile.arn }}</arn>
                    <id>{{ iam_association.iam_instance_profile.id }}</id>
                </iamInstanceProfile>
                <instanceId>{{ iam_association.instance.id }}</instanceId>
                <state>{{ iam_association.state }}</state>
            </item>
        {% endfor %}
    </iamInstanceProfileAssociationSet>
</DescribeIamInstanceProfileAssociationsResponse>
"""
