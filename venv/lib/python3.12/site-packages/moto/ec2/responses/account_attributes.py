from moto.core.responses import BaseResponse


class AccountAttributes(BaseResponse):
    def describe_account_attributes(self) -> str:
        template = self.response_template(DESCRIBE_ACCOUNT_ATTRIBUTES_RESULT)
        return template.render()


DESCRIBE_ACCOUNT_ATTRIBUTES_RESULT = """
<DescribeAccountAttributesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <accountAttributeSet>
    <item>
      <attributeName>vpc-max-security-groups-per-interface</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>5</attributeValue>
        </item>
      </attributeValueSet>
    </item>
    <item>
      <attributeName>max-instances</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>20</attributeValue>
        </item>
      </attributeValueSet>
    </item>
    <item>
      <attributeName>supported-platforms</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>EC2</attributeValue>
        </item>
        <item>
          <attributeValue>VPC</attributeValue>
        </item>
      </attributeValueSet>
    </item>
    <item>
      <attributeName>default-vpc</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>none</attributeValue>
        </item>
      </attributeValueSet>
    </item>
    <item>
      <attributeName>max-elastic-ips</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>5</attributeValue>
        </item>
      </attributeValueSet>
    </item>
    <item>
      <attributeName>vpc-max-elastic-ips</attributeName>
      <attributeValueSet>
        <item>
          <attributeValue>5</attributeValue>
        </item>
      </attributeValueSet>
    </item>
  </accountAttributeSet>
</DescribeAccountAttributesResponse>
"""
