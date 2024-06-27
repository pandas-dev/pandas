from ._base_response import EC2BaseResponse


class KeyPairs(EC2BaseResponse):
    def create_key_pair(self) -> str:
        name = self._get_param("KeyName")
        key_type = self._get_param("KeyType")
        tags = self._parse_tag_specification("key-pair").get("key-pair", {})
        self.error_on_dryrun()
        keypair = self.ec2_backend.create_key_pair(name, key_type, tags=tags)
        return self.response_template(CREATE_KEY_PAIR_RESPONSE).render(keypair=keypair)

    def delete_key_pair(self) -> str:
        name = self._get_param("KeyName")
        self.error_on_dryrun()

        self.ec2_backend.delete_key_pair(name)
        return self.response_template(DELETE_KEY_PAIR_RESPONSE).render()

    def describe_key_pairs(self) -> str:
        names = self._get_multi_param("KeyName")
        filters = self._filters_from_querystring()
        keypairs = self.ec2_backend.describe_key_pairs(names, filters)
        template = self.response_template(DESCRIBE_KEY_PAIRS_RESPONSE)
        return template.render(keypairs=keypairs)

    def import_key_pair(self) -> str:
        name = self._get_param("KeyName")
        material = self._get_param("PublicKeyMaterial")
        tags = self._parse_tag_specification("key-pair").get("key-pair", {})
        self.error_on_dryrun()

        keypair = self.ec2_backend.import_key_pair(name, material, tags=tags)
        return self.response_template(IMPORT_KEYPAIR_RESPONSE).render(keypair=keypair)


DESCRIBE_KEY_PAIRS_RESPONSE = """<DescribeKeyPairsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
    <keySet>
    {% for keypair in keypairs %}
      <item>
           <createTime>{{ keypair.created_iso_8601 }}</createTime>
           <keyPairId>{{ keypair.id }}</keyPairId>
           <keyName>{{ keypair.name }}</keyName>
           <keyFingerprint>{{ keypair.fingerprint }}</keyFingerprint>
           {% if keypair.get_tags() %}
             <tagSet>
               {% for tag in keypair.get_tags() %}
                 <item>
                   <key>{{ tag.key }}</key>
                   <value>{{ tag.value }}</value>
                 </item>
               {% endfor %}
             </tagSet>
           {% endif %}
      </item>
    {% endfor %}
    </keySet>
 </DescribeKeyPairsResponse>"""


CREATE_KEY_PAIR_RESPONSE = """<CreateKeyPairResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <keyPairId>{{ keypair.id }}</keyPairId>
   <keyName>{{ keypair.name }}</keyName>
   <keyFingerprint>{{ keypair.fingerprint }}</keyFingerprint>
   <keyMaterial>{{ keypair.material }}</keyMaterial>
   {% if keypair.get_tags() %}
   <tagSet>
     {% for tag in keypair.get_tags() %}
       <item>
         <key>{{ tag.key }}</key>
         <value>{{ tag.value }}</value>
       </item>
     {% endfor %}
   </tagSet>
   {% endif %}
</CreateKeyPairResponse>"""


DELETE_KEY_PAIR_RESPONSE = """<DeleteKeyPairResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DeleteKeyPairResponse>"""

IMPORT_KEYPAIR_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
  <ImportKeyPairResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>471f9fdd-8fe2-4a84-86b0-bd3d3e350979</requestId>
    <keyPairId>{{ keypair.id }}</keyPairId>
    <keyName>{{ keypair.name }}</keyName>
    <keyFingerprint>{{ keypair.fingerprint }}</keyFingerprint>
    {% if keypair.get_tags() %}
   <tagSet>
     {% for tag in keypair.get_tags() %}
       <item>
         <key>{{ tag.key }}</key>
         <value>{{ tag.value }}</value>
       </item>
     {% endfor %}
   </tagSet>
   {% endif %}
  </ImportKeyPairResponse>"""
