from moto.core.responses import BaseResponse

from .exceptions import STSValidationError
from .models import STSBackend, sts_backends

MAX_FEDERATION_TOKEN_POLICY_LENGTH = 2048


class TokenResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="sts")

    @property
    def backend(self) -> STSBackend:
        return sts_backends[self.current_account][self.partition]

    def _determine_resource(self) -> str:
        if "AssumeRole" in self.querystring.get("Action", []):
            return self.querystring.get("RoleArn")[0]  # type: ignore[index]
        return "*"

    def get_session_token(self) -> str:
        duration = int(self.querystring.get("DurationSeconds", [43200])[0])
        token = self.backend.get_session_token(duration=duration)
        template = self.response_template(GET_SESSION_TOKEN_RESPONSE)
        return template.render(token=token)

    def get_federation_token(self) -> str:
        duration = int(self.querystring.get("DurationSeconds", [43200])[0])
        policy = self.querystring.get("Policy", [None])[0]

        if policy is not None and len(policy) > MAX_FEDERATION_TOKEN_POLICY_LENGTH:
            raise STSValidationError(
                "1 validation error detected: Value "
                '\'{"Version": "2012-10-17", "Statement": [...]}\' '
                "at 'policy' failed to satisfy constraint: Member must have length less than or "
                f" equal to {MAX_FEDERATION_TOKEN_POLICY_LENGTH}"
            )

        name = self.querystring.get("Name")[0]  # type: ignore
        token = self.backend.get_federation_token(duration=duration, name=name)
        template = self.response_template(GET_FEDERATION_TOKEN_RESPONSE)
        return template.render(
            token=token, account_id=self.current_account, partition=self.partition
        )

    def assume_role(self) -> str:
        role_session_name = self.querystring.get("RoleSessionName")[0]  # type: ignore
        role_arn = self.querystring.get("RoleArn")[0]  # type: ignore

        policy = self.querystring.get("Policy", [None])[0]
        duration = int(self.querystring.get("DurationSeconds", [3600])[0])
        external_id = self.querystring.get("ExternalId", [None])[0]

        role = self.backend.assume_role(
            region_name=self.region,
            role_session_name=role_session_name,
            role_arn=role_arn,
            policy=policy,
            duration=duration,
            external_id=external_id,
        )
        template = self.response_template(ASSUME_ROLE_RESPONSE)
        return template.render(role=role)

    def assume_role_with_web_identity(self) -> str:
        role_session_name = self.querystring.get("RoleSessionName")[0]  # type: ignore
        role_arn = self.querystring.get("RoleArn")[0]  # type: ignore

        policy = self.querystring.get("Policy", [None])[0]
        duration = int(self.querystring.get("DurationSeconds", [3600])[0])
        external_id = self.querystring.get("ExternalId", [None])[0]

        role = self.backend.assume_role_with_web_identity(
            region_name=self.region,
            role_session_name=role_session_name,
            role_arn=role_arn,
            policy=policy,
            duration=duration,
            external_id=external_id,
        )
        template = self.response_template(ASSUME_ROLE_WITH_WEB_IDENTITY_RESPONSE)
        return template.render(role=role)

    def assume_role_with_saml(self) -> str:
        role_arn = self.querystring.get("RoleArn")[0]  # type: ignore
        principal_arn = self.querystring.get("PrincipalArn")[0]  # type: ignore
        saml_assertion = self.querystring.get("SAMLAssertion")[0]  # type: ignore

        role = self.backend.assume_role_with_saml(
            role_arn=role_arn,
            principal_arn=principal_arn,
            saml_assertion=saml_assertion,
        )
        template = self.response_template(ASSUME_ROLE_WITH_SAML_RESPONSE)
        return template.render(role=role)

    def get_caller_identity(self) -> str:
        template = self.response_template(GET_CALLER_IDENTITY_RESPONSE)

        access_key_id = self.get_access_key()
        user_id, arn, account_id = self.backend.get_caller_identity(
            access_key_id, self.region
        )

        return template.render(account_id=account_id, user_id=user_id, arn=arn)


GET_SESSION_TOKEN_RESPONSE = """<GetSessionTokenResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <GetSessionTokenResult>
    <Credentials>
      <SessionToken>AQoEXAMPLEH4aoAH0gNCAPyJxz4BlCFFxWNE1OPTgk5TthT+FvwqnKwRcOIfrRh3c/LTo6UDdyJwOOvEVPvLXCrrrUtdnniCEXAMPLE/IvU1dYUg2RVAJBanLiHb4IgRmpRV3zrkuWJOgQs8IZZaIv2BXIa2R4OlgkBN9bkUDNCJiBeb/AXlzBBko7b15fjrBs2+cTQtpZ3CYWFXG8C5zqx37wnOE49mRl/+OtkIKGO7fAE</SessionToken>
      <SecretAccessKey>wJalrXUtnFEMI/K7MDENG/bPxRfiCYzEXAMPLEKEY</SecretAccessKey>
      <Expiration>{{ token.expiration_ISO8601 }}</Expiration>
      <AccessKeyId>AKIAIOSFODNN7EXAMPLE</AccessKeyId>
    </Credentials>
  </GetSessionTokenResult>
  <ResponseMetadata>
    <RequestId>58c5dbae-abef-11e0-8cfe-09039844ac7d</RequestId>
  </ResponseMetadata>
</GetSessionTokenResponse>"""


GET_FEDERATION_TOKEN_RESPONSE = """<GetFederationTokenResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <GetFederationTokenResult>
    <Credentials>
      <SessionToken>AQoDYXdzEPT//////////wEXAMPLEtc764bNrC9SAPBSM22wDOk4x4HIZ8j4FZTwdQWLWsKWHGBuFqwAeMicRXmxfpSPfIeoIYRqTflfKD8YUuwthAx7mSEI/qkPpKPi/kMcGdQrmGdeehM4IC1NtBmUpp2wUE8phUZampKsburEDy0KPkyQDYwT7WZ0wq5VSXDvp75YU9HFvlRd8Tx6q6fE8YQcHNVXAkiY9q6d+xo0rKwT38xVqr7ZD0u0iPPkUL64lIZbqBAz+scqKmlzm8FDrypNC9Yjc8fPOLn9FX9KSYvKTr4rvx3iSIlTJabIQwj2ICCR/oLxBA==</SessionToken>
      <SecretAccessKey>wJalrXUtnFEMI/K7MDENG/bPxRfiCYzEXAMPLEKEY</SecretAccessKey>
      <Expiration>{{ token.expiration_ISO8601 }}</Expiration>
      <AccessKeyId>AKIAIOSFODNN7EXAMPLE</AccessKeyId>
    </Credentials>
    <FederatedUser>
      <Arn>arn:{{ partition }}:sts::{{ account_id }}:federated-user/{{ token.name }}</Arn>
      <FederatedUserId>{{ account_id }}:{{ token.name }}</FederatedUserId>
    </FederatedUser>
    <PackedPolicySize>6</PackedPolicySize>
  </GetFederationTokenResult>
  <ResponseMetadata>
    <RequestId>c6104cbe-af31-11e0-8154-cbc7ccf896c7</RequestId>
  </ResponseMetadata>
</GetFederationTokenResponse>"""


ASSUME_ROLE_RESPONSE = """<AssumeRoleResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <AssumeRoleResult>
    <Credentials>
      <SessionToken>{{ role.session_token }}</SessionToken>
      <SecretAccessKey>{{ role.secret_access_key }}</SecretAccessKey>
      <Expiration>{{ role.expiration_ISO8601 }}</Expiration>
      <AccessKeyId>{{ role.access_key_id }}</AccessKeyId>
    </Credentials>
    <AssumedRoleUser>
      <Arn>{{ role.arn }}</Arn>
      <AssumedRoleId>{{ role.user_id }}</AssumedRoleId>
    </AssumedRoleUser>
    <PackedPolicySize>6</PackedPolicySize>
  </AssumeRoleResult>
  <ResponseMetadata>
    <RequestId>c6104cbe-af31-11e0-8154-cbc7ccf896c7</RequestId>
  </ResponseMetadata>
</AssumeRoleResponse>"""


ASSUME_ROLE_WITH_WEB_IDENTITY_RESPONSE = """<AssumeRoleWithWebIdentityResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <AssumeRoleWithWebIdentityResult>
    <Credentials>
      <SessionToken>{{ role.session_token }}</SessionToken>
      <SecretAccessKey>{{ role.secret_access_key }}</SecretAccessKey>
      <Expiration>{{ role.expiration_ISO8601 }}</Expiration>
      <AccessKeyId>{{ role.access_key_id }}</AccessKeyId>
    </Credentials>
    <AssumedRoleUser>
      <Arn>{{ role.arn }}</Arn>
      <AssumedRoleId>ARO123EXAMPLE123:{{ role.session_name }}</AssumedRoleId>
    </AssumedRoleUser>
    <PackedPolicySize>6</PackedPolicySize>
  </AssumeRoleWithWebIdentityResult>
  <ResponseMetadata>
    <RequestId>c6104cbe-af31-11e0-8154-cbc7ccf896c7</RequestId>
  </ResponseMetadata>
</AssumeRoleWithWebIdentityResponse>"""


ASSUME_ROLE_WITH_SAML_RESPONSE = """<AssumeRoleWithSAMLResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <AssumeRoleWithSAMLResult>
    <Audience>https://signin.aws.amazon.com/saml</Audience>
    <AssumedRoleUser>
      <AssumedRoleId>{{ role.user_id }}</AssumedRoleId>
      <Arn>{{ role.arn }}</Arn>
    </AssumedRoleUser>
    <Credentials>
      <AccessKeyId>{{ role.access_key_id }}</AccessKeyId>
      <SecretAccessKey>{{ role.secret_access_key }}</SecretAccessKey>
      <SessionToken>{{ role.session_token }}</SessionToken>
      <Expiration>{{ role.expiration_ISO8601 }}</Expiration>
    </Credentials>
    <Subject>{{ role.user_id }}</Subject>
    <NameQualifier>B64EncodedStringOfHashOfIssuerAccountIdAndUserId=</NameQualifier>
    <SubjectType>persistent</SubjectType>
    <Issuer>http://localhost:3000/</Issuer>
  </AssumeRoleWithSAMLResult>
  <ResponseMetadata>
    <RequestId>c6104cbe-af31-11e0-8154-cbc7ccf896c7</RequestId>
  </ResponseMetadata>
</AssumeRoleWithSAMLResponse>"""


GET_CALLER_IDENTITY_RESPONSE = """<GetCallerIdentityResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
  <GetCallerIdentityResult>
    <Arn>{{ arn }}</Arn>
    <UserId>{{ user_id }}</UserId>
    <Account>{{ account_id }}</Account>
  </GetCallerIdentityResult>
  <ResponseMetadata>
    <RequestId>c6104cbe-af31-11e0-8154-cbc7ccf896c7</RequestId>
  </ResponseMetadata>
</GetCallerIdentityResponse>
"""
