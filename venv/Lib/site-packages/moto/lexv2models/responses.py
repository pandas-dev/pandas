"""Handles incoming lexv2models requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import LexModelsV2Backend, lexv2models_backends


class LexModelsV2Response(BaseResponse):
    """Handler for LexModelsV2 requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="lexv2-models")

    @property
    def lexv2models_backend(self) -> LexModelsV2Backend:
        """Return backend instance specific for this region."""
        # lexv2models_backends is not yet typed
        # Please modify moto/backends.py to add the appropriate type annotations for this service
        return lexv2models_backends[self.current_account][self.region]

    def create_bot(self) -> str:
        bot_name = self._get_param("botName")
        description = self._get_param("description")
        role_arn = self._get_param("roleArn")
        data_privacy = self._get_param("dataPrivacy")
        idle_session_ttl_in_seconds = self._get_param("idleSessionTTLInSeconds")
        bot_tags = self._get_param("botTags")
        test_bot_alias_tags = self._get_param("testBotAliasTags")
        bot_type = self._get_param("botType")
        bot_members = self._get_param("botMembers")

        resp = self.lexv2models_backend.create_bot(
            bot_name=bot_name,
            description=description,
            role_arn=role_arn,
            data_privacy=data_privacy,
            idle_session_ttl_in_seconds=idle_session_ttl_in_seconds,
            bot_tags=bot_tags,
            test_bot_alias_tags=test_bot_alias_tags,
            bot_type=bot_type,
            bot_members=bot_members,
        )

        return json.dumps(resp)

    def describe_bot(self) -> str:
        bot_id = self._get_param("botId")

        return json.dumps(
            self.lexv2models_backend.describe_bot(
                bot_id=bot_id,
            )
        )

    def update_bot(self) -> str:
        bot_id = self._get_param("botId")
        bot_name = self._get_param("botName")
        description = self._get_param("description")
        role_arn = self._get_param("roleArn")
        data_privacy = self._get_param("dataPrivacy")
        idle_session_ttl_in_seconds = self._get_param("idleSessionTTLInSeconds")
        bot_type = self._get_param("botType")
        bot_members = self._get_param("botMembers")

        return json.dumps(
            self.lexv2models_backend.update_bot(
                bot_id=bot_id,
                bot_name=bot_name,
                description=description,
                role_arn=role_arn,
                data_privacy=data_privacy,
                idle_session_ttl_in_seconds=idle_session_ttl_in_seconds,
                bot_type=bot_type,
                bot_members=bot_members,
            )
        )

    def list_bots(self) -> str:
        bot_summaries = self.lexv2models_backend.list_bots()
        return json.dumps({"botSummaries": bot_summaries, "nextToken": None})

    def delete_bot(self) -> str:
        bot_id = self._get_param("botId")
        skip_resource_in_use_check = self._get_param("skipResourceInUseCheck")
        bot_id, bot_status = self.lexv2models_backend.delete_bot(
            bot_id=bot_id,
            skip_resource_in_use_check=skip_resource_in_use_check,
        )
        return json.dumps({"botId": bot_id, "botStatus": bot_status})

    def create_bot_alias(self) -> str:
        bot_alias_name = self._get_param("botAliasName")
        description = self._get_param("description")
        bot_version = self._get_param("botVersion")
        bot_alias_locale_settings = self._get_param("botAliasLocaleSettings")
        conversation_log_settings = self._get_param("conversationLogSettings")
        sentiment_analysis_settings = self._get_param("sentimentAnalysisSettings")
        bot_id = self._get_param("botId")
        tags = self._get_param("tags")

        return json.dumps(
            self.lexv2models_backend.create_bot_alias(
                bot_alias_name=bot_alias_name,
                description=description,
                bot_version=bot_version,
                bot_alias_locale_settings=bot_alias_locale_settings,
                conversation_log_settings=conversation_log_settings,
                sentiment_analysis_settings=sentiment_analysis_settings,
                bot_id=bot_id,
                tags=tags,
            )
        )

    def describe_bot_alias(self) -> str:
        bot_alias_id = self._get_param("botAliasId")
        bot_id = self._get_param("botId")

        return json.dumps(
            self.lexv2models_backend.describe_bot_alias(
                bot_alias_id=bot_alias_id,
                bot_id=bot_id,
            )
        )

    def update_bot_alias(self) -> str:
        bot_alias_id = self._get_param("botAliasId")
        bot_alias_name = self._get_param("botAliasName")
        description = self._get_param("description")
        bot_version = self._get_param("botVersion")
        bot_alias_locale_settings = self._get_param("botAliasLocaleSettings")
        conversation_log_settings = self._get_param("conversationLogSettings")
        sentiment_analysis_settings = self._get_param("sentimentAnalysisSettings")
        bot_id = self._get_param("botId")

        return json.dumps(
            self.lexv2models_backend.update_bot_alias(
                bot_alias_id=bot_alias_id,
                bot_alias_name=bot_alias_name,
                description=description,
                bot_version=bot_version,
                bot_alias_locale_settings=bot_alias_locale_settings,
                conversation_log_settings=conversation_log_settings,
                sentiment_analysis_settings=sentiment_analysis_settings,
                bot_id=bot_id,
            )
        )

    def list_bot_aliases(self) -> str:
        bot_id = self._get_param("botId")
        max_results = self._get_param("maxResults")
        bot_alias_summaries, bot_id = self.lexv2models_backend.list_bot_aliases(
            bot_id=bot_id, max_results=max_results
        )
        return json.dumps(
            {
                "botAliasSummaries": bot_alias_summaries,
                "nextToken": None,
                "botId": bot_id,
            }
        )

    def delete_bot_alias(self) -> str:
        bot_alias_id = self._get_param("botAliasId")
        bot_id = self._get_param("botId")
        skip_resource_in_use_check = self._get_param("skipResourceInUseCheck")
        bot_alias_id, bot_id, bot_alias_status = (
            self.lexv2models_backend.delete_bot_alias(
                bot_alias_id=bot_alias_id,
                bot_id=bot_id,
                skip_resource_in_use_check=skip_resource_in_use_check,
            )
        )
        return json.dumps(
            {
                "botAliasId": bot_alias_id,
                "botId": bot_id,
                "botAliasStatus": bot_alias_status,
            }
        )

    def create_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        policy = self._get_param("policy")
        resource_arn, revision_id = self.lexv2models_backend.create_resource_policy(
            resource_arn=resource_arn,
            policy=policy,
        )
        return json.dumps({"resourceArn": resource_arn, "revisionId": revision_id})

    def describe_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        resource_arn, policy, revision_id = (
            self.lexv2models_backend.describe_resource_policy(
                resource_arn=resource_arn,
            )
        )
        return json.dumps(
            {"resourceArn": resource_arn, "policy": policy, "revisionId": revision_id}
        )

    def update_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        policy = self._get_param("policy")
        expected_revision_id = self._get_param("expectedRevisionId")
        resource_arn, revision_id = self.lexv2models_backend.update_resource_policy(
            resource_arn=resource_arn,
            policy=policy,
            expected_revision_id=expected_revision_id,
        )
        return json.dumps({"resourceArn": resource_arn, "revisionId": revision_id})

    def delete_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        expected_revision_id = self._get_param("expectedRevisionId")
        resource_arn, revision_id = self.lexv2models_backend.delete_resource_policy(
            resource_arn=resource_arn,
            expected_revision_id=expected_revision_id,
        )
        return json.dumps({"resourceArn": resource_arn, "revisionId": revision_id})

    def tag_resource(self) -> str:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = self._get_param("tags")
        self.lexv2models_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps({})

    def untag_resource(self) -> str:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tag_keys = self.__dict__["data"]["tagKeys"]
        self.lexv2models_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps({})

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = self.lexv2models_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps({"tags": tags})
