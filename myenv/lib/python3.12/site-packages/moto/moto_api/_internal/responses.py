import json
from typing import Any, Dict, List

from botocore.awsrequest import AWSPreparedRequest

from moto import settings
from moto.core import DEFAULT_ACCOUNT_ID
from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import ActionAuthenticatorMixin, BaseResponse


class MotoAPIResponse(BaseResponse):
    def reset_response(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,  # pylint: disable=unused-argument
    ) -> TYPE_RESPONSE:
        if request.method == "POST":
            from .models import moto_api_backend

            moto_api_backend.reset()
            return 200, {}, json.dumps({"status": "ok"})
        return 400, {}, json.dumps({"Error": "Need to POST to reset Moto"})

    def reset_auth_response(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,  # pylint: disable=unused-argument
    ) -> TYPE_RESPONSE:
        if request.method == "POST":
            previous_initial_no_auth_action_count = (
                settings.INITIAL_NO_AUTH_ACTION_COUNT
            )
            settings.INITIAL_NO_AUTH_ACTION_COUNT = float(request.data.decode())
            ActionAuthenticatorMixin.request_count = 0
            return (
                200,
                {},
                json.dumps(
                    {
                        "status": "ok",
                        "PREVIOUS_INITIAL_NO_AUTH_ACTION_COUNT": str(
                            previous_initial_no_auth_action_count
                        ),
                    }
                ),
            )
        return 400, {}, json.dumps({"Error": "Need to POST to reset Moto Auth"})

    def model_data(
        self,
        request: Any,  # pylint: disable=unused-argument
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,  # pylint: disable=unused-argument
    ) -> TYPE_RESPONSE:
        from moto.core.model_instances import model_data

        results: Dict[str, Dict[str, List[Any]]] = {}
        for service in sorted(model_data):
            models = model_data[service]
            results[service] = {}
            for name in sorted(models):
                model = models[name]
                results[service][name] = []
                for instance in model.instances:  # type: ignore[attr-defined]
                    inst_result = {}
                    for attr in dir(instance):
                        if not attr.startswith("_"):
                            try:
                                json.dumps(getattr(instance, attr))
                            except (TypeError, AttributeError, ValueError):
                                pass
                            else:
                                inst_result[attr] = getattr(instance, attr)
                    results[service][name].append(inst_result)
        return 200, {"Content-Type": "application/javascript"}, json.dumps(results)

    def dashboard(
        self,
        request: Any,  # pylint: disable=unused-argument
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,  # pylint: disable=unused-argument
    ) -> str:
        from flask import render_template

        return render_template("dashboard.html")

    def get_transition(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,  # pylint: disable=unused-argument
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        qs_dict = dict(
            x.split("=") for x in request.query_string.decode("utf-8").split("&")
        )
        model_name = qs_dict["model_name"]

        resp = moto_api_backend.get_transition(model_name=model_name)

        return 200, {}, json.dumps(resp)

    def set_transition(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        model_name = body["model_name"]
        transition = body["transition"]

        moto_api_backend.set_transition(model_name, transition)
        return 201, {}, ""

    def unset_transition(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        model_name = body["model_name"]

        moto_api_backend.unset_transition(model_name)
        return 201, {}, ""

    def seed(self, req: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        self.setup_class(req, full_url, headers)
        from . import mock_random

        a = self._get_param("a")
        mock_random.seed(int(a))
        return 200, {}, ""

    def set_athena_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            rows = result["rows"]
            column_info = result.get("column_info", [])
            moto_api_backend.set_athena_result(
                rows=rows,
                column_info=column_info,
                account_id=account_id,
                region=region,
            )
        return 201, {}, ""

    def set_ce_cost_usage_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)

        for result in body.get("results", []):
            moto_api_backend.set_ce_cost_usage(result=result, account_id=account_id)
        return 201, {}, ""

    def set_lambda_simple_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            moto_api_backend.set_lambda_simple_result(
                result=result, account_id=account_id, region=region
            )
        return 201, {}, ""

    def set_resilience_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            moto_api_backend.set_resilience_result(
                result=result, account_id=account_id, region=region
            )
        return 201, {}, ""

    def set_sagemaker_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            body = result["Body"]
            content_type = result.get("ContentType")
            prod_variant = result.get("InvokedProductionVariant")
            custom_attrs = result.get("CustomAttributes")
            moto_api_backend.set_sagemaker_result(
                body=body,
                content_type=content_type,
                prod_variant=prod_variant,
                custom_attrs=custom_attrs,
                account_id=account_id,
                region=region,
            )
        return 201, {}, ""

    def set_rds_data_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            records = result.get("records")
            column_metadata = result.get("columnMetadata")
            nr_of_records_updated = result.get("numberOfRecordsUpdated")
            generated_fields = result.get("generatedFields")
            formatted_records = result.get("formattedRecords")
            moto_api_backend.set_rds_data_result(
                records=records,
                column_metadata=column_metadata,
                nr_of_records_updated=nr_of_records_updated,
                generated_fields=generated_fields,
                formatted_records=formatted_records,
                account_id=account_id,
                region=region,
            )
        return 201, {}, ""

    def set_inspector2_findings_result(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        body = self._get_body(headers, request)
        account_id = body.get("account_id", DEFAULT_ACCOUNT_ID)
        region = body.get("region", "us-east-1")

        for result in body.get("results", []):
            moto_api_backend.set_inspector2_findings_result(
                results=result,
                account_id=account_id,
                region=region,
            )
        return 201, {}, ""

    def set_proxy_passthrough(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        res_headers = {"Content-Type": "application/json"}

        if request.method == "POST":
            body = self._get_body(headers, request)
            http_urls = body.get("http_urls", [])
            https_hosts = body.get("https_hosts", [])
            moto_api_backend.set_proxy_passthrough(http_urls, https_hosts)
        if request.method == "DELETE":
            moto_api_backend.delete_proxy_passthroughs()

        urls, hosts = moto_api_backend.get_proxy_passthrough()
        resp = {"http_urls": list(urls), "https_hosts": list(hosts)}
        return 201, res_headers, json.dumps(resp).encode("utf-8")

    def config(
        self,
        request: Any,
        full_url: str,  # pylint: disable=unused-argument
        headers: Any,
    ) -> TYPE_RESPONSE:
        from .models import moto_api_backend

        res_headers = {"Content-Type": "application/json"}

        if request.method == "POST":
            config = self._get_body(headers, request)
            moto_api_backend.set_config(config)

        config = moto_api_backend.get_config()
        return 201, res_headers, json.dumps(config).encode("utf-8")

    def _get_body(self, headers: Any, request: Any) -> Any:
        if isinstance(request, AWSPreparedRequest):
            return json.loads(request.body)  # type: ignore[arg-type]
        else:
            # Werkzeug request
            request_body_size = int(headers["Content-Length"])
            body = request.environ["wsgi.input"].read(request_body_size).decode("utf-8")
            return json.loads(body)
