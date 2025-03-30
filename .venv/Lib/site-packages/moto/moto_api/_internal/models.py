from typing import Any, Dict, List, Optional, Set, Tuple

from moto.core import DEFAULT_ACCOUNT_ID
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.config import default_user_config
from moto.core.model_instances import reset_model_data


class MotoAPIBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.proxy_urls_to_passthrough: Set[str] = set()
        self.proxy_hosts_to_passthrough: Set[str] = set()

    def reset(self) -> None:
        region_name = self.region_name
        account_id = self.account_id

        BackendDict.reset()
        reset_model_data()
        self.__init__(region_name, account_id)  # type: ignore[misc]

    def get_transition(self, model_name: str) -> Dict[str, Any]:
        from moto.moto_api import state_manager

        return state_manager.get_transition(model_name)

    def set_transition(self, model_name: str, transition: Dict[str, Any]) -> None:
        from moto.moto_api import state_manager

        state_manager.set_transition(model_name, transition)

    def unset_transition(self, model_name: str) -> None:
        from moto.moto_api import state_manager

        state_manager.unset_transition(model_name)

    def set_athena_result(
        self,
        rows: List[Dict[str, Any]],
        column_info: List[Dict[str, str]],
        account_id: str,
        region: str,
    ) -> None:
        from moto.athena.models import QueryResults, athena_backends

        backend = athena_backends[account_id][region]
        results = QueryResults(rows=rows, column_info=column_info)
        backend.query_results_queue.append(results)

    def set_ce_cost_usage(self, result: Dict[str, Any], account_id: str) -> None:
        from moto.ce.models import ce_backends

        backend = ce_backends[account_id]["global"]
        backend.cost_usage_results_queue.append(result)

    def set_lambda_simple_result(
        self, result: str, account_id: str, region: str
    ) -> None:
        from moto.awslambda_simple.models import lambda_simple_backends

        backend = lambda_simple_backends[account_id][region]
        backend.lambda_simple_results_queue.append(result)

    def set_resilience_result(
        self, result: List[Dict[str, Any]], account_id: str, region: str
    ) -> None:
        from moto.resiliencehub.models import resiliencehub_backends

        backend = resiliencehub_backends[account_id][region]
        backend.app_assessments_queue.append(result)

    def set_sagemaker_result(
        self,
        body: str,
        content_type: str,
        prod_variant: str,
        custom_attrs: str,
        account_id: str,
        region: str,
    ) -> None:
        from moto.sagemakerruntime.models import sagemakerruntime_backends

        backend = sagemakerruntime_backends[account_id][region]
        backend.results_queue.append((body, content_type, prod_variant, custom_attrs))

    def set_sagemaker_async_result(
        self,
        is_failure: bool,
        data: str,
        account_id: str,
        region: str,
    ) -> None:
        from moto.sagemakerruntime.models import sagemakerruntime_backends

        backend = sagemakerruntime_backends[account_id][region]
        backend.async_results_queue.append((is_failure, data))

    def set_rds_data_result(
        self,
        records: Optional[List[List[Dict[str, Any]]]],
        column_metadata: Optional[List[Dict[str, Any]]],
        nr_of_records_updated: Optional[int],
        generated_fields: Optional[List[Dict[str, Any]]],
        formatted_records: Optional[str],
        account_id: str,
        region: str,
    ) -> None:
        from moto.rdsdata.models import QueryResults, rdsdata_backends

        backend = rdsdata_backends[account_id][region]
        backend.results_queue.append(
            QueryResults(
                records=records,
                column_metadata=column_metadata,
                number_of_records_updated=nr_of_records_updated,
                generated_fields=generated_fields,
                formatted_records=formatted_records,
            )
        )

    def set_inspector2_findings_result(
        self,
        results: Optional[List[List[Dict[str, Any]]]],
        account_id: str,
        region: str,
    ) -> None:
        from moto.inspector2.models import inspector2_backends

        backend = inspector2_backends[account_id][region]
        backend.findings_queue.append(results)

    def set_timestream_result(
        self,
        query: Optional[str],
        query_results: List[Dict[str, Any]],
        account_id: str,
        region: str,
    ) -> None:
        from moto.timestreamquery.models import (
            TimestreamQueryBackend,
            timestreamquery_backends,
        )

        backend: TimestreamQueryBackend = timestreamquery_backends[account_id][region]
        if query not in backend.query_result_queue:
            backend.query_result_queue[query] = []
        backend.query_result_queue[query].extend(query_results)

    def get_proxy_passthrough(self) -> Tuple[Set[str], Set[str]]:
        return self.proxy_urls_to_passthrough, self.proxy_hosts_to_passthrough

    def set_proxy_passthrough(
        self, http_urls: List[str], https_hosts: List[str]
    ) -> None:
        for url in http_urls:
            self.proxy_urls_to_passthrough.add(url)
        for host in https_hosts:
            self.proxy_hosts_to_passthrough.add(host)

    def delete_proxy_passthroughs(self) -> None:
        self.proxy_urls_to_passthrough.clear()
        self.proxy_hosts_to_passthrough.clear()

    def get_config(self) -> Dict[str, Any]:
        return {
            "batch": default_user_config["batch"],
            "lambda": default_user_config["lambda"],
        }

    def set_config(self, config: Dict[str, Any]) -> None:
        if "batch" in config:
            default_user_config["batch"] = config["batch"]
        if "lambda" in config:
            default_user_config["lambda"] = config["lambda"]
        if "stepfunctions" in config:
            default_user_config["stepfunctions"] = config["stepfunctions"]


moto_api_backend = MotoAPIBackend(region_name="global", account_id=DEFAULT_ACCOUNT_ID)
