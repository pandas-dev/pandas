import json
from typing import Any

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import GlacierBackend, glacier_backends
from .utils import vault_from_glacier_url


class GlacierResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="glacier")

    def setup_class(
        self, request: Any, full_url: str, headers: Any, use_raw_body: bool = False
    ) -> None:
        super().setup_class(request, full_url, headers, use_raw_body=True)

    @property
    def glacier_backend(self) -> GlacierBackend:
        return glacier_backends[self.current_account][self.region]

    def list_vaults(self) -> TYPE_RESPONSE:
        vaults = self.glacier_backend.list_vaults()
        response = json.dumps(
            {"Marker": None, "VaultList": [vault.to_dict() for vault in vaults]}
        )

        headers = {"content-type": "application/json"}
        return 200, headers, response

    def describe_vault(self) -> TYPE_RESPONSE:
        vault_name = vault_from_glacier_url(self.uri)
        vault = self.glacier_backend.describe_vault(vault_name)
        headers = {"content-type": "application/json"}
        return 200, headers, json.dumps(vault.to_dict())

    def create_vault(self) -> TYPE_RESPONSE:
        vault_name = vault_from_glacier_url(self.uri)
        self.glacier_backend.create_vault(vault_name)
        return 201, {"status": 201}, ""

    def delete_vault(self) -> TYPE_RESPONSE:
        vault_name = vault_from_glacier_url(self.uri)
        self.glacier_backend.delete_vault(vault_name)
        return 204, {"status": 204}, ""

    def upload_archive(self) -> TYPE_RESPONSE:
        description = self.headers.get("x-amz-archive-description") or ""
        vault_name = self.parsed_url.path.split("/")[-2]
        vault = self.glacier_backend.upload_archive(vault_name, self.body, description)
        headers = {
            "x-amz-archive-id": vault["archive_id"],
            "x-amz-sha256-tree-hash": vault["sha256"],
            "status": 201,
        }
        return 201, headers, ""

    def delete_archive(self) -> TYPE_RESPONSE:
        vault_name = self.parsed_url.path.split("/")[-3]
        archive_id = self.parsed_url.path.split("/")[-1]

        self.glacier_backend.delete_archive(vault_name, archive_id)
        return 204, {"status": 204}, ""

    def list_jobs(self) -> TYPE_RESPONSE:
        vault_name = self.parsed_url.path.split("/")[-2]
        jobs = self.glacier_backend.list_jobs(vault_name)
        headers = {"content-type": "application/json"}
        response = json.dumps(
            {"JobList": [job.to_dict() for job in jobs], "Marker": None}
        )
        return 200, headers, response

    def initiate_job(self) -> TYPE_RESPONSE:
        account_id = self.uri.split("/")[1]
        vault_name = self.parsed_url.path.split("/")[-2]
        json_body = json.loads(self.body.decode("utf-8"))
        job_type = json_body["Type"]
        archive_id = json_body.get("ArchiveId")
        tier = json_body.get("Tier") or "Standard"
        job_id = self.glacier_backend.initiate_job(
            vault_name, job_type, tier, archive_id
        )
        headers = {
            "x-amz-job-id": job_id,
            "Location": f"/{account_id}/vaults/{vault_name}/jobs/{job_id}",
            "status": 202,
        }
        return 202, headers, ""

    def describe_job(self) -> str:
        vault_name = self.parsed_url.path.split("/")[-3]
        archive_id = self.parsed_url.path.split("/")[-1]

        job = self.glacier_backend.describe_job(vault_name, archive_id)
        return json.dumps(job.to_dict())  # type: ignore

    def get_job_output(self) -> TYPE_RESPONSE:
        vault_name = self.parsed_url.path.split("/")[-4]
        job_id = self.parsed_url.path.split("/")[-2]
        output = self.glacier_backend.get_job_output(vault_name, job_id)
        if output is None:
            return 404, {"status": 404}, "404 Not Found"
        if isinstance(output, dict):
            headers = {"content-type": "application/json"}
            return 200, headers, json.dumps(output)
        else:
            headers = {"content-type": "application/octet-stream"}
            return 200, headers, output
