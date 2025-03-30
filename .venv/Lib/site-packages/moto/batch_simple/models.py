import datetime
from os import getenv
from time import sleep
from typing import Any, Dict, List, Optional, Tuple

from moto.batch.exceptions import ClientException
from moto.batch.models import BatchBackend, Job, batch_backends
from moto.core.base_backend import BackendDict


class BatchSimpleBackend(BatchBackend):
    """
    Implements a Batch-Backend that does not use Docker containers. Submitted Jobs are marked as Success by default.

    Set the environment variable MOTO_SIMPLE_BATCH_FAIL_AFTER=0 to fail jobs immediately, or set this variable to a positive integer to control after how many seconds the job fails.

    Annotate your tests with `@mock_aws(config={"batch": {"use_docker": False}})`-decorator to use this Batch-implementation.
    """

    @property
    def backend(self) -> BatchBackend:
        return batch_backends[self.account_id][self.region_name]

    def __getattribute__(self, name: str) -> Any:
        """
        Magic part that makes this class behave like a wrapper around the regular batch_backend
        We intercept calls to `submit_job` and replace this with our own (non-Docker) implementation
        Every other method call is send through to batch_backend
        """
        if name in [
            "backend",
            "account_id",
            "region_name",
            "urls",
            "_url_module",
            "__class__",
            "url_bases",
        ]:
            return object.__getattribute__(self, name)
        if name in ["submit_job", "_mark_job_as_finished"]:

            def newfunc(*args: Any, **kwargs: Any) -> Any:
                attr = object.__getattribute__(self, name)
                return attr(*args, **kwargs)

            return newfunc
        else:
            return object.__getattribute__(self.backend, name)

    def submit_job(
        self,
        job_name: str,
        job_def_id: str,
        job_queue: str,
        array_properties: Dict[str, Any],
        depends_on: Optional[List[Dict[str, str]]] = None,
        container_overrides: Optional[Dict[str, Any]] = None,
        timeout: Optional[Dict[str, int]] = None,
        parameters: Optional[Dict[str, str]] = None,
    ) -> Tuple[str, str, str]:
        # Look for job definition
        job_def = self.get_job_definition(job_def_id)
        if job_def is None:
            raise ClientException(f"Job definition {job_def_id} does not exist")

        queue = self.get_job_queue(job_queue)
        if queue is None:
            raise ClientException(f"Job queue {job_queue} does not exist")

        job = Job(
            job_name,
            job_def,
            queue,
            log_backend=self.logs_backend,
            container_overrides=container_overrides,
            depends_on=depends_on,
            all_jobs=self._jobs,
            timeout=timeout,
            array_properties=array_properties,
            parameters=parameters,
        )

        if "size" in array_properties:
            child_jobs: List[Job] = []
            for array_index in range(array_properties.get("size", 0)):
                provided_job_id = f"{job.job_id}:{array_index}"
                child_job = Job(
                    job_name,
                    job_def,
                    queue,
                    log_backend=self.logs_backend,
                    container_overrides=container_overrides,
                    depends_on=depends_on,
                    all_jobs=self._jobs,
                    timeout=timeout,
                    array_properties={"statusSummary": {}, "index": array_index},
                    provided_job_id=provided_job_id,
                    parameters=parameters,
                )
                self._mark_job_as_finished(include_start_attempt=True, job=child_job)
                child_jobs.append(child_job)
            self._mark_job_as_finished(include_start_attempt=False, job=job)
            job._child_jobs = child_jobs
        else:
            self._mark_job_as_finished(include_start_attempt=True, job=job)

        return job_name, job.job_id, job.arn

    def _mark_job_as_finished(self, include_start_attempt: bool, job: Job) -> None:
        self.backend._jobs[job.job_id] = job
        job.job_started_at = datetime.datetime.now()
        job.log_stream_name = job._stream_name
        if include_start_attempt:
            job._start_attempt()
        # We don't want to actually run the job - just mark it as succeeded or failed
        # depending on whether env var MOTO_SIMPLE_BATCH_FAIL_AFTER is set
        # if MOTO_SIMPLE_BATCH_FAIL_AFTER is set to an integer then batch will
        # sleep this many seconds
        should_batch_fail = getenv("MOTO_SIMPLE_BATCH_FAIL_AFTER")
        if should_batch_fail:
            try:
                batch_fail_delay = int(should_batch_fail)
                sleep(batch_fail_delay)
            except ValueError:
                # Unable to parse value of MOTO_SIMPLE_BATCH_FAIL_AFTER as an integer
                pass

            # fail the job
            job._mark_stopped(success=False)
        else:
            job._mark_stopped(success=True)


batch_simple_backends = BackendDict(BatchSimpleBackend, "batch")
