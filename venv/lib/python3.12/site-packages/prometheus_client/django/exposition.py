import os

from django.http import HttpResponse
from django.views import View

import prometheus_client
from prometheus_client import multiprocess
from prometheus_client.exposition import _bake_output


class PrometheusDjangoView(View):
    multiprocess_mode: bool = "PROMETHEUS_MULTIPROC_DIR" in os.environ or "prometheus_multiproc_dir" in os.environ
    registry: prometheus_client.CollectorRegistry = None

    def get(self, request, *args, **kwargs):
        if self.registry is None:
            if self.multiprocess_mode:
                self.registry = prometheus_client.CollectorRegistry()
                multiprocess.MultiProcessCollector(self.registry)
            else:
                self.registry = prometheus_client.REGISTRY
        accept_header = request.headers.get("Accept")
        accept_encoding_header = request.headers.get("Accept-Encoding")
        # Bake output
        status, headers, output = _bake_output(
            registry=self.registry,
            accept_header=accept_header,
            accept_encoding_header=accept_encoding_header,
            params=request.GET,
            disable_compression=False,
        )
        status = int(status.split(" ")[0])
        return HttpResponse(
            output,
            status=status,
            headers=headers,
        )

    def options(self, request, *args, **kwargs):
        return HttpResponse(
            status=200,
            headers={"Allow": "OPTIONS,GET"},
        )
