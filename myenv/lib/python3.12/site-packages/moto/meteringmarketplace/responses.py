import json

from moto.core.responses import BaseResponse

from .models import MeteringMarketplaceBackend, meteringmarketplace_backends


class MarketplaceMeteringResponse(BaseResponse):
    @property
    def backend(self) -> MeteringMarketplaceBackend:
        return meteringmarketplace_backends[self.current_account][self.region]

    def batch_meter_usage(self) -> str:
        usage_records = json.loads(self.body)["UsageRecords"]
        product_code = json.loads(self.body)["ProductCode"]
        results = self.backend.batch_meter_usage(product_code, usage_records)
        return json.dumps({"Results": results, "UnprocessedRecords": []})
