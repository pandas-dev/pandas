from braintree.resource import Resource

class SettlementBatchSummary(Resource):
    @staticmethod
    def generate(settlement_date, group_by_custom_field=None): ...
