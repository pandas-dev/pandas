from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class ReconciledPayments(Entity):
    id: Incomplete
    foo: Incomplete
    add_on_id: Incomplete

class ContractAddOns(Entity):
    id: Incomplete
    reconciled_payments: Incomplete

r1: ReconciledPayments
r2: ReconciledPayments
c1: ContractAddOns
old_val: Incomplete
