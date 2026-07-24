from _typeshed import Incomplete

class ExtendedOperationContainer:
    def __init__(self, connection) -> None: ...

class StandardExtendedOperations(ExtendedOperationContainer):
    def who_am_i(self, controls=None): ...
    def modify_password(self, user=None, old_password=None, new_password=None, hash_algorithm=None, salt=None, controls=None): ...
    def paged_search(
        self,
        search_base,
        search_filter,
        search_scope="SUBTREE",
        dereference_aliases="ALWAYS",
        attributes=None,
        size_limit: int = 0,
        time_limit: int = 0,
        types_only: bool = False,
        get_operational_attributes: bool = False,
        controls=None,
        paged_size: int = 100,
        paged_criticality: bool = False,
        generator: bool = True,
    ): ...
    def persistent_search(
        self,
        search_base: str = "",
        search_filter: str = "(objectclass=*)",
        search_scope="SUBTREE",
        dereference_aliases="NEVER",
        attributes="*",
        size_limit: int = 0,
        time_limit: int = 0,
        controls=None,
        changes_only: bool = True,
        show_additions: bool = True,
        show_deletions: bool = True,
        show_modifications: bool = True,
        show_dn_modifications: bool = True,
        notifications: bool = True,
        streaming: bool = True,
        callback=None,
    ): ...
    def funnel_search(
        self,
        search_base: str = "",
        search_filter: str = "",
        search_scope="SUBTREE",
        dereference_aliases="NEVER",
        attributes="*",
        size_limit: int = 0,
        time_limit: int = 0,
        controls=None,
        streaming: bool = False,
        callback=None,
    ): ...

class NovellExtendedOperations(ExtendedOperationContainer):
    def get_bind_dn(self, controls=None): ...
    def get_universal_password(self, user, controls=None): ...
    def set_universal_password(self, user, new_password=None, controls=None): ...
    def list_replicas(self, server_dn, controls=None): ...
    def partition_entry_count(self, partition_dn, controls=None): ...
    def replica_info(self, server_dn, partition_dn, controls=None): ...
    def start_transaction(self, controls=None): ...
    def end_transaction(self, commit: bool = True, controls=None): ...
    def add_members_to_groups(self, members, groups, fix: bool = True, transaction: bool = True): ...
    def remove_members_from_groups(self, members, groups, fix: bool = True, transaction: bool = True): ...
    def check_groups_memberships(self, members, groups, fix: bool = False, transaction: bool = True): ...

class MicrosoftExtendedOperations(ExtendedOperationContainer):
    def dir_sync(
        self,
        sync_base,
        sync_filter: str = "(objectclass=*)",
        attributes="*",
        cookie=None,
        object_security: bool = False,
        ancestors_first: bool = True,
        public_data_only: bool = False,
        incremental_values: bool = True,
        max_length: int = 2147483647,
        hex_guid: bool = False,
    ): ...
    def modify_password(self, user, new_password, old_password=None, controls=None): ...
    def unlock_account(self, user): ...
    def add_members_to_groups(self, members, groups, fix: bool = True): ...
    def remove_members_from_groups(self, members, groups, fix: bool = True): ...
    def persistent_search(
        self, search_base: str = "", search_scope="SUBTREE", attributes="*", streaming: bool = True, callback=None
    ): ...

class ExtendedOperationsRoot(ExtendedOperationContainer):
    standard: Incomplete
    novell: Incomplete
    microsoft: Incomplete
    def __init__(self, connection) -> None: ...
