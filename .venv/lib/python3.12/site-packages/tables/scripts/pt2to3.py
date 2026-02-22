"""Utility to helps the migration from PyTables 2.x APIs to 3.x APIs.

The new API is PEP 8 compliant.

"""

import re
import sys
import argparse
from pathlib import Path

old2newnames = dict(  # noqa: C406
    [
        # from __init__.py
        ("hdf5Version", "hdf5_version"),  # data
        # from array.py
        ("parentNode", "parentnode"),  # kwarg
        ("getEnum", "get_enum"),
        ("_initLoop", "_init_loop"),
        ("_fancySelection", "_fancy_selection"),
        ("_checkShape", "_check_shape"),
        ("_readSlice", "_read_slice"),
        ("_readCoords", "_read_coords"),
        ("_readSelection", "_read_selection"),
        ("_writeSlice", "_write_slice"),
        ("_writeCoords", "_write_coords"),
        ("_writeSelection", "_write_selection"),
        ("_g_copyWithStats", "_g_copy_with_stats"),
        ("_c_classId", "_c_classid"),  # attr
        # from atom.py
        ("_checkBase", "_checkbase"),
        # from attributeset.py
        ("newSet", "newset"),  # kwarg
        ("copyClass", "copyclass"),  # kwarg
        ("_g_updateNodeLocation", "_g_update_node_location"),
        ("_g_logAdd", "_g_log_add"),
        ("_g_delAndLog", "_g_del_and_log"),
        ("_v__nodeFile", "_v__nodefile"),  # attr (private)
        ("_v__nodePath", "_v__nodepath"),  # attr (private)
        # from carray.py
        # ('parentNode', 'parentnode'),                       # kwarg
        # from description.py
        ("_g_setNestedNamesDescr", "_g_set_nested_names_descr"),
        ("_g_setPathNames", "_g_set_path_names"),
        ("_v_colObjects", "_v_colobjects"),  # attr
        ("_v_nestedFormats", "_v_nested_formats"),  # attr
        ("_v_nestedNames", "_v_nested_names"),  # attr
        ("_v_nestedDescr", "_v_nested_descr"),  # attr
        ("getColsInOrder", "get_cols_in_order"),
        ("joinPaths", "join_paths"),
        ("metaIsDescription", "MetaIsDescription"),
        # from earray.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("_checkShapeAppend", "_check_shape_append"),
        # from expression.py
        ("_exprvarsCache", "_exprvars_cache"),  # attr (private)
        ("_requiredExprVars", "_required_expr_vars"),
        ("setInputsRange", "set_inputs_range"),
        ("setOutput", "set_output"),
        ("setOutputRange", "set_output_range"),
        # from file.py
        ("_opToCode", "_op_to_code"),  # data (private)
        ("_codeToOp", "_code_to_op"),  # data (private)
        ("_transVersion", "_trans_version"),  # data (private)
        ("_transGroupParent", "_trans_group_parent"),  # data (private)
        ("_transGroupName", "_trans_group_name"),  # data (private)
        ("_transGroupPath", "_trans_group_path"),  # data (private)
        ("_actionLogParent", "_action_log_parent"),  # data (private)
        ("_actionLogName", "_action_log_name"),  # data (private)
        ("_actionLogPath", "_action_log_path"),  # data (private)
        ("_transParent", "_trans_parent"),  # data (private)
        ("_transName", "_trans_name"),  # data (private)
        ("_transPath", "_trans_path"),  # data (private)
        ("_shadowParent", "_shadow_parent"),  # data (private)
        ("_shadowName", "_shadow_name"),  # data (private)
        ("_shadowPath", "_shadow_path"),  # data (private)
        ("copyFile", "copy_file"),
        ("openFile", "open_file"),
        ("_getValueFromContainer", "_get_value_from_container"),
        ("__getRootGroup", "__get_root_group"),
        ("rootUEP", "root_uep"),  # attr
        ("_getOrCreatePath", "_get_or_create_path"),
        ("_createPath", "_create_path"),
        ("createGroup", "create_group"),
        ("createTable", "create_table"),
        ("createArray", "create_array"),
        ("createCArray", "create_carray"),
        ("createEArray", "create_earray"),
        ("createVLArray", "create_vlarray"),
        ("createHardLink", "create_hard_link"),
        ("createSoftLink", "create_soft_link"),
        ("createExternalLink", "create_external_link"),
        ("_getNode", "_get_node"),
        ("getNode", "get_node"),
        ("isVisibleNode", "is_visible_node"),
        ("renameNode", "rename_node"),
        ("moveNode", "move_node"),
        ("copyNode", "copy_node"),
        ("removeNode", "remove_node"),
        ("getNodeAttr", "get_node_attr"),
        ("setNodeAttr", "set_node_attr"),
        ("delNodeAttr", "del_node_attr"),
        ("copyNodeAttrs", "copy_node_attrs"),
        ("copyChildren", "copy_children"),
        ("listNodes", "list_nodes"),
        ("iterNodes", "iter_nodes"),
        ("walkNodes", "walk_nodes"),
        ("walkGroups", "walk_groups"),
        ("_checkOpen", "_check_open"),
        ("_isWritable", "_iswritable"),
        ("_checkWritable", "_check_writable"),
        ("_checkGroup", "_check_group"),
        ("isUndoEnabled", "is_undo_enabled"),
        ("_checkUndoEnabled", "_check_undo_enabled"),
        ("_createTransactionGroup", "_create_transaction_group"),
        ("_createTransaction", "_create_transaction"),
        ("_createMark", "_create_mark"),
        ("enableUndo", "enable_undo"),
        ("disableUndo", "disable_undo"),
        ("_getMarkID", "_get_mark_id"),
        ("_getFinalAction", "_get_final_action"),
        ("getCurrentMark", "get_current_mark"),
        ("_updateNodeLocations", "_update_node_locations"),
        # from group.py
        # ('parentNode', 'parentnode'),                       # kwarg
        # ('ptFile', 'ptfile'),                               # kwarg
        ("_getValueFromContainer", "_get_value_from_container"),
        ("_g_postInitHook", "_g_post_init_hook"),
        ("_g_getChildGroupClass", "_g_get_child_group_class"),
        ("_g_getChildLeafClass", "_g_get_child_leaf_class"),
        ("_g_addChildrenNames", "_g_add_children_names"),
        ("_g_checkHasChild", "_g_check_has_child"),
        ("_f_walkNodes", "_f_walknodes"),
        ("_g_widthWarning", "_g_width_warning"),
        ("_g_refNode", "_g_refnode"),
        ("_g_unrefNode", "_g_unrefnode"),
        ("_g_copyChildren", "_g_copy_children"),
        ("_f_getChild", "_f_get_child"),
        ("_f_listNodes", "_f_list_nodes"),
        ("_f_iterNodes", "_f_iter_nodes"),
        ("_f_walkGroups", "_f_walk_groups"),
        ("_g_closeDescendents", "_g_close_descendents"),
        ("_f_copyChildren", "_f_copy_children"),
        ("_v_maxGroupWidth", "_v_max_group_width"),  # attr
        ("_v_objectID", "_v_objectid"),  # attr
        ("_g_loadChild", "_g_load_child"),
        ("childName", "childname"),  # ???
        ("_c_shadowNameRE", "_c_shadow_name_re"),  # attr (private)
        # from hdf5extension.p{yx,xd}
        ("hdf5Extension", "hdf5extension"),
        ("_getFileId", "_get_file_id"),
        ("_flushFile", "_flush_file"),
        ("_closeFile", "_close_file"),
        ("_g_listAttr", "_g_list_attr"),
        ("_g_setAttr", "_g_setattr"),
        ("_g_getAttr", "_g_getattr"),
        ("_g_listGroup", "_g_list_group"),
        ("_g_getGChildAttr", "_g_get_gchild_attr"),
        ("_g_getLChildAttr", "_g_get_lchild_attr"),
        ("_g_flushGroup", "_g_flush_group"),
        ("_g_closeGroup", "_g_close_group"),
        ("_g_moveNode", "_g_move_node"),
        ("_convertTime64", "_convert_time64"),
        ("_createArray", "_create_array"),
        ("_createCArray", "_create_carray"),
        ("_openArray", "_open_array"),
        ("_readArray", "_read_array"),
        ("_g_readSlice", "_g_read_slice"),
        ("_g_readCoords", "_g_read_coords"),
        ("_g_readSelection", "_g_read_selection"),
        ("_g_writeSlice", "_g_write_slice"),
        ("_g_writeCoords", "_g_write_coords"),
        ("_g_writeSelection", "_g_write_selection"),
        # from idxutils.py
        ("calcChunksize", "calc_chunksize"),
        ("infinityF", "infinityf"),  # data
        ("infinityMap", "infinitymap"),  # data
        ("infType", "inftype"),
        ("StringNextAfter", "string_next_after"),
        ("IntTypeNextAfter", "int_type_next_after"),
        ("BoolTypeNextAfter", "bool_type_next_after"),
        # from index.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("defaultAutoIndex", "default_auto_index"),  # data
        ("defaultIndexFilters", "default_index_filters"),  # data
        ("_tableColumnPathnameOfIndex", "_table_column_pathname_of_index"),
        ("_is_CSI", "_is_csi"),
        ("is_CSI", "is_csi"),  # property
        ("appendLastRow", "append_last_row"),
        ("read_sliceLR", "read_slice_lr"),
        ("readSorted", "read_sorted"),
        ("readIndices", "read_indices"),
        ("_processRange", "_process_range"),
        ("searchLastRow", "search_last_row"),
        ("getLookupRange", "get_lookup_range"),
        ("_g_checkName", "_g_check_name"),
        # from indexes.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("_searchBin", "_search_bin"),
        # from indexesextension
        ("indexesExtension", "indexesextension"),
        ("initRead", "initread"),
        ("readSlice", "read_slice"),
        ("_readIndexSlice", "_read_index_slice"),
        ("_initSortedSlice", "_init_sorted_slice"),
        ("_g_readSortedSlice", "_g_read_sorted_slice"),
        ("_readSortedSlice", "_read_sorted_slice"),
        ("getLRUbounds", "get_lru_bounds"),
        ("getLRUsorted", "get_lru_sorted"),
        ("_searchBinNA_b", "_search_bin_na_b"),
        ("_searchBinNA_ub", "_search_bin_na_ub"),
        ("_searchBinNA_s", "_search_bin_na_s"),
        ("_searchBinNA_us", "_search_bin_na_us"),
        ("_searchBinNA_i", "_search_bin_na_i"),
        ("_searchBinNA_ui", "_search_bin_na_ui"),
        ("_searchBinNA_ll", "_search_bin_na_ll"),
        ("_searchBinNA_ull", "_search_bin_na_ull"),
        ("_searchBinNA_e", "_search_bin_na_e"),
        ("_searchBinNA_f", "_search_bin_na_f"),
        ("_searchBinNA_d", "_search_bin_na_d"),
        ("_searchBinNA_g", "_search_bin_na_g"),
        # from leaf.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("objectID", "object_id"),  # property
        ("_processRangeRead", "_process_range_read"),
        ("_pointSelection", "_point_selection"),
        ("isVisible", "isvisible"),
        ("getAttr", "get_attr"),
        ("setAttr", "set_attr"),
        ("delAttr", "del_attr"),
        # from link.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("_g_getLinkClass", "_g_get_link_class"),
        # from linkextension
        ("linkExtension", "linkextension"),
        ("_getLinkClass", "_get_link_class"),
        ("_g_createHardLink", "_g_create_hard_link"),
        # from lrucacheextension
        ("lrucacheExtension", "lrucacheextension"),
        # from misc/enum.py
        ("_checkAndSetPair", "_check_and_set_pair"),
        ("_getContainer", "_get_container"),
        # from misc/proxydict.py
        ("containerRef", "containerref"),  # attr
        # from node.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("_g_logCreate", "_g_log_create"),
        ("_g_preKillHook", "_g_pre_kill_hook"),
        ("_g_checkOpen", "_g_check_open"),
        ("_g_setLocation", "_g_set_location"),
        ("_g_updateLocation", "_g_update_location"),
        ("_g_delLocation", "_g_del_location"),
        ("_g_updateDependent", "_g_update_dependent"),
        ("_g_removeAndLog", "_g_remove_and_log"),
        ("_g_logMove", "_g_log_move"),
        ("oldPathname", "oldpathname"),  # ??
        ("_g_copyAsChild", "_g_copy_as_child"),
        ("_f_isVisible", "_f_isvisible"),
        ("_g_checkGroup", "_g_check_group"),
        ("_g_checkNotContains", "_g_check_not_contains"),
        ("_g_maybeRemove", "_g_maybe_remove"),
        ("_f_getAttr", "_f_getattr"),
        ("_f_setAttr", "_f_setattr"),
        ("_f_delAttr", "_f_delattr"),
        ("_v_maxTreeDepth", "_v_maxtreedepth"),  # attr
        # from nodes/filenode.py
        ("newNode", "new_node"),
        ("openNode", "open_node"),
        ("_lineChunkSize", "_line_chunksize"),  # attr (private)
        ("_lineSeparator", "_line_separator"),  # attr (private)
        # ('getLineSeparator', 'get_line_separator'),        # dropped
        # ('setLineSeparator', 'set_line_separator'),        # dropped
        # ('delLineSeparator', 'del_line_separator'),        # dropped
        # ('lineSeparator', 'line_separator'),                # property -- dropped
        ("_notReadableError", "_not_readable_error"),
        ("_appendZeros", "_append_zeros"),
        ("getAttrs", "_get_attrs"),
        ("setAttrs", "_set_attrs"),
        ("delAttrs", "_del_attrs"),
        ("_setAttributes", "_set_attributes"),
        ("_checkAttributes", "_check_attributes"),
        ("_checkNotClosed", "_check_not_closed"),
        ("__allowedInitKwArgs", "__allowed_init_kwargs"),  # attr (private)
        ("_byteShape", "_byte_shape"),  # attr (private)
        ("_sizeToShape", "_size_to_shape"),  # attr (private)
        ("_vType", "_vtype"),  # attr (private)
        ("_vShape", "_vshape"),  # attr (private)
        # from path.py
        ("parentPath", "parentpath"),  # kwarg
        ("_pythonIdRE", "_python_id_re"),  # attr (private)
        ("_reservedIdRE", "_reserved_id_re"),  # attr (private)
        ("_hiddenNameRE", "_hidden_name_re"),  # attr (private)
        ("_hiddenPathRE", "_hidden_path_re"),  # attr (private)
        ("checkNameValidity", "check_name_validity"),
        ("joinPath", "join_path"),
        ("splitPath", "split_path"),
        ("isVisibleName", "isvisiblename"),
        ("isVisiblePath", "isvisiblepath"),
        # from registry.py
        ("className", "classname"),  # kwarg
        ("classNameDict", "class_name_dict"),  # data
        ("classIdDict", "class_id_dict"),  # data
        ("getClassByName", "get_class_by_name"),
        # from scripts/ptdump.py
        ("dumpLeaf", "dump_leaf"),
        ("dumpGroup", "dump_group"),
        # from scripts/ptrepack.py
        ("newdstGroup", "newdst_group"),
        ("recreateIndexes", "recreate_indexes"),
        ("copyLeaf", "copy_leaf"),
        # from table.py
        # ('parentNode', 'parentnode'),                       # kwarg
        ("_nxTypeFromNPType", "_nxtype_from_nptype"),  # data (private)
        ("_npSizeType", "_npsizetype"),  # data (private)
        ("_indexNameOf", "_index_name_of"),
        ("_indexPathnameOf", "_index_pathname_of"),
        ("_indexPathnameOfColumn", "_index_pathname_of_column"),
        ("_indexNameOf_", "_index_name_of_"),
        ("_indexPathnameOf_", "_index_pathname_of_"),
        ("_indexPathnameOfColumn_", "_index_pathname_of_column_"),
        ("_table__setautoIndex", "_table__setautoindex"),
        ("_table__getautoIndex", "_table__getautoindex"),
        ("_table__autoIndex", "_table__autoindex"),  # data (private)
        ("_table__whereIndexed", "_table__where_indexed"),
        ("createIndexesTable", "create_indexes_table"),
        ("createIndexesDescr", "create_indexes_descr"),
        ("_column__createIndex", "_column__create_index"),
        ("_autoIndex", "_autoindex"),  # attr
        ("autoIndex", "autoindex"),  # attr
        ("_useIndex", "_use_index"),
        ("_whereCondition", "_where_condition"),  # attr (private)
        ("_conditionCache", "_condition_cache"),  # attr (private)
        # ('_exprvarsCache', '_exprvars_cache'),
        (
            "_enabledIndexingInQueries",
            "_enabled_indexing_in_queries",
        ),  # attr (private)
        ("_emptyArrayCache", "_empty_array_cache"),  # attr (private)
        ("_getTypeColNames", "_get_type_col_names"),
        ("_getEnumMap", "_get_enum_map"),
        ("_cacheDescriptionData", "_cache_description_data"),
        ("_getColumnInstance", "_get_column_instance"),
        ("_checkColumn", "_check_column"),
        ("_disableIndexingInQueries", "_disable_indexing_in_queries"),
        ("_enableIndexingInQueries", "_enable_indexing_in_queries"),
        # ('_requiredExprVars', '_required_expr_vars'),
        ("_getConditionKey", "_get_condition_key"),
        ("_compileCondition", "_compile_condition"),
        ("willQueryUseIndexing", "will_query_use_indexing"),
        ("readWhere", "read_where"),
        ("whereAppend", "append_where"),
        ("getWhereList", "get_where_list"),
        ("_check_sortby_CSI", "_check_sortby_csi"),
        ("_readCoordinates", "_read_coordinates"),
        ("readCoordinates", "read_coordinates"),
        ("_saveBufferedRows", "_save_buffered_rows"),
        ("modifyCoordinates", "modify_coordinates"),
        ("modifyRows", "modify_rows"),
        ("modifyColumn", "modify_column"),
        ("modifyColumns", "modify_columns"),
        ("flushRowsToIndex", "flush_rows_to_index"),
        ("_addRowsToIndex", "_add_rows_to_index"),
        ("removeRows", "remove_rows"),
        ("_setColumnIndexing", "_set_column_indexing"),
        ("_markColumnsAsDirty", "_mark_columns_as_dirty"),
        ("_reIndex", "_reindex"),
        ("_doReIndex", "_do_reindex"),
        ("reIndex", "reindex"),
        ("reIndexDirty", "reindex_dirty"),
        ("_g_copyRows", "_g_copy_rows"),
        ("_g_copyRows_optim", "_g_copy_rows_optim"),
        ("_g_propIndexes", "_g_prop_indexes"),
        ("_g_updateTableLocation", "_g_update_table_location"),
        ("_tableFile", "_table_file"),  # attr (private)
        ("_tablePath", "_table_path"),  # attr (private)
        ("createIndex", "create_index"),
        ("createCSIndex", "create_csindex"),
        ("removeIndex", "remove_index"),
        # from tableextension
        ("tableExtension", "tableextension"),
        ("getNestedFieldCache", "get_nested_field_cache"),
        ("getNestedType", "get_nested_type"),
        ("_createTable", "_create_table"),
        ("_getInfo", "_get_info"),
        ("indexChunk", "indexchunk"),  # attr
        ("indexValid", "indexvalid"),  # attr
        ("indexValues", "indexvalues"),  # attr
        ("bufcoordsData", "bufcoords_data"),  # attr
        ("indexValuesData", "index_values_data"),  # attr
        ("chunkmapData", "chunkmap_data"),  # attr
        ("indexValidData", "index_valid_data"),  # attr
        ("whereCond", "wherecond"),  # attr
        ("iterseqMaxElements", "iterseq_max_elements"),  # attr
        ("IObuf", "iobuf"),  # attr
        ("IObufcpy", "iobufcpy"),  # attr
        ("_convertTime64_", "_convert_time64_"),
        ("_convertTypes", "_convert_types"),
        ("_newBuffer", "_new_buffer"),
        ("__next__inKernel", "__next__inkernel"),
        ("_fillCol", "_fill_col"),
        ("_flushBufferedRows", "_flush_buffered_rows"),
        ("_getUnsavedNrows", "_get_unsaved_nrows"),
        ("_flushModRows", "_flush_mod_rows"),
        # from undoredo.py
        ("moveToShadow", "move_to_shadow"),
        ("moveFromShadow", "move_from_shadow"),
        ("undoCreate", "undo_create"),
        ("redoCreate", "redo_create"),
        ("undoRemove", "undo_remove"),
        ("redoRemove", "redo_remove"),
        ("undoMove", "undo_move"),
        ("redoMove", "redo_move"),
        ("attrToShadow", "attr_to_shadow"),
        ("attrFromShadow", "attr_from_shadow"),
        ("undoAddAttr", "undo_add_attr"),
        ("redoAddAttr", "redo_add_attr"),
        ("undoDelAttr", "undo_del_attr"),
        ("redoDelAttr", "redo_del_attr"),
        # from utils.py
        ("convertToNPAtom", "convert_to_np_atom"),
        ("convertToNPAtom2", "convert_to_np_atom2"),
        ("checkFileAccess", "check_file_access"),
        ("logInstanceCreation", "log_instance_creation"),
        ("fetchLoggedInstances", "fetch_logged_instances"),
        ("countLoggedInstances", "count_logged_instances"),
        ("listLoggedInstances", "list_logged_instances"),
        ("dumpLoggedInstances", "dump_logged_instances"),
        ("detectNumberOfCores", "detect_number_of_cores"),
        # from utilsextension
        ("utilsExtension", "utilsextension"),
        ("PTTypeToHDF5", "pttype_to_hdf5"),  # data
        ("PTSpecialKinds", "pt_special_kinds"),  # data
        ("NPExtPrefixesToPTKinds", "npext_prefixes_to_ptkinds"),  # data
        ("HDF5ClassToString", "hdf5_class_to_string"),  # data
        ("setBloscMaxThreads", "set_blosc_max_threads"),
        ("silenceHDF5Messages", "silence_hdf5_messages"),
        ("isHDF5File", "is_hdf5_file"),
        ("isPyTablesFile", "is_pytables_file"),
        ("getHDF5Version", "get_hdf5_version"),
        ("getPyTablesVersion", "get_pytables_version"),
        ("whichLibVersion", "which_lib_version"),
        ("whichClass", "which_class"),
        ("getNestedField", "get_nested_field"),
        ("getFilters", "get_filters"),
        ("getTypeEnum", "get_type_enum"),
        ("enumFromHDF5", "enum_from_hdf5"),
        ("enumToHDF5", "enum_to_hdf5"),
        ("AtomToHDF5Type", "atom_to_hdf5_type"),
        ("loadEnum", "load_enum"),
        ("HDF5ToNPNestedType", "hdf5_to_np_nested_type"),
        ("HDF5ToNPExtType", "hdf5_to_np_ext_type"),
        ("AtomFromHDF5Type", "atom_from_hdf5_type"),
        ("createNestedType", "create_nested_type"),
        # from unimlemented.py
        ("_openUnImplemented", "_open_unimplemented"),
        # from vlarray.py
        # ('parentNode', 'parentnode'),                       # kwarg
        # ('expectedsizeinMB', 'expected_mb'),                # --> expectedrows
        # ('_v_expectedsizeinMB', '_v_expected_mb'),          # --> expectedrows
    ]
)

new2oldnames = {v: k for k, v in old2newnames.items()}

# Note that it is tempting to use the ast module here, but then this
# breaks transforming cython files.  So instead we are going to do the
# dumb thing with replace.


def make_subs(ns):
    """Make stubs."""
    names = new2oldnames if ns.reverse else old2newnames
    s = r"(?<=\W)({})(?=\W)".format("|".join(list(names)))
    if ns.ignore_previous:
        s += r"(?!\s*?=\s*?previous_api(_property)?\()"
        s += r"(?!\* to \*\w+\*)"
        s += r"(?!\* parameter has been renamed into \*\w+\*\.)"
        s += r"(?! is pending deprecation, import \w+ instead\.)"
    subs = re.compile(s, flags=re.MULTILINE)

    def repl(m):
        return names.get(m.group(1), m.group(0))

    return subs, repl


def main():
    """Implement the main CLI interface."""
    desc = (
        "PyTables 2.x -> 3.x API transition tool\n\n"
        "This tool displays to standard out, so it is \n"
        "common to pipe this to another file:\n\n"
        "$ pt2to3 oldfile.py > newfile.py"
    )
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        default=False,
        dest="reverse",
        help="reverts changes, going from 3.x -> 2.x.",
    )
    parser.add_argument(
        "-p",
        "--no-ignore-previous",
        action="store_false",
        default=True,
        dest="ignore_previous",
        help="ignores previous_api() calls.",
    )
    parser.add_argument(
        "-o", default=None, dest="output", help="output file to write to."
    )
    parser.add_argument(
        "-i",
        "--inplace",
        action="store_true",
        default=False,
        dest="inplace",
        help="overwrites the file in-place.",
    )
    parser.add_argument("filename", help="path to input file.")
    ns = parser.parse_args()

    if not Path(ns.filename).is_file():
        sys.exit(f"file {ns.filename!r} not found")
    src = Path(ns.filename).read_text()

    subs, repl = make_subs(ns)
    targ = subs.sub(repl, src)

    ns.output = ns.filename if ns.inplace else ns.output
    if ns.output is None:
        sys.stdout.write(targ)
    else:
        Path(ns.output).write_text(targ)


if __name__ == "__main__":
    main()
