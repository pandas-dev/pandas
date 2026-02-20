"""
Functions written in pure python
"""
from copy import deepcopy, copy
import warnings

import narwhals.stable.v2 as nw

# Functions to deal with value labels

def set_value_labels(dataframe, metadata, formats_as_category=True, formats_as_ordered_category=False):
    """
    Changes the values in the dataframe according to the value formats in the metadata.
    It will return a copy of the dataframe. If no appropiate formats were found, the result will be an unchanged copy
    of the original dataframe.

    Parameters
    ----------
        dataframe : pandas dataframe
            resulting from parsing a file
        metadata : dictionary
            resulting from parsing a file
        formats_as_category : bool, optional
            defaults to True. If True the variables having formats will be transformed into pandas categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into pandas ordered categories.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.

    Returns
    -------
        df_copy : pandas dataframe
            a copy of the original dataframe with the values changed, if appropiate formats were found, unaltered
            otherwise
    """

    df_copy = nw.from_native(dataframe).clone()

    if metadata.value_labels and metadata.variable_to_label:
        for var_name, label_name in metadata.variable_to_label.items():
            labels = metadata.value_labels.get(label_name)
            if labels:
                labels = deepcopy(labels)
                if var_name in df_copy.columns:
                    # unique does not work for polars Object
                    if not df_copy.implementation.is_pandas() and df_copy[var_name].dtype == nw.Object:
                        unvals = list(set(df_copy[var_name].to_list()))
                    else:
                        unvals = df_copy[var_name].unique()
                    for uval in unvals:
                        if uval not in labels:
                            labels[uval] = uval
                    # if all values are null, there will be nothing to replace. However we cannot do replace_strict on null dtype, it raises an error
                    if not df_copy.implementation.is_pandas() and (len(df_copy[var_name])==df_copy[var_name].null_count()):
                        continue
                    # replace_strict requires that all the values are in the map. Could not get map_batches or when/then/otherwise to work
                    elif not df_copy.implementation.is_pandas() and (df_copy[var_name].dtype==nw.Object or not all([type(v)==type(list(labels.values())[0]) for v in labels.values() if v is not None])):
                        # polars is very difficult to convince to mix strings and numbers, so we have to do it this way
                        temp = [labels[x] for x in df_copy[var_name]]
                        newser = nw.new_series(name=var_name, values= temp, dtype=nw.Object, backend=df_copy.implementation) 
                        df_copy = df_copy.with_columns(newser.alias(var_name))
                        if formats_as_category or formats_as_ordered_category:
                            msg = f"You requested formats_as_category=True or formats_as_ordered_category=True, but it was not possible to cast variable '{var_name}' to category"
                            warnings.warn(msg, RuntimeWarning)
                            continue
                    # not sure if we get into this situation ever or what would exactly happen, maybe this is not needed?
                    elif not df_copy.implementation.is_pandas() and df_copy[var_name].dtype==nw.Unknown:
                        msg = f"It was not possible to apply value formats to variable '{var_name}' due to unknown/not supported data type"
                        warnings.warn(msg, RuntimeWarning)
                        continue
                    else:
                        df_copy = df_copy.with_columns(nw.col(var_name).replace_strict(labels))
                    if formats_as_ordered_category:
                        categories = list(set(labels.values()))
                        original_values = list(labels.keys())
                        original_values.sort()
                        revdict= dict()
                        for orival in original_values:
                            curcat = labels.get(orival)
                            if not revdict.get(curcat):
                                revdict[curcat] = orival
                        categories.sort(key=revdict.get)
                        df_copy = df_copy.with_columns(nw.col(var_name).cast(nw.Enum(categories)))
                    elif formats_as_category:
                        df_copy = df_copy.with_columns(nw.col(var_name).cast(nw.Categorical))


    return df_copy.to_native()

def set_catalog_to_sas(sas_dataframe, sas_metadata, catalog_metadata, formats_as_category=True,  
                       formats_as_ordered_category=False):
    """
    Changes the values in the dataframe and sas_metadata according to the formats in the catalog.
    It will return a copy of the dataframe and metadata. If no appropriate formats were found, the result will
    be an unchanged copy of the original dataframe.

    Parameters
    ----------
        sas_dataframe : pandas dataframe
            resulting from parsing a sas7bdat file
        sas_metadata : pyreadstat metadata object
            resulting from parsing a sas7bdat file
        catalog_metadata : pyreadstat metadata object
            resulting from parsing a sas7bcat (catalog) file
        formats_as_category : bool, optional
            defaults to True. If True the variables having formats will be transformed into pandas categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into pandas ordered categories.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.

    Returns
    -------
        df_copy : pandas dataframe
            a copy of the original dataframe with the values changed, if appropriate formats were found, unaltered
            otherwise
        metadata : dict
            a copy of the original sas_metadata enriched with catalog information if found, otherwise unaltered
    """

    if catalog_metadata.value_labels and sas_metadata.variable_to_label:
        catalog_metadata_copy = deepcopy(catalog_metadata)
        metadata = deepcopy(sas_metadata)
        metadata.value_labels = catalog_metadata_copy.value_labels
        df_copy = set_value_labels(sas_dataframe, metadata, formats_as_category=formats_as_category,
                                   formats_as_ordered_category=formats_as_ordered_category)

        variable_value_labels = dict()
        for var_name, var_label in metadata.variable_to_label.items():
            current_labels = catalog_metadata_copy.value_labels.get(var_label)
            if current_labels:
                variable_value_labels[var_name] = current_labels
        metadata.variable_value_labels = variable_value_labels

    else:
        #df_copy = sas_dataframe.copy()
        df_copy = nw.from_native(sas_dataframe).clone().to_native()
        metadata = deepcopy(sas_metadata)

    return df_copy, metadata

