(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1871,7634,9531,8701],{

/***/ 40018:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  ZP: () => (/* binding */ lib)
});

// UNUSED EXPORTS: getDefaultRegistry, withTheme

// EXTERNAL MODULE: ../node_modules/react/jsx-runtime.js
var jsx_runtime = __webpack_require__(24246);
// EXTERNAL MODULE: consume shared module (default) react@~18.2.0 (singleton) (fallback: ../node_modules/react/index.js)
var index_js_ = __webpack_require__(78156);
// EXTERNAL MODULE: consume shared module (default) @rjsf/utils@^5.13.4 (strict) (fallback: ../node_modules/@rjsf/utils/lib/index.js)
var lib_index_js_ = __webpack_require__(24885);
// EXTERNAL MODULE: ../node_modules/lodash/get.js
var get = __webpack_require__(99729);
var get_default = /*#__PURE__*/__webpack_require__.n(get);
// EXTERNAL MODULE: ../node_modules/lodash/isEmpty.js
var isEmpty = __webpack_require__(90104);
var isEmpty_default = /*#__PURE__*/__webpack_require__.n(isEmpty);
// EXTERNAL MODULE: ../node_modules/lodash/pick.js
var pick = __webpack_require__(14648);
var pick_default = /*#__PURE__*/__webpack_require__.n(pick);
// EXTERNAL MODULE: ../node_modules/lodash/toPath.js
var toPath = __webpack_require__(40110);
var toPath_default = /*#__PURE__*/__webpack_require__.n(toPath);
// EXTERNAL MODULE: ../node_modules/lodash/cloneDeep.js
var cloneDeep = __webpack_require__(30454);
var cloneDeep_default = /*#__PURE__*/__webpack_require__.n(cloneDeep);
// EXTERNAL MODULE: ../node_modules/lodash/isObject.js
var isObject = __webpack_require__(11611);
var isObject_default = /*#__PURE__*/__webpack_require__.n(isObject);
// EXTERNAL MODULE: ../node_modules/lodash/set.js
var set = __webpack_require__(47215);
var set_default = /*#__PURE__*/__webpack_require__.n(set);
;// CONCATENATED MODULE: ../node_modules/nanoid/index.browser.js
// This file replaces `index.js` in bundlers like webpack or Rollup,
// according to `browser` config in `package.json`.



let random = bytes => crypto.getRandomValues(new Uint8Array(bytes))

let customRandom = (alphabet, defaultSize, getRandom) => {
  // First, a bitmask is necessary to generate the ID. The bitmask makes bytes
  // values closer to the alphabet size. The bitmask calculates the closest
  // `2^31 - 1` number, which exceeds the alphabet size.
  // For example, the bitmask for the alphabet size 30 is 31 (00011111).
  // `Math.clz32` is not used, because it is not available in browsers.
  let mask = (2 << (Math.log(alphabet.length - 1) / Math.LN2)) - 1
  // Though, the bitmask solution is not perfect since the bytes exceeding
  // the alphabet size are refused. Therefore, to reliably generate the ID,
  // the random bytes redundancy has to be satisfied.

  // Note: every hardware random generator call is performance expensive,
  // because the system call for entropy collection takes a lot of time.
  // So, to avoid additional system calls, extra bytes are requested in advance.

  // Next, a step determines how many random bytes to generate.
  // The number of random bytes gets decided upon the ID size, mask,
  // alphabet size, and magic number 1.6 (using 1.6 peaks at performance
  // according to benchmarks).

  // `-~f => Math.ceil(f)` if f is a float
  // `-~i => i + 1` if i is an integer
  let step = -~((1.6 * mask * defaultSize) / alphabet.length)

  return (size = defaultSize) => {
    let id = ''
    while (true) {
      let bytes = getRandom(step)
      // A compact alternative for `for (var i = 0; i < step; i++)`.
      let j = step | 0
      while (j--) {
        // Adding `|| ''` refuses a random byte that exceeds the alphabet size.
        id += alphabet[bytes[j] & mask] || ''
        if (id.length === size) return id
      }
    }
  }
}

let customAlphabet = (alphabet, size = 21) =>
  customRandom(alphabet, size, random)

let nanoid = (size = 21) =>
  crypto.getRandomValues(new Uint8Array(size)).reduce((id, byte) => {
    // It is incorrect to use bytes exceeding the alphabet size.
    // The following mask reduces the random byte in the 0-255 value
    // range to the 0-63 value range. Therefore, adding hacks, such
    // as empty string fallback or magic numbers, is unneccessary because
    // the bitmask trims bytes down to the alphabet size.
    byte &= 63
    if (byte < 36) {
      // `0-9a-z`
      id += byte.toString(36)
    } else if (byte < 62) {
      // `A-Z`
      id += (byte - 26).toString(36).toUpperCase()
    } else if (byte > 62) {
      id += '-'
    } else {
      id += '_'
    }
    return id
  }, '')



;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/ArrayField.js








/** Used to generate a unique ID for an element in a row */
function generateRowId() {
    return nanoid();
}
/** Converts the `formData` into `KeyedFormDataType` data, using the `generateRowId()` function to create the key
 *
 * @param formData - The data for the form
 * @returns - The `formData` converted into a `KeyedFormDataType` element
 */
function generateKeyedFormData(formData) {
    return !Array.isArray(formData)
        ? []
        : formData.map((item) => {
            return {
                key: generateRowId(),
                item,
            };
        });
}
/** Converts `KeyedFormDataType` data into the inner `formData`
 *
 * @param keyedFormData - The `KeyedFormDataType` to be converted
 * @returns - The inner `formData` item(s) in the `keyedFormData`
 */
function keyedToPlainFormData(keyedFormData) {
    if (Array.isArray(keyedFormData)) {
        return keyedFormData.map((keyedItem) => keyedItem.item);
    }
    return [];
}
/** The `ArrayField` component is used to render a field in the schema that is of type `array`. It supports both normal
 * and fixed array, allowing user to add and remove elements from the array data.
 */
class ArrayField extends index_js_.Component {
    /** Constructs an `ArrayField` from the `props`, generating the initial keyed data from the `formData`
     *
     * @param props - The `FieldProps` for this template
     */
    constructor(props) {
        super(props);
        /** Returns the default form information for an item based on the schema for that item. Deals with the possibility
         * that the schema is fixed and allows additional items.
         */
        this._getNewFormDataRow = () => {
            const { schema, registry } = this.props;
            const { schemaUtils } = registry;
            let itemSchema = schema.items;
            if ((0,lib_index_js_.isFixedItems)(schema) && (0,lib_index_js_.allowAdditionalItems)(schema)) {
                itemSchema = schema.additionalItems;
            }
            // Cast this as a T to work around schema utils being for T[] caused by the FieldProps<T[], S, F> call on the class
            return schemaUtils.getDefaultFormState(itemSchema);
        };
        /** Callback handler for when the user clicks on the add button. Creates a new row of keyed form data at the end of
         * the list, adding it into the state, and then returning `onChange()` with the plain form data converted from the
         * keyed data
         *
         * @param event - The event for the click
         */
        this.onAddClick = (event) => {
            this._handleAddClick(event);
        };
        /** Callback handler for when the user clicks on the add button on an existing array element. Creates a new row of
         * keyed form data inserted at the `index`, adding it into the state, and then returning `onChange()` with the plain
         * form data converted from the keyed data
         *
         * @param index - The index at which the add button is clicked
         */
        this.onAddIndexClick = (index) => {
            return (event) => {
                this._handleAddClick(event, index);
            };
        };
        /** Callback handler for when the user clicks on the copy button on an existing array element. Clones the row of
         * keyed form data at the `index` into the next position in the state, and then returning `onChange()` with the plain
         * form data converted from the keyed data
         *
         * @param index - The index at which the copy button is clicked
         */
        this.onCopyIndexClick = (index) => {
            return (event) => {
                if (event) {
                    event.preventDefault();
                }
                const { onChange, errorSchema } = this.props;
                const { keyedFormData } = this.state;
                // refs #195: revalidate to ensure properly reindexing errors
                let newErrorSchema;
                if (errorSchema) {
                    newErrorSchema = {};
                    for (const idx in errorSchema) {
                        const i = parseInt(idx);
                        if (i <= index) {
                            set_default()(newErrorSchema, [i], errorSchema[idx]);
                        }
                        else if (i > index) {
                            set_default()(newErrorSchema, [i + 1], errorSchema[idx]);
                        }
                    }
                }
                const newKeyedFormDataRow = {
                    key: generateRowId(),
                    item: cloneDeep_default()(keyedFormData[index].item),
                };
                const newKeyedFormData = [...keyedFormData];
                if (index !== undefined) {
                    newKeyedFormData.splice(index + 1, 0, newKeyedFormDataRow);
                }
                else {
                    newKeyedFormData.push(newKeyedFormDataRow);
                }
                this.setState({
                    keyedFormData: newKeyedFormData,
                    updatedKeyedFormData: true,
                }, () => onChange(keyedToPlainFormData(newKeyedFormData), newErrorSchema));
            };
        };
        /** Callback handler for when the user clicks on the remove button on an existing array element. Removes the row of
         * keyed form data at the `index` in the state, and then returning `onChange()` with the plain form data converted
         * from the keyed data
         *
         * @param index - The index at which the remove button is clicked
         */
        this.onDropIndexClick = (index) => {
            return (event) => {
                if (event) {
                    event.preventDefault();
                }
                const { onChange, errorSchema } = this.props;
                const { keyedFormData } = this.state;
                // refs #195: revalidate to ensure properly reindexing errors
                let newErrorSchema;
                if (errorSchema) {
                    newErrorSchema = {};
                    for (const idx in errorSchema) {
                        const i = parseInt(idx);
                        if (i < index) {
                            set_default()(newErrorSchema, [i], errorSchema[idx]);
                        }
                        else if (i > index) {
                            set_default()(newErrorSchema, [i - 1], errorSchema[idx]);
                        }
                    }
                }
                const newKeyedFormData = keyedFormData.filter((_, i) => i !== index);
                this.setState({
                    keyedFormData: newKeyedFormData,
                    updatedKeyedFormData: true,
                }, () => onChange(keyedToPlainFormData(newKeyedFormData), newErrorSchema));
            };
        };
        /** Callback handler for when the user clicks on one of the move item buttons on an existing array element. Moves the
         * row of keyed form data at the `index` to the `newIndex` in the state, and then returning `onChange()` with the
         * plain form data converted from the keyed data
         *
         * @param index - The index of the item to move
         * @param newIndex - The index to where the item is to be moved
         */
        this.onReorderClick = (index, newIndex) => {
            return (event) => {
                if (event) {
                    event.preventDefault();
                    event.currentTarget.blur();
                }
                const { onChange, errorSchema } = this.props;
                let newErrorSchema;
                if (errorSchema) {
                    newErrorSchema = {};
                    for (const idx in errorSchema) {
                        const i = parseInt(idx);
                        if (i == index) {
                            set_default()(newErrorSchema, [newIndex], errorSchema[index]);
                        }
                        else if (i == newIndex) {
                            set_default()(newErrorSchema, [index], errorSchema[newIndex]);
                        }
                        else {
                            set_default()(newErrorSchema, [idx], errorSchema[i]);
                        }
                    }
                }
                const { keyedFormData } = this.state;
                function reOrderArray() {
                    // Copy item
                    const _newKeyedFormData = keyedFormData.slice();
                    // Moves item from index to newIndex
                    _newKeyedFormData.splice(index, 1);
                    _newKeyedFormData.splice(newIndex, 0, keyedFormData[index]);
                    return _newKeyedFormData;
                }
                const newKeyedFormData = reOrderArray();
                this.setState({
                    keyedFormData: newKeyedFormData,
                }, () => onChange(keyedToPlainFormData(newKeyedFormData), newErrorSchema));
            };
        };
        /** Callback handler used to deal with changing the value of the data in the array at the `index`. Calls the
         * `onChange` callback with the updated form data
         *
         * @param index - The index of the item being changed
         */
        this.onChangeForIndex = (index) => {
            return (value, newErrorSchema, id) => {
                const { formData, onChange, errorSchema } = this.props;
                const arrayData = Array.isArray(formData) ? formData : [];
                const newFormData = arrayData.map((item, i) => {
                    // We need to treat undefined items as nulls to have validation.
                    // See https://github.com/tdegrunt/jsonschema/issues/206
                    const jsonValue = typeof value === 'undefined' ? null : value;
                    return index === i ? jsonValue : item;
                });
                onChange(newFormData, errorSchema &&
                    errorSchema && {
                    ...errorSchema,
                    [index]: newErrorSchema,
                }, id);
            };
        };
        /** Callback handler used to change the value for a checkbox */
        this.onSelectChange = (value) => {
            const { onChange, idSchema } = this.props;
            onChange(value, undefined, idSchema && idSchema.$id);
        };
        const { formData = [] } = props;
        const keyedFormData = generateKeyedFormData(formData);
        this.state = {
            keyedFormData,
            updatedKeyedFormData: false,
        };
    }
    /** React lifecycle method that is called when the props are about to change allowing the state to be updated. It
     * regenerates the keyed form data and returns it
     *
     * @param nextProps - The next set of props data
     * @param prevState - The previous set of state data
     */
    static getDerivedStateFromProps(nextProps, prevState) {
        // Don't call getDerivedStateFromProps if keyed formdata was just updated.
        if (prevState.updatedKeyedFormData) {
            return {
                updatedKeyedFormData: false,
            };
        }
        const nextFormData = Array.isArray(nextProps.formData) ? nextProps.formData : [];
        const previousKeyedFormData = prevState.keyedFormData || [];
        const newKeyedFormData = nextFormData.length === previousKeyedFormData.length
            ? previousKeyedFormData.map((previousKeyedFormDatum, index) => {
                return {
                    key: previousKeyedFormDatum.key,
                    item: nextFormData[index],
                };
            })
            : generateKeyedFormData(nextFormData);
        return {
            keyedFormData: newKeyedFormData,
        };
    }
    /** Returns the appropriate title for an item by getting first the title from the schema.items, then falling back to
     * the description from the schema.items, and finally the string "Item"
     */
    get itemTitle() {
        const { schema, registry } = this.props;
        const { translateString } = registry;
        return get_default()(schema, [lib_index_js_.ITEMS_KEY, 'title'], get_default()(schema, [lib_index_js_.ITEMS_KEY, 'description'], translateString(lib_index_js_.TranslatableString.ArrayItemTitle)));
    }
    /** Determines whether the item described in the schema is always required, which is determined by whether any item
     * may be null.
     *
     * @param itemSchema - The schema for the item
     * @return - True if the item schema type does not contain the "null" type
     */
    isItemRequired(itemSchema) {
        if (Array.isArray(itemSchema.type)) {
            // While we don't yet support composite/nullable jsonschema types, it's
            // future-proof to check for requirement against these.
            return !itemSchema.type.includes('null');
        }
        // All non-null array item types are inherently required by design
        return itemSchema.type !== 'null';
    }
    /** Determines whether more items can be added to the array. If the uiSchema indicates the array doesn't allow adding
     * then false is returned. Otherwise, if the schema indicates that there are a maximum number of items and the
     * `formData` matches that value, then false is returned, otherwise true is returned.
     *
     * @param formItems - The list of items in the form
     * @returns - True if the item is addable otherwise false
     */
    canAddItem(formItems) {
        const { schema, uiSchema, registry } = this.props;
        let { addable } = (0,lib_index_js_.getUiOptions)(uiSchema, registry.globalUiOptions);
        if (addable !== false) {
            // if ui:options.addable was not explicitly set to false, we can add
            // another item if we have not exceeded maxItems yet
            if (schema.maxItems !== undefined) {
                addable = formItems.length < schema.maxItems;
            }
            else {
                addable = true;
            }
        }
        return addable;
    }
    /** Callback handler for when the user clicks on the add or add at index buttons. Creates a new row of keyed form data
     * either at the end of the list (when index is not specified) or inserted at the `index` when it is, adding it into
     * the state, and then returning `onChange()` with the plain form data converted from the keyed data
     *
     * @param event - The event for the click
     * @param [index] - The optional index at which to add the new data
     */
    _handleAddClick(event, index) {
        if (event) {
            event.preventDefault();
        }
        const { onChange, errorSchema } = this.props;
        const { keyedFormData } = this.state;
        // refs #195: revalidate to ensure properly reindexing errors
        let newErrorSchema;
        if (errorSchema) {
            newErrorSchema = {};
            for (const idx in errorSchema) {
                const i = parseInt(idx);
                if (index === undefined || i < index) {
                    set_default()(newErrorSchema, [i], errorSchema[idx]);
                }
                else if (i >= index) {
                    set_default()(newErrorSchema, [i + 1], errorSchema[idx]);
                }
            }
        }
        const newKeyedFormDataRow = {
            key: generateRowId(),
            item: this._getNewFormDataRow(),
        };
        const newKeyedFormData = [...keyedFormData];
        if (index !== undefined) {
            newKeyedFormData.splice(index, 0, newKeyedFormDataRow);
        }
        else {
            newKeyedFormData.push(newKeyedFormDataRow);
        }
        this.setState({
            keyedFormData: newKeyedFormData,
            updatedKeyedFormData: true,
        }, () => onChange(keyedToPlainFormData(newKeyedFormData), newErrorSchema));
    }
    /** Renders the `ArrayField` depending on the specific needs of the schema and uischema elements
     */
    render() {
        const { schema, uiSchema, idSchema, registry } = this.props;
        const { schemaUtils, translateString } = registry;
        if (!(lib_index_js_.ITEMS_KEY in schema)) {
            const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema);
            const UnsupportedFieldTemplate = (0,lib_index_js_.getTemplate)('UnsupportedFieldTemplate', registry, uiOptions);
            return ((0,jsx_runtime.jsx)(UnsupportedFieldTemplate, { schema: schema, idSchema: idSchema, reason: translateString(lib_index_js_.TranslatableString.MissingItems), registry: registry }));
        }
        if (schemaUtils.isMultiSelect(schema)) {
            // If array has enum or uniqueItems set to true, call renderMultiSelect() to render the default multiselect widget or a custom widget, if specified.
            return this.renderMultiSelect();
        }
        if ((0,lib_index_js_.isCustomWidget)(uiSchema)) {
            return this.renderCustomWidget();
        }
        if ((0,lib_index_js_.isFixedItems)(schema)) {
            return this.renderFixedArray();
        }
        if (schemaUtils.isFilesArray(schema, uiSchema)) {
            return this.renderFiles();
        }
        return this.renderNormalArray();
    }
    /** Renders a normal array without any limitations of length
     */
    renderNormalArray() {
        const { schema, uiSchema = {}, errorSchema, idSchema, name, disabled = false, readonly = false, autofocus = false, required = false, registry, onBlur, onFocus, idPrefix, idSeparator = '_', rawErrors, } = this.props;
        const { keyedFormData } = this.state;
        const title = schema.title === undefined ? name : schema.title;
        const { schemaUtils, formContext } = registry;
        const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema);
        const _schemaItems = isObject_default()(schema.items) ? schema.items : {};
        const itemsSchema = schemaUtils.retrieveSchema(_schemaItems);
        const formData = keyedToPlainFormData(this.state.keyedFormData);
        const canAdd = this.canAddItem(formData);
        const arrayProps = {
            canAdd,
            items: keyedFormData.map((keyedItem, index) => {
                const { key, item } = keyedItem;
                // While we are actually dealing with a single item of type T, the types require a T[], so cast
                const itemCast = item;
                const itemSchema = schemaUtils.retrieveSchema(_schemaItems, itemCast);
                const itemErrorSchema = errorSchema ? errorSchema[index] : undefined;
                const itemIdPrefix = idSchema.$id + idSeparator + index;
                const itemIdSchema = schemaUtils.toIdSchema(itemSchema, itemIdPrefix, itemCast, idPrefix, idSeparator);
                return this.renderArrayFieldItem({
                    key,
                    index,
                    name: name && `${name}-${index}`,
                    canAdd,
                    canMoveUp: index > 0,
                    canMoveDown: index < formData.length - 1,
                    itemSchema,
                    itemIdSchema,
                    itemErrorSchema,
                    itemData: itemCast,
                    itemUiSchema: uiSchema.items,
                    autofocus: autofocus && index === 0,
                    onBlur,
                    onFocus,
                    rawErrors,
                    totalItems: keyedFormData.length,
                });
            }),
            className: `field field-array field-array-of-${itemsSchema.type}`,
            disabled,
            idSchema,
            uiSchema,
            onAddClick: this.onAddClick,
            readonly,
            required,
            schema,
            title,
            formContext,
            formData,
            rawErrors,
            registry,
        };
        const Template = (0,lib_index_js_.getTemplate)('ArrayFieldTemplate', registry, uiOptions);
        return (0,jsx_runtime.jsx)(Template, { ...arrayProps });
    }
    /** Renders an array using the custom widget provided by the user in the `uiSchema`
     */
    renderCustomWidget() {
        var _a;
        const { schema, idSchema, uiSchema, disabled = false, readonly = false, autofocus = false, required = false, hideError, placeholder, onBlur, onFocus, formData: items = [], registry, rawErrors, name, } = this.props;
        const { widgets, formContext, globalUiOptions, schemaUtils } = registry;
        const { widget, title: uiTitle, ...options } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const Widget = (0,lib_index_js_.getWidget)(schema, widget, widgets);
        const label = (_a = uiTitle !== null && uiTitle !== void 0 ? uiTitle : schema.title) !== null && _a !== void 0 ? _a : name;
        const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
        return ((0,jsx_runtime.jsx)(Widget, { id: idSchema.$id, name: name, multiple: true, onChange: this.onSelectChange, onBlur: onBlur, onFocus: onFocus, options: options, schema: schema, uiSchema: uiSchema, registry: registry, value: items, disabled: disabled, readonly: readonly, hideError: hideError, required: required, label: label, hideLabel: !displayLabel, placeholder: placeholder, formContext: formContext, autofocus: autofocus, rawErrors: rawErrors }));
    }
    /** Renders an array as a set of checkboxes
     */
    renderMultiSelect() {
        var _a;
        const { schema, idSchema, uiSchema, formData: items = [], disabled = false, readonly = false, autofocus = false, required = false, placeholder, onBlur, onFocus, registry, rawErrors, name, } = this.props;
        const { widgets, schemaUtils, formContext, globalUiOptions } = registry;
        const itemsSchema = schemaUtils.retrieveSchema(schema.items, items);
        const enumOptions = (0,lib_index_js_.optionsList)(itemsSchema);
        const { widget = 'select', title: uiTitle, ...options } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const Widget = (0,lib_index_js_.getWidget)(schema, widget, widgets);
        const label = (_a = uiTitle !== null && uiTitle !== void 0 ? uiTitle : schema.title) !== null && _a !== void 0 ? _a : name;
        const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
        return ((0,jsx_runtime.jsx)(Widget, { id: idSchema.$id, name: name, multiple: true, onChange: this.onSelectChange, onBlur: onBlur, onFocus: onFocus, options: { ...options, enumOptions }, schema: schema, uiSchema: uiSchema, registry: registry, value: items, disabled: disabled, readonly: readonly, required: required, label: label, hideLabel: !displayLabel, placeholder: placeholder, formContext: formContext, autofocus: autofocus, rawErrors: rawErrors }));
    }
    /** Renders an array of files using the `FileWidget`
     */
    renderFiles() {
        var _a;
        const { schema, uiSchema, idSchema, name, disabled = false, readonly = false, autofocus = false, required = false, onBlur, onFocus, registry, formData: items = [], rawErrors, } = this.props;
        const { widgets, formContext, globalUiOptions, schemaUtils } = registry;
        const { widget = 'files', title: uiTitle, ...options } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const Widget = (0,lib_index_js_.getWidget)(schema, widget, widgets);
        const label = (_a = uiTitle !== null && uiTitle !== void 0 ? uiTitle : schema.title) !== null && _a !== void 0 ? _a : name;
        const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
        return ((0,jsx_runtime.jsx)(Widget, { options: options, id: idSchema.$id, name: name, multiple: true, onChange: this.onSelectChange, onBlur: onBlur, onFocus: onFocus, schema: schema, uiSchema: uiSchema, value: items, disabled: disabled, readonly: readonly, required: required, registry: registry, formContext: formContext, autofocus: autofocus, rawErrors: rawErrors, label: label, hideLabel: !displayLabel }));
    }
    /** Renders an array that has a maximum limit of items
     */
    renderFixedArray() {
        const { schema, uiSchema = {}, formData = [], errorSchema, idPrefix, idSeparator = '_', idSchema, name, disabled = false, readonly = false, autofocus = false, required = false, registry, onBlur, onFocus, rawErrors, } = this.props;
        const { keyedFormData } = this.state;
        let { formData: items = [] } = this.props;
        const title = schema.title || name;
        const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema);
        const { schemaUtils, formContext } = registry;
        const _schemaItems = isObject_default()(schema.items) ? schema.items : [];
        const itemSchemas = _schemaItems.map((item, index) => schemaUtils.retrieveSchema(item, formData[index]));
        const additionalSchema = isObject_default()(schema.additionalItems)
            ? schemaUtils.retrieveSchema(schema.additionalItems, formData)
            : null;
        if (!items || items.length < itemSchemas.length) {
            // to make sure at least all fixed items are generated
            items = items || [];
            items = items.concat(new Array(itemSchemas.length - items.length));
        }
        // These are the props passed into the render function
        const canAdd = this.canAddItem(items) && !!additionalSchema;
        const arrayProps = {
            canAdd,
            className: 'field field-array field-array-fixed-items',
            disabled,
            idSchema,
            formData,
            items: keyedFormData.map((keyedItem, index) => {
                const { key, item } = keyedItem;
                // While we are actually dealing with a single item of type T, the types require a T[], so cast
                const itemCast = item;
                const additional = index >= itemSchemas.length;
                const itemSchema = (additional && isObject_default()(schema.additionalItems)
                    ? schemaUtils.retrieveSchema(schema.additionalItems, itemCast)
                    : itemSchemas[index]) || {};
                const itemIdPrefix = idSchema.$id + idSeparator + index;
                const itemIdSchema = schemaUtils.toIdSchema(itemSchema, itemIdPrefix, itemCast, idPrefix, idSeparator);
                const itemUiSchema = additional
                    ? uiSchema.additionalItems || {}
                    : Array.isArray(uiSchema.items)
                        ? uiSchema.items[index]
                        : uiSchema.items || {};
                const itemErrorSchema = errorSchema ? errorSchema[index] : undefined;
                return this.renderArrayFieldItem({
                    key,
                    index,
                    name: name && `${name}-${index}`,
                    canAdd,
                    canRemove: additional,
                    canMoveUp: index >= itemSchemas.length + 1,
                    canMoveDown: additional && index < items.length - 1,
                    itemSchema,
                    itemData: itemCast,
                    itemUiSchema,
                    itemIdSchema,
                    itemErrorSchema,
                    autofocus: autofocus && index === 0,
                    onBlur,
                    onFocus,
                    rawErrors,
                    totalItems: keyedFormData.length,
                });
            }),
            onAddClick: this.onAddClick,
            readonly,
            required,
            registry,
            schema,
            uiSchema,
            title,
            formContext,
            rawErrors,
        };
        const Template = (0,lib_index_js_.getTemplate)('ArrayFieldTemplate', registry, uiOptions);
        return (0,jsx_runtime.jsx)(Template, { ...arrayProps });
    }
    /** Renders the individual array item using a `SchemaField` along with the additional properties required to be send
     * back to the `ArrayFieldItemTemplate`.
     *
     * @param props - The props for the individual array item to be rendered
     */
    renderArrayFieldItem(props) {
        const { key, index, name, canAdd, canRemove = true, canMoveUp, canMoveDown, itemSchema, itemData, itemUiSchema, itemIdSchema, itemErrorSchema, autofocus, onBlur, onFocus, rawErrors, totalItems, } = props;
        const { disabled, hideError, idPrefix, idSeparator, readonly, uiSchema, registry, formContext } = this.props;
        const { fields: { ArraySchemaField, SchemaField }, globalUiOptions, } = registry;
        const ItemSchemaField = ArraySchemaField || SchemaField;
        const { orderable = true, removable = true, copyable = false } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const has = {
            moveUp: orderable && canMoveUp,
            moveDown: orderable && canMoveDown,
            copy: copyable && canAdd,
            remove: removable && canRemove,
            toolbar: false,
        };
        has.toolbar = Object.keys(has).some((key) => has[key]);
        return {
            children: ((0,jsx_runtime.jsx)(ItemSchemaField, { name: name, index: index, schema: itemSchema, uiSchema: itemUiSchema, formData: itemData, formContext: formContext, errorSchema: itemErrorSchema, idPrefix: idPrefix, idSeparator: idSeparator, idSchema: itemIdSchema, required: this.isItemRequired(itemSchema), onChange: this.onChangeForIndex(index), onBlur: onBlur, onFocus: onFocus, registry: registry, disabled: disabled, readonly: readonly, hideError: hideError, autofocus: autofocus, rawErrors: rawErrors })),
            className: 'array-item',
            disabled,
            canAdd,
            hasCopy: has.copy,
            hasToolbar: has.toolbar,
            hasMoveUp: has.moveUp,
            hasMoveDown: has.moveDown,
            hasRemove: has.remove,
            index,
            totalItems,
            key,
            onAddIndexClick: this.onAddIndexClick,
            onCopyIndexClick: this.onCopyIndexClick,
            onDropIndexClick: this.onDropIndexClick,
            onReorderClick: this.onReorderClick,
            readonly,
            registry,
            schema: itemSchema,
            uiSchema: itemUiSchema,
        };
    }
}
/** `ArrayField` is `React.ComponentType<FieldProps<T[], S, F>>` (necessarily) but the `registry` requires things to be a
 * `Field` which is defined as `React.ComponentType<FieldProps<T, S, F>>`, so cast it to make `registry` happy.
 */
/* harmony default export */ const fields_ArrayField = (ArrayField);
//# sourceMappingURL=ArrayField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/BooleanField.js



/** The `BooleanField` component is used to render a field in the schema is boolean. It constructs `enumOptions` for the
 * two boolean values based on the various alternatives in the schema.
 *
 * @param props - The `FieldProps` for this template
 */
function BooleanField(props) {
    var _a, _b;
    const { schema, name, uiSchema, idSchema, formData, registry, required, disabled, readonly, hideError, autofocus, onChange, onFocus, onBlur, rawErrors, } = props;
    const { title } = schema;
    const { widgets, formContext, translateString, globalUiOptions } = registry;
    const { widget = 'checkbox', title: uiTitle, 
    // Unlike the other fields, don't use `getDisplayLabel()` since it always returns false for the boolean type
    label: displayLabel = true, ...options } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
    const Widget = (0,lib_index_js_.getWidget)(schema, widget, widgets);
    const yes = translateString(lib_index_js_.TranslatableString.YesLabel);
    const no = translateString(lib_index_js_.TranslatableString.NoLabel);
    let enumOptions;
    const label = (_a = uiTitle !== null && uiTitle !== void 0 ? uiTitle : title) !== null && _a !== void 0 ? _a : name;
    if (Array.isArray(schema.oneOf)) {
        enumOptions = (0,lib_index_js_.optionsList)({
            oneOf: schema.oneOf
                .map((option) => {
                if (isObject_default()(option)) {
                    return {
                        ...option,
                        title: option.title || (option.const === true ? yes : no),
                    };
                }
                return undefined;
            })
                .filter((o) => o), // cast away the error that typescript can't grok is fixed
        });
    }
    else {
        // We deprecated enumNames in v5. It's intentionally omitted from RSJFSchema type, so we need to cast here.
        const schemaWithEnumNames = schema;
        const enums = (_b = schema.enum) !== null && _b !== void 0 ? _b : [true, false];
        if (!schemaWithEnumNames.enumNames && enums.length === 2 && enums.every((v) => typeof v === 'boolean')) {
            enumOptions = [
                {
                    value: enums[0],
                    label: enums[0] ? yes : no,
                },
                {
                    value: enums[1],
                    label: enums[1] ? yes : no,
                },
            ];
        }
        else {
            enumOptions = (0,lib_index_js_.optionsList)({
                enum: enums,
                // NOTE: enumNames is deprecated, but still supported for now.
                enumNames: schemaWithEnumNames.enumNames,
            });
        }
    }
    return ((0,jsx_runtime.jsx)(Widget, { options: { ...options, enumOptions }, schema: schema, uiSchema: uiSchema, id: idSchema.$id, name: name, onChange: onChange, onFocus: onFocus, onBlur: onBlur, label: label, hideLabel: !displayLabel, value: formData, required: required, disabled: disabled, readonly: readonly, hideError: hideError, registry: registry, formContext: formContext, autofocus: autofocus, rawErrors: rawErrors }));
}
/* harmony default export */ const fields_BooleanField = (BooleanField);
//# sourceMappingURL=BooleanField.js.map
// EXTERNAL MODULE: ../node_modules/lodash/omit.js
var omit = __webpack_require__(48159);
var omit_default = /*#__PURE__*/__webpack_require__.n(omit);
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/MultiSchemaField.js






/** The `AnyOfField` component is used to render a field in the schema that is an `anyOf`, `allOf` or `oneOf`. It tracks
 * the currently selected option and cleans up any irrelevant data in `formData`.
 *
 * @param props - The `FieldProps` for this template
 */
class AnyOfField extends index_js_.Component {
    /** Constructs an `AnyOfField` with the given `props` to initialize the initially selected option in state
     *
     * @param props - The `FieldProps` for this template
     */
    constructor(props) {
        super(props);
        /** Callback handler to remember what the currently selected option is. In addition to that the `formData` is updated
         * to remove properties that are not part of the newly selected option schema, and then the updated data is passed to
         * the `onChange` handler.
         *
         * @param option - The new option value being selected
         */
        this.onOptionChange = (option) => {
            const { selectedOption, retrievedOptions } = this.state;
            const { formData, onChange, registry } = this.props;
            const { schemaUtils } = registry;
            const intOption = option !== undefined ? parseInt(option, 10) : -1;
            if (intOption === selectedOption) {
                return;
            }
            const newOption = intOption >= 0 ? retrievedOptions[intOption] : undefined;
            const oldOption = selectedOption >= 0 ? retrievedOptions[selectedOption] : undefined;
            let newFormData = schemaUtils.sanitizeDataForNewSchema(newOption, oldOption, formData);
            if (newFormData && newOption) {
                // Call getDefaultFormState to make sure defaults are populated on change. Pass "excludeObjectChildren"
                // so that only the root objects themselves are created without adding undefined children properties
                newFormData = schemaUtils.getDefaultFormState(newOption, newFormData, 'excludeObjectChildren');
            }
            onChange(newFormData, undefined, this.getFieldId());
            this.setState({ selectedOption: intOption });
        };
        const { formData, options, registry: { schemaUtils }, } = this.props;
        // cache the retrieved options in state in case they have $refs to save doing it later
        const retrievedOptions = options.map((opt) => schemaUtils.retrieveSchema(opt, formData));
        this.state = {
            retrievedOptions,
            selectedOption: this.getMatchingOption(0, formData, retrievedOptions),
        };
    }
    /** React lifecycle method that is called when the props and/or state for this component is updated. It recomputes the
     * currently selected option based on the overall `formData`
     *
     * @param prevProps - The previous `FieldProps` for this template
     * @param prevState - The previous `AnyOfFieldState` for this template
     */
    componentDidUpdate(prevProps, prevState) {
        const { formData, options, idSchema } = this.props;
        const { selectedOption } = this.state;
        let newState = this.state;
        if (!(0,lib_index_js_.deepEquals)(prevProps.options, options)) {
            const { registry: { schemaUtils }, } = this.props;
            // re-cache the retrieved options in state in case they have $refs to save doing it later
            const retrievedOptions = options.map((opt) => schemaUtils.retrieveSchema(opt, formData));
            newState = { selectedOption, retrievedOptions };
        }
        if (!(0,lib_index_js_.deepEquals)(formData, prevProps.formData) && idSchema.$id === prevProps.idSchema.$id) {
            const { retrievedOptions } = newState;
            const matchingOption = this.getMatchingOption(selectedOption, formData, retrievedOptions);
            if (prevState && matchingOption !== selectedOption) {
                newState = { selectedOption: matchingOption, retrievedOptions };
            }
        }
        if (newState !== this.state) {
            this.setState(newState);
        }
    }
    /** Determines the best matching option for the given `formData` and `options`.
     *
     * @param formData - The new formData
     * @param options - The list of options to choose from
     * @return - The index of the `option` that best matches the `formData`
     */
    getMatchingOption(selectedOption, formData, options) {
        const { schema, registry: { schemaUtils }, } = this.props;
        const discriminator = (0,lib_index_js_.getDiscriminatorFieldFromSchema)(schema);
        const option = schemaUtils.getClosestMatchingOption(formData, options, selectedOption, discriminator);
        return option;
    }
    getFieldId() {
        const { idSchema, schema } = this.props;
        return `${idSchema.$id}${schema.oneOf ? '__oneof_select' : '__anyof_select'}`;
    }
    /** Renders the `AnyOfField` selector along with a `SchemaField` for the value of the `formData`
     */
    render() {
        const { name, disabled = false, errorSchema = {}, formContext, onBlur, onFocus, registry, schema, uiSchema, } = this.props;
        const { widgets, fields, translateString, globalUiOptions, schemaUtils } = registry;
        const { SchemaField: _SchemaField } = fields;
        const { selectedOption, retrievedOptions } = this.state;
        const { widget = 'select', placeholder, autofocus, autocomplete, title = schema.title, ...uiOptions } = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const Widget = (0,lib_index_js_.getWidget)({ type: 'number' }, widget, widgets);
        const rawErrors = get_default()(errorSchema, lib_index_js_.ERRORS_KEY, []);
        const fieldErrorSchema = omit_default()(errorSchema, [lib_index_js_.ERRORS_KEY]);
        const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
        const option = selectedOption >= 0 ? retrievedOptions[selectedOption] || null : null;
        let optionSchema;
        if (option) {
            // merge top level required field
            const { required } = schema;
            // Merge in all the non-oneOf/anyOf properties and also skip the special ADDITIONAL_PROPERTY_FLAG property
            optionSchema = required ? (0,lib_index_js_.mergeSchemas)({ required }, option) : option;
        }
        const translateEnum = title
            ? lib_index_js_.TranslatableString.TitleOptionPrefix
            : lib_index_js_.TranslatableString.OptionPrefix;
        const translateParams = title ? [title] : [];
        const enumOptions = retrievedOptions.map((opt, index) => ({
            label: opt.title || translateString(translateEnum, translateParams.concat(String(index + 1))),
            value: index,
        }));
        return ((0,jsx_runtime.jsxs)("div", { className: 'panel panel-default panel-body', children: [(0,jsx_runtime.jsx)("div", { className: 'form-group', children: (0,jsx_runtime.jsx)(Widget, { id: this.getFieldId(), name: `${name}${schema.oneOf ? '__oneof_select' : '__anyof_select'}`, schema: { type: 'number', default: 0 }, onChange: this.onOptionChange, onBlur: onBlur, onFocus: onFocus, disabled: disabled || isEmpty_default()(enumOptions), multiple: false, rawErrors: rawErrors, errorSchema: fieldErrorSchema, value: selectedOption >= 0 ? selectedOption : undefined, options: { enumOptions, ...uiOptions }, registry: registry, formContext: formContext, placeholder: placeholder, autocomplete: autocomplete, autofocus: autofocus, label: title !== null && title !== void 0 ? title : name, hideLabel: !displayLabel }) }), option !== null && (0,jsx_runtime.jsx)(_SchemaField, { ...this.props, schema: optionSchema })] }));
    }
}
/* harmony default export */ const MultiSchemaField = (AnyOfField);
//# sourceMappingURL=MultiSchemaField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/NumberField.js



// Matches a string that ends in a . character, optionally followed by a sequence of
// digits followed by any number of 0 characters up until the end of the line.
// Ensuring that there is at least one prefixed character is important so that
// you don't incorrectly match against "0".
const trailingCharMatcherWithPrefix = /\.([0-9]*0)*$/;
// This is used for trimming the trailing 0 and . characters without affecting
// the rest of the string. Its possible to use one RegEx with groups for this
// functionality, but it is fairly complex compared to simply defining two
// different matchers.
const trailingCharMatcher = /[0.]0*$/;
/**
 * The NumberField class has some special handling for dealing with trailing
 * decimal points and/or zeroes. This logic is designed to allow trailing values
 * to be visible in the input element, but not be represented in the
 * corresponding form data.
 *
 * The algorithm is as follows:
 *
 * 1. When the input value changes the value is cached in the component state
 *
 * 2. The value is then normalized, removing trailing decimal points and zeros,
 *    then passed to the "onChange" callback
 *
 * 3. When the component is rendered, the formData value is checked against the
 *    value cached in the state. If it matches the cached value, the cached
 *    value is passed to the input instead of the formData value
 */
function NumberField(props) {
    const { registry, onChange, formData, value: initialValue } = props;
    const [lastValue, setLastValue] = (0,index_js_.useState)(initialValue);
    const { StringField } = registry.fields;
    let value = formData;
    /** Handle the change from the `StringField` to properly convert to a number
     *
     * @param value - The current value for the change occurring
     */
    const handleChange = (0,index_js_.useCallback)((value) => {
        // Cache the original value in component state
        setLastValue(value);
        // Normalize decimals that don't start with a zero character in advance so
        // that the rest of the normalization logic is simpler
        if (`${value}`.charAt(0) === '.') {
            value = `0${value}`;
        }
        // Check that the value is a string (this can happen if the widget used is a
        // <select>, due to an enum declaration etc) then, if the value ends in a
        // trailing decimal point or multiple zeroes, strip the trailing values
        const processed = typeof value === 'string' && value.match(trailingCharMatcherWithPrefix)
            ? (0,lib_index_js_.asNumber)(value.replace(trailingCharMatcher, ''))
            : (0,lib_index_js_.asNumber)(value);
        onChange(processed);
    }, [onChange]);
    if (typeof lastValue === 'string' && typeof value === 'number') {
        // Construct a regular expression that checks for a string that consists
        // of the formData value suffixed with zero or one '.' characters and zero
        // or more '0' characters
        const re = new RegExp(`${value}`.replace('.', '\\.') + '\\.?0*$');
        // If the cached "lastValue" is a match, use that instead of the formData
        // value to prevent the input value from changing in the UI
        if (lastValue.match(re)) {
            value = lastValue;
        }
    }
    return (0,jsx_runtime.jsx)(StringField, { ...props, formData: value, onChange: handleChange });
}
/* harmony default export */ const fields_NumberField = (NumberField);
//# sourceMappingURL=NumberField.js.map
;// CONCATENATED MODULE: ../node_modules/markdown-to-jsx/dist/index.modern.js
function t(){return t=Object.assign?Object.assign.bind():function(e){for(var t=1;t<arguments.length;t++){var n=arguments[t];for(var r in n)Object.prototype.hasOwnProperty.call(n,r)&&(e[r]=n[r])}return e},t.apply(this,arguments)}const n=["children","options"],r={blockQuote:"0",breakLine:"1",breakThematic:"2",codeBlock:"3",codeFenced:"4",codeInline:"5",footnote:"6",footnoteReference:"7",gfmTask:"8",heading:"9",headingSetext:"10",htmlBlock:"11",htmlComment:"12",htmlSelfClosing:"13",image:"14",link:"15",linkAngleBraceStyleDetector:"16",linkBareUrlDetector:"17",linkMailtoDetector:"18",newlineCoalescer:"19",orderedList:"20",paragraph:"21",ref:"22",refImage:"23",refLink:"24",table:"25",tableSeparator:"26",text:"27",textBolded:"28",textEmphasized:"29",textEscaped:"30",textMarked:"31",textStrikethroughed:"32",unorderedList:"33"};var i;!function(e){e[e.MAX=0]="MAX",e[e.HIGH=1]="HIGH",e[e.MED=2]="MED",e[e.LOW=3]="LOW",e[e.MIN=4]="MIN"}(i||(i={}));const l=["allowFullScreen","allowTransparency","autoComplete","autoFocus","autoPlay","cellPadding","cellSpacing","charSet","classId","colSpan","contentEditable","contextMenu","crossOrigin","encType","formAction","formEncType","formMethod","formNoValidate","formTarget","frameBorder","hrefLang","inputMode","keyParams","keyType","marginHeight","marginWidth","maxLength","mediaGroup","minLength","noValidate","radioGroup","readOnly","rowSpan","spellCheck","srcDoc","srcLang","srcSet","tabIndex","useMap"].reduce((e,t)=>(e[t.toLowerCase()]=t,e),{class:"className",for:"htmlFor"}),a={amp:"&",apos:"'",gt:">",lt:"<",nbsp:" ",quot:"“"},o=["style","script"],c=/([-A-Z0-9_:]+)(?:\s*=\s*(?:(?:"((?:\\.|[^"])*)")|(?:'((?:\\.|[^'])*)')|(?:\{((?:\\.|{[^}]*?}|[^}])*)\})))?/gi,s=/mailto:/i,d=/\n{2,}$/,u=/^(\s*>[\s\S]*?)(?=\n\n|$)/,p=/^ *> ?/gm,f=/^(?:\[!([^\]]*)\]\n)?([\s\S]*)/,h=/^ {2,}\n/,m=/^(?:( *[-*_])){3,} *(?:\n *)+\n/,g=/^(?: {1,3})?(`{3,}|~{3,}) *(\S+)? *([^\n]*?)?\n([\s\S]*?)(?:\1\n?|$)/,y=/^(?: {4}[^\n]+\n*)+(?:\n *)+\n?/,k=/^(`+)((?:\\`|[^`])+)\1/,x=/^(?:\n *)*\n/,b=/\r\n?/g,v=/^\[\^([^\]]+)](:(.*)((\n+ {4,}.*)|(\n(?!\[\^).+))*)/,C=/^\[\^([^\]]+)]/,$=/\f/g,S=/^---[ \t]*\n(.|\n)*\n---[ \t]*\n/,w=/^\s*?\[(x|\s)\]/,E=/^ *(#{1,6}) *([^\n]+?)(?: +#*)?(?:\n *)*(?:\n|$)/,z=/^ *(#{1,6}) +([^\n]+?)(?: +#*)?(?:\n *)*(?:\n|$)/,L=/^([^\n]+)\n *(=|-){3,} *(?:\n *)+\n/,A=/^ *(?!<[a-z][^ >/]* ?\/>)<([a-z][^ >/]*) ?((?:[^>]*[^/])?)>\n?(\s*(?:<\1[^>]*?>[\s\S]*?<\/\1>|(?!<\1\b)[\s\S])*?)<\/\1>(?!<\/\1>)\n*/i,T=/&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-fA-F]{1,6});/gi,B=/^<!--[\s\S]*?(?:-->)/,O=/^(data|aria|x)-[a-z_][a-z\d_.-]*$/,M=/^ *<([a-z][a-z0-9:]*)(?:\s+((?:<.*?>|[^>])*))?\/?>(?!<\/\1>)(\s*\n)?/i,R=/^\{.*\}$/,I=/^(https?:\/\/[^\s<]+[^<.,:;"')\]\s])/,U=/^<([^ >]+@[^ >]+)>/,D=/^<([^ >]+:\/[^ >]+)>/,N=/-([a-z])?/gi,j=/^(\|.*)\n(?: *(\|? *[-:]+ *\|[-| :]*)\n((?:.*\|.*\n)*))?\n?/,H=/^\[([^\]]*)\]:\s+<?([^\s>]+)>?\s*("([^"]*)")?/,P=/^!\[([^\]]*)\] ?\[([^\]]*)\]/,_=/^\[([^\]]*)\] ?\[([^\]]*)\]/,F=/(\n|^[-*]\s|^#|^ {2,}|^-{2,}|^>\s)/,G=/\t/g,W=/(^ *\||\| *$)/g,Z=/^ *:-+: *$/,q=/^ *:-+ *$/,Q=/^ *-+: *$/,V="((?:\\[.*?\\][([].*?[)\\]]|<.*?>(?:.*?<.*?>)?|`.*?`|\\\\\\1|[\\s\\S])+?)",X=new RegExp(`^([*_])\\1${V}\\1\\1(?!\\1)`),J=new RegExp(`^([*_])${V}\\1(?!\\1)`),K=new RegExp(`^(==)${V}\\1`),Y=new RegExp(`^(~~)${V}\\1`),ee=/^\\([^0-9A-Za-z\s])/,te=/\\([^0-9A-Za-z\s])/g,ne=/^([\s\S](?:(?!  |[0-9]\.)[^*_~\-\n<`\\\[!])*)/,re=/^\n+/,ie=/^([ \t]*)/,le=/\\([^\\])/g,ae=/(?:^|\n)( *)$/,oe="(?:\\d+\\.)",ce="(?:[*+-])";function se(e){return"( *)("+(1===e?oe:ce)+") +"}const de=se(1),ue=se(2);function pe(e){return new RegExp("^"+(1===e?de:ue))}const fe=pe(1),he=pe(2);function me(e){return new RegExp("^"+(1===e?de:ue)+"[^\\n]*(?:\\n(?!\\1"+(1===e?oe:ce)+" )[^\\n]*)*(\\n|$)","gm")}const ge=me(1),ye=me(2);function ke(e){const t=1===e?oe:ce;return new RegExp("^( *)("+t+") [\\s\\S]+?(?:\\n{2,}(?! )(?!\\1"+t+" (?!"+t+" ))\\n*|\\s*\\n*$)")}const xe=ke(1),be=ke(2);function ve(e,t){const n=1===t,i=n?xe:be,l=n?ge:ye,a=n?fe:he;return{match:Oe(function(e,t){const n=ae.exec(t.prevCapture);return n&&(t.list||!t.inline&&!t.simple)?i.exec(e=n[1]+e):null}),order:1,parse(e,t,r){const i=n?+e[2]:void 0,o=e[0].replace(d,"\n").match(l);let c=!1;return{items:o.map(function(e,n){const i=a.exec(e)[0].length,l=new RegExp("^ {1,"+i+"}","gm"),s=e.replace(l,"").replace(a,""),d=n===o.length-1,u=-1!==s.indexOf("\n\n")||d&&c;c=u;const p=r.inline,f=r.list;let h;r.list=!0,u?(r.inline=!1,h=Ee(s)+"\n\n"):(r.inline=!0,h=Ee(s));const m=t(h,r);return r.inline=p,r.list=f,m}),ordered:n,start:i}},render:(t,n,i)=>e(t.ordered?"ol":"ul",{key:i.key,start:t.type===r.orderedList?t.start:void 0},t.items.map(function(t,r){return e("li",{key:r},n(t,i))}))}}const Ce=new RegExp("^\\[((?:\\[[^\\]]*\\]|[^\\[\\]]|\\](?=[^\\[]*\\]))*)\\]\\(\\s*<?((?:\\([^)]*\\)|[^\\s\\\\]|\\\\.)*?)>?(?:\\s+['\"]([\\s\\S]*?)['\"])?\\s*\\)"),$e=/^!\[(.*?)\]\( *((?:\([^)]*\)|[^() ])*) *"?([^)"]*)?"?\)/,Se=[u,g,y,E,L,z,j,xe,be],we=[...Se,/^[^\n]+(?:  \n|\n{2,})/,A,B,M];function Ee(e){let t=e.length;for(;t>0&&e[t-1]<=" ";)t--;return e.slice(0,t)}function ze(e){return e.replace(/[ÀÁÂÃÄÅàáâãäåæÆ]/g,"a").replace(/[çÇ]/g,"c").replace(/[ðÐ]/g,"d").replace(/[ÈÉÊËéèêë]/g,"e").replace(/[ÏïÎîÍíÌì]/g,"i").replace(/[Ññ]/g,"n").replace(/[øØœŒÕõÔôÓóÒò]/g,"o").replace(/[ÜüÛûÚúÙù]/g,"u").replace(/[ŸÿÝý]/g,"y").replace(/[^a-z0-9- ]/gi,"").replace(/ /gi,"-").toLowerCase()}function Le(e){return Q.test(e)?"right":Z.test(e)?"center":q.test(e)?"left":null}function Ae(e,t,n,r){const i=n.inTable;n.inTable=!0;let l=[[]],a="";function o(){if(!a)return;const e=l[l.length-1];e.push.apply(e,t(a,n)),a=""}return e.trim().split(/(`[^`]*`|\\\||\|)/).filter(Boolean).forEach((e,t,n)=>{"|"===e.trim()&&(o(),r)?0!==t&&t!==n.length-1&&l.push([]):a+=e}),o(),n.inTable=i,l}function Te(e,t,n){n.inline=!0;const i=e[2]?e[2].replace(W,"").split("|").map(Le):[],l=e[3]?function(e,t,n){return e.trim().split("\n").map(function(e){return Ae(e,t,n,!0)})}(e[3],t,n):[],a=Ae(e[1],t,n,!!l.length);return n.inline=!1,l.length?{align:i,cells:l,header:a,type:r.table}:{children:a,type:r.paragraph}}function Be(e,t){return null==e.align[t]?{}:{textAlign:e.align[t]}}function Oe(e){return e.inline=1,e}function Me(e){return Oe(function(t,n){return n.inline?e.exec(t):null})}function Re(e){return Oe(function(t,n){return n.inline||n.simple?e.exec(t):null})}function Ie(e){return function(t,n){return n.inline||n.simple?null:e.exec(t)}}function Ue(e){return Oe(function(t){return e.exec(t)})}function De(e,t){if(t.inline||t.simple)return null;let n="";e.split("\n").every(e=>(e+="\n",!Se.some(t=>t.test(e))&&(n+=e,!!e.trim())));const r=Ee(n);return""==r?null:[n,,r]}function Ne(e){try{if(decodeURIComponent(e).replace(/[^A-Za-z0-9/:]/g,"").match(/^\s*(javascript|vbscript|data(?!:image)):/i))return null}catch(e){return null}return e}function je(e){return e.replace(le,"$1")}function He(e,t,n){const r=n.inline||!1,i=n.simple||!1;n.inline=!0,n.simple=!0;const l=e(t,n);return n.inline=r,n.simple=i,l}function Pe(e,t,n){const r=n.inline||!1,i=n.simple||!1;n.inline=!1,n.simple=!0;const l=e(t,n);return n.inline=r,n.simple=i,l}function _e(e,t,n){const r=n.inline||!1;n.inline=!1;const i=e(t,n);return n.inline=r,i}const Fe=(e,t,n)=>({children:He(t,e[2],n)});function Ge(){return{}}function We(){return null}function Ze(...e){return e.filter(Boolean).join(" ")}function qe(e,t,n){let r=e;const i=t.split(".");for(;i.length&&(r=r[i[0]],void 0!==r);)i.shift();return r||n}function Qe(n="",i={}){function d(e,n,...r){const l=qe(i.overrides,`${e}.props`,{});return i.createElement(function(e,t){const n=qe(t,e);return n?"function"==typeof n||"object"==typeof n&&"render"in n?n:qe(t,`${e}.component`,e):e}(e,i.overrides),t({},n,l,{className:Ze(null==n?void 0:n.className,l.className)||void 0}),...r)}function W(e){e=e.replace(S,"");let t=!1;i.forceInline?t=!0:i.forceBlock||(t=!1===F.test(e));const n=ae(le(t?e:`${Ee(e).replace(re,"")}\n\n`,{inline:t}));for(;"string"==typeof n[n.length-1]&&!n[n.length-1].trim();)n.pop();if(null===i.wrapper)return n;const r=i.wrapper||(t?"span":"div");let l;if(n.length>1||i.forceWrapper)l=n;else{if(1===n.length)return l=n[0],"string"==typeof l?d("span",{key:"outer"},l):l;l=null}return i.createElement(r,{key:"outer"},l)}function Z(e,t){const n=t.match(c);return n?n.reduce(function(t,n){const r=n.indexOf("=");if(-1!==r){const a=function(e){return-1!==e.indexOf("-")&&null===e.match(O)&&(e=e.replace(N,function(e,t){return t.toUpperCase()})),e}(n.slice(0,r)).trim(),o=function(e){const t=e[0];return('"'===t||"'"===t)&&e.length>=2&&e[e.length-1]===t?e.slice(1,-1):e}(n.slice(r+1).trim()),c=l[a]||a;if("ref"===c)return t;const s=t[c]=function(e,t,n,r){return"style"===t?n.split(/;\s?/).reduce(function(e,t){const n=t.slice(0,t.indexOf(":"));return e[n.trim().replace(/(-[a-z])/g,e=>e[1].toUpperCase())]=t.slice(n.length+1).trim(),e},{}):"href"===t||"src"===t?r(n,e,t):(n.match(R)&&(n=n.slice(1,n.length-1)),"true"===n||"false"!==n&&n)}(e,a,o,i.sanitizer);"string"==typeof s&&(A.test(s)||M.test(s))&&(t[c]=W(s.trim()))}else"style"!==n&&(t[l[n]||n]=!0);return t},{}):null}i.overrides=i.overrides||{},i.sanitizer=i.sanitizer||Ne,i.slugify=i.slugify||ze,i.namedCodesToUnicode=i.namedCodesToUnicode?t({},a,i.namedCodesToUnicode):a,i.createElement=i.createElement||index_js_.createElement;const q=[],Q={},V={[r.blockQuote]:{match:Ie(u),order:1,parse(e,t,n){const[,r,i]=e[0].replace(p,"").match(f);return{alert:r,children:t(i,n)}},render(e,t,n){const l={key:n.key};return e.alert&&(l.className="markdown-alert-"+i.slugify(e.alert.toLowerCase(),ze),e.children.unshift({attrs:{},children:[{type:r.text,text:e.alert}],noInnerParse:!0,type:r.htmlBlock,tag:"header"})),d("blockquote",l,t(e.children,n))}},[r.breakLine]:{match:Ue(h),order:1,parse:Ge,render:(e,t,n)=>d("br",{key:n.key})},[r.breakThematic]:{match:Ie(m),order:1,parse:Ge,render:(e,t,n)=>d("hr",{key:n.key})},[r.codeBlock]:{match:Ie(y),order:0,parse:e=>({lang:void 0,text:Ee(e[0].replace(/^ {4}/gm,"")).replace(te,"$1")}),render:(e,n,r)=>d("pre",{key:r.key},d("code",t({},e.attrs,{className:e.lang?`lang-${e.lang}`:""}),e.text))},[r.codeFenced]:{match:Ie(g),order:0,parse:e=>({attrs:Z("code",e[3]||""),lang:e[2]||void 0,text:e[4].replace(te,"$1"),type:r.codeBlock})},[r.codeInline]:{match:Re(k),order:3,parse:e=>({text:e[2].replace(te,"$1")}),render:(e,t,n)=>d("code",{key:n.key},e.text)},[r.footnote]:{match:Ie(v),order:0,parse:e=>(q.push({footnote:e[2],identifier:e[1]}),{}),render:We},[r.footnoteReference]:{match:Me(C),order:1,parse:e=>({target:`#${i.slugify(e[1],ze)}`,text:e[1]}),render:(e,t,n)=>d("a",{key:n.key,href:i.sanitizer(e.target,"a","href")},d("sup",{key:n.key},e.text))},[r.gfmTask]:{match:Me(w),order:1,parse:e=>({completed:"x"===e[1].toLowerCase()}),render:(e,t,n)=>d("input",{checked:e.completed,key:n.key,readOnly:!0,type:"checkbox"})},[r.heading]:{match:Ie(i.enforceAtxHeadings?z:E),order:1,parse:(e,t,n)=>({children:He(t,e[2],n),id:i.slugify(e[2],ze),level:e[1].length}),render:(e,t,n)=>d(`h${e.level}`,{id:e.id,key:n.key},t(e.children,n))},[r.headingSetext]:{match:Ie(L),order:0,parse:(e,t,n)=>({children:He(t,e[1],n),level:"="===e[2]?1:2,type:r.heading})},[r.htmlBlock]:{match:Ue(A),order:1,parse(e,t,n){const[,r]=e[3].match(ie),i=new RegExp(`^${r}`,"gm"),l=e[3].replace(i,""),a=(c=l,we.some(e=>e.test(c))?_e:He);var c;const s=e[1].toLowerCase(),d=-1!==o.indexOf(s),u=(d?s:e[1]).trim(),p={attrs:Z(u,e[2]),noInnerParse:d,tag:u};return n.inAnchor=n.inAnchor||"a"===s,d?p.text=e[3]:p.children=a(t,l,n),n.inAnchor=!1,p},render:(e,n,r)=>d(e.tag,t({key:r.key},e.attrs),e.text||(e.children?n(e.children,r):""))},[r.htmlSelfClosing]:{match:Ue(M),order:1,parse(e){const t=e[1].trim();return{attrs:Z(t,e[2]||""),tag:t}},render:(e,n,r)=>d(e.tag,t({},e.attrs,{key:r.key}))},[r.htmlComment]:{match:Ue(B),order:1,parse:()=>({}),render:We},[r.image]:{match:Re($e),order:1,parse:e=>({alt:e[1],target:je(e[2]),title:e[3]}),render:(e,t,n)=>d("img",{key:n.key,alt:e.alt||void 0,title:e.title||void 0,src:i.sanitizer(e.target,"img","src")})},[r.link]:{match:Me(Ce),order:3,parse:(e,t,n)=>({children:Pe(t,e[1],n),target:je(e[2]),title:e[3]}),render:(e,t,n)=>d("a",{key:n.key,href:i.sanitizer(e.target,"a","href"),title:e.title},t(e.children,n))},[r.linkAngleBraceStyleDetector]:{match:Me(D),order:0,parse:e=>({children:[{text:e[1],type:r.text}],target:e[1],type:r.link})},[r.linkBareUrlDetector]:{match:Oe((e,t)=>t.inAnchor||i.disableAutoLink?null:Me(I)(e,t)),order:0,parse:e=>({children:[{text:e[1],type:r.text}],target:e[1],title:void 0,type:r.link})},[r.linkMailtoDetector]:{match:Me(U),order:0,parse(e){let t=e[1],n=e[1];return s.test(n)||(n="mailto:"+n),{children:[{text:t.replace("mailto:",""),type:r.text}],target:n,type:r.link}}},[r.orderedList]:ve(d,1),[r.unorderedList]:ve(d,2),[r.newlineCoalescer]:{match:Ie(x),order:3,parse:Ge,render:()=>"\n"},[r.paragraph]:{match:Oe(De),order:3,parse:Fe,render:(e,t,n)=>d("p",{key:n.key},t(e.children,n))},[r.ref]:{match:Me(H),order:0,parse:e=>(Q[e[1]]={target:e[2],title:e[4]},{}),render:We},[r.refImage]:{match:Re(P),order:0,parse:e=>({alt:e[1]||void 0,ref:e[2]}),render:(e,t,n)=>Q[e.ref]?d("img",{key:n.key,alt:e.alt,src:i.sanitizer(Q[e.ref].target,"img","src"),title:Q[e.ref].title}):null},[r.refLink]:{match:Me(_),order:0,parse:(e,t,n)=>({children:t(e[1],n),fallbackChildren:e[0],ref:e[2]}),render:(e,t,n)=>Q[e.ref]?d("a",{key:n.key,href:i.sanitizer(Q[e.ref].target,"a","href"),title:Q[e.ref].title},t(e.children,n)):d("span",{key:n.key},e.fallbackChildren)},[r.table]:{match:Ie(j),order:1,parse:Te,render(e,t,n){const r=e;return d("table",{key:n.key},d("thead",null,d("tr",null,r.header.map(function(e,i){return d("th",{key:i,style:Be(r,i)},t(e,n))}))),d("tbody",null,r.cells.map(function(e,i){return d("tr",{key:i},e.map(function(e,i){return d("td",{key:i,style:Be(r,i)},t(e,n))}))})))}},[r.text]:{match:Ue(ne),order:4,parse:e=>({text:e[0].replace(T,(e,t)=>i.namedCodesToUnicode[t]?i.namedCodesToUnicode[t]:e)}),render:e=>e.text},[r.textBolded]:{match:Re(X),order:2,parse:(e,t,n)=>({children:t(e[2],n)}),render:(e,t,n)=>d("strong",{key:n.key},t(e.children,n))},[r.textEmphasized]:{match:Re(J),order:3,parse:(e,t,n)=>({children:t(e[2],n)}),render:(e,t,n)=>d("em",{key:n.key},t(e.children,n))},[r.textEscaped]:{match:Re(ee),order:1,parse:e=>({text:e[1],type:r.text})},[r.textMarked]:{match:Re(K),order:3,parse:Fe,render:(e,t,n)=>d("mark",{key:n.key},t(e.children,n))},[r.textStrikethroughed]:{match:Re(Y),order:3,parse:Fe,render:(e,t,n)=>d("del",{key:n.key},t(e.children,n))}};!0===i.disableParsingRawHTML&&(delete V[r.htmlBlock],delete V[r.htmlSelfClosing]);const le=function(e){let t=Object.keys(e);function n(r,i){let l,a,o=[],c="",s="";for(i.prevCapture=i.prevCapture||"";r;){let d=0;for(;d<t.length;){if(c=t[d],l=e[c],i.inline&&!l.match.inline){d++;continue}const u=l.match(r,i);if(u){s=u[0],i.prevCapture+=s,r=r.substring(s.length),a=l.parse(u,n,i),null==a.type&&(a.type=c),o.push(a);break}d++}}return i.prevCapture="",o}return t.sort(function(t,n){let r=e[t].order,i=e[n].order;return r!==i?r-i:t<n?-1:1}),function(e,t){return n(function(e){return e.replace(b,"\n").replace($,"").replace(G,"    ")}(e),t)}}(V),ae=(oe=function(e,t){return function(n,r,i){const l=e[n.type].render;return t?t(()=>l(n,r,i),n,r,i):l(n,r,i)}}(V,i.renderRule),function e(t,n={}){if(Array.isArray(t)){const r=n.key,i=[];let l=!1;for(let r=0;r<t.length;r++){n.key=r;const a=e(t[r],n),o="string"==typeof a;o&&l?i[i.length-1]+=a:null!==a&&i.push(a),l=o}return n.key=r,i}return oe(t,e,n)});var oe;const ce=W(n);return q.length?d("div",null,ce,d("footer",{key:"footer"},q.map(function(e){return d("div",{id:i.slugify(e.identifier,ze),key:e.identifier},e.identifier,ae(le(e.footnote,{inline:!0})))}))):ce}/* harmony default export */ const index_modern = (t=>{let{children:r="",options:i}=t,l=function(e,t){if(null==e)return{};var n,r,i={},l=Object.keys(e);for(r=0;r<l.length;r++)t.indexOf(n=l[r])>=0||(i[n]=e[n]);return i}(t,n);return index_js_.cloneElement(Qe(r,i),l)});
//# sourceMappingURL=index.modern.js.map

// EXTERNAL MODULE: ../node_modules/lodash/has.js
var has = __webpack_require__(73915);
var has_default = /*#__PURE__*/__webpack_require__.n(has);
// EXTERNAL MODULE: ../node_modules/lodash/unset.js
var unset = __webpack_require__(43551);
var unset_default = /*#__PURE__*/__webpack_require__.n(unset);
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/ObjectField.js









/** The `ObjectField` component is used to render a field in the schema that is of type `object`. It tracks whether an
 * additional property key was modified and what it was modified to
 *
 * @param props - The `FieldProps` for this template
 */
class ObjectField extends index_js_.Component {
    constructor() {
        super(...arguments);
        /** Set up the initial state */
        this.state = {
            wasPropertyKeyModified: false,
            additionalProperties: {},
        };
        /** Returns the `onPropertyChange` handler for the `name` field. Handles the special case where a user is attempting
         * to clear the data for a field added as an additional property. Calls the `onChange()` handler with the updated
         * formData.
         *
         * @param name - The name of the property
         * @param addedByAdditionalProperties - Flag indicating whether this property is an additional property
         * @returns - The onPropertyChange callback for the `name` property
         */
        this.onPropertyChange = (name, addedByAdditionalProperties = false) => {
            return (value, newErrorSchema, id) => {
                const { formData, onChange, errorSchema } = this.props;
                if (value === undefined && addedByAdditionalProperties) {
                    // Don't set value = undefined for fields added by
                    // additionalProperties. Doing so removes them from the
                    // formData, which causes them to completely disappear
                    // (including the input field for the property name). Unlike
                    // fields which are "mandated" by the schema, these fields can
                    // be set to undefined by clicking a "delete field" button, so
                    // set empty values to the empty string.
                    value = '';
                }
                const newFormData = { ...formData, [name]: value };
                onChange(newFormData, errorSchema &&
                    errorSchema && {
                    ...errorSchema,
                    [name]: newErrorSchema,
                }, id);
            };
        };
        /** Returns a callback to handle the onDropPropertyClick event for the given `key` which removes the old `key` data
         * and calls the `onChange` callback with it
         *
         * @param key - The key for which the drop callback is desired
         * @returns - The drop property click callback
         */
        this.onDropPropertyClick = (key) => {
            return (event) => {
                event.preventDefault();
                const { onChange, formData } = this.props;
                const copiedFormData = { ...formData };
                unset_default()(copiedFormData, key);
                onChange(copiedFormData);
            };
        };
        /** Computes the next available key name from the `preferredKey`, indexing through the already existing keys until one
         * that is already not assigned is found.
         *
         * @param preferredKey - The preferred name of a new key
         * @param [formData] - The form data in which to check if the desired key already exists
         * @returns - The name of the next available key from `preferredKey`
         */
        this.getAvailableKey = (preferredKey, formData) => {
            const { uiSchema, registry } = this.props;
            const { duplicateKeySuffixSeparator = '-' } = (0,lib_index_js_.getUiOptions)(uiSchema, registry.globalUiOptions);
            let index = 0;
            let newKey = preferredKey;
            while (has_default()(formData, newKey)) {
                newKey = `${preferredKey}${duplicateKeySuffixSeparator}${++index}`;
            }
            return newKey;
        };
        /** Returns a callback function that deals with the rename of a key for an additional property for a schema. That
         * callback will attempt to rename the key and move the existing data to that key, calling `onChange` when it does.
         *
         * @param oldValue - The old value of a field
         * @returns - The key change callback function
         */
        this.onKeyChange = (oldValue) => {
            return (value, newErrorSchema) => {
                if (oldValue === value) {
                    return;
                }
                const { formData, onChange, errorSchema } = this.props;
                value = this.getAvailableKey(value, formData);
                const newFormData = {
                    ...formData,
                };
                const newKeys = { [oldValue]: value };
                const keyValues = Object.keys(newFormData).map((key) => {
                    const newKey = newKeys[key] || key;
                    return { [newKey]: newFormData[key] };
                });
                const renamedObj = Object.assign({}, ...keyValues);
                this.setState({ wasPropertyKeyModified: true });
                onChange(renamedObj, errorSchema &&
                    errorSchema && {
                    ...errorSchema,
                    [value]: newErrorSchema,
                });
            };
        };
        /** Handles the adding of a new additional property on the given `schema`. Calls the `onChange` callback once the new
         * default data for that field has been added to the formData.
         *
         * @param schema - The schema element to which the new property is being added
         */
        this.handleAddClick = (schema) => () => {
            if (!schema.additionalProperties) {
                return;
            }
            const { formData, onChange, registry } = this.props;
            const newFormData = { ...formData };
            let type = undefined;
            if (isObject_default()(schema.additionalProperties)) {
                type = schema.additionalProperties.type;
                let apSchema = schema.additionalProperties;
                if (lib_index_js_.REF_KEY in apSchema) {
                    const { schemaUtils } = registry;
                    apSchema = schemaUtils.retrieveSchema({ $ref: apSchema[lib_index_js_.REF_KEY] }, formData);
                    type = apSchema.type;
                }
                if (!type && (lib_index_js_.ANY_OF_KEY in apSchema || lib_index_js_.ONE_OF_KEY in apSchema)) {
                    type = 'object';
                }
            }
            const newKey = this.getAvailableKey('newKey', newFormData);
            // Cast this to make the `set` work properly
            set_default()(newFormData, newKey, this.getDefaultValue(type));
            onChange(newFormData);
        };
    }
    /** Returns a flag indicating whether the `name` field is required in the object schema
     *
     * @param name - The name of the field to check for required-ness
     * @returns - True if the field `name` is required, false otherwise
     */
    isRequired(name) {
        const { schema } = this.props;
        return Array.isArray(schema.required) && schema.required.indexOf(name) !== -1;
    }
    /** Returns a default value to be used for a new additional schema property of the given `type`
     *
     * @param type - The type of the new additional schema property
     */
    getDefaultValue(type) {
        const { registry: { translateString }, } = this.props;
        switch (type) {
            case 'array':
                return [];
            case 'boolean':
                return false;
            case 'null':
                return null;
            case 'number':
                return 0;
            case 'object':
                return {};
            case 'string':
            default:
                // We don't have a datatype for some reason (perhaps additionalProperties was true)
                return translateString(lib_index_js_.TranslatableString.NewStringDefault);
        }
    }
    /** Renders the `ObjectField` from the given props
     */
    render() {
        var _a, _b, _c;
        const { schema: rawSchema, uiSchema = {}, formData, errorSchema, idSchema, name, required = false, disabled = false, readonly = false, hideError, idPrefix, idSeparator, onBlur, onFocus, registry, } = this.props;
        const { fields, formContext, schemaUtils, translateString, globalUiOptions } = registry;
        const { SchemaField } = fields;
        const schema = schemaUtils.retrieveSchema(rawSchema, formData);
        const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
        const { properties: schemaProperties = {} } = schema;
        const title = (_b = (_a = uiOptions.title) !== null && _a !== void 0 ? _a : schema.title) !== null && _b !== void 0 ? _b : name;
        const description = (_c = uiOptions.description) !== null && _c !== void 0 ? _c : schema.description;
        let orderedProperties;
        try {
            const properties = Object.keys(schemaProperties);
            orderedProperties = (0,lib_index_js_.orderProperties)(properties, uiOptions.order);
        }
        catch (err) {
            return ((0,jsx_runtime.jsxs)("div", { children: [(0,jsx_runtime.jsx)("p", { className: 'config-error', style: { color: 'red' }, children: (0,jsx_runtime.jsx)(index_modern, { children: translateString(lib_index_js_.TranslatableString.InvalidObjectField, [name || 'root', err.message]) }) }), (0,jsx_runtime.jsx)("pre", { children: JSON.stringify(schema) })] }));
        }
        const Template = (0,lib_index_js_.getTemplate)('ObjectFieldTemplate', registry, uiOptions);
        const templateProps = {
            // getDisplayLabel() always returns false for object types, so just check the `uiOptions.label`
            title: uiOptions.label === false ? '' : title,
            description: uiOptions.label === false ? undefined : description,
            properties: orderedProperties.map((name) => {
                const addedByAdditionalProperties = has_default()(schema, [lib_index_js_.PROPERTIES_KEY, name, lib_index_js_.ADDITIONAL_PROPERTY_FLAG]);
                const fieldUiSchema = addedByAdditionalProperties ? uiSchema.additionalProperties : uiSchema[name];
                const hidden = (0,lib_index_js_.getUiOptions)(fieldUiSchema).widget === 'hidden';
                const fieldIdSchema = get_default()(idSchema, [name], {});
                return {
                    content: ((0,jsx_runtime.jsx)(SchemaField, { name: name, required: this.isRequired(name), schema: get_default()(schema, [lib_index_js_.PROPERTIES_KEY, name], {}), uiSchema: fieldUiSchema, errorSchema: get_default()(errorSchema, name), idSchema: fieldIdSchema, idPrefix: idPrefix, idSeparator: idSeparator, formData: get_default()(formData, name), formContext: formContext, wasPropertyKeyModified: this.state.wasPropertyKeyModified, onKeyChange: this.onKeyChange(name), onChange: this.onPropertyChange(name, addedByAdditionalProperties), onBlur: onBlur, onFocus: onFocus, registry: registry, disabled: disabled, readonly: readonly, hideError: hideError, onDropPropertyClick: this.onDropPropertyClick }, name)),
                    name,
                    readonly,
                    disabled,
                    required,
                    hidden,
                };
            }),
            readonly,
            disabled,
            required,
            idSchema,
            uiSchema,
            errorSchema,
            schema,
            formData,
            formContext,
            registry,
        };
        return (0,jsx_runtime.jsx)(Template, { ...templateProps, onAddClick: this.handleAddClick });
    }
}
/* harmony default export */ const fields_ObjectField = (ObjectField);
//# sourceMappingURL=ObjectField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/SchemaField.js






/** The map of component type to FieldName */
const COMPONENT_TYPES = {
    array: 'ArrayField',
    boolean: 'BooleanField',
    integer: 'NumberField',
    number: 'NumberField',
    object: 'ObjectField',
    string: 'StringField',
    null: 'NullField',
};
/** Computes and returns which `Field` implementation to return in order to render the field represented by the
 * `schema`. The `uiOptions` are used to alter what potential `Field` implementation is actually returned. If no
 * appropriate `Field` implementation can be found then a wrapper around `UnsupportedFieldTemplate` is used.
 *
 * @param schema - The schema from which to obtain the type
 * @param uiOptions - The UI Options that may affect the component decision
 * @param idSchema - The id that is passed to the `UnsupportedFieldTemplate`
 * @param registry - The registry from which fields and templates are obtained
 * @returns - The `Field` component that is used to render the actual field data
 */
function getFieldComponent(schema, uiOptions, idSchema, registry) {
    const field = uiOptions.field;
    const { fields, translateString } = registry;
    if (typeof field === 'function') {
        return field;
    }
    if (typeof field === 'string' && field in fields) {
        return fields[field];
    }
    const schemaType = (0,lib_index_js_.getSchemaType)(schema);
    const type = Array.isArray(schemaType) ? schemaType[0] : schemaType || '';
    const schemaId = schema.$id;
    let componentName = COMPONENT_TYPES[type];
    if (schemaId && schemaId in fields) {
        componentName = schemaId;
    }
    // If the type is not defined and the schema uses 'anyOf' or 'oneOf', don't
    // render a field and let the MultiSchemaField component handle the form display
    if (!componentName && (schema.anyOf || schema.oneOf)) {
        return () => null;
    }
    return componentName in fields
        ? fields[componentName]
        : () => {
            const UnsupportedFieldTemplate = (0,lib_index_js_.getTemplate)('UnsupportedFieldTemplate', registry, uiOptions);
            return ((0,jsx_runtime.jsx)(UnsupportedFieldTemplate, { schema: schema, idSchema: idSchema, reason: translateString(lib_index_js_.TranslatableString.UnknownFieldType, [String(schema.type)]), registry: registry }));
        };
}
/** The `SchemaFieldRender` component is the work-horse of react-jsonschema-form, determining what kind of real field to
 * render based on the `schema`, `uiSchema` and all the other props. It also deals with rendering the `anyOf` and
 * `oneOf` fields.
 *
 * @param props - The `FieldProps` for this component
 */
function SchemaFieldRender(props) {
    const { schema: _schema, idSchema: _idSchema, uiSchema, formData, errorSchema, idPrefix, idSeparator, name, onChange, onKeyChange, onDropPropertyClick, required, registry, wasPropertyKeyModified = false, } = props;
    const { formContext, schemaUtils, globalUiOptions } = registry;
    const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema, globalUiOptions);
    const FieldTemplate = (0,lib_index_js_.getTemplate)('FieldTemplate', registry, uiOptions);
    const DescriptionFieldTemplate = (0,lib_index_js_.getTemplate)('DescriptionFieldTemplate', registry, uiOptions);
    const FieldHelpTemplate = (0,lib_index_js_.getTemplate)('FieldHelpTemplate', registry, uiOptions);
    const FieldErrorTemplate = (0,lib_index_js_.getTemplate)('FieldErrorTemplate', registry, uiOptions);
    const schema = schemaUtils.retrieveSchema(_schema, formData);
    const fieldId = _idSchema[lib_index_js_.ID_KEY];
    const idSchema = (0,lib_index_js_.mergeObjects)(schemaUtils.toIdSchema(schema, fieldId, formData, idPrefix, idSeparator), _idSchema);
    /** Intermediary `onChange` handler for field components that will inject the `id` of the current field into the
     * `onChange` chain if it is not already being provided from a deeper level in the hierarchy
     */
    const handleFieldComponentChange = (0,index_js_.useCallback)((formData, newErrorSchema, id) => {
        const theId = id || fieldId;
        return onChange(formData, newErrorSchema, theId);
    }, [fieldId, onChange]);
    const FieldComponent = getFieldComponent(schema, uiOptions, idSchema, registry);
    const disabled = Boolean(props.disabled || uiOptions.disabled);
    const readonly = Boolean(props.readonly || uiOptions.readonly || props.schema.readOnly || schema.readOnly);
    const uiSchemaHideError = uiOptions.hideError;
    // Set hideError to the value provided in the uiSchema, otherwise stick with the prop to propagate to children
    const hideError = uiSchemaHideError === undefined ? props.hideError : Boolean(uiSchemaHideError);
    const autofocus = Boolean(props.autofocus || uiOptions.autofocus);
    if (Object.keys(schema).length === 0) {
        return null;
    }
    const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
    const { __errors, ...fieldErrorSchema } = errorSchema || {};
    // See #439: uiSchema: Don't pass consumed class names or style to child components
    const fieldUiSchema = omit_default()(uiSchema, ['ui:classNames', 'classNames', 'ui:style']);
    if (lib_index_js_.UI_OPTIONS_KEY in fieldUiSchema) {
        fieldUiSchema[lib_index_js_.UI_OPTIONS_KEY] = omit_default()(fieldUiSchema[lib_index_js_.UI_OPTIONS_KEY], ['classNames', 'style']);
    }
    const field = ((0,jsx_runtime.jsx)(FieldComponent, { ...props, onChange: handleFieldComponentChange, idSchema: idSchema, schema: schema, uiSchema: fieldUiSchema, disabled: disabled, readonly: readonly, hideError: hideError, autofocus: autofocus, errorSchema: fieldErrorSchema, formContext: formContext, rawErrors: __errors }));
    const id = idSchema[lib_index_js_.ID_KEY];
    // If this schema has a title defined, but the user has set a new key/label, retain their input.
    let label;
    if (wasPropertyKeyModified) {
        label = name;
    }
    else {
        label = lib_index_js_.ADDITIONAL_PROPERTY_FLAG in schema ? name : uiOptions.title || props.schema.title || schema.title || name;
    }
    const description = uiOptions.description || props.schema.description || schema.description || '';
    const richDescription = uiOptions.enableMarkdownInDescription ? (0,jsx_runtime.jsx)(index_modern, { children: description }) : description;
    const help = uiOptions.help;
    const hidden = uiOptions.widget === 'hidden';
    const classNames = ['form-group', 'field', `field-${(0,lib_index_js_.getSchemaType)(schema)}`];
    if (!hideError && __errors && __errors.length > 0) {
        classNames.push('field-error has-error has-danger');
    }
    if (uiSchema === null || uiSchema === void 0 ? void 0 : uiSchema.classNames) {
        if (false) {}
        classNames.push(uiSchema.classNames);
    }
    if (uiOptions.classNames) {
        classNames.push(uiOptions.classNames);
    }
    const helpComponent = ((0,jsx_runtime.jsx)(FieldHelpTemplate, { help: help, idSchema: idSchema, schema: schema, uiSchema: uiSchema, hasErrors: !hideError && __errors && __errors.length > 0, registry: registry }));
    /*
     * AnyOf/OneOf errors handled by child schema
     * unless it can be rendered as select control
     */
    const errorsComponent = hideError || ((schema.anyOf || schema.oneOf) && !schemaUtils.isSelect(schema)) ? undefined : ((0,jsx_runtime.jsx)(FieldErrorTemplate, { errors: __errors, errorSchema: errorSchema, idSchema: idSchema, schema: schema, uiSchema: uiSchema, registry: registry }));
    const fieldProps = {
        description: ((0,jsx_runtime.jsx)(DescriptionFieldTemplate, { id: (0,lib_index_js_.descriptionId)(id), description: richDescription, schema: schema, uiSchema: uiSchema, registry: registry })),
        rawDescription: description,
        help: helpComponent,
        rawHelp: typeof help === 'string' ? help : undefined,
        errors: errorsComponent,
        rawErrors: hideError ? undefined : __errors,
        id,
        label,
        hidden,
        onChange,
        onKeyChange,
        onDropPropertyClick,
        required,
        disabled,
        readonly,
        hideError,
        displayLabel,
        classNames: classNames.join(' ').trim(),
        style: uiOptions.style,
        formContext,
        formData,
        schema,
        uiSchema,
        registry,
    };
    const _AnyOfField = registry.fields.AnyOfField;
    const _OneOfField = registry.fields.OneOfField;
    const isReplacingAnyOrOneOf = (uiSchema === null || uiSchema === void 0 ? void 0 : uiSchema['ui:field']) && (uiSchema === null || uiSchema === void 0 ? void 0 : uiSchema['ui:fieldReplacesAnyOrOneOf']) === true;
    return ((0,jsx_runtime.jsx)(FieldTemplate, { ...fieldProps, children: (0,jsx_runtime.jsxs)(jsx_runtime.Fragment, { children: [field, schema.anyOf && !isReplacingAnyOrOneOf && !schemaUtils.isSelect(schema) && ((0,jsx_runtime.jsx)(_AnyOfField, { name: name, disabled: disabled, readonly: readonly, hideError: hideError, errorSchema: errorSchema, formData: formData, formContext: formContext, idPrefix: idPrefix, idSchema: idSchema, idSeparator: idSeparator, onBlur: props.onBlur, onChange: props.onChange, onFocus: props.onFocus, options: schema.anyOf.map((_schema) => schemaUtils.retrieveSchema(isObject_default()(_schema) ? _schema : {}, formData)), registry: registry, schema: schema, uiSchema: uiSchema })), schema.oneOf && !isReplacingAnyOrOneOf && !schemaUtils.isSelect(schema) && ((0,jsx_runtime.jsx)(_OneOfField, { name: name, disabled: disabled, readonly: readonly, hideError: hideError, errorSchema: errorSchema, formData: formData, formContext: formContext, idPrefix: idPrefix, idSchema: idSchema, idSeparator: idSeparator, onBlur: props.onBlur, onChange: props.onChange, onFocus: props.onFocus, options: schema.oneOf.map((_schema) => schemaUtils.retrieveSchema(isObject_default()(_schema) ? _schema : {}, formData)), registry: registry, schema: schema, uiSchema: uiSchema }))] }) }));
}
/** The `SchemaField` component determines whether it is necessary to rerender the component based on any props changes
 * and if so, calls the `SchemaFieldRender` component with the props.
 */
class SchemaField extends index_js_.Component {
    shouldComponentUpdate(nextProps) {
        return !(0,lib_index_js_.deepEquals)(this.props, nextProps);
    }
    render() {
        return (0,jsx_runtime.jsx)(SchemaFieldRender, { ...this.props });
    }
}
/* harmony default export */ const fields_SchemaField = (SchemaField);
//# sourceMappingURL=SchemaField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/StringField.js


/** The `StringField` component is used to render a schema field that represents a string type
 *
 * @param props - The `FieldProps` for this template
 */
function StringField(props) {
    var _a;
    const { schema, name, uiSchema, idSchema, formData, required, disabled = false, readonly = false, autofocus = false, onChange, onBlur, onFocus, registry, rawErrors, hideError, } = props;
    const { title, format } = schema;
    const { widgets, formContext, schemaUtils, globalUiOptions } = registry;
    const enumOptions = schemaUtils.isSelect(schema) ? (0,lib_index_js_.optionsList)(schema) : undefined;
    let defaultWidget = enumOptions ? 'select' : 'text';
    if (format && (0,lib_index_js_.hasWidget)(schema, format, widgets)) {
        defaultWidget = format;
    }
    const { widget = defaultWidget, placeholder = '', title: uiTitle, ...options } = (0,lib_index_js_.getUiOptions)(uiSchema);
    const displayLabel = schemaUtils.getDisplayLabel(schema, uiSchema, globalUiOptions);
    const label = (_a = uiTitle !== null && uiTitle !== void 0 ? uiTitle : title) !== null && _a !== void 0 ? _a : name;
    const Widget = (0,lib_index_js_.getWidget)(schema, widget, widgets);
    return ((0,jsx_runtime.jsx)(Widget, { options: { ...options, enumOptions }, schema: schema, uiSchema: uiSchema, id: idSchema.$id, name: name, label: label, hideLabel: !displayLabel, hideError: hideError, value: formData, onChange: onChange, onBlur: onBlur, onFocus: onFocus, required: required, disabled: disabled, readonly: readonly, formContext: formContext, autofocus: autofocus, registry: registry, placeholder: placeholder, rawErrors: rawErrors }));
}
/* harmony default export */ const fields_StringField = (StringField);
//# sourceMappingURL=StringField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/NullField.js

/** The `NullField` component is used to render a field in the schema is null. It also ensures that the `formData` is
 * also set to null if it has no value.
 *
 * @param props - The `FieldProps` for this template
 */
function NullField(props) {
    const { formData, onChange } = props;
    (0,index_js_.useEffect)(() => {
        if (formData === undefined) {
            onChange(null);
        }
    }, [formData, onChange]);
    return null;
}
/* harmony default export */ const fields_NullField = (NullField);
//# sourceMappingURL=NullField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/fields/index.js








function fields() {
    return {
        AnyOfField: MultiSchemaField,
        ArrayField: fields_ArrayField,
        // ArrayField falls back to SchemaField if ArraySchemaField is not defined, which it isn't by default
        BooleanField: fields_BooleanField,
        NumberField: fields_NumberField,
        ObjectField: fields_ObjectField,
        OneOfField: MultiSchemaField,
        SchemaField: fields_SchemaField,
        StringField: fields_StringField,
        NullField: fields_NullField,
    };
}
/* harmony default export */ const components_fields = (fields);
//# sourceMappingURL=index.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ArrayFieldDescriptionTemplate.js


/** The `ArrayFieldDescriptionTemplate` component renders a `DescriptionFieldTemplate` with an `id` derived from
 * the `idSchema`.
 *
 * @param props - The `ArrayFieldDescriptionProps` for the component
 */
function ArrayFieldDescriptionTemplate(props) {
    const { idSchema, description, registry, schema, uiSchema } = props;
    const options = (0,lib_index_js_.getUiOptions)(uiSchema, registry.globalUiOptions);
    const { label: displayLabel = true } = options;
    if (!description || !displayLabel) {
        return null;
    }
    const DescriptionFieldTemplate = (0,lib_index_js_.getTemplate)('DescriptionFieldTemplate', registry, options);
    return ((0,jsx_runtime.jsx)(DescriptionFieldTemplate, { id: (0,lib_index_js_.descriptionId)(idSchema), description: description, schema: schema, uiSchema: uiSchema, registry: registry }));
}
//# sourceMappingURL=ArrayFieldDescriptionTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ArrayFieldItemTemplate.js

/** The `ArrayFieldItemTemplate` component is the template used to render an items of an array.
 *
 * @param props - The `ArrayFieldTemplateItemType` props for the component
 */
function ArrayFieldItemTemplate(props) {
    const { children, className, disabled, hasToolbar, hasMoveDown, hasMoveUp, hasRemove, hasCopy, index, onCopyIndexClick, onDropIndexClick, onReorderClick, readonly, registry, uiSchema, } = props;
    const { CopyButton, MoveDownButton, MoveUpButton, RemoveButton } = registry.templates.ButtonTemplates;
    const btnStyle = {
        flex: 1,
        paddingLeft: 6,
        paddingRight: 6,
        fontWeight: 'bold',
    };
    return ((0,jsx_runtime.jsxs)("div", { className: className, children: [(0,jsx_runtime.jsx)("div", { className: hasToolbar ? 'col-xs-9' : 'col-xs-12', children: children }), hasToolbar && ((0,jsx_runtime.jsx)("div", { className: 'col-xs-3 array-item-toolbox', children: (0,jsx_runtime.jsxs)("div", { className: 'btn-group', style: {
                        display: 'flex',
                        justifyContent: 'space-around',
                    }, children: [(hasMoveUp || hasMoveDown) && ((0,jsx_runtime.jsx)(MoveUpButton, { style: btnStyle, disabled: disabled || readonly || !hasMoveUp, onClick: onReorderClick(index, index - 1), uiSchema: uiSchema, registry: registry })), (hasMoveUp || hasMoveDown) && ((0,jsx_runtime.jsx)(MoveDownButton, { style: btnStyle, disabled: disabled || readonly || !hasMoveDown, onClick: onReorderClick(index, index + 1), uiSchema: uiSchema, registry: registry })), hasCopy && ((0,jsx_runtime.jsx)(CopyButton, { style: btnStyle, disabled: disabled || readonly, onClick: onCopyIndexClick(index), uiSchema: uiSchema, registry: registry })), hasRemove && ((0,jsx_runtime.jsx)(RemoveButton, { style: btnStyle, disabled: disabled || readonly, onClick: onDropIndexClick(index), uiSchema: uiSchema, registry: registry }))] }) }))] }));
}
//# sourceMappingURL=ArrayFieldItemTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ArrayFieldTemplate.js


/** The `ArrayFieldTemplate` component is the template used to render all items in an array.
 *
 * @param props - The `ArrayFieldTemplateItemType` props for the component
 */
function ArrayFieldTemplate(props) {
    const { canAdd, className, disabled, idSchema, uiSchema, items, onAddClick, readonly, registry, required, schema, title, } = props;
    const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema);
    const ArrayFieldDescriptionTemplate = (0,lib_index_js_.getTemplate)('ArrayFieldDescriptionTemplate', registry, uiOptions);
    const ArrayFieldItemTemplate = (0,lib_index_js_.getTemplate)('ArrayFieldItemTemplate', registry, uiOptions);
    const ArrayFieldTitleTemplate = (0,lib_index_js_.getTemplate)('ArrayFieldTitleTemplate', registry, uiOptions);
    // Button templates are not overridden in the uiSchema
    const { ButtonTemplates: { AddButton }, } = registry.templates;
    return ((0,jsx_runtime.jsxs)("fieldset", { className: className, id: idSchema.$id, children: [(0,jsx_runtime.jsx)(ArrayFieldTitleTemplate, { idSchema: idSchema, title: uiOptions.title || title, required: required, schema: schema, uiSchema: uiSchema, registry: registry }), (0,jsx_runtime.jsx)(ArrayFieldDescriptionTemplate, { idSchema: idSchema, description: uiOptions.description || schema.description, schema: schema, uiSchema: uiSchema, registry: registry }), (0,jsx_runtime.jsx)("div", { className: 'row array-item-list', children: items &&
                    items.map(({ key, ...itemProps }) => ((0,jsx_runtime.jsx)(ArrayFieldItemTemplate, { ...itemProps }, key))) }), canAdd && ((0,jsx_runtime.jsx)(AddButton, { className: 'array-item-add', onClick: onAddClick, disabled: disabled || readonly, uiSchema: uiSchema, registry: registry }))] }));
}
//# sourceMappingURL=ArrayFieldTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ArrayFieldTitleTemplate.js


/** The `ArrayFieldTitleTemplate` component renders a `TitleFieldTemplate` with an `id` derived from
 * the `idSchema`.
 *
 * @param props - The `ArrayFieldTitleProps` for the component
 */
function ArrayFieldTitleTemplate(props) {
    const { idSchema, title, schema, uiSchema, required, registry } = props;
    const options = (0,lib_index_js_.getUiOptions)(uiSchema, registry.globalUiOptions);
    const { label: displayLabel = true } = options;
    if (!title || !displayLabel) {
        return null;
    }
    const TitleFieldTemplate = (0,lib_index_js_.getTemplate)('TitleFieldTemplate', registry, options);
    return ((0,jsx_runtime.jsx)(TitleFieldTemplate, { id: (0,lib_index_js_.titleId)(idSchema), title: title, required: required, schema: schema, uiSchema: uiSchema, registry: registry }));
}
//# sourceMappingURL=ArrayFieldTitleTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/BaseInputTemplate.js



/** The `BaseInputTemplate` is the template to use to render the basic `<input>` component for the `core` theme.
 * It is used as the template for rendering many of the <input> based widgets that differ by `type` and callbacks only.
 * It can be customized/overridden for other themes or individual implementations as needed.
 *
 * @param props - The `WidgetProps` for this template
 */
function BaseInputTemplate(props) {
    const { id, name, // remove this from ...rest
    value, readonly, disabled, autofocus, onBlur, onFocus, onChange, onChangeOverride, options, schema, uiSchema, formContext, registry, rawErrors, type, hideLabel, // remove this from ...rest
    hideError, // remove this from ...rest
    ...rest } = props;
    // Note: since React 15.2.0 we can't forward unknown element attributes, so we
    // exclude the "options" and "schema" ones here.
    if (!id) {
        console.log('No id for', props);
        throw new Error(`no id for props ${JSON.stringify(props)}`);
    }
    const inputProps = {
        ...rest,
        ...(0,lib_index_js_.getInputProps)(schema, type, options),
    };
    let inputValue;
    if (inputProps.type === 'number' || inputProps.type === 'integer') {
        inputValue = value || value === 0 ? value : '';
    }
    else {
        inputValue = value == null ? '' : value;
    }
    const _onChange = (0,index_js_.useCallback)(({ target: { value } }) => onChange(value === '' ? options.emptyValue : value), [onChange, options]);
    const _onBlur = (0,index_js_.useCallback)(({ target: { value } }) => onBlur(id, value), [onBlur, id]);
    const _onFocus = (0,index_js_.useCallback)(({ target: { value } }) => onFocus(id, value), [onFocus, id]);
    return ((0,jsx_runtime.jsxs)(jsx_runtime.Fragment, { children: [(0,jsx_runtime.jsx)("input", { id: id, name: id, className: 'form-control', readOnly: readonly, disabled: disabled, autoFocus: autofocus, value: inputValue, ...inputProps, list: schema.examples ? (0,lib_index_js_.examplesId)(id) : undefined, onChange: onChangeOverride || _onChange, onBlur: _onBlur, onFocus: _onFocus, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id, !!schema.examples) }), Array.isArray(schema.examples) && ((0,jsx_runtime.jsx)("datalist", { id: (0,lib_index_js_.examplesId)(id), children: schema.examples
                    .concat(schema.default && !schema.examples.includes(schema.default) ? [schema.default] : [])
                    .map((example) => {
                    return (0,jsx_runtime.jsx)("option", { value: example }, example);
                }) }, `datalist_${id}`))] }));
}
//# sourceMappingURL=BaseInputTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ButtonTemplates/SubmitButton.js


/** The `SubmitButton` renders a button that represent the `Submit` action on a form
 */
function SubmitButton({ uiSchema }) {
    const { submitText, norender, props: submitButtonProps = {} } = (0,lib_index_js_.getSubmitButtonOptions)(uiSchema);
    if (norender) {
        return null;
    }
    return ((0,jsx_runtime.jsx)("div", { children: (0,jsx_runtime.jsx)("button", { type: 'submit', ...submitButtonProps, className: `btn btn-info ${submitButtonProps.className || ''}`, children: submitText }) }));
}
//# sourceMappingURL=SubmitButton.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ButtonTemplates/IconButton.js


function IconButton(props) {
    const { iconType = 'default', icon, className, uiSchema, registry, ...otherProps } = props;
    return ((0,jsx_runtime.jsx)("button", { type: 'button', className: `btn btn-${iconType} ${className}`, ...otherProps, children: (0,jsx_runtime.jsx)("i", { className: `glyphicon glyphicon-${icon}` }) }));
}
function CopyButton(props) {
    const { registry: { translateString }, } = props;
    return ((0,jsx_runtime.jsx)(IconButton, { title: translateString(lib_index_js_.TranslatableString.CopyButton), className: 'array-item-copy', ...props, icon: 'copy' }));
}
function MoveDownButton(props) {
    const { registry: { translateString }, } = props;
    return ((0,jsx_runtime.jsx)(IconButton, { title: translateString(lib_index_js_.TranslatableString.MoveDownButton), className: 'array-item-move-down', ...props, icon: 'arrow-down' }));
}
function MoveUpButton(props) {
    const { registry: { translateString }, } = props;
    return ((0,jsx_runtime.jsx)(IconButton, { title: translateString(lib_index_js_.TranslatableString.MoveUpButton), className: 'array-item-move-up', ...props, icon: 'arrow-up' }));
}
function RemoveButton(props) {
    const { registry: { translateString }, } = props;
    return ((0,jsx_runtime.jsx)(IconButton, { title: translateString(lib_index_js_.TranslatableString.RemoveButton), className: 'array-item-remove', ...props, iconType: 'danger', icon: 'remove' }));
}
//# sourceMappingURL=IconButton.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ButtonTemplates/AddButton.js



/** The `AddButton` renders a button that represent the `Add` action on a form
 */
function AddButton({ className, onClick, disabled, registry, }) {
    const { translateString } = registry;
    return ((0,jsx_runtime.jsx)("div", { className: 'row', children: (0,jsx_runtime.jsx)("p", { className: `col-xs-3 col-xs-offset-9 text-right ${className}`, children: (0,jsx_runtime.jsx)(IconButton, { iconType: 'info', icon: 'plus', className: 'btn-add col-xs-12', title: translateString(lib_index_js_.TranslatableString.AddButton), onClick: onClick, disabled: disabled, registry: registry }) }) }));
}
//# sourceMappingURL=AddButton.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ButtonTemplates/index.js



function buttonTemplates() {
    return {
        SubmitButton: SubmitButton,
        AddButton: AddButton,
        CopyButton: CopyButton,
        MoveDownButton: MoveDownButton,
        MoveUpButton: MoveUpButton,
        RemoveButton: RemoveButton,
    };
}
/* harmony default export */ const ButtonTemplates = (buttonTemplates);
//# sourceMappingURL=index.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/DescriptionField.js

/** The `DescriptionField` is the template to use to render the description of a field
 *
 * @param props - The `DescriptionFieldProps` for this component
 */
function DescriptionField(props) {
    const { id, description } = props;
    if (!description) {
        return null;
    }
    if (typeof description === 'string') {
        return ((0,jsx_runtime.jsx)("p", { id: id, className: 'field-description', children: description }));
    }
    else {
        return ((0,jsx_runtime.jsx)("div", { id: id, className: 'field-description', children: description }));
    }
}
//# sourceMappingURL=DescriptionField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ErrorList.js


/** The `ErrorList` component is the template that renders the all the errors associated with the fields in the `Form`
 *
 * @param props - The `ErrorListProps` for this component
 */
function ErrorList({ errors, registry, }) {
    const { translateString } = registry;
    return ((0,jsx_runtime.jsxs)("div", { className: 'panel panel-danger errors', children: [(0,jsx_runtime.jsx)("div", { className: 'panel-heading', children: (0,jsx_runtime.jsx)("h3", { className: 'panel-title', children: translateString(lib_index_js_.TranslatableString.ErrorsLabel) }) }), (0,jsx_runtime.jsx)("ul", { className: 'list-group', children: errors.map((error, i) => {
                    return ((0,jsx_runtime.jsx)("li", { className: 'list-group-item text-danger', children: error.stack }, i));
                }) })] }));
}
//# sourceMappingURL=ErrorList.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/FieldTemplate/Label.js

const REQUIRED_FIELD_SYMBOL = '*';
/** Renders a label for a field
 *
 * @param props - The `LabelProps` for this component
 */
function Label(props) {
    const { label, required, id } = props;
    if (!label) {
        return null;
    }
    return ((0,jsx_runtime.jsxs)("label", { className: 'control-label', htmlFor: id, children: [label, required && (0,jsx_runtime.jsx)("span", { className: 'required', children: REQUIRED_FIELD_SYMBOL })] }));
}
//# sourceMappingURL=Label.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/FieldTemplate/FieldTemplate.js



/** The `FieldTemplate` component is the template used by `SchemaField` to render any field. It renders the field
 * content, (label, description, children, errors and help) inside of a `WrapIfAdditional` component.
 *
 * @param props - The `FieldTemplateProps` for this component
 */
function FieldTemplate(props) {
    const { id, label, children, errors, help, description, hidden, required, displayLabel, registry, uiSchema } = props;
    const uiOptions = (0,lib_index_js_.getUiOptions)(uiSchema);
    const WrapIfAdditionalTemplate = (0,lib_index_js_.getTemplate)('WrapIfAdditionalTemplate', registry, uiOptions);
    if (hidden) {
        return (0,jsx_runtime.jsx)("div", { className: 'hidden', children: children });
    }
    return ((0,jsx_runtime.jsxs)(WrapIfAdditionalTemplate, { ...props, children: [displayLabel && (0,jsx_runtime.jsx)(Label, { label: label, required: required, id: id }), displayLabel && description ? description : null, children, errors, help] }));
}
//# sourceMappingURL=FieldTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/FieldTemplate/index.js

/* harmony default export */ const templates_FieldTemplate = (FieldTemplate);
//# sourceMappingURL=index.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/FieldErrorTemplate.js


/** The `FieldErrorTemplate` component renders the errors local to the particular field
 *
 * @param props - The `FieldErrorProps` for the errors being rendered
 */
function FieldErrorTemplate(props) {
    const { errors = [], idSchema } = props;
    if (errors.length === 0) {
        return null;
    }
    const id = (0,lib_index_js_.errorId)(idSchema);
    return ((0,jsx_runtime.jsx)("div", { children: (0,jsx_runtime.jsx)("ul", { id: id, className: 'error-detail bs-callout bs-callout-info', children: errors
                .filter((elem) => !!elem)
                .map((error, index) => {
                return ((0,jsx_runtime.jsx)("li", { className: 'text-danger', children: error }, index));
            }) }) }));
}
//# sourceMappingURL=FieldErrorTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/FieldHelpTemplate.js


/** The `FieldHelpTemplate` component renders any help desired for a field
 *
 * @param props - The `FieldHelpProps` to be rendered
 */
function FieldHelpTemplate(props) {
    const { idSchema, help } = props;
    if (!help) {
        return null;
    }
    const id = (0,lib_index_js_.helpId)(idSchema);
    if (typeof help === 'string') {
        return ((0,jsx_runtime.jsx)("p", { id: id, className: 'help-block', children: help }));
    }
    return ((0,jsx_runtime.jsx)("div", { id: id, className: 'help-block', children: help }));
}
//# sourceMappingURL=FieldHelpTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/ObjectFieldTemplate.js


/** The `ObjectFieldTemplate` is the template to use to render all the inner properties of an object along with the
 * title and description if available. If the object is expandable, then an `AddButton` is also rendered after all
 * the properties.
 *
 * @param props - The `ObjectFieldTemplateProps` for this component
 */
function ObjectFieldTemplate(props) {
    const { description, disabled, formData, idSchema, onAddClick, properties, readonly, registry, required, schema, title, uiSchema, } = props;
    const options = (0,lib_index_js_.getUiOptions)(uiSchema);
    const TitleFieldTemplate = (0,lib_index_js_.getTemplate)('TitleFieldTemplate', registry, options);
    const DescriptionFieldTemplate = (0,lib_index_js_.getTemplate)('DescriptionFieldTemplate', registry, options);
    // Button templates are not overridden in the uiSchema
    const { ButtonTemplates: { AddButton }, } = registry.templates;
    return ((0,jsx_runtime.jsxs)("fieldset", { id: idSchema.$id, children: [title && ((0,jsx_runtime.jsx)(TitleFieldTemplate, { id: (0,lib_index_js_.titleId)(idSchema), title: title, required: required, schema: schema, uiSchema: uiSchema, registry: registry })), description && ((0,jsx_runtime.jsx)(DescriptionFieldTemplate, { id: (0,lib_index_js_.descriptionId)(idSchema), description: description, schema: schema, uiSchema: uiSchema, registry: registry })), properties.map((prop) => prop.content), (0,lib_index_js_.canExpand)(schema, uiSchema, formData) && ((0,jsx_runtime.jsx)(AddButton, { className: 'object-property-expand', onClick: onAddClick(schema), disabled: disabled || readonly, uiSchema: uiSchema, registry: registry }))] }));
}
//# sourceMappingURL=ObjectFieldTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/TitleField.js

const TitleField_REQUIRED_FIELD_SYMBOL = '*';
/** The `TitleField` is the template to use to render the title of a field
 *
 * @param props - The `TitleFieldProps` for this component
 */
function TitleField(props) {
    const { id, title, required } = props;
    return ((0,jsx_runtime.jsxs)("legend", { id: id, children: [title, required && (0,jsx_runtime.jsx)("span", { className: 'required', children: TitleField_REQUIRED_FIELD_SYMBOL })] }));
}
//# sourceMappingURL=TitleField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/UnsupportedField.js



/** The `UnsupportedField` component is used to render a field in the schema is one that is not supported by
 * react-jsonschema-form.
 *
 * @param props - The `FieldProps` for this template
 */
function UnsupportedField(props) {
    const { schema, idSchema, reason, registry } = props;
    const { translateString } = registry;
    let translateEnum = lib_index_js_.TranslatableString.UnsupportedField;
    const translateParams = [];
    if (idSchema && idSchema.$id) {
        translateEnum = lib_index_js_.TranslatableString.UnsupportedFieldWithId;
        translateParams.push(idSchema.$id);
    }
    if (reason) {
        translateEnum =
            translateEnum === lib_index_js_.TranslatableString.UnsupportedField
                ? lib_index_js_.TranslatableString.UnsupportedFieldWithReason
                : lib_index_js_.TranslatableString.UnsupportedFieldWithIdAndReason;
        translateParams.push(reason);
    }
    return ((0,jsx_runtime.jsxs)("div", { className: 'unsupported-field', children: [(0,jsx_runtime.jsx)("p", { children: (0,jsx_runtime.jsx)(index_modern, { children: translateString(translateEnum, translateParams) }) }), schema && (0,jsx_runtime.jsx)("pre", { children: JSON.stringify(schema, null, 2) })] }));
}
/* harmony default export */ const templates_UnsupportedField = (UnsupportedField);
//# sourceMappingURL=UnsupportedField.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/WrapIfAdditionalTemplate.js



/** The `WrapIfAdditional` component is used by the `FieldTemplate` to rename, or remove properties that are
 * part of an `additionalProperties` part of a schema.
 *
 * @param props - The `WrapIfAdditionalProps` for this component
 */
function WrapIfAdditionalTemplate(props) {
    const { id, classNames, style, disabled, label, onKeyChange, onDropPropertyClick, readonly, required, schema, children, uiSchema, registry, } = props;
    const { templates, translateString } = registry;
    // Button templates are not overridden in the uiSchema
    const { RemoveButton } = templates.ButtonTemplates;
    const keyLabel = translateString(lib_index_js_.TranslatableString.KeyLabel, [label]);
    const additional = lib_index_js_.ADDITIONAL_PROPERTY_FLAG in schema;
    if (!additional) {
        return ((0,jsx_runtime.jsx)("div", { className: classNames, style: style, children: children }));
    }
    return ((0,jsx_runtime.jsx)("div", { className: classNames, style: style, children: (0,jsx_runtime.jsxs)("div", { className: 'row', children: [(0,jsx_runtime.jsx)("div", { className: 'col-xs-5 form-additional', children: (0,jsx_runtime.jsxs)("div", { className: 'form-group', children: [(0,jsx_runtime.jsx)(Label, { label: keyLabel, required: required, id: `${id}-key` }), (0,jsx_runtime.jsx)("input", { className: 'form-control', type: 'text', id: `${id}-key`, onBlur: (event) => onKeyChange(event.target.value), defaultValue: label })] }) }), (0,jsx_runtime.jsx)("div", { className: 'form-additional form-group col-xs-5', children: children }), (0,jsx_runtime.jsx)("div", { className: 'col-xs-2', children: (0,jsx_runtime.jsx)(RemoveButton, { className: 'array-item-remove btn-block', style: { border: '0' }, disabled: disabled || readonly, onClick: onDropPropertyClick(label), uiSchema: uiSchema, registry: registry }) })] }) }));
}
//# sourceMappingURL=WrapIfAdditionalTemplate.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/templates/index.js















function templates() {
    return {
        ArrayFieldDescriptionTemplate: ArrayFieldDescriptionTemplate,
        ArrayFieldItemTemplate: ArrayFieldItemTemplate,
        ArrayFieldTemplate: ArrayFieldTemplate,
        ArrayFieldTitleTemplate: ArrayFieldTitleTemplate,
        ButtonTemplates: ButtonTemplates(),
        BaseInputTemplate: BaseInputTemplate,
        DescriptionFieldTemplate: DescriptionField,
        ErrorListTemplate: ErrorList,
        FieldTemplate: templates_FieldTemplate,
        FieldErrorTemplate: FieldErrorTemplate,
        FieldHelpTemplate: FieldHelpTemplate,
        ObjectFieldTemplate: ObjectFieldTemplate,
        TitleFieldTemplate: TitleField,
        UnsupportedFieldTemplate: templates_UnsupportedField,
        WrapIfAdditionalTemplate: WrapIfAdditionalTemplate,
    };
}
/* harmony default export */ const components_templates = (templates);
//# sourceMappingURL=index.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/AltDateWidget.js



function rangeOptions(start, stop) {
    const options = [];
    for (let i = start; i <= stop; i++) {
        options.push({ value: i, label: (0,lib_index_js_.pad)(i, 2) });
    }
    return options;
}
function readyForChange(state) {
    return Object.values(state).every((value) => value !== -1);
}
function dateElementProps(state, time, yearsRange = [1900, new Date().getFullYear() + 2]) {
    const { year, month, day, hour, minute, second } = state;
    const data = [
        {
            type: 'year',
            range: yearsRange,
            value: year,
        },
        { type: 'month', range: [1, 12], value: month },
        { type: 'day', range: [1, 31], value: day },
    ];
    if (time) {
        data.push({ type: 'hour', range: [0, 23], value: hour }, { type: 'minute', range: [0, 59], value: minute }, { type: 'second', range: [0, 59], value: second });
    }
    return data;
}
function DateElement({ type, range, value, select, rootId, name, disabled, readonly, autofocus, registry, onBlur, onFocus, }) {
    const id = rootId + '_' + type;
    const { SelectWidget } = registry.widgets;
    return ((0,jsx_runtime.jsx)(SelectWidget, { schema: { type: 'integer' }, id: id, name: name, className: 'form-control', options: { enumOptions: rangeOptions(range[0], range[1]) }, placeholder: type, value: value, disabled: disabled, readonly: readonly, autofocus: autofocus, onChange: (value) => select(type, value), onBlur: onBlur, onFocus: onFocus, registry: registry, label: '', "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(rootId) }));
}
/** The `AltDateWidget` is an alternative widget for rendering date properties.
 * @param props - The `WidgetProps` for this component
 */
function AltDateWidget({ time = false, disabled = false, readonly = false, autofocus = false, options, id, name, registry, onBlur, onFocus, onChange, value, }) {
    const { translateString } = registry;
    const [lastValue, setLastValue] = (0,index_js_.useState)(value);
    const [state, setState] = (0,index_js_.useReducer)((state, action) => {
        return { ...state, ...action };
    }, (0,lib_index_js_.parseDateString)(value, time));
    (0,index_js_.useEffect)(() => {
        const stateValue = (0,lib_index_js_.toDateString)(state, time);
        if (readyForChange(state) && stateValue !== value) {
            // The user changed the date to a new valid data via the comboboxes, so call onChange
            onChange(stateValue);
        }
        else if (lastValue !== value) {
            // We got a new value in the props
            setLastValue(value);
            setState((0,lib_index_js_.parseDateString)(value, time));
        }
    }, [time, value, onChange, state, lastValue]);
    const handleChange = (0,index_js_.useCallback)((property, value) => {
        setState({ [property]: value });
    }, []);
    const handleSetNow = (0,index_js_.useCallback)((event) => {
        event.preventDefault();
        if (disabled || readonly) {
            return;
        }
        const nextState = (0,lib_index_js_.parseDateString)(new Date().toJSON(), time);
        onChange((0,lib_index_js_.toDateString)(nextState, time));
    }, [disabled, readonly, time]);
    const handleClear = (0,index_js_.useCallback)((event) => {
        event.preventDefault();
        if (disabled || readonly) {
            return;
        }
        onChange(undefined);
    }, [disabled, readonly, onChange]);
    return ((0,jsx_runtime.jsxs)("ul", { className: 'list-inline', children: [dateElementProps(state, time, options.yearsRange).map((elemProps, i) => ((0,jsx_runtime.jsx)("li", { className: 'list-inline-item', children: (0,jsx_runtime.jsx)(DateElement, { rootId: id, name: name, select: handleChange, ...elemProps, disabled: disabled, readonly: readonly, registry: registry, onBlur: onBlur, onFocus: onFocus, autofocus: autofocus && i === 0 }) }, i))), (options.hideNowButton !== 'undefined' ? !options.hideNowButton : true) && ((0,jsx_runtime.jsx)("li", { className: 'list-inline-item', children: (0,jsx_runtime.jsx)("a", { href: '#', className: 'btn btn-info btn-now', onClick: handleSetNow, children: translateString(lib_index_js_.TranslatableString.NowLabel) }) })), (options.hideClearButton !== 'undefined' ? !options.hideClearButton : true) && ((0,jsx_runtime.jsx)("li", { className: 'list-inline-item', children: (0,jsx_runtime.jsx)("a", { href: '#', className: 'btn btn-warning btn-clear', onClick: handleClear, children: translateString(lib_index_js_.TranslatableString.ClearLabel) }) }))] }));
}
/* harmony default export */ const widgets_AltDateWidget = (AltDateWidget);
//# sourceMappingURL=AltDateWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/AltDateTimeWidget.js

/** The `AltDateTimeWidget` is an alternative widget for rendering datetime properties.
 *  It uses the AltDateWidget for rendering, with the `time` prop set to true by default.
 *
 * @param props - The `WidgetProps` for this component
 */
function AltDateTimeWidget({ time = true, ...props }) {
    const { AltDateWidget } = props.registry.widgets;
    return (0,jsx_runtime.jsx)(AltDateWidget, { time: time, ...props });
}
/* harmony default export */ const widgets_AltDateTimeWidget = (AltDateTimeWidget);
//# sourceMappingURL=AltDateTimeWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/CheckboxWidget.js



/** The `CheckBoxWidget` is a widget for rendering boolean properties.
 *  It is typically used to represent a boolean.
 *
 * @param props - The `WidgetProps` for this component
 */
function CheckboxWidget({ schema, uiSchema, options, id, value, disabled, readonly, label, hideLabel, autofocus = false, onBlur, onFocus, onChange, registry, }) {
    var _a;
    const DescriptionFieldTemplate = (0,lib_index_js_.getTemplate)('DescriptionFieldTemplate', registry, options);
    // Because an unchecked checkbox will cause html5 validation to fail, only add
    // the "required" attribute if the field value must be "true", due to the
    // "const" or "enum" keywords
    const required = (0,lib_index_js_.schemaRequiresTrueValue)(schema);
    const handleChange = (0,index_js_.useCallback)((event) => onChange(event.target.checked), [onChange]);
    const handleBlur = (0,index_js_.useCallback)((event) => onBlur(id, event.target.checked), [onBlur, id]);
    const handleFocus = (0,index_js_.useCallback)((event) => onFocus(id, event.target.checked), [onFocus, id]);
    const description = (_a = options.description) !== null && _a !== void 0 ? _a : schema.description;
    return ((0,jsx_runtime.jsxs)("div", { className: `checkbox ${disabled || readonly ? 'disabled' : ''}`, children: [!hideLabel && !!description && ((0,jsx_runtime.jsx)(DescriptionFieldTemplate, { id: (0,lib_index_js_.descriptionId)(id), description: description, schema: schema, uiSchema: uiSchema, registry: registry })), (0,jsx_runtime.jsxs)("label", { children: [(0,jsx_runtime.jsx)("input", { type: 'checkbox', id: id, name: id, checked: typeof value === 'undefined' ? false : value, required: required, disabled: disabled || readonly, autoFocus: autofocus, onChange: handleChange, onBlur: handleBlur, onFocus: handleFocus, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id) }), (0,lib_index_js_.labelValue)((0,jsx_runtime.jsx)("span", { children: label }), hideLabel)] })] }));
}
/* harmony default export */ const widgets_CheckboxWidget = (CheckboxWidget);
//# sourceMappingURL=CheckboxWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/CheckboxesWidget.js



/** The `CheckboxesWidget` is a widget for rendering checkbox groups.
 *  It is typically used to represent an array of enums.
 *
 * @param props - The `WidgetProps` for this component
 */
function CheckboxesWidget({ id, disabled, options: { inline = false, enumOptions, enumDisabled, emptyValue }, value, autofocus = false, readonly, onChange, onBlur, onFocus, }) {
    const checkboxesValues = Array.isArray(value) ? value : [value];
    const handleBlur = (0,index_js_.useCallback)(({ target: { value } }) => onBlur(id, (0,lib_index_js_.enumOptionsValueForIndex)(value, enumOptions, emptyValue)), [onBlur, id]);
    const handleFocus = (0,index_js_.useCallback)(({ target: { value } }) => onFocus(id, (0,lib_index_js_.enumOptionsValueForIndex)(value, enumOptions, emptyValue)), [onFocus, id]);
    return ((0,jsx_runtime.jsx)("div", { className: 'checkboxes', id: id, children: Array.isArray(enumOptions) &&
            enumOptions.map((option, index) => {
                const checked = (0,lib_index_js_.enumOptionsIsSelected)(option.value, checkboxesValues);
                const itemDisabled = Array.isArray(enumDisabled) && enumDisabled.indexOf(option.value) !== -1;
                const disabledCls = disabled || itemDisabled || readonly ? 'disabled' : '';
                const handleChange = (event) => {
                    if (event.target.checked) {
                        onChange((0,lib_index_js_.enumOptionsSelectValue)(index, checkboxesValues, enumOptions));
                    }
                    else {
                        onChange((0,lib_index_js_.enumOptionsDeselectValue)(index, checkboxesValues, enumOptions));
                    }
                };
                const checkbox = ((0,jsx_runtime.jsxs)("span", { children: [(0,jsx_runtime.jsx)("input", { type: 'checkbox', id: (0,lib_index_js_.optionId)(id, index), name: id, checked: checked, value: String(index), disabled: disabled || itemDisabled || readonly, autoFocus: autofocus && index === 0, onChange: handleChange, onBlur: handleBlur, onFocus: handleFocus, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id) }), (0,jsx_runtime.jsx)("span", { children: option.label })] }));
                return inline ? ((0,jsx_runtime.jsx)("label", { className: `checkbox-inline ${disabledCls}`, children: checkbox }, index)) : ((0,jsx_runtime.jsx)("div", { className: `checkbox ${disabledCls}`, children: (0,jsx_runtime.jsx)("label", { children: checkbox }) }, index));
            }) }));
}
/* harmony default export */ const widgets_CheckboxesWidget = (CheckboxesWidget);
//# sourceMappingURL=CheckboxesWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/ColorWidget.js


/** The `ColorWidget` component uses the `BaseInputTemplate` changing the type to `color` and disables it when it is
 * either disabled or readonly.
 *
 * @param props - The `WidgetProps` for this component
 */
function ColorWidget(props) {
    const { disabled, readonly, options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'color', ...props, disabled: disabled || readonly });
}
//# sourceMappingURL=ColorWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/DateWidget.js



/** The `DateWidget` component uses the `BaseInputTemplate` changing the type to `date` and transforms
 * the value to undefined when it is falsy during the `onChange` handling.
 *
 * @param props - The `WidgetProps` for this component
 */
function DateWidget(props) {
    const { onChange, options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    const handleChange = (0,index_js_.useCallback)((value) => onChange(value || undefined), [onChange]);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'date', ...props, onChange: handleChange });
}
//# sourceMappingURL=DateWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/DateTimeWidget.js


/** The `DateTimeWidget` component uses the `BaseInputTemplate` changing the type to `datetime-local` and transforms
 * the value to/from utc using the appropriate utility functions.
 *
 * @param props - The `WidgetProps` for this component
 */
function DateTimeWidget(props) {
    const { onChange, value, options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return ((0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'datetime-local', ...props, value: (0,lib_index_js_.utcToLocal)(value), onChange: (value) => onChange((0,lib_index_js_.localToUTC)(value)) }));
}
//# sourceMappingURL=DateTimeWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/EmailWidget.js


/** The `EmailWidget` component uses the `BaseInputTemplate` changing the type to `email`.
 *
 * @param props - The `WidgetProps` for this component
 */
function EmailWidget(props) {
    const { options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'email', ...props });
}
//# sourceMappingURL=EmailWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/FileWidget.js




function addNameToDataURL(dataURL, name) {
    if (dataURL === null) {
        return null;
    }
    return dataURL.replace(';base64', `;name=${encodeURIComponent(name)};base64`);
}
function processFile(file) {
    const { name, size, type } = file;
    return new Promise((resolve, reject) => {
        const reader = new window.FileReader();
        reader.onerror = reject;
        reader.onload = (event) => {
            var _a;
            if (typeof ((_a = event.target) === null || _a === void 0 ? void 0 : _a.result) === 'string') {
                resolve({
                    dataURL: addNameToDataURL(event.target.result, name),
                    name,
                    size,
                    type,
                });
            }
            else {
                resolve({
                    dataURL: null,
                    name,
                    size,
                    type,
                });
            }
        };
        reader.readAsDataURL(file);
    });
}
function processFiles(files) {
    return Promise.all(Array.from(files).map(processFile));
}
function FileInfoPreview({ fileInfo, registry, }) {
    const { translateString } = registry;
    const { dataURL, type, name } = fileInfo;
    if (!dataURL) {
        return null;
    }
    if (type.indexOf('image') !== -1) {
        return (0,jsx_runtime.jsx)("img", { src: dataURL, style: { maxWidth: '100%' }, className: 'file-preview' });
    }
    return ((0,jsx_runtime.jsxs)(jsx_runtime.Fragment, { children: [' ', (0,jsx_runtime.jsx)("a", { download: `preview-${name}`, href: dataURL, className: 'file-download', children: translateString(lib_index_js_.TranslatableString.PreviewLabel) })] }));
}
function FilesInfo({ filesInfo, registry, preview, }) {
    if (filesInfo.length === 0) {
        return null;
    }
    const { translateString } = registry;
    return ((0,jsx_runtime.jsx)("ul", { className: 'file-info', children: filesInfo.map((fileInfo, key) => {
            const { name, size, type } = fileInfo;
            return ((0,jsx_runtime.jsxs)("li", { children: [(0,jsx_runtime.jsx)(index_modern, { children: translateString(lib_index_js_.TranslatableString.FilesInfo, [name, type, String(size)]) }), preview && (0,jsx_runtime.jsx)(FileInfoPreview, { fileInfo: fileInfo, registry: registry })] }, key));
        }) }));
}
function extractFileInfo(dataURLs) {
    return dataURLs
        .filter((dataURL) => dataURL)
        .map((dataURL) => {
        const { blob, name } = (0,lib_index_js_.dataURItoBlob)(dataURL);
        return {
            dataURL,
            name: name,
            size: blob.size,
            type: blob.type,
        };
    });
}
/**
 *  The `FileWidget` is a widget for rendering file upload fields.
 *  It is typically used with a string property with data-url format.
 */
function FileWidget(props) {
    const { disabled, readonly, required, multiple, onChange, value, options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    const [filesInfo, setFilesInfo] = (0,index_js_.useState)(Array.isArray(value) ? extractFileInfo(value) : extractFileInfo([value]));
    const handleChange = (0,index_js_.useCallback)((event) => {
        if (!event.target.files) {
            return;
        }
        // Due to variances in themes, dealing with multiple files for the array case now happens one file at a time.
        // This is because we don't pass `multiple` into the `BaseInputTemplate` anymore. Instead, we deal with the single
        // file in each event and concatenate them together ourselves
        processFiles(event.target.files).then((filesInfoEvent) => {
            const newValue = filesInfoEvent.map((fileInfo) => fileInfo.dataURL);
            if (multiple) {
                setFilesInfo(filesInfo.concat(filesInfoEvent[0]));
                onChange(value.concat(newValue[0]));
            }
            else {
                setFilesInfo(filesInfoEvent);
                onChange(newValue[0]);
            }
        });
    }, [multiple, value, filesInfo, onChange]);
    return ((0,jsx_runtime.jsxs)("div", { children: [(0,jsx_runtime.jsx)(BaseInputTemplate, { ...props, disabled: disabled || readonly, type: 'file', required: value ? false : required, onChangeOverride: handleChange, value: '', accept: options.accept ? String(options.accept) : undefined }), (0,jsx_runtime.jsx)(FilesInfo, { filesInfo: filesInfo, registry: registry, preview: options.filePreview })] }));
}
/* harmony default export */ const widgets_FileWidget = (FileWidget);
//# sourceMappingURL=FileWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/HiddenWidget.js

/** The `HiddenWidget` is a widget for rendering a hidden input field.
 *  It is typically used by setting type to "hidden".
 *
 * @param props - The `WidgetProps` for this component
 */
function HiddenWidget({ id, value, }) {
    return (0,jsx_runtime.jsx)("input", { type: 'hidden', id: id, name: id, value: typeof value === 'undefined' ? '' : value });
}
/* harmony default export */ const widgets_HiddenWidget = (HiddenWidget);
//# sourceMappingURL=HiddenWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/PasswordWidget.js


/** The `PasswordWidget` component uses the `BaseInputTemplate` changing the type to `password`.
 *
 * @param props - The `WidgetProps` for this component
 */
function PasswordWidget(props) {
    const { options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'password', ...props });
}
//# sourceMappingURL=PasswordWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/RadioWidget.js



/** The `RadioWidget` is a widget for rendering a radio group.
 *  It is typically used with a string property constrained with enum options.
 *
 * @param props - The `WidgetProps` for this component
 */
function RadioWidget({ options, value, required, disabled, readonly, autofocus = false, onBlur, onFocus, onChange, id, }) {
    const { enumOptions, enumDisabled, inline, emptyValue } = options;
    const handleBlur = (0,index_js_.useCallback)(({ target: { value } }) => onBlur(id, (0,lib_index_js_.enumOptionsValueForIndex)(value, enumOptions, emptyValue)), [onBlur, id]);
    const handleFocus = (0,index_js_.useCallback)(({ target: { value } }) => onFocus(id, (0,lib_index_js_.enumOptionsValueForIndex)(value, enumOptions, emptyValue)), [onFocus, id]);
    return ((0,jsx_runtime.jsx)("div", { className: 'field-radio-group', id: id, children: Array.isArray(enumOptions) &&
            enumOptions.map((option, i) => {
                const checked = (0,lib_index_js_.enumOptionsIsSelected)(option.value, value);
                const itemDisabled = Array.isArray(enumDisabled) && enumDisabled.indexOf(option.value) !== -1;
                const disabledCls = disabled || itemDisabled || readonly ? 'disabled' : '';
                const handleChange = () => onChange(option.value);
                const radio = ((0,jsx_runtime.jsxs)("span", { children: [(0,jsx_runtime.jsx)("input", { type: 'radio', id: (0,lib_index_js_.optionId)(id, i), checked: checked, name: id, required: required, value: String(i), disabled: disabled || itemDisabled || readonly, autoFocus: autofocus && i === 0, onChange: handleChange, onBlur: handleBlur, onFocus: handleFocus, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id) }), (0,jsx_runtime.jsx)("span", { children: option.label })] }));
                return inline ? ((0,jsx_runtime.jsx)("label", { className: `radio-inline ${disabledCls}`, children: radio }, i)) : ((0,jsx_runtime.jsx)("div", { className: `radio ${disabledCls}`, children: (0,jsx_runtime.jsx)("label", { children: radio }) }, i));
            }) }));
}
/* harmony default export */ const widgets_RadioWidget = (RadioWidget);
//# sourceMappingURL=RadioWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/RangeWidget.js

/** The `RangeWidget` component uses the `BaseInputTemplate` changing the type to `range` and wrapping the result
 * in a div, with the value along side it.
 *
 * @param props - The `WidgetProps` for this component
 */
function RangeWidget(props) {
    const { value, registry: { templates: { BaseInputTemplate }, }, } = props;
    return ((0,jsx_runtime.jsxs)("div", { className: 'field-range-wrapper', children: [(0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'range', ...props }), (0,jsx_runtime.jsx)("span", { className: 'range-view', children: value })] }));
}
//# sourceMappingURL=RangeWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/SelectWidget.js



function getValue(event, multiple) {
    if (multiple) {
        return Array.from(event.target.options)
            .slice()
            .filter((o) => o.selected)
            .map((o) => o.value);
    }
    return event.target.value;
}
/** The `SelectWidget` is a widget for rendering dropdowns.
 *  It is typically used with string properties constrained with enum options.
 *
 * @param props - The `WidgetProps` for this component
 */
function SelectWidget({ schema, id, options, value, required, disabled, readonly, multiple = false, autofocus = false, onChange, onBlur, onFocus, placeholder, }) {
    const { enumOptions, enumDisabled, emptyValue: optEmptyVal } = options;
    const emptyValue = multiple ? [] : '';
    const handleFocus = (0,index_js_.useCallback)((event) => {
        const newValue = getValue(event, multiple);
        return onFocus(id, (0,lib_index_js_.enumOptionsValueForIndex)(newValue, enumOptions, optEmptyVal));
    }, [onFocus, id, schema, multiple, options]);
    const handleBlur = (0,index_js_.useCallback)((event) => {
        const newValue = getValue(event, multiple);
        return onBlur(id, (0,lib_index_js_.enumOptionsValueForIndex)(newValue, enumOptions, optEmptyVal));
    }, [onBlur, id, schema, multiple, options]);
    const handleChange = (0,index_js_.useCallback)((event) => {
        const newValue = getValue(event, multiple);
        return onChange((0,lib_index_js_.enumOptionsValueForIndex)(newValue, enumOptions, optEmptyVal));
    }, [onChange, schema, multiple, options]);
    const selectedIndexes = (0,lib_index_js_.enumOptionsIndexForValue)(value, enumOptions, multiple);
    return ((0,jsx_runtime.jsxs)("select", { id: id, name: id, multiple: multiple, className: 'form-control', value: typeof selectedIndexes === 'undefined' ? emptyValue : selectedIndexes, required: required, disabled: disabled || readonly, autoFocus: autofocus, onBlur: handleBlur, onFocus: handleFocus, onChange: handleChange, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id), children: [!multiple && schema.default === undefined && (0,jsx_runtime.jsx)("option", { value: '', children: placeholder }), Array.isArray(enumOptions) &&
                enumOptions.map(({ value, label }, i) => {
                    const disabled = enumDisabled && enumDisabled.indexOf(value) !== -1;
                    return ((0,jsx_runtime.jsx)("option", { value: String(i), disabled: disabled, children: label }, i));
                })] }));
}
/* harmony default export */ const widgets_SelectWidget = (SelectWidget);
//# sourceMappingURL=SelectWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/TextareaWidget.js



/** The `TextareaWidget` is a widget for rendering input fields as textarea.
 *
 * @param props - The `WidgetProps` for this component
 */
function TextareaWidget({ id, options = {}, placeholder, value, required, disabled, readonly, autofocus = false, onChange, onBlur, onFocus, }) {
    const handleChange = (0,index_js_.useCallback)(({ target: { value } }) => onChange(value === '' ? options.emptyValue : value), [onChange, options.emptyValue]);
    const handleBlur = (0,index_js_.useCallback)(({ target: { value } }) => onBlur(id, value), [onBlur, id]);
    const handleFocus = (0,index_js_.useCallback)(({ target: { value } }) => onFocus(id, value), [id, onFocus]);
    return ((0,jsx_runtime.jsx)("textarea", { id: id, name: id, className: 'form-control', value: value ? value : '', placeholder: placeholder, required: required, disabled: disabled, readOnly: readonly, autoFocus: autofocus, rows: options.rows, onBlur: handleBlur, onFocus: handleFocus, onChange: handleChange, "aria-describedby": (0,lib_index_js_.ariaDescribedByIds)(id) }));
}
TextareaWidget.defaultProps = {
    autofocus: false,
    options: {},
};
/* harmony default export */ const widgets_TextareaWidget = (TextareaWidget);
//# sourceMappingURL=TextareaWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/TextWidget.js


/** The `TextWidget` component uses the `BaseInputTemplate`.
 *
 * @param props - The `WidgetProps` for this component
 */
function TextWidget(props) {
    const { options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { ...props });
}
//# sourceMappingURL=TextWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/TimeWidget.js



/** The `TimeWidget` component uses the `BaseInputTemplate` changing the type to `time` and transforms
 * the value to undefined when it is falsy during the `onChange` handling.
 *
 * @param props - The `WidgetProps` for this component
 */
function TimeWidget(props) {
    const { onChange, options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    const handleChange = (0,index_js_.useCallback)((value) => onChange(value ? `${value}:00` : undefined), [onChange]);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'time', ...props, onChange: handleChange });
}
//# sourceMappingURL=TimeWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/URLWidget.js


/** The `URLWidget` component uses the `BaseInputTemplate` changing the type to `url`.
 *
 * @param props - The `WidgetProps` for this component
 */
function URLWidget(props) {
    const { options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'url', ...props });
}
//# sourceMappingURL=URLWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/UpDownWidget.js


/** The `UpDownWidget` component uses the `BaseInputTemplate` changing the type to `number`.
 *
 * @param props - The `WidgetProps` for this component
 */
function UpDownWidget(props) {
    const { options, registry } = props;
    const BaseInputTemplate = (0,lib_index_js_.getTemplate)('BaseInputTemplate', registry, options);
    return (0,jsx_runtime.jsx)(BaseInputTemplate, { type: 'number', ...props });
}
//# sourceMappingURL=UpDownWidget.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/widgets/index.js



















function widgets() {
    return {
        AltDateWidget: widgets_AltDateWidget,
        AltDateTimeWidget: widgets_AltDateTimeWidget,
        CheckboxWidget: widgets_CheckboxWidget,
        CheckboxesWidget: widgets_CheckboxesWidget,
        ColorWidget: ColorWidget,
        DateWidget: DateWidget,
        DateTimeWidget: DateTimeWidget,
        EmailWidget: EmailWidget,
        FileWidget: widgets_FileWidget,
        HiddenWidget: widgets_HiddenWidget,
        PasswordWidget: PasswordWidget,
        RadioWidget: widgets_RadioWidget,
        RangeWidget: RangeWidget,
        SelectWidget: widgets_SelectWidget,
        TextWidget: TextWidget,
        TextareaWidget: widgets_TextareaWidget,
        TimeWidget: TimeWidget,
        UpDownWidget: UpDownWidget,
        URLWidget: URLWidget,
    };
}
/* harmony default export */ const components_widgets = (widgets);
//# sourceMappingURL=index.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/getDefaultRegistry.js




/** The default registry consists of all the fields, templates and widgets provided in the core implementation,
 * plus an empty `rootSchema` and `formContext. We omit schemaUtils here because it cannot be defaulted without a
 * rootSchema and validator. It will be added into the computed registry later in the Form.
 */
function getDefaultRegistry() {
    return {
        fields: components_fields(),
        templates: components_templates(),
        widgets: components_widgets(),
        rootSchema: {},
        formContext: {},
        translateString: lib_index_js_.englishStringTranslator,
    };
}
//# sourceMappingURL=getDefaultRegistry.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/components/Form.js








/** The `Form` component renders the outer form and all the fields defined in the `schema` */
class Form_Form extends index_js_.Component {
    /** Constructs the `Form` from the `props`. Will setup the initial state from the props. It will also call the
     * `onChange` handler if the initially provided `formData` is modified to add missing default values as part of the
     * state construction.
     *
     * @param props - The initial props for the `Form`
     */
    constructor(props) {
        super(props);
        /** Returns the `formData` with only the elements specified in the `fields` list
         *
         * @param formData - The data for the `Form`
         * @param fields - The fields to keep while filtering
         */
        this.getUsedFormData = (formData, fields) => {
            // For the case of a single input form
            if (fields.length === 0 && typeof formData !== 'object') {
                return formData;
            }
            // _pick has incorrect type definition, it works with string[][], because lodash/hasIn supports it
            const data = pick_default()(formData, fields);
            if (Array.isArray(formData)) {
                return Object.keys(data).map((key) => data[key]);
            }
            return data;
        };
        /** Returns the list of field names from inspecting the `pathSchema` as well as using the `formData`
         *
         * @param pathSchema - The `PathSchema` object for the form
         * @param [formData] - The form data to use while checking for empty objects/arrays
         */
        this.getFieldNames = (pathSchema, formData) => {
            const getAllPaths = (_obj, acc = [], paths = [[]]) => {
                Object.keys(_obj).forEach((key) => {
                    if (typeof _obj[key] === 'object') {
                        const newPaths = paths.map((path) => [...path, key]);
                        // If an object is marked with additionalProperties, all its keys are valid
                        if (_obj[key][lib_index_js_.RJSF_ADDITONAL_PROPERTIES_FLAG] && _obj[key][lib_index_js_.NAME_KEY] !== '') {
                            acc.push(_obj[key][lib_index_js_.NAME_KEY]);
                        }
                        else {
                            getAllPaths(_obj[key], acc, newPaths);
                        }
                    }
                    else if (key === lib_index_js_.NAME_KEY && _obj[key] !== '') {
                        paths.forEach((path) => {
                            const formValue = get_default()(formData, path);
                            // adds path to fieldNames if it points to a value
                            // or an empty object/array
                            if (typeof formValue !== 'object' ||
                                isEmpty_default()(formValue) ||
                                (Array.isArray(formValue) && formValue.every((val) => typeof val !== 'object'))) {
                                acc.push(path);
                            }
                        });
                    }
                });
                return acc;
            };
            return getAllPaths(pathSchema);
        };
        /** Function to handle changes made to a field in the `Form`. This handler receives an entirely new copy of the
         * `formData` along with a new `ErrorSchema`. It will first update the `formData` with any missing default fields and
         * then, if `omitExtraData` and `liveOmit` are turned on, the `formData` will be filterer to remove any extra data not
         * in a form field. Then, the resulting formData will be validated if required. The state will be updated with the new
         * updated (potentially filtered) `formData`, any errors that resulted from validation. Finally the `onChange`
         * callback will be called if specified with the updated state.
         *
         * @param formData - The new form data from a change to a field
         * @param newErrorSchema - The new `ErrorSchema` based on the field change
         * @param id - The id of the field that caused the change
         */
        this.onChange = (formData, newErrorSchema, id) => {
            const { extraErrors, omitExtraData, liveOmit, noValidate, liveValidate, onChange } = this.props;
            const { schemaUtils, schema, retrievedSchema } = this.state;
            if ((0,lib_index_js_.isObject)(formData) || Array.isArray(formData)) {
                const newState = this.getStateFromProps(this.props, formData, retrievedSchema);
                formData = newState.formData;
            }
            const mustValidate = !noValidate && liveValidate;
            let state = { formData, schema };
            let newFormData = formData;
            let _retrievedSchema;
            if (omitExtraData === true && liveOmit === true) {
                _retrievedSchema = schemaUtils.retrieveSchema(schema, formData);
                const pathSchema = schemaUtils.toPathSchema(_retrievedSchema, '', formData);
                const fieldNames = this.getFieldNames(pathSchema, formData);
                newFormData = this.getUsedFormData(formData, fieldNames);
                state = {
                    formData: newFormData,
                };
            }
            if (mustValidate) {
                const schemaValidation = this.validate(newFormData, schema, schemaUtils, retrievedSchema);
                let errors = schemaValidation.errors;
                let errorSchema = schemaValidation.errorSchema;
                const schemaValidationErrors = errors;
                const schemaValidationErrorSchema = errorSchema;
                if (extraErrors) {
                    const merged = (0,lib_index_js_.validationDataMerge)(schemaValidation, extraErrors);
                    errorSchema = merged.errorSchema;
                    errors = merged.errors;
                }
                state = {
                    formData: newFormData,
                    errors,
                    errorSchema,
                    schemaValidationErrors,
                    schemaValidationErrorSchema,
                };
            }
            else if (!noValidate && newErrorSchema) {
                const errorSchema = extraErrors
                    ? (0,lib_index_js_.mergeObjects)(newErrorSchema, extraErrors, 'preventDuplicates')
                    : newErrorSchema;
                state = {
                    formData: newFormData,
                    errorSchema: errorSchema,
                    errors: (0,lib_index_js_.toErrorList)(errorSchema),
                };
            }
            if (_retrievedSchema) {
                state.retrievedSchema = _retrievedSchema;
            }
            this.setState(state, () => onChange && onChange({ ...this.state, ...state }, id));
        };
        /**
         * Callback function to handle reset form data.
         * - Reset all fields with default values.
         * - Reset validations and errors
         *
         */
        this.reset = () => {
            const { onChange } = this.props;
            const newState = this.getStateFromProps(this.props, undefined);
            const newFormData = newState.formData;
            const state = {
                formData: newFormData,
                errorSchema: {},
                errors: [],
                schemaValidationErrors: [],
                schemaValidationErrorSchema: {},
            };
            this.setState(state, () => onChange && onChange({ ...this.state, ...state }));
        };
        /** Callback function to handle when a field on the form is blurred. Calls the `onBlur` callback for the `Form` if it
         * was provided.
         *
         * @param id - The unique `id` of the field that was blurred
         * @param data - The data associated with the field that was blurred
         */
        this.onBlur = (id, data) => {
            const { onBlur } = this.props;
            if (onBlur) {
                onBlur(id, data);
            }
        };
        /** Callback function to handle when a field on the form is focused. Calls the `onFocus` callback for the `Form` if it
         * was provided.
         *
         * @param id - The unique `id` of the field that was focused
         * @param data - The data associated with the field that was focused
         */
        this.onFocus = (id, data) => {
            const { onFocus } = this.props;
            if (onFocus) {
                onFocus(id, data);
            }
        };
        /** Callback function to handle when the form is submitted. First, it prevents the default event behavior. Nothing
         * happens if the target and currentTarget of the event are not the same. It will omit any extra data in the
         * `formData` in the state if `omitExtraData` is true. It will validate the resulting `formData`, reporting errors
         * via the `onError()` callback unless validation is disabled. Finally, it will add in any `extraErrors` and then call
         * back the `onSubmit` callback if it was provided.
         *
         * @param event - The submit HTML form event
         */
        this.onSubmit = (event) => {
            event.preventDefault();
            if (event.target !== event.currentTarget) {
                return;
            }
            event.persist();
            const { omitExtraData, extraErrors, noValidate, onSubmit } = this.props;
            let { formData: newFormData } = this.state;
            const { schema, schemaUtils } = this.state;
            if (omitExtraData === true) {
                const retrievedSchema = schemaUtils.retrieveSchema(schema, newFormData);
                const pathSchema = schemaUtils.toPathSchema(retrievedSchema, '', newFormData);
                const fieldNames = this.getFieldNames(pathSchema, newFormData);
                newFormData = this.getUsedFormData(newFormData, fieldNames);
            }
            if (noValidate || this.validateForm()) {
                // There are no errors generated through schema validation.
                // Check for user provided errors and update state accordingly.
                const errorSchema = extraErrors || {};
                const errors = extraErrors ? (0,lib_index_js_.toErrorList)(extraErrors) : [];
                this.setState({
                    formData: newFormData,
                    errors,
                    errorSchema,
                    schemaValidationErrors: [],
                    schemaValidationErrorSchema: {},
                }, () => {
                    if (onSubmit) {
                        onSubmit({ ...this.state, formData: newFormData, status: 'submitted' }, event);
                    }
                });
            }
        };
        if (!props.validator) {
            throw new Error('A validator is required for Form functionality to work');
        }
        this.state = this.getStateFromProps(props, props.formData);
        if (this.props.onChange && !(0,lib_index_js_.deepEquals)(this.state.formData, this.props.formData)) {
            this.props.onChange(this.state);
        }
        this.formElement = (0,index_js_.createRef)();
    }
    /**
     * `getSnapshotBeforeUpdate` is a React lifecycle method that is invoked right before the most recently rendered
     * output is committed to the DOM. It enables your component to capture current values (e.g., scroll position) before
     * they are potentially changed.
     *
     * In this case, it checks if the props have changed since the last render. If they have, it computes the next state
     * of the component using `getStateFromProps` method and returns it along with a `shouldUpdate` flag set to `true` IF
     * the `nextState` and `prevState` are different, otherwise `false`. This ensures that we have the most up-to-date
     * state ready to be applied in `componentDidUpdate`.
     *
     * If `formData` hasn't changed, it simply returns an object with `shouldUpdate` set to `false`, indicating that a
     * state update is not necessary.
     *
     * @param prevProps - The previous set of props before the update.
     * @param prevState - The previous state before the update.
     * @returns Either an object containing the next state and a flag indicating that an update should occur, or an object
     *        with a flag indicating that an update is not necessary.
     */
    getSnapshotBeforeUpdate(prevProps, prevState) {
        if (!(0,lib_index_js_.deepEquals)(this.props, prevProps)) {
            const nextState = this.getStateFromProps(this.props, this.props.formData, prevProps.schema !== this.props.schema ? undefined : this.state.retrievedSchema);
            const shouldUpdate = !(0,lib_index_js_.deepEquals)(nextState, prevState);
            return { nextState, shouldUpdate };
        }
        return { shouldUpdate: false };
    }
    /**
     * `componentDidUpdate` is a React lifecycle method that is invoked immediately after updating occurs. This method is
     * not called for the initial render.
     *
     * Here, it checks if an update is necessary based on the `shouldUpdate` flag received from `getSnapshotBeforeUpdate`.
     * If an update is required, it applies the next state and, if needed, triggers the `onChange` handler to inform about
     * changes.
     *
     * This method effectively replaces the deprecated `UNSAFE_componentWillReceiveProps`, providing a safer alternative
     * to handle prop changes and state updates.
     *
     * @param _ - The previous set of props.
     * @param prevState - The previous state of the component before the update.
     * @param snapshot - The value returned from `getSnapshotBeforeUpdate`.
     */
    componentDidUpdate(_, prevState, snapshot) {
        if (snapshot.shouldUpdate) {
            const { nextState } = snapshot;
            if (!(0,lib_index_js_.deepEquals)(nextState.formData, this.props.formData) &&
                !(0,lib_index_js_.deepEquals)(nextState.formData, prevState.formData) &&
                this.props.onChange) {
                this.props.onChange(nextState);
            }
            this.setState(nextState);
        }
    }
    /** Extracts the updated state from the given `props` and `inputFormData`. As part of this process, the
     * `inputFormData` is first processed to add any missing required defaults. After that, the data is run through the
     * validation process IF required by the `props`.
     *
     * @param props - The props passed to the `Form`
     * @param inputFormData - The new or current data for the `Form`
     * @returns - The new state for the `Form`
     */
    getStateFromProps(props, inputFormData, retrievedSchema) {
        const state = this.state || {};
        const schema = 'schema' in props ? props.schema : this.props.schema;
        const uiSchema = ('uiSchema' in props ? props.uiSchema : this.props.uiSchema) || {};
        const edit = typeof inputFormData !== 'undefined';
        const liveValidate = 'liveValidate' in props ? props.liveValidate : this.props.liveValidate;
        const mustValidate = edit && !props.noValidate && liveValidate;
        const rootSchema = schema;
        const experimental_defaultFormStateBehavior = 'experimental_defaultFormStateBehavior' in props
            ? props.experimental_defaultFormStateBehavior
            : this.props.experimental_defaultFormStateBehavior;
        let schemaUtils = state.schemaUtils;
        if (!schemaUtils ||
            schemaUtils.doesSchemaUtilsDiffer(props.validator, rootSchema, experimental_defaultFormStateBehavior)) {
            schemaUtils = (0,lib_index_js_.createSchemaUtils)(props.validator, rootSchema, experimental_defaultFormStateBehavior);
        }
        const formData = schemaUtils.getDefaultFormState(schema, inputFormData);
        const _retrievedSchema = retrievedSchema !== null && retrievedSchema !== void 0 ? retrievedSchema : schemaUtils.retrieveSchema(schema, formData);
        const getCurrentErrors = () => {
            if (props.noValidate) {
                return { errors: [], errorSchema: {} };
            }
            else if (!props.liveValidate) {
                return {
                    errors: state.schemaValidationErrors || [],
                    errorSchema: state.schemaValidationErrorSchema || {},
                };
            }
            return {
                errors: state.errors || [],
                errorSchema: state.errorSchema || {},
            };
        };
        let errors;
        let errorSchema;
        let schemaValidationErrors = state.schemaValidationErrors;
        let schemaValidationErrorSchema = state.schemaValidationErrorSchema;
        if (mustValidate) {
            const schemaValidation = this.validate(formData, schema, schemaUtils, _retrievedSchema);
            errors = schemaValidation.errors;
            errorSchema = schemaValidation.errorSchema;
            schemaValidationErrors = errors;
            schemaValidationErrorSchema = errorSchema;
        }
        else {
            const currentErrors = getCurrentErrors();
            errors = currentErrors.errors;
            errorSchema = currentErrors.errorSchema;
        }
        if (props.extraErrors) {
            const merged = (0,lib_index_js_.validationDataMerge)({ errorSchema, errors }, props.extraErrors);
            errorSchema = merged.errorSchema;
            errors = merged.errors;
        }
        const idSchema = schemaUtils.toIdSchema(_retrievedSchema, uiSchema['ui:rootFieldId'], formData, props.idPrefix, props.idSeparator);
        const nextState = {
            schemaUtils,
            schema,
            uiSchema,
            idSchema,
            formData,
            edit,
            errors,
            errorSchema,
            schemaValidationErrors,
            schemaValidationErrorSchema,
            retrievedSchema: _retrievedSchema,
        };
        return nextState;
    }
    /** React lifecycle method that is used to determine whether component should be updated.
     *
     * @param nextProps - The next version of the props
     * @param nextState - The next version of the state
     * @returns - True if the component should be updated, false otherwise
     */
    shouldComponentUpdate(nextProps, nextState) {
        return (0,lib_index_js_.shouldRender)(this, nextProps, nextState);
    }
    /** Validates the `formData` against the `schema` using the `altSchemaUtils` (if provided otherwise it uses the
     * `schemaUtils` in the state), returning the results.
     *
     * @param formData - The new form data to validate
     * @param schema - The schema used to validate against
     * @param altSchemaUtils - The alternate schemaUtils to use for validation
     */
    validate(formData, schema = this.props.schema, altSchemaUtils, retrievedSchema) {
        const schemaUtils = altSchemaUtils ? altSchemaUtils : this.state.schemaUtils;
        const { customValidate, transformErrors, uiSchema } = this.props;
        const resolvedSchema = retrievedSchema !== null && retrievedSchema !== void 0 ? retrievedSchema : schemaUtils.retrieveSchema(schema, formData);
        return schemaUtils
            .getValidator()
            .validateFormData(formData, resolvedSchema, customValidate, transformErrors, uiSchema);
    }
    /** Renders any errors contained in the `state` in using the `ErrorList`, if not disabled by `showErrorList`. */
    renderErrors(registry) {
        const { errors, errorSchema, schema, uiSchema } = this.state;
        const { formContext } = this.props;
        const options = (0,lib_index_js_.getUiOptions)(uiSchema);
        const ErrorListTemplate = (0,lib_index_js_.getTemplate)('ErrorListTemplate', registry, options);
        if (errors && errors.length) {
            return ((0,jsx_runtime.jsx)(ErrorListTemplate, { errors: errors, errorSchema: errorSchema || {}, schema: schema, uiSchema: uiSchema, formContext: formContext, registry: registry }));
        }
        return null;
    }
    /** Returns the registry for the form */
    getRegistry() {
        var _a;
        const { translateString: customTranslateString, uiSchema = {} } = this.props;
        const { schemaUtils } = this.state;
        const { fields, templates, widgets, formContext, translateString } = getDefaultRegistry();
        return {
            fields: { ...fields, ...this.props.fields },
            templates: {
                ...templates,
                ...this.props.templates,
                ButtonTemplates: {
                    ...templates.ButtonTemplates,
                    ...(_a = this.props.templates) === null || _a === void 0 ? void 0 : _a.ButtonTemplates,
                },
            },
            widgets: { ...widgets, ...this.props.widgets },
            rootSchema: this.props.schema,
            formContext: this.props.formContext || formContext,
            schemaUtils,
            translateString: customTranslateString || translateString,
            globalUiOptions: uiSchema[lib_index_js_.UI_GLOBAL_OPTIONS_KEY],
        };
    }
    /** Provides a function that can be used to programmatically submit the `Form` */
    submit() {
        if (this.formElement.current) {
            this.formElement.current.dispatchEvent(new CustomEvent('submit', {
                cancelable: true,
            }));
            this.formElement.current.requestSubmit();
        }
    }
    /** Attempts to focus on the field associated with the `error`. Uses the `property` field to compute path of the error
     * field, then, using the `idPrefix` and `idSeparator` converts that path into an id. Then the input element with that
     * id is attempted to be found using the `formElement` ref. If it is located, then it is focused.
     *
     * @param error - The error on which to focus
     */
    focusOnError(error) {
        const { idPrefix = 'root', idSeparator = '_' } = this.props;
        const { property } = error;
        const path = toPath_default()(property);
        if (path[0] === '') {
            // Most of the time the `.foo` property results in the first element being empty, so replace it with the idPrefix
            path[0] = idPrefix;
        }
        else {
            // Otherwise insert the idPrefix into the first location using unshift
            path.unshift(idPrefix);
        }
        const elementId = path.join(idSeparator);
        let field = this.formElement.current.elements[elementId];
        if (!field) {
            // if not an exact match, try finding an input starting with the element id (like radio buttons or checkboxes)
            field = this.formElement.current.querySelector(`input[id^=${elementId}`);
        }
        if (field && field.length) {
            // If we got a list with length > 0
            field = field[0];
        }
        if (field) {
            field.focus();
        }
    }
    /** Programmatically validate the form. If `onError` is provided, then it will be called with the list of errors the
     * same way as would happen on form submission.
     *
     * @returns - True if the form is valid, false otherwise.
     */
    validateForm() {
        const { extraErrors, extraErrorsBlockSubmit, focusOnFirstError, onError } = this.props;
        const { formData, errors: prevErrors } = this.state;
        const schemaValidation = this.validate(formData);
        let errors = schemaValidation.errors;
        let errorSchema = schemaValidation.errorSchema;
        const schemaValidationErrors = errors;
        const schemaValidationErrorSchema = errorSchema;
        const hasError = errors.length > 0 || (extraErrors && extraErrorsBlockSubmit);
        if (hasError) {
            if (extraErrors) {
                const merged = (0,lib_index_js_.validationDataMerge)(schemaValidation, extraErrors);
                errorSchema = merged.errorSchema;
                errors = merged.errors;
            }
            if (focusOnFirstError) {
                if (typeof focusOnFirstError === 'function') {
                    focusOnFirstError(errors[0]);
                }
                else {
                    this.focusOnError(errors[0]);
                }
            }
            this.setState({
                errors,
                errorSchema,
                schemaValidationErrors,
                schemaValidationErrorSchema,
            }, () => {
                if (onError) {
                    onError(errors);
                }
                else {
                    console.error('Form validation failed', errors);
                }
            });
        }
        else if (prevErrors.length > 0) {
            this.setState({
                errors: [],
                errorSchema: {},
                schemaValidationErrors: [],
                schemaValidationErrorSchema: {},
            });
        }
        return !hasError;
    }
    /** Renders the `Form` fields inside the <form> | `tagName` or `_internalFormWrapper`, rendering any errors if
     * needed along with the submit button or any children of the form.
     */
    render() {
        const { children, id, idPrefix, idSeparator, className = '', tagName, name, method, target, action, autoComplete, enctype, acceptcharset, noHtml5Validate = false, disabled = false, readonly = false, formContext, showErrorList = 'top', _internalFormWrapper, } = this.props;
        const { schema, uiSchema, formData, errorSchema, idSchema } = this.state;
        const registry = this.getRegistry();
        const { SchemaField: _SchemaField } = registry.fields;
        const { SubmitButton } = registry.templates.ButtonTemplates;
        // The `semantic-ui` and `material-ui` themes have `_internalFormWrapper`s that take an `as` prop that is the
        // PropTypes.elementType to use for the inner tag, so we'll need to pass `tagName` along if it is provided.
        // NOTE, the `as` prop is native to `semantic-ui` and is emulated in the `material-ui` theme
        const as = _internalFormWrapper ? tagName : undefined;
        const FormTag = _internalFormWrapper || tagName || 'form';
        let { [lib_index_js_.SUBMIT_BTN_OPTIONS_KEY]: submitOptions = {} } = (0,lib_index_js_.getUiOptions)(uiSchema);
        if (disabled) {
            submitOptions = { ...submitOptions, props: { ...submitOptions.props, disabled: true } };
        }
        const submitUiSchema = { [lib_index_js_.UI_OPTIONS_KEY]: { [lib_index_js_.SUBMIT_BTN_OPTIONS_KEY]: submitOptions } };
        return ((0,jsx_runtime.jsxs)(FormTag, { className: className ? className : 'rjsf', id: id, name: name, method: method, target: target, action: action, autoComplete: autoComplete, encType: enctype, acceptCharset: acceptcharset, noValidate: noHtml5Validate, onSubmit: this.onSubmit, as: as, ref: this.formElement, children: [showErrorList === 'top' && this.renderErrors(registry), (0,jsx_runtime.jsx)(_SchemaField, { name: '', schema: schema, uiSchema: uiSchema, errorSchema: errorSchema, idSchema: idSchema, idPrefix: idPrefix, idSeparator: idSeparator, formContext: formContext, formData: formData, onChange: this.onChange, onBlur: this.onBlur, onFocus: this.onFocus, registry: registry, disabled: disabled, readonly: readonly }), children ? children : (0,jsx_runtime.jsx)(SubmitButton, { uiSchema: submitUiSchema, registry: registry }), showErrorList === 'bottom' && this.renderErrors(registry)] }));
    }
}
//# sourceMappingURL=Form.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/withTheme.js



/** A Higher-Order component that creates a wrapper around a `Form` with the overrides from the `WithThemeProps` */
function withTheme(themeProps) {
    return forwardRef(({ fields, widgets, templates, ...directProps }, ref) => {
        var _a;
        fields = { ...themeProps === null || themeProps === void 0 ? void 0 : themeProps.fields, ...fields };
        widgets = { ...themeProps === null || themeProps === void 0 ? void 0 : themeProps.widgets, ...widgets };
        templates = {
            ...themeProps === null || themeProps === void 0 ? void 0 : themeProps.templates,
            ...templates,
            ButtonTemplates: {
                ...(_a = themeProps === null || themeProps === void 0 ? void 0 : themeProps.templates) === null || _a === void 0 ? void 0 : _a.ButtonTemplates,
                ...templates === null || templates === void 0 ? void 0 : templates.ButtonTemplates,
            },
        };
        return (_jsx(Form, { ...themeProps, ...directProps, fields: fields, widgets: widgets, templates: templates, ref: ref }));
    });
}
//# sourceMappingURL=withTheme.js.map
;// CONCATENATED MODULE: ../node_modules/@rjsf/core/lib/index.js




/* harmony default export */ const lib = (Form_Form);
//# sourceMappingURL=index.js.map

/***/ }),

/***/ 8843:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Cache: () => (/* binding */ Cache),
/* harmony export */   FreeStyle: () => (/* binding */ FreeStyle),
/* harmony export */   Rule: () => (/* binding */ Rule),
/* harmony export */   Selector: () => (/* binding */ Selector),
/* harmony export */   Style: () => (/* binding */ Style),
/* harmony export */   create: () => (/* binding */ create)
/* harmony export */ });
/**
 * The unique id is used for unique hashes.
 */
let uniqueId = 0;
/**
 * Quick dictionary lookup for unit-less numbers.
 */
const CSS_NUMBER = Object.create(null);
/**
 * CSS properties that are valid unit-less numbers.
 *
 * Ref: https://github.com/facebook/react/blob/master/packages/react-dom/src/shared/CSSProperty.js
 */
const CSS_NUMBER_KEYS = [
    "animation-iteration-count",
    "border-image-outset",
    "border-image-slice",
    "border-image-width",
    "box-flex",
    "box-flex-group",
    "box-ordinal-group",
    "column-count",
    "columns",
    "counter-increment",
    "counter-reset",
    "flex",
    "flex-grow",
    "flex-positive",
    "flex-shrink",
    "flex-negative",
    "flex-order",
    "font-weight",
    "grid-area",
    "grid-column",
    "grid-column-end",
    "grid-column-span",
    "grid-column-start",
    "grid-row",
    "grid-row-end",
    "grid-row-span",
    "grid-row-start",
    "line-clamp",
    "line-height",
    "opacity",
    "order",
    "orphans",
    "tab-size",
    "widows",
    "z-index",
    "zoom",
    // SVG properties.
    "fill-opacity",
    "flood-opacity",
    "stop-opacity",
    "stroke-dasharray",
    "stroke-dashoffset",
    "stroke-miterlimit",
    "stroke-opacity",
    "stroke-width"
];
// Add vendor prefixes to all unit-less properties.
for (const property of CSS_NUMBER_KEYS) {
    for (const prefix of ["-webkit-", "-ms-", "-moz-", "-o-", ""]) {
        CSS_NUMBER[prefix + property] = true;
    }
}
/**
 * Escape a CSS class name.
 */
function escape(str) {
    return str.replace(/[ !#$%&()*+,./;<=>?@[\]^`{|}~"'\\]/g, "\\$&");
}
/**
 * Transform a JavaScript property into a CSS property.
 */
function hyphenate(propertyName) {
    return propertyName
        .replace(/[A-Z]/g, (m) => `-${m.toLowerCase()}`)
        .replace(/^ms-/, "-ms-"); // Internet Explorer vendor prefix.
}
/**
 * Generate a hash value from a string.
 */
function stringHash(str) {
    let value = 5381;
    let len = str.length;
    while (len--)
        value = (value * 33) ^ str.charCodeAt(len);
    return (value >>> 0).toString(36);
}
/**
 * Transform a style string to a CSS string.
 */
function styleToString(key, value) {
    if (value && typeof value === "number" && !CSS_NUMBER[key]) {
        return `${key}:${value}px`;
    }
    return `${key}:${value}`;
}
/**
 * Sort an array of tuples by first value.
 */
function sortTuples(value) {
    return value.sort((a, b) => (a[0] > b[0] ? 1 : -1));
}
/**
 * Categorize user styles.
 */
function parseStyles(styles, hasNestedStyles) {
    const properties = [];
    const nestedStyles = [];
    // Sort keys before adding to styles.
    for (const key of Object.keys(styles)) {
        const name = key.trim();
        const value = styles[key];
        if (name.charCodeAt(0) !== 36 /* $ */ && value != null) {
            if (typeof value === "object" && !Array.isArray(value)) {
                nestedStyles.push([name, value]);
            }
            else {
                properties.push([hyphenate(name), value]);
            }
        }
    }
    return {
        style: stringifyProperties(sortTuples(properties)),
        nested: hasNestedStyles ? nestedStyles : sortTuples(nestedStyles),
        isUnique: !!styles.$unique
    };
}
/**
 * Stringify an array of property tuples.
 */
function stringifyProperties(properties) {
    return properties
        .map(([name, value]) => {
        if (!Array.isArray(value))
            return styleToString(name, value);
        return value.map(x => styleToString(name, x)).join(";");
    })
        .join(";");
}
/**
 * Interpolate CSS selectors.
 */
function interpolate(selector, parent) {
    if (selector.indexOf("&") === -1)
        return `${parent} ${selector}`;
    return selector.replace(/&/g, parent);
}
/**
 * Recursive loop building styles with deferred selectors.
 */
function stylize(selector, styles, rulesList, stylesList, parent) {
    const { style, nested, isUnique } = parseStyles(styles, selector !== "");
    let pid = style;
    if (selector.charCodeAt(0) === 64 /* @ */) {
        const child = {
            selector,
            styles: [],
            rules: [],
            style: parent ? "" : style
        };
        rulesList.push(child);
        // Nested styles support (e.g. `.foo > @media > .bar`).
        if (style && parent) {
            child.styles.push({ selector: parent, style, isUnique });
        }
        for (const [name, value] of nested) {
            pid += name + stylize(name, value, child.rules, child.styles, parent);
        }
    }
    else {
        const key = parent ? interpolate(selector, parent) : selector;
        if (style)
            stylesList.push({ selector: key, style, isUnique });
        for (const [name, value] of nested) {
            pid += name + stylize(name, value, rulesList, stylesList, key);
        }
    }
    return pid;
}
/**
 * Transform `stylize` tree into style objects.
 */
function composeStylize(cache, pid, rulesList, stylesList, className, isStyle) {
    for (const { selector, style, isUnique } of stylesList) {
        const key = isStyle ? interpolate(selector, className) : selector;
        const id = isUnique
            ? `u\0${(++uniqueId).toString(36)}`
            : `s\0${pid}\0${style}`;
        const item = new Style(style, id);
        item.add(new Selector(key, `k\0${pid}\0${key}`));
        cache.add(item);
    }
    for (const { selector, style, rules, styles } of rulesList) {
        const item = new Rule(selector, style, `r\0${pid}\0${selector}\0${style}`);
        composeStylize(item, pid, rules, styles, className, isStyle);
        cache.add(item);
    }
}
/**
 * Cache to list to styles.
 */
function join(arr) {
    let res = "";
    for (let i = 0; i < arr.length; i++)
        res += arr[i];
    return res;
}
/**
 * Noop changes.
 */
const noopChanges = {
    add: () => undefined,
    change: () => undefined,
    remove: () => undefined
};
/**
 * Implement a cache/event emitter.
 */
class Cache {
    constructor(changes = noopChanges) {
        this.changes = changes;
        this.sheet = [];
        this.changeId = 0;
        this._keys = [];
        this._children = Object.create(null);
        this._counters = Object.create(null);
    }
    add(style) {
        const count = this._counters[style.id] || 0;
        const item = this._children[style.id] || style.clone();
        this._counters[style.id] = count + 1;
        if (count === 0) {
            this._children[item.id] = item;
            this._keys.push(item.id);
            this.sheet.push(item.getStyles());
            this.changeId++;
            this.changes.add(item, this._keys.length - 1);
        }
        else if (item instanceof Cache && style instanceof Cache) {
            const curIndex = this._keys.indexOf(style.id);
            const prevItemChangeId = item.changeId;
            item.merge(style);
            if (item.changeId !== prevItemChangeId) {
                this.sheet.splice(curIndex, 1, item.getStyles());
                this.changeId++;
                this.changes.change(item, curIndex, curIndex);
            }
        }
    }
    remove(style) {
        const count = this._counters[style.id];
        if (count) {
            this._counters[style.id] = count - 1;
            const item = this._children[style.id];
            const index = this._keys.indexOf(item.id);
            if (count === 1) {
                delete this._counters[style.id];
                delete this._children[style.id];
                this._keys.splice(index, 1);
                this.sheet.splice(index, 1);
                this.changeId++;
                this.changes.remove(item, index);
            }
            else if (item instanceof Cache && style instanceof Cache) {
                const prevChangeId = item.changeId;
                item.unmerge(style);
                if (item.changeId !== prevChangeId) {
                    this.sheet.splice(index, 1, item.getStyles());
                    this.changeId++;
                    this.changes.change(item, index, index);
                }
            }
        }
    }
    values() {
        return this._keys.map(key => this._children[key]);
    }
    merge(cache) {
        for (const item of cache.values())
            this.add(item);
        return this;
    }
    unmerge(cache) {
        for (const item of cache.values())
            this.remove(item);
        return this;
    }
    clone() {
        return new Cache().merge(this);
    }
}
/**
 * Selector is a dumb class made to represent nested CSS selectors.
 */
class Selector {
    constructor(selector, id) {
        this.selector = selector;
        this.id = id;
    }
    getStyles() {
        return this.selector;
    }
    clone() {
        return this;
    }
}
/**
 * The style container registers a style string with selectors.
 */
class Style extends Cache {
    constructor(style, id) {
        super();
        this.style = style;
        this.id = id;
    }
    getStyles() {
        return `${this.sheet.join(",")}{${this.style}}`;
    }
    clone() {
        return new Style(this.style, this.id).merge(this);
    }
}
/**
 * Implement rule logic for style output.
 */
class Rule extends Cache {
    constructor(rule, style, id) {
        super();
        this.rule = rule;
        this.style = style;
        this.id = id;
    }
    getStyles() {
        return `${this.rule}{${this.style}${join(this.sheet)}}`;
    }
    clone() {
        return new Rule(this.rule, this.style, this.id).merge(this);
    }
}
function key(pid, styles) {
    const key = `f${stringHash(pid)}`;
    if (true)
        return key;
    return `${styles.$displayName}_${key}`;
}
/**
 * The FreeStyle class implements the API for everything else.
 */
class FreeStyle extends Cache {
    constructor(id, changes) {
        super(changes);
        this.id = id;
    }
    registerStyle(styles) {
        const rulesList = [];
        const stylesList = [];
        const pid = stylize("&", styles, rulesList, stylesList);
        const id = key(pid, styles);
        const selector = `.${ true ? id : 0}`;
        composeStylize(this, pid, rulesList, stylesList, selector, true);
        return id;
    }
    registerKeyframes(keyframes) {
        return this.registerHashRule("@keyframes", keyframes);
    }
    registerHashRule(prefix, styles) {
        const rulesList = [];
        const stylesList = [];
        const pid = stylize("", styles, rulesList, stylesList);
        const id = key(pid, styles);
        const selector = `${prefix} ${ true ? id : 0}`;
        const rule = new Rule(selector, "", `h\0${pid}\0${prefix}`);
        composeStylize(rule, pid, rulesList, stylesList, "", false);
        this.add(rule);
        return id;
    }
    registerRule(rule, styles) {
        const rulesList = [];
        const stylesList = [];
        const pid = stylize(rule, styles, rulesList, stylesList);
        composeStylize(this, pid, rulesList, stylesList, "", false);
    }
    registerCss(styles) {
        return this.registerRule("", styles);
    }
    getStyles() {
        return join(this.sheet);
    }
    clone() {
        return new FreeStyle(this.id, this.changes).merge(this);
    }
}
/**
 * Exports a simple function to create a new instance.
 */
function create(changes) {
    return new FreeStyle(`f${(++uniqueId).toString(36)}`, changes);
}
//# sourceMappingURL=index.js.map

/***/ }),

/***/ 63005:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

var basePickBy = __webpack_require__(10228),
    hasIn = __webpack_require__(79749);

/**
 * The base implementation of `_.pick` without support for individual
 * property identifiers.
 *
 * @private
 * @param {Object} object The source object.
 * @param {string[]} paths The property paths to pick.
 * @returns {Object} Returns the new object.
 */
function basePick(object, paths) {
  return basePickBy(object, paths, function(value, path) {
    return hasIn(object, path);
  });
}

module.exports = basePick;


/***/ }),

/***/ 10228:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

var baseGet = __webpack_require__(79867),
    baseSet = __webpack_require__(78859),
    castPath = __webpack_require__(76747);

/**
 * The base implementation of  `_.pickBy` without support for iteratee shorthands.
 *
 * @private
 * @param {Object} object The source object.
 * @param {string[]} paths The property paths to pick.
 * @param {Function} predicate The function invoked per property.
 * @returns {Object} Returns the new object.
 */
function basePickBy(object, paths, predicate) {
  var index = -1,
      length = paths.length,
      result = {};

  while (++index < length) {
    var path = paths[index],
        value = baseGet(object, path);

    if (predicate(value, path)) {
      baseSet(result, castPath(path, object), value);
    }
  }
  return result;
}

module.exports = basePickBy;


/***/ }),

/***/ 14648:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

var basePick = __webpack_require__(63005),
    flatRest = __webpack_require__(24288);

/**
 * Creates an object composed of the picked `object` properties.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Object
 * @param {Object} object The source object.
 * @param {...(string|string[])} [paths] The property paths to pick.
 * @returns {Object} Returns the new object.
 * @example
 *
 * var object = { 'a': 1, 'b': '2', 'c': 3 };
 *
 * _.pick(object, ['a', 'c']);
 * // => { 'a': 1, 'c': 3 }
 */
var pick = flatRest(function(object, paths) {
  return object == null ? {} : basePick(object, paths);
});

module.exports = pick;


/***/ }),

/***/ 43551:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

var baseUnset = __webpack_require__(70830);

/**
 * Removes the property at `path` of `object`.
 *
 * **Note:** This method mutates `object`.
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Object
 * @param {Object} object The object to modify.
 * @param {Array|string} path The path of the property to unset.
 * @returns {boolean} Returns `true` if the property is deleted, else `false`.
 * @example
 *
 * var object = { 'a': [{ 'b': { 'c': 7 } }] };
 * _.unset(object, 'a[0].b.c');
 * // => true
 *
 * console.log(object);
 * // => { 'a': [{ 'b': {} }] };
 *
 * _.unset(object, ['a', '0', 'b', 'c']);
 * // => true
 *
 * console.log(object);
 * // => { 'a': [{ 'b': {} }] };
 */
function unset(object, path) {
  return object == null ? true : baseUnset(object, path);
}

module.exports = unset;


/***/ }),

/***/ 37634:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";
var __webpack_unused_export__;


var m = __webpack_require__(38005);
if (true) {
  exports.s = m.createRoot;
  __webpack_unused_export__ = m.hydrateRoot;
} else { var i; }


/***/ }),

/***/ 73062:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";
var __webpack_unused_export__;

__webpack_unused_export__ = ({ value: true });
var typestyle_1 = __webpack_require__(53861);
__webpack_unused_export__ = typestyle_1.TypeStyle;
/**
 * All the CSS types in the 'types' namespace
 */
var types = __webpack_require__(66720);
__webpack_unused_export__ = types;
/**
 * Export certain utilities
 */
var utilities_1 = __webpack_require__(51833);
__webpack_unused_export__ = utilities_1.extend;
__webpack_unused_export__ = utilities_1.classes;
__webpack_unused_export__ = utilities_1.media;
/** Zero configuration, default instance of TypeStyle */
var ts = new typestyle_1.TypeStyle({ autoGenerateTag: true });
/** Sets the target tag where we write the css on style updates */
__webpack_unused_export__ = ts.setStylesTarget;
/**
 * Insert `raw` CSS as a string. This is useful for e.g.
 * - third party CSS that you are customizing with template strings
 * - generating raw CSS in JavaScript
 * - reset libraries like normalize.css that you can use without loaders
 */
__webpack_unused_export__ = ts.cssRaw;
/**
 * Takes CSSProperties and registers it to a global selector (body, html, etc.)
 */
__webpack_unused_export__ = ts.cssRule;
/**
 * Renders styles to the singleton tag imediately
 * NOTE: You should only call it on initial render to prevent any non CSS flash.
 * After that it is kept sync using `requestAnimationFrame` and we haven't noticed any bad flashes.
 **/
__webpack_unused_export__ = ts.forceRenderStyles;
/**
 * Utility function to register an @font-face
 */
__webpack_unused_export__ = ts.fontFace;
/**
 * Allows use to use the stylesheet in a node.js environment
 */
__webpack_unused_export__ = ts.getStyles;
/**
 * Takes keyframes and returns a generated animationName
 */
__webpack_unused_export__ = ts.keyframes;
/**
 * Helps with testing. Reinitializes FreeStyle + raw
 */
__webpack_unused_export__ = ts.reinit;
/**
 * Takes CSSProperties and return a generated className you can use on your component
 */
exports.oB = ts.style;
/**
 * Takes an object where property names are ideal class names and property values are CSSProperties, and
 * returns an object where property names are the same ideal class names and the property values are
 * the actual generated class names using the ideal class name as the $debugName
 */
__webpack_unused_export__ = ts.stylesheet;
/**
 * Creates a new instance of TypeStyle separate from the default instance.
 *
 * - Use this for creating a different typestyle instance for a shadow dom component.
 * - Use this if you don't want an auto tag generated and you just want to collect the CSS.
 *
 * NOTE: styles aren't shared between different instances.
 */
function createTypeStyle(target) {
    var instance = new typestyle_1.TypeStyle({ autoGenerateTag: false });
    if (target) {
        instance.setStylesTarget(target);
    }
    return instance;
}
__webpack_unused_export__ = createTypeStyle;


/***/ }),

/***/ 62034:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
/**
 * We need to do the following to *our* objects before passing to freestyle:
 * - For any `$nest` directive move up to FreeStyle style nesting
 * - For any `$unique` directive map to FreeStyle Unique
 * - For any `$debugName` directive return the debug name
 */
function convertToStyles(object) {
    /** The final result we will return */
    var styles = {};
    for (var key in object) {
        /** Grab the value upfront */
        var val = object[key];
        /** TypeStyle configuration options */
        if (key === '$nest') {
            var nested = val;
            for (var selector in nested) {
                var subproperties = nested[selector];
                styles[selector] = convertToStyles(subproperties);
            }
        }
        else if (key === '$debugName') {
            styles.$displayName = val;
        }
        else {
            styles[key] = val;
        }
    }
    return styles;
}
exports.convertToStyles = convertToStyles;
// todo: better name here
function convertToKeyframes(frames) {
    var result = {};
    for (var offset in frames) {
        if (offset !== '$debugName') {
            result[offset] = frames[offset];
        }
    }
    if (frames.$debugName) {
        result.$displayName = frames.$debugName;
    }
    return result;
}
exports.convertToKeyframes = convertToKeyframes;


/***/ }),

/***/ 53861:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
var FreeStyle = __webpack_require__(8843);
var formatting_1 = __webpack_require__(62034);
var utilities_1 = __webpack_require__(51833);
/**
 * Creates an instance of free style with our options
 */
var createFreeStyle = function () { return FreeStyle.create(); };
/**
 * Maintains a single stylesheet and keeps it in sync with requested styles
 */
var TypeStyle = /** @class */ (function () {
    function TypeStyle(_a) {
        var _this = this;
        var autoGenerateTag = _a.autoGenerateTag;
        /**
         * Insert `raw` CSS as a string. This is useful for e.g.
         * - third party CSS that you are customizing with template strings
         * - generating raw CSS in JavaScript
         * - reset libraries like normalize.css that you can use without loaders
         */
        this.cssRaw = function (mustBeValidCSS) {
            if (!mustBeValidCSS) {
                return;
            }
            _this._raw += mustBeValidCSS || '';
            _this._pendingRawChange = true;
            _this._styleUpdated();
        };
        /**
         * Takes CSSProperties and registers it to a global selector (body, html, etc.)
         */
        this.cssRule = function (selector) {
            var objects = [];
            for (var _i = 1; _i < arguments.length; _i++) {
                objects[_i - 1] = arguments[_i];
            }
            var styles = formatting_1.convertToStyles(utilities_1.extend.apply(void 0, objects));
            _this._freeStyle.registerRule(selector, styles);
            _this._styleUpdated();
            return;
        };
        /**
         * Renders styles to the singleton tag imediately
         * NOTE: You should only call it on initial render to prevent any non CSS flash.
         * After that it is kept sync using `requestAnimationFrame` and we haven't noticed any bad flashes.
         **/
        this.forceRenderStyles = function () {
            var target = _this._getTag();
            if (!target) {
                return;
            }
            target.textContent = _this.getStyles();
        };
        /**
         * Utility function to register an @font-face
         */
        this.fontFace = function () {
            var fontFace = [];
            for (var _i = 0; _i < arguments.length; _i++) {
                fontFace[_i] = arguments[_i];
            }
            var freeStyle = _this._freeStyle;
            for (var _a = 0, _b = fontFace; _a < _b.length; _a++) {
                var face = _b[_a];
                freeStyle.registerRule('@font-face', face);
            }
            _this._styleUpdated();
            return;
        };
        /**
         * Allows use to use the stylesheet in a node.js environment
         */
        this.getStyles = function () {
            return (_this._raw || '') + _this._freeStyle.getStyles();
        };
        /**
         * Takes keyframes and returns a generated animationName
         */
        this.keyframes = function (frames) {
            var keyframes = formatting_1.convertToKeyframes(frames);
            // TODO: replace $debugName with display name
            var animationName = _this._freeStyle.registerKeyframes(keyframes);
            _this._styleUpdated();
            return animationName;
        };
        /**
         * Helps with testing. Reinitializes FreeStyle + raw
         */
        this.reinit = function () {
            /** reinit freestyle */
            var freeStyle = createFreeStyle();
            _this._freeStyle = freeStyle;
            _this._lastFreeStyleChangeId = freeStyle.changeId;
            /** reinit raw */
            _this._raw = '';
            _this._pendingRawChange = false;
            /** Clear any styles that were flushed */
            var target = _this._getTag();
            if (target) {
                target.textContent = '';
            }
        };
        /** Sets the target tag where we write the css on style updates */
        this.setStylesTarget = function (tag) {
            /** Clear any data in any previous tag */
            if (_this._tag) {
                _this._tag.textContent = '';
            }
            _this._tag = tag;
            /** This special time buffer immediately */
            _this.forceRenderStyles();
        };
        /**
         * Takes an object where property names are ideal class names and property values are CSSProperties, and
         * returns an object where property names are the same ideal class names and the property values are
         * the actual generated class names using the ideal class name as the $debugName
         */
        this.stylesheet = function (classes) {
            var classNames = Object.getOwnPropertyNames(classes);
            var result = {};
            for (var _i = 0, classNames_1 = classNames; _i < classNames_1.length; _i++) {
                var className = classNames_1[_i];
                var classDef = classes[className];
                if (classDef) {
                    classDef.$debugName = className;
                    result[className] = _this.style(classDef);
                }
            }
            return result;
        };
        var freeStyle = createFreeStyle();
        this._autoGenerateTag = autoGenerateTag;
        this._freeStyle = freeStyle;
        this._lastFreeStyleChangeId = freeStyle.changeId;
        this._pending = 0;
        this._pendingRawChange = false;
        this._raw = '';
        this._tag = undefined;
        // rebind prototype to TypeStyle.  It might be better to do a function() { return this.style.apply(this, arguments)}
        this.style = this.style.bind(this);
    }
    /**
     * Only calls cb all sync operations settle
     */
    TypeStyle.prototype._afterAllSync = function (cb) {
        var _this = this;
        this._pending++;
        var pending = this._pending;
        utilities_1.raf(function () {
            if (pending !== _this._pending) {
                return;
            }
            cb();
        });
    };
    TypeStyle.prototype._getTag = function () {
        if (this._tag) {
            return this._tag;
        }
        if (this._autoGenerateTag) {
            var tag = typeof window === 'undefined'
                ? { textContent: '' }
                : document.createElement('style');
            if (typeof document !== 'undefined') {
                document.head.appendChild(tag);
            }
            this._tag = tag;
            return tag;
        }
        return undefined;
    };
    /** Checks if the style tag needs updating and if so queues up the change */
    TypeStyle.prototype._styleUpdated = function () {
        var _this = this;
        var changeId = this._freeStyle.changeId;
        var lastChangeId = this._lastFreeStyleChangeId;
        if (!this._pendingRawChange && changeId === lastChangeId) {
            return;
        }
        this._lastFreeStyleChangeId = changeId;
        this._pendingRawChange = false;
        this._afterAllSync(function () { return _this.forceRenderStyles(); });
    };
    TypeStyle.prototype.style = function () {
        var className = this._freeStyle.registerStyle(formatting_1.convertToStyles(utilities_1.extend.apply(undefined, arguments)));
        this._styleUpdated();
        return className;
    };
    return TypeStyle;
}());
exports.TypeStyle = TypeStyle;


/***/ }),

/***/ 51833:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
/** Raf for node + browser */
exports.raf = typeof requestAnimationFrame === 'undefined'
    /**
     * Make sure setTimeout is always invoked with
     * `this` set to `window` or `global` automatically
     **/
    ? function (cb) { return setTimeout(cb); }
    /**
     * Make sure window.requestAnimationFrame is always invoked with `this` window
     * We might have raf without window in case of `raf/polyfill` (recommended by React)
     **/
    : typeof window === 'undefined'
        ? requestAnimationFrame
        : requestAnimationFrame.bind(window);
/**
 * Utility to join classes conditionally
 */
function classes() {
    var classes = [];
    for (var _i = 0; _i < arguments.length; _i++) {
        classes[_i] = arguments[_i];
    }
    return classes
        .map(function (c) { return c && typeof c === 'object' ? Object.keys(c).map(function (key) { return !!c[key] && key; }) : [c]; })
        .reduce(function (flattened, c) { return flattened.concat(c); }, [])
        .filter(function (c) { return !!c; })
        .join(' ');
}
exports.classes = classes;
/**
 * Merges various styles into a single style object.
 * Note: if two objects have the same property the last one wins
 */
function extend() {
    var objects = [];
    for (var _i = 0; _i < arguments.length; _i++) {
        objects[_i] = arguments[_i];
    }
    /** The final result we will return */
    var result = {};
    for (var _a = 0, objects_1 = objects; _a < objects_1.length; _a++) {
        var object = objects_1[_a];
        if (object == null || object === false) {
            continue;
        }
        for (var key in object) {
            /** Falsy values except a explicit 0 is ignored */
            var val = object[key];
            if (!val && val !== 0) {
                continue;
            }
            /** if nested media or pseudo selector */
            if (key === '$nest' && val) {
                result[key] = result['$nest'] ? extend(result['$nest'], val) : val;
            }
            /** if freestyle sub key that needs merging. We come here due to our recursive calls */
            else if ((key.indexOf('&') !== -1 || key.indexOf('@media') === 0)) {
                result[key] = result[key] ? extend(result[key], val) : val;
            }
            else {
                result[key] = val;
            }
        }
    }
    return result;
}
exports.extend = extend;
/**
 * Utility to help customize styles with media queries. e.g.
 * ```
 * style(
 *  media({maxWidth:500}, {color:'red'})
 * )
 * ```
 */
exports.media = function (mediaQuery) {
    var _a;
    var objects = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        objects[_i - 1] = arguments[_i];
    }
    var mediaQuerySections = [];
    if (mediaQuery.type)
        mediaQuerySections.push(mediaQuery.type);
    if (mediaQuery.orientation)
        mediaQuerySections.push("(orientation: " + mediaQuery.orientation + ")");
    if (mediaQuery.minWidth)
        mediaQuerySections.push("(min-width: " + mediaLength(mediaQuery.minWidth) + ")");
    if (mediaQuery.maxWidth)
        mediaQuerySections.push("(max-width: " + mediaLength(mediaQuery.maxWidth) + ")");
    if (mediaQuery.minHeight)
        mediaQuerySections.push("(min-height: " + mediaLength(mediaQuery.minHeight) + ")");
    if (mediaQuery.maxHeight)
        mediaQuerySections.push("(max-height: " + mediaLength(mediaQuery.maxHeight) + ")");
    if (mediaQuery.prefersColorScheme)
        mediaQuerySections.push("(prefers-color-scheme: " + mediaQuery.prefersColorScheme + ")");
    var stringMediaQuery = "@media " + mediaQuerySections.join(' and ');
    var object = {
        $nest: (_a = {},
            _a[stringMediaQuery] = extend.apply(void 0, objects),
            _a)
    };
    return object;
};
var mediaLength = function (value) {
    return typeof value === 'string' ? value : value + "px";
};


/***/ }),

/***/ 66720:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));


/***/ })

}]);
//# sourceMappingURL=1871.c375ee093b7e51966390.js.map?v=c375ee093b7e51966390