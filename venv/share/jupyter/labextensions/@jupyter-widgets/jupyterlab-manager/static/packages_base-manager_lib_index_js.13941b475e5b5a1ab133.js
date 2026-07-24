(self["webpackChunk_jupyter_widgets_jupyterlab_manager"] = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] || []).push([["packages_base-manager_lib_index_js"],{

/***/ "../../packages/base-manager/lib/index.js":
/*!************************************************!*\
  !*** ../../packages/base-manager/lib/index.js ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   CONTROL_COMM_PROTOCOL_VERSION: () => (/* reexport safe */ _manager_base__WEBPACK_IMPORTED_MODULE_0__.CONTROL_COMM_PROTOCOL_VERSION),
/* harmony export */   CONTROL_COMM_TARGET: () => (/* reexport safe */ _manager_base__WEBPACK_IMPORTED_MODULE_0__.CONTROL_COMM_TARGET),
/* harmony export */   CONTROL_COMM_TIMEOUT: () => (/* reexport safe */ _manager_base__WEBPACK_IMPORTED_MODULE_0__.CONTROL_COMM_TIMEOUT),
/* harmony export */   ManagerBase: () => (/* reexport safe */ _manager_base__WEBPACK_IMPORTED_MODULE_0__.ManagerBase),
/* harmony export */   base64ToBuffer: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_1__.base64ToBuffer),
/* harmony export */   bufferToBase64: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_1__.bufferToBase64),
/* harmony export */   bufferToHex: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_1__.bufferToHex),
/* harmony export */   hexToBuffer: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_1__.hexToBuffer),
/* harmony export */   serialize_state: () => (/* reexport safe */ _manager_base__WEBPACK_IMPORTED_MODULE_0__.serialize_state)
/* harmony export */ });
/* harmony import */ var _manager_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./manager-base */ "../../packages/base-manager/lib/manager-base.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./utils */ "../../packages/base-manager/lib/utils.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




/***/ }),

/***/ "../../packages/base-manager/lib/latex.js":
/*!************************************************!*\
  !*** ../../packages/base-manager/lib/latex.js ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   removeMath: () => (/* binding */ removeMath),
/* harmony export */   replaceMath: () => (/* binding */ replaceMath)
/* harmony export */ });
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/
// Some magic for deferring mathematical expressions to MathJax
// by hiding them from the Markdown parser.
// Some of the code here is adapted with permission from Davide Cervone
// under the terms of the Apache2 license governing the MathJax project.
// Other minor modifications are also due to StackExchange and are used with
// permission.
const inline = '$'; // the inline math delimiter
// MATHSPLIT contains the pattern for math delimiters and special symbols
// needed for searching for math in the text input.
const MATHSPLIT = /(\$\$?|\\(?:begin|end)\{[a-z]*\*?\}|\\[{}$]|[{}]|(?:\n\s*)+|@@\d+@@|\\\\(?:\(|\)|\[|\]))/i;
/**
 *  Break up the text into its component parts and search
 *    through them for math delimiters, braces, linebreaks, etc.
 *  Math delimiters must match and braces must balance.
 *  Don't allow math to pass through a double linebreak
 *    (which will be a paragraph).
 */
function removeMath(text) {
    const math = []; // stores math strings for later
    let start = null;
    let end = null;
    let last = null;
    let braces = 0;
    let deTilde;
    // Except for extreme edge cases, this should catch precisely those pieces of the markdown
    // source that will later be turned into code spans. While MathJax will not TeXify code spans,
    // we still have to consider them at this point; the following issue has happened several times:
    //
    //     `$foo` and `$bar` are variables.  -->  <code>$foo ` and `$bar</code> are variables.
    const hasCodeSpans = /`/.test(text);
    if (hasCodeSpans) {
        text = text
            .replace(/~/g, '~T')
            .replace(/(^|[^\\])(`+)([^\n]*?[^`\n])\2(?!`)/gm, (wholematch) => wholematch.replace(/\$/g, '~D'));
        deTilde = (text) => {
            return text.replace(/~([TD])/g, (wholematch, character) => character === 'T' ? '~' : inline);
        };
    }
    else {
        deTilde = (text) => {
            return text;
        };
    }
    let blocks = text.replace(/\r\n?/g, '\n').split(MATHSPLIT);
    for (let i = 1, m = blocks.length; i < m; i += 2) {
        const block = blocks[i];
        if (block.charAt(0) === '@') {
            //
            //  Things that look like our math markers will get
            //  stored and then retrieved along with the math.
            //
            blocks[i] = '@@' + math.length + '@@';
            math.push(block);
        }
        else if (start !== null) {
            //
            //  If we are in math, look for the end delimiter,
            //    but don't go past double line breaks, and
            //    and balance braces within the math.
            //
            if (block === end) {
                if (braces) {
                    last = i;
                }
                else {
                    blocks = processMath(start, i, deTilde, math, blocks);
                    start = null;
                    end = null;
                    last = null;
                }
            }
            else if (block.match(/\n.*\n/)) {
                if (last !== null) {
                    i = last;
                    blocks = processMath(start, i, deTilde, math, blocks);
                }
                start = null;
                end = null;
                last = null;
                braces = 0;
            }
            else if (block === '{') {
                braces++;
            }
            else if (block === '}' && braces) {
                braces--;
            }
        }
        else {
            //
            //  Look for math start delimiters and when
            //    found, set up the end delimiter.
            //
            if (block === inline || block === '$$') {
                start = i;
                end = block;
                braces = 0;
            }
            else if (block === '\\\\(' || block === '\\\\[') {
                start = i;
                end = block.slice(-1) === '(' ? '\\\\)' : '\\\\]';
                braces = 0;
            }
            else if (block.substr(1, 5) === 'begin') {
                start = i;
                end = '\\end' + block.substr(6);
                braces = 0;
            }
        }
    }
    if (start !== null && last !== null) {
        blocks = processMath(start, last, deTilde, math, blocks);
        start = null;
        end = null;
        last = null;
    }
    return { text: deTilde(blocks.join('')), math };
}
/**
 * Put back the math strings that were saved,
 * and clear the math array (no need to keep it around).
 */
function replaceMath(text, math) {
    /**
     * Replace a math placeholder with its corresponding group.
     * The math delimiters "\\(", "\\[", "\\)" and "\\]" are replaced
     * removing one backslash in order to be interpreted correctly by MathJax.
     */
    const process = (match, n) => {
        let group = math[n];
        if (group.substr(0, 3) === '\\\\(' &&
            group.substr(group.length - 3) === '\\\\)') {
            group = '\\(' + group.substring(3, group.length - 3) + '\\)';
        }
        else if (group.substr(0, 3) === '\\\\[' &&
            group.substr(group.length - 3) === '\\\\]') {
            group = '\\[' + group.substring(3, group.length - 3) + '\\]';
        }
        return group;
    };
    // Replace all the math group placeholders in the text
    // with the saved strings.
    return text.replace(/@@(\d+)@@/g, process);
}
/**
 * Process math blocks.
 *
 * The math is in blocks i through j, so
 *   collect it into one block and clear the others.
 *  Replace &, <, and > by named entities.
 *  For IE, put <br> at the ends of comments since IE removes \n.
 *  Clear the current math positions and store the index of the
 *   math, then push the math string onto the storage array.
 *  The preProcess function is called on all blocks if it has been passed in
 */
function processMath(i, j, preProcess, math, blocks) {
    let block = blocks
        .slice(i, j + 1)
        .join('')
        .replace(/&/g, '&amp;') // use HTML entity for &
        .replace(/</g, '&lt;') // use HTML entity for <
        .replace(/>/g, '&gt;'); // use HTML entity for >
    if (navigator && navigator.appName === 'Microsoft Internet Explorer') {
        block = block.replace(/(%[^\n]*)\n/g, '$1<br/>\n');
    }
    while (j > i) {
        blocks[j] = '';
        j--;
    }
    blocks[i] = '@@' + math.length + '@@'; // replace the current block text with a unique tag to find later
    if (preProcess) {
        block = preProcess(block);
    }
    math.push(block);
    return blocks;
}


/***/ }),

/***/ "../../packages/base-manager/lib/manager-base.js":
/*!*******************************************************!*\
  !*** ../../packages/base-manager/lib/manager-base.js ***!
  \*******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   CONTROL_COMM_PROTOCOL_VERSION: () => (/* binding */ CONTROL_COMM_PROTOCOL_VERSION),
/* harmony export */   CONTROL_COMM_TARGET: () => (/* binding */ CONTROL_COMM_TARGET),
/* harmony export */   CONTROL_COMM_TIMEOUT: () => (/* binding */ CONTROL_COMM_TIMEOUT),
/* harmony export */   ManagerBase: () => (/* binding */ ManagerBase),
/* harmony export */   serialize_state: () => (/* binding */ serialize_state)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/base-manager/lib/utils.js");
/* harmony import */ var _latex__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./latex */ "../../packages/base-manager/lib/latex.js");
/* harmony import */ var sanitize_html__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! sanitize-html */ "../../node_modules/sanitize-html/index.js");
/* harmony import */ var sanitize_html__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(sanitize_html__WEBPACK_IMPORTED_MODULE_4__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.






const PROTOCOL_MAJOR_VERSION = _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.PROTOCOL_VERSION.split('.', 1)[0];
/**
 * The control comm target name.
 */
const CONTROL_COMM_TARGET = 'jupyter.widget.control';
/**
 * The supported version for the control comm channel.
 */
const CONTROL_COMM_PROTOCOL_VERSION = '1.0.0';
/**
 * Time (in ms) after which we consider the control comm target not responding.
 */
const CONTROL_COMM_TIMEOUT = 4000;
/**
 * Sanitize HTML-formatted descriptions.
 */
function default_inline_sanitize(s) {
    const allowedTags = [
        'a',
        'abbr',
        'b',
        'code',
        'em',
        'i',
        'img',
        'li',
        'ol',
        'span',
        'strong',
        'ul',
    ];
    const allowedAttributes = {
        '*': ['aria-*', 'class', 'style', 'title'],
        a: ['href'],
        img: ['src'],
        style: ['media', 'type'],
    };
    return sanitize_html__WEBPACK_IMPORTED_MODULE_4___default()(s, {
        allowedTags: allowedTags,
        allowedAttributes: allowedAttributes,
    });
}
/**
 * Manager abstract base class
 */
class ManagerBase {
    constructor() {
        /**
         * The comm target name to register
         */
        this.comm_target_name = 'jupyter.widget';
        /**
         * Dictionary of model ids and model instance promises
         */
        this._models = Object.create(null);
    }
    /**
     * Modifies view options. Generally overloaded in custom widget manager
     * implementations.
     */
    setViewOptions(options = {}) {
        return options;
    }
    create_view(model, options = {}) {
        const id = (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.uuid)();
        const viewPromise = (model.state_change = model.state_change.then(async () => {
            const _view_name = model.get('_view_name');
            const _view_module = model.get('_view_module');
            try {
                const ViewType = (await this.loadViewClass(_view_name, _view_module, model.get('_view_module_version')));
                const view = new ViewType({
                    model: model,
                    options: this.setViewOptions(options),
                });
                view.listenTo(model, 'destroy', view.remove);
                await view.render();
                // This presumes the view is added to the list of model views below
                view.once('remove', () => {
                    if (model.views) {
                        delete model.views[id];
                    }
                });
                return view;
            }
            catch (e) {
                console.error(`Could not create a view for model id ${model.model_id}`);
                const msg = `Failed to create view for '${_view_name}' from module '${_view_module}' with model '${model.name}' from module '${model.module}'`;
                const ModelCls = _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.createErrorWidgetModel(e, msg);
                const errorModel = new ModelCls();
                const view = new _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.ErrorWidgetView({
                    model: errorModel,
                    options: this.setViewOptions(options),
                });
                await view.render();
                return view;
            }
        }));
        if (model.views) {
            model.views[id] = viewPromise;
        }
        return viewPromise;
    }
    /**
     * callback handlers specific to a view
     */
    callbacks(view) {
        return {};
    }
    /**
     * Get a promise for a model by model id.
     *
     * #### Notes
     * If the model is not found, the returned Promise object is rejected.
     *
     * If you would like to synchronously test if a model exists, use .has_model().
     */
    async get_model(model_id) {
        const modelPromise = this._models[model_id];
        if (modelPromise === undefined) {
            throw new Error('widget model not found');
        }
        return modelPromise;
    }
    /**
     * Returns true if the given model is registered, otherwise false.
     *
     * #### Notes
     * This is a synchronous way to check if a model is registered.
     */
    has_model(model_id) {
        return this._models[model_id] !== undefined;
    }
    /**
     * Handle when a comm is opened.
     */
    handle_comm_open(comm, msg) {
        const protocolVersion = (msg.metadata || {})['version'] || '';
        if (protocolVersion.split('.', 1)[0] !== PROTOCOL_MAJOR_VERSION) {
            const error = `Wrong widget protocol version: received protocol version '${protocolVersion}', but was expecting major version '${PROTOCOL_MAJOR_VERSION}'`;
            console.error(error);
            return Promise.reject(error);
        }
        const data = msg.content.data;
        const buffer_paths = data.buffer_paths || [];
        const buffers = msg.buffers || [];
        (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.put_buffers)(data.state, buffer_paths, buffers);
        return this.new_model({
            model_name: data.state['_model_name'],
            model_module: data.state['_model_module'],
            model_module_version: data.state['_model_module_version'],
            comm: comm,
        }, data.state).catch((0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.reject)('Could not create a model.', true));
    }
    /**
     * Create a comm and new widget model.
     * @param  options - same options as new_model but comm is not
     *                          required and additional options are available.
     * @param  serialized_state - serialized model attributes.
     */
    new_widget(options, serialized_state = {}) {
        let commPromise;
        // we check to make sure the view information is provided, to help catch
        // backwards incompatibility errors.
        if (options.view_name === undefined ||
            options.view_module === undefined ||
            options.view_module_version === undefined) {
            return Promise.reject('new_widget(...) must be given view information in the options.');
        }
        // If no comm is provided, a new comm is opened for the jupyter.widget
        // target.
        if (options.comm) {
            commPromise = Promise.resolve(options.comm);
        }
        else {
            commPromise = this._create_comm(this.comm_target_name, options.model_id, {
                state: {
                    _model_module: options.model_module,
                    _model_module_version: options.model_module_version,
                    _model_name: options.model_name,
                    _view_module: options.view_module,
                    _view_module_version: options.view_module_version,
                    _view_name: options.view_name,
                },
            }, { version: _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.PROTOCOL_VERSION });
        }
        // The options dictionary is copied since data will be added to it.
        const options_clone = Object.assign({}, options);
        // Create the model. In the case where the comm promise is rejected a
        // comm-less model is still created with the required model id.
        return commPromise.then((comm) => {
            // Comm Promise Resolved.
            options_clone.comm = comm;
            const widget_model = this.new_model(options_clone, serialized_state);
            return widget_model.then((model) => {
                model.sync('create', model);
                return model;
            });
        }, () => {
            // Comm Promise Rejected.
            if (!options_clone.model_id) {
                options_clone.model_id = (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.uuid)();
            }
            return this.new_model(options_clone, serialized_state);
        });
    }
    register_model(model_id, modelPromise) {
        this._models[model_id] = modelPromise;
        modelPromise.then((model) => {
            model.once('comm:close', () => {
                delete this._models[model_id];
            });
        });
    }
    /**
     * Create and return a promise for a new widget model
     *
     * @param options - the options for creating the model.
     * @param serialized_state - attribute values for the model.
     *
     * @example
     * widget_manager.new_model({
     *      model_name: 'IntSlider',
     *      model_module: '@jupyter-widgets/controls',
     *      model_module_version: '1.0.0',
     *      model_id: 'u-u-i-d'
     * }).then((model) => { console.log('Create success!', model); },
     *  (err) => {console.error(err)});
     *
     */
    async new_model(options, serialized_state = {}) {
        var _a, _b;
        const model_id = (_a = options.model_id) !== null && _a !== void 0 ? _a : (_b = options.comm) === null || _b === void 0 ? void 0 : _b.comm_id;
        if (!model_id) {
            throw new Error('Neither comm nor model_id provided in options object. At least one must exist.');
        }
        options.model_id = model_id;
        const modelPromise = this._make_model(options, serialized_state);
        // this call needs to happen before the first `await`, see note in `set_state`:
        this.register_model(model_id, modelPromise);
        return await modelPromise;
    }
    /**
     * Fetch all widgets states from the kernel using the control comm channel
     * If this fails (control comm handler not implemented kernel side),
     * it will fall back to `_loadFromKernelModels`.
     *
     * This is a utility function that can be used in subclasses.
     */
    async _loadFromKernel() {
        // Try fetching all widget states through the control comm
        let data;
        let buffers;
        try {
            const initComm = await this._create_comm(CONTROL_COMM_TARGET, (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.uuid)(), {}, { version: CONTROL_COMM_PROTOCOL_VERSION });
            await new Promise((resolve, reject) => {
                initComm.on_msg((msg) => {
                    data = msg['content']['data'];
                    if (data.method !== 'update_states') {
                        console.warn(`
              Unknown ${data.method} message on the Control channel
            `);
                        return;
                    }
                    buffers = (msg.buffers || []).map((b) => {
                        if (b instanceof DataView) {
                            return b;
                        }
                        else {
                            return new DataView(b instanceof ArrayBuffer ? b : b.buffer);
                        }
                    });
                    resolve(null);
                });
                initComm.on_close(() => reject('Control comm was closed too early'));
                // Send a states request msg
                initComm.send({ method: 'request_states' }, {});
                // Reject if we didn't get a response in time
                setTimeout(() => reject('Control comm did not respond in time'), CONTROL_COMM_TIMEOUT);
            });
            initComm.close();
        }
        catch (error) {
            // Fall back to the old implementation for old ipywidgets backend versions (ipywidgets<=7.6)
            return this._loadFromKernelModels();
        }
        const states = data.states;
        const bufferPaths = {};
        const bufferGroups = {};
        // Group buffers and buffer paths by widget id
        for (let i = 0; i < data.buffer_paths.length; i++) {
            const [widget_id, ...path] = data.buffer_paths[i];
            const b = buffers[i];
            if (!bufferPaths[widget_id]) {
                bufferPaths[widget_id] = [];
                bufferGroups[widget_id] = [];
            }
            bufferPaths[widget_id].push(path);
            bufferGroups[widget_id].push(b);
        }
        // Create comms for all new widgets.
        const widget_comms = await Promise.all(Object.keys(states).map(async (widget_id) => {
            const comm = this.has_model(widget_id)
                ? undefined
                : await this._create_comm('jupyter.widget', widget_id);
            return { widget_id, comm };
        }));
        await Promise.all(widget_comms.map(async ({ widget_id, comm }) => {
            const state = states[widget_id];
            // Put binary buffers
            if (widget_id in bufferPaths) {
                (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.put_buffers)(state, bufferPaths[widget_id], bufferGroups[widget_id]);
            }
            try {
                if (comm) {
                    // This must be the first await in the code path that
                    // reaches here so that registering the model promise in
                    // new_model can register the widget promise before it may
                    // be required by other widgets.
                    await this.new_model({
                        model_name: state.model_name,
                        model_module: state.model_module,
                        model_module_version: state.model_module_version,
                        model_id: widget_id,
                        comm: comm,
                    }, state.state);
                }
                else {
                    // model already exists here
                    const model = await this.get_model(widget_id);
                    const deserializedState = await model.constructor._deserialize_state(state.state, this);
                    model.set_state(deserializedState);
                }
            }
            catch (error) {
                // Failed to create a widget model, we continue creating other models so that
                // other widgets can render
                console.error(error);
            }
        }));
    }
    /**
     * Old implementation of fetching widget models one by one using
     * the request_state message on each comm.
     *
     * This is a utility function that can be used in subclasses.
     */
    async _loadFromKernelModels() {
        const comm_ids = await this._get_comm_info();
        // For each comm id that we do not know about, create the comm, and request the state.
        const widgets_info = await Promise.all(Object.keys(comm_ids).map(async (comm_id) => {
            if (this.has_model(comm_id)) {
                return;
            }
            const comm = await this._create_comm(this.comm_target_name, comm_id);
            let msg_id = '';
            const info = new _lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__.PromiseDelegate();
            comm.on_msg((msg) => {
                if (msg.parent_header.msg_id === msg_id &&
                    msg.header.msg_type === 'comm_msg' &&
                    msg.content.data.method === 'update') {
                    const data = msg.content.data;
                    const buffer_paths = data.buffer_paths || [];
                    const buffers = msg.buffers || [];
                    (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.put_buffers)(data.state, buffer_paths, buffers);
                    info.resolve({ comm, msg });
                }
            });
            msg_id = comm.send({
                method: 'request_state',
            }, this.callbacks(undefined));
            return info.promise;
        }));
        // We put in a synchronization barrier here so that we don't have to
        // topologically sort the restored widgets. `new_model` synchronously
        // registers the widget ids before reconstructing their state
        // asynchronously, so promises to every widget reference should be available
        // by the time they are used.
        await Promise.all(widgets_info.map(async (widget_info) => {
            if (!widget_info) {
                return;
            }
            const content = widget_info.msg.content;
            await this.new_model({
                model_name: content.data.state._model_name,
                model_module: content.data.state._model_module,
                model_module_version: content.data.state._model_module_version,
                comm: widget_info.comm,
            }, content.data.state);
        }));
    }
    async _make_model(options, serialized_state = {}) {
        const model_id = options.model_id;
        const model_promise = this.loadModelClass(options.model_name, options.model_module, options.model_module_version);
        let ModelType;
        const makeErrorModel = (error, msg) => {
            const Cls = _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.createErrorWidgetModel(error, msg);
            const widget_model = new Cls();
            return widget_model;
        };
        try {
            ModelType = await model_promise;
        }
        catch (error) {
            const msg = 'Could not instantiate widget';
            console.error(msg);
            return makeErrorModel(error, msg);
        }
        if (!ModelType) {
            const msg = 'Could not instantiate widget';
            console.error(msg);
            const error = new Error(`Cannot find model module ${options.model_module}@${options.model_module_version}, ${options.model_name}`);
            return makeErrorModel(error, msg);
        }
        let widget_model;
        try {
            const attributes = await ModelType._deserialize_state(serialized_state, this);
            const modelOptions = {
                widget_manager: this,
                model_id: model_id,
                comm: options.comm,
            };
            widget_model = new ModelType(attributes, modelOptions);
        }
        catch (error) {
            console.error(error);
            const msg = `Model class '${options.model_name}' from module '${options.model_module}' is loaded but can not be instantiated`;
            widget_model = makeErrorModel(error, msg);
        }
        widget_model.name = options.model_name;
        widget_model.module = options.model_module;
        return widget_model;
    }
    /**
     * Close all widgets and empty the widget state.
     * @return Promise that resolves when the widget state is cleared.
     */
    clear_state() {
        return (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.resolvePromisesDict)(this._models).then((models) => {
            Object.keys(models).forEach((id) => models[id].close());
            this._models = Object.create(null);
        });
    }
    /**
     * Asynchronously get the state of the widget manager.
     *
     * This includes all of the widget models, and follows the format given in
     * the @jupyter-widgets/schema package.
     *
     * @param options - The options for what state to return.
     * @returns Promise for a state dictionary
     */
    get_state(options = {}) {
        const modelPromises = Object.keys(this._models).map((id) => this._models[id]);
        return Promise.all(modelPromises).then((models) => {
            return serialize_state(models, options);
        });
    }
    /**
     * Set the widget manager state.
     *
     * @param state - a Javascript object conforming to the application/vnd.jupyter.widget-state+json spec.
     *
     * Reconstructs all of the widget models in the state, merges that with the
     * current manager state, and then attempts to redisplay the widgets in the
     * state.
     */
    set_state(state) {
        // Check to make sure that it's the same version we are parsing.
        if (!(state.version_major && state.version_major <= 2)) {
            throw 'Unsupported widget state format';
        }
        const models = state.state;
        // Recreate all the widget models for the given widget manager state.
        const all_models = this._get_comm_info().then((live_comms) => {
            /* Note: It is currently safe to just loop over the models in any order,
                     given that the following holds (does at the time of writing):
                     1: any call to `new_model` with state registers the model promise (e.g. with `register_model`)
                        synchronously (before it's first `await` statement).
                     2: any calls to a model constructor or the `set_state` method on a model,
                        happens asynchronously (in a `then` clause, or after an `await` statement).
      
                    Without these assumptions, one risks trying to set model state with a reference
                    to another model that doesn't exist yet!
                  */
            return Promise.all(Object.keys(models).map((model_id) => {
                // First put back the binary buffers
                const decode = {
                    base64: _utils__WEBPACK_IMPORTED_MODULE_2__.base64ToBuffer,
                    hex: _utils__WEBPACK_IMPORTED_MODULE_2__.hexToBuffer,
                };
                const model = models[model_id];
                const modelState = model.state;
                if (model.buffers) {
                    const bufferPaths = model.buffers.map((b) => b.path);
                    // put_buffers expects buffers to be DataViews
                    const buffers = model.buffers.map((b) => new DataView(decode[b.encoding](b.data)));
                    (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.put_buffers)(model.state, bufferPaths, buffers);
                }
                // If the model has already been created, set its state and then
                // return it.
                if (this.has_model(model_id)) {
                    return this.get_model(model_id).then((model) => {
                        // deserialize state
                        return model.constructor
                            ._deserialize_state(modelState || {}, this)
                            .then((attributes) => {
                            model.set_state(attributes); // case 2
                            return model;
                        });
                    });
                }
                const modelCreate = {
                    model_id: model_id,
                    model_name: model.model_name,
                    model_module: model.model_module,
                    model_module_version: model.model_module_version,
                };
                if (Object.prototype.hasOwnProperty.call(live_comms, 'model_id')) {
                    // live comm
                    // This connects to an existing comm if it exists, and
                    // should *not* send a comm open message.
                    return this._create_comm(this.comm_target_name, model_id).then((comm) => {
                        modelCreate.comm = comm;
                        return this.new_model(modelCreate); // No state, so safe wrt. case 1
                    });
                }
                else {
                    return this.new_model(modelCreate, modelState); // case 1
                }
            }));
        });
        return all_models;
    }
    /**
     * Disconnect the widget manager from the kernel, setting each model's comm
     * as dead.
     */
    disconnect() {
        Object.keys(this._models).forEach((i) => {
            this._models[i].then((model) => {
                model.comm_live = false;
            });
        });
    }
    /**
     * Resolve a URL relative to the current notebook location.
     *
     * The default implementation just returns the original url.
     */
    resolveUrl(url) {
        return Promise.resolve(url);
    }
    inline_sanitize(source) {
        const parts = (0,_latex__WEBPACK_IMPORTED_MODULE_3__.removeMath)(source);
        // Sanitize tags for inline output.
        const sanitized = default_inline_sanitize(parts['text']);
        return (0,_latex__WEBPACK_IMPORTED_MODULE_3__.replaceMath)(sanitized, parts['math']);
    }
    async loadModelClass(className, moduleName, moduleVersion) {
        try {
            const promise = this.loadClass(className, moduleName, moduleVersion);
            await promise;
            return promise;
        }
        catch (error) {
            console.error(error);
            const msg = `Failed to load model class '${className}' from module '${moduleName}'`;
            return _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.createErrorWidgetModel(error, msg);
        }
    }
    async loadViewClass(className, moduleName, moduleVersion) {
        try {
            const promise = this.loadClass(className, moduleName, moduleVersion);
            await promise;
            return promise;
        }
        catch (error) {
            console.error(error);
            const msg = `Failed to load view class '${className}' from module '${moduleName}'`;
            return _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.createErrorWidgetView(error, msg);
        }
    }
    /**
     * Filter serialized widget state to remove any ID's already present in manager.
     *
     * @param {*} state Serialized state to filter
     *
     * @returns {*} A copy of the state, with its 'state' attribute filtered
     */
    filterExistingModelState(serialized_state) {
        let models = serialized_state.state;
        models = Object.keys(models)
            .filter((model_id) => !this.has_model(model_id))
            .reduce((res, model_id) => {
            res[model_id] = models[model_id];
            return res;
        }, {});
        return Object.assign(Object.assign({}, serialized_state), { state: models });
    }
}
/**
 * Serialize an array of widget models
 *
 * #### Notes
 * The return value follows the format given in the
 * @jupyter-widgets/schema package.
 */
function serialize_state(models, options = {}) {
    const state = {};
    models.forEach((model) => {
        const model_id = model.model_id;
        const split = (0,_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.remove_buffers)(model.serialize(model.get_state(options.drop_defaults)));
        const buffers = split.buffers.map((buffer, index) => {
            return {
                data: (0,_utils__WEBPACK_IMPORTED_MODULE_2__.bufferToBase64)(buffer),
                path: split.buffer_paths[index],
                encoding: 'base64',
            };
        });
        state[model_id] = {
            model_name: model.name,
            model_module: model.module,
            model_module_version: model.get('_model_module_version'),
            state: split.state,
        };
        // To save space, only include the buffers key if we have buffers
        if (buffers.length > 0) {
            state[model_id].buffers = buffers;
        }
    });
    return { version_major: 2, version_minor: 0, state: state };
}


/***/ }),

/***/ "../../packages/base-manager/lib/utils.js":
/*!************************************************!*\
  !*** ../../packages/base-manager/lib/utils.js ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   base64ToBuffer: () => (/* binding */ base64ToBuffer),
/* harmony export */   bufferToBase64: () => (/* binding */ bufferToBase64),
/* harmony export */   bufferToHex: () => (/* binding */ bufferToHex),
/* harmony export */   hexToBuffer: () => (/* binding */ hexToBuffer)
/* harmony export */ });
/* harmony import */ var base64_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! base64-js */ "../../node_modules/base64-js/index.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

const hexTable = [
    '00',
    '01',
    '02',
    '03',
    '04',
    '05',
    '06',
    '07',
    '08',
    '09',
    '0A',
    '0B',
    '0C',
    '0D',
    '0E',
    '0F',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18',
    '19',
    '1A',
    '1B',
    '1C',
    '1D',
    '1E',
    '1F',
    '20',
    '21',
    '22',
    '23',
    '24',
    '25',
    '26',
    '27',
    '28',
    '29',
    '2A',
    '2B',
    '2C',
    '2D',
    '2E',
    '2F',
    '30',
    '31',
    '32',
    '33',
    '34',
    '35',
    '36',
    '37',
    '38',
    '39',
    '3A',
    '3B',
    '3C',
    '3D',
    '3E',
    '3F',
    '40',
    '41',
    '42',
    '43',
    '44',
    '45',
    '46',
    '47',
    '48',
    '49',
    '4A',
    '4B',
    '4C',
    '4D',
    '4E',
    '4F',
    '50',
    '51',
    '52',
    '53',
    '54',
    '55',
    '56',
    '57',
    '58',
    '59',
    '5A',
    '5B',
    '5C',
    '5D',
    '5E',
    '5F',
    '60',
    '61',
    '62',
    '63',
    '64',
    '65',
    '66',
    '67',
    '68',
    '69',
    '6A',
    '6B',
    '6C',
    '6D',
    '6E',
    '6F',
    '70',
    '71',
    '72',
    '73',
    '74',
    '75',
    '76',
    '77',
    '78',
    '79',
    '7A',
    '7B',
    '7C',
    '7D',
    '7E',
    '7F',
    '80',
    '81',
    '82',
    '83',
    '84',
    '85',
    '86',
    '87',
    '88',
    '89',
    '8A',
    '8B',
    '8C',
    '8D',
    '8E',
    '8F',
    '90',
    '91',
    '92',
    '93',
    '94',
    '95',
    '96',
    '97',
    '98',
    '99',
    '9A',
    '9B',
    '9C',
    '9D',
    '9E',
    '9F',
    'A0',
    'A1',
    'A2',
    'A3',
    'A4',
    'A5',
    'A6',
    'A7',
    'A8',
    'A9',
    'AA',
    'AB',
    'AC',
    'AD',
    'AE',
    'AF',
    'B0',
    'B1',
    'B2',
    'B3',
    'B4',
    'B5',
    'B6',
    'B7',
    'B8',
    'B9',
    'BA',
    'BB',
    'BC',
    'BD',
    'BE',
    'BF',
    'C0',
    'C1',
    'C2',
    'C3',
    'C4',
    'C5',
    'C6',
    'C7',
    'C8',
    'C9',
    'CA',
    'CB',
    'CC',
    'CD',
    'CE',
    'CF',
    'D0',
    'D1',
    'D2',
    'D3',
    'D4',
    'D5',
    'D6',
    'D7',
    'D8',
    'D9',
    'DA',
    'DB',
    'DC',
    'DD',
    'DE',
    'DF',
    'E0',
    'E1',
    'E2',
    'E3',
    'E4',
    'E5',
    'E6',
    'E7',
    'E8',
    'E9',
    'EA',
    'EB',
    'EC',
    'ED',
    'EE',
    'EF',
    'F0',
    'F1',
    'F2',
    'F3',
    'F4',
    'F5',
    'F6',
    'F7',
    'F8',
    'F9',
    'FA',
    'FB',
    'FC',
    'FD',
    'FE',
    'FF',
];
/**
 * Convert an ArrayBuffer to a hex string.
 */
function bufferToHex(buffer) {
    const x = new Uint8Array(buffer);
    const s = [];
    for (let i = 0; i < x.length; i++) {
        s.push(hexTable[x[i]]);
    }
    return s.join('');
}
/**
 * Convert a hex string to an ArrayBuffer.
 */
function hexToBuffer(hex) {
    const x = new Uint8Array(hex.length / 2);
    for (let i = 0; i < hex.length; i += 2) {
        x[i / 2] = parseInt(hex.slice(i, i + 2), 16);
    }
    return x.buffer;
}
/**
 * Convert an ArrayBuffer to a base64 string.
 */
function bufferToBase64(buffer) {
    return (0,base64_js__WEBPACK_IMPORTED_MODULE_0__.fromByteArray)(new Uint8Array(buffer));
}
/**
 * Convert a base64 string to an ArrayBuffer.
 */
function base64ToBuffer(base64) {
    return (0,base64_js__WEBPACK_IMPORTED_MODULE_0__.toByteArray)(base64).buffer;
}


/***/ }),

/***/ "?c11c":
/*!**************************************!*\
  !*** ./terminal-highlight (ignored) ***!
  \**************************************/
/***/ (() => {

/* (ignored) */

/***/ }),

/***/ "?783f":
/*!********************!*\
  !*** fs (ignored) ***!
  \********************/
/***/ (() => {

/* (ignored) */

/***/ }),

/***/ "?c929":
/*!**********************!*\
  !*** path (ignored) ***!
  \**********************/
/***/ (() => {

/* (ignored) */

/***/ }),

/***/ "?4a0c":
/*!*******************************!*\
  !*** source-map-js (ignored) ***!
  \*******************************/
/***/ (() => {

/* (ignored) */

/***/ }),

/***/ "?e350":
/*!*********************!*\
  !*** url (ignored) ***!
  \*********************/
/***/ (() => {

/* (ignored) */

/***/ })

}]);
//# sourceMappingURL=packages_base-manager_lib_index_js.13941b475e5b5a1ab133.js.map