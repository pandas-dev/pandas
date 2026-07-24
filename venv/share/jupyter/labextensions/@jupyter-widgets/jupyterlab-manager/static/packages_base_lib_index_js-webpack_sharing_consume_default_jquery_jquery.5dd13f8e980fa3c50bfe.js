"use strict";
(self["webpackChunk_jupyter_widgets_jupyterlab_manager"] = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] || []).push([["packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery"],{

/***/ "../../packages/base/lib/backbone-patch.js":
/*!*************************************************!*\
  !*** ../../packages/base/lib/backbone-patch.js ***!
  \*************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   set: () => (/* binding */ set)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__);
// This file contains a modified version of the set function from the Backbone
// (see
// https://github.com/jashkenas/backbone/blob/05fde9e201f7e2137796663081105cd6dad12a98/backbone.js#L460,
// with changes below marked with an EDIT comment). This file in Backbone has the following license.
//     (c) 2010-2015 Jeremy Ashkenas, DocumentCloud and Investigative Reporters & Editors
//     Backbone may be freely distributed under the MIT license.
//     For all details and documentation:
//     http://backbonejs.org
// Backbone's full license is below (from https://github.com/jashkenas/backbone/blob/05fde9e201f7e2137796663081105cd6dad12a98/LICENSE)
/*
Copyright (c) 2010-2015 Jeremy Ashkenas, DocumentCloud

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/


// Set a hash of model attributes on the object, firing `"change"`. This is
// the core primitive operation of a model, updating the data and notifying
// anyone who needs to know about the change in state. The heart of the beast.
// This *MUST* be called with the model as the `this` context.
function set(key, val, options) {
    if (key == null) {
        return this;
    }
    // Handle both `"key", value` and `{key: value}` -style arguments.
    let attrs;
    if (_lumino_coreutils__WEBPACK_IMPORTED_MODULE_1__.JSONExt.isObject(key)) {
        attrs = key;
        options = val;
    }
    else {
        (attrs = {})[key] = val;
    }
    options || (options = {});
    // Run validation.
    if (!this._validate(attrs, options)) {
        return false;
    }
    // Extract attributes and options.
    const unset = options.unset;
    const silent = options.silent;
    const changes = [];
    const changing = this._changing;
    this._changing = true;
    try {
        if (!changing) {
            // EDIT: changed to use object spread instead of _.clone
            this._previousAttributes = Object.assign({}, this.attributes);
            this.changed = {};
        }
        const current = this.attributes;
        const changed = this.changed;
        const prev = this._previousAttributes;
        // For each `set` attribute, update or delete the current value.
        for (const attr in attrs) {
            val = attrs[attr];
            // EDIT: the following two lines use our isEqual instead of _.isEqual
            if (!_utils__WEBPACK_IMPORTED_MODULE_0__.isEqual(current[attr], val)) {
                changes.push(attr);
            }
            if (!_utils__WEBPACK_IMPORTED_MODULE_0__.isEqual(prev[attr], val)) {
                changed[attr] = val;
            }
            else {
                delete changed[attr];
            }
            unset ? delete current[attr] : (current[attr] = val);
        }
        // Update the `id`.
        this.id = this.get(this.idAttribute);
        // Trigger all relevant attribute changes.
        if (!silent) {
            if (changes.length) {
                this._pending = options;
            }
            for (let i = 0; i < changes.length; i++) {
                this.trigger('change:' + changes[i], this, current[changes[i]], options);
            }
        }
        // You might be wondering why there's a `while` loop here. Changes can
        // be recursively nested within `"change"` events.
        if (changing) {
            return this;
        }
        if (!silent) {
            while (this._pending) {
                options = this._pending;
                this._pending = false;
                this.trigger('change', this, options);
            }
        }
    }
    finally {
        this._pending = false;
        this._changing = false;
    }
    return this;
}


/***/ }),

/***/ "../../packages/base/lib/errorwidget.js":
/*!**********************************************!*\
  !*** ../../packages/base/lib/errorwidget.js ***!
  \**********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ErrorWidgetView: () => (/* binding */ ErrorWidgetView),
/* harmony export */   createErrorWidgetModel: () => (/* binding */ createErrorWidgetModel),
/* harmony export */   createErrorWidgetView: () => (/* binding */ createErrorWidgetView)
/* harmony export */ });
/* harmony import */ var _widget__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget */ "../../packages/base/lib/widget.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./version */ "../../packages/base/lib/version.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");



// create a Widget Model that captures an error object
function createErrorWidgetModel(error, msg) {
    class ErrorWidget extends _widget__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetModel {
        constructor(attributes, options) {
            attributes = Object.assign(Object.assign({}, attributes), { _view_name: 'ErrorWidgetView', _view_module: '@jupyter-widgets/base', _model_module_version: _version__WEBPACK_IMPORTED_MODULE_1__.JUPYTER_WIDGETS_VERSION, _view_module_version: _version__WEBPACK_IMPORTED_MODULE_1__.JUPYTER_WIDGETS_VERSION, msg: msg, error: error });
            super(attributes, options);
            this.comm_live = true;
        }
    }
    return ErrorWidget;
}
class ErrorWidgetView extends _widget__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
    generateErrorMessage() {
        return {
            msg: this.model.get('msg'),
            stack: String(this.model.get('error').stack),
        };
    }
    render() {
        const { msg, stack } = this.generateErrorMessage();
        this.el.classList.add('jupyter-widgets');
        const content = document.createElement('div');
        content.classList.add('jupyter-widgets-error-widget', 'icon-error');
        content.innerHTML = _utils__WEBPACK_IMPORTED_MODULE_2__.BROKEN_FILE_SVG_ICON;
        const text = document.createElement('pre');
        text.style.textAlign = 'center';
        text.innerText = 'Click to show javascript error.';
        content.append(text);
        this.el.appendChild(content);
        let width;
        let height;
        this.el.onclick = () => {
            if (content.classList.contains('icon-error')) {
                height = height || content.clientHeight;
                width = width || content.clientWidth;
                content.classList.remove('icon-error');
                content.innerHTML = `
        <pre>[Open Browser Console for more detailed log - Double click to close this message]\n${msg}\n${stack}</pre>
        `;
                content.style.height = `${height}px`;
                content.style.width = `${width}px`;
                content.classList.add('text-error');
            }
        };
        this.el.ondblclick = () => {
            if (content.classList.contains('text-error')) {
                content.classList.remove('text-error');
                content.innerHTML = _utils__WEBPACK_IMPORTED_MODULE_2__.BROKEN_FILE_SVG_ICON;
                content.append(text);
                content.classList.add('icon-error');
            }
        };
    }
}
function createErrorWidgetView(error, msg) {
    return class InnerErrorWidgetView extends ErrorWidgetView {
        generateErrorMessage() {
            return {
                msg,
                stack: String(error instanceof Error ? error.stack : error),
            };
        }
    };
}


/***/ }),

/***/ "../../packages/base/lib/index.js":
/*!****************************************!*\
  !*** ../../packages/base/lib/index.js ***!
  \****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BROKEN_FILE_SVG_ICON: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.BROKEN_FILE_SVG_ICON),
/* harmony export */   DOMWidgetModel: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetModel),
/* harmony export */   DOMWidgetView: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView),
/* harmony export */   ErrorWidgetView: () => (/* reexport safe */ _errorwidget__WEBPACK_IMPORTED_MODULE_9__.ErrorWidgetView),
/* harmony export */   IJupyterWidgetRegistry: () => (/* reexport safe */ _registry__WEBPACK_IMPORTED_MODULE_8__.IJupyterWidgetRegistry),
/* harmony export */   JUPYTER_WIDGETS_VERSION: () => (/* reexport safe */ _version__WEBPACK_IMPORTED_MODULE_6__.JUPYTER_WIDGETS_VERSION),
/* harmony export */   JupyterLuminoPanelWidget: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.JupyterLuminoPanelWidget),
/* harmony export */   JupyterLuminoWidget: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.JupyterLuminoWidget),
/* harmony export */   JupyterPhosphorPanelWidget: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.JupyterPhosphorPanelWidget),
/* harmony export */   JupyterPhosphorWidget: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.JupyterPhosphorWidget),
/* harmony export */   LayoutModel: () => (/* reexport safe */ _widget_layout__WEBPACK_IMPORTED_MODULE_2__.LayoutModel),
/* harmony export */   LayoutView: () => (/* reexport safe */ _widget_layout__WEBPACK_IMPORTED_MODULE_2__.LayoutView),
/* harmony export */   PROTOCOL_VERSION: () => (/* reexport safe */ _version__WEBPACK_IMPORTED_MODULE_6__.PROTOCOL_VERSION),
/* harmony export */   StyleModel: () => (/* reexport safe */ _widget_style__WEBPACK_IMPORTED_MODULE_3__.StyleModel),
/* harmony export */   StyleView: () => (/* reexport safe */ _widget_style__WEBPACK_IMPORTED_MODULE_3__.StyleView),
/* harmony export */   ViewList: () => (/* reexport safe */ _viewlist__WEBPACK_IMPORTED_MODULE_5__.ViewList),
/* harmony export */   WidgetModel: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.WidgetModel),
/* harmony export */   WidgetView: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.WidgetView),
/* harmony export */   assign: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.assign),
/* harmony export */   createErrorWidgetModel: () => (/* reexport safe */ _errorwidget__WEBPACK_IMPORTED_MODULE_9__.createErrorWidgetModel),
/* harmony export */   createErrorWidgetView: () => (/* reexport safe */ _errorwidget__WEBPACK_IMPORTED_MODULE_9__.createErrorWidgetView),
/* harmony export */   difference: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.difference),
/* harmony export */   isEqual: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.isEqual),
/* harmony export */   isObject: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.isObject),
/* harmony export */   isSerializable: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.isSerializable),
/* harmony export */   pack_models: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.pack_models),
/* harmony export */   put_buffers: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.put_buffers),
/* harmony export */   reject: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.reject),
/* harmony export */   remove_buffers: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.remove_buffers),
/* harmony export */   resolvePromisesDict: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.resolvePromisesDict),
/* harmony export */   shims: () => (/* reexport safe */ _services_shim__WEBPACK_IMPORTED_MODULE_4__.shims),
/* harmony export */   unpack_models: () => (/* reexport safe */ _widget__WEBPACK_IMPORTED_MODULE_0__.unpack_models),
/* harmony export */   uuid: () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_7__.uuid)
/* harmony export */ });
/* harmony import */ var _widget__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./widget */ "../../packages/base/lib/widget.js");
/* harmony import */ var _manager__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./manager */ "../../packages/base/lib/manager.js");
/* harmony import */ var _widget_layout__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./widget_layout */ "../../packages/base/lib/widget_layout.js");
/* harmony import */ var _widget_style__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./widget_style */ "../../packages/base/lib/widget_style.js");
/* harmony import */ var _services_shim__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./services-shim */ "../../packages/base/lib/services-shim.js");
/* harmony import */ var _viewlist__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./viewlist */ "../../packages/base/lib/viewlist.js");
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./version */ "../../packages/base/lib/version.js");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");
/* harmony import */ var _registry__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./registry */ "../../packages/base/lib/registry.js");
/* harmony import */ var _errorwidget__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./errorwidget */ "../../packages/base/lib/errorwidget.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.












/***/ }),

/***/ "../../packages/base/lib/manager.js":
/*!******************************************!*\
  !*** ../../packages/base/lib/manager.js ***!
  \******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.



/***/ }),

/***/ "../../packages/base/lib/nativeview.js":
/*!*********************************************!*\
  !*** ../../packages/base/lib/nativeview.js ***!
  \*********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   NativeView: () => (/* binding */ NativeView)
/* harmony export */ });
/* harmony import */ var backbone__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! backbone */ "../../node_modules/backbone/backbone.js");
/* harmony import */ var backbone__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(backbone__WEBPACK_IMPORTED_MODULE_0__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
/*
 This file contains substantial portions of https://github.com/akre54/Backbone.NativeView/blob/521188d9554b53d95d70ed34f878d8ac9fc10df2/backbone.nativeview.js, which has the following license:

(c) 2015 Adam Krebs, Jimmy Yuen Ho Wong
Backbone.NativeView may be freely distributed under the MIT license.

Copyright (c) 2014 Adam Krebs

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

*/

// Caches a local reference to `Element.prototype` for faster access.
const ElementProto = typeof Element !== 'undefined' ? Element.prototype : undefined;
// Find the right `Element#matches` for IE>=9 and modern browsers.
function matchesFallback(selector) {
    const matches = (this.document || this.ownerDocument).querySelectorAll(selector);
    let i = matches.length;
    while (--i >= 0 && matches.item(i) !== this) {
        continue;
    }
    return i > -1;
}
const matchesSelector = ElementProto
    ? ElementProto.matches ||
        ElementProto['webkitMatchesSelector'] ||
        ElementProto['mozMatchesSelector'] ||
        ElementProto['msMatchesSelector'] ||
        ElementProto['oMatchesSelector'] ||
        matchesFallback
    : matchesFallback;
class NativeView extends backbone__WEBPACK_IMPORTED_MODULE_0__.View {
    _removeElement() {
        this.undelegateEvents();
        if (this.el.parentNode) {
            this.el.parentNode.removeChild(this.el);
        }
    }
    // Apply the `element` to the view.
    _setElement(element) {
        this.el = element;
    }
    // Set a hash of attributes to the view's `el`. We use the "prop" version
    // if available, falling back to `setAttribute` for the catch-all.
    _setAttributes(attrs) {
        for (const attr in attrs) {
            attr in this.el
                ? (this.el[attr] = attrs[attr])
                : this.el.setAttribute(attr, attrs[attr]);
        }
    }
    delegate(eventName, selector, listener) {
        if (typeof selector !== 'string') {
            listener = selector;
            selector = null;
        }
        // We have to initialize this here, instead of in the constructor, because the
        // super constructor eventually calls this method before we get a chance to initialize
        // this._domEvents to an empty list.
        if (this._domEvents === void 0) {
            this._domEvents = [];
        }
        const root = this.el;
        const handler = selector
            ? function (e) {
                let node = e.target || e.srcElement;
                for (; node && node !== root; node = node.parentNode) {
                    if (matchesSelector.call(node, selector)) {
                        e.delegateTarget = node;
                        if (listener.handleEvent) {
                            return listener.handleEvent(e);
                        }
                        else {
                            return listener(e);
                        }
                    }
                }
            }
            : listener;
        this.el.addEventListener(eventName, handler, false);
        this._domEvents.push({ eventName, handler, listener, selector });
        return handler;
    }
    undelegate(eventName, selector, listener) {
        if (typeof selector === 'function') {
            listener = selector;
            selector = null;
        }
        if (this.el && this._domEvents) {
            const handlers = this._domEvents.slice();
            let i = handlers.length;
            while (i--) {
                const item = handlers[i];
                const match = item.eventName === eventName &&
                    (listener ? item.listener === listener : true) &&
                    (selector ? item.selector === selector : true);
                if (!match) {
                    continue;
                }
                this.el.removeEventListener(item.eventName, item.handler, false);
                this._domEvents.splice(i, 1);
            }
        }
        return this;
    }
    // Remove all events created with `delegate` from `el`
    undelegateEvents() {
        if (this.el && this._domEvents) {
            const len = this._domEvents.length;
            for (let i = 0; i < len; i++) {
                const item = this._domEvents[i];
                this.el.removeEventListener(item.eventName, item.handler, false);
            }
            this._domEvents.length = 0;
        }
        return this;
    }
}


/***/ }),

/***/ "../../packages/base/lib/registry.js":
/*!*******************************************!*\
  !*** ../../packages/base/lib/registry.js ***!
  \*******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   IJupyterWidgetRegistry: () => (/* binding */ IJupyterWidgetRegistry)
/* harmony export */ });
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

/**
 * A runtime interface token for a widget registry.
 */
const IJupyterWidgetRegistry = new _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__.Token('jupyter.extensions.jupyterWidgetRegistry');


/***/ }),

/***/ "../../packages/base/lib/services-shim.js":
/*!************************************************!*\
  !*** ../../packages/base/lib/services-shim.js ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   shims: () => (/* binding */ shims)
/* harmony export */ });
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
var shims;
(function (shims) {
    let services;
    (function (services) {
        /**
         * Public constructor
         * @param jsServicesKernel - @jupyterlab/services Kernel.IKernel instance
         */
        class CommManager {
            constructor(jsServicesKernel) {
                this.targets = Object.create(null);
                this.comms = Object.create(null);
                this.init_kernel(jsServicesKernel);
            }
            /**
             * Hookup kernel events.
             * @param  {Kernel.IKernel} jsServicesKernel - @jupyterlab/services Kernel.IKernel instance
             */
            init_kernel(jsServicesKernel) {
                this.kernel = jsServicesKernel; // These aren't really the same.
                this.jsServicesKernel = jsServicesKernel;
            }
            /**
             * Creates a new connected comm
             */
            async new_comm(target_name, data, callbacks, metadata, comm_id, buffers) {
                const c = this.jsServicesKernel.createComm(target_name, comm_id);
                const comm = new Comm(c);
                this.register_comm(comm);
                comm.open(data, callbacks, metadata, buffers);
                return comm;
            }
            /**
             * Register a comm target
             * @param  {string} target_name
             * @param  {(Comm, object) => void} f - callback that is called when the
             *                         comm is made.  Signature of f(comm, msg).
             */
            register_target(target_name, f) {
                const handle = this.jsServicesKernel.registerCommTarget(target_name, (jsServicesComm, msg) => {
                    // Create the comm.
                    const comm = new Comm(jsServicesComm);
                    this.register_comm(comm);
                    // Call the callback for the comm.
                    try {
                        return f(comm, msg);
                    }
                    catch (e) {
                        comm.close();
                        console.error(e);
                        console.error(new Error('Exception opening new comm'));
                    }
                });
                this.targets[target_name] = handle;
            }
            /**
             * Unregisters a comm target
             * @param  {string} target_name
             */
            unregister_target(target_name, f) {
                const handle = this.targets[target_name];
                handle.dispose();
                delete this.targets[target_name];
            }
            /**
             * Register a comm in the mapping
             */
            register_comm(comm) {
                this.comms[comm.comm_id] = Promise.resolve(comm);
                comm.kernel = this.kernel;
                return comm.comm_id;
            }
        }
        services.CommManager = CommManager;
        /**
         * Public constructor
         * @param  {IComm} jsServicesComm - @jupyterlab/services IComm instance
         */
        class Comm {
            constructor(jsServicesComm) {
                this.jsServicesComm = jsServicesComm;
            }
            /**
             * Comm id
             * @return {string}
             */
            get comm_id() {
                return this.jsServicesComm.commId;
            }
            /**
             * Target name
             * @return {string}
             */
            get target_name() {
                return this.jsServicesComm.targetName;
            }
            /**
             * Opens a sibling comm in the backend
             * @param  data
             * @param  callbacks
             * @param  metadata
             * @return msg id
             */
            open(data, callbacks, metadata, buffers) {
                const future = this.jsServicesComm.open(data, metadata, buffers);
                this._hookupCallbacks(future, callbacks);
                return future.msg.header.msg_id;
            }
            /**
             * Sends a message to the sibling comm in the backend
             * @param  data
             * @param  callbacks
             * @param  metadata
             * @param  buffers
             * @return message id
             */
            send(data, callbacks, metadata, buffers) {
                const future = this.jsServicesComm.send(data, metadata, buffers);
                this._hookupCallbacks(future, callbacks);
                return future.msg.header.msg_id;
            }
            /**
             * Closes the sibling comm in the backend
             * @param  data
             * @param  callbacks
             * @param  metadata
             * @return msg id
             */
            close(data, callbacks, metadata, buffers) {
                const future = this.jsServicesComm.close(data, metadata, buffers);
                this._hookupCallbacks(future, callbacks);
                return future.msg.header.msg_id;
            }
            /**
             * Register a message handler
             * @param  callback, which is given a message
             */
            on_msg(callback) {
                this.jsServicesComm.onMsg = callback.bind(this);
            }
            /**
             * Register a handler for when the comm is closed by the backend
             * @param  callback, which is given a message
             */
            on_close(callback) {
                this.jsServicesComm.onClose = callback.bind(this);
            }
            /**
             * Hooks callback object up with @jupyterlab/services IKernelFuture
             * @param  @jupyterlab/services IKernelFuture instance
             * @param  callbacks
             */
            _hookupCallbacks(future, callbacks) {
                if (callbacks) {
                    future.onReply = function (msg) {
                        if (callbacks.shell && callbacks.shell.reply) {
                            callbacks.shell.reply(msg);
                        }
                    };
                    future.onStdin = function (msg) {
                        if (callbacks.input) {
                            callbacks.input(msg);
                        }
                    };
                    future.onIOPub = function (msg) {
                        if (callbacks.iopub) {
                            if (callbacks.iopub.status && msg.header.msg_type === 'status') {
                                callbacks.iopub.status(msg);
                            }
                            else if (callbacks.iopub.clear_output &&
                                msg.header.msg_type === 'clear_output') {
                                callbacks.iopub.clear_output(msg);
                            }
                            else if (callbacks.iopub.output) {
                                switch (msg.header.msg_type) {
                                    case 'display_data':
                                    case 'execute_result':
                                    case 'stream':
                                    case 'error':
                                        callbacks.iopub.output(msg);
                                        break;
                                    default:
                                        break;
                                }
                            }
                        }
                    };
                }
            }
        }
        services.Comm = Comm;
    })(services = shims.services || (shims.services = {}));
})(shims || (shims = {}));


/***/ }),

/***/ "../../packages/base/lib/utils.js":
/*!****************************************!*\
  !*** ../../packages/base/lib/utils.js ***!
  \****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   BROKEN_FILE_SVG_ICON: () => (/* binding */ BROKEN_FILE_SVG_ICON),
/* harmony export */   assign: () => (/* binding */ assign),
/* harmony export */   difference: () => (/* binding */ difference),
/* harmony export */   isEqual: () => (/* binding */ isEqual),
/* harmony export */   isObject: () => (/* binding */ isObject),
/* harmony export */   isSerializable: () => (/* binding */ isSerializable),
/* harmony export */   put_buffers: () => (/* binding */ put_buffers),
/* harmony export */   reject: () => (/* binding */ reject),
/* harmony export */   remove_buffers: () => (/* binding */ remove_buffers),
/* harmony export */   resolvePromisesDict: () => (/* binding */ resolvePromisesDict),
/* harmony export */   uuid: () => (/* binding */ uuid)
/* harmony export */ });
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var lodash_isEqual__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! lodash/isEqual */ "../../node_modules/lodash/isEqual.js");
/* harmony import */ var lodash_isEqual__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(lodash_isEqual__WEBPACK_IMPORTED_MODULE_1__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


/**
 * Find all strings in the first argument that are not in the second.
 */
function difference(a, b) {
    return a.filter((v) => b.indexOf(v) === -1);
}
/**
 * Compare two objects deeply to see if they are equal.
 */
function isEqual(a, b) {
    return lodash_isEqual__WEBPACK_IMPORTED_MODULE_1___default()(a, b);
}
/**
 * A polyfill for Object.assign
 *
 * This is from code that Typescript 2.4 generates for a polyfill.
 */
const assign = Object.assign ||
    function (t, ...args) {
        for (let i = 1; i < args.length; i++) {
            const s = args[i];
            for (const p in s) {
                if (Object.prototype.hasOwnProperty.call(s, p)) {
                    t[p] = s[p];
                }
            }
        }
        return t;
    };
/**
 * Generate a UUID
 *
 * http://www.ietf.org/rfc/rfc4122.txt
 */
function uuid() {
    return _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__.UUID.uuid4();
}
/**
 * Resolve a promiseful dictionary.
 * Returns a single Promise.
 */
function resolvePromisesDict(d) {
    const keys = Object.keys(d);
    const values = [];
    keys.forEach(function (key) {
        values.push(d[key]);
    });
    return Promise.all(values).then((v) => {
        const d = {};
        for (let i = 0; i < keys.length; i++) {
            d[keys[i]] = v[i];
        }
        return d;
    });
}
/**
 * Creates a wrappable Promise rejection function.
 *
 * Creates a function that logs an error message before rethrowing
 * the original error that caused the promise to reject.
 */
function reject(message, log) {
    return function promiseRejection(error) {
        if (log) {
            console.error(new Error(message));
        }
        throw error;
    };
}
/**
 * Takes an object 'state' and fills in buffer[i] at 'path' buffer_paths[i]
 * where buffer_paths[i] is a list indicating where in the object buffer[i] should
 * be placed
 * Example: state = {a: 1, b: {}, c: [0, null]}
 * buffers = [array1, array2]
 * buffer_paths = [['b', 'data'], ['c', 1]]
 * Will lead to {a: 1, b: {data: array1}, c: [0, array2]}
 */
function put_buffers(state, buffer_paths, buffers) {
    for (let i = 0; i < buffer_paths.length; i++) {
        const buffer_path = buffer_paths[i];
        // make sure the buffers are DataViews
        let buffer = buffers[i];
        if (!(buffer instanceof DataView)) {
            buffer = new DataView(buffer instanceof ArrayBuffer ? buffer : buffer.buffer);
        }
        // say we want to set state[x][y][z] = buffer
        let obj = state;
        // we first get obj = state[x][y]
        for (let j = 0; j < buffer_path.length - 1; j++) {
            obj = obj[buffer_path[j]];
        }
        // and then set: obj[z] = buffer
        obj[buffer_path[buffer_path.length - 1]] = buffer;
    }
}
function isSerializable(object) {
    var _a;
    return (_a = (typeof object === 'object' && object && 'toJSON' in object)) !== null && _a !== void 0 ? _a : false;
}
function isObject(data) {
    return _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__.JSONExt.isObject(data);
}
/**
 * The inverse of put_buffers, return an objects with the new state where all buffers(ArrayBuffer)
 * are removed. If a buffer is a member of an object, that object is cloned, and the key removed. If a buffer
 * is an element of an array, that array is cloned, and the element is set to null.
 * See put_buffers for the meaning of buffer_paths
 * Returns an object with the new state (.state) an array with paths to the buffers (.buffer_paths),
 * and the buffers associated to those paths (.buffers).
 */
function remove_buffers(state) {
    const buffers = [];
    const buffer_paths = [];
    // if we need to remove an object from a list, we need to clone that list, otherwise we may modify
    // the internal state of the widget model
    // however, we do not want to clone everything, for performance
    function remove(obj, path) {
        if (isSerializable(obj)) {
            // We need to get the JSON form of the object before recursing.
            // See https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/JSON/stringify#toJSON()_behavior
            obj = obj.toJSON();
        }
        if (Array.isArray(obj)) {
            let is_cloned = false;
            for (let i = 0; i < obj.length; i++) {
                const value = obj[i];
                if (value) {
                    if (value instanceof ArrayBuffer || ArrayBuffer.isView(value)) {
                        if (!is_cloned) {
                            obj = obj.slice();
                            is_cloned = true;
                        }
                        buffers.push(ArrayBuffer.isView(value) ? value.buffer : value);
                        buffer_paths.push(path.concat([i]));
                        // easier to just keep the array, but clear the entry, otherwise we have to think
                        // about array length, much easier this way
                        obj[i] = null;
                    }
                    else {
                        const new_value = remove(value, path.concat([i]));
                        // only assigned when the value changes, we may serialize objects that don't support assignment
                        if (new_value !== value) {
                            if (!is_cloned) {
                                obj = obj.slice();
                                is_cloned = true;
                            }
                            obj[i] = new_value;
                        }
                    }
                }
            }
        }
        else if (isObject(obj)) {
            for (const key in obj) {
                let is_cloned = false;
                if (Object.prototype.hasOwnProperty.call(obj, key)) {
                    const value = obj[key];
                    if (value) {
                        if (value instanceof ArrayBuffer || ArrayBuffer.isView(value)) {
                            if (!is_cloned) {
                                obj = Object.assign({}, obj);
                                is_cloned = true;
                            }
                            buffers.push(ArrayBuffer.isView(value) ? value.buffer : value);
                            buffer_paths.push(path.concat([key]));
                            delete obj[key]; // for objects/dicts we just delete them
                        }
                        else {
                            const new_value = remove(value, path.concat([key]));
                            // only assigned when the value changes, we may serialize objects that don't support assignment
                            if (new_value !== value) {
                                if (!is_cloned) {
                                    obj = Object.assign({}, obj);
                                    is_cloned = true;
                                }
                                obj[key] = new_value;
                            }
                        }
                    }
                }
            }
        }
        return obj;
    }
    const new_state = remove(state, []);
    return { state: new_state, buffers: buffers, buffer_paths: buffer_paths };
}
const BROKEN_FILE_SVG_ICON = `<svg style="height:50%;max-height: 50px;" role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 48 48">
<g >
  <g transform="translate(0.24520123,0.93464292)">
    <path  d="M 8.2494641,21.074514 V 5.6225142 c 0,-0.314 0.254,-0.567 0.57,-0.567 H 29.978464 c 2.388,0 9.268,5.8269998 9.268,8.3029998 v 5.5835 l -3.585749,4.407396 -2.772971,-3.535534 -5.126524,3.414213 -5.944543,-3.237436 -5.722718,3.06066 z m 30.9969999,3.8675 v 15.5835 c 0,0.314 -0.254,0.567 -0.57,0.567 H 8.8194641 c -0.315,0.002 -0.57,-0.251 -0.57,-0.566 v -15.452 l 7.8444949,2.628449 5.656854,-2.65165 4.24264,3.005204 5.833631,-3.237437 3.712311,3.944543 z" style="fill:url(#linearGradient3448);stroke:#888a85"  />
    <path d="m 30.383464,12.110514 c 4.108,0.159 7.304,-0.978 8.867,1.446 0.304,-3.9679998 -7.254,-8.8279998 -9.285,-8.4979998 0.813,0.498 0.418,7.0519998 0.418,7.0519998 z" style="fill:url(#linearGradient3445);stroke:#868a84" />
    <path enable-background="new" d="m 31.443464,11.086514 c 2.754,-0.019 4.106,-0.49 5.702,0.19 -1.299,-1.8809998 -4.358,-3.3439998 -5.728,-4.0279998 0.188,0.775 0.026,3.8379998 0.026,3.8379998 z" style="opacity:0.36930003;fill:none;stroke:url(#linearGradient3442)" />
  </g>
</g>
</svg>`;


/***/ }),

/***/ "../../packages/base/lib/version.js":
/*!******************************************!*\
  !*** ../../packages/base/lib/version.js ***!
  \******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   JUPYTER_WIDGETS_VERSION: () => (/* binding */ JUPYTER_WIDGETS_VERSION),
/* harmony export */   PROTOCOL_VERSION: () => (/* binding */ PROTOCOL_VERSION)
/* harmony export */ });
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
const JUPYTER_WIDGETS_VERSION = '2.0.0';
const PROTOCOL_VERSION = '2.1.0';


/***/ }),

/***/ "../../packages/base/lib/viewlist.js":
/*!*******************************************!*\
  !*** ../../packages/base/lib/viewlist.js ***!
  \*******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ViewList: () => (/* binding */ ViewList)
/* harmony export */ });
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
/**
 * - create_view and remove_view are default functions called when adding or removing views
 * - create_view takes a model and an index and returns a view or a promise for a view for that model
 * - remove_view takes a view and destroys it (including calling `view.remove()`)
 * - each time the update() function is called with a new list, the create and remove
 *   callbacks will be called in an order so that if you append the views created in the
 *   create callback and remove the views in the remove callback, you will duplicate
 *   the order of the list.
 * - the remove callback defaults to just removing the view (e.g., pass in null for the second parameter)
 * - the context defaults to the created ViewList.  If you pass another context, the create and remove
 *   will be called in that context.
 */
class ViewList {
    constructor(create_view, remove_view, context) {
        this.initialize(create_view, remove_view, context);
    }
    initialize(create_view, remove_view, context) {
        this._handler_context = context || this;
        this._models = [];
        this.views = []; // list of promises for views
        this._create_view = create_view;
        this._remove_view =
            remove_view ||
                function (view) {
                    view.remove();
                };
    }
    /**
     * the create_view, remove_view, and context arguments override the defaults
     * specified when the list is created.
     * after this function, the .views attribute is a list of promises for views
     * if you want to perform some action on the list of views, do something like
     * `Promise.all(myviewlist.views).then(function(views) {...});`
     */
    update(new_models, create_view, remove_view, context) {
        const remove = remove_view || this._remove_view;
        const create = create_view || this._create_view;
        context = context || this._handler_context;
        let i = 0;
        // first, skip past the beginning of the lists if they are identical
        for (; i < new_models.length; i++) {
            if (i >= this._models.length || new_models[i] !== this._models[i]) {
                break;
            }
        }
        const first_removed = i;
        // Remove the non-matching items from the old list.
        const removed = this.views.splice(first_removed, this.views.length - first_removed);
        for (let j = 0; j < removed.length; j++) {
            removed[j].then(function (view) {
                remove.call(context, view);
            });
        }
        // Add the rest of the new list items.
        for (; i < new_models.length; i++) {
            this.views.push(Promise.resolve(create.call(context, new_models[i], i)));
        }
        // make a copy of the input array
        this._models = new_models.slice();
        // return a promise that resolves to all of the resolved views
        return Promise.all(this.views);
    }
    /**
     * removes every view in the list; convenience function for `.update([])`
     * that should be faster
     * returns a promise that resolves after this removal is done
     */
    remove() {
        return Promise.all(this.views).then((views) => {
            views.forEach((value) => this._remove_view.call(this._handler_context, value));
            this.views = [];
            this._models = [];
        });
    }
    /**
     * Dispose this viewlist.
     *
     * A synchronous function which just deletes references to child views. This
     * function does not call .remove() on child views because that is
     * asynchronous. Use this in cases where child views will be removed in
     * another way.
     */
    dispose() {
        this.views = null;
        this._models = null;
    }
}


/***/ }),

/***/ "../../packages/base/lib/widget.js":
/*!*****************************************!*\
  !*** ../../packages/base/lib/widget.js ***!
  \*****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DOMWidgetModel: () => (/* binding */ DOMWidgetModel),
/* harmony export */   DOMWidgetView: () => (/* binding */ DOMWidgetView),
/* harmony export */   JupyterLuminoPanelWidget: () => (/* binding */ JupyterLuminoPanelWidget),
/* harmony export */   JupyterLuminoWidget: () => (/* binding */ JupyterLuminoWidget),
/* harmony export */   JupyterPhosphorPanelWidget: () => (/* binding */ JupyterPhosphorPanelWidget),
/* harmony export */   JupyterPhosphorWidget: () => (/* binding */ JupyterPhosphorWidget),
/* harmony export */   WidgetModel: () => (/* binding */ WidgetModel),
/* harmony export */   WidgetView: () => (/* binding */ WidgetView),
/* harmony export */   pack_models: () => (/* binding */ pack_models),
/* harmony export */   unpack_models: () => (/* binding */ unpack_models)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");
/* harmony import */ var _backbone_patch__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./backbone-patch */ "../../packages/base/lib/backbone-patch.js");
/* harmony import */ var backbone__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! backbone */ "../../node_modules/backbone/backbone.js");
/* harmony import */ var backbone__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(backbone__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! jquery */ "webpack/sharing/consume/default/jquery/jquery?123b");
/* harmony import */ var jquery__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(jquery__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _nativeview__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./nativeview */ "../../packages/base/lib/nativeview.js");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_5__);
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @lumino/messaging */ "webpack/sharing/consume/default/@lumino/messaging");
/* harmony import */ var _lumino_messaging__WEBPACK_IMPORTED_MODULE_6___default = /*#__PURE__*/__webpack_require__.n(_lumino_messaging__WEBPACK_IMPORTED_MODULE_6__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! @lumino/widgets */ "webpack/sharing/consume/default/@lumino/widgets");
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_7___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_7__);
/* harmony import */ var _version__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./version */ "../../packages/base/lib/version.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.









/**
 * The magic key used in the widget graph serialization.
 */
const IPY_MODEL_ = 'IPY_MODEL_';
/**
 * Replace model ids with models recursively.
 */
function unpack_models(value, manager // actually required, but typed to be compatible with ISerializers
) {
    if (Array.isArray(value)) {
        const unpacked = [];
        for (const sub_value of value) {
            unpacked.push(unpack_models(sub_value, manager));
        }
        return Promise.all(unpacked);
    }
    else if (value instanceof Object && typeof value !== 'string') {
        const unpacked = {};
        for (const [key, sub_value] of Object.entries(value)) {
            unpacked[key] = unpack_models(sub_value, manager);
        }
        return _utils__WEBPACK_IMPORTED_MODULE_0__.resolvePromisesDict(unpacked);
    }
    else if (typeof value === 'string' && value.slice(0, 10) === IPY_MODEL_) {
        // get_model returns a promise already
        return manager.get_model(value.slice(10, value.length));
    }
    else {
        return Promise.resolve(value);
    }
}
/** Replace models with ids recursively.
 *
 * If the commonly-used `unpack_models` is given as the `deseralize` method,
 * pack_models would be the appropriate `serialize`.
 * However, the default serialize method will have the same effect, when
 * `unpack_models` is used as the deserialize method.
 * This is to ensure backwards compatibility, see:
 *   https://github.com/jupyter-widgets/ipywidgets/pull/3738/commits/f9e27328bb631eb5247a7a6563595d3e655492c7#diff-efb19099381ae8911dd7f69b015a0138d08da7164512c1ee112aa75100bc9be2
 */
function pack_models(value, widget) {
    if (Array.isArray(value)) {
        const model_ids = [];
        for (const model of value) {
            model_ids.push(pack_models(model, widget));
        }
        return model_ids;
    }
    else if (value instanceof WidgetModel) {
        return `${IPY_MODEL_}${value.model_id}`;
    }
    else if (value instanceof Object && typeof value !== 'string') {
        const packed = {};
        for (const [key, sub_value] of Object.entries(value)) {
            packed[key] = pack_models(sub_value, widget);
        }
        return packed;
    }
    else {
        return value;
    }
}
class WidgetModel extends backbone__WEBPACK_IMPORTED_MODULE_2__.Model {
    /**
     * The default attributes.
     */
    defaults() {
        return {
            _model_module: '@jupyter-widgets/base',
            _model_name: 'WidgetModel',
            _model_module_version: _version__WEBPACK_IMPORTED_MODULE_8__.JUPYTER_WIDGETS_VERSION,
            _view_module: '@jupyter-widgets/base',
            _view_name: null,
            _view_module_version: _version__WEBPACK_IMPORTED_MODULE_8__.JUPYTER_WIDGETS_VERSION,
            _view_count: null,
        };
    }
    /**
     * Test to see if the model has been synced with the server.
     *
     * #### Notes
     * As of backbone 1.1, backbone ignores `patch` if it thinks the
     * model has never been pushed.
     */
    isNew() {
        return false;
    }
    /**
     * Constructor
     *
     * Initializes a WidgetModel instance. Called by the Backbone constructor.
     *
     * Parameters
     * ----------
     * widget_manager : WidgetManager instance
     * model_id : string
     *      An ID unique to this model.
     * comm : Comm instance (optional)
     */
    initialize(attributes, options) {
        this._expectedEchoMsgIds = new Map();
        this._attrsToUpdate = new Set();
        super.initialize(attributes, options);
        // Attributes should be initialized here, since user initialization may depend on it
        this.widget_manager = options.widget_manager;
        this.model_id = options.model_id;
        const comm = options.comm;
        this.views = Object.create(null);
        this.state_change = Promise.resolve();
        this._closed = false;
        this._state_lock = null;
        this._msg_buffer = null;
        this._msg_buffer_callbacks = null;
        this._pending_msgs = 0;
        // _buffered_state_diff must be created *after* the super.initialize
        // call above. See the note in the set() method below.
        this._buffered_state_diff = {};
        if (comm) {
            // Remember comm associated with the model.
            this.comm = comm;
            // Hook comm messages up to model.
            comm.on_close(this._handle_comm_closed.bind(this));
            comm.on_msg(this._handle_comm_msg.bind(this));
            this.comm_live = true;
        }
        else {
            this.comm_live = false;
        }
    }
    get comm_live() {
        return this._comm_live;
    }
    set comm_live(x) {
        this._comm_live = x;
        this.trigger('comm_live_update');
    }
    /**
     * Send a custom msg over the comm.
     */
    send(content, callbacks, buffers) {
        if (this.comm !== undefined) {
            const data = { method: 'custom', content: content };
            this.comm.send(data, callbacks, {}, buffers);
        }
    }
    /**
     * Close model
     *
     * @param comm_closed - true if the comm is already being closed. If false, the comm will be closed.
     *
     * @returns - a promise that is fulfilled when all the associated views have been removed.
     */
    close(comm_closed = false) {
        // can only be closed once.
        if (this._closed) {
            return Promise.resolve();
        }
        this._closed = true;
        if (this.comm && !comm_closed) {
            this.comm.close();
        }
        this.stopListening();
        this.trigger('destroy', this);
        if (this.comm) {
            delete this.comm;
        }
        // Delete all views of this model
        if (this.views) {
            const views = Object.keys(this.views).map((id) => {
                return this.views[id].then((view) => view.remove());
            });
            delete this.views;
            return Promise.all(views).then(() => {
                return;
            });
        }
        return Promise.resolve();
    }
    /**
     * Handle when a widget comm is closed.
     */
    _handle_comm_closed(msg) {
        this.trigger('comm:close');
        this.close(true);
    }
    /**
     * Handle incoming comm msg.
     */
    _handle_comm_msg(msg) {
        const data = msg.content.data;
        const method = data.method;
        switch (method) {
            case 'update':
            case 'echo_update':
                this.state_change = this.state_change
                    .then(() => {
                    var _a, _b, _c;
                    const state = data.state;
                    const buffer_paths = (_a = data.buffer_paths) !== null && _a !== void 0 ? _a : [];
                    const buffers = (_c = (_b = msg.buffers) === null || _b === void 0 ? void 0 : _b.slice(0, buffer_paths.length)) !== null && _c !== void 0 ? _c : [];
                    _utils__WEBPACK_IMPORTED_MODULE_0__.put_buffers(state, buffer_paths, buffers);
                    if (msg.parent_header && method === 'echo_update') {
                        const msgId = msg.parent_header.msg_id;
                        // we may have echos coming from other clients, we only care about
                        // dropping echos for which we expected a reply
                        const expectedEcho = Object.keys(state).filter((attrName) => this._expectedEchoMsgIds.has(attrName));
                        expectedEcho.forEach((attrName) => {
                            // Skip echo messages until we get the reply we are expecting.
                            const isOldMessage = this._expectedEchoMsgIds.get(attrName) !== msgId;
                            if (isOldMessage) {
                                // Ignore an echo update that comes before our echo.
                                delete state[attrName];
                            }
                            else {
                                // we got our echo confirmation, so stop looking for it
                                this._expectedEchoMsgIds.delete(attrName);
                                // Start accepting echo updates unless we plan to send out a new state soon
                                if (this._msg_buffer !== null &&
                                    Object.prototype.hasOwnProperty.call(this._msg_buffer, attrName)) {
                                    delete state[attrName];
                                }
                            }
                        });
                    }
                    return this.constructor._deserialize_state(
                    // Combine the state updates, with preference for kernel updates
                    state, this.widget_manager);
                })
                    .then((state) => {
                    this.set_state(state);
                })
                    .catch(_utils__WEBPACK_IMPORTED_MODULE_0__.reject(`Could not process update msg for model id: ${this.model_id}`, true));
                return this.state_change;
            case 'custom':
                this.trigger('msg:custom', data.content, msg.buffers);
                return Promise.resolve();
        }
        return Promise.resolve();
    }
    /**
     * Handle when a widget is updated from the backend.
     *
     * This function is meant for internal use only. Values set here will not be propagated on a sync.
     */
    set_state(state) {
        this._state_lock = state;
        try {
            this.set(state);
        }
        catch (e) {
            console.error(`Error setting state: ${e instanceof Error ? e.message : e}`);
        }
        finally {
            this._state_lock = null;
        }
    }
    /**
     * Get the serializable state of the model.
     *
     * If drop_default is truthy, attributes that are equal to their default
     * values are dropped.
     */
    get_state(drop_defaults) {
        const fullState = this.attributes;
        if (drop_defaults) {
            // if defaults is a function, call it
            const d = this.defaults;
            const defaults = typeof d === 'function' ? d.call(this) : d;
            const state = {};
            Object.keys(fullState).forEach((key) => {
                if (!_utils__WEBPACK_IMPORTED_MODULE_0__.isEqual(fullState[key], defaults[key])) {
                    state[key] = fullState[key];
                }
            });
            return state;
        }
        else {
            return Object.assign({}, fullState);
        }
    }
    /**
     * Handle status msgs.
     *
     * execution_state : ('busy', 'idle', 'starting')
     */
    _handle_status(msg) {
        if (this.comm !== void 0) {
            if (msg.content.execution_state === 'idle') {
                this._pending_msgs--;
                // Sanity check for logic errors that may push this below zero.
                if (this._pending_msgs < 0) {
                    console.error(`Jupyter Widgets message throttle: Pending messages < 0 (=${this._pending_msgs}), which is unexpected. Resetting to 0 to continue.`);
                    this._pending_msgs = 0; // do not break message throttling in case of unexpected errors
                }
                // Send buffer if one is waiting and we are below the throttle.
                if (this._msg_buffer !== null && this._pending_msgs < 1) {
                    const msgId = this.send_sync_message(this._msg_buffer, this._msg_buffer_callbacks);
                    this.rememberLastUpdateFor(msgId);
                    this._msg_buffer = null;
                    this._msg_buffer_callbacks = null;
                }
            }
        }
    }
    /**
     * Create msg callbacks for a comm msg.
     */
    callbacks(view) {
        return this.widget_manager.callbacks(view);
    }
    /**
     * Set one or more values.
     *
     * We just call the super method, in which val and options are optional.
     * Handles both "key", value and {key: value} -style arguments.
     */
    set(key, val, options) {
        // Call our patched backbone set. See #1642 and #1643.
        const return_value = _backbone_patch__WEBPACK_IMPORTED_MODULE_1__.set.call(this, key, val, options);
        // Backbone only remembers the diff of the most recent set()
        // operation.  Calling set multiple times in a row results in a
        // loss of change information.  Here we keep our own running diff.
        //
        // We don't buffer the state set in the constructor (including
        // defaults), so we first check to see if we've initialized _buffered_state_diff.
        // which happens after the constructor sets attributes at creation.
        if (this._buffered_state_diff !== void 0) {
            const attrs = this.changedAttributes() || {};
            // The state_lock lists attributes that are currently being changed
            // right now from a kernel message. We don't want to send these
            // non-changes back to the kernel, so we delete them out of attrs if
            // they haven't changed from their state_lock value.
            // The state lock could be null or undefined (if set is being called from
            // the initializer).
            if (this._state_lock) {
                for (const key of Object.keys(this._state_lock)) {
                    if (attrs[key] === this._state_lock[key]) {
                        delete attrs[key];
                    }
                }
            }
            // _buffered_state_diff_synced lists things that have already been sent to the kernel during a top-level call to .set(), so we don't need to buffer these things either.
            if (this._buffered_state_diff_synced) {
                for (const key of Object.keys(this._buffered_state_diff_synced)) {
                    if (attrs[key] === this._buffered_state_diff_synced[key]) {
                        delete attrs[key];
                    }
                }
            }
            this._buffered_state_diff = _utils__WEBPACK_IMPORTED_MODULE_0__.assign(this._buffered_state_diff, attrs);
        }
        // If this ended a top-level call to .set, then reset _buffered_state_diff_synced
        if (this._changing === false) {
            this._buffered_state_diff_synced = {};
        }
        return return_value;
    }
    /**
     * Handle sync to the back-end.  Called when a model.save() is called.
     *
     * Make sure a comm exists.
     *
     * Parameters
     * ----------
     * method : create, update, patch, delete, read
     *   create/update always send the full attribute set
     *   patch - only send attributes listed in options.attrs, and if we
     *   are queuing up messages, combine with previous messages that have
     *   not been sent yet
     * model : the model we are syncing
     *   will normally be the same as `this`
     * options : dict
     *   the `attrs` key, if it exists, gives an {attr: value} dict that
     *   should be synced, otherwise, sync all attributes.
     *
     */
    sync(method, model, options = {}) {
        // the typing is to return `any` since the super.sync method returns a JqXHR, but we just return false if there is an error.
        if (this.comm === undefined) {
            throw 'Syncing error: no comm channel defined';
        }
        const attrs = method === 'patch'
            ? options.attrs
            : model.get_state(options.drop_defaults);
        // The state_lock lists attributes that are currently being changed
        // right now from a kernel message. We don't want to send these
        // non-changes back to the kernel, so we delete them out of attrs if
        // they haven't changed from their state_lock value.
        // The state lock could be null or undefined (if this is triggered
        // from the initializer).
        if (this._state_lock) {
            for (const key of Object.keys(this._state_lock)) {
                if (attrs[key] === this._state_lock[key]) {
                    delete attrs[key];
                }
            }
        }
        Object.keys(attrs).forEach((attrName) => {
            this._attrsToUpdate.add(attrName);
        });
        const msgState = this.serialize(attrs);
        if (Object.keys(msgState).length > 0) {
            // If this message was sent via backbone itself, it will not
            // have any callbacks.  It's important that we create callbacks
            // so we can listen for status messages, etc...
            const callbacks = options.callbacks || this.callbacks();
            // Check throttle.
            if (this._pending_msgs >= 1) {
                // The throttle has been exceeded, buffer the current msg so
                // it can be sent once the kernel has finished processing
                // some of the existing messages.
                // Combine updates if it is a 'patch' sync, otherwise replace updates
                switch (method) {
                    case 'patch':
                        this._msg_buffer = _utils__WEBPACK_IMPORTED_MODULE_0__.assign(this._msg_buffer || {}, msgState);
                        break;
                    case 'update':
                    case 'create':
                        this._msg_buffer = msgState;
                        break;
                    default:
                        throw 'unrecognized syncing method';
                }
                this._msg_buffer_callbacks = callbacks;
            }
            else {
                // We haven't exceeded the throttle, send the message like
                // normal.
                const msgId = this.send_sync_message(attrs, callbacks);
                this.rememberLastUpdateFor(msgId);
                // Since the comm is a one-way communication, assume the message
                // arrived and was processed successfully.
                // Don't call options.success since we don't have a model back from
                // the server. Note that this means we don't have the Backbone
                // 'sync' event.
            }
        }
    }
    rememberLastUpdateFor(msgId) {
        this._attrsToUpdate.forEach((attrName) => {
            this._expectedEchoMsgIds.set(attrName, msgId);
        });
        this._attrsToUpdate = new Set();
    }
    /**
     * Serialize widget state.
     *
     * A serializer is a function which takes in a state attribute and a widget,
     * and synchronously returns a JSONable object. The returned object will
     * have toJSON called if possible, and the final result should be a
     * primitive object that is a snapshot of the widget state that may have
     * binary array buffers.
     */
    serialize(state) {
        const serializers = this.constructor.serializers ||
            _lumino_coreutils__WEBPACK_IMPORTED_MODULE_5__.JSONExt.emptyObject;
        for (const k of Object.keys(state)) {
            try {
                if (serializers[k] && serializers[k].serialize) {
                    state[k] = serializers[k].serialize(state[k], this);
                }
                else {
                    // the default serializer just deep-copies the object
                    state[k] = JSON.parse(JSON.stringify(state[k]));
                }
                if (state[k] && state[k].toJSON) {
                    state[k] = state[k].toJSON();
                }
            }
            catch (e) {
                console.error('Error serializing widget state attribute: ', k);
                throw e;
            }
        }
        return state;
    }
    /**
     * Send a sync message to the kernel.
     *
     * If a message is sent successfully, this returns the message ID of that
     * message. Otherwise it returns an empty string
     */
    send_sync_message(state, callbacks = {}) {
        if (!this.comm) {
            return '';
        }
        try {
            // Make a 2-deep copy so we don't modify the caller's callbacks object.
            callbacks = {
                shell: Object.assign({}, callbacks.shell),
                iopub: Object.assign({}, callbacks.iopub),
                input: callbacks.input,
            };
            // Save the caller's status callback so we can call it after we handle the message.
            const statuscb = callbacks.iopub.status;
            callbacks.iopub.status = (msg) => {
                this._handle_status(msg);
                if (statuscb) {
                    statuscb(msg);
                }
            };
            // split out the binary buffers
            const split = _utils__WEBPACK_IMPORTED_MODULE_0__.remove_buffers(state);
            const msgId = this.comm.send({
                method: 'update',
                state: split.state,
                buffer_paths: split.buffer_paths,
            }, callbacks, {}, split.buffers);
            this._pending_msgs++;
            return msgId;
        }
        catch (e) {
            console.error('Could not send widget sync message', e);
        }
        return '';
    }
    /**
     * Push this model's state to the back-end
     *
     * This invokes a Backbone.Sync.
     */
    save_changes(callbacks) {
        if (this.comm_live) {
            const options = { patch: true };
            if (callbacks) {
                options.callbacks = callbacks;
            }
            this.save(this._buffered_state_diff, options);
            // If we are currently in a .set() call, save what state we have synced
            // to the kernel so we don't buffer it again as we come out of the .set call.
            if (this._changing) {
                _utils__WEBPACK_IMPORTED_MODULE_0__.assign(this._buffered_state_diff_synced, this._buffered_state_diff);
            }
            this._buffered_state_diff = {};
        }
    }
    /**
     * on_some_change(['key1', 'key2'], foo, context) differs from
     * on('change:key1 change:key2', foo, context).
     * If the widget attributes key1 and key2 are both modified,
     * the second form will result in foo being called twice
     * while the first will call foo only once.
     */
    on_some_change(keys, callback, context) {
        this.on('change', (...args) => {
            if (keys.some(this.hasChanged, this)) {
                callback.apply(context, args);
            }
        }, this);
    }
    /**
     * Serialize the model.  See the deserialization function at the top of this file
     * and the kernel-side serializer/deserializer.
     */
    toJSON(options) {
        return `IPY_MODEL_${this.model_id}`;
    }
    /**
     * Returns a promise for the deserialized state. The second argument
     * is an instance of widget manager, which is required for the
     * deserialization of widget models.
     */
    static _deserialize_state(state, manager) {
        const serializers = this.serializers;
        let deserialized;
        if (serializers) {
            deserialized = {};
            for (const k in state) {
                if (serializers[k] && serializers[k].deserialize) {
                    deserialized[k] = serializers[k].deserialize(state[k], manager);
                }
                else {
                    deserialized[k] = state[k];
                }
            }
        }
        else {
            deserialized = state;
        }
        return _utils__WEBPACK_IMPORTED_MODULE_0__.resolvePromisesDict(deserialized);
    }
}
class DOMWidgetModel extends WidgetModel {
    defaults() {
        return _utils__WEBPACK_IMPORTED_MODULE_0__.assign(super.defaults(), {
            _dom_classes: [],
            tabbable: null,
            tooltip: null,
            // We do not declare defaults for the layout and style attributes.
            // Those defaults are constructed on the kernel side and synced here
            // as needed, and our code here copes with those attributes being
            // undefined. See
            // https://github.com/jupyter-widgets/ipywidgets/issues/1620 and
            // https://github.com/jupyter-widgets/ipywidgets/pull/1621
        });
    }
}
DOMWidgetModel.serializers = Object.assign(Object.assign({}, WidgetModel.serializers), { layout: { deserialize: unpack_models }, style: { deserialize: unpack_models } });
class WidgetView extends _nativeview__WEBPACK_IMPORTED_MODULE_4__.NativeView {
    /**
     * Public constructor.
     */
    constructor(options) {
        super(options);
    }
    /**
     * Initializer, called at the end of the constructor.
     */
    initialize(parameters) {
        this.listenTo(this.model, 'change', (model, options) => {
            const changed = Object.keys(this.model.changedAttributes() || {});
            if (changed[0] === '_view_count' && changed.length === 1) {
                // Just the view count was updated
                return;
            }
            this.update(options);
        });
        this.options = parameters.options;
        this.once('remove', () => {
            if (typeof this.model.get('_view_count') === 'number') {
                this.model.set('_view_count', this.model.get('_view_count') - 1);
                this.model.save_changes();
            }
        });
        this.once('displayed', () => {
            if (typeof this.model.get('_view_count') === 'number') {
                this.model.set('_view_count', this.model.get('_view_count') + 1);
                this.model.save_changes();
            }
        });
        this.displayed = new Promise((resolve, reject) => {
            this.once('displayed', resolve);
            this.model.on('msg:custom', this.handle_message.bind(this));
        });
    }
    /**
     * Handle message sent to the front end.
     *
     * Used to focus or blur the widget.
     */
    handle_message(content) {
        if (content.do === 'focus') {
            this.el.focus();
        }
        else if (content.do === 'blur') {
            this.el.blur();
        }
    }
    /**
     * Triggered on model change.
     *
     * Update view to be consistent with this.model
     */
    update(options) {
        return;
    }
    /**
     * Render a view
     *
     * @returns the view or a promise to the view.
     */
    render() {
        return;
    }
    create_child_view(child_model, options = {}) {
        options = Object.assign({ parent: this }, options);
        return this.model.widget_manager
            .create_view(child_model, options)
            .catch(_utils__WEBPACK_IMPORTED_MODULE_0__.reject('Could not create child view', true));
    }
    /**
     * Create msg callbacks for a comm msg.
     */
    callbacks() {
        return this.model.callbacks(this);
    }
    /**
     * Send a custom msg associated with this view.
     */
    send(content, buffers) {
        this.model.send(content, this.callbacks(), buffers);
    }
    touch() {
        this.model.save_changes(this.callbacks());
    }
    remove() {
        // Raise a remove event when the view is removed.
        super.remove();
        this.trigger('remove');
        return this;
    }
}
class JupyterLuminoWidget extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_7__.Widget {
    constructor(options) {
        const view = options.view;
        // Cast as any since we cannot delete a mandatory value
        delete options.view;
        super(options);
        this._view = view;
    }
    /**
     * Dispose the widget.
     *
     * This causes the view to be destroyed as well with 'remove'
     */
    dispose() {
        if (this.isDisposed) {
            return;
        }
        super.dispose();
        this._view.remove();
        this._view = null;
    }
    /**
     * Process the Lumino message.
     *
     * Any custom Lumino widget used inside a Jupyter widget should override
     * the processMessage function like this.
     */
    processMessage(msg) {
        super.processMessage(msg);
        this._view.processLuminoMessage(msg);
    }
}
/**
 * @deprecated Use {@link JupyterLuminoWidget} instead (Since 8.0).
 */
const JupyterPhosphorWidget = JupyterLuminoWidget;
class JupyterLuminoPanelWidget extends _lumino_widgets__WEBPACK_IMPORTED_MODULE_7__.Panel {
    constructor(options) {
        const view = options.view;
        delete options.view;
        super(options);
        this._view = view;
    }
    /**
     * Process the Lumino message.
     *
     * Any custom Lumino widget used inside a Jupyter widget should override
     * the processMessage function like this.
     */
    processMessage(msg) {
        super.processMessage(msg);
        this._view.processLuminoMessage(msg);
    }
    /**
     * Dispose the widget.
     *
     * This causes the view to be destroyed as well with 'remove'
     */
    dispose() {
        var _a;
        if (this.isDisposed) {
            return;
        }
        super.dispose();
        (_a = this._view) === null || _a === void 0 ? void 0 : _a.remove();
        this._view = null;
    }
}
/**
 * @deprecated Use {@link JupyterLuminoPanelWidget} instead (Since 8.0).
 */
const JupyterPhosphorPanelWidget = JupyterLuminoPanelWidget;
class DOMWidgetView extends WidgetView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        super.initialize(parameters);
        this.listenTo(this.model, 'change:_dom_classes', (model, new_classes) => {
            const old_classes = model.previous('_dom_classes');
            this.update_classes(old_classes, new_classes);
        });
        this.layoutPromise = Promise.resolve();
        this.listenTo(this.model, 'change:layout', (model, value) => {
            this.setLayout(value, model.previous('layout'));
        });
        this.stylePromise = Promise.resolve();
        this.listenTo(this.model, 'change:style', (model, value) => {
            this.setStyle(value, model.previous('style'));
        });
        this.displayed.then(() => {
            this.update_classes([], this.model.get('_dom_classes'));
            this.setLayout(this.model.get('layout'));
            this.setStyle(this.model.get('style'));
        });
        this._comm_live_update();
        this.listenTo(this.model, 'comm_live_update', () => {
            this._comm_live_update();
        });
        this.listenTo(this.model, 'change:tooltip', this.updateTooltip);
        this.updateTooltip();
    }
    setLayout(layout, oldLayout) {
        if (layout) {
            this.layoutPromise = this.layoutPromise.then((oldLayoutView) => {
                if (oldLayoutView) {
                    oldLayoutView.unlayout();
                    this.stopListening(oldLayoutView.model);
                    oldLayoutView.remove();
                }
                return this.create_child_view(layout)
                    .then((view) => {
                    // Trigger the displayed event of the child view.
                    return this.displayed.then(() => {
                        view.trigger('displayed');
                        this.listenTo(view.model, 'change', () => {
                            // Post (asynchronous) so layout changes can take
                            // effect first.
                            _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__.MessageLoop.postMessage(this.luminoWidget, _lumino_widgets__WEBPACK_IMPORTED_MODULE_7__.Widget.ResizeMessage.UnknownSize);
                        });
                        _lumino_messaging__WEBPACK_IMPORTED_MODULE_6__.MessageLoop.postMessage(this.luminoWidget, _lumino_widgets__WEBPACK_IMPORTED_MODULE_7__.Widget.ResizeMessage.UnknownSize);
                        this.trigger('layout-changed');
                        return view;
                    });
                })
                    .catch(_utils__WEBPACK_IMPORTED_MODULE_0__.reject('Could not add LayoutView to DOMWidgetView', true));
            });
        }
    }
    setStyle(style, oldStyle) {
        if (style) {
            this.stylePromise = this.stylePromise.then((oldStyleView) => {
                if (oldStyleView) {
                    oldStyleView.unstyle();
                    this.stopListening(oldStyleView.model);
                    oldStyleView.remove();
                }
                return this.create_child_view(style)
                    .then((view) => {
                    // Trigger the displayed event of the child view.
                    return this.displayed.then(() => {
                        view.trigger('displayed');
                        this.trigger('style-changed');
                        // Unlike for the layout attribute, style changes don't
                        // trigger Lumino resize messages.
                        return view;
                    });
                })
                    .catch(_utils__WEBPACK_IMPORTED_MODULE_0__.reject('Could not add styleView to DOMWidgetView', true));
            });
        }
    }
    updateTooltip() {
        const title = this.model.get('tooltip');
        if (!title) {
            this.el.removeAttribute('title');
        }
        else if (this.model.get('description').length === 0) {
            this.el.setAttribute('title', title);
        }
    }
    /**
     * Update the DOM classes applied to an element, default to this.el.
     */
    update_classes(old_classes, new_classes, el) {
        if (el === undefined) {
            el = this.el;
        }
        _utils__WEBPACK_IMPORTED_MODULE_0__.difference(old_classes, new_classes).map(function (c) {
            if (el.classList) {
                // classList is not supported by IE for svg elements
                el.classList.remove(c);
            }
            else {
                el.setAttribute('class', el.getAttribute('class').replace(c, ''));
            }
        });
        _utils__WEBPACK_IMPORTED_MODULE_0__.difference(new_classes, old_classes).map(function (c) {
            if (el.classList) {
                // classList is not supported by IE for svg elements
                el.classList.add(c);
            }
            else {
                el.setAttribute('class', el.getAttribute('class').concat(' ', c));
            }
        });
    }
    /**
     * Update the DOM classes applied to the widget based on a single
     * trait's value.
     *
     * Given a trait value classes map, this function automatically
     * handles applying the appropriate classes to the widget element
     * and removing classes that are no longer valid.
     *
     * Parameters
     * ----------
     * class_map: dictionary
     *  Dictionary of trait values to class lists.
     *  Example:
     *      {
     *          success: ['alert', 'alert-success'],
     *          info: ['alert', 'alert-info'],
     *          warning: ['alert', 'alert-warning'],
     *          danger: ['alert', 'alert-danger']
     *      };
     * trait_name: string
     *  Name of the trait to check the value of.
     * el: optional DOM element handle, defaults to this.el
     *  Element that the classes are applied to.
     */
    update_mapped_classes(class_map, trait_name, el) {
        let key = this.model.previous(trait_name);
        const old_classes = class_map[key] ? class_map[key] : [];
        key = this.model.get(trait_name);
        const new_classes = class_map[key] ? class_map[key] : [];
        this.update_classes(old_classes, new_classes, el || this.el);
    }
    set_mapped_classes(class_map, trait_name, el) {
        const key = this.model.get(trait_name);
        const new_classes = class_map[key] ? class_map[key] : [];
        this.update_classes([], new_classes, el || this.el);
    }
    _setElement(el) {
        if (this.luminoWidget) {
            this.luminoWidget.dispose();
        }
        this.$el = el instanceof (jquery__WEBPACK_IMPORTED_MODULE_3___default()) ? el : jquery__WEBPACK_IMPORTED_MODULE_3___default()(el);
        this.el = this.$el[0];
        this.luminoWidget = new JupyterLuminoWidget({
            node: el,
            view: this,
        });
    }
    remove() {
        if (this.luminoWidget) {
            this.luminoWidget.dispose();
        }
        return super.remove();
    }
    /**
     * @deprecated Use {@link processLuminoMessage} instead (Since 8.0).
     */
    processPhosphorMessage(msg) {
        this.processLuminoMessage(msg);
    }
    processLuminoMessage(msg) {
        switch (msg.type) {
            case 'after-attach':
                this.trigger('displayed');
                break;
            case 'show':
                this.trigger('shown');
                break;
        }
    }
    _comm_live_update() {
        if (this.model.comm_live) {
            this.luminoWidget.removeClass('jupyter-widgets-disconnected');
        }
        else {
            this.luminoWidget.addClass('jupyter-widgets-disconnected');
        }
    }
    updateTabindex() {
        const tabbable = this.model.get('tabbable');
        if (tabbable === true) {
            this.el.setAttribute('tabIndex', '0');
        }
        else if (tabbable === false) {
            this.el.setAttribute('tabIndex', '-1');
        }
        else if (tabbable === null) {
            this.el.removeAttribute('tabIndex');
        }
    }
    /**
     * @deprecated Use {@link luminoWidget} instead (Since 8.0).
     */
    get pWidget() {
        return this.luminoWidget;
    }
    /**
     * @deprecated Use {@link luminoWidget} instead (Since 8.0).
     */
    set pWidget(value) {
        this.luminoWidget = value;
    }
}


/***/ }),

/***/ "../../packages/base/lib/widget_layout.js":
/*!************************************************!*\
  !*** ../../packages/base/lib/widget_layout.js ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   LayoutModel: () => (/* binding */ LayoutModel),
/* harmony export */   LayoutView: () => (/* binding */ LayoutView)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");
/* harmony import */ var _widget__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget */ "../../packages/base/lib/widget.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


/**
 * css properties exposed by the layout widget with their default values.
 */
const css_properties = {
    align_content: null,
    align_items: null,
    align_self: null,
    border_top: null,
    border_right: null,
    border_bottom: null,
    border_left: null,
    bottom: null,
    display: null,
    flex: null,
    flex_flow: null,
    height: null,
    justify_content: null,
    justify_items: null,
    left: null,
    margin: null,
    max_height: null,
    max_width: null,
    min_height: null,
    min_width: null,
    overflow: null,
    order: null,
    padding: null,
    right: null,
    top: null,
    visibility: null,
    width: null,
    // image-specific
    object_fit: null,
    object_position: null,
    // container
    grid_auto_columns: null,
    grid_auto_flow: null,
    grid_auto_rows: null,
    grid_gap: null,
    grid_template_rows: null,
    grid_template_columns: null,
    grid_template_areas: null,
    // items
    grid_row: null,
    grid_column: null,
    grid_area: null,
};
class LayoutModel extends _widget__WEBPACK_IMPORTED_MODULE_1__.WidgetModel {
    defaults() {
        return (0,_utils__WEBPACK_IMPORTED_MODULE_0__.assign)(super.defaults(), {
            _model_name: 'LayoutModel',
            _view_name: 'LayoutView',
        }, css_properties);
    }
}
class LayoutView extends _widget__WEBPACK_IMPORTED_MODULE_1__.WidgetView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        this._traitNames = [];
        super.initialize(parameters);
        // Register the traits that live on the Python side
        for (const key of Object.keys(css_properties)) {
            this.registerTrait(key);
        }
    }
    /**
     * Register a CSS trait that is known by the model
     * @param trait
     */
    registerTrait(trait) {
        this._traitNames.push(trait);
        // Listen to changes, and set the value on change.
        this.listenTo(this.model, 'change:' + trait, (model, value) => {
            this.handleChange(trait, value);
        });
        // Set the initial value on display.
        this.handleChange(trait, this.model.get(trait));
    }
    /**
     * Get the the name of the css property from the trait name
     * @param  model attribute name
     * @return css property name
     */
    css_name(trait) {
        return trait.replace(/_/g, '-');
    }
    /**
     * Handles when a trait value changes
     */
    handleChange(trait, value) {
        // should be synchronous so that we can measure later.
        const parent = this.options.parent;
        if (parent) {
            if (value === null) {
                parent.el.style.removeProperty(this.css_name(trait));
            }
            else {
                parent.el.style.setProperty(this.css_name(trait), value);
            }
        }
        else {
            console.warn('Style not applied because a parent view does not exist');
        }
    }
    /**
     * Remove the styling from the parent view.
     */
    unlayout() {
        const parent = this.options.parent;
        this._traitNames.forEach((trait) => {
            if (parent) {
                parent.el.style.removeProperty(this.css_name(trait));
            }
            else {
                console.warn('Style not removed because a parent view does not exist');
            }
        }, this);
    }
}


/***/ }),

/***/ "../../packages/base/lib/widget_style.js":
/*!***********************************************!*\
  !*** ../../packages/base/lib/widget_style.js ***!
  \***********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   StyleModel: () => (/* binding */ StyleModel),
/* harmony export */   StyleView: () => (/* binding */ StyleView)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "../../packages/base/lib/utils.js");
/* harmony import */ var _widget__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./widget */ "../../packages/base/lib/widget.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


/**
 * Three functions to deal with some CSS attributes
 * to make them easier to use.
 */
class StyleModel extends _widget__WEBPACK_IMPORTED_MODULE_1__.WidgetModel {
    defaults() {
        const Derived = this.constructor;
        return (0,_utils__WEBPACK_IMPORTED_MODULE_0__.assign)(super.defaults(), {
            _model_name: 'StyleModel',
            _view_name: 'StyleView',
        }, Object.keys(Derived.styleProperties).reduce((obj, key) => {
            obj[key] = Derived.styleProperties[key].default;
            return obj;
        }, {}));
    }
}
StyleModel.styleProperties = {};
class StyleView extends _widget__WEBPACK_IMPORTED_MODULE_1__.WidgetView {
    /**
     * Public constructor
     */
    initialize(parameters) {
        this._traitNames = [];
        super.initialize(parameters);
        // Register the traits that live on the Python side
        const ModelType = this.model.constructor;
        for (const key of Object.keys(ModelType.styleProperties)) {
            this.registerTrait(key);
        }
        // Set the initial styles
        this.style();
    }
    /**
     * Register a CSS trait that is known by the model
     * @param trait
     */
    registerTrait(trait) {
        this._traitNames.push(trait);
        // Listen to changes, and set the value on change.
        this.listenTo(this.model, 'change:' + trait, (model, value) => {
            this.handleChange(trait, value);
        });
    }
    /**
     * Handles when a trait value changes
     */
    handleChange(trait, value) {
        // should be synchronous so that we can measure later.
        const parent = this.options.parent;
        if (parent) {
            const ModelType = this.model.constructor;
            const styleProperties = ModelType.styleProperties;
            const attribute = styleProperties[trait].attribute;
            const selector = styleProperties[trait].selector;
            const elements = selector
                ? parent.el.querySelectorAll(selector)
                : [parent.el];
            if (value === null) {
                for (let i = 0; i !== elements.length; ++i) {
                    elements[i].style.removeProperty(attribute);
                }
            }
            else {
                for (let i = 0; i !== elements.length; ++i) {
                    elements[i].style.setProperty(attribute, value);
                }
            }
        }
        else {
            console.warn('Style not applied because a parent view does not exist');
        }
    }
    /**
     * Apply styles for all registered traits
     */
    style() {
        for (const trait of this._traitNames) {
            this.handleChange(trait, this.model.get(trait));
        }
    }
    /**
     * Remove the styling from the parent view.
     */
    unstyle() {
        const parent = this.options.parent;
        const ModelType = this.model.constructor;
        const styleProperties = ModelType.styleProperties;
        this._traitNames.forEach((trait) => {
            if (parent) {
                const attribute = styleProperties[trait].attribute;
                const selector = styleProperties[trait].selector;
                const elements = selector
                    ? parent.el.querySelectorAll(selector)
                    : [parent.el];
                for (let i = 0; i !== elements.length; ++i) {
                    elements[i].style.removeProperty(attribute);
                }
            }
            else {
                console.warn('Style not removed because a parent view does not exist');
            }
        }, this);
    }
}


/***/ })

}]);
//# sourceMappingURL=packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery.5dd13f8e980fa3c50bfe.js.map