"use strict";
(self["webpackChunk_jupyterlab_galata_extension"] = self["webpackChunk_jupyterlab_galata_extension"] || []).push([["lib_extension_index_js"],{

/***/ "../lib/extension/global.js":
/*!**********************************!*\
  !*** ../lib/extension/global.js ***!
  \**********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   GalataInpage: () => (/* binding */ GalataInpage)
/* harmony export */ });
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/algorithm */ "webpack/sharing/consume/default/@lumino/algorithm");
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! @lumino/signaling */ "webpack/sharing/consume/default/@lumino/signaling");
/* harmony import */ var _lumino_signaling__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_lumino_signaling__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _tokens__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./tokens */ "../lib/extension/tokens.js");
/* eslint-disable @typescript-eslint/ban-types */
/*
 * Copyright (c) Jupyter Development Team.
 * Distributed under the terms of the Modified BSD License.
 * Copyright (c) Bloomberg Finance LP.
 */



const PLUGIN_ID_DOC_MANAGER = '@jupyterlab/docmanager-extension:manager';
const PLUGIN_ID_ROUTER = '@jupyterlab/application-extension:router';
/**
 * In-Page Galata helpers
 */
class GalataInpage {
    constructor() {
        this.listeners = new WeakMap();
    }
    /**
     * Get an application plugin
     *
     * @param pluginId Plugin ID
     * @returns Application plugin
     */
    async getPlugin(pluginId) {
        return new Promise((resolve, reject) => {
            const app = this._app;
            const hasPlugin = app.hasPlugin(pluginId);
            if (hasPlugin) {
                try {
                    const appAny = app;
                    const plugin = appAny.pluginRegistry._plugins
                        ? appAny.pluginRegistry._plugins.get(pluginId)
                        : undefined;
                    if (plugin.activated) {
                        resolve(plugin.service);
                    }
                    else {
                        void app.activatePlugin(pluginId).then(response => {
                            resolve(plugin.service);
                        });
                    }
                }
                catch (error) {
                    console.error('Failed to get plugin', error);
                }
            }
        });
    }
    /**
     * Get the Jupyter notifications
     *
     * @returns Jupyter Notifications
     */
    async getNotifications() {
        var _a;
        const plugin = await this.getPlugin(_tokens__WEBPACK_IMPORTED_MODULE_2__.PLUGIN_ID_GALATA_HELPERS);
        return (_a = plugin === null || plugin === void 0 ? void 0 : plugin.notifications.notifications) !== null && _a !== void 0 ? _a : [];
    }
    off(event, listener) {
        const callback = this.listeners.get(listener);
        if (callback) {
            this.getPlugin(_tokens__WEBPACK_IMPORTED_MODULE_2__.PLUGIN_ID_GALATA_HELPERS)
                .then(plugin => {
                if (!plugin) {
                    return;
                }
                switch (event) {
                    case 'dialog':
                        plugin.dialogs.currentChanged.disconnect(callback);
                        break;
                    case 'notification':
                        plugin.notifications.changed.disconnect(callback);
                        break;
                }
                this.listeners.delete(listener);
            })
                .catch(reason => {
                console.log(`Failed to disconnect listener for '${event}' event.`, reason);
            });
        }
    }
    /**
     * Connect a listener to new or updated Jupyter notification events.
     *
     * @param event Event type
     * @param listener Event listener
     */
    on(event, listener) {
        this.getPlugin(_tokens__WEBPACK_IMPORTED_MODULE_2__.PLUGIN_ID_GALATA_HELPERS)
            .then(plugin => {
            switch (event) {
                case 'dialog':
                    {
                        const callback = (tracker, dialog) => {
                            listener(dialog);
                        };
                        this.listeners.set(listener, callback);
                        plugin === null || plugin === void 0 ? void 0 : plugin.dialogs.currentChanged.connect(callback);
                    }
                    break;
                case 'notification':
                    {
                        const callback = (manager, notification) => {
                            if (notification.type !== 'removed') {
                                listener(notification.notification);
                            }
                        };
                        this.listeners.set(listener, callback);
                        plugin === null || plugin === void 0 ? void 0 : plugin.notifications.changed.connect(callback);
                    }
                    break;
            }
        })
            .catch(reason => {
            console.error(`Failed to add listener to JupyterLab dialog event:\n${reason}`);
        });
    }
    once(event, listener) {
        const onceListener = (arg) => {
            try {
                listener(arg);
            }
            finally {
                this.off(event, onceListener);
            }
        };
        this.on(event, onceListener);
    }
    /**
     * Wait for a function to finish for max. timeout milliseconds
     *
     * @param fn Function
     * @param timeout Timeout
     */
    async waitForFunction(fn, timeout) {
        return new Promise((resolve, reject) => {
            let checkTimer = null;
            let timeoutTimer = null;
            const check = async () => {
                checkTimer = null;
                if (await Promise.resolve(fn())) {
                    if (timeoutTimer) {
                        clearTimeout(timeoutTimer);
                    }
                    resolve();
                }
                else {
                    checkTimer = window.setTimeout(check, 200);
                }
            };
            void check();
            if (timeout) {
                timeoutTimer = window.setTimeout(() => {
                    timeoutTimer = null;
                    if (checkTimer) {
                        clearTimeout(checkTimer);
                    }
                    reject(new Error('Timed out waiting for condition to be fulfilled.'));
                }, timeout);
            }
        });
    }
    /**
     * Waits for the given `timeout` in milliseconds.
     *
     * @param timeout A timeout to wait for
     */
    async waitForTimeout(timeout) {
        return new Promise(resolve => {
            setTimeout(() => {
                resolve();
            }, timeout);
        });
    }
    /**
     * Wait for a condition to fulfill or for a certain time
     *
     * @param condition Condition or timeout to wait for
     * @param timeout Timeout
     */
    async waitFor(condition, timeout) {
        const conditionType = typeof condition;
        if (conditionType === 'function') {
            return this.waitForFunction(condition, timeout);
        }
        else if (conditionType === 'number') {
            return this.waitForTimeout(condition);
        }
    }
    /**
     * Wait for the route to be on path and close all documents
     *
     * @param path Path to monitor
     */
    async waitForLaunch(path = '/lab') {
        let resolver;
        const delegate = new Promise(resolve => {
            resolver = resolve;
        });
        const router = await this.getPlugin(PLUGIN_ID_ROUTER);
        const docManager = await this.getPlugin(PLUGIN_ID_DOC_MANAGER);
        const listener = async (sender, args) => {
            if (args.path === path) {
                router === null || router === void 0 ? void 0 : router.routed.disconnect(listener);
                await (docManager === null || docManager === void 0 ? void 0 : docManager.closeAll());
                resolver();
            }
        };
        router === null || router === void 0 ? void 0 : router.routed.connect(listener);
        return delegate;
    }
    /**
     * Wait for an element to be found from a CSS selector
     *
     * @param selector CSS selector
     * @param node Element
     * @param options Options
     * @returns Selected element
     */
    async waitForSelector(selector, node, options) {
        const waitForHidden = options && options.hidden;
        return new Promise((resolve, reject) => {
            const timer = setInterval(() => {
                const parent = node || document;
                const found = parent.querySelector(selector);
                if (waitForHidden) {
                    if (!found) {
                        clearInterval(timer);
                        resolve(null);
                    }
                }
                else if (found) {
                    clearInterval(timer);
                    resolve(found);
                }
            }, 200);
        });
    }
    /**
     * Wait for an element to be found from a XPath
     *
     * @param selector CSS selector
     * @param node Element
     * @param options Options
     * @returns Selected element
     */
    async waitForXPath(selector, node, options) {
        const waitForHidden = options && options.hidden;
        return new Promise((resolve, reject) => {
            const timer = setInterval(() => {
                const parent = node || document;
                const iterator = document.evaluate(selector, parent, null, XPathResult.UNORDERED_NODE_ITERATOR_TYPE, null);
                const found = iterator && iterator.iterateNext();
                if (waitForHidden) {
                    if (!found) {
                        clearInterval(timer);
                        resolve(null);
                    }
                }
                else if (found) {
                    clearInterval(timer);
                    resolve(found);
                }
            }, 200);
        });
    }
    /**
     * Delete all cells of the active notebook
     */
    async deleteNotebookCells() {
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        void this._app.commands.execute('notebook:delete-cell');
        nb.update();
    }
    /**
     * Add a cell to the active notebook
     *
     * @param cellType Cell type
     * @param source Cell input source
     * @returns Action success result
     */
    addNotebookCell(cellType, source) {
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        if (nb !== null) {
            void this._app.commands.execute('notebook:insert-cell-below');
            if (nb.model) {
                const sharedModel = nb.model.sharedModel;
                sharedModel.insertCell(sharedModel.cells.length, {
                    cell_type: cellType,
                    source
                });
            }
            nb.update();
        }
        else {
            return false;
        }
        return true;
    }
    /**
     * Reset execution count of one or all cells.
     *
     * @param cellIndex Cell index
     */
    resetExecutionCount(cellIndex) {
        var _a;
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        if (nb) {
            if (nb.model) {
                if (cellIndex === undefined) {
                    for (const cell of nb.model.cells) {
                        switch (cell.type) {
                            case 'code':
                                cell.executionCount = null;
                                break;
                        }
                    }
                }
                else if (((_a = nb.model.cells.get(cellIndex)) === null || _a === void 0 ? void 0 : _a.type) === 'code') {
                    nb.model.cells.get(cellIndex).executionCount =
                        null;
                }
            }
        }
    }
    /**
     * Set the type and content of a cell in the active notebook
     *
     * @param cellIndex Cell index
     * @param cellType Cell type
     * @param source Cell input source
     * @returns Action success status
     */
    setNotebookCell(cellIndex, cellType, source) {
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        if (nb) {
            const numCells = nb.widgets.length;
            if (cellIndex < 0 || cellIndex >= numCells) {
                return false;
            }
            if (nb.model) {
                const sharedModel = nb.model.sharedModel;
                sharedModel.transact(() => {
                    sharedModel.deleteCell(cellIndex);
                    sharedModel.insertCell(cellIndex, { cell_type: cellType, source });
                });
            }
            nb.update();
        }
        else {
            return false;
        }
        return true;
    }
    /**
     * Test if a cell is selected in the active notebook
     *
     * @param cellIndex Cell index
     * @returns Whether the cell is selected or not
     */
    isNotebookCellSelected(cellIndex) {
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        const numCells = nb.widgets.length;
        if (cellIndex < 0 || cellIndex >= numCells) {
            return false;
        }
        return nb.isSelected(nb.widgets[cellIndex]);
    }
    /**
     * Save the active notebook
     */
    async saveActiveNotebook() {
        const nbPanel = this._app.shell.currentWidget;
        await nbPanel.context.save();
    }
    /**
     * Wait for a Markdown cell to be rendered
     *
     * @param cell Cell
     */
    async waitForMarkdownCellRendered(cell) {
        if (!cell.inViewport) {
            return;
        }
        await cell.ready;
        if (cell.rendered) {
            return;
        }
        let resolver;
        const delegate = new Promise(resolve => {
            resolver = resolve;
        });
        const onRenderedChanged = () => {
            _lumino_signaling__WEBPACK_IMPORTED_MODULE_1__.Signal.disconnectReceiver(onRenderedChanged);
            resolver();
        };
        cell.renderedChanged.connect(onRenderedChanged);
        return delegate;
    }
    /**
     * Wait for a cell to be run
     *
     * @param cell Cell
     * @param timeout Timeout
     */
    async waitForCellRun(cell, timeout = 2000) {
        const model = cell.model;
        const code = model.sharedModel.getSource();
        if (!code.trim()) {
            return;
        }
        const cellType = cell.model.type;
        if (cellType === 'markdown') {
            await this.waitForMarkdownCellRendered(cell);
            return;
        }
        else if (cellType === 'raw') {
            return;
        }
        else {
            // 'code'
            let resolver;
            const delegate = new Promise(resolve => {
                resolver = resolve;
            });
            let numTries = 0;
            let timer = null;
            let timeoutTimer = null;
            const clearAndResolve = () => {
                clearInterval(timer);
                timer = null;
                clearTimeout(timeoutTimer);
                timeoutTimer = null;
                resolver();
            };
            const startTimeout = () => {
                if (!timeoutTimer) {
                    timeoutTimer = setTimeout(() => {
                        clearAndResolve();
                    }, timeout);
                }
            };
            const checkIfDone = () => {
                if (cell.model.executionCount !== null) {
                    if (!cell.inViewport) {
                        clearAndResolve();
                        return;
                    }
                    const output = cell.node.querySelector('.jp-Cell-outputArea .jp-OutputArea-output');
                    if (output) {
                        if (output.textContent === 'Loading widget...') {
                            startTimeout();
                        }
                        else {
                            clearAndResolve();
                        }
                    }
                    else {
                        if (numTries > 0) {
                            clearAndResolve();
                        }
                    }
                }
                else {
                    startTimeout();
                }
                numTries++;
            };
            checkIfDone();
            timer = setInterval(() => {
                checkIfDone();
            }, 200);
            return delegate;
        }
    }
    /**
     * Run the active notebook cell by cell
     * and execute the callback after each cell execution
     *
     * @param callback Callback
     */
    async runActiveNotebookCellByCell(callback) {
        const nbPanel = this._app.shell.currentWidget;
        const notebook = nbPanel.content;
        if (!notebook.widgets) {
            console.error('NOTEBOOK CELL PROBLEM', notebook);
        }
        const numCells = notebook.widgets.length;
        if (numCells === 0) {
            return;
        }
        void this._app.commands.execute('notebook:deselect-all');
        for (let i = 0; i < numCells; ++i) {
            const cell = notebook.widgets[i];
            notebook.activeCellIndex = i;
            notebook.select(cell);
            await this._app.commands.execute('notebook:run-cell');
            await this.waitForCellRun(cell);
            if (callback === null || callback === void 0 ? void 0 : callback.onAfterCellRun) {
                await callback.onAfterCellRun(i);
            }
            if (callback === null || callback === void 0 ? void 0 : callback.onBeforeScroll) {
                await callback.onBeforeScroll();
            }
            await notebook.scrollToItem(i).catch(reason => {
                // no-op
            });
            if (callback === null || callback === void 0 ? void 0 : callback.onAfterScroll) {
                await callback.onAfterScroll();
            }
        }
    }
    /**
     * Test if one or all cells have an execution number.
     *
     * @param cellIndex Cell index
     * @returns Whether the cell was executed or not
     *
     * ### Notes
     * It checks that no cells have a `null` execution count.
     */
    haveBeenExecuted(cellIndex) {
        const nbPanel = this._app.shell.currentWidget;
        const nb = nbPanel.content;
        let counter = 0;
        if (nb) {
            if (nb.model) {
                if (cellIndex === undefined) {
                    for (const cell of nb.model.cells) {
                        if (cell.type === 'code') {
                            counter +=
                                cell.sharedModel.getSource().length > 0 &&
                                    cell.executionCount === null
                                    ? 1
                                    : 0;
                        }
                    }
                }
                else {
                    const cell = nb.model.cells.get(cellIndex);
                    if ((cell === null || cell === void 0 ? void 0 : cell.type) === 'code') {
                        counter +=
                            cell.sharedModel.getSource().length > 0 &&
                                cell.executionCount === null
                                ? 1
                                : 0;
                    }
                }
            }
        }
        return counter === 0;
    }
    /**
     * Get the index of a toolbar item
     *
     * @param itemName Item name
     * @returns Index
     */
    getNotebookToolbarItemIndex(itemName) {
        const nbPanel = this._app.shell.currentWidget;
        return (0,_lumino_algorithm__WEBPACK_IMPORTED_MODULE_0__.findIndex)(nbPanel.toolbar.names(), name => name === itemName);
    }
    /**
     * Test if a element is visible or not
     *
     * @param el Element
     * @returns Test result
     */
    isElementVisible(el) {
        return !!(el.offsetWidth || el.offsetHeight || el.getClientRects().length);
    }
    /**
     * Set the application theme
     *
     * @param themeName Theme name
     */
    async setTheme(themeName) {
        await this._app.commands.execute('apputils:change-theme', {
            theme: themeName
        });
        await this.waitFor(async () => {
            return document.body.dataset.jpThemeName === themeName;
        });
    }
    /**
     * Application object
     */
    get app() {
        if (!this._app) {
            this._app = window.jupyterapp;
        }
        return this._app;
    }
}


/***/ }),

/***/ "../lib/extension/index.js":
/*!*********************************!*\
  !*** ../lib/extension/index.js ***!
  \*********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   IGalataHelpers: () => (/* reexport safe */ _tokens__WEBPACK_IMPORTED_MODULE_2__.IGalataHelpers),
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyterlab/apputils */ "webpack/sharing/consume/default/@jupyterlab/apputils");
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _global__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./global */ "../lib/extension/global.js");
/* harmony import */ var _tokens__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./tokens */ "../lib/extension/tokens.js");
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.




/**
 * Galata in-page object
 */
window.galata = new _global__WEBPACK_IMPORTED_MODULE_1__.GalataInpage();
/**
 * Galata in-page object
 * @deprecated
 */
window.galataip = window.galata;
const galataPlugin = {
    id: _tokens__WEBPACK_IMPORTED_MODULE_2__.PLUGIN_ID_GALATA_HELPERS,
    autoStart: true,
    activate: (app) => {
        return Object.freeze({
            notifications: _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__.Notification.manager,
            dialogs: _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__.Dialog.tracker
        });
    }
};
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (galataPlugin);


/***/ }),

/***/ "../lib/extension/tokens.js":
/*!**********************************!*\
  !*** ../lib/extension/tokens.js ***!
  \**********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   IGalataHelpers: () => (/* binding */ IGalataHelpers),
/* harmony export */   PLUGIN_ID_GALATA_HELPERS: () => (/* binding */ PLUGIN_ID_GALATA_HELPERS)
/* harmony export */ });
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @lumino/coreutils */ "webpack/sharing/consume/default/@lumino/coreutils");
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__);
/*
 * Copyright (c) Jupyter Development Team.
 * Distributed under the terms of the Modified BSD License.
 * Copyright (c) Bloomberg Finance LP.
 */

/**
 * Galata in-page extension helpers.
 */
const PLUGIN_ID_GALATA_HELPERS = '@jupyterlab/galata-extension:helpers';
/**
 * Test token exposing some static JupyterLab objects.
 */
const IGalataHelpers = new _lumino_coreutils__WEBPACK_IMPORTED_MODULE_0__.Token('@jupyterlab/galata-extension:IGalataHelpers');


/***/ })

}]);
//# sourceMappingURL=lib_extension_index_js.cc9baf28438ab1b21a09.js.map