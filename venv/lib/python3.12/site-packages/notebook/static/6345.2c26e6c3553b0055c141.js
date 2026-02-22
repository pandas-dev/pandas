"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[6345,4645],{

/***/ 94645:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(21217);
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _jupyterlab_console__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11798);
/* harmony import */ var _jupyterlab_console__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_console__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(52733);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(83228);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(78839);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.





/**
 * A plugin to open consoles in a new tab
 */
const opener = {
    id: '@jupyter-notebook/console-extension:opener',
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter],
    autoStart: true,
    description: 'A plugin to open consoles in a new tab',
    activate: (app, router) => {
        const { commands } = app;
        const consolePattern = new RegExp('/consoles/(.*)');
        const command = 'router:console';
        commands.addCommand(command, {
            execute: (args) => {
                const parsed = args;
                const matches = parsed.path.match(consolePattern);
                if (!matches) {
                    return;
                }
                const [, match] = matches;
                if (!match) {
                    return;
                }
                const path = decodeURIComponent(match);
                commands.execute('console:create', { path });
            },
        });
        router.register({ command, pattern: consolePattern });
    },
};
/**
 * Open consoles in a new tab.
 */
const redirect = {
    id: '@jupyter-notebook/console-extension:redirect',
    requires: [_jupyterlab_console__WEBPACK_IMPORTED_MODULE_1__.IConsoleTracker],
    optional: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__.INotebookPathOpener],
    autoStart: true,
    description: 'Open consoles in a new tab',
    activate: (app, tracker, notebookPathOpener) => {
        const baseUrl = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_2__.PageConfig.getBaseUrl();
        const opener = notebookPathOpener !== null && notebookPathOpener !== void 0 ? notebookPathOpener : _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__.defaultNotebookPathOpener;
        tracker.widgetAdded.connect(async (send, console) => {
            const { sessionContext } = console;
            await sessionContext.ready;
            const widget = (0,_lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__.find)(app.shell.widgets('main'), (w) => w.id === console.id);
            if (widget) {
                // bail if the console is already added to the main area
                return;
            }
            opener.open({
                prefix: _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_2__.URLExt.join(baseUrl, 'consoles'),
                path: sessionContext.path,
                target: '_blank',
            });
            // the widget is not needed anymore
            console.dispose();
        });
    },
};
/**
 * Export the plugins as default.
 */
const plugins = [opener, redirect];
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (plugins);


/***/ })

}]);
//# sourceMappingURL=6345.2c26e6c3553b0055c141.js.map?v=2c26e6c3553b0055c141