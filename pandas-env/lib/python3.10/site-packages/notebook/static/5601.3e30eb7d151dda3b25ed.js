"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5601,1684],{

/***/ 95601:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(38075);
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(76107);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _jupyterlab_terminal__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(27368);
/* harmony import */ var _jupyterlab_terminal__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_terminal__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(87362);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(14931);
/* harmony import */ var _lumino_algorithm__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(_lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.





/**
 * A plugin to open terminals in a new tab
 */
const opener = {
    id: '@jupyter-notebook/terminal-extension:opener',
    description: 'A plugin to open terminals in a new tab.',
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter, _jupyterlab_terminal__WEBPACK_IMPORTED_MODULE_2__.ITerminalTracker],
    autoStart: true,
    activate: (app, router, tracker) => {
        const { commands } = app;
        const terminalPattern = new RegExp('/terminals/(.*)');
        const command = 'router:terminal';
        commands.addCommand(command, {
            execute: (args) => {
                const parsed = args;
                const matches = parsed.path.match(terminalPattern);
                if (!matches) {
                    return;
                }
                const [, name] = matches;
                if (!name) {
                    return;
                }
                tracker.widgetAdded.connect((send, terminal) => {
                    terminal.content.setOption('closeOnExit', false);
                });
                commands.execute('terminal:open', { name });
            },
        });
        router.register({ command, pattern: terminalPattern });
    },
};
/**
 * Open terminals in a new tab.
 */
const redirect = {
    id: '@jupyter-notebook/terminal-extension:redirect',
    description: 'Open terminals in a new tab.',
    requires: [_jupyterlab_terminal__WEBPACK_IMPORTED_MODULE_2__.ITerminalTracker],
    optional: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__.INotebookPathOpener],
    autoStart: true,
    activate: (app, tracker, notebookPathOpener) => {
        const baseUrl = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_1__.PageConfig.getBaseUrl();
        const opener = notebookPathOpener !== null && notebookPathOpener !== void 0 ? notebookPathOpener : _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_3__.defaultNotebookPathOpener;
        tracker.widgetAdded.connect((send, terminal) => {
            const widget = (0,_lumino_algorithm__WEBPACK_IMPORTED_MODULE_4__.find)(app.shell.widgets('main'), (w) => w.id === terminal.id);
            if (widget) {
                // bail if the terminal is already added to the main area
                return;
            }
            const name = terminal.content.session.name;
            opener.open({
                prefix: _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_1__.URLExt.join(baseUrl, 'terminals'),
                path: name,
                target: '_blank',
            });
            // dispose the widget since it is not used on this page
            terminal.dispose();
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
//# sourceMappingURL=5601.3e30eb7d151dda3b25ed.js.map?v=3e30eb7d151dda3b25ed