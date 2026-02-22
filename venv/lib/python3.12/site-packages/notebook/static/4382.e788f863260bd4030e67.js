"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4382,7906],{

/***/ 54382:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_documentsearch__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(79331);
/* harmony import */ var _jupyterlab_documentsearch__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_documentsearch__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(83228);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_1__);


const SEARCHABLE_CLASS = 'jp-mod-searchable';
/**
 * A plugin to add document search functionalities.
 */
const notebookShellWidgetListener = {
    id: '@jupyter-notebook/documentsearch-extension:notebookShellWidgetListener',
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_1__.INotebookShell, _jupyterlab_documentsearch__WEBPACK_IMPORTED_MODULE_0__.ISearchProviderRegistry],
    autoStart: true,
    description: 'A plugin to add document search functionalities',
    activate: (app, notebookShell, registry) => {
        // If a given widget is searchable, apply the searchable class.
        // If it's not searchable, remove the class.
        const transformWidgetSearchability = (widget) => {
            if (!widget) {
                return;
            }
            if (registry.hasProvider(widget)) {
                widget.addClass(SEARCHABLE_CLASS);
            }
            else {
                widget.removeClass(SEARCHABLE_CLASS);
            }
        };
        // Update searchability of the active widget when the registry
        // changes, in case a provider for the current widget was added
        // or removed
        registry.changed.connect(() => transformWidgetSearchability(notebookShell.currentWidget));
        // Apply the searchable class only to the active widget if it is actually
        // searchable. Remove the searchable class from a widget when it's
        // no longer active.
        notebookShell.currentChanged.connect((_, args) => {
            if (notebookShell.currentWidget) {
                transformWidgetSearchability(notebookShell.currentWidget);
            }
        });
    },
};
/**
 * Export the plugins as default.
 */
const plugins = [notebookShellWidgetListener];
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (plugins);


/***/ })

}]);
//# sourceMappingURL=4382.e788f863260bd4030e67.js.map?v=e788f863260bd4030e67