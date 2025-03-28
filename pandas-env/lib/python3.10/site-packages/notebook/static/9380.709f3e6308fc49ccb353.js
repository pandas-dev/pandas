"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[9380],{

/***/ 19380:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(11363);
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(55985);
/* harmony import */ var _jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(7801);
/* harmony import */ var _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_translation__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(4174);
/* harmony import */ var _jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var react__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(78156);
/* harmony import */ var react__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_4__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.





/**
 * A list of resources to show in the help menu.
 */
const RESOURCES = [
    {
        text: 'About Jupyter',
        url: 'https://jupyter.org',
    },
    {
        text: 'Markdown Reference',
        url: 'https://commonmark.org/help/',
    },
    {
        text: 'Documentation',
        url: 'https://jupyter-notebook.readthedocs.io/en/stable/',
    },
];
/**
 * The command IDs used by the help plugin.
 */
var CommandIDs;
(function (CommandIDs) {
    CommandIDs.open = 'help:open';
    CommandIDs.about = 'help:about';
})(CommandIDs || (CommandIDs = {}));
/**
 * A plugin to open the about section with resources.
 */
const open = {
    id: '@jupyter-notebook/help-extension:open',
    autoStart: true,
    description: 'A plugin to open the about section with resources',
    activate: (app) => {
        const { commands } = app;
        commands.addCommand(CommandIDs.open, {
            label: (args) => args['text'],
            execute: (args) => {
                const url = args['url'];
                window.open(url);
            },
        });
    },
};
/**
 * Plugin to add a command to show an About Jupyter Notebook and Markdown Reference.
 */
const about = {
    id: '@jupyter-notebook/help-extension:about',
    autoStart: true,
    requires: [_jupyterlab_translation__WEBPACK_IMPORTED_MODULE_2__.ITranslator],
    optional: [_jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_1__.IMainMenu, _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__.ICommandPalette],
    description: 'Plugin to add a command to show an About Jupyter Notebook and Markdown Reference',
    activate: (app, translator, menu, palette) => {
        const { commands } = app;
        const trans = translator.load('notebook');
        const category = trans.__('Help');
        commands.addCommand(CommandIDs.about, {
            label: trans.__('About %1', app.name),
            execute: () => {
                const title = (react__WEBPACK_IMPORTED_MODULE_4__.createElement(react__WEBPACK_IMPORTED_MODULE_4__.Fragment, null,
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("span", { className: "jp-AboutNotebook-header" },
                        react__WEBPACK_IMPORTED_MODULE_4__.createElement(_jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_3__.jupyterIcon.react, { width: "196px", height: "auto" }))));
                const notebookURL = 'https://github.com/jupyter/notebook';
                const contributorURL = 'https://github.com/jupyter/notebook/pulse';
                const aboutJupyter = trans.__('JUPYTER NOTEBOOK ON GITHUB');
                const contributorList = trans.__('CONTRIBUTOR LIST');
                const externalLinks = (react__WEBPACK_IMPORTED_MODULE_4__.createElement("span", null,
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("a", { href: notebookURL, target: "_blank", rel: "noopener noreferrer", className: "jp-Button-flat jp-AboutNotebook-about-externalLinks" }, aboutJupyter),
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("a", { href: contributorURL, target: "_blank", rel: "noopener noreferrer", className: "jp-Button-flat jp-AboutNotebook-about-externalLinks" }, contributorList)));
                const version = trans.__('Version: %1', app.version);
                const copyright = trans.__('© 2021-2023 Jupyter Notebook Contributors');
                const body = (react__WEBPACK_IMPORTED_MODULE_4__.createElement(react__WEBPACK_IMPORTED_MODULE_4__.Fragment, null,
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("span", { className: "jp-AboutNotebook-version" }, version),
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("div", null, externalLinks),
                    react__WEBPACK_IMPORTED_MODULE_4__.createElement("span", { className: "jp-AboutNotebook-about-copyright" }, copyright)));
                const dialog = new _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__.Dialog({
                    title,
                    body,
                    buttons: [
                        _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_0__.Dialog.createButton({
                            label: trans.__('Dismiss'),
                            className: 'jp-AboutNotebook-about-button jp-mod-reject jp-mod-styled',
                        }),
                    ],
                });
                dialog.addClass('jp-AboutNotebook');
                void dialog.launch();
            },
        });
        if (palette) {
            palette.addItem({ command: CommandIDs.about, category });
        }
        const resourcesGroup = RESOURCES.map((args) => ({
            args,
            command: CommandIDs.open,
        }));
        if (menu) {
            menu.helpMenu.addGroup(resourcesGroup, 30);
        }
    },
};
const plugins = [open, about];
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (plugins);


/***/ })

}]);
//# sourceMappingURL=9380.709f3e6308fc49ccb353.js.map?v=709f3e6308fc49ccb353