"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4837,3146],{

/***/ 73146:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  INotebookTree: () => (/* reexport */ INotebookTree),
  NotebookTreeWidget: () => (/* reexport */ NotebookTreeWidget)
});

// EXTERNAL MODULE: consume shared module (default) @jupyterlab/ui-components@~4.3.6 (singleton) (fallback: ../node_modules/@jupyterlab/ui-components/lib/index.js)
var index_js_ = __webpack_require__(13077);
// EXTERNAL MODULE: consume shared module (default) @lumino/widgets@~2.5.0 (singleton) (fallback: ../node_modules/@lumino/widgets/dist/index.es6.js)
var index_es6_js_ = __webpack_require__(2260);
;// CONCATENATED MODULE: ../packages/tree/lib/notebook-tree.js


/**
 * The widget added in main area of the tree view.
 */
class NotebookTreeWidget extends index_es6_js_.TabPanel {
    /**
     * Constructor of the NotebookTreeWidget.
     */
    constructor() {
        super({
            tabPlacement: 'top',
            tabsMovable: true,
            renderer: index_js_.TabBarSvg.defaultRenderer,
        });
        this.addClass('jp-TreePanel');
    }
}

// EXTERNAL MODULE: consume shared module (default) @lumino/coreutils@~2.2.0 (singleton) (fallback: ../node_modules/@lumino/coreutils/dist/index.js)
var dist_index_js_ = __webpack_require__(21961);
;// CONCATENATED MODULE: ../packages/tree/lib/token.js

/**
 * The INotebookTree token.
 */
const INotebookTree = new dist_index_js_.Token('@jupyter-notebook/tree:INotebookTree');

;// CONCATENATED MODULE: ../packages/tree/lib/index.js




/***/ })

}]);
//# sourceMappingURL=4837.84a75547bf8a4a8e3759.js.map?v=84a75547bf8a4a8e3759