"use strict";
(self["webpackChunk_jupyter_widgets_jupyterlab_manager"] = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] || []).push([["packages_output_lib_index_js"],{

/***/ "../../packages/output/lib/index.js":
/*!******************************************!*\
  !*** ../../packages/output/lib/index.js ***!
  \******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   OUTPUT_WIDGET_VERSION: () => (/* reexport safe */ _output__WEBPACK_IMPORTED_MODULE_0__.OUTPUT_WIDGET_VERSION),
/* harmony export */   OutputModel: () => (/* reexport safe */ _output__WEBPACK_IMPORTED_MODULE_0__.OutputModel),
/* harmony export */   OutputView: () => (/* reexport safe */ _output__WEBPACK_IMPORTED_MODULE_0__.OutputView)
/* harmony export */ });
/* harmony import */ var _output__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./output */ "../../packages/output/lib/output.js");



/***/ }),

/***/ "../../packages/output/lib/output.js":
/*!*******************************************!*\
  !*** ../../packages/output/lib/output.js ***!
  \*******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   OUTPUT_WIDGET_VERSION: () => (/* binding */ OUTPUT_WIDGET_VERSION),
/* harmony export */   OutputModel: () => (/* binding */ OutputModel),
/* harmony export */   OutputView: () => (/* binding */ OutputView)
/* harmony export */ });
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @jupyter-widgets/base */ "webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base");
/* harmony import */ var _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__);

const OUTPUT_WIDGET_VERSION = '1.0.0';
class OutputModel extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetModel {
    defaults() {
        return Object.assign(Object.assign({}, super.defaults()), { _model_name: 'OutputModel', _view_name: 'OutputView', _model_module: '@jupyter-widgets/output', _view_module: '@jupyter-widgets/output', _model_module_version: OUTPUT_WIDGET_VERSION, _view_module_version: OUTPUT_WIDGET_VERSION });
    }
}
class OutputView extends _jupyter_widgets_base__WEBPACK_IMPORTED_MODULE_0__.DOMWidgetView {
}


/***/ })

}]);
//# sourceMappingURL=packages_output_lib_index_js.49c9e4037a3b9e9e3b18.js.map