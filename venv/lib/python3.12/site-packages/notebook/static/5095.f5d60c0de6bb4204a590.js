"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5095],{

/***/ 26167:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   e: () => (/* binding */ __name)
/* harmony export */ });
var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });




/***/ }),

/***/ 5095:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (/* binding */ layouts_default)
/* harmony export */ });
/* harmony import */ var _chunks_mermaid_layout_elk_core_chunk_ZW26E7AF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(26167);


// src/layouts.ts
var loader = /* @__PURE__ */ (0,_chunks_mermaid_layout_elk_core_chunk_ZW26E7AF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .e)(async () => await Promise.all(/* import() */[__webpack_require__.e(6577), __webpack_require__.e(5321), __webpack_require__.e(4915)]).then(__webpack_require__.bind(__webpack_require__, 14915)), "loader");
var algos = ["elk.stress", "elk.force", "elk.mrtree", "elk.sporeOverlap"];
var layouts = [
  {
    name: "elk",
    loader,
    algorithm: "elk.layered"
  },
  ...algos.map((algo) => ({
    name: algo,
    loader,
    algorithm: algo
  }))
];
var layouts_default = layouts;



/***/ })

}]);
//# sourceMappingURL=5095.f5d60c0de6bb4204a590.js.map?v=f5d60c0de6bb4204a590