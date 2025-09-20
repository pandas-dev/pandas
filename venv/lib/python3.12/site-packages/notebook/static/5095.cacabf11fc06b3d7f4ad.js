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
var loader = /* @__PURE__ */ (0,_chunks_mermaid_layout_elk_core_chunk_ZW26E7AF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .e)(async () => await Promise.all(/* import() */[__webpack_require__.e(1715), __webpack_require__.e(3619), __webpack_require__.e(9310)]).then(__webpack_require__.bind(__webpack_require__, 79310)), "loader");
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
//# sourceMappingURL=5095.cacabf11fc06b3d7f4ad.js.map?v=cacabf11fc06b3d7f4ad