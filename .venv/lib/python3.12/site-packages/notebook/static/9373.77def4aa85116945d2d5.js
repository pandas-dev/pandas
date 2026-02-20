"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[9373],{

/***/ 99373:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_6MN3ZHY7_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(50692);
/* harmony import */ var _chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(98154);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(74999);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(13197);





// src/diagrams/info/infoParser.ts

var parser = {
  parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__/* .parse */ .Qc)("info", input);
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(ast);
  }, "parse")
};

// src/diagrams/info/infoDb.ts
var DEFAULT_INFO_DB = {
  version: _chunk_6MN3ZHY7_mjs__WEBPACK_IMPORTED_MODULE_0__/* .package_default */ .X.version + ( true ? "" : 0)
};
var getVersion = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => DEFAULT_INFO_DB.version, "getVersion");
var db = {
  getVersion
};

// src/diagrams/info/infoRenderer.ts
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((text, id, version) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("rendering info diagram\n" + text);
  const svg = (0,_chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_1__/* .selectSvgElement */ .P)(id);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .configureSvgSize */ .v2)(svg, 100, 400, true);
  const group = svg.append("g");
  group.append("text").attr("x", 100).attr("y", 40).attr("class", "version").attr("font-size", 32).style("text-anchor", "middle").text(`v${version}`);
}, "draw");
var renderer = { draw };

// src/diagrams/info/infoDiagram.ts
var diagram = {
  parser,
  db,
  renderer
};



/***/ })

}]);
//# sourceMappingURL=9373.77def4aa85116945d2d5.js.map?v=77def4aa85116945d2d5