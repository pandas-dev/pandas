"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1616],{

/***/ 18556:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ populateCommonDb)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


// src/diagrams/common/populateCommonDb.ts
function populateCommonDb(ast, db) {
  if (ast.accDescr) {
    db.setAccDescription?.(ast.accDescr);
  }
  if (ast.accTitle) {
    db.setAccTitle?.(ast.accTitle);
  }
  if (ast.title) {
    db.setDiagramTitle?.(ast.title);
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(populateCommonDb, "populateCommonDb");




/***/ }),

/***/ 1616:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(98154);
/* harmony import */ var _chunk_4BX2VUAB_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(18556);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(74999);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(13197);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(35321);






// src/diagrams/pie/pieParser.ts


// src/diagrams/pie/pieDb.ts
var DEFAULT_PIE_CONFIG = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .defaultConfig_default */ .vZ.pie;
var DEFAULT_PIE_DB = {
  sections: /* @__PURE__ */ new Map(),
  showData: false,
  config: DEFAULT_PIE_CONFIG
};
var sections = DEFAULT_PIE_DB.sections;
var showData = DEFAULT_PIE_DB.showData;
var config = structuredClone(DEFAULT_PIE_CONFIG);
var getConfig2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(() => structuredClone(config), "getConfig");
var clear2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(() => {
  sections = /* @__PURE__ */ new Map();
  showData = DEFAULT_PIE_DB.showData;
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
}, "clear");
var addSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(({ label, value }) => {
  if (value < 0) {
    throw new Error(
      `"${label}" has invalid value: ${value}. Negative values are not allowed in pie charts. All slice values must be >= 0.`
    );
  }
  if (!sections.has(label)) {
    sections.set(label, value);
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.debug(`added new section: ${label}, with value: ${value}`);
  }
}, "addSection");
var getSections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(() => sections, "getSections");
var setShowData = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((toggle) => {
  showData = toggle;
}, "setShowData");
var getShowData = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(() => showData, "getShowData");
var db = {
  getConfig: getConfig2,
  clear: clear2,
  setDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr,
  setAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu,
  setAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$,
  getAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx,
  addSection,
  getSections,
  setShowData,
  getShowData
};

// src/diagrams/pie/pieParser.ts
var populateDb = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((ast, db2) => {
  (0,_chunk_4BX2VUAB_mjs__WEBPACK_IMPORTED_MODULE_1__/* .populateCommonDb */ .A)(ast, db2);
  db2.setShowData(ast.showData);
  ast.sections.map(db2.addSection);
}, "populateDb");
var parser = {
  parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_5__/* .parse */ .Qc)("pie", input);
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.debug(ast);
    populateDb(ast, db);
  }, "parse")
};

// src/diagrams/pie/pieStyles.ts
var getStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((options) => `
  .pieCircle{
    stroke: ${options.pieStrokeColor};
    stroke-width : ${options.pieStrokeWidth};
    opacity : ${options.pieOpacity};
  }
  .pieOuterCircle{
    stroke: ${options.pieOuterStrokeColor};
    stroke-width: ${options.pieOuterStrokeWidth};
    fill: none;
  }
  .pieTitleText {
    text-anchor: middle;
    font-size: ${options.pieTitleTextSize};
    fill: ${options.pieTitleTextColor};
    font-family: ${options.fontFamily};
  }
  .slice {
    font-family: ${options.fontFamily};
    fill: ${options.pieSectionTextColor};
    font-size:${options.pieSectionTextSize};
    // fill: white;
  }
  .legend text {
    fill: ${options.pieLegendTextColor};
    font-family: ${options.fontFamily};
    font-size: ${options.pieLegendTextSize};
  }
`, "getStyles");
var pieStyles_default = getStyles;

// src/diagrams/pie/pieRenderer.ts

var createPieArcs = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((sections2) => {
  const sum = [...sections2.values()].reduce((acc, val) => acc + val, 0);
  const pieData = [...sections2.entries()].map(([label, value]) => ({ label, value })).filter((d) => d.value / sum * 100 >= 1).sort((a, b) => b.value - a.value);
  const pie = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .pie */ .ve8)().value((d) => d.value);
  return pie(pieData);
}, "createPieArcs");
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((text, id, _version, diagObj) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.debug("rendering pie chart\n" + text);
  const db2 = diagObj.db;
  const globalConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
  const pieConfig = (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__/* .cleanAndMerge */ .Rb)(db2.getConfig(), globalConfig.pie);
  const MARGIN = 40;
  const LEGEND_RECT_SIZE = 18;
  const LEGEND_SPACING = 4;
  const height = 450;
  const pieWidth = height;
  const svg = (0,_chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__/* .selectSvgElement */ .P)(id);
  const group = svg.append("g");
  group.attr("transform", "translate(" + pieWidth / 2 + "," + height / 2 + ")");
  const { themeVariables } = globalConfig;
  let [outerStrokeWidth] = (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__/* .parseFontSize */ .VG)(themeVariables.pieOuterStrokeWidth);
  outerStrokeWidth ??= 2;
  const textPosition = pieConfig.textPosition;
  const radius = Math.min(pieWidth, height) / 2 - MARGIN;
  const arcGenerator = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .arc */ .Nb1)().innerRadius(0).outerRadius(radius);
  const labelArcGenerator = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .arc */ .Nb1)().innerRadius(radius * textPosition).outerRadius(radius * textPosition);
  group.append("circle").attr("cx", 0).attr("cy", 0).attr("r", radius + outerStrokeWidth / 2).attr("class", "pieOuterCircle");
  const sections2 = db2.getSections();
  const arcs = createPieArcs(sections2);
  const myGeneratedColors = [
    themeVariables.pie1,
    themeVariables.pie2,
    themeVariables.pie3,
    themeVariables.pie4,
    themeVariables.pie5,
    themeVariables.pie6,
    themeVariables.pie7,
    themeVariables.pie8,
    themeVariables.pie9,
    themeVariables.pie10,
    themeVariables.pie11,
    themeVariables.pie12
  ];
  let sum = 0;
  sections2.forEach((section) => {
    sum += section;
  });
  const filteredArcs = arcs.filter((datum) => (datum.data.value / sum * 100).toFixed(0) !== "0");
  const color = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .scaleOrdinal */ .PKp)(myGeneratedColors);
  group.selectAll("mySlices").data(filteredArcs).enter().append("path").attr("d", arcGenerator).attr("fill", (datum) => {
    return color(datum.data.label);
  }).attr("class", "pieCircle");
  group.selectAll("mySlices").data(filteredArcs).enter().append("text").text((datum) => {
    return (datum.data.value / sum * 100).toFixed(0) + "%";
  }).attr("transform", (datum) => {
    return "translate(" + labelArcGenerator.centroid(datum) + ")";
  }).style("text-anchor", "middle").attr("class", "slice");
  group.append("text").text(db2.getDiagramTitle()).attr("x", 0).attr("y", -(height - 50) / 2).attr("class", "pieTitleText");
  const allSectionData = [...sections2.entries()].map(([label, value]) => ({
    label,
    value
  }));
  const legend = group.selectAll(".legend").data(allSectionData).enter().append("g").attr("class", "legend").attr("transform", (_datum, index) => {
    const height2 = LEGEND_RECT_SIZE + LEGEND_SPACING;
    const offset = height2 * allSectionData.length / 2;
    const horizontal = 12 * LEGEND_RECT_SIZE;
    const vertical = index * height2 - offset;
    return "translate(" + horizontal + "," + vertical + ")";
  });
  legend.append("rect").attr("width", LEGEND_RECT_SIZE).attr("height", LEGEND_RECT_SIZE).style("fill", (d) => color(d.label)).style("stroke", (d) => color(d.label));
  legend.append("text").attr("x", LEGEND_RECT_SIZE + LEGEND_SPACING).attr("y", LEGEND_RECT_SIZE - LEGEND_SPACING).text((d) => {
    if (db2.getShowData()) {
      return `${d.label} [${d.value}]`;
    }
    return d.label;
  });
  const longestTextWidth = Math.max(
    ...legend.selectAll("text").nodes().map((node) => node?.getBoundingClientRect().width ?? 0)
  );
  const totalWidth = pieWidth + MARGIN + LEGEND_RECT_SIZE + LEGEND_SPACING + longestTextWidth;
  svg.attr("viewBox", `0 0 ${totalWidth} ${height}`);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .configureSvgSize */ .v2)(svg, height, totalWidth, pieConfig.useMaxWidth);
}, "draw");
var renderer = { draw };

// src/diagrams/pie/pieDiagram.ts
var diagram = {
  parser,
  db,
  renderer,
  styles: pieStyles_default
};



/***/ })

}]);
//# sourceMappingURL=1616.ee161d92c1ef1d77afcc.js.map?v=ee161d92c1ef1d77afcc