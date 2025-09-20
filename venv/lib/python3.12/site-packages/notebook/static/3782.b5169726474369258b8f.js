"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3782],{

/***/ 80910:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ populateCommonDb)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);


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
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(populateCommonDb, "populateCommonDb");




/***/ }),

/***/ 3782:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(80910);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(42626);
/* harmony import */ var _chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(30103);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86906);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(13197);





// src/diagrams/radar/db.ts
var defaultOptions = {
  showLegend: true,
  ticks: 5,
  max: null,
  min: 0,
  graticule: "circle"
};
var defaultRadarData = {
  axes: [],
  curves: [],
  options: defaultOptions
};
var data = structuredClone(defaultRadarData);
var DEFAULT_RADAR_CONFIG = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .defaultConfig_default */ .vZ.radar;
var getConfig2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => {
  const config = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)({
    ...DEFAULT_RADAR_CONFIG,
    ...(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)().radar
  });
  return config;
}, "getConfig");
var getAxes = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => data.axes, "getAxes");
var getCurves = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => data.curves, "getCurves");
var getOptions = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => data.options, "getOptions");
var setAxes = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((axes) => {
  data.axes = axes.map((axis) => {
    return {
      name: axis.name,
      label: axis.label ?? axis.name
    };
  });
}, "setAxes");
var setCurves = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((curves) => {
  data.curves = curves.map((curve) => {
    return {
      name: curve.name,
      label: curve.label ?? curve.name,
      entries: computeCurveEntries(curve.entries)
    };
  });
}, "setCurves");
var computeCurveEntries = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((entries) => {
  if (entries[0].axis == void 0) {
    return entries.map((entry) => entry.value);
  }
  const axes = getAxes();
  if (axes.length === 0) {
    throw new Error("Axes must be populated before curves for reference entries");
  }
  return axes.map((axis) => {
    const entry = entries.find((entry2) => entry2.axis?.$refText === axis.name);
    if (entry === void 0) {
      throw new Error("Missing entry for axis " + axis.label);
    }
    return entry.value;
  });
}, "computeCurveEntries");
var setOptions = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((options) => {
  const optionMap = options.reduce(
    (acc, option) => {
      acc[option.name] = option;
      return acc;
    },
    {}
  );
  data.options = {
    showLegend: optionMap.showLegend?.value ?? defaultOptions.showLegend,
    ticks: optionMap.ticks?.value ?? defaultOptions.ticks,
    max: optionMap.max?.value ?? defaultOptions.max,
    min: optionMap.min?.value ?? defaultOptions.min,
    graticule: optionMap.graticule?.value ?? defaultOptions.graticule
  };
}, "setOptions");
var clear2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => {
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
  data = structuredClone(defaultRadarData);
}, "clear");
var db = {
  getAxes,
  getCurves,
  getOptions,
  setAxes,
  setCurves,
  setOptions,
  getConfig: getConfig2,
  clear: clear2,
  setAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu,
  setDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr,
  getAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx,
  setAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$
};

// src/diagrams/radar/parser.ts

var populate = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((ast) => {
  (0,_chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .populateCommonDb */ .A)(ast, db);
  const { axes, curves, options } = ast;
  db.setAxes(axes);
  db.setCurves(curves);
  db.setOptions(options);
}, "populate");
var parser = {
  parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__/* .parse */ .Qc)("radar", input);
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(ast);
    populate(ast);
  }, "parse")
};

// src/diagrams/radar/renderer.ts
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((_text, id, _version, diagram2) => {
  const db2 = diagram2.db;
  const axes = db2.getAxes();
  const curves = db2.getCurves();
  const options = db2.getOptions();
  const config = db2.getConfig();
  const title = db2.getDiagramTitle();
  const svg = (0,_chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__/* .selectSvgElement */ .P)(id);
  const g = drawFrame(svg, config);
  const maxValue = options.max ?? Math.max(...curves.map((curve) => Math.max(...curve.entries)));
  const minValue = options.min;
  const radius = Math.min(config.width, config.height) / 2;
  drawGraticule(g, axes, radius, options.ticks, options.graticule);
  drawAxes(g, axes, radius, config);
  drawCurves(g, axes, curves, minValue, maxValue, options.graticule, config);
  drawLegend(g, curves, options.showLegend, config);
  g.append("text").attr("class", "radarTitle").text(title).attr("x", 0).attr("y", -config.height / 2 - config.marginTop);
}, "draw");
var drawFrame = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((svg, config) => {
  const totalWidth = config.width + config.marginLeft + config.marginRight;
  const totalHeight = config.height + config.marginTop + config.marginBottom;
  const center = {
    x: config.marginLeft + config.width / 2,
    y: config.marginTop + config.height / 2
  };
  svg.attr("viewbox", `0 0 ${totalWidth} ${totalHeight}`).attr("width", totalWidth).attr("height", totalHeight);
  return svg.append("g").attr("transform", `translate(${center.x}, ${center.y})`);
}, "drawFrame");
var drawGraticule = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((g, axes, radius, ticks, graticule) => {
  if (graticule === "circle") {
    for (let i = 0; i < ticks; i++) {
      const r = radius * (i + 1) / ticks;
      g.append("circle").attr("r", r).attr("class", "radarGraticule");
    }
  } else if (graticule === "polygon") {
    const numAxes = axes.length;
    for (let i = 0; i < ticks; i++) {
      const r = radius * (i + 1) / ticks;
      const points = axes.map((_, j) => {
        const angle = 2 * j * Math.PI / numAxes - Math.PI / 2;
        const x = r * Math.cos(angle);
        const y = r * Math.sin(angle);
        return `${x},${y}`;
      }).join(" ");
      g.append("polygon").attr("points", points).attr("class", "radarGraticule");
    }
  }
}, "drawGraticule");
var drawAxes = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((g, axes, radius, config) => {
  const numAxes = axes.length;
  for (let i = 0; i < numAxes; i++) {
    const label = axes[i].label;
    const angle = 2 * i * Math.PI / numAxes - Math.PI / 2;
    g.append("line").attr("x1", 0).attr("y1", 0).attr("x2", radius * config.axisScaleFactor * Math.cos(angle)).attr("y2", radius * config.axisScaleFactor * Math.sin(angle)).attr("class", "radarAxisLine");
    g.append("text").text(label).attr("x", radius * config.axisLabelFactor * Math.cos(angle)).attr("y", radius * config.axisLabelFactor * Math.sin(angle)).attr("class", "radarAxisLabel");
  }
}, "drawAxes");
function drawCurves(g, axes, curves, minValue, maxValue, graticule, config) {
  const numAxes = axes.length;
  const radius = Math.min(config.width, config.height) / 2;
  curves.forEach((curve, index) => {
    if (curve.entries.length !== numAxes) {
      return;
    }
    const points = curve.entries.map((entry, i) => {
      const angle = 2 * Math.PI * i / numAxes - Math.PI / 2;
      const r = relativeRadius(entry, minValue, maxValue, radius);
      const x = r * Math.cos(angle);
      const y = r * Math.sin(angle);
      return { x, y };
    });
    if (graticule === "circle") {
      g.append("path").attr("d", closedRoundCurve(points, config.curveTension)).attr("class", `radarCurve-${index}`);
    } else if (graticule === "polygon") {
      g.append("polygon").attr("points", points.map((p) => `${p.x},${p.y}`).join(" ")).attr("class", `radarCurve-${index}`);
    }
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(drawCurves, "drawCurves");
function relativeRadius(value, minValue, maxValue, radius) {
  const clippedValue = Math.min(Math.max(value, minValue), maxValue);
  return radius * (clippedValue - minValue) / (maxValue - minValue);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(relativeRadius, "relativeRadius");
function closedRoundCurve(points, tension) {
  const numPoints = points.length;
  let d = `M${points[0].x},${points[0].y}`;
  for (let i = 0; i < numPoints; i++) {
    const p0 = points[(i - 1 + numPoints) % numPoints];
    const p1 = points[i];
    const p2 = points[(i + 1) % numPoints];
    const p3 = points[(i + 2) % numPoints];
    const cp1 = {
      x: p1.x + (p2.x - p0.x) * tension,
      y: p1.y + (p2.y - p0.y) * tension
    };
    const cp2 = {
      x: p2.x - (p3.x - p1.x) * tension,
      y: p2.y - (p3.y - p1.y) * tension
    };
    d += ` C${cp1.x},${cp1.y} ${cp2.x},${cp2.y} ${p2.x},${p2.y}`;
  }
  return `${d} Z`;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(closedRoundCurve, "closedRoundCurve");
function drawLegend(g, curves, showLegend, config) {
  if (!showLegend) {
    return;
  }
  const legendX = (config.width / 2 + config.marginRight) * 3 / 4;
  const legendY = -(config.height / 2 + config.marginTop) * 3 / 4;
  const lineHeight = 20;
  curves.forEach((curve, index) => {
    const itemGroup = g.append("g").attr("transform", `translate(${legendX}, ${legendY + index * lineHeight})`);
    itemGroup.append("rect").attr("width", 12).attr("height", 12).attr("class", `radarLegendBox-${index}`);
    itemGroup.append("text").attr("x", 16).attr("y", 0).attr("class", "radarLegendText").text(curve.label);
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(drawLegend, "drawLegend");
var renderer = { draw };

// src/diagrams/radar/styles.ts
var genIndexStyles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((themeVariables, radarOptions) => {
  let sections = "";
  for (let i = 0; i < themeVariables.THEME_COLOR_LIMIT; i++) {
    const indexColor = themeVariables[`cScale${i}`];
    sections += `
		.radarCurve-${i} {
			color: ${indexColor};
			fill: ${indexColor};
			fill-opacity: ${radarOptions.curveOpacity};
			stroke: ${indexColor};
			stroke-width: ${radarOptions.curveStrokeWidth};
		}
		.radarLegendBox-${i} {
			fill: ${indexColor};
			fill-opacity: ${radarOptions.curveOpacity};
			stroke: ${indexColor};
		}
		`;
  }
  return sections;
}, "genIndexStyles");
var buildRadarStyleOptions = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((radar) => {
  const defaultThemeVariables = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getThemeVariables */ .xN)();
  const currentConfig = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)();
  const themeVariables = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)(defaultThemeVariables, currentConfig.themeVariables);
  const radarOptions = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)(themeVariables.radar, radar);
  return { themeVariables, radarOptions };
}, "buildRadarStyleOptions");
var styles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(({ radar } = {}) => {
  const { themeVariables, radarOptions } = buildRadarStyleOptions(radar);
  return `
	.radarTitle {
		font-size: ${themeVariables.fontSize};
		color: ${themeVariables.titleColor};
		dominant-baseline: hanging;
		text-anchor: middle;
	}
	.radarAxisLine {
		stroke: ${radarOptions.axisColor};
		stroke-width: ${radarOptions.axisStrokeWidth};
	}
	.radarAxisLabel {
		dominant-baseline: middle;
		text-anchor: middle;
		font-size: ${radarOptions.axisLabelFontSize}px;
		color: ${radarOptions.axisColor};
	}
	.radarGraticule {
		fill: ${radarOptions.graticuleColor};
		fill-opacity: ${radarOptions.graticuleOpacity};
		stroke: ${radarOptions.graticuleColor};
		stroke-width: ${radarOptions.graticuleStrokeWidth};
	}
	.radarLegendText {
		text-anchor: start;
		font-size: ${radarOptions.legendFontSize}px;
		dominant-baseline: hanging;
	}
	${genIndexStyles(themeVariables, radarOptions)}
	`;
}, "styles");

// src/diagrams/radar/diagram.ts
var diagram = {
  parser,
  db,
  renderer,
  styles
};



/***/ })

}]);
//# sourceMappingURL=3782.b5169726474369258b8f.js.map?v=b5169726474369258b8f