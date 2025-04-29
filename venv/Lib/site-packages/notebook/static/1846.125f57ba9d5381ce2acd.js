"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1846],{

/***/ 11846:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _styles_9dd40fb9_js__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(17010);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(23617);
/* harmony import */ var dagre_d3_es_src_graphlib_index_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(67406);
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(24028);
/* harmony import */ var _index_0980fb80_js__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(86281);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(27693);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(7608);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(31699);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(dompurify__WEBPACK_IMPORTED_MODULE_4__);
/* harmony import */ var dagre_d3_es_src_dagre_index_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(7259);
/* harmony import */ var dagre_d3_es_src_graphlib_json_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(81779);



















const sanitizeText = (txt) => _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.e.sanitizeText(txt, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)());
let conf = {
  dividerMargin: 10,
  padding: 5,
  textHeight: 10,
  curve: void 0
};
const addNamespaces = function(namespaces, g, _id, diagObj) {
  const keys = Object.keys(namespaces);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("keys:", keys);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info(namespaces);
  keys.forEach(function(id) {
    var _a, _b;
    const vertex = namespaces[id];
    const shape = "rect";
    const node = {
      shape,
      id: vertex.id,
      domId: vertex.domId,
      labelText: sanitizeText(vertex.id),
      labelStyle: "",
      style: "fill: none; stroke: black",
      // TODO V10: Flowchart ? Keeping flowchart for backwards compatibility. Remove in next major release
      padding: ((_a = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart) == null ? void 0 : _a.padding) ?? ((_b = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().class) == null ? void 0 : _b.padding)
    };
    g.setNode(vertex.id, node);
    addClasses(vertex.classes, g, _id, diagObj, vertex.id);
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("setNode", node);
  });
};
const addClasses = function(classes, g, _id, diagObj, parent) {
  const keys = Object.keys(classes);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("keys:", keys);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info(classes);
  keys.filter((id) => classes[id].parent == parent).forEach(function(id) {
    var _a, _b;
    const vertex = classes[id];
    const cssClassStr = vertex.cssClasses.join(" ");
    const styles2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.k)(vertex.styles);
    const vertexText = vertex.label ?? vertex.id;
    const radius = 0;
    const shape = "class_box";
    const node = {
      labelStyle: styles2.labelStyle,
      shape,
      labelText: sanitizeText(vertexText),
      classData: vertex,
      rx: radius,
      ry: radius,
      class: cssClassStr,
      style: styles2.style,
      id: vertex.id,
      domId: vertex.domId,
      tooltip: diagObj.db.getTooltip(vertex.id, parent) || "",
      haveCallback: vertex.haveCallback,
      link: vertex.link,
      width: vertex.type === "group" ? 500 : void 0,
      type: vertex.type,
      // TODO V10: Flowchart ? Keeping flowchart for backwards compatibility. Remove in next major release
      padding: ((_a = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart) == null ? void 0 : _a.padding) ?? ((_b = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().class) == null ? void 0 : _b.padding)
    };
    g.setNode(vertex.id, node);
    if (parent) {
      g.setParent(vertex.id, parent);
    }
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("setNode", node);
  });
};
const addNotes = function(notes, g, startEdgeId, classes) {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info(notes);
  notes.forEach(function(note, i) {
    var _a, _b;
    const vertex = note;
    const cssNoteStr = "";
    const styles2 = { labelStyle: "", style: "" };
    const vertexText = vertex.text;
    const radius = 0;
    const shape = "note";
    const node = {
      labelStyle: styles2.labelStyle,
      shape,
      labelText: sanitizeText(vertexText),
      noteData: vertex,
      rx: radius,
      ry: radius,
      class: cssNoteStr,
      style: styles2.style,
      id: vertex.id,
      domId: vertex.id,
      tooltip: "",
      type: "note",
      // TODO V10: Flowchart ? Keeping flowchart for backwards compatibility. Remove in next major release
      padding: ((_a = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart) == null ? void 0 : _a.padding) ?? ((_b = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().class) == null ? void 0 : _b.padding)
    };
    g.setNode(vertex.id, node);
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("setNode", node);
    if (!vertex.class || !(vertex.class in classes)) {
      return;
    }
    const edgeId = startEdgeId + i;
    const edgeData = {
      id: `edgeNote${edgeId}`,
      //Set relationship style and line type
      classes: "relation",
      pattern: "dotted",
      // Set link type for rendering
      arrowhead: "none",
      //Set edge extra labels
      startLabelRight: "",
      endLabelLeft: "",
      //Set relation arrow types
      arrowTypeStart: "none",
      arrowTypeEnd: "none",
      style: "fill:none",
      labelStyle: "",
      curve: (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.n)(conf.curve, d3__WEBPACK_IMPORTED_MODULE_0__/* .curveLinear */ .c_6)
    };
    g.setEdge(vertex.id, vertex.class, edgeData, edgeId);
  });
};
const addRelations = function(relations, g) {
  const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart;
  let cnt = 0;
  relations.forEach(function(edge) {
    var _a;
    cnt++;
    const edgeData = {
      //Set relationship style and line type
      classes: "relation",
      pattern: edge.relation.lineType == 1 ? "dashed" : "solid",
      id: `id_${edge.id1}_${edge.id2}_${cnt}`,
      // Set link type for rendering
      arrowhead: edge.type === "arrow_open" ? "none" : "normal",
      //Set edge extra labels
      startLabelRight: edge.relationTitle1 === "none" ? "" : edge.relationTitle1,
      endLabelLeft: edge.relationTitle2 === "none" ? "" : edge.relationTitle2,
      //Set relation arrow types
      arrowTypeStart: getArrowMarker(edge.relation.type1),
      arrowTypeEnd: getArrowMarker(edge.relation.type2),
      style: "fill:none",
      labelStyle: "",
      curve: (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.n)(conf2 == null ? void 0 : conf2.curve, d3__WEBPACK_IMPORTED_MODULE_0__/* .curveLinear */ .c_6)
    };
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info(edgeData, edge);
    if (edge.style !== void 0) {
      const styles2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.k)(edge.style);
      edgeData.style = styles2.style;
      edgeData.labelStyle = styles2.labelStyle;
    }
    edge.text = edge.title;
    if (edge.text === void 0) {
      if (edge.style !== void 0) {
        edgeData.arrowheadStyle = "fill: #333";
      }
    } else {
      edgeData.arrowheadStyle = "fill: #333";
      edgeData.labelpos = "c";
      if (((_a = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart) == null ? void 0 : _a.htmlLabels) ?? (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().htmlLabels) {
        edgeData.labelType = "html";
        edgeData.label = '<span class="edgeLabel">' + edge.text + "</span>";
      } else {
        edgeData.labelType = "text";
        edgeData.label = edge.text.replace(_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.e.lineBreakRegex, "\n");
        if (edge.style === void 0) {
          edgeData.style = edgeData.style || "stroke: #333; stroke-width: 1.5px;fill:none";
        }
        edgeData.labelStyle = edgeData.labelStyle.replace("color:", "fill:");
      }
    }
    g.setEdge(edge.id1, edge.id2, edgeData, cnt);
  });
};
const setConf = function(cnf) {
  conf = {
    ...conf,
    ...cnf
  };
};
const draw = async function(text, id, _version, diagObj) {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("Drawing class - ", id);
  const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().flowchart ?? (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().class;
  const securityLevel = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.c)().securityLevel;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info("config:", conf2);
  const nodeSpacing = (conf2 == null ? void 0 : conf2.nodeSpacing) ?? 50;
  const rankSpacing = (conf2 == null ? void 0 : conf2.rankSpacing) ?? 50;
  const g = new dagre_d3_es_src_graphlib_index_js__WEBPACK_IMPORTED_MODULE_1__/* .Graph */ .k({
    multigraph: true,
    compound: true
  }).setGraph({
    rankdir: diagObj.db.getDirection(),
    nodesep: nodeSpacing,
    ranksep: rankSpacing,
    marginx: 8,
    marginy: 8
  }).setDefaultEdgeLabel(function() {
    return {};
  });
  const namespaces = diagObj.db.getNamespaces();
  const classes = diagObj.db.getClasses();
  const relations = diagObj.db.getRelations();
  const notes = diagObj.db.getNotes();
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.l.info(relations);
  addNamespaces(namespaces, g, id, diagObj);
  addClasses(classes, g, id, diagObj);
  addRelations(relations, g);
  addNotes(notes, g, relations.length + 1, classes);
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("body");
  const svg = root.select(`[id="${id}"]`);
  const element = root.select("#" + id + " g");
  await (0,_index_0980fb80_js__WEBPACK_IMPORTED_MODULE_8__.r)(
    element,
    g,
    ["aggregation", "extension", "composition", "dependency", "lollipop"],
    "classDiagram",
    id
  );
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.u.insertTitle(svg, "classTitleText", (conf2 == null ? void 0 : conf2.titleTopMargin) ?? 5, diagObj.db.getDiagramTitle());
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_7__.o)(g, svg, conf2 == null ? void 0 : conf2.diagramPadding, conf2 == null ? void 0 : conf2.useMaxWidth);
  if (!(conf2 == null ? void 0 : conf2.htmlLabels)) {
    const doc = securityLevel === "sandbox" ? sandboxElement.nodes()[0].contentDocument : document;
    const labels = doc.querySelectorAll('[id="' + id + '"] .edgeLabel .label');
    for (const label of labels) {
      const dim = label.getBBox();
      const rect = doc.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("rx", 0);
      rect.setAttribute("ry", 0);
      rect.setAttribute("width", dim.width);
      rect.setAttribute("height", dim.height);
      label.insertBefore(rect, label.firstChild);
    }
  }
};
function getArrowMarker(type) {
  let marker;
  switch (type) {
    case 0:
      marker = "aggregation";
      break;
    case 1:
      marker = "extension";
      break;
    case 2:
      marker = "composition";
      break;
    case 3:
      marker = "dependency";
      break;
    case 4:
      marker = "lollipop";
      break;
    default:
      marker = "none";
  }
  return marker;
}
const renderer = {
  setConf,
  draw
};
const diagram = {
  parser: _styles_9dd40fb9_js__WEBPACK_IMPORTED_MODULE_9__.p,
  db: _styles_9dd40fb9_js__WEBPACK_IMPORTED_MODULE_9__.d,
  renderer,
  styles: _styles_9dd40fb9_js__WEBPACK_IMPORTED_MODULE_9__.s,
  init: (cnf) => {
    if (!cnf.class) {
      cnf.class = {};
    }
    cnf.class.arrowMarkerAbsolute = cnf.arrowMarkerAbsolute;
    _styles_9dd40fb9_js__WEBPACK_IMPORTED_MODULE_9__.d.clear();
  }
};



/***/ })

}]);
//# sourceMappingURL=1846.125f57ba9d5381ce2acd.js.map?v=125f57ba9d5381ce2acd