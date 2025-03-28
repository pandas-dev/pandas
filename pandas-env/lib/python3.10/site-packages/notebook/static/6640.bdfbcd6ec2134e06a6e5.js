"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[6640],{

/***/ 6640:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(75343);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(23617);
/* harmony import */ var dagre_d3_es_src_dagre_index_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(7259);
/* harmony import */ var dagre_d3_es_src_graphlib_index_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(67406);
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(24028);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(27693);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(7608);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(31699);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(dompurify__WEBPACK_IMPORTED_MODULE_5__);














const idCache = {};
const set = (key, val) => {
  idCache[key] = val;
};
const get = (k) => idCache[k];
const keys = () => Object.keys(idCache);
const size = () => keys().length;
const idCache$1 = {
  get,
  set,
  keys,
  size
};
const drawStartState = (g) => g.append("circle").attr("class", "start-state").attr("r", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit).attr("cx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit).attr("cy", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit);
const drawDivider = (g) => g.append("line").style("stroke", "grey").style("stroke-dasharray", "3").attr("x1", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight).attr("class", "divider").attr("x2", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight * 2).attr("y1", 0).attr("y2", 0);
const drawSimpleState = (g, stateDef) => {
  const state = g.append("text").attr("x", 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("font-size", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.fontSize).attr("class", "state-title").text(stateDef.id);
  const classBox = state.node().getBBox();
  g.insert("rect", ":first-child").attr("x", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("width", classBox.width + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("height", classBox.height + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("rx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.radius);
  return state;
};
const drawDescrState = (g, stateDef) => {
  const addTspan = function(textEl, txt, isFirst2) {
    const tSpan = textEl.append("tspan").attr("x", 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).text(txt);
    if (!isFirst2) {
      tSpan.attr("dy", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight);
    }
  };
  const title = g.append("text").attr("x", 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight + 1.3 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("font-size", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.fontSize).attr("class", "state-title").text(stateDef.descriptions[0]);
  const titleBox = title.node().getBBox();
  const titleHeight = titleBox.height;
  const description = g.append("text").attr("x", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr(
    "y",
    titleHeight + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding * 0.4 + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.dividerMargin + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight
  ).attr("class", "state-description");
  let isFirst = true;
  let isSecond = true;
  stateDef.descriptions.forEach(function(descr) {
    if (!isFirst) {
      addTspan(description, descr, isSecond);
      isSecond = false;
    }
    isFirst = false;
  });
  const descrLine = g.append("line").attr("x1", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y1", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + titleHeight + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.dividerMargin / 2).attr("y2", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + titleHeight + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.dividerMargin / 2).attr("class", "descr-divider");
  const descrBox = description.node().getBBox();
  const width = Math.max(descrBox.width, titleBox.width);
  descrLine.attr("x2", width + 3 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding);
  g.insert("rect", ":first-child").attr("x", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("width", width + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("height", descrBox.height + titleHeight + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("rx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.radius);
  return g;
};
const addTitleAndBox = (g, stateDef, altBkg) => {
  const pad = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding;
  const dblPad = 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding;
  const orgBox = g.node().getBBox();
  const orgWidth = orgBox.width;
  const orgX = orgBox.x;
  const title = g.append("text").attr("x", 0).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.titleShift).attr("font-size", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.fontSize).attr("class", "state-title").text(stateDef.id);
  const titleBox = title.node().getBBox();
  const titleWidth = titleBox.width + dblPad;
  let width = Math.max(titleWidth, orgWidth);
  if (width === orgWidth) {
    width = width + dblPad;
  }
  let startX;
  const graphBox = g.node().getBBox();
  if (stateDef.doc)
    ;
  startX = orgX - pad;
  if (titleWidth > orgWidth) {
    startX = (orgWidth - width) / 2 + pad;
  }
  if (Math.abs(orgX - graphBox.x) < pad && titleWidth > orgWidth) {
    startX = orgX - (titleWidth - orgWidth) / 2;
  }
  const lineY = 1 - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight;
  g.insert("rect", ":first-child").attr("x", startX).attr("y", lineY).attr("class", altBkg ? "alt-composit" : "composit").attr("width", width).attr(
    "height",
    graphBox.height + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.titleShift + 1
  ).attr("rx", "0");
  title.attr("x", startX + pad);
  if (titleWidth <= orgWidth) {
    title.attr("x", orgX + (width - dblPad) / 2 - titleWidth / 2 + pad);
  }
  g.insert("rect", ":first-child").attr("x", startX).attr(
    "y",
    (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.titleShift - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding
  ).attr("width", width).attr("height", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight * 3).attr("rx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.radius);
  g.insert("rect", ":first-child").attr("x", startX).attr(
    "y",
    (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.titleShift - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding
  ).attr("width", width).attr("height", graphBox.height + 3 + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.textHeight).attr("rx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.radius);
  return g;
};
const drawEndState = (g) => {
  g.append("circle").attr("class", "end-state-outer").attr("r", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.miniPadding).attr(
    "cx",
    (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.miniPadding
  ).attr(
    "cy",
    (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.miniPadding
  );
  return g.append("circle").attr("class", "end-state-inner").attr("r", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit).attr("cx", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit + 2).attr("cy", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.sizeUnit + 2);
};
const drawForkJoinState = (g, stateDef) => {
  let width = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.forkWidth;
  let height = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.forkHeight;
  if (stateDef.parentId) {
    let tmp = width;
    width = height;
    height = tmp;
  }
  return g.append("rect").style("stroke", "black").style("fill", "black").attr("width", width).attr("height", height).attr("x", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding);
};
const _drawLongText = (_text, x, y, g) => {
  let textHeight = 0;
  const textElem = g.append("text");
  textElem.style("text-anchor", "start");
  textElem.attr("class", "noteText");
  let text = _text.replace(/\r\n/g, "<br/>");
  text = text.replace(/\n/g, "<br/>");
  const lines = text.split(_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.e.lineBreakRegex);
  let tHeight = 1.25 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.noteMargin;
  for (const line2 of lines) {
    const txt = line2.trim();
    if (txt.length > 0) {
      const span = textElem.append("tspan");
      span.text(txt);
      if (tHeight === 0) {
        const textBounds = span.node().getBBox();
        tHeight += textBounds.height;
      }
      textHeight += tHeight;
      span.attr("x", x + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.noteMargin);
      span.attr("y", y + textHeight + 1.25 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.noteMargin);
    }
  }
  return { textWidth: textElem.node().getBBox().width, textHeight };
};
const drawNote = (text, g) => {
  g.attr("class", "state-note");
  const note = g.append("rect").attr("x", 0).attr("y", (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding);
  const rectElem = g.append("g");
  const { textWidth, textHeight } = _drawLongText(text, 0, 0, rectElem);
  note.attr("height", textHeight + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.noteMargin);
  note.attr("width", textWidth + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.noteMargin * 2);
  return note;
};
const drawState = function(elem, stateDef) {
  const id = stateDef.id;
  const stateInfo = {
    id,
    label: stateDef.id,
    width: 0,
    height: 0
  };
  const g = elem.append("g").attr("id", id).attr("class", "stateGroup");
  if (stateDef.type === "start") {
    drawStartState(g);
  }
  if (stateDef.type === "end") {
    drawEndState(g);
  }
  if (stateDef.type === "fork" || stateDef.type === "join") {
    drawForkJoinState(g, stateDef);
  }
  if (stateDef.type === "note") {
    drawNote(stateDef.note.text, g);
  }
  if (stateDef.type === "divider") {
    drawDivider(g);
  }
  if (stateDef.type === "default" && stateDef.descriptions.length === 0) {
    drawSimpleState(g, stateDef);
  }
  if (stateDef.type === "default" && stateDef.descriptions.length > 0) {
    drawDescrState(g, stateDef);
  }
  const stateBox = g.node().getBBox();
  stateInfo.width = stateBox.width + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding;
  stateInfo.height = stateBox.height + 2 * (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding;
  idCache$1.set(id, stateInfo);
  return stateInfo;
};
let edgeCount = 0;
const drawEdge = function(elem, path, relation) {
  const getRelationType = function(type) {
    switch (type) {
      case _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.relationType.AGGREGATION:
        return "aggregation";
      case _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.relationType.EXTENSION:
        return "extension";
      case _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.relationType.COMPOSITION:
        return "composition";
      case _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.relationType.DEPENDENCY:
        return "dependency";
    }
  };
  path.points = path.points.filter((p) => !Number.isNaN(p.y));
  const lineData = path.points;
  const lineFunction = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .line */ .jvg)().x(function(d) {
    return d.x;
  }).y(function(d) {
    return d.y;
  }).curve(d3__WEBPACK_IMPORTED_MODULE_0__/* .curveBasis */ .$0Z);
  const svgPath = elem.append("path").attr("d", lineFunction(lineData)).attr("id", "edge" + edgeCount).attr("class", "transition");
  let url = "";
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.arrowMarkerAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = url.replace(/\(/g, "\\(");
    url = url.replace(/\)/g, "\\)");
  }
  svgPath.attr(
    "marker-end",
    "url(" + url + "#" + getRelationType(_styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.relationType.DEPENDENCY) + "End)"
  );
  if (relation.title !== void 0) {
    const label = elem.append("g").attr("class", "stateLabel");
    const { x, y } = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.u.calcLabelPosition(path.points);
    const rows = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.e.getRows(relation.title);
    let titleHeight = 0;
    const titleRows = [];
    let maxWidth = 0;
    let minX = 0;
    for (let i = 0; i <= rows.length; i++) {
      const title = label.append("text").attr("text-anchor", "middle").text(rows[i]).attr("x", x).attr("y", y + titleHeight);
      const boundstmp = title.node().getBBox();
      maxWidth = Math.max(maxWidth, boundstmp.width);
      minX = Math.min(minX, boundstmp.x);
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.info(boundstmp.x, x, y + titleHeight);
      if (titleHeight === 0) {
        const titleBox = title.node().getBBox();
        titleHeight = titleBox.height;
        _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.info("Title height", titleHeight, y);
      }
      titleRows.push(title);
    }
    let boxHeight = titleHeight * rows.length;
    if (rows.length > 1) {
      const heightAdj = (rows.length - 1) * titleHeight * 0.5;
      titleRows.forEach((title, i) => title.attr("y", y + i * titleHeight - heightAdj));
      boxHeight = titleHeight * rows.length;
    }
    const bounds = label.node().getBBox();
    label.insert("rect", ":first-child").attr("class", "box").attr("x", x - maxWidth / 2 - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding / 2).attr("y", y - boxHeight / 2 - (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding / 2 - 3.5).attr("width", maxWidth + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding).attr("height", boxHeight + (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state.padding);
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.info(bounds);
  }
  edgeCount++;
};
let conf;
const transformationLog = {};
const setConf = function() {
};
const insertMarkers = function(elem) {
  elem.append("defs").append("marker").attr("id", "dependencyEnd").attr("refX", 19).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 19,7 L9,13 L14,7 L9,1 Z");
};
const draw = function(text, id, _version, diagObj) {
  conf = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().state;
  const securityLevel = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.c)().securityLevel;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("body");
  const doc = securityLevel === "sandbox" ? sandboxElement.nodes()[0].contentDocument : document;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Rendering diagram " + text);
  const diagram2 = root.select(`[id='${id}']`);
  insertMarkers(diagram2);
  const rootDoc = diagObj.db.getRootDoc();
  renderDoc(rootDoc, diagram2, void 0, false, root, doc, diagObj);
  const padding = conf.padding;
  const bounds = diagram2.node().getBBox();
  const width = bounds.width + padding * 2;
  const height = bounds.height + padding * 2;
  const svgWidth = width * 1.75;
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.i)(diagram2, height, svgWidth, conf.useMaxWidth);
  diagram2.attr(
    "viewBox",
    `${bounds.x - conf.padding}  ${bounds.y - conf.padding} ` + width + " " + height
  );
};
const getLabelWidth = (text) => {
  return text ? text.length * conf.fontSizeFactor : 1;
};
const renderDoc = (doc, diagram2, parentId, altBkg, root, domDocument, diagObj) => {
  const graph = new dagre_d3_es_src_graphlib_index_js__WEBPACK_IMPORTED_MODULE_2__/* .Graph */ .k({
    compound: true,
    multigraph: true
  });
  let i;
  let edgeFreeDoc = true;
  for (i = 0; i < doc.length; i++) {
    if (doc[i].stmt === "relation") {
      edgeFreeDoc = false;
      break;
    }
  }
  if (parentId) {
    graph.setGraph({
      rankdir: "LR",
      multigraph: true,
      compound: true,
      // acyclicer: 'greedy',
      ranker: "tight-tree",
      ranksep: edgeFreeDoc ? 1 : conf.edgeLengthFactor,
      nodeSep: edgeFreeDoc ? 1 : 50,
      isMultiGraph: true
      // ranksep: 5,
      // nodesep: 1
    });
  } else {
    graph.setGraph({
      rankdir: "TB",
      multigraph: true,
      compound: true,
      // isCompound: true,
      // acyclicer: 'greedy',
      // ranker: 'longest-path'
      ranksep: edgeFreeDoc ? 1 : conf.edgeLengthFactor,
      nodeSep: edgeFreeDoc ? 1 : 50,
      ranker: "tight-tree",
      // ranker: 'network-simplex'
      isMultiGraph: true
    });
  }
  graph.setDefaultEdgeLabel(function() {
    return {};
  });
  diagObj.db.extract(doc);
  const states = diagObj.db.getStates();
  const relations = diagObj.db.getRelations();
  const keys2 = Object.keys(states);
  for (const key of keys2) {
    const stateDef = states[key];
    if (parentId) {
      stateDef.parentId = parentId;
    }
    let node;
    if (stateDef.doc) {
      let sub = diagram2.append("g").attr("id", stateDef.id).attr("class", "stateGroup");
      node = renderDoc(stateDef.doc, sub, stateDef.id, !altBkg, root, domDocument, diagObj);
      {
        sub = addTitleAndBox(sub, stateDef, altBkg);
        let boxBounds = sub.node().getBBox();
        node.width = boxBounds.width;
        node.height = boxBounds.height + conf.padding / 2;
        transformationLog[stateDef.id] = { y: conf.compositTitleSize };
      }
    } else {
      node = drawState(diagram2, stateDef);
    }
    if (stateDef.note) {
      const noteDef = {
        descriptions: [],
        id: stateDef.id + "-note",
        note: stateDef.note,
        type: "note"
      };
      const note = drawState(diagram2, noteDef);
      if (stateDef.note.position === "left of") {
        graph.setNode(node.id + "-note", note);
        graph.setNode(node.id, node);
      } else {
        graph.setNode(node.id, node);
        graph.setNode(node.id + "-note", note);
      }
      graph.setParent(node.id, node.id + "-group");
      graph.setParent(node.id + "-note", node.id + "-group");
    } else {
      graph.setNode(node.id, node);
    }
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Count=", graph.nodeCount(), graph);
  let cnt = 0;
  relations.forEach(function(relation) {
    cnt++;
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Setting edge", relation);
    graph.setEdge(
      relation.id1,
      relation.id2,
      {
        relation,
        width: getLabelWidth(relation.title),
        height: conf.labelHeight * _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.e.getRows(relation.title).length,
        labelpos: "c"
      },
      "id" + cnt
    );
  });
  (0,dagre_d3_es_src_dagre_index_js__WEBPACK_IMPORTED_MODULE_1__/* .layout */ .bK)(graph);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Graph after layout", graph.nodes());
  const svgElem = diagram2.node();
  graph.nodes().forEach(function(v) {
    if (v !== void 0 && graph.node(v) !== void 0) {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.warn("Node " + v + ": " + JSON.stringify(graph.node(v)));
      root.select("#" + svgElem.id + " #" + v).attr(
        "transform",
        "translate(" + (graph.node(v).x - graph.node(v).width / 2) + "," + (graph.node(v).y + (transformationLog[v] ? transformationLog[v].y : 0) - graph.node(v).height / 2) + " )"
      );
      root.select("#" + svgElem.id + " #" + v).attr("data-x-shift", graph.node(v).x - graph.node(v).width / 2);
      const dividers = domDocument.querySelectorAll("#" + svgElem.id + " #" + v + " .divider");
      dividers.forEach((divider) => {
        const parent = divider.parentElement;
        let pWidth = 0;
        let pShift = 0;
        if (parent) {
          if (parent.parentElement) {
            pWidth = parent.parentElement.getBBox().width;
          }
          pShift = parseInt(parent.getAttribute("data-x-shift"), 10);
          if (Number.isNaN(pShift)) {
            pShift = 0;
          }
        }
        divider.setAttribute("x1", 0 - pShift + 8);
        divider.setAttribute("x2", pWidth - pShift - 8);
      });
    } else {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("No Node " + v + ": " + JSON.stringify(graph.node(v)));
    }
  });
  let stateBox = svgElem.getBBox();
  graph.edges().forEach(function(e) {
    if (e !== void 0 && graph.edge(e) !== void 0) {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Edge " + e.v + " -> " + e.w + ": " + JSON.stringify(graph.edge(e)));
      drawEdge(diagram2, graph.edge(e), graph.edge(e).relation);
    }
  });
  stateBox = svgElem.getBBox();
  const stateInfo = {
    id: parentId ? parentId : "root",
    label: parentId ? parentId : "root",
    width: 0,
    height: 0
  };
  stateInfo.width = stateBox.width + 2 * conf.padding;
  stateInfo.height = stateBox.height + 2 * conf.padding;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_6__.l.debug("Doc rendered", stateInfo, graph);
  return stateInfo;
};
const renderer = {
  setConf,
  draw
};
const diagram = {
  parser: _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.p,
  db: _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d,
  renderer,
  styles: _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.s,
  init: (cnf) => {
    if (!cnf.state) {
      cnf.state = {};
    }
    cnf.state.arrowMarkerAbsolute = cnf.arrowMarkerAbsolute;
    _styles_401768ef_js__WEBPACK_IMPORTED_MODULE_7__.d.clear();
  }
};



/***/ })

}]);
//# sourceMappingURL=6640.bdfbcd6ec2134e06a6e5.js.map?v=bdfbcd6ec2134e06a6e5