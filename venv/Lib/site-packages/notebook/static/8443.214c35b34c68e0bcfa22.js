"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8443],{

/***/ 98443:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   a: () => (/* binding */ insertMarkers$1),
/* harmony export */   b: () => (/* binding */ clear$1),
/* harmony export */   c: () => (/* binding */ createLabel$1),
/* harmony export */   d: () => (/* binding */ clear),
/* harmony export */   e: () => (/* binding */ insertNode),
/* harmony export */   f: () => (/* binding */ insertEdgeLabel),
/* harmony export */   g: () => (/* binding */ insertEdge),
/* harmony export */   h: () => (/* binding */ positionEdgeLabel),
/* harmony export */   i: () => (/* binding */ intersectRect$1),
/* harmony export */   j: () => (/* binding */ getLineFunctionsWithOffset),
/* harmony export */   k: () => (/* binding */ addEdgeMarkers),
/* harmony export */   l: () => (/* binding */ labelHelper),
/* harmony export */   p: () => (/* binding */ positionNode),
/* harmony export */   s: () => (/* binding */ setNodeElem),
/* harmony export */   u: () => (/* binding */ updateNodeBounds)
/* harmony export */ });
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(24028);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(23617);
/* harmony import */ var _createText_a48a4c90_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(14780);



const insertMarkers = (elem, markerArray, type, id) => {
  markerArray.forEach((markerName) => {
    markers[markerName](elem, type, id);
  });
};
const extension = (elem, type, id) => {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.trace("Making markers for ", id);
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-extensionStart").attr("class", "marker extension " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 1,7 L18,13 V 1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-extensionEnd").attr("class", "marker extension " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 1,1 V 13 L18,7 Z");
};
const composition = (elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-compositionStart").attr("class", "marker composition " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-compositionEnd").attr("class", "marker composition " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
};
const aggregation = (elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-aggregationStart").attr("class", "marker aggregation " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-aggregationEnd").attr("class", "marker aggregation " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
};
const dependency = (elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-dependencyStart").attr("class", "marker dependency " + type).attr("refX", 6).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 5,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-dependencyEnd").attr("class", "marker dependency " + type).attr("refX", 13).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L14,7 L9,1 Z");
};
const lollipop = (elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-lollipopStart").attr("class", "marker lollipop " + type).attr("refX", 13).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("circle").attr("stroke", "black").attr("fill", "transparent").attr("cx", 7).attr("cy", 7).attr("r", 6);
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-lollipopEnd").attr("class", "marker lollipop " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("circle").attr("stroke", "black").attr("fill", "transparent").attr("cx", 7).attr("cy", 7).attr("r", 6);
};
const point = (elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-pointEnd").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 6).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 12).attr("markerHeight", 12).attr("orient", "auto").append("path").attr("d", "M 0 0 L 10 5 L 0 10 z").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-pointStart").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 4.5).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 12).attr("markerHeight", 12).attr("orient", "auto").append("path").attr("d", "M 0 5 L 10 10 L 10 0 z").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
};
const circle$1 = (elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-circleEnd").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 11).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("circle").attr("cx", "5").attr("cy", "5").attr("r", "5").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-circleStart").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", -1).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("circle").attr("cx", "5").attr("cy", "5").attr("r", "5").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
};
const cross = (elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-crossEnd").attr("class", "marker cross " + type).attr("viewBox", "0 0 11 11").attr("refX", 12).attr("refY", 5.2).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("path").attr("d", "M 1,1 l 9,9 M 10,1 l -9,9").attr("class", "arrowMarkerPath").style("stroke-width", 2).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-crossStart").attr("class", "marker cross " + type).attr("viewBox", "0 0 11 11").attr("refX", -1).attr("refY", 5.2).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("path").attr("d", "M 1,1 l 9,9 M 10,1 l -9,9").attr("class", "arrowMarkerPath").style("stroke-width", 2).style("stroke-dasharray", "1,0");
};
const barb = (elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-barbEnd").attr("refX", 19).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 14).attr("markerUnits", "strokeWidth").attr("orient", "auto").append("path").attr("d", "M 19,7 L9,13 L14,7 L9,1 Z");
};
const markers = {
  extension,
  composition,
  aggregation,
  dependency,
  lollipop,
  point,
  circle: circle$1,
  cross,
  barb
};
const insertMarkers$1 = insertMarkers;
function applyStyle(dom, styleFn) {
  if (styleFn) {
    dom.attr("style", styleFn);
  }
}
function addHtmlLabel(node) {
  const fo = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(document.createElementNS("http://www.w3.org/2000/svg", "foreignObject"));
  const div = fo.append("xhtml:div");
  const label = node.label;
  const labelClass = node.isNode ? "nodeLabel" : "edgeLabel";
  div.html(
    '<span class="' + labelClass + '" ' + (node.labelStyle ? 'style="' + node.labelStyle + '"' : "") + ">" + label + "</span>"
  );
  applyStyle(div, node.labelStyle);
  div.style("display", "inline-block");
  div.style("white-space", "nowrap");
  div.attr("xmlns", "http://www.w3.org/1999/xhtml");
  return fo.node();
}
const createLabel = (_vertexText, style, isTitle, isNode) => {
  let vertexText = _vertexText || "";
  if (typeof vertexText === "object") {
    vertexText = vertexText[0];
  }
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    vertexText = vertexText.replace(/\\n|\n/g, "<br />");
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("vertexText" + vertexText);
    const node = {
      isNode,
      label: (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.J)(vertexText).replace(
        /fa[blrs]?:fa-[\w-]+/g,
        (s) => `<i class='${s.replace(":", " ")}'></i>`
      ),
      labelStyle: style.replace("fill:", "color:")
    };
    let vertexNode = addHtmlLabel(node);
    return vertexNode;
  } else {
    const svgLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
    svgLabel.setAttribute("style", style.replace("color:", "fill:"));
    let rows = [];
    if (typeof vertexText === "string") {
      rows = vertexText.split(/\\n|\n|<br\s*\/?>/gi);
    } else if (Array.isArray(vertexText)) {
      rows = vertexText;
    } else {
      rows = [];
    }
    for (const row of rows) {
      const tspan = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
      tspan.setAttributeNS("http://www.w3.org/XML/1998/namespace", "xml:space", "preserve");
      tspan.setAttribute("dy", "1em");
      tspan.setAttribute("x", "0");
      if (isTitle) {
        tspan.setAttribute("class", "title-row");
      } else {
        tspan.setAttribute("class", "row");
      }
      tspan.textContent = row.trim();
      svgLabel.appendChild(tspan);
    }
    return svgLabel;
  }
};
const createLabel$1 = createLabel;
const labelHelper = async (parent, node, _classes, isNode) => {
  let classes;
  const useHtmlLabels = node.useHtmlLabels || (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels);
  if (!_classes) {
    classes = "node default";
  } else {
    classes = _classes;
  }
  const shapeSvg = parent.insert("g").attr("class", classes).attr("id", node.domId || node.id);
  const label = shapeSvg.insert("g").attr("class", "label").attr("style", node.labelStyle);
  let labelText;
  if (node.labelText === void 0) {
    labelText = "";
  } else {
    labelText = typeof node.labelText === "string" ? node.labelText : node.labelText[0];
  }
  const textNode = label.node();
  let text;
  if (node.labelType === "markdown") {
    text = (0,_createText_a48a4c90_js__WEBPACK_IMPORTED_MODULE_2__.a)(label, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.d)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.J)(labelText), (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)()), {
      useHtmlLabels,
      width: node.width || (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.wrappingWidth,
      classes: "markdown-node-label"
    });
  } else {
    text = textNode.appendChild(
      createLabel$1(
        (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.d)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.J)(labelText), (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)()),
        node.labelStyle,
        false,
        isNode
      )
    );
  }
  let bbox = text.getBBox();
  const halfPadding = node.padding / 2;
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    const div = text.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(text);
    const images = div.getElementsByTagName("img");
    if (images) {
      const noImgText = labelText.replace(/<img[^>]*>/g, "").trim() === "";
      await Promise.all(
        [...images].map(
          (img) => new Promise((res) => {
            function setupImage() {
              img.style.display = "flex";
              img.style.flexDirection = "column";
              if (noImgText) {
                const bodyFontSize = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().fontSize ? (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().fontSize : window.getComputedStyle(document.body).fontSize;
                const enlargingFactor = 5;
                const width = parseInt(bodyFontSize, 10) * enlargingFactor + "px";
                img.style.minWidth = width;
                img.style.maxWidth = width;
              } else {
                img.style.width = "100%";
              }
              res(img);
            }
            setTimeout(() => {
              if (img.complete) {
                setupImage();
              }
            });
            img.addEventListener("error", setupImage);
            img.addEventListener("load", setupImage);
          })
        )
      );
    }
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  if (useHtmlLabels) {
    label.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  } else {
    label.attr("transform", "translate(0, " + -bbox.height / 2 + ")");
  }
  if (node.centerLabel) {
    label.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  }
  label.insert("rect", ":first-child");
  return { shapeSvg, bbox, halfPadding, label };
};
const updateNodeBounds = (node, element) => {
  const bbox = element.node().getBBox();
  node.width = bbox.width;
  node.height = bbox.height;
};
function insertPolygonShape(parent, w, h, points) {
  return parent.insert("polygon", ":first-child").attr(
    "points",
    points.map(function(d) {
      return d.x + "," + d.y;
    }).join(" ")
  ).attr("class", "label-container").attr("transform", "translate(" + -w / 2 + "," + h / 2 + ")");
}
function intersectNode(node, point2) {
  return node.intersect(point2);
}
function intersectEllipse(node, rx, ry, point2) {
  var cx = node.x;
  var cy = node.y;
  var px = cx - point2.x;
  var py = cy - point2.y;
  var det = Math.sqrt(rx * rx * py * py + ry * ry * px * px);
  var dx = Math.abs(rx * ry * px / det);
  if (point2.x < cx) {
    dx = -dx;
  }
  var dy = Math.abs(rx * ry * py / det);
  if (point2.y < cy) {
    dy = -dy;
  }
  return { x: cx + dx, y: cy + dy };
}
function intersectCircle(node, rx, point2) {
  return intersectEllipse(node, rx, rx, point2);
}
function intersectLine(p1, p2, q1, q2) {
  var a1, a2, b1, b2, c1, c2;
  var r1, r2, r3, r4;
  var denom, offset, num;
  var x, y;
  a1 = p2.y - p1.y;
  b1 = p1.x - p2.x;
  c1 = p2.x * p1.y - p1.x * p2.y;
  r3 = a1 * q1.x + b1 * q1.y + c1;
  r4 = a1 * q2.x + b1 * q2.y + c1;
  if (r3 !== 0 && r4 !== 0 && sameSign(r3, r4)) {
    return;
  }
  a2 = q2.y - q1.y;
  b2 = q1.x - q2.x;
  c2 = q2.x * q1.y - q1.x * q2.y;
  r1 = a2 * p1.x + b2 * p1.y + c2;
  r2 = a2 * p2.x + b2 * p2.y + c2;
  if (r1 !== 0 && r2 !== 0 && sameSign(r1, r2)) {
    return;
  }
  denom = a1 * b2 - a2 * b1;
  if (denom === 0) {
    return;
  }
  offset = Math.abs(denom / 2);
  num = b1 * c2 - b2 * c1;
  x = num < 0 ? (num - offset) / denom : (num + offset) / denom;
  num = a2 * c1 - a1 * c2;
  y = num < 0 ? (num - offset) / denom : (num + offset) / denom;
  return { x, y };
}
function sameSign(r1, r2) {
  return r1 * r2 > 0;
}
function intersectPolygon(node, polyPoints, point2) {
  var x1 = node.x;
  var y1 = node.y;
  var intersections = [];
  var minX = Number.POSITIVE_INFINITY;
  var minY = Number.POSITIVE_INFINITY;
  if (typeof polyPoints.forEach === "function") {
    polyPoints.forEach(function(entry) {
      minX = Math.min(minX, entry.x);
      minY = Math.min(minY, entry.y);
    });
  } else {
    minX = Math.min(minX, polyPoints.x);
    minY = Math.min(minY, polyPoints.y);
  }
  var left = x1 - node.width / 2 - minX;
  var top = y1 - node.height / 2 - minY;
  for (var i = 0; i < polyPoints.length; i++) {
    var p1 = polyPoints[i];
    var p2 = polyPoints[i < polyPoints.length - 1 ? i + 1 : 0];
    var intersect2 = intersectLine(
      node,
      point2,
      { x: left + p1.x, y: top + p1.y },
      { x: left + p2.x, y: top + p2.y }
    );
    if (intersect2) {
      intersections.push(intersect2);
    }
  }
  if (!intersections.length) {
    return node;
  }
  if (intersections.length > 1) {
    intersections.sort(function(p, q) {
      var pdx = p.x - point2.x;
      var pdy = p.y - point2.y;
      var distp = Math.sqrt(pdx * pdx + pdy * pdy);
      var qdx = q.x - point2.x;
      var qdy = q.y - point2.y;
      var distq = Math.sqrt(qdx * qdx + qdy * qdy);
      return distp < distq ? -1 : distp === distq ? 0 : 1;
    });
  }
  return intersections[0];
}
const intersectRect = (node, point2) => {
  var x = node.x;
  var y = node.y;
  var dx = point2.x - x;
  var dy = point2.y - y;
  var w = node.width / 2;
  var h = node.height / 2;
  var sx, sy;
  if (Math.abs(dy) * w > Math.abs(dx) * h) {
    if (dy < 0) {
      h = -h;
    }
    sx = dy === 0 ? 0 : h * dx / dy;
    sy = h;
  } else {
    if (dx < 0) {
      w = -w;
    }
    sx = w;
    sy = dx === 0 ? 0 : w * dy / dx;
  }
  return { x: x + sx, y: y + sy };
};
const intersectRect$1 = intersectRect;
const intersect = {
  node: intersectNode,
  circle: intersectCircle,
  ellipse: intersectEllipse,
  polygon: intersectPolygon,
  rect: intersectRect$1
};
const note = async (parent, node) => {
  const useHtmlLabels = node.useHtmlLabels || (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels;
  if (!useHtmlLabels) {
    node.centerLabel = true;
  }
  const { shapeSvg, bbox, halfPadding } = await labelHelper(
    parent,
    node,
    "node " + node.classes,
    true
  );
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Classes = ", node.classes);
  const rect2 = shapeSvg.insert("rect", ":first-child");
  rect2.attr("rx", node.rx).attr("ry", node.ry).attr("x", -bbox.width / 2 - halfPadding).attr("y", -bbox.height / 2 - halfPadding).attr("width", bbox.width + node.padding).attr("height", bbox.height + node.padding);
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const note$1 = note;
const formatClass = (str) => {
  if (str) {
    return " " + str;
  }
  return "";
};
const getClassesFromNode = (node, otherClasses) => {
  return `${otherClasses ? otherClasses : "node default"}${formatClass(node.classes)} ${formatClass(
    node.class
  )}`;
};
const question = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const s = w + h;
  const points = [
    { x: s / 2, y: 0 },
    { x: s, y: -s / 2 },
    { x: s / 2, y: -s },
    { x: 0, y: -s / 2 }
  ];
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Question main (Circle)");
  const questionElem = insertPolygonShape(shapeSvg, s, s, points);
  questionElem.attr("style", node.style);
  updateNodeBounds(node, questionElem);
  node.intersect = function(point2) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("Intersect called");
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const choice = (parent, node) => {
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  const s = 28;
  const points = [
    { x: 0, y: s / 2 },
    { x: s / 2, y: 0 },
    { x: 0, y: -s / 2 },
    { x: -s / 2, y: 0 }
  ];
  const choice2 = shapeSvg.insert("polygon", ":first-child").attr(
    "points",
    points.map(function(d) {
      return d.x + "," + d.y;
    }).join(" ")
  );
  choice2.attr("class", "state-start").attr("r", 7).attr("width", 28).attr("height", 28);
  node.width = 28;
  node.height = 28;
  node.intersect = function(point2) {
    return intersect.circle(node, 14, point2);
  };
  return shapeSvg;
};
const hexagon = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const f = 4;
  const h = bbox.height + node.padding;
  const m = h / f;
  const w = bbox.width + 2 * m + node.padding;
  const points = [
    { x: m, y: 0 },
    { x: w - m, y: 0 },
    { x: w, y: -h / 2 },
    { x: w - m, y: -h },
    { x: m, y: -h },
    { x: 0, y: -h / 2 }
  ];
  const hex = insertPolygonShape(shapeSvg, w, h, points);
  hex.attr("style", node.style);
  updateNodeBounds(node, hex);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const rect_left_inv_arrow = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: -h / 2, y: 0 },
    { x: w, y: 0 },
    { x: w, y: -h },
    { x: -h / 2, y: -h },
    { x: 0, y: -h / 2 }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  node.width = w + h;
  node.height = h;
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const lean_right = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(parent, node, getClassesFromNode(node), true);
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: -2 * h / 6, y: 0 },
    { x: w - h / 6, y: 0 },
    { x: w + 2 * h / 6, y: -h },
    { x: h / 6, y: -h }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const lean_left = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: 2 * h / 6, y: 0 },
    { x: w + h / 6, y: 0 },
    { x: w - 2 * h / 6, y: -h },
    { x: -h / 6, y: -h }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const trapezoid = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: -2 * h / 6, y: 0 },
    { x: w + 2 * h / 6, y: 0 },
    { x: w - h / 6, y: -h },
    { x: h / 6, y: -h }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const inv_trapezoid = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: h / 6, y: 0 },
    { x: w - h / 6, y: 0 },
    { x: w + 2 * h / 6, y: -h },
    { x: -2 * h / 6, y: -h }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const rect_right_inv_arrow = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: 0, y: 0 },
    { x: w + h / 2, y: 0 },
    { x: w, y: -h / 2 },
    { x: w + h / 2, y: -h },
    { x: 0, y: -h }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const cylinder = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const rx = w / 2;
  const ry = rx / (2.5 + w / 50);
  const h = bbox.height + ry + node.padding;
  const shape = "M 0," + ry + " a " + rx + "," + ry + " 0,0,0 " + w + " 0 a " + rx + "," + ry + " 0,0,0 " + -w + " 0 l 0," + h + " a " + rx + "," + ry + " 0,0,0 " + w + " 0 l 0," + -h;
  const el = shapeSvg.attr("label-offset-y", ry).insert("path", ":first-child").attr("style", node.style).attr("d", shape).attr("transform", "translate(" + -w / 2 + "," + -(h / 2 + ry) + ")");
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    const pos = intersect.rect(node, point2);
    const x = pos.x - node.x;
    if (rx != 0 && (Math.abs(x) < node.width / 2 || Math.abs(x) == node.width / 2 && Math.abs(pos.y - node.y) > node.height / 2 - ry)) {
      let y = ry * ry * (1 - x * x / (rx * rx));
      if (y != 0) {
        y = Math.sqrt(y);
      }
      y = ry - y;
      if (point2.y - node.y > 0) {
        y = -y;
      }
      pos.y += y;
    }
    return pos;
  };
  return shapeSvg;
};
const rect = async (parent, node) => {
  const { shapeSvg, bbox, halfPadding } = await labelHelper(
    parent,
    node,
    "node " + node.classes + " " + node.class,
    true
  );
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const totalWidth = bbox.width + node.padding;
  const totalHeight = bbox.height + node.padding;
  rect2.attr("class", "basic label-container").attr("style", node.style).attr("rx", node.rx).attr("ry", node.ry).attr("x", -bbox.width / 2 - halfPadding).attr("y", -bbox.height / 2 - halfPadding).attr("width", totalWidth).attr("height", totalHeight);
  if (node.props) {
    const propKeys = new Set(Object.keys(node.props));
    if (node.props.borders) {
      applyNodePropertyBorders(rect2, node.props.borders, totalWidth, totalHeight);
      propKeys.delete("borders");
    }
    propKeys.forEach((propKey) => {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`Unknown node property ${propKey}`);
    });
  }
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const labelRect = async (parent, node) => {
  const { shapeSvg } = await labelHelper(parent, node, "label", true);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.trace("Classes = ", node.class);
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const totalWidth = 0;
  const totalHeight = 0;
  rect2.attr("width", totalWidth).attr("height", totalHeight);
  shapeSvg.attr("class", "label edgeLabel");
  if (node.props) {
    const propKeys = new Set(Object.keys(node.props));
    if (node.props.borders) {
      applyNodePropertyBorders(rect2, node.props.borders, totalWidth, totalHeight);
      propKeys.delete("borders");
    }
    propKeys.forEach((propKey) => {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`Unknown node property ${propKey}`);
    });
  }
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
function applyNodePropertyBorders(rect2, borders, totalWidth, totalHeight) {
  const strokeDashArray = [];
  const addBorder = (length) => {
    strokeDashArray.push(length, 0);
  };
  const skipBorder = (length) => {
    strokeDashArray.push(0, length);
  };
  if (borders.includes("t")) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.debug("add top border");
    addBorder(totalWidth);
  } else {
    skipBorder(totalWidth);
  }
  if (borders.includes("r")) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.debug("add right border");
    addBorder(totalHeight);
  } else {
    skipBorder(totalHeight);
  }
  if (borders.includes("b")) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.debug("add bottom border");
    addBorder(totalWidth);
  } else {
    skipBorder(totalWidth);
  }
  if (borders.includes("l")) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.debug("add left border");
    addBorder(totalHeight);
  } else {
    skipBorder(totalHeight);
  }
  rect2.attr("stroke-dasharray", strokeDashArray.join(" "));
}
const rectWithTitle = (parent, node) => {
  let classes;
  if (!node.classes) {
    classes = "node default";
  } else {
    classes = "node " + node.classes;
  }
  const shapeSvg = parent.insert("g").attr("class", classes).attr("id", node.domId || node.id);
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const innerLine = shapeSvg.insert("line");
  const label = shapeSvg.insert("g").attr("class", "label");
  const text2 = node.labelText.flat ? node.labelText.flat() : node.labelText;
  let title = "";
  if (typeof text2 === "object") {
    title = text2[0];
  } else {
    title = text2;
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Label text abc79", title, text2, typeof text2 === "object");
  const text = label.node().appendChild(createLabel$1(title, node.labelStyle, true, true));
  let bbox = { width: 0, height: 0 };
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    const div = text.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(text);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Text 2", text2);
  const textRows = text2.slice(1, text2.length);
  let titleBox = text.getBBox();
  const descr = label.node().appendChild(
    createLabel$1(textRows.join ? textRows.join("<br/>") : textRows, node.labelStyle, true, true)
  );
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    const div = descr.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(descr);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  const halfPadding = node.padding / 2;
  (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(descr).attr(
    "transform",
    "translate( " + // (titleBox.width - bbox.width) / 2 +
    (bbox.width > titleBox.width ? 0 : (titleBox.width - bbox.width) / 2) + ", " + (titleBox.height + halfPadding + 5) + ")"
  );
  (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(text).attr(
    "transform",
    "translate( " + // (titleBox.width - bbox.width) / 2 +
    (bbox.width < titleBox.width ? 0 : -(titleBox.width - bbox.width) / 2) + ", 0)"
  );
  bbox = label.node().getBBox();
  label.attr(
    "transform",
    "translate(" + -bbox.width / 2 + ", " + (-bbox.height / 2 - halfPadding + 3) + ")"
  );
  rect2.attr("class", "outer title-state").attr("x", -bbox.width / 2 - halfPadding).attr("y", -bbox.height / 2 - halfPadding).attr("width", bbox.width + node.padding).attr("height", bbox.height + node.padding);
  innerLine.attr("class", "divider").attr("x1", -bbox.width / 2 - halfPadding).attr("x2", bbox.width / 2 + halfPadding).attr("y1", -bbox.height / 2 - halfPadding + titleBox.height + halfPadding).attr("y2", -bbox.height / 2 - halfPadding + titleBox.height + halfPadding);
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const stadium = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const h = bbox.height + node.padding;
  const w = bbox.width + h / 4 + node.padding;
  const rect2 = shapeSvg.insert("rect", ":first-child").attr("style", node.style).attr("rx", h / 2).attr("ry", h / 2).attr("x", -w / 2).attr("y", -h / 2).attr("width", w).attr("height", h);
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const circle = async (parent, node) => {
  const { shapeSvg, bbox, halfPadding } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const circle2 = shapeSvg.insert("circle", ":first-child");
  circle2.attr("style", node.style).attr("rx", node.rx).attr("ry", node.ry).attr("r", bbox.width / 2 + halfPadding).attr("width", bbox.width + node.padding).attr("height", bbox.height + node.padding);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Circle main");
  updateNodeBounds(node, circle2);
  node.intersect = function(point2) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Circle intersect", node, bbox.width / 2 + halfPadding, point2);
    return intersect.circle(node, bbox.width / 2 + halfPadding, point2);
  };
  return shapeSvg;
};
const doublecircle = async (parent, node) => {
  const { shapeSvg, bbox, halfPadding } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const gap = 5;
  const circleGroup = shapeSvg.insert("g", ":first-child");
  const outerCircle = circleGroup.insert("circle");
  const innerCircle = circleGroup.insert("circle");
  circleGroup.attr("class", node.class);
  outerCircle.attr("style", node.style).attr("rx", node.rx).attr("ry", node.ry).attr("r", bbox.width / 2 + halfPadding + gap).attr("width", bbox.width + node.padding + gap * 2).attr("height", bbox.height + node.padding + gap * 2);
  innerCircle.attr("style", node.style).attr("rx", node.rx).attr("ry", node.ry).attr("r", bbox.width / 2 + halfPadding).attr("width", bbox.width + node.padding).attr("height", bbox.height + node.padding);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("DoubleCircle main");
  updateNodeBounds(node, outerCircle);
  node.intersect = function(point2) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("DoubleCircle intersect", node, bbox.width / 2 + halfPadding + gap, point2);
    return intersect.circle(node, bbox.width / 2 + halfPadding + gap, point2);
  };
  return shapeSvg;
};
const subroutine = async (parent, node) => {
  const { shapeSvg, bbox } = await labelHelper(
    parent,
    node,
    getClassesFromNode(node, void 0),
    true
  );
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: 0, y: 0 },
    { x: w, y: 0 },
    { x: w, y: -h },
    { x: 0, y: -h },
    { x: 0, y: 0 },
    { x: -8, y: 0 },
    { x: w + 8, y: 0 },
    { x: w + 8, y: -h },
    { x: -8, y: -h },
    { x: -8, y: 0 }
  ];
  const el = insertPolygonShape(shapeSvg, w, h, points);
  el.attr("style", node.style);
  updateNodeBounds(node, el);
  node.intersect = function(point2) {
    return intersect.polygon(node, points, point2);
  };
  return shapeSvg;
};
const start = (parent, node) => {
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  const circle2 = shapeSvg.insert("circle", ":first-child");
  circle2.attr("class", "state-start").attr("r", 7).attr("width", 14).attr("height", 14);
  updateNodeBounds(node, circle2);
  node.intersect = function(point2) {
    return intersect.circle(node, 7, point2);
  };
  return shapeSvg;
};
const forkJoin = (parent, node, dir) => {
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  let width = 70;
  let height = 10;
  if (dir === "LR") {
    width = 10;
    height = 70;
  }
  const shape = shapeSvg.append("rect").attr("x", -1 * width / 2).attr("y", -1 * height / 2).attr("width", width).attr("height", height).attr("class", "fork-join");
  updateNodeBounds(node, shape);
  node.height = node.height + node.padding / 2;
  node.width = node.width + node.padding / 2;
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const end = (parent, node) => {
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  const innerCircle = shapeSvg.insert("circle", ":first-child");
  const circle2 = shapeSvg.insert("circle", ":first-child");
  circle2.attr("class", "state-start").attr("r", 7).attr("width", 14).attr("height", 14);
  innerCircle.attr("class", "state-end").attr("r", 5).attr("width", 10).attr("height", 10);
  updateNodeBounds(node, circle2);
  node.intersect = function(point2) {
    return intersect.circle(node, 7, point2);
  };
  return shapeSvg;
};
const class_box = (parent, node) => {
  const halfPadding = node.padding / 2;
  const rowPadding = 4;
  const lineHeight = 8;
  let classes;
  if (!node.classes) {
    classes = "node default";
  } else {
    classes = "node " + node.classes;
  }
  const shapeSvg = parent.insert("g").attr("class", classes).attr("id", node.domId || node.id);
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const topLine = shapeSvg.insert("line");
  const bottomLine = shapeSvg.insert("line");
  let maxWidth = 0;
  let maxHeight = rowPadding;
  const labelContainer = shapeSvg.insert("g").attr("class", "label");
  let verticalPos = 0;
  const hasInterface = node.classData.annotations && node.classData.annotations[0];
  const interfaceLabelText = node.classData.annotations[0] ? "«" + node.classData.annotations[0] + "»" : "";
  const interfaceLabel = labelContainer.node().appendChild(createLabel$1(interfaceLabelText, node.labelStyle, true, true));
  let interfaceBBox = interfaceLabel.getBBox();
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    const div = interfaceLabel.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(interfaceLabel);
    interfaceBBox = div.getBoundingClientRect();
    dv.attr("width", interfaceBBox.width);
    dv.attr("height", interfaceBBox.height);
  }
  if (node.classData.annotations[0]) {
    maxHeight += interfaceBBox.height + rowPadding;
    maxWidth += interfaceBBox.width;
  }
  let classTitleString = node.classData.label;
  if (node.classData.type !== void 0 && node.classData.type !== "") {
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels) {
      classTitleString += "&lt;" + node.classData.type + "&gt;";
    } else {
      classTitleString += "<" + node.classData.type + ">";
    }
  }
  const classTitleLabel = labelContainer.node().appendChild(createLabel$1(classTitleString, node.labelStyle, true, true));
  (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(classTitleLabel).attr("class", "classTitle");
  let classTitleBBox = classTitleLabel.getBBox();
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
    const div = classTitleLabel.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(classTitleLabel);
    classTitleBBox = div.getBoundingClientRect();
    dv.attr("width", classTitleBBox.width);
    dv.attr("height", classTitleBBox.height);
  }
  maxHeight += classTitleBBox.height + rowPadding;
  if (classTitleBBox.width > maxWidth) {
    maxWidth = classTitleBBox.width;
  }
  const classAttributes = [];
  node.classData.members.forEach((member) => {
    const parsedInfo = member.getDisplayDetails();
    let parsedText = parsedInfo.displayText;
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels) {
      parsedText = parsedText.replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }
    const lbl = labelContainer.node().appendChild(
      createLabel$1(
        parsedText,
        parsedInfo.cssStyle ? parsedInfo.cssStyle : node.labelStyle,
        true,
        true
      )
    );
    let bbox = lbl.getBBox();
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
      const div = lbl.children[0];
      const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(lbl);
      bbox = div.getBoundingClientRect();
      dv.attr("width", bbox.width);
      dv.attr("height", bbox.height);
    }
    if (bbox.width > maxWidth) {
      maxWidth = bbox.width;
    }
    maxHeight += bbox.height + rowPadding;
    classAttributes.push(lbl);
  });
  maxHeight += lineHeight;
  const classMethods = [];
  node.classData.methods.forEach((member) => {
    const parsedInfo = member.getDisplayDetails();
    let displayText = parsedInfo.displayText;
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels) {
      displayText = displayText.replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }
    const lbl = labelContainer.node().appendChild(
      createLabel$1(
        displayText,
        parsedInfo.cssStyle ? parsedInfo.cssStyle : node.labelStyle,
        true,
        true
      )
    );
    let bbox = lbl.getBBox();
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels)) {
      const div = lbl.children[0];
      const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(lbl);
      bbox = div.getBoundingClientRect();
      dv.attr("width", bbox.width);
      dv.attr("height", bbox.height);
    }
    if (bbox.width > maxWidth) {
      maxWidth = bbox.width;
    }
    maxHeight += bbox.height + rowPadding;
    classMethods.push(lbl);
  });
  maxHeight += lineHeight;
  if (hasInterface) {
    let diffX2 = (maxWidth - interfaceBBox.width) / 2;
    (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(interfaceLabel).attr(
      "transform",
      "translate( " + (-1 * maxWidth / 2 + diffX2) + ", " + -1 * maxHeight / 2 + ")"
    );
    verticalPos = interfaceBBox.height + rowPadding;
  }
  let diffX = (maxWidth - classTitleBBox.width) / 2;
  (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(classTitleLabel).attr(
    "transform",
    "translate( " + (-1 * maxWidth / 2 + diffX) + ", " + (-1 * maxHeight / 2 + verticalPos) + ")"
  );
  verticalPos += classTitleBBox.height + rowPadding;
  topLine.attr("class", "divider").attr("x1", -maxWidth / 2 - halfPadding).attr("x2", maxWidth / 2 + halfPadding).attr("y1", -maxHeight / 2 - halfPadding + lineHeight + verticalPos).attr("y2", -maxHeight / 2 - halfPadding + lineHeight + verticalPos);
  verticalPos += lineHeight;
  classAttributes.forEach((lbl) => {
    (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(lbl).attr(
      "transform",
      "translate( " + -maxWidth / 2 + ", " + (-1 * maxHeight / 2 + verticalPos + lineHeight / 2) + ")"
    );
    const memberBBox = lbl == null ? void 0 : lbl.getBBox();
    verticalPos += ((memberBBox == null ? void 0 : memberBBox.height) ?? 0) + rowPadding;
  });
  verticalPos += lineHeight;
  bottomLine.attr("class", "divider").attr("x1", -maxWidth / 2 - halfPadding).attr("x2", maxWidth / 2 + halfPadding).attr("y1", -maxHeight / 2 - halfPadding + lineHeight + verticalPos).attr("y2", -maxHeight / 2 - halfPadding + lineHeight + verticalPos);
  verticalPos += lineHeight;
  classMethods.forEach((lbl) => {
    (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(lbl).attr(
      "transform",
      "translate( " + -maxWidth / 2 + ", " + (-1 * maxHeight / 2 + verticalPos) + ")"
    );
    const memberBBox = lbl == null ? void 0 : lbl.getBBox();
    verticalPos += ((memberBBox == null ? void 0 : memberBBox.height) ?? 0) + rowPadding;
  });
  rect2.attr("style", node.style).attr("class", "outer title-state").attr("x", -maxWidth / 2 - halfPadding).attr("y", -(maxHeight / 2) - halfPadding).attr("width", maxWidth + node.padding).attr("height", maxHeight + node.padding);
  updateNodeBounds(node, rect2);
  node.intersect = function(point2) {
    return intersect.rect(node, point2);
  };
  return shapeSvg;
};
const shapes = {
  rhombus: question,
  question,
  rect,
  labelRect,
  rectWithTitle,
  choice,
  circle,
  doublecircle,
  stadium,
  hexagon,
  rect_left_inv_arrow,
  lean_right,
  lean_left,
  trapezoid,
  inv_trapezoid,
  rect_right_inv_arrow,
  cylinder,
  start,
  end,
  note: note$1,
  subroutine,
  fork: forkJoin,
  join: forkJoin,
  class_box
};
let nodeElems = {};
const insertNode = async (elem, node, dir) => {
  let newEl;
  let el;
  if (node.link) {
    let target;
    if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().securityLevel === "sandbox") {
      target = "_top";
    } else if (node.linkTarget) {
      target = node.linkTarget || "_blank";
    }
    newEl = elem.insert("svg:a").attr("xlink:href", node.link).attr("target", target);
    el = await shapes[node.shape](newEl, node, dir);
  } else {
    el = await shapes[node.shape](elem, node, dir);
    newEl = el;
  }
  if (node.tooltip) {
    el.attr("title", node.tooltip);
  }
  if (node.class) {
    el.attr("class", "node default " + node.class);
  }
  nodeElems[node.id] = newEl;
  if (node.haveCallback) {
    nodeElems[node.id].attr("class", nodeElems[node.id].attr("class") + " clickable");
  }
  return newEl;
};
const setNodeElem = (elem, node) => {
  nodeElems[node.id] = elem;
};
const clear$1 = () => {
  nodeElems = {};
};
const positionNode = (node) => {
  const el = nodeElems[node.id];
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.trace(
    "Transforming node",
    node.diff,
    node,
    "translate(" + (node.x - node.width / 2 - 5) + ", " + node.width / 2 + ")"
  );
  const padding = 8;
  const diff = node.diff || 0;
  if (node.clusterNode) {
    el.attr(
      "transform",
      "translate(" + (node.x + diff - node.width / 2) + ", " + (node.y - node.height / 2 - padding) + ")"
    );
  } else {
    el.attr("transform", "translate(" + node.x + ", " + node.y + ")");
  }
  return diff;
};
const markerOffsets = {
  aggregation: 18,
  extension: 18,
  composition: 18,
  dependency: 6,
  lollipop: 13.5,
  arrow_point: 5.3
};
function calculateDeltaAndAngle(point1, point2) {
  if (point1 === void 0 || point2 === void 0) {
    return { angle: 0, deltaX: 0, deltaY: 0 };
  }
  point1 = pointTransformer(point1);
  point2 = pointTransformer(point2);
  const [x1, y1] = [point1.x, point1.y];
  const [x2, y2] = [point2.x, point2.y];
  const deltaX = x2 - x1;
  const deltaY = y2 - y1;
  return { angle: Math.atan(deltaY / deltaX), deltaX, deltaY };
}
const pointTransformer = (data) => {
  if (Array.isArray(data)) {
    return { x: data[0], y: data[1] };
  }
  return data;
};
const getLineFunctionsWithOffset = (edge) => {
  return {
    x: function(d, i, data) {
      let offset = 0;
      if (i === 0 && Object.hasOwn(markerOffsets, edge.arrowTypeStart)) {
        const { angle, deltaX } = calculateDeltaAndAngle(data[0], data[1]);
        offset = markerOffsets[edge.arrowTypeStart] * Math.cos(angle) * (deltaX >= 0 ? 1 : -1);
      } else if (i === data.length - 1 && Object.hasOwn(markerOffsets, edge.arrowTypeEnd)) {
        const { angle, deltaX } = calculateDeltaAndAngle(
          data[data.length - 1],
          data[data.length - 2]
        );
        offset = markerOffsets[edge.arrowTypeEnd] * Math.cos(angle) * (deltaX >= 0 ? 1 : -1);
      }
      return pointTransformer(d).x + offset;
    },
    y: function(d, i, data) {
      let offset = 0;
      if (i === 0 && Object.hasOwn(markerOffsets, edge.arrowTypeStart)) {
        const { angle, deltaY } = calculateDeltaAndAngle(data[0], data[1]);
        offset = markerOffsets[edge.arrowTypeStart] * Math.abs(Math.sin(angle)) * (deltaY >= 0 ? 1 : -1);
      } else if (i === data.length - 1 && Object.hasOwn(markerOffsets, edge.arrowTypeEnd)) {
        const { angle, deltaY } = calculateDeltaAndAngle(
          data[data.length - 1],
          data[data.length - 2]
        );
        offset = markerOffsets[edge.arrowTypeEnd] * Math.abs(Math.sin(angle)) * (deltaY >= 0 ? 1 : -1);
      }
      return pointTransformer(d).y + offset;
    }
  };
};
const addEdgeMarkers = (svgPath, edge, url, id, diagramType) => {
  if (edge.arrowTypeStart) {
    addEdgeMarker(svgPath, "start", edge.arrowTypeStart, url, id, diagramType);
  }
  if (edge.arrowTypeEnd) {
    addEdgeMarker(svgPath, "end", edge.arrowTypeEnd, url, id, diagramType);
  }
};
const arrowTypesMap = {
  arrow_cross: "cross",
  arrow_point: "point",
  arrow_barb: "barb",
  arrow_circle: "circle",
  aggregation: "aggregation",
  extension: "extension",
  composition: "composition",
  dependency: "dependency",
  lollipop: "lollipop"
};
const addEdgeMarker = (svgPath, position, arrowType, url, id, diagramType) => {
  const endMarkerType = arrowTypesMap[arrowType];
  if (!endMarkerType) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`Unknown arrow type: ${arrowType}`);
    return;
  }
  const suffix = position === "start" ? "Start" : "End";
  svgPath.attr(`marker-${position}`, `url(${url}#${id}_${diagramType}-${endMarkerType}${suffix})`);
};
let edgeLabels = {};
let terminalLabels = {};
const clear = () => {
  edgeLabels = {};
  terminalLabels = {};
};
const insertEdgeLabel = (elem, edge) => {
  const useHtmlLabels = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.m)((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels);
  const labelElement = edge.labelType === "markdown" ? (0,_createText_a48a4c90_js__WEBPACK_IMPORTED_MODULE_2__.a)(elem, edge.label, {
    style: edge.labelStyle,
    useHtmlLabels,
    addSvgBackground: true
  }) : createLabel$1(edge.label, edge.labelStyle);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("abc82", edge, edge.labelType);
  const edgeLabel = elem.insert("g").attr("class", "edgeLabel");
  const label = edgeLabel.insert("g").attr("class", "label");
  label.node().appendChild(labelElement);
  let bbox = labelElement.getBBox();
  if (useHtmlLabels) {
    const div = labelElement.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(labelElement);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  label.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  edgeLabels[edge.id] = edgeLabel;
  edge.width = bbox.width;
  edge.height = bbox.height;
  let fo;
  if (edge.startLabelLeft) {
    const startLabelElement = createLabel$1(edge.startLabelLeft, edge.labelStyle);
    const startEdgeLabelLeft = elem.insert("g").attr("class", "edgeTerminals");
    const inner = startEdgeLabelLeft.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(startLabelElement);
    const slBox = startLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    if (!terminalLabels[edge.id]) {
      terminalLabels[edge.id] = {};
    }
    terminalLabels[edge.id].startLeft = startEdgeLabelLeft;
    setTerminalWidth(fo, edge.startLabelLeft);
  }
  if (edge.startLabelRight) {
    const startLabelElement = createLabel$1(edge.startLabelRight, edge.labelStyle);
    const startEdgeLabelRight = elem.insert("g").attr("class", "edgeTerminals");
    const inner = startEdgeLabelRight.insert("g").attr("class", "inner");
    fo = startEdgeLabelRight.node().appendChild(startLabelElement);
    inner.node().appendChild(startLabelElement);
    const slBox = startLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    if (!terminalLabels[edge.id]) {
      terminalLabels[edge.id] = {};
    }
    terminalLabels[edge.id].startRight = startEdgeLabelRight;
    setTerminalWidth(fo, edge.startLabelRight);
  }
  if (edge.endLabelLeft) {
    const endLabelElement = createLabel$1(edge.endLabelLeft, edge.labelStyle);
    const endEdgeLabelLeft = elem.insert("g").attr("class", "edgeTerminals");
    const inner = endEdgeLabelLeft.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(endLabelElement);
    const slBox = endLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    endEdgeLabelLeft.node().appendChild(endLabelElement);
    if (!terminalLabels[edge.id]) {
      terminalLabels[edge.id] = {};
    }
    terminalLabels[edge.id].endLeft = endEdgeLabelLeft;
    setTerminalWidth(fo, edge.endLabelLeft);
  }
  if (edge.endLabelRight) {
    const endLabelElement = createLabel$1(edge.endLabelRight, edge.labelStyle);
    const endEdgeLabelRight = elem.insert("g").attr("class", "edgeTerminals");
    const inner = endEdgeLabelRight.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(endLabelElement);
    const slBox = endLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    endEdgeLabelRight.node().appendChild(endLabelElement);
    if (!terminalLabels[edge.id]) {
      terminalLabels[edge.id] = {};
    }
    terminalLabels[edge.id].endRight = endEdgeLabelRight;
    setTerminalWidth(fo, edge.endLabelRight);
  }
  return labelElement;
};
function setTerminalWidth(fo, value) {
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.htmlLabels && fo) {
    fo.style.width = value.length * 9 + "px";
    fo.style.height = "12px";
  }
}
const positionEdgeLabel = (edge, paths) => {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("Moving label abc78 ", edge.id, edge.label, edgeLabels[edge.id]);
  let path = paths.updatedPath ? paths.updatedPath : paths.originalPath;
  if (edge.label) {
    const el = edgeLabels[edge.id];
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.u.calcLabelPosition(path);
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info(
        "Moving label " + edge.label + " from (",
        x,
        ",",
        y,
        ") to (",
        pos.x,
        ",",
        pos.y,
        ") abc78"
      );
      if (paths.updatedPath) {
        x = pos.x;
        y = pos.y;
      }
    }
    el.attr("transform", "translate(" + x + ", " + y + ")");
  }
  if (edge.startLabelLeft) {
    const el = terminalLabels[edge.id].startLeft;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.u.calcTerminalLabelPosition(edge.arrowTypeStart ? 10 : 0, "start_left", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", "translate(" + x + ", " + y + ")");
  }
  if (edge.startLabelRight) {
    const el = terminalLabels[edge.id].startRight;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.u.calcTerminalLabelPosition(
        edge.arrowTypeStart ? 10 : 0,
        "start_right",
        path
      );
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", "translate(" + x + ", " + y + ")");
  }
  if (edge.endLabelLeft) {
    const el = terminalLabels[edge.id].endLeft;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.u.calcTerminalLabelPosition(edge.arrowTypeEnd ? 10 : 0, "end_left", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", "translate(" + x + ", " + y + ")");
  }
  if (edge.endLabelRight) {
    const el = terminalLabels[edge.id].endRight;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.u.calcTerminalLabelPosition(edge.arrowTypeEnd ? 10 : 0, "end_right", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", "translate(" + x + ", " + y + ")");
  }
};
const outsideNode = (node, point2) => {
  const x = node.x;
  const y = node.y;
  const dx = Math.abs(point2.x - x);
  const dy = Math.abs(point2.y - y);
  const w = node.width / 2;
  const h = node.height / 2;
  if (dx >= w || dy >= h) {
    return true;
  }
  return false;
};
const intersection = (node, outsidePoint, insidePoint) => {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`intersection calc abc89:
  outsidePoint: ${JSON.stringify(outsidePoint)}
  insidePoint : ${JSON.stringify(insidePoint)}
  node        : x:${node.x} y:${node.y} w:${node.width} h:${node.height}`);
  const x = node.x;
  const y = node.y;
  const dx = Math.abs(x - insidePoint.x);
  const w = node.width / 2;
  let r = insidePoint.x < outsidePoint.x ? w - dx : w + dx;
  const h = node.height / 2;
  const Q = Math.abs(outsidePoint.y - insidePoint.y);
  const R = Math.abs(outsidePoint.x - insidePoint.x);
  if (Math.abs(y - outsidePoint.y) * w > Math.abs(x - outsidePoint.x) * h) {
    let q = insidePoint.y < outsidePoint.y ? outsidePoint.y - h - y : y - h - outsidePoint.y;
    r = R * q / Q;
    const res = {
      x: insidePoint.x < outsidePoint.x ? insidePoint.x + r : insidePoint.x - R + r,
      y: insidePoint.y < outsidePoint.y ? insidePoint.y + Q - q : insidePoint.y - Q + q
    };
    if (r === 0) {
      res.x = outsidePoint.x;
      res.y = outsidePoint.y;
    }
    if (R === 0) {
      res.x = outsidePoint.x;
    }
    if (Q === 0) {
      res.y = outsidePoint.y;
    }
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`abc89 topp/bott calc, Q ${Q}, q ${q}, R ${R}, r ${r}`, res);
    return res;
  } else {
    if (insidePoint.x < outsidePoint.x) {
      r = outsidePoint.x - w - x;
    } else {
      r = x - w - outsidePoint.x;
    }
    let q = Q * r / R;
    let _x = insidePoint.x < outsidePoint.x ? insidePoint.x + R - r : insidePoint.x - R + r;
    let _y = insidePoint.y < outsidePoint.y ? insidePoint.y + q : insidePoint.y - q;
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn(`sides calc abc89, Q ${Q}, q ${q}, R ${R}, r ${r}`, { _x, _y });
    if (r === 0) {
      _x = outsidePoint.x;
      _y = outsidePoint.y;
    }
    if (R === 0) {
      _x = outsidePoint.x;
    }
    if (Q === 0) {
      _y = outsidePoint.y;
    }
    return { x: _x, y: _y };
  }
};
const cutPathAtIntersect = (_points, boundryNode) => {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 cutPathAtIntersect", _points, boundryNode);
  let points = [];
  let lastPointOutside = _points[0];
  let isInside = false;
  _points.forEach((point2) => {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("abc88 checking point", point2, boundryNode);
    if (!outsideNode(boundryNode, point2) && !isInside) {
      const inter = intersection(boundryNode, lastPointOutside, point2);
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 inside", point2, lastPointOutside, inter);
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 intersection", inter);
      let pointPresent = false;
      points.forEach((p) => {
        pointPresent = pointPresent || p.x === inter.x && p.y === inter.y;
      });
      if (!points.some((e) => e.x === inter.x && e.y === inter.y)) {
        points.push(inter);
      } else {
        _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 no intersect", inter, points);
      }
      isInside = true;
    } else {
      _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 outside", point2, lastPointOutside);
      lastPointOutside = point2;
      if (!isInside) {
        points.push(point2);
      }
    }
  });
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.warn("abc88 returning points", points);
  return points;
};
const insertEdge = function(elem, e, edge, clusterDb, diagramType, graph, id) {
  let points = edge.points;
  let pointsHasChanged = false;
  const tail = graph.node(e.v);
  var head = graph.node(e.w);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("abc88 InsertEdge: ", edge);
  if (head.intersect && tail.intersect) {
    points = points.slice(1, edge.points.length - 1);
    points.unshift(tail.intersect(points[0]));
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info(
      "Last point",
      points[points.length - 1],
      head,
      head.intersect(points[points.length - 1])
    );
    points.push(head.intersect(points[points.length - 1]));
  }
  if (edge.toCluster) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("to cluster abc88", clusterDb[edge.toCluster]);
    points = cutPathAtIntersect(edge.points, clusterDb[edge.toCluster].node);
    pointsHasChanged = true;
  }
  if (edge.fromCluster) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("from cluster abc88", clusterDb[edge.fromCluster]);
    points = cutPathAtIntersect(points.reverse(), clusterDb[edge.fromCluster].node).reverse();
    pointsHasChanged = true;
  }
  const lineData = points.filter((p) => !Number.isNaN(p.y));
  let curve = d3__WEBPACK_IMPORTED_MODULE_0__/* .curveBasis */ .$0Z;
  if (edge.curve && (diagramType === "graph" || diagramType === "flowchart")) {
    curve = edge.curve;
  }
  const { x, y } = getLineFunctionsWithOffset(edge);
  const lineFunction = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .line */ .jvg)().x(x).y(y).curve(curve);
  let strokeClasses;
  switch (edge.thickness) {
    case "normal":
      strokeClasses = "edge-thickness-normal";
      break;
    case "thick":
      strokeClasses = "edge-thickness-thick";
      break;
    case "invisible":
      strokeClasses = "edge-thickness-thick";
      break;
    default:
      strokeClasses = "";
  }
  switch (edge.pattern) {
    case "solid":
      strokeClasses += " edge-pattern-solid";
      break;
    case "dotted":
      strokeClasses += " edge-pattern-dotted";
      break;
    case "dashed":
      strokeClasses += " edge-pattern-dashed";
      break;
  }
  const svgPath = elem.append("path").attr("d", lineFunction(lineData)).attr("id", edge.id).attr("class", " " + strokeClasses + (edge.classes ? " " + edge.classes : "")).attr("style", edge.style);
  let url = "";
  if ((0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().flowchart.arrowMarkerAbsolute || (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.c)().state.arrowMarkerAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = url.replace(/\(/g, "\\(");
    url = url.replace(/\)/g, "\\)");
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("arrowTypeStart", edge.arrowTypeStart);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.l.info("arrowTypeEnd", edge.arrowTypeEnd);
  addEdgeMarkers(svgPath, edge, url, id, diagramType);
  let paths = {};
  if (pointsHasChanged) {
    paths.updatedPath = points;
  }
  paths.originalPath = edge.points;
  return paths;
};



/***/ })

}]);
//# sourceMappingURL=8443.214c35b34c68e0bcfa22.js.map?v=214c35b34c68e0bcfa22