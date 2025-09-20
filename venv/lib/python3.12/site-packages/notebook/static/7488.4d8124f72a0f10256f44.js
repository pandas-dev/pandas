"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7488],{

/***/ 26616:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _utils_index_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(90267);
/* harmony import */ var _color_index_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(42528);
/* IMPORT */


/* MAIN */
const channel = (color, channel) => {
    return _utils_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.lang.round(_color_index_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z.parse(color)[channel]);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (channel);


/***/ }),

/***/ 49479:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   j: () => (/* binding */ setupViewPortForSVG),
/* harmony export */   q: () => (/* binding */ getDiagramElement)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(83619);


// src/rendering-util/insertElementsForSize.js

var getDiagramElement = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((id, securityLevel) => {
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("body");
  const svg = root.select(`[id="${id}"]`);
  return svg;
}, "getDiagramElement");

// src/rendering-util/setupViewPortForSVG.ts
var setupViewPortForSVG = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((svg, padding, cssDiagram, useMaxWidth) => {
  svg.attr("class", cssDiagram);
  const { width, height, x, y } = calculateDimensionsWithPadding(svg, padding);
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .configureSvgSize */ .v2)(svg, height, width, useMaxWidth);
  const viewBox = createViewBox(x, y, width, height, padding);
  svg.attr("viewBox", viewBox);
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .log */ .cM.debug(`viewBox configured: ${viewBox} with padding: ${padding}`);
}, "setupViewPortForSVG");
var calculateDimensionsWithPadding = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((svg, padding) => {
  const bounds = svg.node()?.getBBox() || { width: 0, height: 0, x: 0, y: 0 };
  return {
    width: bounds.width + padding * 2,
    height: bounds.height + padding * 2,
    x: bounds.x,
    y: bounds.y
  };
}, "calculateDimensionsWithPadding");
var createViewBox = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((x, y, width, height, padding) => {
  return `${x - padding} ${y - padding} ${width} ${height}`;
}, "createViewBox");




/***/ }),

/***/ 57488:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(49479);
/* harmony import */ var _chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(83696);
/* harmony import */ var _chunk_IIMUDSI4_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(52017);
/* harmony import */ var _chunk_VV3M67IP_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(83133);
/* harmony import */ var _chunk_HRU6DDCH_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(41921);
/* harmony import */ var _chunk_K557N5IZ_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(81861);
/* harmony import */ var _chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(62072);
/* harmony import */ var _chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(39769);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(42626);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(83619);
/* harmony import */ var khroma__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(26616);
/* harmony import */ var khroma__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(14728);











// src/diagrams/er/parser/erDiagram.jison
var parser = function() {
  var o = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [6, 8, 10, 22, 24, 26, 28, 33, 34, 35, 36, 37, 40, 43, 44, 50], $V1 = [1, 10], $V2 = [1, 11], $V3 = [1, 12], $V4 = [1, 13], $V5 = [1, 20], $V6 = [1, 21], $V7 = [1, 22], $V8 = [1, 23], $V9 = [1, 24], $Va = [1, 19], $Vb = [1, 25], $Vc = [1, 26], $Vd = [1, 18], $Ve = [1, 33], $Vf = [1, 34], $Vg = [1, 35], $Vh = [1, 36], $Vi = [1, 37], $Vj = [6, 8, 10, 13, 15, 17, 20, 21, 22, 24, 26, 28, 33, 34, 35, 36, 37, 40, 43, 44, 50, 63, 64, 65, 66, 67], $Vk = [1, 42], $Vl = [1, 43], $Vm = [1, 52], $Vn = [40, 50, 68, 69], $Vo = [1, 63], $Vp = [1, 61], $Vq = [1, 58], $Vr = [1, 62], $Vs = [1, 64], $Vt = [6, 8, 10, 13, 17, 22, 24, 26, 28, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 48, 49, 50, 63, 64, 65, 66, 67], $Vu = [63, 64, 65, 66, 67], $Vv = [1, 81], $Vw = [1, 80], $Vx = [1, 78], $Vy = [1, 79], $Vz = [6, 10, 42, 47], $VA = [6, 10, 13, 41, 42, 47, 48, 49], $VB = [1, 89], $VC = [1, 88], $VD = [1, 87], $VE = [19, 56], $VF = [1, 98], $VG = [1, 97], $VH = [19, 56, 58, 60];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "ER_DIAGRAM": 4, "document": 5, "EOF": 6, "line": 7, "SPACE": 8, "statement": 9, "NEWLINE": 10, "entityName": 11, "relSpec": 12, "COLON": 13, "role": 14, "STYLE_SEPARATOR": 15, "idList": 16, "BLOCK_START": 17, "attributes": 18, "BLOCK_STOP": 19, "SQS": 20, "SQE": 21, "title": 22, "title_value": 23, "acc_title": 24, "acc_title_value": 25, "acc_descr": 26, "acc_descr_value": 27, "acc_descr_multiline_value": 28, "direction": 29, "classDefStatement": 30, "classStatement": 31, "styleStatement": 32, "direction_tb": 33, "direction_bt": 34, "direction_rl": 35, "direction_lr": 36, "CLASSDEF": 37, "stylesOpt": 38, "separator": 39, "UNICODE_TEXT": 40, "STYLE_TEXT": 41, "COMMA": 42, "CLASS": 43, "STYLE": 44, "style": 45, "styleComponent": 46, "SEMI": 47, "NUM": 48, "BRKT": 49, "ENTITY_NAME": 50, "attribute": 51, "attributeType": 52, "attributeName": 53, "attributeKeyTypeList": 54, "attributeComment": 55, "ATTRIBUTE_WORD": 56, "attributeKeyType": 57, ",": 58, "ATTRIBUTE_KEY": 59, "COMMENT": 60, "cardinality": 61, "relType": 62, "ZERO_OR_ONE": 63, "ZERO_OR_MORE": 64, "ONE_OR_MORE": 65, "ONLY_ONE": 66, "MD_PARENT": 67, "NON_IDENTIFYING": 68, "IDENTIFYING": 69, "WORD": 70, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "ER_DIAGRAM", 6: "EOF", 8: "SPACE", 10: "NEWLINE", 13: "COLON", 15: "STYLE_SEPARATOR", 17: "BLOCK_START", 19: "BLOCK_STOP", 20: "SQS", 21: "SQE", 22: "title", 23: "title_value", 24: "acc_title", 25: "acc_title_value", 26: "acc_descr", 27: "acc_descr_value", 28: "acc_descr_multiline_value", 33: "direction_tb", 34: "direction_bt", 35: "direction_rl", 36: "direction_lr", 37: "CLASSDEF", 40: "UNICODE_TEXT", 41: "STYLE_TEXT", 42: "COMMA", 43: "CLASS", 44: "STYLE", 47: "SEMI", 48: "NUM", 49: "BRKT", 50: "ENTITY_NAME", 56: "ATTRIBUTE_WORD", 58: ",", 59: "ATTRIBUTE_KEY", 60: "COMMENT", 63: "ZERO_OR_ONE", 64: "ZERO_OR_MORE", 65: "ONE_OR_MORE", 66: "ONLY_ONE", 67: "MD_PARENT", 68: "NON_IDENTIFYING", 69: "IDENTIFYING", 70: "WORD" },
    productions_: [0, [3, 3], [5, 0], [5, 2], [7, 2], [7, 1], [7, 1], [7, 1], [9, 5], [9, 9], [9, 7], [9, 7], [9, 4], [9, 6], [9, 3], [9, 5], [9, 1], [9, 3], [9, 7], [9, 9], [9, 6], [9, 8], [9, 4], [9, 6], [9, 2], [9, 2], [9, 2], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [29, 1], [29, 1], [29, 1], [29, 1], [30, 4], [16, 1], [16, 1], [16, 3], [16, 3], [31, 3], [32, 4], [38, 1], [38, 3], [45, 1], [45, 2], [39, 1], [39, 1], [39, 1], [46, 1], [46, 1], [46, 1], [46, 1], [11, 1], [11, 1], [18, 1], [18, 2], [51, 2], [51, 3], [51, 3], [51, 4], [52, 1], [53, 1], [54, 1], [54, 3], [57, 1], [55, 1], [12, 3], [61, 1], [61, 1], [61, 1], [61, 1], [61, 1], [62, 1], [62, 1], [14, 1], [14, 1], [14, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 1:
          break;
        case 2:
          this.$ = [];
          break;
        case 3:
          $$[$0 - 1].push($$[$0]);
          this.$ = $$[$0 - 1];
          break;
        case 4:
        case 5:
          this.$ = $$[$0];
          break;
        case 6:
        case 7:
          this.$ = [];
          break;
        case 8:
          yy.addEntity($$[$0 - 4]);
          yy.addEntity($$[$0 - 2]);
          yy.addRelationship($$[$0 - 4], $$[$0], $$[$0 - 2], $$[$0 - 3]);
          break;
        case 9:
          yy.addEntity($$[$0 - 8]);
          yy.addEntity($$[$0 - 4]);
          yy.addRelationship($$[$0 - 8], $$[$0], $$[$0 - 4], $$[$0 - 5]);
          yy.setClass([$$[$0 - 8]], $$[$0 - 6]);
          yy.setClass([$$[$0 - 4]], $$[$0 - 2]);
          break;
        case 10:
          yy.addEntity($$[$0 - 6]);
          yy.addEntity($$[$0 - 2]);
          yy.addRelationship($$[$0 - 6], $$[$0], $$[$0 - 2], $$[$0 - 3]);
          yy.setClass([$$[$0 - 6]], $$[$0 - 4]);
          break;
        case 11:
          yy.addEntity($$[$0 - 6]);
          yy.addEntity($$[$0 - 4]);
          yy.addRelationship($$[$0 - 6], $$[$0], $$[$0 - 4], $$[$0 - 5]);
          yy.setClass([$$[$0 - 4]], $$[$0 - 2]);
          break;
        case 12:
          yy.addEntity($$[$0 - 3]);
          yy.addAttributes($$[$0 - 3], $$[$0 - 1]);
          break;
        case 13:
          yy.addEntity($$[$0 - 5]);
          yy.addAttributes($$[$0 - 5], $$[$0 - 1]);
          yy.setClass([$$[$0 - 5]], $$[$0 - 3]);
          break;
        case 14:
          yy.addEntity($$[$0 - 2]);
          break;
        case 15:
          yy.addEntity($$[$0 - 4]);
          yy.setClass([$$[$0 - 4]], $$[$0 - 2]);
          break;
        case 16:
          yy.addEntity($$[$0]);
          break;
        case 17:
          yy.addEntity($$[$0 - 2]);
          yy.setClass([$$[$0 - 2]], $$[$0]);
          break;
        case 18:
          yy.addEntity($$[$0 - 6], $$[$0 - 4]);
          yy.addAttributes($$[$0 - 6], $$[$0 - 1]);
          break;
        case 19:
          yy.addEntity($$[$0 - 8], $$[$0 - 6]);
          yy.addAttributes($$[$0 - 8], $$[$0 - 1]);
          yy.setClass([$$[$0 - 8]], $$[$0 - 3]);
          break;
        case 20:
          yy.addEntity($$[$0 - 5], $$[$0 - 3]);
          break;
        case 21:
          yy.addEntity($$[$0 - 7], $$[$0 - 5]);
          yy.setClass([$$[$0 - 7]], $$[$0 - 2]);
          break;
        case 22:
          yy.addEntity($$[$0 - 3], $$[$0 - 1]);
          break;
        case 23:
          yy.addEntity($$[$0 - 5], $$[$0 - 3]);
          yy.setClass([$$[$0 - 5]], $$[$0]);
          break;
        case 24:
        case 25:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 26:
        case 27:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 32:
          yy.setDirection("TB");
          break;
        case 33:
          yy.setDirection("BT");
          break;
        case 34:
          yy.setDirection("RL");
          break;
        case 35:
          yy.setDirection("LR");
          break;
        case 36:
          this.$ = $$[$0 - 3];
          yy.addClass($$[$0 - 2], $$[$0 - 1]);
          break;
        case 37:
        case 38:
        case 56:
        case 64:
          this.$ = [$$[$0]];
          break;
        case 39:
        case 40:
          this.$ = $$[$0 - 2].concat([$$[$0]]);
          break;
        case 41:
          this.$ = $$[$0 - 2];
          yy.setClass($$[$0 - 1], $$[$0]);
          break;
        case 42:
          ;
          this.$ = $$[$0 - 3];
          yy.addCssStyles($$[$0 - 2], $$[$0 - 1]);
          break;
        case 43:
          this.$ = [$$[$0]];
          break;
        case 44:
          $$[$0 - 2].push($$[$0]);
          this.$ = $$[$0 - 2];
          break;
        case 46:
          this.$ = $$[$0 - 1] + $$[$0];
          break;
        case 54:
        case 76:
        case 77:
          this.$ = $$[$0].replace(/"/g, "");
          break;
        case 55:
        case 78:
          this.$ = $$[$0];
          break;
        case 57:
          $$[$0].push($$[$0 - 1]);
          this.$ = $$[$0];
          break;
        case 58:
          this.$ = { type: $$[$0 - 1], name: $$[$0] };
          break;
        case 59:
          this.$ = { type: $$[$0 - 2], name: $$[$0 - 1], keys: $$[$0] };
          break;
        case 60:
          this.$ = { type: $$[$0 - 2], name: $$[$0 - 1], comment: $$[$0] };
          break;
        case 61:
          this.$ = { type: $$[$0 - 3], name: $$[$0 - 2], keys: $$[$0 - 1], comment: $$[$0] };
          break;
        case 62:
        case 63:
        case 66:
          this.$ = $$[$0];
          break;
        case 65:
          $$[$0 - 2].push($$[$0]);
          this.$ = $$[$0 - 2];
          break;
        case 67:
          this.$ = $$[$0].replace(/"/g, "");
          break;
        case 68:
          this.$ = { cardA: $$[$0], relType: $$[$0 - 1], cardB: $$[$0 - 2] };
          break;
        case 69:
          this.$ = yy.Cardinality.ZERO_OR_ONE;
          break;
        case 70:
          this.$ = yy.Cardinality.ZERO_OR_MORE;
          break;
        case 71:
          this.$ = yy.Cardinality.ONE_OR_MORE;
          break;
        case 72:
          this.$ = yy.Cardinality.ONLY_ONE;
          break;
        case 73:
          this.$ = yy.Cardinality.MD_PARENT;
          break;
        case 74:
          this.$ = yy.Identification.NON_IDENTIFYING;
          break;
        case 75:
          this.$ = yy.Identification.IDENTIFYING;
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: [1, 2] }, { 1: [3] }, o($V0, [2, 2], { 5: 3 }), { 6: [1, 4], 7: 5, 8: [1, 6], 9: 7, 10: [1, 8], 11: 9, 22: $V1, 24: $V2, 26: $V3, 28: $V4, 29: 14, 30: 15, 31: 16, 32: 17, 33: $V5, 34: $V6, 35: $V7, 36: $V8, 37: $V9, 40: $Va, 43: $Vb, 44: $Vc, 50: $Vd }, o($V0, [2, 7], { 1: [2, 1] }), o($V0, [2, 3]), { 9: 27, 11: 9, 22: $V1, 24: $V2, 26: $V3, 28: $V4, 29: 14, 30: 15, 31: 16, 32: 17, 33: $V5, 34: $V6, 35: $V7, 36: $V8, 37: $V9, 40: $Va, 43: $Vb, 44: $Vc, 50: $Vd }, o($V0, [2, 5]), o($V0, [2, 6]), o($V0, [2, 16], { 12: 28, 61: 32, 15: [1, 29], 17: [1, 30], 20: [1, 31], 63: $Ve, 64: $Vf, 65: $Vg, 66: $Vh, 67: $Vi }), { 23: [1, 38] }, { 25: [1, 39] }, { 27: [1, 40] }, o($V0, [2, 27]), o($V0, [2, 28]), o($V0, [2, 29]), o($V0, [2, 30]), o($V0, [2, 31]), o($Vj, [2, 54]), o($Vj, [2, 55]), o($V0, [2, 32]), o($V0, [2, 33]), o($V0, [2, 34]), o($V0, [2, 35]), { 16: 41, 40: $Vk, 41: $Vl }, { 16: 44, 40: $Vk, 41: $Vl }, { 16: 45, 40: $Vk, 41: $Vl }, o($V0, [2, 4]), { 11: 46, 40: $Va, 50: $Vd }, { 16: 47, 40: $Vk, 41: $Vl }, { 18: 48, 19: [1, 49], 51: 50, 52: 51, 56: $Vm }, { 11: 53, 40: $Va, 50: $Vd }, { 62: 54, 68: [1, 55], 69: [1, 56] }, o($Vn, [2, 69]), o($Vn, [2, 70]), o($Vn, [2, 71]), o($Vn, [2, 72]), o($Vn, [2, 73]), o($V0, [2, 24]), o($V0, [2, 25]), o($V0, [2, 26]), { 13: $Vo, 38: 57, 41: $Vp, 42: $Vq, 45: 59, 46: 60, 48: $Vr, 49: $Vs }, o($Vt, [2, 37]), o($Vt, [2, 38]), { 16: 65, 40: $Vk, 41: $Vl, 42: $Vq }, { 13: $Vo, 38: 66, 41: $Vp, 42: $Vq, 45: 59, 46: 60, 48: $Vr, 49: $Vs }, { 13: [1, 67], 15: [1, 68] }, o($V0, [2, 17], { 61: 32, 12: 69, 17: [1, 70], 42: $Vq, 63: $Ve, 64: $Vf, 65: $Vg, 66: $Vh, 67: $Vi }), { 19: [1, 71] }, o($V0, [2, 14]), { 18: 72, 19: [2, 56], 51: 50, 52: 51, 56: $Vm }, { 53: 73, 56: [1, 74] }, { 56: [2, 62] }, { 21: [1, 75] }, { 61: 76, 63: $Ve, 64: $Vf, 65: $Vg, 66: $Vh, 67: $Vi }, o($Vu, [2, 74]), o($Vu, [2, 75]), { 6: $Vv, 10: $Vw, 39: 77, 42: $Vx, 47: $Vy }, { 40: [1, 82], 41: [1, 83] }, o($Vz, [2, 43], { 46: 84, 13: $Vo, 41: $Vp, 48: $Vr, 49: $Vs }), o($VA, [2, 45]), o($VA, [2, 50]), o($VA, [2, 51]), o($VA, [2, 52]), o($VA, [2, 53]), o($V0, [2, 41], { 42: $Vq }), { 6: $Vv, 10: $Vw, 39: 85, 42: $Vx, 47: $Vy }, { 14: 86, 40: $VB, 50: $VC, 70: $VD }, { 16: 90, 40: $Vk, 41: $Vl }, { 11: 91, 40: $Va, 50: $Vd }, { 18: 92, 19: [1, 93], 51: 50, 52: 51, 56: $Vm }, o($V0, [2, 12]), { 19: [2, 57] }, o($VE, [2, 58], { 54: 94, 55: 95, 57: 96, 59: $VF, 60: $VG }), o([19, 56, 59, 60], [2, 63]), o($V0, [2, 22], { 15: [1, 100], 17: [1, 99] }), o([40, 50], [2, 68]), o($V0, [2, 36]), { 13: $Vo, 41: $Vp, 45: 101, 46: 60, 48: $Vr, 49: $Vs }, o($V0, [2, 47]), o($V0, [2, 48]), o($V0, [2, 49]), o($Vt, [2, 39]), o($Vt, [2, 40]), o($VA, [2, 46]), o($V0, [2, 42]), o($V0, [2, 8]), o($V0, [2, 76]), o($V0, [2, 77]), o($V0, [2, 78]), { 13: [1, 102], 42: $Vq }, { 13: [1, 104], 15: [1, 103] }, { 19: [1, 105] }, o($V0, [2, 15]), o($VE, [2, 59], { 55: 106, 58: [1, 107], 60: $VG }), o($VE, [2, 60]), o($VH, [2, 64]), o($VE, [2, 67]), o($VH, [2, 66]), { 18: 108, 19: [1, 109], 51: 50, 52: 51, 56: $Vm }, { 16: 110, 40: $Vk, 41: $Vl }, o($Vz, [2, 44], { 46: 84, 13: $Vo, 41: $Vp, 48: $Vr, 49: $Vs }), { 14: 111, 40: $VB, 50: $VC, 70: $VD }, { 16: 112, 40: $Vk, 41: $Vl }, { 14: 113, 40: $VB, 50: $VC, 70: $VD }, o($V0, [2, 13]), o($VE, [2, 61]), { 57: 114, 59: $VF }, { 19: [1, 115] }, o($V0, [2, 20]), o($V0, [2, 23], { 17: [1, 116], 42: $Vq }), o($V0, [2, 11]), { 13: [1, 117], 42: $Vq }, o($V0, [2, 10]), o($VH, [2, 65]), o($V0, [2, 18]), { 18: 118, 19: [1, 119], 51: 50, 52: 51, 56: $Vm }, { 14: 120, 40: $VB, 50: $VC, 70: $VD }, { 19: [1, 121] }, o($V0, [2, 21]), o($V0, [2, 9]), o($V0, [2, 19])],
    defaultActions: { 52: [2, 62], 72: [2, 57] },
    parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parse(input) {
      var self = this, stack = [0], tstack = [], vstack = [null], lstack = [], table = this.table, yytext = "", yylineno = 0, yyleng = 0, recovering = 0, TERROR = 2, EOF = 1;
      var args = lstack.slice.call(arguments, 1);
      var lexer2 = Object.create(this.lexer);
      var sharedState = { yy: {} };
      for (var k in this.yy) {
        if (Object.prototype.hasOwnProperty.call(this.yy, k)) {
          sharedState.yy[k] = this.yy[k];
        }
      }
      lexer2.setInput(input, sharedState.yy);
      sharedState.yy.lexer = lexer2;
      sharedState.yy.parser = this;
      if (typeof lexer2.yylloc == "undefined") {
        lexer2.yylloc = {};
      }
      var yyloc = lexer2.yylloc;
      lstack.push(yyloc);
      var ranges = lexer2.options && lexer2.options.ranges;
      if (typeof sharedState.yy.parseError === "function") {
        this.parseError = sharedState.yy.parseError;
      } else {
        this.parseError = Object.getPrototypeOf(this).parseError;
      }
      function popStack(n) {
        stack.length = stack.length - 2 * n;
        vstack.length = vstack.length - n;
        lstack.length = lstack.length - n;
      }
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(popStack, "popStack");
      function lex() {
        var token;
        token = tstack.pop() || lexer2.lex() || EOF;
        if (typeof token !== "number") {
          if (token instanceof Array) {
            tstack = token;
            token = tstack.pop();
          }
          token = self.symbols_[token] || token;
        }
        return token;
      }
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(lex, "lex");
      var symbol, preErrorSymbol, state, action, a, r, yyval = {}, p, len, newState, expected;
      while (true) {
        state = stack[stack.length - 1];
        if (this.defaultActions[state]) {
          action = this.defaultActions[state];
        } else {
          if (symbol === null || typeof symbol == "undefined") {
            symbol = lex();
          }
          action = table[state] && table[state][symbol];
        }
        if (typeof action === "undefined" || !action.length || !action[0]) {
          var errStr = "";
          expected = [];
          for (p in table[state]) {
            if (this.terminals_[p] && p > TERROR) {
              expected.push("'" + this.terminals_[p] + "'");
            }
          }
          if (lexer2.showPosition) {
            errStr = "Parse error on line " + (yylineno + 1) + ":\n" + lexer2.showPosition() + "\nExpecting " + expected.join(", ") + ", got '" + (this.terminals_[symbol] || symbol) + "'";
          } else {
            errStr = "Parse error on line " + (yylineno + 1) + ": Unexpected " + (symbol == EOF ? "end of input" : "'" + (this.terminals_[symbol] || symbol) + "'");
          }
          this.parseError(errStr, {
            text: lexer2.match,
            token: this.terminals_[symbol] || symbol,
            line: lexer2.yylineno,
            loc: yyloc,
            expected
          });
        }
        if (action[0] instanceof Array && action.length > 1) {
          throw new Error("Parse Error: multiple actions possible at state: " + state + ", token: " + symbol);
        }
        switch (action[0]) {
          case 1:
            stack.push(symbol);
            vstack.push(lexer2.yytext);
            lstack.push(lexer2.yylloc);
            stack.push(action[1]);
            symbol = null;
            if (!preErrorSymbol) {
              yyleng = lexer2.yyleng;
              yytext = lexer2.yytext;
              yylineno = lexer2.yylineno;
              yyloc = lexer2.yylloc;
              if (recovering > 0) {
                recovering--;
              }
            } else {
              symbol = preErrorSymbol;
              preErrorSymbol = null;
            }
            break;
          case 2:
            len = this.productions_[action[1]][1];
            yyval.$ = vstack[vstack.length - len];
            yyval._$ = {
              first_line: lstack[lstack.length - (len || 1)].first_line,
              last_line: lstack[lstack.length - 1].last_line,
              first_column: lstack[lstack.length - (len || 1)].first_column,
              last_column: lstack[lstack.length - 1].last_column
            };
            if (ranges) {
              yyval._$.range = [
                lstack[lstack.length - (len || 1)].range[0],
                lstack[lstack.length - 1].range[1]
              ];
            }
            r = this.performAction.apply(yyval, [
              yytext,
              yyleng,
              yylineno,
              sharedState.yy,
              action[1],
              vstack,
              lstack
            ].concat(args));
            if (typeof r !== "undefined") {
              return r;
            }
            if (len) {
              stack = stack.slice(0, -1 * len * 2);
              vstack = vstack.slice(0, -1 * len);
              lstack = lstack.slice(0, -1 * len);
            }
            stack.push(this.productions_[action[1]][0]);
            vstack.push(yyval.$);
            lstack.push(yyval._$);
            newState = table[stack[stack.length - 2]][stack[stack.length - 1]];
            stack.push(newState);
            break;
          case 3:
            return true;
        }
      }
      return true;
    }, "parse")
  };
  var lexer = /* @__PURE__ */ function() {
    var lexer2 = {
      EOF: 1,
      parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(input, yy) {
        this.yy = yy || this.yy || {};
        this._input = input;
        this._more = this._backtrack = this.done = false;
        this.yylineno = this.yyleng = 0;
        this.yytext = this.matched = this.match = "";
        this.conditionStack = ["INITIAL"];
        this.yylloc = {
          first_line: 1,
          first_column: 0,
          last_line: 1,
          last_column: 0
        };
        if (this.options.ranges) {
          this.yylloc.range = [0, 0];
        }
        this.offset = 0;
        return this;
      }, "setInput"),
      // consumes and returns one char from the input
      input: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var ch = this._input[0];
        this.yytext += ch;
        this.yyleng++;
        this.offset++;
        this.match += ch;
        this.matched += ch;
        var lines = ch.match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno++;
          this.yylloc.last_line++;
        } else {
          this.yylloc.last_column++;
        }
        if (this.options.ranges) {
          this.yylloc.range[1]++;
        }
        this._input = this._input.slice(1);
        return ch;
      }, "input"),
      // unshifts one char (or a string) into the input
      unput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(ch) {
        var len = ch.length;
        var lines = ch.split(/(?:\r\n?|\n)/g);
        this._input = ch + this._input;
        this.yytext = this.yytext.substr(0, this.yytext.length - len);
        this.offset -= len;
        var oldLines = this.match.split(/(?:\r\n?|\n)/g);
        this.match = this.match.substr(0, this.match.length - 1);
        this.matched = this.matched.substr(0, this.matched.length - 1);
        if (lines.length - 1) {
          this.yylineno -= lines.length - 1;
        }
        var r = this.yylloc.range;
        this.yylloc = {
          first_line: this.yylloc.first_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.first_column,
          last_column: lines ? (lines.length === oldLines.length ? this.yylloc.first_column : 0) + oldLines[oldLines.length - lines.length].length - lines[0].length : this.yylloc.first_column - len
        };
        if (this.options.ranges) {
          this.yylloc.range = [r[0], r[0] + this.yyleng - len];
        }
        this.yyleng = this.yytext.length;
        return this;
      }, "unput"),
      // When called from action, caches matched text and appends it on next action
      more: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        if (this.options.backtrack_lexer) {
          this._backtrack = true;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". You can only invoke reject() in the lexer when the lexer is of the backtracking persuasion (options.backtrack_lexer = true).\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
        return this;
      }, "reject"),
      // retain first n characters of the match
      less: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(match, indexed_rule) {
        var token, lines, backup;
        if (this.options.backtrack_lexer) {
          backup = {
            yylineno: this.yylineno,
            yylloc: {
              first_line: this.yylloc.first_line,
              last_line: this.last_line,
              first_column: this.yylloc.first_column,
              last_column: this.yylloc.last_column
            },
            yytext: this.yytext,
            match: this.match,
            matches: this.matches,
            matched: this.matched,
            yyleng: this.yyleng,
            offset: this.offset,
            _more: this._more,
            _input: this._input,
            yy: this.yy,
            conditionStack: this.conditionStack.slice(0),
            done: this.done
          };
          if (this.options.ranges) {
            backup.yylloc.range = this.yylloc.range.slice(0);
          }
        }
        lines = match[0].match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno += lines.length;
        }
        this.yylloc = {
          first_line: this.yylloc.last_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.last_column,
          last_column: lines ? lines[lines.length - 1].length - lines[lines.length - 1].match(/\r?\n?/)[0].length : this.yylloc.last_column + match[0].length
        };
        this.yytext += match[0];
        this.match += match[0];
        this.matches = match;
        this.yyleng = this.yytext.length;
        if (this.options.ranges) {
          this.yylloc.range = [this.offset, this.offset += this.yyleng];
        }
        this._more = false;
        this._backtrack = false;
        this._input = this._input.slice(match[0].length);
        this.matched += match[0];
        token = this.performAction.call(this, this.yy, this, indexed_rule, this.conditionStack[this.conditionStack.length - 1]);
        if (this.done && this._input) {
          this.done = false;
        }
        if (token) {
          return token;
        } else if (this._backtrack) {
          for (var k in backup) {
            this[k] = backup[k];
          }
          return false;
        }
        return false;
      }, "test_match"),
      // return next match in input
      next: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        if (this.done) {
          return this.EOF;
        }
        if (!this._input) {
          this.done = true;
        }
        var token, match, tempMatch, index;
        if (!this._more) {
          this.yytext = "";
          this.match = "";
        }
        var rules = this._currentRules();
        for (var i = 0; i < rules.length; i++) {
          tempMatch = this._input.match(this.rules[rules[i]]);
          if (tempMatch && (!match || tempMatch[0].length > match[0].length)) {
            match = tempMatch;
            index = i;
            if (this.options.backtrack_lexer) {
              token = this.test_match(tempMatch, rules[i]);
              if (token !== false) {
                return token;
              } else if (this._backtrack) {
                match = false;
                continue;
              } else {
                return false;
              }
            } else if (!this.options.flex) {
              break;
            }
          }
        }
        if (match) {
          token = this.test_match(match, rules[index]);
          if (token !== false) {
            return token;
          }
          return false;
        }
        if (this._input === "") {
          return this.EOF;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". Unrecognized text.\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
      }, "next"),
      // return next match that has a token
      lex: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            this.begin("acc_title");
            return 24;
            break;
          case 1:
            this.popState();
            return "acc_title_value";
            break;
          case 2:
            this.begin("acc_descr");
            return 26;
            break;
          case 3:
            this.popState();
            return "acc_descr_value";
            break;
          case 4:
            this.begin("acc_descr_multiline");
            break;
          case 5:
            this.popState();
            break;
          case 6:
            return "acc_descr_multiline_value";
            break;
          case 7:
            return 33;
            break;
          case 8:
            return 34;
            break;
          case 9:
            return 35;
            break;
          case 10:
            return 36;
            break;
          case 11:
            return 10;
            break;
          case 12:
            break;
          case 13:
            return 8;
            break;
          case 14:
            return 50;
            break;
          case 15:
            return 70;
            break;
          case 16:
            return 4;
            break;
          case 17:
            this.begin("block");
            return 17;
            break;
          case 18:
            return 49;
            break;
          case 19:
            return 49;
            break;
          case 20:
            return 42;
            break;
          case 21:
            return 15;
            break;
          case 22:
            return 13;
            break;
          case 23:
            break;
          case 24:
            return 59;
            break;
          case 25:
            return 56;
            break;
          case 26:
            return 56;
            break;
          case 27:
            return 60;
            break;
          case 28:
            break;
          case 29:
            this.popState();
            return 19;
            break;
          case 30:
            return yy_.yytext[0];
            break;
          case 31:
            return 20;
            break;
          case 32:
            return 21;
            break;
          case 33:
            this.begin("style");
            return 44;
            break;
          case 34:
            this.popState();
            return 10;
            break;
          case 35:
            break;
          case 36:
            return 13;
            break;
          case 37:
            return 42;
            break;
          case 38:
            return 49;
            break;
          case 39:
            this.begin("style");
            return 37;
            break;
          case 40:
            return 43;
            break;
          case 41:
            return 63;
            break;
          case 42:
            return 65;
            break;
          case 43:
            return 65;
            break;
          case 44:
            return 65;
            break;
          case 45:
            return 63;
            break;
          case 46:
            return 63;
            break;
          case 47:
            return 64;
            break;
          case 48:
            return 64;
            break;
          case 49:
            return 64;
            break;
          case 50:
            return 64;
            break;
          case 51:
            return 64;
            break;
          case 52:
            return 65;
            break;
          case 53:
            return 64;
            break;
          case 54:
            return 65;
            break;
          case 55:
            return 66;
            break;
          case 56:
            return 66;
            break;
          case 57:
            return 66;
            break;
          case 58:
            return 66;
            break;
          case 59:
            return 63;
            break;
          case 60:
            return 64;
            break;
          case 61:
            return 65;
            break;
          case 62:
            return 67;
            break;
          case 63:
            return 68;
            break;
          case 64:
            return 69;
            break;
          case 65:
            return 69;
            break;
          case 66:
            return 68;
            break;
          case 67:
            return 68;
            break;
          case 68:
            return 68;
            break;
          case 69:
            return 41;
            break;
          case 70:
            return 47;
            break;
          case 71:
            return 40;
            break;
          case 72:
            return 48;
            break;
          case 73:
            return yy_.yytext[0];
            break;
          case 74:
            return 6;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:.*direction\s+TB[^\n]*)/i, /^(?:.*direction\s+BT[^\n]*)/i, /^(?:.*direction\s+RL[^\n]*)/i, /^(?:.*direction\s+LR[^\n]*)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:[\s]+)/i, /^(?:"[^"%\r\n\v\b\\]+")/i, /^(?:"[^"]*")/i, /^(?:erDiagram\b)/i, /^(?:\{)/i, /^(?:#)/i, /^(?:#)/i, /^(?:,)/i, /^(?::::)/i, /^(?::)/i, /^(?:\s+)/i, /^(?:\b((?:PK)|(?:FK)|(?:UK))\b)/i, /^(?:([^\s]*)[~].*[~]([^\s]*))/i, /^(?:([\*A-Za-z_\u00C0-\uFFFF][A-Za-z0-9\-\_\[\]\(\)\u00C0-\uFFFF\*]*))/i, /^(?:"[^"]*")/i, /^(?:[\n]+)/i, /^(?:\})/i, /^(?:.)/i, /^(?:\[)/i, /^(?:\])/i, /^(?:style\b)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?::)/i, /^(?:,)/i, /^(?:#)/i, /^(?:classDef\b)/i, /^(?:class\b)/i, /^(?:one or zero\b)/i, /^(?:one or more\b)/i, /^(?:one or many\b)/i, /^(?:1\+)/i, /^(?:\|o\b)/i, /^(?:zero or one\b)/i, /^(?:zero or more\b)/i, /^(?:zero or many\b)/i, /^(?:0\+)/i, /^(?:\}o\b)/i, /^(?:many\(0\))/i, /^(?:many\(1\))/i, /^(?:many\b)/i, /^(?:\}\|)/i, /^(?:one\b)/i, /^(?:only one\b)/i, /^(?:1\b)/i, /^(?:\|\|)/i, /^(?:o\|)/i, /^(?:o\{)/i, /^(?:\|\{)/i, /^(?:\s*u\b)/i, /^(?:\.\.)/i, /^(?:--)/i, /^(?:to\b)/i, /^(?:optionally to\b)/i, /^(?:\.-)/i, /^(?:-\.)/i, /^(?:([^\x00-\x7F]|\w|-|\*)+)/i, /^(?:;)/i, /^(?:([^\x00-\x7F]|\w|-|\*)+)/i, /^(?:[0-9])/i, /^(?:.)/i, /^(?:$)/i],
      conditions: { "style": { "rules": [34, 35, 36, 37, 38, 69, 70], "inclusive": false }, "acc_descr_multiline": { "rules": [5, 6], "inclusive": false }, "acc_descr": { "rules": [3], "inclusive": false }, "acc_title": { "rules": [1], "inclusive": false }, "block": { "rules": [23, 24, 25, 26, 27, 28, 29, 30], "inclusive": false }, "INITIAL": { "rules": [0, 2, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
var erDiagram_default = parser;

// src/diagrams/er/erDb.ts
var ErDB = class {
  constructor() {
    this.entities = /* @__PURE__ */ new Map();
    this.relationships = [];
    this.classes = /* @__PURE__ */ new Map();
    this.direction = "TB";
    this.Cardinality = {
      ZERO_OR_ONE: "ZERO_OR_ONE",
      ZERO_OR_MORE: "ZERO_OR_MORE",
      ONE_OR_MORE: "ONE_OR_MORE",
      ONLY_ONE: "ONLY_ONE",
      MD_PARENT: "MD_PARENT"
    };
    this.Identification = {
      NON_IDENTIFYING: "NON_IDENTIFYING",
      IDENTIFYING: "IDENTIFYING"
    };
    this.setAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .setAccTitle */ .GN;
    this.getAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getAccTitle */ .eu;
    this.setAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .setAccDescription */ .U$;
    this.getAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getAccDescription */ .Mx;
    this.setDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .setDiagramTitle */ .g2;
    this.getDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getDiagramTitle */ .Kr;
    this.getConfig = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(() => (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getConfig2 */ .nV)().er, "getConfig");
    this.clear();
    this.addEntity = this.addEntity.bind(this);
    this.addAttributes = this.addAttributes.bind(this);
    this.addRelationship = this.addRelationship.bind(this);
    this.setDirection = this.setDirection.bind(this);
    this.addCssStyles = this.addCssStyles.bind(this);
    this.addClass = this.addClass.bind(this);
    this.setClass = this.setClass.bind(this);
    this.setAccTitle = this.setAccTitle.bind(this);
    this.setAccDescription = this.setAccDescription.bind(this);
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(this, "ErDB");
  }
  /**
   * Add entity
   * @param name - The name of the entity
   * @param alias - The alias of the entity
   */
  addEntity(name, alias = "") {
    if (!this.entities.has(name)) {
      this.entities.set(name, {
        id: `entity-${name}-${this.entities.size}`,
        label: name,
        attributes: [],
        alias,
        shape: "erBox",
        look: (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getConfig2 */ .nV)().look ?? "default",
        cssClasses: "default",
        cssStyles: []
      });
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.info("Added new entity :", name);
    } else if (!this.entities.get(name)?.alias && alias) {
      this.entities.get(name).alias = alias;
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.info(`Add alias '${alias}' to entity '${name}'`);
    }
    return this.entities.get(name);
  }
  getEntity(name) {
    return this.entities.get(name);
  }
  getEntities() {
    return this.entities;
  }
  getClasses() {
    return this.classes;
  }
  addAttributes(entityName, attribs) {
    const entity = this.addEntity(entityName);
    let i;
    for (i = attribs.length - 1; i >= 0; i--) {
      if (!attribs[i].keys) {
        attribs[i].keys = [];
      }
      if (!attribs[i].comment) {
        attribs[i].comment = "";
      }
      entity.attributes.push(attribs[i]);
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.debug("Added attribute ", attribs[i].name);
    }
  }
  /**
   * Add a relationship
   *
   * @param entA - The first entity in the relationship
   * @param rolA - The role played by the first entity in relation to the second
   * @param entB - The second entity in the relationship
   * @param rSpec - The details of the relationship between the two entities
   */
  addRelationship(entA, rolA, entB, rSpec) {
    const entityA = this.entities.get(entA);
    const entityB = this.entities.get(entB);
    if (!entityA || !entityB) {
      return;
    }
    const rel = {
      entityA: entityA.id,
      roleA: rolA,
      entityB: entityB.id,
      relSpec: rSpec
    };
    this.relationships.push(rel);
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.debug("Added new relationship :", rel);
  }
  getRelationships() {
    return this.relationships;
  }
  getDirection() {
    return this.direction;
  }
  setDirection(dir) {
    this.direction = dir;
  }
  getCompiledStyles(classDefs) {
    let compiledStyles = [];
    for (const customClass of classDefs) {
      const cssClass = this.classes.get(customClass);
      if (cssClass?.styles) {
        compiledStyles = [...compiledStyles, ...cssClass.styles ?? []].map((s) => s.trim());
      }
      if (cssClass?.textStyles) {
        compiledStyles = [...compiledStyles, ...cssClass.textStyles ?? []].map((s) => s.trim());
      }
    }
    return compiledStyles;
  }
  addCssStyles(ids, styles) {
    for (const id of ids) {
      const entity = this.entities.get(id);
      if (!styles || !entity) {
        return;
      }
      for (const style of styles) {
        entity.cssStyles.push(style);
      }
    }
  }
  addClass(ids, style) {
    ids.forEach((id) => {
      let classNode = this.classes.get(id);
      if (classNode === void 0) {
        classNode = { id, styles: [], textStyles: [] };
        this.classes.set(id, classNode);
      }
      if (style) {
        style.forEach(function(s) {
          if (/color/.exec(s)) {
            const newStyle = s.replace("fill", "bgFill");
            classNode.textStyles.push(newStyle);
          }
          classNode.styles.push(s);
        });
      }
    });
  }
  setClass(ids, classNames) {
    for (const id of ids) {
      const entity = this.entities.get(id);
      if (entity) {
        for (const className of classNames) {
          entity.cssClasses += " " + className;
        }
      }
    }
  }
  clear() {
    this.entities = /* @__PURE__ */ new Map();
    this.classes = /* @__PURE__ */ new Map();
    this.relationships = [];
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .clear */ .ZH)();
  }
  getData() {
    const nodes = [];
    const edges = [];
    const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getConfig2 */ .nV)();
    for (const entityKey of this.entities.keys()) {
      const entityNode = this.entities.get(entityKey);
      if (entityNode) {
        entityNode.cssCompiledStyles = this.getCompiledStyles(entityNode.cssClasses.split(" "));
        nodes.push(entityNode);
      }
    }
    let count = 0;
    for (const relationship of this.relationships) {
      const edge = {
        id: (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getEdgeId */ .Ln)(relationship.entityA, relationship.entityB, {
          prefix: "id",
          counter: count++
        }),
        type: "normal",
        curve: "basis",
        start: relationship.entityA,
        end: relationship.entityB,
        label: relationship.roleA,
        labelpos: "c",
        thickness: "normal",
        classes: "relationshipLine",
        arrowTypeStart: relationship.relSpec.cardB.toLowerCase(),
        arrowTypeEnd: relationship.relSpec.cardA.toLowerCase(),
        pattern: relationship.relSpec.relType == "IDENTIFYING" ? "solid" : "dashed",
        look: config.look
      };
      edges.push(edge);
    }
    return { nodes, edges, other: {}, config, direction: "TB" };
  }
};

// src/diagrams/er/erRenderer-unified.ts
var erRenderer_unified_exports = {};
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__export */ .r2)(erRenderer_unified_exports, {
  draw: () => draw
});

var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(async function(text, id, _version, diag) {
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.info("REF0:");
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.info("Drawing er diagram (unified)", id);
  const { securityLevel, er: conf, layout } = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .getConfig2 */ .nV)();
  const data4Layout = diag.db.getData();
  const svg = (0,_chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getDiagramElement */ .q)(id, securityLevel);
  data4Layout.type = diag.type;
  data4Layout.layoutAlgorithm = (0,_chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getRegisteredLayoutAlgorithm */ ._b)(layout);
  data4Layout.config.flowchart.nodeSpacing = conf?.nodeSpacing || 140;
  data4Layout.config.flowchart.rankSpacing = conf?.rankSpacing || 80;
  data4Layout.direction = diag.db.getDirection();
  data4Layout.markers = ["only_one", "zero_or_one", "one_or_more", "zero_or_more"];
  data4Layout.diagramId = id;
  await (0,_chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__/* .render */ .sY)(data4Layout, svg);
  if (data4Layout.layoutAlgorithm === "elk") {
    svg.select(".edges").lower();
  }
  const backgroundNodes = svg.selectAll('[id*="-background"]');
  if (Array.from(backgroundNodes).length > 0) {
    backgroundNodes.each(function() {
      const backgroundNode = (0,d3__WEBPACK_IMPORTED_MODULE_10__/* .select */ .Ys)(this);
      const backgroundId = backgroundNode.attr("id");
      const nonBackgroundId = backgroundId.replace("-background", "");
      const nonBackgroundNode = svg.select(`#${CSS.escape(nonBackgroundId)}`);
      if (!nonBackgroundNode.empty()) {
        const transform = nonBackgroundNode.attr("transform");
        backgroundNode.attr("transform", transform);
      }
    });
  }
  const padding = 8;
  _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_8__/* .utils_default */ .w8.insertTitle(
    svg,
    "erDiagramTitleText",
    conf?.titleTopMargin ?? 25,
    diag.db.getDiagramTitle()
  );
  (0,_chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__/* .setupViewPortForSVG */ .j)(svg, padding, "erDiagram", conf?.useMaxWidth ?? true);
}, "draw");

// src/diagrams/er/styles.ts

var fade = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((color, opacity) => {
  const channel2 = khroma__WEBPACK_IMPORTED_MODULE_11__/* ["default"] */ .Z;
  const r = channel2(color, "r");
  const g = channel2(color, "g");
  const b = channel2(color, "b");
  return khroma__WEBPACK_IMPORTED_MODULE_12__/* ["default"] */ .Z(r, g, b, opacity);
}, "fade");
var getStyles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((options) => `
  .entityBox {
    fill: ${options.mainBkg};
    stroke: ${options.nodeBorder};
  }

  .relationshipLabelBox {
    fill: ${options.tertiaryColor};
    opacity: 0.7;
    background-color: ${options.tertiaryColor};
      rect {
        opacity: 0.5;
      }
  }

  .labelBkg {
    background-color: ${fade(options.tertiaryColor, 0.5)};
  }

  .edgeLabel .label {
    fill: ${options.nodeBorder};
    font-size: 14px;
  }

  .label {
    font-family: ${options.fontFamily};
    color: ${options.nodeTextColor || options.textColor};
  }

  .edge-pattern-dashed {
    stroke-dasharray: 8,8;
  }

  .node rect,
  .node circle,
  .node ellipse,
  .node polygon
  {
    fill: ${options.mainBkg};
    stroke: ${options.nodeBorder};
    stroke-width: 1px;
  }

  .relationshipLine {
    stroke: ${options.lineColor};
    stroke-width: 1;
    fill: none;
  }

  .marker {
    fill: none !important;
    stroke: ${options.lineColor} !important;
    stroke-width: 1;
  }
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/er/erDiagram.ts
var diagram = {
  parser: erDiagram_default,
  get db() {
    return new ErDB();
  },
  renderer: erRenderer_unified_exports,
  styles: styles_default
};



/***/ })

}]);
//# sourceMappingURL=7488.4d8124f72a0f10256f44.js.map?v=4d8124f72a0f10256f44