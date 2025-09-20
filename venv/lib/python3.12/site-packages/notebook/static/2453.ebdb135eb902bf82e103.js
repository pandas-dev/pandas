"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2453],{

/***/ 2520:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   AD: () => (/* binding */ getTextObj),
/* harmony export */   AE: () => (/* binding */ drawImage),
/* harmony export */   Mu: () => (/* binding */ drawRect),
/* harmony export */   O: () => (/* binding */ drawBackgroundRect),
/* harmony export */   kc: () => (/* binding */ getNoteRect),
/* harmony export */   rB: () => (/* binding */ drawEmbeddedImage),
/* harmony export */   yU: () => (/* binding */ drawText)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(7608);


// src/diagrams/common/svgDrawCommon.ts

var drawRect = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((element, rectData) => {
  const rectElement = element.append("rect");
  rectElement.attr("x", rectData.x);
  rectElement.attr("y", rectData.y);
  rectElement.attr("fill", rectData.fill);
  rectElement.attr("stroke", rectData.stroke);
  rectElement.attr("width", rectData.width);
  rectElement.attr("height", rectData.height);
  if (rectData.name) {
    rectElement.attr("name", rectData.name);
  }
  if (rectData.rx) {
    rectElement.attr("rx", rectData.rx);
  }
  if (rectData.ry) {
    rectElement.attr("ry", rectData.ry);
  }
  if (rectData.attrs !== void 0) {
    for (const attrKey in rectData.attrs) {
      rectElement.attr(attrKey, rectData.attrs[attrKey]);
    }
  }
  if (rectData.class) {
    rectElement.attr("class", rectData.class);
  }
  return rectElement;
}, "drawRect");
var drawBackgroundRect = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((element, bounds) => {
  const rectData = {
    x: bounds.startx,
    y: bounds.starty,
    width: bounds.stopx - bounds.startx,
    height: bounds.stopy - bounds.starty,
    fill: bounds.fill,
    stroke: bounds.stroke,
    class: "rect"
  };
  const rectElement = drawRect(element, rectData);
  rectElement.lower();
}, "drawBackgroundRect");
var drawText = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((element, textData) => {
  const nText = textData.text.replace(_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .lineBreakRegex */ .Vw, " ");
  const textElem = element.append("text");
  textElem.attr("x", textData.x);
  textElem.attr("y", textData.y);
  textElem.attr("class", "legend");
  textElem.style("text-anchor", textData.anchor);
  if (textData.class) {
    textElem.attr("class", textData.class);
  }
  const tspan = textElem.append("tspan");
  tspan.attr("x", textData.x + textData.textMargin * 2);
  tspan.text(nText);
  return textElem;
}, "drawText");
var drawImage = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((elem, x, y, link) => {
  const imageElement = elem.append("image");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_1__/* .sanitizeUrl */ .N)(link);
  imageElement.attr("xlink:href", sanitizedLink);
}, "drawImage");
var drawEmbeddedImage = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((element, x, y, link) => {
  const imageElement = element.append("use");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_1__/* .sanitizeUrl */ .N)(link);
  imageElement.attr("xlink:href", `#${sanitizedLink}`);
}, "drawEmbeddedImage");
var getNoteRect = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => {
  const noteRectData = {
    x: 0,
    y: 0,
    width: 100,
    height: 100,
    fill: "#EDF2AE",
    stroke: "#666",
    anchor: "start",
    rx: 0,
    ry: 0
  };
  return noteRectData;
}, "getNoteRect");
var getTextObj = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => {
  const testObject = {
    x: 0,
    y: 0,
    width: 100,
    height: 100,
    "text-anchor": "start",
    style: "#666",
    textMargin: 0,
    rx: 0,
    ry: 0,
    tspan: true
  };
  return testObject;
}, "getTextObj");




/***/ }),

/***/ 82002:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ ImperativeState)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);


// src/utils/imperativeState.ts
var ImperativeState = class {
  /**
   * @param init - Function that creates the default state.
   */
  constructor(init) {
    this.init = init;
    this.records = this.init();
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "ImperativeState");
  }
  reset() {
    this.records = this.init();
  }
};




/***/ }),

/***/ 62453:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(2520);
/* harmony import */ var _chunk_XZIHB7SX_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(82002);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(42626);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(83619);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(7608);





// src/diagrams/sequence/parser/sequenceDiagram.jison
var parser = function() {
  var o = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 2], $V1 = [1, 3], $V2 = [1, 4], $V3 = [2, 4], $V4 = [1, 9], $V5 = [1, 11], $V6 = [1, 13], $V7 = [1, 14], $V8 = [1, 16], $V9 = [1, 17], $Va = [1, 18], $Vb = [1, 24], $Vc = [1, 25], $Vd = [1, 26], $Ve = [1, 27], $Vf = [1, 28], $Vg = [1, 29], $Vh = [1, 30], $Vi = [1, 31], $Vj = [1, 32], $Vk = [1, 33], $Vl = [1, 34], $Vm = [1, 35], $Vn = [1, 36], $Vo = [1, 37], $Vp = [1, 38], $Vq = [1, 39], $Vr = [1, 41], $Vs = [1, 42], $Vt = [1, 43], $Vu = [1, 44], $Vv = [1, 45], $Vw = [1, 46], $Vx = [1, 4, 5, 13, 14, 16, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 47, 48, 49, 50, 52, 53, 54, 59, 60, 61, 62, 70], $Vy = [4, 5, 16, 50, 52, 53], $Vz = [4, 5, 13, 14, 16, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 50, 52, 53, 54, 59, 60, 61, 62, 70], $VA = [4, 5, 13, 14, 16, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 49, 50, 52, 53, 54, 59, 60, 61, 62, 70], $VB = [4, 5, 13, 14, 16, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 48, 50, 52, 53, 54, 59, 60, 61, 62, 70], $VC = [4, 5, 13, 14, 16, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 47, 50, 52, 53, 54, 59, 60, 61, 62, 70], $VD = [68, 69, 70], $VE = [1, 122];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "SPACE": 4, "NEWLINE": 5, "SD": 6, "document": 7, "line": 8, "statement": 9, "box_section": 10, "box_line": 11, "participant_statement": 12, "create": 13, "box": 14, "restOfLine": 15, "end": 16, "signal": 17, "autonumber": 18, "NUM": 19, "off": 20, "activate": 21, "actor": 22, "deactivate": 23, "note_statement": 24, "links_statement": 25, "link_statement": 26, "properties_statement": 27, "details_statement": 28, "title": 29, "legacy_title": 30, "acc_title": 31, "acc_title_value": 32, "acc_descr": 33, "acc_descr_value": 34, "acc_descr_multiline_value": 35, "loop": 36, "rect": 37, "opt": 38, "alt": 39, "else_sections": 40, "par": 41, "par_sections": 42, "par_over": 43, "critical": 44, "option_sections": 45, "break": 46, "option": 47, "and": 48, "else": 49, "participant": 50, "AS": 51, "participant_actor": 52, "destroy": 53, "note": 54, "placement": 55, "text2": 56, "over": 57, "actor_pair": 58, "links": 59, "link": 60, "properties": 61, "details": 62, "spaceList": 63, ",": 64, "left_of": 65, "right_of": 66, "signaltype": 67, "+": 68, "-": 69, "ACTOR": 70, "SOLID_OPEN_ARROW": 71, "DOTTED_OPEN_ARROW": 72, "SOLID_ARROW": 73, "BIDIRECTIONAL_SOLID_ARROW": 74, "DOTTED_ARROW": 75, "BIDIRECTIONAL_DOTTED_ARROW": 76, "SOLID_CROSS": 77, "DOTTED_CROSS": 78, "SOLID_POINT": 79, "DOTTED_POINT": 80, "TXT": 81, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "SPACE", 5: "NEWLINE", 6: "SD", 13: "create", 14: "box", 15: "restOfLine", 16: "end", 18: "autonumber", 19: "NUM", 20: "off", 21: "activate", 23: "deactivate", 29: "title", 30: "legacy_title", 31: "acc_title", 32: "acc_title_value", 33: "acc_descr", 34: "acc_descr_value", 35: "acc_descr_multiline_value", 36: "loop", 37: "rect", 38: "opt", 39: "alt", 41: "par", 43: "par_over", 44: "critical", 46: "break", 47: "option", 48: "and", 49: "else", 50: "participant", 51: "AS", 52: "participant_actor", 53: "destroy", 54: "note", 57: "over", 59: "links", 60: "link", 61: "properties", 62: "details", 64: ",", 65: "left_of", 66: "right_of", 68: "+", 69: "-", 70: "ACTOR", 71: "SOLID_OPEN_ARROW", 72: "DOTTED_OPEN_ARROW", 73: "SOLID_ARROW", 74: "BIDIRECTIONAL_SOLID_ARROW", 75: "DOTTED_ARROW", 76: "BIDIRECTIONAL_DOTTED_ARROW", 77: "SOLID_CROSS", 78: "DOTTED_CROSS", 79: "SOLID_POINT", 80: "DOTTED_POINT", 81: "TXT" },
    productions_: [0, [3, 2], [3, 2], [3, 2], [7, 0], [7, 2], [8, 2], [8, 1], [8, 1], [10, 0], [10, 2], [11, 2], [11, 1], [11, 1], [9, 1], [9, 2], [9, 4], [9, 2], [9, 4], [9, 3], [9, 3], [9, 2], [9, 3], [9, 3], [9, 2], [9, 2], [9, 2], [9, 2], [9, 2], [9, 1], [9, 1], [9, 2], [9, 2], [9, 1], [9, 4], [9, 4], [9, 4], [9, 4], [9, 4], [9, 4], [9, 4], [9, 4], [45, 1], [45, 4], [42, 1], [42, 4], [40, 1], [40, 4], [12, 5], [12, 3], [12, 5], [12, 3], [12, 3], [24, 4], [24, 4], [25, 3], [26, 3], [27, 3], [28, 3], [63, 2], [63, 1], [58, 3], [58, 1], [55, 1], [55, 1], [17, 5], [17, 5], [17, 4], [22, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [67, 1], [56, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 3:
          yy.apply($$[$0]);
          return $$[$0];
          break;
        case 4:
        case 9:
          this.$ = [];
          break;
        case 5:
        case 10:
          $$[$0 - 1].push($$[$0]);
          this.$ = $$[$0 - 1];
          break;
        case 6:
        case 7:
        case 11:
        case 12:
          this.$ = $$[$0];
          break;
        case 8:
        case 13:
          this.$ = [];
          break;
        case 15:
          $$[$0].type = "createParticipant";
          this.$ = $$[$0];
          break;
        case 16:
          $$[$0 - 1].unshift({ type: "boxStart", boxData: yy.parseBoxData($$[$0 - 2]) });
          $$[$0 - 1].push({ type: "boxEnd", boxText: $$[$0 - 2] });
          this.$ = $$[$0 - 1];
          break;
        case 18:
          this.$ = { type: "sequenceIndex", sequenceIndex: Number($$[$0 - 2]), sequenceIndexStep: Number($$[$0 - 1]), sequenceVisible: true, signalType: yy.LINETYPE.AUTONUMBER };
          break;
        case 19:
          this.$ = { type: "sequenceIndex", sequenceIndex: Number($$[$0 - 1]), sequenceIndexStep: 1, sequenceVisible: true, signalType: yy.LINETYPE.AUTONUMBER };
          break;
        case 20:
          this.$ = { type: "sequenceIndex", sequenceVisible: false, signalType: yy.LINETYPE.AUTONUMBER };
          break;
        case 21:
          this.$ = { type: "sequenceIndex", sequenceVisible: true, signalType: yy.LINETYPE.AUTONUMBER };
          break;
        case 22:
          this.$ = { type: "activeStart", signalType: yy.LINETYPE.ACTIVE_START, actor: $$[$0 - 1].actor };
          break;
        case 23:
          this.$ = { type: "activeEnd", signalType: yy.LINETYPE.ACTIVE_END, actor: $$[$0 - 1].actor };
          break;
        case 29:
          yy.setDiagramTitle($$[$0].substring(6));
          this.$ = $$[$0].substring(6);
          break;
        case 30:
          yy.setDiagramTitle($$[$0].substring(7));
          this.$ = $$[$0].substring(7);
          break;
        case 31:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 32:
        case 33:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 34:
          $$[$0 - 1].unshift({ type: "loopStart", loopText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.LOOP_START });
          $$[$0 - 1].push({ type: "loopEnd", loopText: $$[$0 - 2], signalType: yy.LINETYPE.LOOP_END });
          this.$ = $$[$0 - 1];
          break;
        case 35:
          $$[$0 - 1].unshift({ type: "rectStart", color: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.RECT_START });
          $$[$0 - 1].push({ type: "rectEnd", color: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.RECT_END });
          this.$ = $$[$0 - 1];
          break;
        case 36:
          $$[$0 - 1].unshift({ type: "optStart", optText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.OPT_START });
          $$[$0 - 1].push({ type: "optEnd", optText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.OPT_END });
          this.$ = $$[$0 - 1];
          break;
        case 37:
          $$[$0 - 1].unshift({ type: "altStart", altText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.ALT_START });
          $$[$0 - 1].push({ type: "altEnd", signalType: yy.LINETYPE.ALT_END });
          this.$ = $$[$0 - 1];
          break;
        case 38:
          $$[$0 - 1].unshift({ type: "parStart", parText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.PAR_START });
          $$[$0 - 1].push({ type: "parEnd", signalType: yy.LINETYPE.PAR_END });
          this.$ = $$[$0 - 1];
          break;
        case 39:
          $$[$0 - 1].unshift({ type: "parStart", parText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.PAR_OVER_START });
          $$[$0 - 1].push({ type: "parEnd", signalType: yy.LINETYPE.PAR_END });
          this.$ = $$[$0 - 1];
          break;
        case 40:
          $$[$0 - 1].unshift({ type: "criticalStart", criticalText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.CRITICAL_START });
          $$[$0 - 1].push({ type: "criticalEnd", signalType: yy.LINETYPE.CRITICAL_END });
          this.$ = $$[$0 - 1];
          break;
        case 41:
          $$[$0 - 1].unshift({ type: "breakStart", breakText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.BREAK_START });
          $$[$0 - 1].push({ type: "breakEnd", optText: yy.parseMessage($$[$0 - 2]), signalType: yy.LINETYPE.BREAK_END });
          this.$ = $$[$0 - 1];
          break;
        case 43:
          this.$ = $$[$0 - 3].concat([{ type: "option", optionText: yy.parseMessage($$[$0 - 1]), signalType: yy.LINETYPE.CRITICAL_OPTION }, $$[$0]]);
          break;
        case 45:
          this.$ = $$[$0 - 3].concat([{ type: "and", parText: yy.parseMessage($$[$0 - 1]), signalType: yy.LINETYPE.PAR_AND }, $$[$0]]);
          break;
        case 47:
          this.$ = $$[$0 - 3].concat([{ type: "else", altText: yy.parseMessage($$[$0 - 1]), signalType: yy.LINETYPE.ALT_ELSE }, $$[$0]]);
          break;
        case 48:
          $$[$0 - 3].draw = "participant";
          $$[$0 - 3].type = "addParticipant";
          $$[$0 - 3].description = yy.parseMessage($$[$0 - 1]);
          this.$ = $$[$0 - 3];
          break;
        case 49:
          $$[$0 - 1].draw = "participant";
          $$[$0 - 1].type = "addParticipant";
          this.$ = $$[$0 - 1];
          break;
        case 50:
          $$[$0 - 3].draw = "actor";
          $$[$0 - 3].type = "addParticipant";
          $$[$0 - 3].description = yy.parseMessage($$[$0 - 1]);
          this.$ = $$[$0 - 3];
          break;
        case 51:
          $$[$0 - 1].draw = "actor";
          $$[$0 - 1].type = "addParticipant";
          this.$ = $$[$0 - 1];
          break;
        case 52:
          $$[$0 - 1].type = "destroyParticipant";
          this.$ = $$[$0 - 1];
          break;
        case 53:
          this.$ = [$$[$0 - 1], { type: "addNote", placement: $$[$0 - 2], actor: $$[$0 - 1].actor, text: $$[$0] }];
          break;
        case 54:
          $$[$0 - 2] = [].concat($$[$0 - 1], $$[$0 - 1]).slice(0, 2);
          $$[$0 - 2][0] = $$[$0 - 2][0].actor;
          $$[$0 - 2][1] = $$[$0 - 2][1].actor;
          this.$ = [$$[$0 - 1], { type: "addNote", placement: yy.PLACEMENT.OVER, actor: $$[$0 - 2].slice(0, 2), text: $$[$0] }];
          break;
        case 55:
          this.$ = [$$[$0 - 1], { type: "addLinks", actor: $$[$0 - 1].actor, text: $$[$0] }];
          break;
        case 56:
          this.$ = [$$[$0 - 1], { type: "addALink", actor: $$[$0 - 1].actor, text: $$[$0] }];
          break;
        case 57:
          this.$ = [$$[$0 - 1], { type: "addProperties", actor: $$[$0 - 1].actor, text: $$[$0] }];
          break;
        case 58:
          this.$ = [$$[$0 - 1], { type: "addDetails", actor: $$[$0 - 1].actor, text: $$[$0] }];
          break;
        case 61:
          this.$ = [$$[$0 - 2], $$[$0]];
          break;
        case 62:
          this.$ = $$[$0];
          break;
        case 63:
          this.$ = yy.PLACEMENT.LEFTOF;
          break;
        case 64:
          this.$ = yy.PLACEMENT.RIGHTOF;
          break;
        case 65:
          this.$ = [
            $$[$0 - 4],
            $$[$0 - 1],
            { type: "addMessage", from: $$[$0 - 4].actor, to: $$[$0 - 1].actor, signalType: $$[$0 - 3], msg: $$[$0], activate: true },
            { type: "activeStart", signalType: yy.LINETYPE.ACTIVE_START, actor: $$[$0 - 1].actor }
          ];
          break;
        case 66:
          this.$ = [
            $$[$0 - 4],
            $$[$0 - 1],
            { type: "addMessage", from: $$[$0 - 4].actor, to: $$[$0 - 1].actor, signalType: $$[$0 - 3], msg: $$[$0] },
            { type: "activeEnd", signalType: yy.LINETYPE.ACTIVE_END, actor: $$[$0 - 4].actor }
          ];
          break;
        case 67:
          this.$ = [$$[$0 - 3], $$[$0 - 1], { type: "addMessage", from: $$[$0 - 3].actor, to: $$[$0 - 1].actor, signalType: $$[$0 - 2], msg: $$[$0] }];
          break;
        case 68:
          this.$ = { type: "addParticipant", actor: $$[$0] };
          break;
        case 69:
          this.$ = yy.LINETYPE.SOLID_OPEN;
          break;
        case 70:
          this.$ = yy.LINETYPE.DOTTED_OPEN;
          break;
        case 71:
          this.$ = yy.LINETYPE.SOLID;
          break;
        case 72:
          this.$ = yy.LINETYPE.BIDIRECTIONAL_SOLID;
          break;
        case 73:
          this.$ = yy.LINETYPE.DOTTED;
          break;
        case 74:
          this.$ = yy.LINETYPE.BIDIRECTIONAL_DOTTED;
          break;
        case 75:
          this.$ = yy.LINETYPE.SOLID_CROSS;
          break;
        case 76:
          this.$ = yy.LINETYPE.DOTTED_CROSS;
          break;
        case 77:
          this.$ = yy.LINETYPE.SOLID_POINT;
          break;
        case 78:
          this.$ = yy.LINETYPE.DOTTED_POINT;
          break;
        case 79:
          this.$ = yy.parseMessage($$[$0].trim().substring(1));
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: $V0, 5: $V1, 6: $V2 }, { 1: [3] }, { 3: 5, 4: $V0, 5: $V1, 6: $V2 }, { 3: 6, 4: $V0, 5: $V1, 6: $V2 }, o([1, 4, 5, 13, 14, 18, 21, 23, 29, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43, 44, 46, 50, 52, 53, 54, 59, 60, 61, 62, 70], $V3, { 7: 7 }), { 1: [2, 1] }, { 1: [2, 2] }, { 1: [2, 3], 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, o($Vx, [2, 5]), { 9: 47, 12: 12, 13: $V6, 14: $V7, 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, o($Vx, [2, 7]), o($Vx, [2, 8]), o($Vx, [2, 14]), { 12: 48, 50: $Vo, 52: $Vp, 53: $Vq }, { 15: [1, 49] }, { 5: [1, 50] }, { 5: [1, 53], 19: [1, 51], 20: [1, 52] }, { 22: 54, 70: $Vw }, { 22: 55, 70: $Vw }, { 5: [1, 56] }, { 5: [1, 57] }, { 5: [1, 58] }, { 5: [1, 59] }, { 5: [1, 60] }, o($Vx, [2, 29]), o($Vx, [2, 30]), { 32: [1, 61] }, { 34: [1, 62] }, o($Vx, [2, 33]), { 15: [1, 63] }, { 15: [1, 64] }, { 15: [1, 65] }, { 15: [1, 66] }, { 15: [1, 67] }, { 15: [1, 68] }, { 15: [1, 69] }, { 15: [1, 70] }, { 22: 71, 70: $Vw }, { 22: 72, 70: $Vw }, { 22: 73, 70: $Vw }, { 67: 74, 71: [1, 75], 72: [1, 76], 73: [1, 77], 74: [1, 78], 75: [1, 79], 76: [1, 80], 77: [1, 81], 78: [1, 82], 79: [1, 83], 80: [1, 84] }, { 55: 85, 57: [1, 86], 65: [1, 87], 66: [1, 88] }, { 22: 89, 70: $Vw }, { 22: 90, 70: $Vw }, { 22: 91, 70: $Vw }, { 22: 92, 70: $Vw }, o([5, 51, 64, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81], [2, 68]), o($Vx, [2, 6]), o($Vx, [2, 15]), o($Vy, [2, 9], { 10: 93 }), o($Vx, [2, 17]), { 5: [1, 95], 19: [1, 94] }, { 5: [1, 96] }, o($Vx, [2, 21]), { 5: [1, 97] }, { 5: [1, 98] }, o($Vx, [2, 24]), o($Vx, [2, 25]), o($Vx, [2, 26]), o($Vx, [2, 27]), o($Vx, [2, 28]), o($Vx, [2, 31]), o($Vx, [2, 32]), o($Vz, $V3, { 7: 99 }), o($Vz, $V3, { 7: 100 }), o($Vz, $V3, { 7: 101 }), o($VA, $V3, { 40: 102, 7: 103 }), o($VB, $V3, { 42: 104, 7: 105 }), o($VB, $V3, { 7: 105, 42: 106 }), o($VC, $V3, { 45: 107, 7: 108 }), o($Vz, $V3, { 7: 109 }), { 5: [1, 111], 51: [1, 110] }, { 5: [1, 113], 51: [1, 112] }, { 5: [1, 114] }, { 22: 117, 68: [1, 115], 69: [1, 116], 70: $Vw }, o($VD, [2, 69]), o($VD, [2, 70]), o($VD, [2, 71]), o($VD, [2, 72]), o($VD, [2, 73]), o($VD, [2, 74]), o($VD, [2, 75]), o($VD, [2, 76]), o($VD, [2, 77]), o($VD, [2, 78]), { 22: 118, 70: $Vw }, { 22: 120, 58: 119, 70: $Vw }, { 70: [2, 63] }, { 70: [2, 64] }, { 56: 121, 81: $VE }, { 56: 123, 81: $VE }, { 56: 124, 81: $VE }, { 56: 125, 81: $VE }, { 4: [1, 128], 5: [1, 130], 11: 127, 12: 129, 16: [1, 126], 50: $Vo, 52: $Vp, 53: $Vq }, { 5: [1, 131] }, o($Vx, [2, 19]), o($Vx, [2, 20]), o($Vx, [2, 22]), o($Vx, [2, 23]), { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [1, 132], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [1, 133], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [1, 134], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 16: [1, 135] }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [2, 46], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 49: [1, 136], 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 16: [1, 137] }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [2, 44], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 48: [1, 138], 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 16: [1, 139] }, { 16: [1, 140] }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [2, 42], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 47: [1, 141], 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 4: $V4, 5: $V5, 8: 8, 9: 10, 12: 12, 13: $V6, 14: $V7, 16: [1, 142], 17: 15, 18: $V8, 21: $V9, 22: 40, 23: $Va, 24: 19, 25: 20, 26: 21, 27: 22, 28: 23, 29: $Vb, 30: $Vc, 31: $Vd, 33: $Ve, 35: $Vf, 36: $Vg, 37: $Vh, 38: $Vi, 39: $Vj, 41: $Vk, 43: $Vl, 44: $Vm, 46: $Vn, 50: $Vo, 52: $Vp, 53: $Vq, 54: $Vr, 59: $Vs, 60: $Vt, 61: $Vu, 62: $Vv, 70: $Vw }, { 15: [1, 143] }, o($Vx, [2, 49]), { 15: [1, 144] }, o($Vx, [2, 51]), o($Vx, [2, 52]), { 22: 145, 70: $Vw }, { 22: 146, 70: $Vw }, { 56: 147, 81: $VE }, { 56: 148, 81: $VE }, { 56: 149, 81: $VE }, { 64: [1, 150], 81: [2, 62] }, { 5: [2, 55] }, { 5: [2, 79] }, { 5: [2, 56] }, { 5: [2, 57] }, { 5: [2, 58] }, o($Vx, [2, 16]), o($Vy, [2, 10]), { 12: 151, 50: $Vo, 52: $Vp, 53: $Vq }, o($Vy, [2, 12]), o($Vy, [2, 13]), o($Vx, [2, 18]), o($Vx, [2, 34]), o($Vx, [2, 35]), o($Vx, [2, 36]), o($Vx, [2, 37]), { 15: [1, 152] }, o($Vx, [2, 38]), { 15: [1, 153] }, o($Vx, [2, 39]), o($Vx, [2, 40]), { 15: [1, 154] }, o($Vx, [2, 41]), { 5: [1, 155] }, { 5: [1, 156] }, { 56: 157, 81: $VE }, { 56: 158, 81: $VE }, { 5: [2, 67] }, { 5: [2, 53] }, { 5: [2, 54] }, { 22: 159, 70: $Vw }, o($Vy, [2, 11]), o($VA, $V3, { 7: 103, 40: 160 }), o($VB, $V3, { 7: 105, 42: 161 }), o($VC, $V3, { 7: 108, 45: 162 }), o($Vx, [2, 48]), o($Vx, [2, 50]), { 5: [2, 65] }, { 5: [2, 66] }, { 81: [2, 61] }, { 16: [2, 47] }, { 16: [2, 45] }, { 16: [2, 43] }],
    defaultActions: { 5: [2, 1], 6: [2, 2], 87: [2, 63], 88: [2, 64], 121: [2, 55], 122: [2, 79], 123: [2, 56], 124: [2, 57], 125: [2, 58], 147: [2, 67], 148: [2, 53], 149: [2, 54], 157: [2, 65], 158: [2, 66], 159: [2, 61], 160: [2, 47], 161: [2, 45], 162: [2, 43] },
    parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(lex, "lex");
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
      parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            return 5;
            break;
          case 1:
            break;
          case 2:
            break;
          case 3:
            break;
          case 4:
            break;
          case 5:
            break;
          case 6:
            return 19;
            break;
          case 7:
            this.begin("LINE");
            return 14;
            break;
          case 8:
            this.begin("ID");
            return 50;
            break;
          case 9:
            this.begin("ID");
            return 52;
            break;
          case 10:
            return 13;
            break;
          case 11:
            this.begin("ID");
            return 53;
            break;
          case 12:
            yy_.yytext = yy_.yytext.trim();
            this.begin("ALIAS");
            return 70;
            break;
          case 13:
            this.popState();
            this.popState();
            this.begin("LINE");
            return 51;
            break;
          case 14:
            this.popState();
            this.popState();
            return 5;
            break;
          case 15:
            this.begin("LINE");
            return 36;
            break;
          case 16:
            this.begin("LINE");
            return 37;
            break;
          case 17:
            this.begin("LINE");
            return 38;
            break;
          case 18:
            this.begin("LINE");
            return 39;
            break;
          case 19:
            this.begin("LINE");
            return 49;
            break;
          case 20:
            this.begin("LINE");
            return 41;
            break;
          case 21:
            this.begin("LINE");
            return 43;
            break;
          case 22:
            this.begin("LINE");
            return 48;
            break;
          case 23:
            this.begin("LINE");
            return 44;
            break;
          case 24:
            this.begin("LINE");
            return 47;
            break;
          case 25:
            this.begin("LINE");
            return 46;
            break;
          case 26:
            this.popState();
            return 15;
            break;
          case 27:
            return 16;
            break;
          case 28:
            return 65;
            break;
          case 29:
            return 66;
            break;
          case 30:
            return 59;
            break;
          case 31:
            return 60;
            break;
          case 32:
            return 61;
            break;
          case 33:
            return 62;
            break;
          case 34:
            return 57;
            break;
          case 35:
            return 54;
            break;
          case 36:
            this.begin("ID");
            return 21;
            break;
          case 37:
            this.begin("ID");
            return 23;
            break;
          case 38:
            return 29;
            break;
          case 39:
            return 30;
            break;
          case 40:
            this.begin("acc_title");
            return 31;
            break;
          case 41:
            this.popState();
            return "acc_title_value";
            break;
          case 42:
            this.begin("acc_descr");
            return 33;
            break;
          case 43:
            this.popState();
            return "acc_descr_value";
            break;
          case 44:
            this.begin("acc_descr_multiline");
            break;
          case 45:
            this.popState();
            break;
          case 46:
            return "acc_descr_multiline_value";
            break;
          case 47:
            return 6;
            break;
          case 48:
            return 18;
            break;
          case 49:
            return 20;
            break;
          case 50:
            return 64;
            break;
          case 51:
            return 5;
            break;
          case 52:
            yy_.yytext = yy_.yytext.trim();
            return 70;
            break;
          case 53:
            return 73;
            break;
          case 54:
            return 74;
            break;
          case 55:
            return 75;
            break;
          case 56:
            return 76;
            break;
          case 57:
            return 71;
            break;
          case 58:
            return 72;
            break;
          case 59:
            return 77;
            break;
          case 60:
            return 78;
            break;
          case 61:
            return 79;
            break;
          case 62:
            return 80;
            break;
          case 63:
            return 81;
            break;
          case 64:
            return 68;
            break;
          case 65:
            return 69;
            break;
          case 66:
            return 5;
            break;
          case 67:
            return "INVALID";
            break;
        }
      }, "anonymous"),
      rules: [/^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:((?!\n)\s)+)/i, /^(?:#[^\n]*)/i, /^(?:%(?!\{)[^\n]*)/i, /^(?:[^\}]%%[^\n]*)/i, /^(?:[0-9]+(?=[ \n]+))/i, /^(?:box\b)/i, /^(?:participant\b)/i, /^(?:actor\b)/i, /^(?:create\b)/i, /^(?:destroy\b)/i, /^(?:[^\<->\->:\n,;]+?([\-]*[^\<->\->:\n,;]+?)*?(?=((?!\n)\s)+as(?!\n)\s|[#\n;]|$))/i, /^(?:as\b)/i, /^(?:(?:))/i, /^(?:loop\b)/i, /^(?:rect\b)/i, /^(?:opt\b)/i, /^(?:alt\b)/i, /^(?:else\b)/i, /^(?:par\b)/i, /^(?:par_over\b)/i, /^(?:and\b)/i, /^(?:critical\b)/i, /^(?:option\b)/i, /^(?:break\b)/i, /^(?:(?:[:]?(?:no)?wrap)?[^#\n;]*)/i, /^(?:end\b)/i, /^(?:left of\b)/i, /^(?:right of\b)/i, /^(?:links\b)/i, /^(?:link\b)/i, /^(?:properties\b)/i, /^(?:details\b)/i, /^(?:over\b)/i, /^(?:note\b)/i, /^(?:activate\b)/i, /^(?:deactivate\b)/i, /^(?:title\s[^#\n;]+)/i, /^(?:title:\s[^#\n;]+)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:sequenceDiagram\b)/i, /^(?:autonumber\b)/i, /^(?:off\b)/i, /^(?:,)/i, /^(?:;)/i, /^(?:[^\+\<->\->:\n,;]+((?!(-x|--x|-\)|--\)))[\-]*[^\+\<->\->:\n,;]+)*)/i, /^(?:->>)/i, /^(?:<<->>)/i, /^(?:-->>)/i, /^(?:<<-->>)/i, /^(?:->)/i, /^(?:-->)/i, /^(?:-[x])/i, /^(?:--[x])/i, /^(?:-[\)])/i, /^(?:--[\)])/i, /^(?::(?:(?:no)?wrap)?[^#\n;]+)/i, /^(?:\+)/i, /^(?:-)/i, /^(?:$)/i, /^(?:.)/i],
      conditions: { "acc_descr_multiline": { "rules": [45, 46], "inclusive": false }, "acc_descr": { "rules": [43], "inclusive": false }, "acc_title": { "rules": [41], "inclusive": false }, "ID": { "rules": [2, 3, 12], "inclusive": false }, "ALIAS": { "rules": [2, 3, 13, 14], "inclusive": false }, "LINE": { "rules": [2, 3, 26], "inclusive": false }, "INITIAL": { "rules": [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 44, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
var sequenceDiagram_default = parser;

// src/diagrams/sequence/sequenceDb.ts
var LINETYPE = {
  SOLID: 0,
  DOTTED: 1,
  NOTE: 2,
  SOLID_CROSS: 3,
  DOTTED_CROSS: 4,
  SOLID_OPEN: 5,
  DOTTED_OPEN: 6,
  LOOP_START: 10,
  LOOP_END: 11,
  ALT_START: 12,
  ALT_ELSE: 13,
  ALT_END: 14,
  OPT_START: 15,
  OPT_END: 16,
  ACTIVE_START: 17,
  ACTIVE_END: 18,
  PAR_START: 19,
  PAR_AND: 20,
  PAR_END: 21,
  RECT_START: 22,
  RECT_END: 23,
  SOLID_POINT: 24,
  DOTTED_POINT: 25,
  AUTONUMBER: 26,
  CRITICAL_START: 27,
  CRITICAL_OPTION: 28,
  CRITICAL_END: 29,
  BREAK_START: 30,
  BREAK_END: 31,
  PAR_OVER_START: 32,
  BIDIRECTIONAL_SOLID: 33,
  BIDIRECTIONAL_DOTTED: 34
};
var ARROWTYPE = {
  FILLED: 0,
  OPEN: 1
};
var PLACEMENT = {
  LEFTOF: 0,
  RIGHTOF: 1,
  OVER: 2
};
var SequenceDB = class {
  constructor() {
    this.state = new _chunk_XZIHB7SX_mjs__WEBPACK_IMPORTED_MODULE_1__/* .ImperativeState */ .A(() => ({
      prevActor: void 0,
      actors: /* @__PURE__ */ new Map(),
      createdActors: /* @__PURE__ */ new Map(),
      destroyedActors: /* @__PURE__ */ new Map(),
      boxes: [],
      messages: [],
      notes: [],
      sequenceNumbersEnabled: false,
      wrapEnabled: void 0,
      currentBox: void 0,
      lastCreated: void 0,
      lastDestroyed: void 0
    }));
    this.setAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN;
    this.setAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$;
    this.setDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2;
    this.getAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu;
    this.getAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx;
    this.getDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr;
    this.apply = this.apply.bind(this);
    this.parseBoxData = this.parseBoxData.bind(this);
    this.parseMessage = this.parseMessage.bind(this);
    this.clear();
    this.setWrap((0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)().wrap);
    this.LINETYPE = LINETYPE;
    this.ARROWTYPE = ARROWTYPE;
    this.PLACEMENT = PLACEMENT;
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "SequenceDB");
  }
  addBox(data) {
    this.state.records.boxes.push({
      name: data.text,
      wrap: data.wrap ?? this.autoWrap(),
      fill: data.color,
      actorKeys: []
    });
    this.state.records.currentBox = this.state.records.boxes.slice(-1)[0];
  }
  addActor(id, name, description, type) {
    let assignedBox = this.state.records.currentBox;
    const old = this.state.records.actors.get(id);
    if (old) {
      if (this.state.records.currentBox && old.box && this.state.records.currentBox !== old.box) {
        throw new Error(
          `A same participant should only be defined in one Box: ${old.name} can't be in '${old.box.name}' and in '${this.state.records.currentBox.name}' at the same time.`
        );
      }
      assignedBox = old.box ? old.box : this.state.records.currentBox;
      old.box = assignedBox;
      if (old && name === old.name && description == null) {
        return;
      }
    }
    if (description?.text == null) {
      description = { text: name, type };
    }
    if (type == null || description.text == null) {
      description = { text: name, type };
    }
    this.state.records.actors.set(id, {
      box: assignedBox,
      name,
      description: description.text,
      wrap: description.wrap ?? this.autoWrap(),
      prevActor: this.state.records.prevActor,
      links: {},
      properties: {},
      actorCnt: null,
      rectData: null,
      type: type ?? "participant"
    });
    if (this.state.records.prevActor) {
      const prevActorInRecords = this.state.records.actors.get(this.state.records.prevActor);
      if (prevActorInRecords) {
        prevActorInRecords.nextActor = id;
      }
    }
    if (this.state.records.currentBox) {
      this.state.records.currentBox.actorKeys.push(id);
    }
    this.state.records.prevActor = id;
  }
  activationCount(part) {
    let i;
    let count = 0;
    if (!part) {
      return 0;
    }
    for (i = 0; i < this.state.records.messages.length; i++) {
      if (this.state.records.messages[i].type === this.LINETYPE.ACTIVE_START && this.state.records.messages[i].from === part) {
        count++;
      }
      if (this.state.records.messages[i].type === this.LINETYPE.ACTIVE_END && this.state.records.messages[i].from === part) {
        count--;
      }
    }
    return count;
  }
  addMessage(idFrom, idTo, message, answer) {
    this.state.records.messages.push({
      id: this.state.records.messages.length.toString(),
      from: idFrom,
      to: idTo,
      message: message.text,
      wrap: message.wrap ?? this.autoWrap(),
      answer
    });
  }
  addSignal(idFrom, idTo, message, messageType, activate = false) {
    if (messageType === this.LINETYPE.ACTIVE_END) {
      const cnt = this.activationCount(idFrom ?? "");
      if (cnt < 1) {
        const error = new Error("Trying to inactivate an inactive participant (" + idFrom + ")");
        error.hash = {
          text: "->>-",
          token: "->>-",
          line: "1",
          loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
          expected: ["'ACTIVE_PARTICIPANT'"]
        };
        throw error;
      }
    }
    this.state.records.messages.push({
      id: this.state.records.messages.length.toString(),
      from: idFrom,
      to: idTo,
      message: message?.text ?? "",
      wrap: message?.wrap ?? this.autoWrap(),
      type: messageType,
      activate
    });
    return true;
  }
  hasAtLeastOneBox() {
    return this.state.records.boxes.length > 0;
  }
  hasAtLeastOneBoxWithTitle() {
    return this.state.records.boxes.some((b) => b.name);
  }
  getMessages() {
    return this.state.records.messages;
  }
  getBoxes() {
    return this.state.records.boxes;
  }
  getActors() {
    return this.state.records.actors;
  }
  getCreatedActors() {
    return this.state.records.createdActors;
  }
  getDestroyedActors() {
    return this.state.records.destroyedActors;
  }
  getActor(id) {
    return this.state.records.actors.get(id);
  }
  getActorKeys() {
    return [...this.state.records.actors.keys()];
  }
  enableSequenceNumbers() {
    this.state.records.sequenceNumbersEnabled = true;
  }
  disableSequenceNumbers() {
    this.state.records.sequenceNumbersEnabled = false;
  }
  showSequenceNumbers() {
    return this.state.records.sequenceNumbersEnabled;
  }
  setWrap(wrapSetting) {
    this.state.records.wrapEnabled = wrapSetting;
  }
  extractWrap(text) {
    if (text === void 0) {
      return {};
    }
    text = text.trim();
    const wrap = /^:?wrap:/.exec(text) !== null ? true : /^:?nowrap:/.exec(text) !== null ? false : void 0;
    const cleanedText = (wrap === void 0 ? text : text.replace(/^:?(?:no)?wrap:/, "")).trim();
    return { cleanedText, wrap };
  }
  autoWrap() {
    if (this.state.records.wrapEnabled !== void 0) {
      return this.state.records.wrapEnabled;
    }
    return (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)().sequence?.wrap ?? false;
  }
  clear() {
    this.state.reset();
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
  }
  parseMessage(str) {
    const trimmedStr = str.trim();
    const { wrap, cleanedText } = this.extractWrap(trimmedStr);
    const message = {
      text: cleanedText,
      wrap
    };
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(`parseMessage: ${JSON.stringify(message)}`);
    return message;
  }
  // We expect the box statement to be color first then description
  // The color can be rgb,rgba,hsl,hsla, or css code names  #hex codes are not supported for now because of the way the char # is handled
  // We extract first segment as color, the rest of the line is considered as text
  parseBoxData(str) {
    const match = /^((?:rgba?|hsla?)\s*\(.*\)|\w*)(.*)$/.exec(str);
    let color = match?.[1] ? match[1].trim() : "transparent";
    let title = match?.[2] ? match[2].trim() : void 0;
    if (window?.CSS) {
      if (!window.CSS.supports("color", color)) {
        color = "transparent";
        title = str.trim();
      }
    } else {
      const style = new Option().style;
      style.color = color;
      if (style.color !== color) {
        color = "transparent";
        title = str.trim();
      }
    }
    const { wrap, cleanedText } = this.extractWrap(title);
    return {
      text: cleanedText ? (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(cleanedText, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()) : void 0,
      color,
      wrap
    };
  }
  addNote(actor, placement, message) {
    const note = {
      actor,
      placement,
      message: message.text,
      wrap: message.wrap ?? this.autoWrap()
    };
    const actors = [].concat(actor, actor);
    this.state.records.notes.push(note);
    this.state.records.messages.push({
      id: this.state.records.messages.length.toString(),
      from: actors[0],
      to: actors[1],
      message: message.text,
      wrap: message.wrap ?? this.autoWrap(),
      type: this.LINETYPE.NOTE,
      placement
    });
  }
  addLinks(actorId, text) {
    const actor = this.getActor(actorId);
    try {
      let sanitizedText = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(text.text, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
      sanitizedText = sanitizedText.replace(/&equals;/g, "=");
      sanitizedText = sanitizedText.replace(/&amp;/g, "&");
      const links = JSON.parse(sanitizedText);
      this.insertLinks(actor, links);
    } catch (e) {
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.error("error while parsing actor link text", e);
    }
  }
  addALink(actorId, text) {
    const actor = this.getActor(actorId);
    try {
      const links = {};
      let sanitizedText = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(text.text, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
      const sep = sanitizedText.indexOf("@");
      sanitizedText = sanitizedText.replace(/&equals;/g, "=");
      sanitizedText = sanitizedText.replace(/&amp;/g, "&");
      const label = sanitizedText.slice(0, sep - 1).trim();
      const link = sanitizedText.slice(sep + 1).trim();
      links[label] = link;
      this.insertLinks(actor, links);
    } catch (e) {
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.error("error while parsing actor link text", e);
    }
  }
  insertLinks(actor, links) {
    if (actor.links == null) {
      actor.links = links;
    } else {
      for (const key in links) {
        actor.links[key] = links[key];
      }
    }
  }
  addProperties(actorId, text) {
    const actor = this.getActor(actorId);
    try {
      const sanitizedText = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(text.text, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
      const properties = JSON.parse(sanitizedText);
      this.insertProperties(actor, properties);
    } catch (e) {
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.error("error while parsing actor properties text", e);
    }
  }
  insertProperties(actor, properties) {
    if (actor.properties == null) {
      actor.properties = properties;
    } else {
      for (const key in properties) {
        actor.properties[key] = properties[key];
      }
    }
  }
  boxEnd() {
    this.state.records.currentBox = void 0;
  }
  addDetails(actorId, text) {
    const actor = this.getActor(actorId);
    const elem = document.getElementById(text.text);
    try {
      const text2 = elem.innerHTML;
      const details = JSON.parse(text2);
      if (details.properties) {
        this.insertProperties(actor, details.properties);
      }
      if (details.links) {
        this.insertLinks(actor, details.links);
      }
    } catch (e) {
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.error("error while parsing actor details text", e);
    }
  }
  getActorProperty(actor, key) {
    if (actor?.properties !== void 0) {
      return actor.properties[key];
    }
    return void 0;
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any, @typescript-eslint/no-redundant-type-constituents
  apply(param) {
    if (Array.isArray(param)) {
      param.forEach((item) => {
        this.apply(item);
      });
    } else {
      switch (param.type) {
        case "sequenceIndex":
          this.state.records.messages.push({
            id: this.state.records.messages.length.toString(),
            from: void 0,
            to: void 0,
            message: {
              start: param.sequenceIndex,
              step: param.sequenceIndexStep,
              visible: param.sequenceVisible
            },
            wrap: false,
            type: param.signalType
          });
          break;
        case "addParticipant":
          this.addActor(param.actor, param.actor, param.description, param.draw);
          break;
        case "createParticipant":
          if (this.state.records.actors.has(param.actor)) {
            throw new Error(
              "It is not possible to have actors with the same id, even if one is destroyed before the next is created. Use 'AS' aliases to simulate the behavior"
            );
          }
          this.state.records.lastCreated = param.actor;
          this.addActor(param.actor, param.actor, param.description, param.draw);
          this.state.records.createdActors.set(param.actor, this.state.records.messages.length);
          break;
        case "destroyParticipant":
          this.state.records.lastDestroyed = param.actor;
          this.state.records.destroyedActors.set(param.actor, this.state.records.messages.length);
          break;
        case "activeStart":
          this.addSignal(param.actor, void 0, void 0, param.signalType);
          break;
        case "activeEnd":
          this.addSignal(param.actor, void 0, void 0, param.signalType);
          break;
        case "addNote":
          this.addNote(param.actor, param.placement, param.text);
          break;
        case "addLinks":
          this.addLinks(param.actor, param.text);
          break;
        case "addALink":
          this.addALink(param.actor, param.text);
          break;
        case "addProperties":
          this.addProperties(param.actor, param.text);
          break;
        case "addDetails":
          this.addDetails(param.actor, param.text);
          break;
        case "addMessage":
          if (this.state.records.lastCreated) {
            if (param.to !== this.state.records.lastCreated) {
              throw new Error(
                "The created participant " + this.state.records.lastCreated.name + " does not have an associated creating message after its declaration. Please check the sequence diagram."
              );
            } else {
              this.state.records.lastCreated = void 0;
            }
          } else if (this.state.records.lastDestroyed) {
            if (param.to !== this.state.records.lastDestroyed && param.from !== this.state.records.lastDestroyed) {
              throw new Error(
                "The destroyed participant " + this.state.records.lastDestroyed.name + " does not have an associated destroying message after its declaration. Please check the sequence diagram."
              );
            } else {
              this.state.records.lastDestroyed = void 0;
            }
          }
          this.addSignal(param.from, param.to, param.msg, param.signalType, param.activate);
          break;
        case "boxStart":
          this.addBox(param.boxData);
          break;
        case "boxEnd":
          this.boxEnd();
          break;
        case "loopStart":
          this.addSignal(void 0, void 0, param.loopText, param.signalType);
          break;
        case "loopEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "rectStart":
          this.addSignal(void 0, void 0, param.color, param.signalType);
          break;
        case "rectEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "optStart":
          this.addSignal(void 0, void 0, param.optText, param.signalType);
          break;
        case "optEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "altStart":
          this.addSignal(void 0, void 0, param.altText, param.signalType);
          break;
        case "else":
          this.addSignal(void 0, void 0, param.altText, param.signalType);
          break;
        case "altEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "setAccTitle":
          (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN)(param.text);
          break;
        case "parStart":
          this.addSignal(void 0, void 0, param.parText, param.signalType);
          break;
        case "and":
          this.addSignal(void 0, void 0, param.parText, param.signalType);
          break;
        case "parEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "criticalStart":
          this.addSignal(void 0, void 0, param.criticalText, param.signalType);
          break;
        case "option":
          this.addSignal(void 0, void 0, param.optionText, param.signalType);
          break;
        case "criticalEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
        case "breakStart":
          this.addSignal(void 0, void 0, param.breakText, param.signalType);
          break;
        case "breakEnd":
          this.addSignal(void 0, void 0, void 0, param.signalType);
          break;
      }
    }
  }
  getConfig() {
    return (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)().sequence;
  }
};

// src/diagrams/sequence/styles.js
var getStyles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((options) => `.actor {
    stroke: ${options.actorBorder};
    fill: ${options.actorBkg};
  }

  text.actor > tspan {
    fill: ${options.actorTextColor};
    stroke: none;
  }

  .actor-line {
    stroke: ${options.actorLineColor};
  }

  .messageLine0 {
    stroke-width: 1.5;
    stroke-dasharray: none;
    stroke: ${options.signalColor};
  }

  .messageLine1 {
    stroke-width: 1.5;
    stroke-dasharray: 2, 2;
    stroke: ${options.signalColor};
  }

  #arrowhead path {
    fill: ${options.signalColor};
    stroke: ${options.signalColor};
  }

  .sequenceNumber {
    fill: ${options.sequenceNumberColor};
  }

  #sequencenumber {
    fill: ${options.signalColor};
  }

  #crosshead path {
    fill: ${options.signalColor};
    stroke: ${options.signalColor};
  }

  .messageText {
    fill: ${options.signalTextColor};
    stroke: none;
  }

  .labelBox {
    stroke: ${options.labelBoxBorderColor};
    fill: ${options.labelBoxBkgColor};
  }

  .labelText, .labelText > tspan {
    fill: ${options.labelTextColor};
    stroke: none;
  }

  .loopText, .loopText > tspan {
    fill: ${options.loopTextColor};
    stroke: none;
  }

  .loopLine {
    stroke-width: 2px;
    stroke-dasharray: 2, 2;
    stroke: ${options.labelBoxBorderColor};
    fill: ${options.labelBoxBorderColor};
  }

  .note {
    //stroke: #decc93;
    stroke: ${options.noteBorderColor};
    fill: ${options.noteBkgColor};
  }

  .noteText, .noteText > tspan {
    fill: ${options.noteTextColor};
    stroke: none;
  }

  .activation0 {
    fill: ${options.activationBkgColor};
    stroke: ${options.activationBorderColor};
  }

  .activation1 {
    fill: ${options.activationBkgColor};
    stroke: ${options.activationBorderColor};
  }

  .activation2 {
    fill: ${options.activationBkgColor};
    stroke: ${options.activationBorderColor};
  }

  .actorPopupMenu {
    position: absolute;
  }

  .actorPopupMenuPanel {
    position: absolute;
    fill: ${options.actorBkg};
    box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
    filter: drop-shadow(3px 5px 2px rgb(0 0 0 / 0.4));
}
  .actor-man line {
    stroke: ${options.actorBorder};
    fill: ${options.actorBkg};
  }
  .actor-man circle, line {
    stroke: ${options.actorBorder};
    fill: ${options.actorBkg};
    stroke-width: 2px;
  }
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/sequence/sequenceRenderer.ts


// src/diagrams/sequence/svgDraw.js

var ACTOR_TYPE_WIDTH = 18 * 2;
var TOP_ACTOR_CLASS = "actor-top";
var BOTTOM_ACTOR_CLASS = "actor-bottom";
var ACTOR_BOX_CLASS = "actor-box";
var ACTOR_MAN_FIGURE_CLASS = "actor-man";
var drawRect2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, rectData) {
  return (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawRect */ .Mu)(elem, rectData);
}, "drawRect");
var drawPopup = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, actor, minMenuWidth, textAttrs, forceMenus) {
  if (actor.links === void 0 || actor.links === null || Object.keys(actor.links).length === 0) {
    return { height: 0, width: 0 };
  }
  const links = actor.links;
  const actorCnt2 = actor.actorCnt;
  const rectData = actor.rectData;
  var displayValue = "none";
  if (forceMenus) {
    displayValue = "block !important";
  }
  const g = elem.append("g");
  g.attr("id", "actor" + actorCnt2 + "_popup");
  g.attr("class", "actorPopupMenu");
  g.attr("display", displayValue);
  var actorClass = "";
  if (rectData.class !== void 0) {
    actorClass = " " + rectData.class;
  }
  let menuWidth = rectData.width > minMenuWidth ? rectData.width : minMenuWidth;
  const rectElem = g.append("rect");
  rectElem.attr("class", "actorPopupMenuPanel" + actorClass);
  rectElem.attr("x", rectData.x);
  rectElem.attr("y", rectData.height);
  rectElem.attr("fill", rectData.fill);
  rectElem.attr("stroke", rectData.stroke);
  rectElem.attr("width", menuWidth);
  rectElem.attr("height", rectData.height);
  rectElem.attr("rx", rectData.rx);
  rectElem.attr("ry", rectData.ry);
  if (links != null) {
    var linkY = 20;
    for (let key in links) {
      var linkElem = g.append("a");
      var sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_5__/* .sanitizeUrl */ .N)(links[key]);
      linkElem.attr("xlink:href", sanitizedLink);
      linkElem.attr("target", "_blank");
      _drawMenuItemTextCandidateFunc(textAttrs)(
        key,
        linkElem,
        rectData.x + 10,
        rectData.height + linkY,
        menuWidth,
        20,
        { class: "actor" },
        textAttrs
      );
      linkY += 30;
    }
  }
  rectElem.attr("height", linkY);
  return { height: rectData.height + linkY, width: menuWidth };
}, "drawPopup");
var popupMenuToggle = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(popId) {
  return "var pu = document.getElementById('" + popId + "'); if (pu != null) { pu.style.display = pu.style.display == 'block' ? 'none' : 'block'; }";
}, "popupMenuToggle");
var drawKatex = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(elem, textData, msgModel = null) {
  let textElem = elem.append("foreignObject");
  const lines = await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .renderKatex */ .uT)(textData.text, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)());
  const divElem = textElem.append("xhtml:div").attr("style", "width: fit-content;").attr("xmlns", "http://www.w3.org/1999/xhtml").html(lines);
  const dim = divElem.node().getBoundingClientRect();
  textElem.attr("height", Math.round(dim.height)).attr("width", Math.round(dim.width));
  if (textData.class === "noteText") {
    const rectElem = elem.node().firstChild;
    rectElem.setAttribute("height", dim.height + 2 * textData.textMargin);
    const rectDim = rectElem.getBBox();
    textElem.attr("x", Math.round(rectDim.x + rectDim.width / 2 - dim.width / 2)).attr("y", Math.round(rectDim.y + rectDim.height / 2 - dim.height / 2));
  } else if (msgModel) {
    let { startx, stopx, starty } = msgModel;
    if (startx > stopx) {
      const temp = startx;
      startx = stopx;
      stopx = temp;
    }
    textElem.attr("x", Math.round(startx + Math.abs(startx - stopx) / 2 - dim.width / 2));
    if (textData.class === "loopText") {
      textElem.attr("y", Math.round(starty));
    } else {
      textElem.attr("y", Math.round(starty - dim.height));
    }
  }
  return [textElem];
}, "drawKatex");
var drawText = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, textData) {
  let prevTextHeight = 0;
  let textHeight = 0;
  const lines = textData.text.split(_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.lineBreakRegex);
  const [_textFontSize, _textFontSizePx] = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .parseFontSize */ .VG)(textData.fontSize);
  let textElems = [];
  let dy = 0;
  let yfunc = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => textData.y, "yfunc");
  if (textData.valign !== void 0 && textData.textMargin !== void 0 && textData.textMargin > 0) {
    switch (textData.valign) {
      case "top":
      case "start":
        yfunc = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => Math.round(textData.y + textData.textMargin), "yfunc");
        break;
      case "middle":
      case "center":
        yfunc = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => Math.round(textData.y + (prevTextHeight + textHeight + textData.textMargin) / 2), "yfunc");
        break;
      case "bottom":
      case "end":
        yfunc = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => Math.round(
          textData.y + (prevTextHeight + textHeight + 2 * textData.textMargin) - textData.textMargin
        ), "yfunc");
        break;
    }
  }
  if (textData.anchor !== void 0 && textData.textMargin !== void 0 && textData.width !== void 0) {
    switch (textData.anchor) {
      case "left":
      case "start":
        textData.x = Math.round(textData.x + textData.textMargin);
        textData.anchor = "start";
        textData.dominantBaseline = "middle";
        textData.alignmentBaseline = "middle";
        break;
      case "middle":
      case "center":
        textData.x = Math.round(textData.x + textData.width / 2);
        textData.anchor = "middle";
        textData.dominantBaseline = "middle";
        textData.alignmentBaseline = "middle";
        break;
      case "right":
      case "end":
        textData.x = Math.round(textData.x + textData.width - textData.textMargin);
        textData.anchor = "end";
        textData.dominantBaseline = "middle";
        textData.alignmentBaseline = "middle";
        break;
    }
  }
  for (let [i, line] of lines.entries()) {
    if (textData.textMargin !== void 0 && textData.textMargin === 0 && _textFontSize !== void 0) {
      dy = i * _textFontSize;
    }
    const textElem = elem.append("text");
    textElem.attr("x", textData.x);
    textElem.attr("y", yfunc());
    if (textData.anchor !== void 0) {
      textElem.attr("text-anchor", textData.anchor).attr("dominant-baseline", textData.dominantBaseline).attr("alignment-baseline", textData.alignmentBaseline);
    }
    if (textData.fontFamily !== void 0) {
      textElem.style("font-family", textData.fontFamily);
    }
    if (_textFontSizePx !== void 0) {
      textElem.style("font-size", _textFontSizePx);
    }
    if (textData.fontWeight !== void 0) {
      textElem.style("font-weight", textData.fontWeight);
    }
    if (textData.fill !== void 0) {
      textElem.attr("fill", textData.fill);
    }
    if (textData.class !== void 0) {
      textElem.attr("class", textData.class);
    }
    if (textData.dy !== void 0) {
      textElem.attr("dy", textData.dy);
    } else if (dy !== 0) {
      textElem.attr("dy", dy);
    }
    const text = line || _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .ZERO_WIDTH_SPACE */ .$m;
    if (textData.tspan) {
      const span = textElem.append("tspan");
      span.attr("x", textData.x);
      if (textData.fill !== void 0) {
        span.attr("fill", textData.fill);
      }
      span.text(text);
    } else {
      textElem.text(text);
    }
    if (textData.valign !== void 0 && textData.textMargin !== void 0 && textData.textMargin > 0) {
      textHeight += (textElem._groups || textElem)[0][0].getBBox().height;
      prevTextHeight = textHeight;
    }
    textElems.push(textElem);
  }
  return textElems;
}, "drawText");
var drawLabel = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, txtObject) {
  function genPoints(x, y, width, height, cut) {
    return x + "," + y + " " + (x + width) + "," + y + " " + (x + width) + "," + (y + height - cut) + " " + (x + width - cut * 1.2) + "," + (y + height) + " " + x + "," + (y + height);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(genPoints, "genPoints");
  const polygon = elem.append("polygon");
  polygon.attr("points", genPoints(txtObject.x, txtObject.y, txtObject.width, txtObject.height, 7));
  polygon.attr("class", "labelBox");
  txtObject.y = txtObject.y + txtObject.height / 2;
  drawText(elem, txtObject);
  return polygon;
}, "drawLabel");
var actorCnt = -1;
var fixLifeLineHeights = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((diagram2, actors, actorKeys, conf2) => {
  if (!diagram2.select) {
    return;
  }
  actorKeys.forEach((actorKey) => {
    const actor = actors.get(actorKey);
    const actorDOM = diagram2.select("#actor" + actor.actorCnt);
    if (!conf2.mirrorActors && actor.stopy) {
      actorDOM.attr("y2", actor.stopy + actor.height / 2);
    } else if (conf2.mirrorActors) {
      actorDOM.attr("y2", actor.stopy);
    }
  });
}, "fixLifeLineHeights");
var drawActorTypeParticipant = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, actor, conf2, isFooter) {
  const actorY = isFooter ? actor.stopy : actor.starty;
  const center = actor.x + actor.width / 2;
  const centerY = actorY + actor.height;
  const boxplusLineGroup = elem.append("g").lower();
  var g = boxplusLineGroup;
  if (!isFooter) {
    actorCnt++;
    if (Object.keys(actor.links || {}).length && !conf2.forceMenus) {
      g.attr("onclick", popupMenuToggle(`actor${actorCnt}_popup`)).attr("cursor", "pointer");
    }
    g.append("line").attr("id", "actor" + actorCnt).attr("x1", center).attr("y1", centerY).attr("x2", center).attr("y2", 2e3).attr("class", "actor-line 200").attr("stroke-width", "0.5px").attr("stroke", "#999").attr("name", actor.name);
    g = boxplusLineGroup.append("g");
    actor.actorCnt = actorCnt;
    if (actor.links != null) {
      g.attr("id", "root-" + actorCnt);
    }
  }
  const rect = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  var cssclass = "actor";
  if (actor.properties?.class) {
    cssclass = actor.properties.class;
  } else {
    rect.fill = "#eaeaea";
  }
  if (isFooter) {
    cssclass += ` ${BOTTOM_ACTOR_CLASS}`;
  } else {
    cssclass += ` ${TOP_ACTOR_CLASS}`;
  }
  rect.x = actor.x;
  rect.y = actorY;
  rect.width = actor.width;
  rect.height = actor.height;
  rect.class = cssclass;
  rect.rx = 3;
  rect.ry = 3;
  rect.name = actor.name;
  const rectElem = drawRect2(g, rect);
  actor.rectData = rect;
  if (actor.properties?.icon) {
    const iconSrc = actor.properties.icon.trim();
    if (iconSrc.charAt(0) === "@") {
      (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawEmbeddedImage */ .rB)(g, rect.x + rect.width - 20, rect.y + 10, iconSrc.substr(1));
    } else {
      (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawImage */ .AE)(g, rect.x + rect.width - 20, rect.y + 10, iconSrc);
    }
  }
  _drawTextCandidateFunc(conf2, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(actor.description))(
    actor.description,
    g,
    rect.x,
    rect.y,
    rect.width,
    rect.height,
    { class: `actor ${ACTOR_BOX_CLASS}` },
    conf2
  );
  let height = actor.height;
  if (rectElem.node) {
    const bounds2 = rectElem.node().getBBox();
    actor.height = bounds2.height;
    height = bounds2.height;
  }
  return height;
}, "drawActorTypeParticipant");
var drawActorTypeActor = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, actor, conf2, isFooter) {
  const actorY = isFooter ? actor.stopy : actor.starty;
  const center = actor.x + actor.width / 2;
  const centerY = actorY + 80;
  const line = elem.append("g").lower();
  if (!isFooter) {
    actorCnt++;
    line.append("line").attr("id", "actor" + actorCnt).attr("x1", center).attr("y1", centerY).attr("x2", center).attr("y2", 2e3).attr("class", "actor-line 200").attr("stroke-width", "0.5px").attr("stroke", "#999").attr("name", actor.name);
    actor.actorCnt = actorCnt;
  }
  const actElem = elem.append("g");
  let cssClass = ACTOR_MAN_FIGURE_CLASS;
  if (isFooter) {
    cssClass += ` ${BOTTOM_ACTOR_CLASS}`;
  } else {
    cssClass += ` ${TOP_ACTOR_CLASS}`;
  }
  actElem.attr("class", cssClass);
  actElem.attr("name", actor.name);
  const rect = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  rect.x = actor.x;
  rect.y = actorY;
  rect.fill = "#eaeaea";
  rect.width = actor.width;
  rect.height = actor.height;
  rect.class = "actor";
  rect.rx = 3;
  rect.ry = 3;
  actElem.append("line").attr("id", "actor-man-torso" + actorCnt).attr("x1", center).attr("y1", actorY + 25).attr("x2", center).attr("y2", actorY + 45);
  actElem.append("line").attr("id", "actor-man-arms" + actorCnt).attr("x1", center - ACTOR_TYPE_WIDTH / 2).attr("y1", actorY + 33).attr("x2", center + ACTOR_TYPE_WIDTH / 2).attr("y2", actorY + 33);
  actElem.append("line").attr("x1", center - ACTOR_TYPE_WIDTH / 2).attr("y1", actorY + 60).attr("x2", center).attr("y2", actorY + 45);
  actElem.append("line").attr("x1", center).attr("y1", actorY + 45).attr("x2", center + ACTOR_TYPE_WIDTH / 2 - 2).attr("y2", actorY + 60);
  const circle = actElem.append("circle");
  circle.attr("cx", actor.x + actor.width / 2);
  circle.attr("cy", actorY + 10);
  circle.attr("r", 15);
  circle.attr("width", actor.width);
  circle.attr("height", actor.height);
  const bounds2 = actElem.node().getBBox();
  actor.height = bounds2.height;
  _drawTextCandidateFunc(conf2, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(actor.description))(
    actor.description,
    actElem,
    rect.x,
    rect.y + 35,
    rect.width,
    rect.height,
    { class: `actor ${ACTOR_MAN_FIGURE_CLASS}` },
    conf2
  );
  return actor.height;
}, "drawActorTypeActor");
var drawActor = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(elem, actor, conf2, isFooter) {
  switch (actor.type) {
    case "actor":
      return await drawActorTypeActor(elem, actor, conf2, isFooter);
    case "participant":
      return await drawActorTypeParticipant(elem, actor, conf2, isFooter);
  }
}, "drawActor");
var drawBox = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, box, conf2) {
  const boxplusTextGroup = elem.append("g");
  const g = boxplusTextGroup;
  drawBackgroundRect2(g, box);
  if (box.name) {
    _drawTextCandidateFunc(conf2)(
      box.name,
      g,
      box.x,
      box.y + (box.textMaxHeight || 0) / 2,
      box.width,
      0,
      { class: "text" },
      conf2
    );
  }
  g.lower();
}, "drawBox");
var anchorElement = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  return elem.append("g");
}, "anchorElement");
var drawActivation = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, bounds2, verticalPos, conf2, actorActivations2) {
  const rect = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  const g = bounds2.anchored;
  rect.x = bounds2.startx;
  rect.y = bounds2.starty;
  rect.class = "activation" + actorActivations2 % 3;
  rect.width = bounds2.stopx - bounds2.startx;
  rect.height = verticalPos - bounds2.starty;
  drawRect2(g, rect);
}, "drawActivation");
var drawLoop = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(elem, loopModel, labelText, conf2) {
  const {
    boxMargin,
    boxTextMargin,
    labelBoxHeight,
    labelBoxWidth,
    messageFontFamily: fontFamily,
    messageFontSize: fontSize,
    messageFontWeight: fontWeight
  } = conf2;
  const g = elem.append("g");
  const drawLoopLine = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(startx, starty, stopx, stopy) {
    return g.append("line").attr("x1", startx).attr("y1", starty).attr("x2", stopx).attr("y2", stopy).attr("class", "loopLine");
  }, "drawLoopLine");
  drawLoopLine(loopModel.startx, loopModel.starty, loopModel.stopx, loopModel.starty);
  drawLoopLine(loopModel.stopx, loopModel.starty, loopModel.stopx, loopModel.stopy);
  drawLoopLine(loopModel.startx, loopModel.stopy, loopModel.stopx, loopModel.stopy);
  drawLoopLine(loopModel.startx, loopModel.starty, loopModel.startx, loopModel.stopy);
  if (loopModel.sections !== void 0) {
    loopModel.sections.forEach(function(item) {
      drawLoopLine(loopModel.startx, item.y, loopModel.stopx, item.y).style(
        "stroke-dasharray",
        "3, 3"
      );
    });
  }
  let txt = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getTextObj */ .AD)();
  txt.text = labelText;
  txt.x = loopModel.startx;
  txt.y = loopModel.starty;
  txt.fontFamily = fontFamily;
  txt.fontSize = fontSize;
  txt.fontWeight = fontWeight;
  txt.anchor = "middle";
  txt.valign = "middle";
  txt.tspan = false;
  txt.width = labelBoxWidth || 50;
  txt.height = labelBoxHeight || 20;
  txt.textMargin = boxTextMargin;
  txt.class = "labelText";
  drawLabel(g, txt);
  txt = getTextObj2();
  txt.text = loopModel.title;
  txt.x = loopModel.startx + labelBoxWidth / 2 + (loopModel.stopx - loopModel.startx) / 2;
  txt.y = loopModel.starty + boxMargin + boxTextMargin;
  txt.anchor = "middle";
  txt.valign = "middle";
  txt.textMargin = boxTextMargin;
  txt.class = "loopText";
  txt.fontFamily = fontFamily;
  txt.fontSize = fontSize;
  txt.fontWeight = fontWeight;
  txt.wrap = true;
  let textElem = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(txt.text) ? await drawKatex(g, txt, loopModel) : drawText(g, txt);
  if (loopModel.sectionTitles !== void 0) {
    for (const [idx, item] of Object.entries(loopModel.sectionTitles)) {
      if (item.message) {
        txt.text = item.message;
        txt.x = loopModel.startx + (loopModel.stopx - loopModel.startx) / 2;
        txt.y = loopModel.sections[idx].y + boxMargin + boxTextMargin;
        txt.class = "loopText";
        txt.anchor = "middle";
        txt.valign = "middle";
        txt.tspan = false;
        txt.fontFamily = fontFamily;
        txt.fontSize = fontSize;
        txt.fontWeight = fontWeight;
        txt.wrap = loopModel.wrap;
        if ((0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(txt.text)) {
          loopModel.starty = loopModel.sections[idx].y;
          await drawKatex(g, txt, loopModel);
        } else {
          drawText(g, txt);
        }
        let sectionHeight = Math.round(
          textElem.map((te) => (te._groups || te)[0][0].getBBox().height).reduce((acc, curr) => acc + curr)
        );
        loopModel.sections[idx].height += sectionHeight - (boxMargin + boxTextMargin);
      }
    }
  }
  loopModel.height = Math.round(loopModel.stopy - loopModel.starty);
  return g;
}, "drawLoop");
var drawBackgroundRect2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, bounds2) {
  (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawBackgroundRect */ .O)(elem, bounds2);
}, "drawBackgroundRect");
var insertDatabaseIcon = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("symbol").attr("id", "database").attr("fill-rule", "evenodd").attr("clip-rule", "evenodd").append("path").attr("transform", "scale(.5)").attr(
    "d",
    "M12.258.001l.256.004.255.005.253.008.251.01.249.012.247.015.246.016.242.019.241.02.239.023.236.024.233.027.231.028.229.031.225.032.223.034.22.036.217.038.214.04.211.041.208.043.205.045.201.046.198.048.194.05.191.051.187.053.183.054.18.056.175.057.172.059.168.06.163.061.16.063.155.064.15.066.074.033.073.033.071.034.07.034.069.035.068.035.067.035.066.035.064.036.064.036.062.036.06.036.06.037.058.037.058.037.055.038.055.038.053.038.052.038.051.039.05.039.048.039.047.039.045.04.044.04.043.04.041.04.04.041.039.041.037.041.036.041.034.041.033.042.032.042.03.042.029.042.027.042.026.043.024.043.023.043.021.043.02.043.018.044.017.043.015.044.013.044.012.044.011.045.009.044.007.045.006.045.004.045.002.045.001.045v17l-.001.045-.002.045-.004.045-.006.045-.007.045-.009.044-.011.045-.012.044-.013.044-.015.044-.017.043-.018.044-.02.043-.021.043-.023.043-.024.043-.026.043-.027.042-.029.042-.03.042-.032.042-.033.042-.034.041-.036.041-.037.041-.039.041-.04.041-.041.04-.043.04-.044.04-.045.04-.047.039-.048.039-.05.039-.051.039-.052.038-.053.038-.055.038-.055.038-.058.037-.058.037-.06.037-.06.036-.062.036-.064.036-.064.036-.066.035-.067.035-.068.035-.069.035-.07.034-.071.034-.073.033-.074.033-.15.066-.155.064-.16.063-.163.061-.168.06-.172.059-.175.057-.18.056-.183.054-.187.053-.191.051-.194.05-.198.048-.201.046-.205.045-.208.043-.211.041-.214.04-.217.038-.22.036-.223.034-.225.032-.229.031-.231.028-.233.027-.236.024-.239.023-.241.02-.242.019-.246.016-.247.015-.249.012-.251.01-.253.008-.255.005-.256.004-.258.001-.258-.001-.256-.004-.255-.005-.253-.008-.251-.01-.249-.012-.247-.015-.245-.016-.243-.019-.241-.02-.238-.023-.236-.024-.234-.027-.231-.028-.228-.031-.226-.032-.223-.034-.22-.036-.217-.038-.214-.04-.211-.041-.208-.043-.204-.045-.201-.046-.198-.048-.195-.05-.19-.051-.187-.053-.184-.054-.179-.056-.176-.057-.172-.059-.167-.06-.164-.061-.159-.063-.155-.064-.151-.066-.074-.033-.072-.033-.072-.034-.07-.034-.069-.035-.068-.035-.067-.035-.066-.035-.064-.036-.063-.036-.062-.036-.061-.036-.06-.037-.058-.037-.057-.037-.056-.038-.055-.038-.053-.038-.052-.038-.051-.039-.049-.039-.049-.039-.046-.039-.046-.04-.044-.04-.043-.04-.041-.04-.04-.041-.039-.041-.037-.041-.036-.041-.034-.041-.033-.042-.032-.042-.03-.042-.029-.042-.027-.042-.026-.043-.024-.043-.023-.043-.021-.043-.02-.043-.018-.044-.017-.043-.015-.044-.013-.044-.012-.044-.011-.045-.009-.044-.007-.045-.006-.045-.004-.045-.002-.045-.001-.045v-17l.001-.045.002-.045.004-.045.006-.045.007-.045.009-.044.011-.045.012-.044.013-.044.015-.044.017-.043.018-.044.02-.043.021-.043.023-.043.024-.043.026-.043.027-.042.029-.042.03-.042.032-.042.033-.042.034-.041.036-.041.037-.041.039-.041.04-.041.041-.04.043-.04.044-.04.046-.04.046-.039.049-.039.049-.039.051-.039.052-.038.053-.038.055-.038.056-.038.057-.037.058-.037.06-.037.061-.036.062-.036.063-.036.064-.036.066-.035.067-.035.068-.035.069-.035.07-.034.072-.034.072-.033.074-.033.151-.066.155-.064.159-.063.164-.061.167-.06.172-.059.176-.057.179-.056.184-.054.187-.053.19-.051.195-.05.198-.048.201-.046.204-.045.208-.043.211-.041.214-.04.217-.038.22-.036.223-.034.226-.032.228-.031.231-.028.234-.027.236-.024.238-.023.241-.02.243-.019.245-.016.247-.015.249-.012.251-.01.253-.008.255-.005.256-.004.258-.001.258.001zm-9.258 20.499v.01l.001.021.003.021.004.022.005.021.006.022.007.022.009.023.01.022.011.023.012.023.013.023.015.023.016.024.017.023.018.024.019.024.021.024.022.025.023.024.024.025.052.049.056.05.061.051.066.051.07.051.075.051.079.052.084.052.088.052.092.052.097.052.102.051.105.052.11.052.114.051.119.051.123.051.127.05.131.05.135.05.139.048.144.049.147.047.152.047.155.047.16.045.163.045.167.043.171.043.176.041.178.041.183.039.187.039.19.037.194.035.197.035.202.033.204.031.209.03.212.029.216.027.219.025.222.024.226.021.23.02.233.018.236.016.24.015.243.012.246.01.249.008.253.005.256.004.259.001.26-.001.257-.004.254-.005.25-.008.247-.011.244-.012.241-.014.237-.016.233-.018.231-.021.226-.021.224-.024.22-.026.216-.027.212-.028.21-.031.205-.031.202-.034.198-.034.194-.036.191-.037.187-.039.183-.04.179-.04.175-.042.172-.043.168-.044.163-.045.16-.046.155-.046.152-.047.148-.048.143-.049.139-.049.136-.05.131-.05.126-.05.123-.051.118-.052.114-.051.11-.052.106-.052.101-.052.096-.052.092-.052.088-.053.083-.051.079-.052.074-.052.07-.051.065-.051.06-.051.056-.05.051-.05.023-.024.023-.025.021-.024.02-.024.019-.024.018-.024.017-.024.015-.023.014-.024.013-.023.012-.023.01-.023.01-.022.008-.022.006-.022.006-.022.004-.022.004-.021.001-.021.001-.021v-4.127l-.077.055-.08.053-.083.054-.085.053-.087.052-.09.052-.093.051-.095.05-.097.05-.1.049-.102.049-.105.048-.106.047-.109.047-.111.046-.114.045-.115.045-.118.044-.12.043-.122.042-.124.042-.126.041-.128.04-.13.04-.132.038-.134.038-.135.037-.138.037-.139.035-.142.035-.143.034-.144.033-.147.032-.148.031-.15.03-.151.03-.153.029-.154.027-.156.027-.158.026-.159.025-.161.024-.162.023-.163.022-.165.021-.166.02-.167.019-.169.018-.169.017-.171.016-.173.015-.173.014-.175.013-.175.012-.177.011-.178.01-.179.008-.179.008-.181.006-.182.005-.182.004-.184.003-.184.002h-.37l-.184-.002-.184-.003-.182-.004-.182-.005-.181-.006-.179-.008-.179-.008-.178-.01-.176-.011-.176-.012-.175-.013-.173-.014-.172-.015-.171-.016-.17-.017-.169-.018-.167-.019-.166-.02-.165-.021-.163-.022-.162-.023-.161-.024-.159-.025-.157-.026-.156-.027-.155-.027-.153-.029-.151-.03-.15-.03-.148-.031-.146-.032-.145-.033-.143-.034-.141-.035-.14-.035-.137-.037-.136-.037-.134-.038-.132-.038-.13-.04-.128-.04-.126-.041-.124-.042-.122-.042-.12-.044-.117-.043-.116-.045-.113-.045-.112-.046-.109-.047-.106-.047-.105-.048-.102-.049-.1-.049-.097-.05-.095-.05-.093-.052-.09-.051-.087-.052-.085-.053-.083-.054-.08-.054-.077-.054v4.127zm0-5.654v.011l.001.021.003.021.004.021.005.022.006.022.007.022.009.022.01.022.011.023.012.023.013.023.015.024.016.023.017.024.018.024.019.024.021.024.022.024.023.025.024.024.052.05.056.05.061.05.066.051.07.051.075.052.079.051.084.052.088.052.092.052.097.052.102.052.105.052.11.051.114.051.119.052.123.05.127.051.131.05.135.049.139.049.144.048.147.048.152.047.155.046.16.045.163.045.167.044.171.042.176.042.178.04.183.04.187.038.19.037.194.036.197.034.202.033.204.032.209.03.212.028.216.027.219.025.222.024.226.022.23.02.233.018.236.016.24.014.243.012.246.01.249.008.253.006.256.003.259.001.26-.001.257-.003.254-.006.25-.008.247-.01.244-.012.241-.015.237-.016.233-.018.231-.02.226-.022.224-.024.22-.025.216-.027.212-.029.21-.03.205-.032.202-.033.198-.035.194-.036.191-.037.187-.039.183-.039.179-.041.175-.042.172-.043.168-.044.163-.045.16-.045.155-.047.152-.047.148-.048.143-.048.139-.05.136-.049.131-.05.126-.051.123-.051.118-.051.114-.052.11-.052.106-.052.101-.052.096-.052.092-.052.088-.052.083-.052.079-.052.074-.051.07-.052.065-.051.06-.05.056-.051.051-.049.023-.025.023-.024.021-.025.02-.024.019-.024.018-.024.017-.024.015-.023.014-.023.013-.024.012-.022.01-.023.01-.023.008-.022.006-.022.006-.022.004-.021.004-.022.001-.021.001-.021v-4.139l-.077.054-.08.054-.083.054-.085.052-.087.053-.09.051-.093.051-.095.051-.097.05-.1.049-.102.049-.105.048-.106.047-.109.047-.111.046-.114.045-.115.044-.118.044-.12.044-.122.042-.124.042-.126.041-.128.04-.13.039-.132.039-.134.038-.135.037-.138.036-.139.036-.142.035-.143.033-.144.033-.147.033-.148.031-.15.03-.151.03-.153.028-.154.028-.156.027-.158.026-.159.025-.161.024-.162.023-.163.022-.165.021-.166.02-.167.019-.169.018-.169.017-.171.016-.173.015-.173.014-.175.013-.175.012-.177.011-.178.009-.179.009-.179.007-.181.007-.182.005-.182.004-.184.003-.184.002h-.37l-.184-.002-.184-.003-.182-.004-.182-.005-.181-.007-.179-.007-.179-.009-.178-.009-.176-.011-.176-.012-.175-.013-.173-.014-.172-.015-.171-.016-.17-.017-.169-.018-.167-.019-.166-.02-.165-.021-.163-.022-.162-.023-.161-.024-.159-.025-.157-.026-.156-.027-.155-.028-.153-.028-.151-.03-.15-.03-.148-.031-.146-.033-.145-.033-.143-.033-.141-.035-.14-.036-.137-.036-.136-.037-.134-.038-.132-.039-.13-.039-.128-.04-.126-.041-.124-.042-.122-.043-.12-.043-.117-.044-.116-.044-.113-.046-.112-.046-.109-.046-.106-.047-.105-.048-.102-.049-.1-.049-.097-.05-.095-.051-.093-.051-.09-.051-.087-.053-.085-.052-.083-.054-.08-.054-.077-.054v4.139zm0-5.666v.011l.001.02.003.022.004.021.005.022.006.021.007.022.009.023.01.022.011.023.012.023.013.023.015.023.016.024.017.024.018.023.019.024.021.025.022.024.023.024.024.025.052.05.056.05.061.05.066.051.07.051.075.052.079.051.084.052.088.052.092.052.097.052.102.052.105.051.11.052.114.051.119.051.123.051.127.05.131.05.135.05.139.049.144.048.147.048.152.047.155.046.16.045.163.045.167.043.171.043.176.042.178.04.183.04.187.038.19.037.194.036.197.034.202.033.204.032.209.03.212.028.216.027.219.025.222.024.226.021.23.02.233.018.236.017.24.014.243.012.246.01.249.008.253.006.256.003.259.001.26-.001.257-.003.254-.006.25-.008.247-.01.244-.013.241-.014.237-.016.233-.018.231-.02.226-.022.224-.024.22-.025.216-.027.212-.029.21-.03.205-.032.202-.033.198-.035.194-.036.191-.037.187-.039.183-.039.179-.041.175-.042.172-.043.168-.044.163-.045.16-.045.155-.047.152-.047.148-.048.143-.049.139-.049.136-.049.131-.051.126-.05.123-.051.118-.052.114-.051.11-.052.106-.052.101-.052.096-.052.092-.052.088-.052.083-.052.079-.052.074-.052.07-.051.065-.051.06-.051.056-.05.051-.049.023-.025.023-.025.021-.024.02-.024.019-.024.018-.024.017-.024.015-.023.014-.024.013-.023.012-.023.01-.022.01-.023.008-.022.006-.022.006-.022.004-.022.004-.021.001-.021.001-.021v-4.153l-.077.054-.08.054-.083.053-.085.053-.087.053-.09.051-.093.051-.095.051-.097.05-.1.049-.102.048-.105.048-.106.048-.109.046-.111.046-.114.046-.115.044-.118.044-.12.043-.122.043-.124.042-.126.041-.128.04-.13.039-.132.039-.134.038-.135.037-.138.036-.139.036-.142.034-.143.034-.144.033-.147.032-.148.032-.15.03-.151.03-.153.028-.154.028-.156.027-.158.026-.159.024-.161.024-.162.023-.163.023-.165.021-.166.02-.167.019-.169.018-.169.017-.171.016-.173.015-.173.014-.175.013-.175.012-.177.01-.178.01-.179.009-.179.007-.181.006-.182.006-.182.004-.184.003-.184.001-.185.001-.185-.001-.184-.001-.184-.003-.182-.004-.182-.006-.181-.006-.179-.007-.179-.009-.178-.01-.176-.01-.176-.012-.175-.013-.173-.014-.172-.015-.171-.016-.17-.017-.169-.018-.167-.019-.166-.02-.165-.021-.163-.023-.162-.023-.161-.024-.159-.024-.157-.026-.156-.027-.155-.028-.153-.028-.151-.03-.15-.03-.148-.032-.146-.032-.145-.033-.143-.034-.141-.034-.14-.036-.137-.036-.136-.037-.134-.038-.132-.039-.13-.039-.128-.041-.126-.041-.124-.041-.122-.043-.12-.043-.117-.044-.116-.044-.113-.046-.112-.046-.109-.046-.106-.048-.105-.048-.102-.048-.1-.05-.097-.049-.095-.051-.093-.051-.09-.052-.087-.052-.085-.053-.083-.053-.08-.054-.077-.054v4.153zm8.74-8.179l-.257.004-.254.005-.25.008-.247.011-.244.012-.241.014-.237.016-.233.018-.231.021-.226.022-.224.023-.22.026-.216.027-.212.028-.21.031-.205.032-.202.033-.198.034-.194.036-.191.038-.187.038-.183.04-.179.041-.175.042-.172.043-.168.043-.163.045-.16.046-.155.046-.152.048-.148.048-.143.048-.139.049-.136.05-.131.05-.126.051-.123.051-.118.051-.114.052-.11.052-.106.052-.101.052-.096.052-.092.052-.088.052-.083.052-.079.052-.074.051-.07.052-.065.051-.06.05-.056.05-.051.05-.023.025-.023.024-.021.024-.02.025-.019.024-.018.024-.017.023-.015.024-.014.023-.013.023-.012.023-.01.023-.01.022-.008.022-.006.023-.006.021-.004.022-.004.021-.001.021-.001.021.001.021.001.021.004.021.004.022.006.021.006.023.008.022.01.022.01.023.012.023.013.023.014.023.015.024.017.023.018.024.019.024.02.025.021.024.023.024.023.025.051.05.056.05.06.05.065.051.07.052.074.051.079.052.083.052.088.052.092.052.096.052.101.052.106.052.11.052.114.052.118.051.123.051.126.051.131.05.136.05.139.049.143.048.148.048.152.048.155.046.16.046.163.045.168.043.172.043.175.042.179.041.183.04.187.038.191.038.194.036.198.034.202.033.205.032.21.031.212.028.216.027.22.026.224.023.226.022.231.021.233.018.237.016.241.014.244.012.247.011.25.008.254.005.257.004.26.001.26-.001.257-.004.254-.005.25-.008.247-.011.244-.012.241-.014.237-.016.233-.018.231-.021.226-.022.224-.023.22-.026.216-.027.212-.028.21-.031.205-.032.202-.033.198-.034.194-.036.191-.038.187-.038.183-.04.179-.041.175-.042.172-.043.168-.043.163-.045.16-.046.155-.046.152-.048.148-.048.143-.048.139-.049.136-.05.131-.05.126-.051.123-.051.118-.051.114-.052.11-.052.106-.052.101-.052.096-.052.092-.052.088-.052.083-.052.079-.052.074-.051.07-.052.065-.051.06-.05.056-.05.051-.05.023-.025.023-.024.021-.024.02-.025.019-.024.018-.024.017-.023.015-.024.014-.023.013-.023.012-.023.01-.023.01-.022.008-.022.006-.023.006-.021.004-.022.004-.021.001-.021.001-.021-.001-.021-.001-.021-.004-.021-.004-.022-.006-.021-.006-.023-.008-.022-.01-.022-.01-.023-.012-.023-.013-.023-.014-.023-.015-.024-.017-.023-.018-.024-.019-.024-.02-.025-.021-.024-.023-.024-.023-.025-.051-.05-.056-.05-.06-.05-.065-.051-.07-.052-.074-.051-.079-.052-.083-.052-.088-.052-.092-.052-.096-.052-.101-.052-.106-.052-.11-.052-.114-.052-.118-.051-.123-.051-.126-.051-.131-.05-.136-.05-.139-.049-.143-.048-.148-.048-.152-.048-.155-.046-.16-.046-.163-.045-.168-.043-.172-.043-.175-.042-.179-.041-.183-.04-.187-.038-.191-.038-.194-.036-.198-.034-.202-.033-.205-.032-.21-.031-.212-.028-.216-.027-.22-.026-.224-.023-.226-.022-.231-.021-.233-.018-.237-.016-.241-.014-.244-.012-.247-.011-.25-.008-.254-.005-.257-.004-.26-.001-.26.001z"
  );
}, "insertDatabaseIcon");
var insertComputerIcon = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("symbol").attr("id", "computer").attr("width", "24").attr("height", "24").append("path").attr("transform", "scale(.5)").attr(
    "d",
    "M2 2v13h20v-13h-20zm18 11h-16v-9h16v9zm-10.228 6l.466-1h3.524l.467 1h-4.457zm14.228 3h-24l2-6h2.104l-1.33 4h18.45l-1.297-4h2.073l2 6zm-5-10h-14v-7h14v7z"
  );
}, "insertComputerIcon");
var insertClockIcon = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("symbol").attr("id", "clock").attr("width", "24").attr("height", "24").append("path").attr("transform", "scale(.5)").attr(
    "d",
    "M12 2c5.514 0 10 4.486 10 10s-4.486 10-10 10-10-4.486-10-10 4.486-10 10-10zm0-2c-6.627 0-12 5.373-12 12s5.373 12 12 12 12-5.373 12-12-5.373-12-12-12zm5.848 12.459c.202.038.202.333.001.372-1.907.361-6.045 1.111-6.547 1.111-.719 0-1.301-.582-1.301-1.301 0-.512.77-5.447 1.125-7.445.034-.192.312-.181.343.014l.985 6.238 5.394 1.011z"
  );
}, "insertClockIcon");
var insertArrowHead = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("marker").attr("id", "arrowhead").attr("refX", 7.9).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 12).attr("markerHeight", 12).attr("orient", "auto-start-reverse").append("path").attr("d", "M -1 0 L 10 5 L 0 10 z");
}, "insertArrowHead");
var insertArrowFilledHead = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("marker").attr("id", "filled-head").attr("refX", 15.5).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L14,7 L9,1 Z");
}, "insertArrowFilledHead");
var insertSequenceNumber = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  elem.append("defs").append("marker").attr("id", "sequencenumber").attr("refX", 15).attr("refY", 15).attr("markerWidth", 60).attr("markerHeight", 40).attr("orient", "auto").append("circle").attr("cx", 15).attr("cy", 15).attr("r", 6);
}, "insertSequenceNumber");
var insertArrowCrossHead = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem) {
  const defs = elem.append("defs");
  const marker = defs.append("marker").attr("id", "crosshead").attr("markerWidth", 15).attr("markerHeight", 8).attr("orient", "auto").attr("refX", 4).attr("refY", 4.5);
  marker.append("path").attr("fill", "none").attr("stroke", "#000000").style("stroke-dasharray", "0, 0").attr("stroke-width", "1pt").attr("d", "M 1,2 L 6,7 M 6,2 L 1,7");
}, "insertArrowCrossHead");
var getTextObj2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  return {
    x: 0,
    y: 0,
    fill: void 0,
    anchor: void 0,
    style: "#666",
    width: void 0,
    height: void 0,
    textMargin: 0,
    rx: 0,
    ry: 0,
    tspan: true,
    valign: void 0
  };
}, "getTextObj");
var getNoteRect2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  return {
    x: 0,
    y: 0,
    fill: "#EDF2AE",
    stroke: "#666",
    width: 100,
    anchor: "start",
    height: 100,
    rx: 0,
    ry: 0
  };
}, "getNoteRect");
var _drawTextCandidateFunc = /* @__PURE__ */ function() {
  function byText(content, g, x, y, width, height, textAttrs) {
    const text = g.append("text").attr("x", x + width / 2).attr("y", y + height / 2 + 5).style("text-anchor", "middle").text(content);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byText, "byText");
  function byTspan(content, g, x, y, width, height, textAttrs, conf2) {
    const { actorFontSize, actorFontFamily, actorFontWeight } = conf2;
    const [_actorFontSize, _actorFontSizePx] = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .parseFontSize */ .VG)(actorFontSize);
    const lines = content.split(_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.lineBreakRegex);
    for (let i = 0; i < lines.length; i++) {
      const dy = i * _actorFontSize - _actorFontSize * (lines.length - 1) / 2;
      const text = g.append("text").attr("x", x + width / 2).attr("y", y).style("text-anchor", "middle").style("font-size", _actorFontSizePx).style("font-weight", actorFontWeight).style("font-family", actorFontFamily);
      text.append("tspan").attr("x", x + width / 2).attr("dy", dy).text(lines[i]);
      text.attr("y", y + height / 2).attr("dominant-baseline", "central").attr("alignment-baseline", "central");
      _setTextAttrs(text, textAttrs);
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byTspan, "byTspan");
  function byFo(content, g, x, y, width, height, textAttrs, conf2) {
    const s = g.append("switch");
    const f = s.append("foreignObject").attr("x", x).attr("y", y).attr("width", width).attr("height", height);
    const text = f.append("xhtml:div").style("display", "table").style("height", "100%").style("width", "100%");
    text.append("div").style("display", "table-cell").style("text-align", "center").style("vertical-align", "middle").text(content);
    byTspan(content, s, x, y, width, height, textAttrs, conf2);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byFo, "byFo");
  async function byKatex(content, g, x, y, width, height, textAttrs, conf2) {
    const dim = await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateMathMLDimensions */ .nH)(content, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)());
    const s = g.append("switch");
    const f = s.append("foreignObject").attr("x", x + width / 2 - dim.width / 2).attr("y", y + height / 2 - dim.height / 2).attr("width", dim.width).attr("height", dim.height);
    const text = f.append("xhtml:div").style("height", "100%").style("width", "100%");
    text.append("div").style("text-align", "center").style("vertical-align", "middle").html(await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .renderKatex */ .uT)(content, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)()));
    byTspan(content, s, x, y, width, height, textAttrs, conf2);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byKatex, "byKatex");
  function _setTextAttrs(toText, fromTextAttrsDict) {
    for (const key in fromTextAttrsDict) {
      if (fromTextAttrsDict.hasOwnProperty(key)) {
        toText.attr(key, fromTextAttrsDict[key]);
      }
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(_setTextAttrs, "_setTextAttrs");
  return function(conf2, hasKatex2 = false) {
    if (hasKatex2) {
      return byKatex;
    }
    return conf2.textPlacement === "fo" ? byFo : conf2.textPlacement === "old" ? byText : byTspan;
  };
}();
var _drawMenuItemTextCandidateFunc = /* @__PURE__ */ function() {
  function byText(content, g, x, y, width, height, textAttrs) {
    const text = g.append("text").attr("x", x).attr("y", y).style("text-anchor", "start").text(content);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byText, "byText");
  function byTspan(content, g, x, y, width, height, textAttrs, conf2) {
    const { actorFontSize, actorFontFamily, actorFontWeight } = conf2;
    const lines = content.split(_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.lineBreakRegex);
    for (let i = 0; i < lines.length; i++) {
      const dy = i * actorFontSize - actorFontSize * (lines.length - 1) / 2;
      const text = g.append("text").attr("x", x).attr("y", y).style("text-anchor", "start").style("font-size", actorFontSize).style("font-weight", actorFontWeight).style("font-family", actorFontFamily);
      text.append("tspan").attr("x", x).attr("dy", dy).text(lines[i]);
      text.attr("y", y + height / 2).attr("dominant-baseline", "central").attr("alignment-baseline", "central");
      _setTextAttrs(text, textAttrs);
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byTspan, "byTspan");
  function byFo(content, g, x, y, width, height, textAttrs, conf2) {
    const s = g.append("switch");
    const f = s.append("foreignObject").attr("x", x).attr("y", y).attr("width", width).attr("height", height);
    const text = f.append("xhtml:div").style("display", "table").style("height", "100%").style("width", "100%");
    text.append("div").style("display", "table-cell").style("text-align", "center").style("vertical-align", "middle").text(content);
    byTspan(content, s, x, y, width, height, textAttrs, conf2);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byFo, "byFo");
  function _setTextAttrs(toText, fromTextAttrsDict) {
    for (const key in fromTextAttrsDict) {
      if (fromTextAttrsDict.hasOwnProperty(key)) {
        toText.attr(key, fromTextAttrsDict[key]);
      }
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(_setTextAttrs, "_setTextAttrs");
  return function(conf2) {
    return conf2.textPlacement === "fo" ? byFo : conf2.textPlacement === "old" ? byText : byTspan;
  };
}();
var svgDraw_default = {
  drawRect: drawRect2,
  drawText,
  drawLabel,
  drawActor,
  drawBox,
  drawPopup,
  anchorElement,
  drawActivation,
  drawLoop,
  drawBackgroundRect: drawBackgroundRect2,
  insertArrowHead,
  insertArrowFilledHead,
  insertSequenceNumber,
  insertArrowCrossHead,
  insertDatabaseIcon,
  insertComputerIcon,
  insertClockIcon,
  getTextObj: getTextObj2,
  getNoteRect: getNoteRect2,
  fixLifeLineHeights,
  sanitizeUrl: _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_5__/* .sanitizeUrl */ .N
};

// src/diagrams/sequence/sequenceRenderer.ts
var conf = {};
var bounds = {
  data: {
    startx: void 0,
    stopx: void 0,
    starty: void 0,
    stopy: void 0
  },
  verticalPos: 0,
  sequenceItems: [],
  activations: [],
  models: {
    getHeight: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      return Math.max.apply(
        null,
        this.actors.length === 0 ? [0] : this.actors.map((actor) => actor.height || 0)
      ) + (this.loops.length === 0 ? 0 : this.loops.map((it) => it.height || 0).reduce((acc, h) => acc + h)) + (this.messages.length === 0 ? 0 : this.messages.map((it) => it.height || 0).reduce((acc, h) => acc + h)) + (this.notes.length === 0 ? 0 : this.notes.map((it) => it.height || 0).reduce((acc, h) => acc + h));
    }, "getHeight"),
    clear: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      this.actors = [];
      this.boxes = [];
      this.loops = [];
      this.messages = [];
      this.notes = [];
    }, "clear"),
    addBox: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(boxModel) {
      this.boxes.push(boxModel);
    }, "addBox"),
    addActor: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(actorModel) {
      this.actors.push(actorModel);
    }, "addActor"),
    addLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(loopModel) {
      this.loops.push(loopModel);
    }, "addLoop"),
    addMessage: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(msgModel) {
      this.messages.push(msgModel);
    }, "addMessage"),
    addNote: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(noteModel) {
      this.notes.push(noteModel);
    }, "addNote"),
    lastActor: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      return this.actors[this.actors.length - 1];
    }, "lastActor"),
    lastLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      return this.loops[this.loops.length - 1];
    }, "lastLoop"),
    lastMessage: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      return this.messages[this.messages.length - 1];
    }, "lastMessage"),
    lastNote: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
      return this.notes[this.notes.length - 1];
    }, "lastNote"),
    actors: [],
    boxes: [],
    loops: [],
    messages: [],
    notes: []
  },
  init: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    this.sequenceItems = [];
    this.activations = [];
    this.models.clear();
    this.data = {
      startx: void 0,
      stopx: void 0,
      starty: void 0,
      stopy: void 0
    };
    this.verticalPos = 0;
    setConf((0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
  }, "init"),
  updateVal: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(obj, key, val, fun) {
    if (obj[key] === void 0) {
      obj[key] = val;
    } else {
      obj[key] = fun(val, obj[key]);
    }
  }, "updateVal"),
  updateBounds: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(startx, starty, stopx, stopy) {
    const _self = this;
    let cnt = 0;
    function updateFn(type) {
      return /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function updateItemBounds(item) {
        cnt++;
        const n = _self.sequenceItems.length - cnt + 1;
        _self.updateVal(item, "starty", starty - n * conf.boxMargin, Math.min);
        _self.updateVal(item, "stopy", stopy + n * conf.boxMargin, Math.max);
        _self.updateVal(bounds.data, "startx", startx - n * conf.boxMargin, Math.min);
        _self.updateVal(bounds.data, "stopx", stopx + n * conf.boxMargin, Math.max);
        if (!(type === "activation")) {
          _self.updateVal(item, "startx", startx - n * conf.boxMargin, Math.min);
          _self.updateVal(item, "stopx", stopx + n * conf.boxMargin, Math.max);
          _self.updateVal(bounds.data, "starty", starty - n * conf.boxMargin, Math.min);
          _self.updateVal(bounds.data, "stopy", stopy + n * conf.boxMargin, Math.max);
        }
      }, "updateItemBounds");
    }
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(updateFn, "updateFn");
    this.sequenceItems.forEach(updateFn());
    this.activations.forEach(updateFn("activation"));
  }, "updateBounds"),
  insert: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(startx, starty, stopx, stopy) {
    const _startx = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(startx, stopx);
    const _stopx = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(startx, stopx);
    const _starty = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(starty, stopy);
    const _stopy = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(starty, stopy);
    this.updateVal(bounds.data, "startx", _startx, Math.min);
    this.updateVal(bounds.data, "starty", _starty, Math.min);
    this.updateVal(bounds.data, "stopx", _stopx, Math.max);
    this.updateVal(bounds.data, "stopy", _stopy, Math.max);
    this.updateBounds(_startx, _starty, _stopx, _stopy);
  }, "insert"),
  newActivation: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(message, diagram2, actors) {
    const actorRect = actors.get(message.from);
    const stackedSize = actorActivations(message.from).length || 0;
    const x = actorRect.x + actorRect.width / 2 + (stackedSize - 1) * conf.activationWidth / 2;
    this.activations.push({
      startx: x,
      starty: this.verticalPos + 2,
      stopx: x + conf.activationWidth,
      stopy: void 0,
      actor: message.from,
      anchored: svgDraw_default.anchorElement(diagram2)
    });
  }, "newActivation"),
  endActivation: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(message) {
    const lastActorActivationIdx = this.activations.map(function(activation) {
      return activation.actor;
    }).lastIndexOf(message.from);
    return this.activations.splice(lastActorActivationIdx, 1)[0];
  }, "endActivation"),
  createLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(title = { message: void 0, wrap: false, width: void 0 }, fill) {
    return {
      startx: void 0,
      starty: this.verticalPos,
      stopx: void 0,
      stopy: void 0,
      title: title.message,
      wrap: title.wrap,
      width: title.width,
      height: 0,
      fill
    };
  }, "createLoop"),
  newLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(title = { message: void 0, wrap: false, width: void 0 }, fill) {
    this.sequenceItems.push(this.createLoop(title, fill));
  }, "newLoop"),
  endLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return this.sequenceItems.pop();
  }, "endLoop"),
  isLoopOverlap: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return this.sequenceItems.length ? this.sequenceItems[this.sequenceItems.length - 1].overlap : false;
  }, "isLoopOverlap"),
  addSectionToLoop: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(message) {
    const loop = this.sequenceItems.pop();
    loop.sections = loop.sections || [];
    loop.sectionTitles = loop.sectionTitles || [];
    loop.sections.push({ y: bounds.getVerticalPos(), height: 0 });
    loop.sectionTitles.push(message);
    this.sequenceItems.push(loop);
  }, "addSectionToLoop"),
  saveVerticalPos: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    if (this.isLoopOverlap()) {
      this.savedVerticalPos = this.verticalPos;
    }
  }, "saveVerticalPos"),
  resetVerticalPos: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    if (this.isLoopOverlap()) {
      this.verticalPos = this.savedVerticalPos;
    }
  }, "resetVerticalPos"),
  bumpVerticalPos: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(bump) {
    this.verticalPos = this.verticalPos + bump;
    this.data.stopy = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(this.data.stopy, this.verticalPos);
  }, "bumpVerticalPos"),
  getVerticalPos: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return this.verticalPos;
  }, "getVerticalPos"),
  getBounds: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return { bounds: this.data, models: this.models };
  }, "getBounds")
};
var drawNote = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(elem, noteModel) {
  bounds.bumpVerticalPos(conf.boxMargin);
  noteModel.height = conf.boxMargin;
  noteModel.starty = bounds.getVerticalPos();
  const rect = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  rect.x = noteModel.startx;
  rect.y = noteModel.starty;
  rect.width = noteModel.width || conf.width;
  rect.class = "note";
  const g = elem.append("g");
  const rectElem = svgDraw_default.drawRect(g, rect);
  const textObj = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getTextObj */ .AD)();
  textObj.x = noteModel.startx;
  textObj.y = noteModel.starty;
  textObj.width = rect.width;
  textObj.dy = "1em";
  textObj.text = noteModel.message;
  textObj.class = "noteText";
  textObj.fontFamily = conf.noteFontFamily;
  textObj.fontSize = conf.noteFontSize;
  textObj.fontWeight = conf.noteFontWeight;
  textObj.anchor = conf.noteAlign;
  textObj.textMargin = conf.noteMargin;
  textObj.valign = "center";
  const textElem = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(textObj.text) ? await drawKatex(g, textObj) : drawText(g, textObj);
  const textHeight = Math.round(
    textElem.map((te) => (te._groups || te)[0][0].getBBox().height).reduce((acc, curr) => acc + curr)
  );
  rectElem.attr("height", textHeight + 2 * conf.noteMargin);
  noteModel.height += textHeight + 2 * conf.noteMargin;
  bounds.bumpVerticalPos(textHeight + 2 * conf.noteMargin);
  noteModel.stopy = noteModel.starty + textHeight + 2 * conf.noteMargin;
  noteModel.stopx = noteModel.startx + rect.width;
  bounds.insert(noteModel.startx, noteModel.starty, noteModel.stopx, noteModel.stopy);
  bounds.models.addNote(noteModel);
}, "drawNote");
var messageFont = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((cnf) => {
  return {
    fontFamily: cnf.messageFontFamily,
    fontSize: cnf.messageFontSize,
    fontWeight: cnf.messageFontWeight
  };
}, "messageFont");
var noteFont = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((cnf) => {
  return {
    fontFamily: cnf.noteFontFamily,
    fontSize: cnf.noteFontSize,
    fontWeight: cnf.noteFontWeight
  };
}, "noteFont");
var actorFont = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((cnf) => {
  return {
    fontFamily: cnf.actorFontFamily,
    fontSize: cnf.actorFontSize,
    fontWeight: cnf.actorFontWeight
  };
}, "actorFont");
async function boundMessage(_diagram, msgModel) {
  bounds.bumpVerticalPos(10);
  const { startx, stopx, message } = msgModel;
  const lines = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.splitBreaks(message).length;
  const isKatexMsg = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(message);
  const textDims = isKatexMsg ? await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateMathMLDimensions */ .nH)(message, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()) : _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(message, messageFont(conf));
  if (!isKatexMsg) {
    const lineHeight = textDims.height / lines;
    msgModel.height += lineHeight;
    bounds.bumpVerticalPos(lineHeight);
  }
  let lineStartY;
  let totalOffset = textDims.height - 10;
  const textWidth = textDims.width;
  if (startx === stopx) {
    lineStartY = bounds.getVerticalPos() + totalOffset;
    if (!conf.rightAngles) {
      totalOffset += conf.boxMargin;
      lineStartY = bounds.getVerticalPos() + totalOffset;
    }
    totalOffset += 30;
    const dx = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(textWidth / 2, conf.width / 2);
    bounds.insert(
      startx - dx,
      bounds.getVerticalPos() - 10 + totalOffset,
      stopx + dx,
      bounds.getVerticalPos() + 30 + totalOffset
    );
  } else {
    totalOffset += conf.boxMargin;
    lineStartY = bounds.getVerticalPos() + totalOffset;
    bounds.insert(startx, lineStartY - 10, stopx, lineStartY);
  }
  bounds.bumpVerticalPos(totalOffset);
  msgModel.height += totalOffset;
  msgModel.stopy = msgModel.starty + msgModel.height;
  bounds.insert(msgModel.fromBounds, msgModel.starty, msgModel.toBounds, msgModel.stopy);
  return lineStartY;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(boundMessage, "boundMessage");
var drawMessage = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(diagram2, msgModel, lineStartY, diagObj) {
  const { startx, stopx, starty, message, type, sequenceIndex, sequenceVisible } = msgModel;
  const textDims = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(message, messageFont(conf));
  const textObj = (0,_chunk_D6G4REZN_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getTextObj */ .AD)();
  textObj.x = startx;
  textObj.y = starty + 10;
  textObj.width = stopx - startx;
  textObj.class = "messageText";
  textObj.dy = "1em";
  textObj.text = message;
  textObj.fontFamily = conf.messageFontFamily;
  textObj.fontSize = conf.messageFontSize;
  textObj.fontWeight = conf.messageFontWeight;
  textObj.anchor = conf.messageAlign;
  textObj.valign = "center";
  textObj.textMargin = conf.wrapPadding;
  textObj.tspan = false;
  if ((0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(textObj.text)) {
    await drawKatex(diagram2, textObj, { startx, stopx, starty: lineStartY });
  } else {
    drawText(diagram2, textObj);
  }
  const textWidth = textDims.width;
  let line;
  if (startx === stopx) {
    if (conf.rightAngles) {
      line = diagram2.append("path").attr(
        "d",
        `M  ${startx},${lineStartY} H ${startx + _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width / 2, textWidth / 2)} V ${lineStartY + 25} H ${startx}`
      );
    } else {
      line = diagram2.append("path").attr(
        "d",
        "M " + startx + "," + lineStartY + " C " + (startx + 60) + "," + (lineStartY - 10) + " " + (startx + 60) + "," + (lineStartY + 30) + " " + startx + "," + (lineStartY + 20)
      );
    }
  } else {
    line = diagram2.append("line");
    line.attr("x1", startx);
    line.attr("y1", lineStartY);
    line.attr("x2", stopx);
    line.attr("y2", lineStartY);
  }
  if (type === diagObj.db.LINETYPE.DOTTED || type === diagObj.db.LINETYPE.DOTTED_CROSS || type === diagObj.db.LINETYPE.DOTTED_POINT || type === diagObj.db.LINETYPE.DOTTED_OPEN || type === diagObj.db.LINETYPE.BIDIRECTIONAL_DOTTED) {
    line.style("stroke-dasharray", "3, 3");
    line.attr("class", "messageLine1");
  } else {
    line.attr("class", "messageLine0");
  }
  let url = "";
  if (conf.arrowMarkerAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = url.replace(/\(/g, "\\(");
    url = url.replace(/\)/g, "\\)");
  }
  line.attr("stroke-width", 2);
  line.attr("stroke", "none");
  line.style("fill", "none");
  if (type === diagObj.db.LINETYPE.SOLID || type === diagObj.db.LINETYPE.DOTTED) {
    line.attr("marker-end", "url(" + url + "#arrowhead)");
  }
  if (type === diagObj.db.LINETYPE.BIDIRECTIONAL_SOLID || type === diagObj.db.LINETYPE.BIDIRECTIONAL_DOTTED) {
    line.attr("marker-start", "url(" + url + "#arrowhead)");
    line.attr("marker-end", "url(" + url + "#arrowhead)");
  }
  if (type === diagObj.db.LINETYPE.SOLID_POINT || type === diagObj.db.LINETYPE.DOTTED_POINT) {
    line.attr("marker-end", "url(" + url + "#filled-head)");
  }
  if (type === diagObj.db.LINETYPE.SOLID_CROSS || type === diagObj.db.LINETYPE.DOTTED_CROSS) {
    line.attr("marker-end", "url(" + url + "#crosshead)");
  }
  if (sequenceVisible || conf.showSequenceNumbers) {
    line.attr("marker-start", "url(" + url + "#sequencenumber)");
    diagram2.append("text").attr("x", startx).attr("y", lineStartY + 4).attr("font-family", "sans-serif").attr("font-size", "12px").attr("text-anchor", "middle").attr("class", "sequenceNumber").text(sequenceIndex);
  }
}, "drawMessage");
var addActorRenderingData = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(diagram2, actors, createdActors, actorKeys, verticalPos, messages, isFooter) {
  let prevWidth = 0;
  let prevMargin = 0;
  let prevBox = void 0;
  let maxHeight = 0;
  for (const actorKey of actorKeys) {
    const actor = actors.get(actorKey);
    const box = actor.box;
    if (prevBox && prevBox != box) {
      if (!isFooter) {
        bounds.models.addBox(prevBox);
      }
      prevMargin += conf.boxMargin + prevBox.margin;
    }
    if (box && box != prevBox) {
      if (!isFooter) {
        box.x = prevWidth + prevMargin;
        box.y = verticalPos;
      }
      prevMargin += box.margin;
    }
    actor.width = actor.width || conf.width;
    actor.height = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(actor.height || conf.height, conf.height);
    actor.margin = actor.margin || conf.actorMargin;
    maxHeight = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(maxHeight, actor.height);
    if (createdActors.get(actor.name)) {
      prevMargin += actor.width / 2;
    }
    actor.x = prevWidth + prevMargin;
    actor.starty = bounds.getVerticalPos();
    bounds.insert(actor.x, verticalPos, actor.x + actor.width, actor.height);
    prevWidth += actor.width + prevMargin;
    if (actor.box) {
      actor.box.width = prevWidth + box.margin - actor.box.x;
    }
    prevMargin = actor.margin;
    prevBox = actor.box;
    bounds.models.addActor(actor);
  }
  if (prevBox && !isFooter) {
    bounds.models.addBox(prevBox);
  }
  bounds.bumpVerticalPos(maxHeight);
}, "addActorRenderingData");
var drawActors = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(diagram2, actors, actorKeys, isFooter) {
  if (!isFooter) {
    for (const actorKey of actorKeys) {
      const actor = actors.get(actorKey);
      await svgDraw_default.drawActor(diagram2, actor, conf, false);
    }
  } else {
    let maxHeight = 0;
    bounds.bumpVerticalPos(conf.boxMargin * 2);
    for (const actorKey of actorKeys) {
      const actor = actors.get(actorKey);
      if (!actor.stopy) {
        actor.stopy = bounds.getVerticalPos();
      }
      const height = await svgDraw_default.drawActor(diagram2, actor, conf, true);
      maxHeight = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(maxHeight, height);
    }
    bounds.bumpVerticalPos(maxHeight + conf.boxMargin);
  }
}, "drawActors");
var drawActorsPopup = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(diagram2, actors, actorKeys, doc) {
  let maxHeight = 0;
  let maxWidth = 0;
  for (const actorKey of actorKeys) {
    const actor = actors.get(actorKey);
    const minMenuWidth = getRequiredPopupWidth(actor);
    const menuDimensions = svgDraw_default.drawPopup(
      diagram2,
      actor,
      minMenuWidth,
      conf,
      conf.forceMenus,
      doc
    );
    if (menuDimensions.height > maxHeight) {
      maxHeight = menuDimensions.height;
    }
    if (menuDimensions.width + actor.x > maxWidth) {
      maxWidth = menuDimensions.width + actor.x;
    }
  }
  return { maxHeight, maxWidth };
}, "drawActorsPopup");
var setConf = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(cnf) {
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .assignWithDepth_default */ .Yc)(conf, cnf);
  if (cnf.fontFamily) {
    conf.actorFontFamily = conf.noteFontFamily = conf.messageFontFamily = cnf.fontFamily;
  }
  if (cnf.fontSize) {
    conf.actorFontSize = conf.noteFontSize = conf.messageFontSize = cnf.fontSize;
  }
  if (cnf.fontWeight) {
    conf.actorFontWeight = conf.noteFontWeight = conf.messageFontWeight = cnf.fontWeight;
  }
}, "setConf");
var actorActivations = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(actor) {
  return bounds.activations.filter(function(activation) {
    return activation.actor === actor;
  });
}, "actorActivations");
var activationBounds = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(actor, actors) {
  const actorObj = actors.get(actor);
  const activations = actorActivations(actor);
  const left = activations.reduce(
    function(acc, activation) {
      return _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(acc, activation.startx);
    },
    actorObj.x + actorObj.width / 2 - 1
  );
  const right = activations.reduce(
    function(acc, activation) {
      return _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(acc, activation.stopx);
    },
    actorObj.x + actorObj.width / 2 + 1
  );
  return [left, right];
}, "activationBounds");
function adjustLoopHeightForWrap(loopWidths, msg, preMargin, postMargin, addLoopFn) {
  bounds.bumpVerticalPos(preMargin);
  let heightAdjust = postMargin;
  if (msg.id && msg.message && loopWidths[msg.id]) {
    const loopWidth = loopWidths[msg.id].width;
    const textConf = messageFont(conf);
    msg.message = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(`[${msg.message}]`, loopWidth - 2 * conf.wrapPadding, textConf);
    msg.width = loopWidth;
    msg.wrap = true;
    const textDims = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(msg.message, textConf);
    const totalOffset = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(textDims.height, conf.labelBoxHeight);
    heightAdjust = postMargin + totalOffset;
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(`${totalOffset} - ${msg.message}`);
  }
  addLoopFn(msg);
  bounds.bumpVerticalPos(heightAdjust);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(adjustLoopHeightForWrap, "adjustLoopHeightForWrap");
function adjustCreatedDestroyedData(msg, msgModel, lineStartY, index, actors, createdActors, destroyedActors) {
  function receiverAdjustment(actor, adjustment) {
    if (actor.x < actors.get(msg.from).x) {
      bounds.insert(
        msgModel.stopx - adjustment,
        msgModel.starty,
        msgModel.startx,
        msgModel.stopy + actor.height / 2 + conf.noteMargin
      );
      msgModel.stopx = msgModel.stopx + adjustment;
    } else {
      bounds.insert(
        msgModel.startx,
        msgModel.starty,
        msgModel.stopx + adjustment,
        msgModel.stopy + actor.height / 2 + conf.noteMargin
      );
      msgModel.stopx = msgModel.stopx - adjustment;
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(receiverAdjustment, "receiverAdjustment");
  function senderAdjustment(actor, adjustment) {
    if (actor.x < actors.get(msg.to).x) {
      bounds.insert(
        msgModel.startx - adjustment,
        msgModel.starty,
        msgModel.stopx,
        msgModel.stopy + actor.height / 2 + conf.noteMargin
      );
      msgModel.startx = msgModel.startx + adjustment;
    } else {
      bounds.insert(
        msgModel.stopx,
        msgModel.starty,
        msgModel.startx + adjustment,
        msgModel.stopy + actor.height / 2 + conf.noteMargin
      );
      msgModel.startx = msgModel.startx - adjustment;
    }
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(senderAdjustment, "senderAdjustment");
  if (createdActors.get(msg.to) == index) {
    const actor = actors.get(msg.to);
    const adjustment = actor.type == "actor" ? ACTOR_TYPE_WIDTH / 2 + 3 : actor.width / 2 + 3;
    receiverAdjustment(actor, adjustment);
    actor.starty = lineStartY - actor.height / 2;
    bounds.bumpVerticalPos(actor.height / 2);
  } else if (destroyedActors.get(msg.from) == index) {
    const actor = actors.get(msg.from);
    if (conf.mirrorActors) {
      const adjustment = actor.type == "actor" ? ACTOR_TYPE_WIDTH / 2 : actor.width / 2;
      senderAdjustment(actor, adjustment);
    }
    actor.stopy = lineStartY - actor.height / 2;
    bounds.bumpVerticalPos(actor.height / 2);
  } else if (destroyedActors.get(msg.to) == index) {
    const actor = actors.get(msg.to);
    if (conf.mirrorActors) {
      const adjustment = actor.type == "actor" ? ACTOR_TYPE_WIDTH / 2 + 3 : actor.width / 2 + 3;
      receiverAdjustment(actor, adjustment);
    }
    actor.stopy = lineStartY - actor.height / 2;
    bounds.bumpVerticalPos(actor.height / 2);
  }
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(adjustCreatedDestroyedData, "adjustCreatedDestroyedData");
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(_text, id, _version, diagObj) {
  const { securityLevel, sequence } = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
  conf = sequence;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)("body");
  const doc = securityLevel === "sandbox" ? sandboxElement.nodes()[0].contentDocument : document;
  bounds.init();
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(diagObj.db);
  const diagram2 = securityLevel === "sandbox" ? root.select(`[id="${id}"]`) : (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(`[id="${id}"]`);
  const actors = diagObj.db.getActors();
  const createdActors = diagObj.db.getCreatedActors();
  const destroyedActors = diagObj.db.getDestroyedActors();
  const boxes = diagObj.db.getBoxes();
  let actorKeys = diagObj.db.getActorKeys();
  const messages = diagObj.db.getMessages();
  const title = diagObj.db.getDiagramTitle();
  const hasBoxes = diagObj.db.hasAtLeastOneBox();
  const hasBoxTitles = diagObj.db.hasAtLeastOneBoxWithTitle();
  const maxMessageWidthPerActor = await getMaxMessageWidthPerActor(actors, messages, diagObj);
  conf.height = await calculateActorMargins(actors, maxMessageWidthPerActor, boxes);
  svgDraw_default.insertComputerIcon(diagram2);
  svgDraw_default.insertDatabaseIcon(diagram2);
  svgDraw_default.insertClockIcon(diagram2);
  if (hasBoxes) {
    bounds.bumpVerticalPos(conf.boxMargin);
    if (hasBoxTitles) {
      bounds.bumpVerticalPos(boxes[0].textMaxHeight);
    }
  }
  if (conf.hideUnusedParticipants === true) {
    const newActors = /* @__PURE__ */ new Set();
    messages.forEach((message) => {
      newActors.add(message.from);
      newActors.add(message.to);
    });
    actorKeys = actorKeys.filter((actorKey) => newActors.has(actorKey));
  }
  addActorRenderingData(diagram2, actors, createdActors, actorKeys, 0, messages, false);
  const loopWidths = await calculateLoopBounds(messages, actors, maxMessageWidthPerActor, diagObj);
  svgDraw_default.insertArrowHead(diagram2);
  svgDraw_default.insertArrowCrossHead(diagram2);
  svgDraw_default.insertArrowFilledHead(diagram2);
  svgDraw_default.insertSequenceNumber(diagram2);
  function activeEnd(msg, verticalPos) {
    const activationData = bounds.endActivation(msg);
    if (activationData.starty + 18 > verticalPos) {
      activationData.starty = verticalPos - 6;
      verticalPos += 12;
    }
    svgDraw_default.drawActivation(
      diagram2,
      activationData,
      verticalPos,
      conf,
      actorActivations(msg.from).length
    );
    bounds.insert(activationData.startx, verticalPos - 10, activationData.stopx, verticalPos);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(activeEnd, "activeEnd");
  let sequenceIndex = 1;
  let sequenceIndexStep = 1;
  const messagesToDraw = [];
  const backgrounds = [];
  let index = 0;
  for (const msg of messages) {
    let loopModel, noteModel, msgModel;
    switch (msg.type) {
      case diagObj.db.LINETYPE.NOTE:
        bounds.resetVerticalPos();
        noteModel = msg.noteModel;
        await drawNote(diagram2, noteModel);
        break;
      case diagObj.db.LINETYPE.ACTIVE_START:
        bounds.newActivation(msg, diagram2, actors);
        break;
      case diagObj.db.LINETYPE.ACTIVE_END:
        activeEnd(msg, bounds.getVerticalPos());
        break;
      case diagObj.db.LINETYPE.LOOP_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.LOOP_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "loop", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      case diagObj.db.LINETYPE.RECT_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin,
          (message) => bounds.newLoop(void 0, message.message)
        );
        break;
      case diagObj.db.LINETYPE.RECT_END:
        loopModel = bounds.endLoop();
        backgrounds.push(loopModel);
        bounds.models.addLoop(loopModel);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        break;
      case diagObj.db.LINETYPE.OPT_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.OPT_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "opt", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      case diagObj.db.LINETYPE.ALT_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.ALT_ELSE:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin + conf.boxTextMargin,
          conf.boxMargin,
          (message) => bounds.addSectionToLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.ALT_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "alt", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      case diagObj.db.LINETYPE.PAR_START:
      case diagObj.db.LINETYPE.PAR_OVER_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        bounds.saveVerticalPos();
        break;
      case diagObj.db.LINETYPE.PAR_AND:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin + conf.boxTextMargin,
          conf.boxMargin,
          (message) => bounds.addSectionToLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.PAR_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "par", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      case diagObj.db.LINETYPE.AUTONUMBER:
        sequenceIndex = msg.message.start || sequenceIndex;
        sequenceIndexStep = msg.message.step || sequenceIndexStep;
        if (msg.message.visible) {
          diagObj.db.enableSequenceNumbers();
        } else {
          diagObj.db.disableSequenceNumbers();
        }
        break;
      case diagObj.db.LINETYPE.CRITICAL_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.CRITICAL_OPTION:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin + conf.boxTextMargin,
          conf.boxMargin,
          (message) => bounds.addSectionToLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.CRITICAL_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "critical", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      case diagObj.db.LINETYPE.BREAK_START:
        adjustLoopHeightForWrap(
          loopWidths,
          msg,
          conf.boxMargin,
          conf.boxMargin + conf.boxTextMargin,
          (message) => bounds.newLoop(message)
        );
        break;
      case diagObj.db.LINETYPE.BREAK_END:
        loopModel = bounds.endLoop();
        await svgDraw_default.drawLoop(diagram2, loopModel, "break", conf);
        bounds.bumpVerticalPos(loopModel.stopy - bounds.getVerticalPos());
        bounds.models.addLoop(loopModel);
        break;
      default:
        try {
          msgModel = msg.msgModel;
          msgModel.starty = bounds.getVerticalPos();
          msgModel.sequenceIndex = sequenceIndex;
          msgModel.sequenceVisible = diagObj.db.showSequenceNumbers();
          const lineStartY = await boundMessage(diagram2, msgModel);
          adjustCreatedDestroyedData(
            msg,
            msgModel,
            lineStartY,
            index,
            actors,
            createdActors,
            destroyedActors
          );
          messagesToDraw.push({ messageModel: msgModel, lineStartY });
          bounds.models.addMessage(msgModel);
        } catch (e) {
          _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.error("error while drawing message", e);
        }
    }
    if ([
      diagObj.db.LINETYPE.SOLID_OPEN,
      diagObj.db.LINETYPE.DOTTED_OPEN,
      diagObj.db.LINETYPE.SOLID,
      diagObj.db.LINETYPE.DOTTED,
      diagObj.db.LINETYPE.SOLID_CROSS,
      diagObj.db.LINETYPE.DOTTED_CROSS,
      diagObj.db.LINETYPE.SOLID_POINT,
      diagObj.db.LINETYPE.DOTTED_POINT,
      diagObj.db.LINETYPE.BIDIRECTIONAL_SOLID,
      diagObj.db.LINETYPE.BIDIRECTIONAL_DOTTED
    ].includes(msg.type)) {
      sequenceIndex = sequenceIndex + sequenceIndexStep;
    }
    index++;
  }
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("createdActors", createdActors);
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("destroyedActors", destroyedActors);
  await drawActors(diagram2, actors, actorKeys, false);
  for (const e of messagesToDraw) {
    await drawMessage(diagram2, e.messageModel, e.lineStartY, diagObj);
  }
  if (conf.mirrorActors) {
    await drawActors(diagram2, actors, actorKeys, true);
  }
  backgrounds.forEach((e) => svgDraw_default.drawBackgroundRect(diagram2, e));
  fixLifeLineHeights(diagram2, actors, actorKeys, conf);
  for (const box2 of bounds.models.boxes) {
    box2.height = bounds.getVerticalPos() - box2.y;
    bounds.insert(box2.x, box2.y, box2.x + box2.width, box2.height);
    box2.startx = box2.x;
    box2.starty = box2.y;
    box2.stopx = box2.startx + box2.width;
    box2.stopy = box2.starty + box2.height;
    box2.stroke = "rgb(0,0,0, 0.5)";
    svgDraw_default.drawBox(diagram2, box2, conf);
  }
  if (hasBoxes) {
    bounds.bumpVerticalPos(conf.boxMargin);
  }
  const requiredBoxSize = drawActorsPopup(diagram2, actors, actorKeys, doc);
  const { bounds: box } = bounds.getBounds();
  if (box.startx === void 0) {
    box.startx = 0;
  }
  if (box.starty === void 0) {
    box.starty = 0;
  }
  if (box.stopx === void 0) {
    box.stopx = 0;
  }
  if (box.stopy === void 0) {
    box.stopy = 0;
  }
  let boxHeight = box.stopy - box.starty;
  if (boxHeight < requiredBoxSize.maxHeight) {
    boxHeight = requiredBoxSize.maxHeight;
  }
  let height = boxHeight + 2 * conf.diagramMarginY;
  if (conf.mirrorActors) {
    height = height - conf.boxMargin + conf.bottomMarginAdj;
  }
  let boxWidth = box.stopx - box.startx;
  if (boxWidth < requiredBoxSize.maxWidth) {
    boxWidth = requiredBoxSize.maxWidth;
  }
  const width = boxWidth + 2 * conf.diagramMarginX;
  if (title) {
    diagram2.append("text").text(title).attr("x", (box.stopx - box.startx) / 2 - 2 * conf.diagramMarginX).attr("y", -25);
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .configureSvgSize */ .v2)(diagram2, height, width, conf.useMaxWidth);
  const extraVertForTitle = title ? 40 : 0;
  diagram2.attr(
    "viewBox",
    box.startx - conf.diagramMarginX + " -" + (conf.diagramMarginY + extraVertForTitle) + " " + width + " " + (height + extraVertForTitle)
  );
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(`models:`, bounds.models);
}, "draw");
async function getMaxMessageWidthPerActor(actors, messages, diagObj) {
  const maxMessageWidthPerActor = {};
  for (const msg of messages) {
    if (actors.get(msg.to) && actors.get(msg.from)) {
      const actor = actors.get(msg.to);
      if (msg.placement === diagObj.db.PLACEMENT.LEFTOF && !actor.prevActor) {
        continue;
      }
      if (msg.placement === diagObj.db.PLACEMENT.RIGHTOF && !actor.nextActor) {
        continue;
      }
      const isNote = msg.placement !== void 0;
      const isMessage = !isNote;
      const textFont = isNote ? noteFont(conf) : messageFont(conf);
      const wrappedMessage = msg.wrap ? _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(msg.message, conf.width - 2 * conf.wrapPadding, textFont) : msg.message;
      const messageDimensions = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(wrappedMessage) ? await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateMathMLDimensions */ .nH)(msg.message, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()) : _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(wrappedMessage, textFont);
      const messageWidth = messageDimensions.width + 2 * conf.wrapPadding;
      if (isMessage && msg.from === actor.nextActor) {
        maxMessageWidthPerActor[msg.to] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[msg.to] || 0,
          messageWidth
        );
      } else if (isMessage && msg.from === actor.prevActor) {
        maxMessageWidthPerActor[msg.from] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[msg.from] || 0,
          messageWidth
        );
      } else if (isMessage && msg.from === msg.to) {
        maxMessageWidthPerActor[msg.from] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[msg.from] || 0,
          messageWidth / 2
        );
        maxMessageWidthPerActor[msg.to] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[msg.to] || 0,
          messageWidth / 2
        );
      } else if (msg.placement === diagObj.db.PLACEMENT.RIGHTOF) {
        maxMessageWidthPerActor[msg.from] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[msg.from] || 0,
          messageWidth
        );
      } else if (msg.placement === diagObj.db.PLACEMENT.LEFTOF) {
        maxMessageWidthPerActor[actor.prevActor] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
          maxMessageWidthPerActor[actor.prevActor] || 0,
          messageWidth
        );
      } else if (msg.placement === diagObj.db.PLACEMENT.OVER) {
        if (actor.prevActor) {
          maxMessageWidthPerActor[actor.prevActor] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
            maxMessageWidthPerActor[actor.prevActor] || 0,
            messageWidth / 2
          );
        }
        if (actor.nextActor) {
          maxMessageWidthPerActor[msg.from] = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
            maxMessageWidthPerActor[msg.from] || 0,
            messageWidth / 2
          );
        }
      }
    }
  }
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("maxMessageWidthPerActor:", maxMessageWidthPerActor);
  return maxMessageWidthPerActor;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getMaxMessageWidthPerActor, "getMaxMessageWidthPerActor");
var getRequiredPopupWidth = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(actor) {
  let requiredPopupWidth = 0;
  const textFont = actorFont(conf);
  for (const key in actor.links) {
    const labelDimensions = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(key, textFont);
    const labelWidth = labelDimensions.width + 2 * conf.wrapPadding + 2 * conf.boxMargin;
    if (requiredPopupWidth < labelWidth) {
      requiredPopupWidth = labelWidth;
    }
  }
  return requiredPopupWidth;
}, "getRequiredPopupWidth");
async function calculateActorMargins(actors, actorToMessageWidth, boxes) {
  let maxHeight = 0;
  for (const prop of actors.keys()) {
    const actor = actors.get(prop);
    if (actor.wrap) {
      actor.description = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(
        actor.description,
        conf.width - 2 * conf.wrapPadding,
        actorFont(conf)
      );
    }
    const actDims = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(actor.description) ? await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateMathMLDimensions */ .nH)(actor.description, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()) : _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(actor.description, actorFont(conf));
    actor.width = actor.wrap ? conf.width : _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, actDims.width + 2 * conf.wrapPadding);
    actor.height = actor.wrap ? _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(actDims.height, conf.height) : conf.height;
    maxHeight = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(maxHeight, actor.height);
  }
  for (const actorKey in actorToMessageWidth) {
    const actor = actors.get(actorKey);
    if (!actor) {
      continue;
    }
    const nextActor = actors.get(actor.nextActor);
    if (!nextActor) {
      const messageWidth2 = actorToMessageWidth[actorKey];
      const actorWidth2 = messageWidth2 + conf.actorMargin - actor.width / 2;
      actor.margin = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(actorWidth2, conf.actorMargin);
      continue;
    }
    const messageWidth = actorToMessageWidth[actorKey];
    const actorWidth = messageWidth + conf.actorMargin - actor.width / 2 - nextActor.width / 2;
    actor.margin = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(actorWidth, conf.actorMargin);
  }
  let maxBoxHeight = 0;
  boxes.forEach((box) => {
    const textFont = messageFont(conf);
    let totalWidth = box.actorKeys.reduce((total, aKey) => {
      return total += actors.get(aKey).width + (actors.get(aKey).margin || 0);
    }, 0);
    totalWidth -= 2 * conf.boxTextMargin;
    if (box.wrap) {
      box.name = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(box.name, totalWidth - 2 * conf.wrapPadding, textFont);
    }
    const boxMsgDimensions = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(box.name, textFont);
    maxBoxHeight = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(boxMsgDimensions.height, maxBoxHeight);
    const minWidth = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(totalWidth, boxMsgDimensions.width + 2 * conf.wrapPadding);
    box.margin = conf.boxTextMargin;
    if (totalWidth < minWidth) {
      const missing = (minWidth - totalWidth) / 2;
      box.margin += missing;
    }
  });
  boxes.forEach((box) => box.textMaxHeight = maxBoxHeight);
  return _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(maxHeight, conf.height);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(calculateActorMargins, "calculateActorMargins");
var buildNoteModel = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(msg, actors, diagObj) {
  const fromActor = actors.get(msg.from);
  const toActor = actors.get(msg.to);
  const startx = fromActor.x;
  const stopx = toActor.x;
  const shouldWrap = msg.wrap && msg.message;
  let textDimensions = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .hasKatex */ .l0)(msg.message) ? await (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateMathMLDimensions */ .nH)(msg.message, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()) : _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(
    shouldWrap ? _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(msg.message, conf.width, noteFont(conf)) : msg.message,
    noteFont(conf)
  );
  const noteModel = {
    width: shouldWrap ? conf.width : _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, textDimensions.width + 2 * conf.noteMargin),
    height: 0,
    startx: fromActor.x,
    stopx: 0,
    starty: 0,
    stopy: 0,
    message: msg.message
  };
  if (msg.placement === diagObj.db.PLACEMENT.RIGHTOF) {
    noteModel.width = shouldWrap ? _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, textDimensions.width) : _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
      fromActor.width / 2 + toActor.width / 2,
      textDimensions.width + 2 * conf.noteMargin
    );
    noteModel.startx = startx + (fromActor.width + conf.actorMargin) / 2;
  } else if (msg.placement === diagObj.db.PLACEMENT.LEFTOF) {
    noteModel.width = shouldWrap ? _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, textDimensions.width + 2 * conf.noteMargin) : _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
      fromActor.width / 2 + toActor.width / 2,
      textDimensions.width + 2 * conf.noteMargin
    );
    noteModel.startx = startx - noteModel.width + (fromActor.width - conf.actorMargin) / 2;
  } else if (msg.to === msg.from) {
    textDimensions = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(
      shouldWrap ? _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(msg.message, _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, fromActor.width), noteFont(conf)) : msg.message,
      noteFont(conf)
    );
    noteModel.width = shouldWrap ? _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(conf.width, fromActor.width) : _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(fromActor.width, conf.width, textDimensions.width + 2 * conf.noteMargin);
    noteModel.startx = startx + (fromActor.width - noteModel.width) / 2;
  } else {
    noteModel.width = Math.abs(startx + fromActor.width / 2 - (stopx + toActor.width / 2)) + conf.actorMargin;
    noteModel.startx = startx < stopx ? startx + fromActor.width / 2 - conf.actorMargin / 2 : stopx + toActor.width / 2 - conf.actorMargin / 2;
  }
  if (shouldWrap) {
    noteModel.message = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(
      msg.message,
      noteModel.width - 2 * conf.wrapPadding,
      noteFont(conf)
    );
  }
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(
    `NM:[${noteModel.startx},${noteModel.stopx},${noteModel.starty},${noteModel.stopy}:${noteModel.width},${noteModel.height}=${msg.message}]`
  );
  return noteModel;
}, "buildNoteModel");
var buildMessageModel = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(msg, actors, diagObj) {
  if (![
    diagObj.db.LINETYPE.SOLID_OPEN,
    diagObj.db.LINETYPE.DOTTED_OPEN,
    diagObj.db.LINETYPE.SOLID,
    diagObj.db.LINETYPE.DOTTED,
    diagObj.db.LINETYPE.SOLID_CROSS,
    diagObj.db.LINETYPE.DOTTED_CROSS,
    diagObj.db.LINETYPE.SOLID_POINT,
    diagObj.db.LINETYPE.DOTTED_POINT,
    diagObj.db.LINETYPE.BIDIRECTIONAL_SOLID,
    diagObj.db.LINETYPE.BIDIRECTIONAL_DOTTED
  ].includes(msg.type)) {
    return {};
  }
  const [fromLeft, fromRight] = activationBounds(msg.from, actors);
  const [toLeft, toRight] = activationBounds(msg.to, actors);
  const isArrowToRight = fromLeft <= toLeft;
  let startx = isArrowToRight ? fromRight : fromLeft;
  let stopx = isArrowToRight ? toLeft : toRight;
  const isArrowToActivation = Math.abs(toLeft - toRight) > 2;
  const adjustValue = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((value) => {
    return isArrowToRight ? -value : value;
  }, "adjustValue");
  if (msg.from === msg.to) {
    stopx = startx;
  } else {
    if (msg.activate && !isArrowToActivation) {
      stopx += adjustValue(conf.activationWidth / 2 - 1);
    }
    if (![diagObj.db.LINETYPE.SOLID_OPEN, diagObj.db.LINETYPE.DOTTED_OPEN].includes(msg.type)) {
      stopx += adjustValue(3);
    }
    if ([diagObj.db.LINETYPE.BIDIRECTIONAL_SOLID, diagObj.db.LINETYPE.BIDIRECTIONAL_DOTTED].includes(
      msg.type
    )) {
      startx -= adjustValue(3);
    }
  }
  const allBounds = [fromLeft, fromRight, toLeft, toRight];
  const boundedWidth = Math.abs(startx - stopx);
  if (msg.wrap && msg.message) {
    msg.message = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.wrapLabel(
      msg.message,
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(boundedWidth + 2 * conf.wrapPadding, conf.width),
      messageFont(conf)
    );
  }
  const msgDims = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.calculateTextDimensions(msg.message, messageFont(conf));
  return {
    width: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
      msg.wrap ? 0 : msgDims.width + 2 * conf.wrapPadding,
      boundedWidth + 2 * conf.wrapPadding,
      conf.width
    ),
    height: 0,
    startx,
    stopx,
    starty: 0,
    stopy: 0,
    message: msg.message,
    type: msg.type,
    wrap: msg.wrap,
    fromBounds: Math.min.apply(null, allBounds),
    toBounds: Math.max.apply(null, allBounds)
  };
}, "buildMessageModel");
var calculateLoopBounds = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(messages, actors, _maxWidthPerActor, diagObj) {
  const loops = {};
  const stack = [];
  let current, noteModel, msgModel;
  for (const msg of messages) {
    switch (msg.type) {
      case diagObj.db.LINETYPE.LOOP_START:
      case diagObj.db.LINETYPE.ALT_START:
      case diagObj.db.LINETYPE.OPT_START:
      case diagObj.db.LINETYPE.PAR_START:
      case diagObj.db.LINETYPE.PAR_OVER_START:
      case diagObj.db.LINETYPE.CRITICAL_START:
      case diagObj.db.LINETYPE.BREAK_START:
        stack.push({
          id: msg.id,
          msg: msg.message,
          from: Number.MAX_SAFE_INTEGER,
          to: Number.MIN_SAFE_INTEGER,
          width: 0
        });
        break;
      case diagObj.db.LINETYPE.ALT_ELSE:
      case diagObj.db.LINETYPE.PAR_AND:
      case diagObj.db.LINETYPE.CRITICAL_OPTION:
        if (msg.message) {
          current = stack.pop();
          loops[current.id] = current;
          loops[msg.id] = current;
          stack.push(current);
        }
        break;
      case diagObj.db.LINETYPE.LOOP_END:
      case diagObj.db.LINETYPE.ALT_END:
      case diagObj.db.LINETYPE.OPT_END:
      case diagObj.db.LINETYPE.PAR_END:
      case diagObj.db.LINETYPE.CRITICAL_END:
      case diagObj.db.LINETYPE.BREAK_END:
        current = stack.pop();
        loops[current.id] = current;
        break;
      case diagObj.db.LINETYPE.ACTIVE_START:
        {
          const actorRect = actors.get(msg.from ? msg.from : msg.to.actor);
          const stackedSize = actorActivations(msg.from ? msg.from : msg.to.actor).length;
          const x = actorRect.x + actorRect.width / 2 + (stackedSize - 1) * conf.activationWidth / 2;
          const toAdd = {
            startx: x,
            stopx: x + conf.activationWidth,
            actor: msg.from,
            enabled: true
          };
          bounds.activations.push(toAdd);
        }
        break;
      case diagObj.db.LINETYPE.ACTIVE_END:
        {
          const lastActorActivationIdx = bounds.activations.map((a) => a.actor).lastIndexOf(msg.from);
          bounds.activations.splice(lastActorActivationIdx, 1).splice(0, 1);
        }
        break;
    }
    const isNote = msg.placement !== void 0;
    if (isNote) {
      noteModel = await buildNoteModel(msg, actors, diagObj);
      msg.noteModel = noteModel;
      stack.forEach((stk) => {
        current = stk;
        current.from = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(current.from, noteModel.startx);
        current.to = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(current.to, noteModel.startx + noteModel.width);
        current.width = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(current.width, Math.abs(current.from - current.to)) - conf.labelBoxWidth;
      });
    } else {
      msgModel = buildMessageModel(msg, actors, diagObj);
      msg.msgModel = msgModel;
      if (msgModel.startx && msgModel.stopx && stack.length > 0) {
        stack.forEach((stk) => {
          current = stk;
          if (msgModel.startx === msgModel.stopx) {
            const from = actors.get(msg.from);
            const to = actors.get(msg.to);
            current.from = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(
              from.x - msgModel.width / 2,
              from.x - from.width / 2,
              current.from
            );
            current.to = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(
              to.x + msgModel.width / 2,
              to.x + from.width / 2,
              current.to
            );
            current.width = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(current.width, Math.abs(current.to - current.from)) - conf.labelBoxWidth;
          } else {
            current.from = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMin(msgModel.startx, current.from);
            current.to = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(msgModel.stopx, current.to);
            current.width = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.getMax(current.width, msgModel.width) - conf.labelBoxWidth;
          }
        });
      }
    }
  }
  bounds.activations = [];
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("Loop type widths:", loops);
  return loops;
}, "calculateLoopBounds");
var sequenceRenderer_default = {
  bounds,
  drawActors,
  drawActorsPopup,
  setConf,
  draw
};

// src/diagrams/sequence/sequenceDiagram.ts
var diagram = {
  parser: sequenceDiagram_default,
  get db() {
    return new SequenceDB();
  },
  renderer: sequenceRenderer_default,
  styles: styles_default,
  init: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((cnf) => {
    if (!cnf.sequence) {
      cnf.sequence = {};
    }
    if (cnf.wrap) {
      cnf.sequence.wrap = cnf.wrap;
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setConfig2 */ .Y4)({ sequence: { wrap: cnf.wrap } });
    }
  }, "init")
};



/***/ })

}]);
//# sourceMappingURL=2453.ebdb135eb902bf82e103.js.map?v=ebdb135eb902bf82e103