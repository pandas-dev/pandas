(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3197,7061],{

/***/ 77647:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  e: () => (/* reexport */ BaseRegExpVisitor),
  O: () => (/* reexport */ RegExpParser)
});

;// CONCATENATED MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/utils.js
function cc(char) {
    return char.charCodeAt(0);
}
function insertToSet(item, set) {
    if (Array.isArray(item)) {
        item.forEach(function (subItem) {
            set.push(subItem);
        });
    }
    else {
        set.push(item);
    }
}
function addFlag(flagObj, flagKey) {
    if (flagObj[flagKey] === true) {
        throw "duplicate flag " + flagKey;
    }
    const x = flagObj[flagKey];
    flagObj[flagKey] = true;
}
function ASSERT_EXISTS(obj) {
    // istanbul ignore next
    if (obj === undefined) {
        throw Error("Internal Error - Should never get here!");
    }
    return true;
}
// istanbul ignore next
function ASSERT_NEVER_REACH_HERE() {
    throw Error("Internal Error - Should never get here!");
}
function isCharacter(obj) {
    return obj["type"] === "Character";
}
//# sourceMappingURL=utils.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/character-classes.js

const digitsCharCodes = [];
for (let i = cc("0"); i <= cc("9"); i++) {
    digitsCharCodes.push(i);
}
const wordCharCodes = [cc("_")].concat(digitsCharCodes);
for (let i = cc("a"); i <= cc("z"); i++) {
    wordCharCodes.push(i);
}
for (let i = cc("A"); i <= cc("Z"); i++) {
    wordCharCodes.push(i);
}
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/RegExp#character-classes
const whitespaceCodes = [
    cc(" "),
    cc("\f"),
    cc("\n"),
    cc("\r"),
    cc("\t"),
    cc("\v"),
    cc("\t"),
    cc("\u00a0"),
    cc("\u1680"),
    cc("\u2000"),
    cc("\u2001"),
    cc("\u2002"),
    cc("\u2003"),
    cc("\u2004"),
    cc("\u2005"),
    cc("\u2006"),
    cc("\u2007"),
    cc("\u2008"),
    cc("\u2009"),
    cc("\u200a"),
    cc("\u2028"),
    cc("\u2029"),
    cc("\u202f"),
    cc("\u205f"),
    cc("\u3000"),
    cc("\ufeff"),
];
//# sourceMappingURL=character-classes.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/regexp-parser.js


// consts and utilities
const hexDigitPattern = /[0-9a-fA-F]/;
const decimalPattern = /[0-9]/;
const decimalPatternNoZero = /[1-9]/;
// https://hackernoon.com/the-madness-of-parsing-real-world-javascript-regexps-d9ee336df983
// https://www.ecma-international.org/ecma-262/8.0/index.html#prod-Pattern
class RegExpParser {
    constructor() {
        this.idx = 0;
        this.input = "";
        this.groupIdx = 0;
    }
    saveState() {
        return {
            idx: this.idx,
            input: this.input,
            groupIdx: this.groupIdx,
        };
    }
    restoreState(newState) {
        this.idx = newState.idx;
        this.input = newState.input;
        this.groupIdx = newState.groupIdx;
    }
    pattern(input) {
        // parser state
        this.idx = 0;
        this.input = input;
        this.groupIdx = 0;
        this.consumeChar("/");
        const value = this.disjunction();
        this.consumeChar("/");
        const flags = {
            type: "Flags",
            loc: { begin: this.idx, end: input.length },
            global: false,
            ignoreCase: false,
            multiLine: false,
            unicode: false,
            sticky: false,
        };
        while (this.isRegExpFlag()) {
            switch (this.popChar()) {
                case "g":
                    addFlag(flags, "global");
                    break;
                case "i":
                    addFlag(flags, "ignoreCase");
                    break;
                case "m":
                    addFlag(flags, "multiLine");
                    break;
                case "u":
                    addFlag(flags, "unicode");
                    break;
                case "y":
                    addFlag(flags, "sticky");
                    break;
            }
        }
        if (this.idx !== this.input.length) {
            throw Error("Redundant input: " + this.input.substring(this.idx));
        }
        return {
            type: "Pattern",
            flags: flags,
            value: value,
            loc: this.loc(0),
        };
    }
    disjunction() {
        const alts = [];
        const begin = this.idx;
        alts.push(this.alternative());
        while (this.peekChar() === "|") {
            this.consumeChar("|");
            alts.push(this.alternative());
        }
        return { type: "Disjunction", value: alts, loc: this.loc(begin) };
    }
    alternative() {
        const terms = [];
        const begin = this.idx;
        while (this.isTerm()) {
            terms.push(this.term());
        }
        return { type: "Alternative", value: terms, loc: this.loc(begin) };
    }
    term() {
        if (this.isAssertion()) {
            return this.assertion();
        }
        else {
            return this.atom();
        }
    }
    assertion() {
        const begin = this.idx;
        switch (this.popChar()) {
            case "^":
                return {
                    type: "StartAnchor",
                    loc: this.loc(begin),
                };
            case "$":
                return { type: "EndAnchor", loc: this.loc(begin) };
            // '\b' or '\B'
            case "\\":
                switch (this.popChar()) {
                    case "b":
                        return {
                            type: "WordBoundary",
                            loc: this.loc(begin),
                        };
                    case "B":
                        return {
                            type: "NonWordBoundary",
                            loc: this.loc(begin),
                        };
                }
                // istanbul ignore next
                throw Error("Invalid Assertion Escape");
            // '(?=' or '(?!'
            case "(":
                this.consumeChar("?");
                let type;
                switch (this.popChar()) {
                    case "=":
                        type = "Lookahead";
                        break;
                    case "!":
                        type = "NegativeLookahead";
                        break;
                }
                ASSERT_EXISTS(type);
                const disjunction = this.disjunction();
                this.consumeChar(")");
                return {
                    type: type,
                    value: disjunction,
                    loc: this.loc(begin),
                };
        }
        // istanbul ignore next
        return ASSERT_NEVER_REACH_HERE();
    }
    quantifier(isBacktracking = false) {
        let range = undefined;
        const begin = this.idx;
        switch (this.popChar()) {
            case "*":
                range = {
                    atLeast: 0,
                    atMost: Infinity,
                };
                break;
            case "+":
                range = {
                    atLeast: 1,
                    atMost: Infinity,
                };
                break;
            case "?":
                range = {
                    atLeast: 0,
                    atMost: 1,
                };
                break;
            case "{":
                const atLeast = this.integerIncludingZero();
                switch (this.popChar()) {
                    case "}":
                        range = {
                            atLeast: atLeast,
                            atMost: atLeast,
                        };
                        break;
                    case ",":
                        let atMost;
                        if (this.isDigit()) {
                            atMost = this.integerIncludingZero();
                            range = {
                                atLeast: atLeast,
                                atMost: atMost,
                            };
                        }
                        else {
                            range = {
                                atLeast: atLeast,
                                atMost: Infinity,
                            };
                        }
                        this.consumeChar("}");
                        break;
                }
                // throwing exceptions from "ASSERT_EXISTS" during backtracking
                // causes severe performance degradations
                if (isBacktracking === true && range === undefined) {
                    return undefined;
                }
                ASSERT_EXISTS(range);
                break;
        }
        // throwing exceptions from "ASSERT_EXISTS" during backtracking
        // causes severe performance degradations
        if (isBacktracking === true && range === undefined) {
            return undefined;
        }
        // istanbul ignore else
        if (ASSERT_EXISTS(range)) {
            if (this.peekChar(0) === "?") {
                this.consumeChar("?");
                range.greedy = false;
            }
            else {
                range.greedy = true;
            }
            range.type = "Quantifier";
            range.loc = this.loc(begin);
            return range;
        }
    }
    atom() {
        let atom;
        const begin = this.idx;
        switch (this.peekChar()) {
            case ".":
                atom = this.dotAll();
                break;
            case "\\":
                atom = this.atomEscape();
                break;
            case "[":
                atom = this.characterClass();
                break;
            case "(":
                atom = this.group();
                break;
        }
        if (atom === undefined && this.isPatternCharacter()) {
            atom = this.patternCharacter();
        }
        // istanbul ignore else
        if (ASSERT_EXISTS(atom)) {
            atom.loc = this.loc(begin);
            if (this.isQuantifier()) {
                atom.quantifier = this.quantifier();
            }
            return atom;
        }
        // istanbul ignore next
        return ASSERT_NEVER_REACH_HERE();
    }
    dotAll() {
        this.consumeChar(".");
        return {
            type: "Set",
            complement: true,
            value: [cc("\n"), cc("\r"), cc("\u2028"), cc("\u2029")],
        };
    }
    atomEscape() {
        this.consumeChar("\\");
        switch (this.peekChar()) {
            case "1":
            case "2":
            case "3":
            case "4":
            case "5":
            case "6":
            case "7":
            case "8":
            case "9":
                return this.decimalEscapeAtom();
            case "d":
            case "D":
            case "s":
            case "S":
            case "w":
            case "W":
                return this.characterClassEscape();
            case "f":
            case "n":
            case "r":
            case "t":
            case "v":
                return this.controlEscapeAtom();
            case "c":
                return this.controlLetterEscapeAtom();
            case "0":
                return this.nulCharacterAtom();
            case "x":
                return this.hexEscapeSequenceAtom();
            case "u":
                return this.regExpUnicodeEscapeSequenceAtom();
            default:
                return this.identityEscapeAtom();
        }
    }
    decimalEscapeAtom() {
        const value = this.positiveInteger();
        return { type: "GroupBackReference", value: value };
    }
    characterClassEscape() {
        let set;
        let complement = false;
        switch (this.popChar()) {
            case "d":
                set = digitsCharCodes;
                break;
            case "D":
                set = digitsCharCodes;
                complement = true;
                break;
            case "s":
                set = whitespaceCodes;
                break;
            case "S":
                set = whitespaceCodes;
                complement = true;
                break;
            case "w":
                set = wordCharCodes;
                break;
            case "W":
                set = wordCharCodes;
                complement = true;
                break;
        }
        // istanbul ignore else
        if (ASSERT_EXISTS(set)) {
            return { type: "Set", value: set, complement: complement };
        }
        // istanbul ignore next
        return ASSERT_NEVER_REACH_HERE();
    }
    controlEscapeAtom() {
        let escapeCode;
        switch (this.popChar()) {
            case "f":
                escapeCode = cc("\f");
                break;
            case "n":
                escapeCode = cc("\n");
                break;
            case "r":
                escapeCode = cc("\r");
                break;
            case "t":
                escapeCode = cc("\t");
                break;
            case "v":
                escapeCode = cc("\v");
                break;
        }
        // istanbul ignore else
        if (ASSERT_EXISTS(escapeCode)) {
            return { type: "Character", value: escapeCode };
        }
        // istanbul ignore next
        return ASSERT_NEVER_REACH_HERE();
    }
    controlLetterEscapeAtom() {
        this.consumeChar("c");
        const letter = this.popChar();
        if (/[a-zA-Z]/.test(letter) === false) {
            throw Error("Invalid ");
        }
        const letterCode = letter.toUpperCase().charCodeAt(0) - 64;
        return { type: "Character", value: letterCode };
    }
    nulCharacterAtom() {
        // TODO implement '[lookahead ∉ DecimalDigit]'
        // TODO: for the deprecated octal escape sequence
        this.consumeChar("0");
        return { type: "Character", value: cc("\0") };
    }
    hexEscapeSequenceAtom() {
        this.consumeChar("x");
        return this.parseHexDigits(2);
    }
    regExpUnicodeEscapeSequenceAtom() {
        this.consumeChar("u");
        return this.parseHexDigits(4);
    }
    identityEscapeAtom() {
        // TODO: implement "SourceCharacter but not UnicodeIDContinue"
        // // http://unicode.org/reports/tr31/#Specific_Character_Adjustments
        const escapedChar = this.popChar();
        return { type: "Character", value: cc(escapedChar) };
    }
    classPatternCharacterAtom() {
        switch (this.peekChar()) {
            // istanbul ignore next
            case "\n":
            // istanbul ignore next
            case "\r":
            // istanbul ignore next
            case "\u2028":
            // istanbul ignore next
            case "\u2029":
            // istanbul ignore next
            case "\\":
            // istanbul ignore next
            case "]":
                throw Error("TBD");
            default:
                const nextChar = this.popChar();
                return { type: "Character", value: cc(nextChar) };
        }
    }
    characterClass() {
        const set = [];
        let complement = false;
        this.consumeChar("[");
        if (this.peekChar(0) === "^") {
            this.consumeChar("^");
            complement = true;
        }
        while (this.isClassAtom()) {
            const from = this.classAtom();
            const isFromSingleChar = from.type === "Character";
            if (isCharacter(from) && this.isRangeDash()) {
                this.consumeChar("-");
                const to = this.classAtom();
                const isToSingleChar = to.type === "Character";
                // a range can only be used when both sides are single characters
                if (isCharacter(to)) {
                    if (to.value < from.value) {
                        throw Error("Range out of order in character class");
                    }
                    set.push({ from: from.value, to: to.value });
                }
                else {
                    // literal dash
                    insertToSet(from.value, set);
                    set.push(cc("-"));
                    insertToSet(to.value, set);
                }
            }
            else {
                insertToSet(from.value, set);
            }
        }
        this.consumeChar("]");
        return { type: "Set", complement: complement, value: set };
    }
    classAtom() {
        switch (this.peekChar()) {
            // istanbul ignore next
            case "]":
            // istanbul ignore next
            case "\n":
            // istanbul ignore next
            case "\r":
            // istanbul ignore next
            case "\u2028":
            // istanbul ignore next
            case "\u2029":
                throw Error("TBD");
            case "\\":
                return this.classEscape();
            default:
                return this.classPatternCharacterAtom();
        }
    }
    classEscape() {
        this.consumeChar("\\");
        switch (this.peekChar()) {
            // Matches a backspace.
            // (Not to be confused with \b word boundary outside characterClass)
            case "b":
                this.consumeChar("b");
                return { type: "Character", value: cc("\u0008") };
            case "d":
            case "D":
            case "s":
            case "S":
            case "w":
            case "W":
                return this.characterClassEscape();
            case "f":
            case "n":
            case "r":
            case "t":
            case "v":
                return this.controlEscapeAtom();
            case "c":
                return this.controlLetterEscapeAtom();
            case "0":
                return this.nulCharacterAtom();
            case "x":
                return this.hexEscapeSequenceAtom();
            case "u":
                return this.regExpUnicodeEscapeSequenceAtom();
            default:
                return this.identityEscapeAtom();
        }
    }
    group() {
        let capturing = true;
        this.consumeChar("(");
        switch (this.peekChar(0)) {
            case "?":
                this.consumeChar("?");
                this.consumeChar(":");
                capturing = false;
                break;
            default:
                this.groupIdx++;
                break;
        }
        const value = this.disjunction();
        this.consumeChar(")");
        const groupAst = {
            type: "Group",
            capturing: capturing,
            value: value,
        };
        if (capturing) {
            groupAst["idx"] = this.groupIdx;
        }
        return groupAst;
    }
    positiveInteger() {
        let number = this.popChar();
        // istanbul ignore next - can't ever get here due to previous lookahead checks
        // still implementing this error checking in case this ever changes.
        if (decimalPatternNoZero.test(number) === false) {
            throw Error("Expecting a positive integer");
        }
        while (decimalPattern.test(this.peekChar(0))) {
            number += this.popChar();
        }
        return parseInt(number, 10);
    }
    integerIncludingZero() {
        let number = this.popChar();
        if (decimalPattern.test(number) === false) {
            throw Error("Expecting an integer");
        }
        while (decimalPattern.test(this.peekChar(0))) {
            number += this.popChar();
        }
        return parseInt(number, 10);
    }
    patternCharacter() {
        const nextChar = this.popChar();
        switch (nextChar) {
            // istanbul ignore next
            case "\n":
            // istanbul ignore next
            case "\r":
            // istanbul ignore next
            case "\u2028":
            // istanbul ignore next
            case "\u2029":
            // istanbul ignore next
            case "^":
            // istanbul ignore next
            case "$":
            // istanbul ignore next
            case "\\":
            // istanbul ignore next
            case ".":
            // istanbul ignore next
            case "*":
            // istanbul ignore next
            case "+":
            // istanbul ignore next
            case "?":
            // istanbul ignore next
            case "(":
            // istanbul ignore next
            case ")":
            // istanbul ignore next
            case "[":
            // istanbul ignore next
            case "|":
                // istanbul ignore next
                throw Error("TBD");
            default:
                return { type: "Character", value: cc(nextChar) };
        }
    }
    isRegExpFlag() {
        switch (this.peekChar(0)) {
            case "g":
            case "i":
            case "m":
            case "u":
            case "y":
                return true;
            default:
                return false;
        }
    }
    isRangeDash() {
        return this.peekChar() === "-" && this.isClassAtom(1);
    }
    isDigit() {
        return decimalPattern.test(this.peekChar(0));
    }
    isClassAtom(howMuch = 0) {
        switch (this.peekChar(howMuch)) {
            case "]":
            case "\n":
            case "\r":
            case "\u2028":
            case "\u2029":
                return false;
            default:
                return true;
        }
    }
    isTerm() {
        return this.isAtom() || this.isAssertion();
    }
    isAtom() {
        if (this.isPatternCharacter()) {
            return true;
        }
        switch (this.peekChar(0)) {
            case ".":
            case "\\": // atomEscape
            case "[": // characterClass
            // TODO: isAtom must be called before isAssertion - disambiguate
            case "(": // group
                return true;
            default:
                return false;
        }
    }
    isAssertion() {
        switch (this.peekChar(0)) {
            case "^":
            case "$":
                return true;
            // '\b' or '\B'
            case "\\":
                switch (this.peekChar(1)) {
                    case "b":
                    case "B":
                        return true;
                    default:
                        return false;
                }
            // '(?=' or '(?!'
            case "(":
                return (this.peekChar(1) === "?" &&
                    (this.peekChar(2) === "=" || this.peekChar(2) === "!"));
            default:
                return false;
        }
    }
    isQuantifier() {
        const prevState = this.saveState();
        try {
            return this.quantifier(true) !== undefined;
        }
        catch (e) {
            return false;
        }
        finally {
            this.restoreState(prevState);
        }
    }
    isPatternCharacter() {
        switch (this.peekChar()) {
            case "^":
            case "$":
            case "\\":
            case ".":
            case "*":
            case "+":
            case "?":
            case "(":
            case ")":
            case "[":
            case "|":
            case "/":
            case "\n":
            case "\r":
            case "\u2028":
            case "\u2029":
                return false;
            default:
                return true;
        }
    }
    parseHexDigits(howMany) {
        let hexString = "";
        for (let i = 0; i < howMany; i++) {
            const hexChar = this.popChar();
            if (hexDigitPattern.test(hexChar) === false) {
                throw Error("Expecting a HexDecimal digits");
            }
            hexString += hexChar;
        }
        const charCode = parseInt(hexString, 16);
        return { type: "Character", value: charCode };
    }
    peekChar(howMuch = 0) {
        return this.input[this.idx + howMuch];
    }
    popChar() {
        const nextChar = this.peekChar(0);
        this.consumeChar(undefined);
        return nextChar;
    }
    consumeChar(char) {
        if (char !== undefined && this.input[this.idx] !== char) {
            throw Error("Expected: '" +
                char +
                "' but found: '" +
                this.input[this.idx] +
                "' at offset: " +
                this.idx);
        }
        if (this.idx >= this.input.length) {
            throw Error("Unexpected end of input");
        }
        this.idx++;
    }
    loc(begin) {
        return { begin: begin, end: this.idx };
    }
}
//# sourceMappingURL=regexp-parser.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/base-regexp-visitor.js
class BaseRegExpVisitor {
    visitChildren(node) {
        for (const key in node) {
            const child = node[key];
            /* istanbul ignore else */
            if (node.hasOwnProperty(key)) {
                if (child.type !== undefined) {
                    this.visit(child);
                }
                else if (Array.isArray(child)) {
                    child.forEach((subChild) => {
                        this.visit(subChild);
                    }, this);
                }
            }
        }
    }
    visit(node) {
        switch (node.type) {
            case "Pattern":
                this.visitPattern(node);
                break;
            case "Flags":
                this.visitFlags(node);
                break;
            case "Disjunction":
                this.visitDisjunction(node);
                break;
            case "Alternative":
                this.visitAlternative(node);
                break;
            case "StartAnchor":
                this.visitStartAnchor(node);
                break;
            case "EndAnchor":
                this.visitEndAnchor(node);
                break;
            case "WordBoundary":
                this.visitWordBoundary(node);
                break;
            case "NonWordBoundary":
                this.visitNonWordBoundary(node);
                break;
            case "Lookahead":
                this.visitLookahead(node);
                break;
            case "NegativeLookahead":
                this.visitNegativeLookahead(node);
                break;
            case "Character":
                this.visitCharacter(node);
                break;
            case "Set":
                this.visitSet(node);
                break;
            case "Group":
                this.visitGroup(node);
                break;
            case "GroupBackReference":
                this.visitGroupBackReference(node);
                break;
            case "Quantifier":
                this.visitQuantifier(node);
                break;
        }
        this.visitChildren(node);
    }
    visitPattern(node) { }
    visitFlags(node) { }
    visitDisjunction(node) { }
    visitAlternative(node) { }
    // Assertion
    visitStartAnchor(node) { }
    visitEndAnchor(node) { }
    visitWordBoundary(node) { }
    visitNonWordBoundary(node) { }
    visitLookahead(node) { }
    visitNegativeLookahead(node) { }
    // atoms
    visitCharacter(node) { }
    visitSet(node) { }
    visitGroup(node) { }
    visitGroupBackReference(node) { }
    visitQuantifier(node) { }
}
//# sourceMappingURL=base-regexp-visitor.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/api.js


//# sourceMappingURL=api.js.map

/***/ }),

/***/ 20078:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   z: () => (/* binding */ createGitGraphServices)
/* harmony export */ });
/* unused harmony export GitGraphModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/gitGraph/module.ts


// src/language/gitGraph/tokenBuilder.ts
var GitGraphTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "GitGraphTokenBuilder");
  }
  constructor() {
    super(["gitGraph"]);
  }
};

// src/language/gitGraph/module.ts
var GitGraphModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new GitGraphTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .CommonValueConverter */ .nr(), "ValueConverter")
  }
};
function createGitGraphServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const GitGraph = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .GitGraphGeneratedModule */ .vn,
    GitGraphModule
  );
  shared.ServiceRegistry.register(GitGraph);
  return { shared, GitGraph };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createGitGraphServices, "createGitGraphServices");




/***/ }),

/***/ 51991:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   T: () => (/* binding */ createRadarServices)
/* harmony export */ });
/* unused harmony export RadarModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/radar/module.ts


// src/language/radar/tokenBuilder.ts
var RadarTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "RadarTokenBuilder");
  }
  constructor() {
    super(["radar-beta"]);
  }
};

// src/language/radar/module.ts
var RadarModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new RadarTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .CommonValueConverter */ .nr(), "ValueConverter")
  }
};
function createRadarServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const Radar = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .RadarGeneratedModule */ .gB,
    RadarModule
  );
  shared.ServiceRegistry.register(Radar);
  return { shared, Radar };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createRadarServices, "createRadarServices");




/***/ }),

/***/ 60960:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  T7: () => (/* binding */ AbstractMermaidTokenBuilder),
  kb: () => (/* binding */ AbstractMermaidValueConverter),
  Qr: () => (/* binding */ ArchitectureGeneratedModule),
  nr: () => (/* binding */ CommonValueConverter),
  vn: () => (/* binding */ GitGraphGeneratedModule),
  F_: () => (/* binding */ InfoGeneratedModule),
  GS: () => (/* binding */ MermaidGeneratedSharedModule),
  bb: () => (/* binding */ PacketGeneratedModule),
  WH: () => (/* binding */ PieGeneratedModule),
  gB: () => (/* binding */ RadarGeneratedModule),
  eW: () => (/* binding */ __name)
});

// UNUSED EXPORTS: Architecture, Branch, Commit, CommonTokenBuilder, GitGraph, Info, Merge, Packet, PacketBlock, Pie, PieSection, Radar, Statement, isArchitecture, isBranch, isCommit, isCommon, isGitGraph, isInfo, isMerge, isPacket, isPacketBlock, isPie, isPieSection

// EXTERNAL MODULE: ../node_modules/langium/lib/syntax-tree.js
var syntax_tree = __webpack_require__(91303);
// EXTERNAL MODULE: ../node_modules/langium/lib/default-module.js + 42 modules
var default_module = __webpack_require__(73001);
// EXTERNAL MODULE: ../node_modules/langium/lib/dependency-injection.js
var dependency_injection = __webpack_require__(81210);
// EXTERNAL MODULE: ../node_modules/langium/lib/languages/generated/ast.js
var ast = __webpack_require__(34905);
// EXTERNAL MODULE: ../node_modules/langium/lib/workspace/file-system-provider.js
var file_system_provider = __webpack_require__(44014);
// EXTERNAL MODULE: ../node_modules/vscode-uri/lib/esm/index.mjs
var esm = __webpack_require__(37943);
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/grammar-loader.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/





const minimalGrammarModule = {
    Grammar: () => undefined,
    LanguageMetaData: () => ({
        caseInsensitive: false,
        fileExtensions: ['.langium'],
        languageId: 'langium'
    })
};
const minimalSharedGrammarModule = {
    AstReflection: () => new ast/* LangiumGrammarAstReflection */.SV()
};
function createMinimalGrammarServices() {
    const shared = (0,dependency_injection/* inject */.f3)((0,default_module/* createDefaultSharedCoreModule */.T)(file_system_provider/* EmptyFileSystem */.u), minimalSharedGrammarModule);
    const grammar = (0,dependency_injection/* inject */.f3)((0,default_module/* createDefaultCoreModule */.Q)({ shared }), minimalGrammarModule);
    shared.ServiceRegistry.register(grammar);
    return grammar;
}
/**
 * Load a Langium grammar for your language from a JSON string. This is used by several services,
 * most notably the parser builder which interprets the grammar to create a parser.
 */
function loadGrammarFromJson(json) {
    var _a;
    const services = createMinimalGrammarServices();
    const astNode = services.serializer.JsonSerializer.deserialize(json);
    services.shared.workspace.LangiumDocumentFactory.fromModel(astNode, esm/* URI */.o.parse(`memory://${(_a = astNode.name) !== null && _a !== void 0 ? _a : 'grammar'}.langium`));
    return astNode;
}
//# sourceMappingURL=grammar-loader.js.map
// EXTERNAL MODULE: ../node_modules/langium/lib/parser/value-converter.js
var value_converter = __webpack_require__(46174);
// EXTERNAL MODULE: ../node_modules/langium/lib/parser/token-builder.js
var token_builder = __webpack_require__(35481);
;// CONCATENATED MODULE: ../node_modules/@mermaid-js/parser/dist/chunks/mermaid-parser.core/chunk-7PKI6E2E.mjs
var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });

// src/language/generated/ast.ts

var Statement = "Statement";
var Architecture = "Architecture";
function isArchitecture(item) {
  return reflection.isInstance(item, Architecture);
}
__name(isArchitecture, "isArchitecture");
var Axis = "Axis";
var Branch = "Branch";
function isBranch(item) {
  return reflection.isInstance(item, Branch);
}
__name(isBranch, "isBranch");
var Checkout = "Checkout";
var CherryPicking = "CherryPicking";
var Commit = "Commit";
function isCommit(item) {
  return reflection.isInstance(item, Commit);
}
__name(isCommit, "isCommit");
var Common = "Common";
function isCommon(item) {
  return reflection.isInstance(item, Common);
}
__name(isCommon, "isCommon");
var Curve = "Curve";
var Edge = "Edge";
var Entry = "Entry";
var GitGraph = "GitGraph";
function isGitGraph(item) {
  return reflection.isInstance(item, GitGraph);
}
__name(isGitGraph, "isGitGraph");
var Group = "Group";
var Info = "Info";
function isInfo(item) {
  return reflection.isInstance(item, Info);
}
__name(isInfo, "isInfo");
var Junction = "Junction";
var Merge = "Merge";
function isMerge(item) {
  return reflection.isInstance(item, Merge);
}
__name(isMerge, "isMerge");
var Option = "Option";
var Packet = "Packet";
function isPacket(item) {
  return reflection.isInstance(item, Packet);
}
__name(isPacket, "isPacket");
var PacketBlock = "PacketBlock";
function isPacketBlock(item) {
  return reflection.isInstance(item, PacketBlock);
}
__name(isPacketBlock, "isPacketBlock");
var Pie = "Pie";
function isPie(item) {
  return reflection.isInstance(item, Pie);
}
__name(isPie, "isPie");
var PieSection = "PieSection";
function isPieSection(item) {
  return reflection.isInstance(item, PieSection);
}
__name(isPieSection, "isPieSection");
var Radar = "Radar";
var Service = "Service";
var Direction = "Direction";
var MermaidAstReflection = class extends syntax_tree/* AbstractAstReflection */.$v {
  static {
    __name(this, "MermaidAstReflection");
  }
  getAllTypes() {
    return [Architecture, Axis, Branch, Checkout, CherryPicking, Commit, Common, Curve, Direction, Edge, Entry, GitGraph, Group, Info, Junction, Merge, Option, Packet, PacketBlock, Pie, PieSection, Radar, Service, Statement];
  }
  computeIsSubtype(subtype, supertype) {
    switch (subtype) {
      case Branch:
      case Checkout:
      case CherryPicking:
      case Commit:
      case Merge: {
        return this.isSubtype(Statement, supertype);
      }
      case Direction: {
        return this.isSubtype(GitGraph, supertype);
      }
      default: {
        return false;
      }
    }
  }
  getReferenceType(refInfo) {
    const referenceId = `${refInfo.container.$type}:${refInfo.property}`;
    switch (referenceId) {
      case "Entry:axis": {
        return Axis;
      }
      default: {
        throw new Error(`${referenceId} is not a valid reference id.`);
      }
    }
  }
  getTypeMetaData(type) {
    switch (type) {
      case Architecture: {
        return {
          name: Architecture,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "edges", defaultValue: [] },
            { name: "groups", defaultValue: [] },
            { name: "junctions", defaultValue: [] },
            { name: "services", defaultValue: [] },
            { name: "title" }
          ]
        };
      }
      case Axis: {
        return {
          name: Axis,
          properties: [
            { name: "label" },
            { name: "name" }
          ]
        };
      }
      case Branch: {
        return {
          name: Branch,
          properties: [
            { name: "name" },
            { name: "order" }
          ]
        };
      }
      case Checkout: {
        return {
          name: Checkout,
          properties: [
            { name: "branch" }
          ]
        };
      }
      case CherryPicking: {
        return {
          name: CherryPicking,
          properties: [
            { name: "id" },
            { name: "parent" },
            { name: "tags", defaultValue: [] }
          ]
        };
      }
      case Commit: {
        return {
          name: Commit,
          properties: [
            { name: "id" },
            { name: "message" },
            { name: "tags", defaultValue: [] },
            { name: "type" }
          ]
        };
      }
      case Common: {
        return {
          name: Common,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "title" }
          ]
        };
      }
      case Curve: {
        return {
          name: Curve,
          properties: [
            { name: "entries", defaultValue: [] },
            { name: "label" },
            { name: "name" }
          ]
        };
      }
      case Edge: {
        return {
          name: Edge,
          properties: [
            { name: "lhsDir" },
            { name: "lhsGroup", defaultValue: false },
            { name: "lhsId" },
            { name: "lhsInto", defaultValue: false },
            { name: "rhsDir" },
            { name: "rhsGroup", defaultValue: false },
            { name: "rhsId" },
            { name: "rhsInto", defaultValue: false },
            { name: "title" }
          ]
        };
      }
      case Entry: {
        return {
          name: Entry,
          properties: [
            { name: "axis" },
            { name: "value" }
          ]
        };
      }
      case GitGraph: {
        return {
          name: GitGraph,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "statements", defaultValue: [] },
            { name: "title" }
          ]
        };
      }
      case Group: {
        return {
          name: Group,
          properties: [
            { name: "icon" },
            { name: "id" },
            { name: "in" },
            { name: "title" }
          ]
        };
      }
      case Info: {
        return {
          name: Info,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "title" }
          ]
        };
      }
      case Junction: {
        return {
          name: Junction,
          properties: [
            { name: "id" },
            { name: "in" }
          ]
        };
      }
      case Merge: {
        return {
          name: Merge,
          properties: [
            { name: "branch" },
            { name: "id" },
            { name: "tags", defaultValue: [] },
            { name: "type" }
          ]
        };
      }
      case Option: {
        return {
          name: Option,
          properties: [
            { name: "name" },
            { name: "value", defaultValue: false }
          ]
        };
      }
      case Packet: {
        return {
          name: Packet,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "blocks", defaultValue: [] },
            { name: "title" }
          ]
        };
      }
      case PacketBlock: {
        return {
          name: PacketBlock,
          properties: [
            { name: "end" },
            { name: "label" },
            { name: "start" }
          ]
        };
      }
      case Pie: {
        return {
          name: Pie,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "sections", defaultValue: [] },
            { name: "showData", defaultValue: false },
            { name: "title" }
          ]
        };
      }
      case PieSection: {
        return {
          name: PieSection,
          properties: [
            { name: "label" },
            { name: "value" }
          ]
        };
      }
      case Radar: {
        return {
          name: Radar,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "axes", defaultValue: [] },
            { name: "curves", defaultValue: [] },
            { name: "options", defaultValue: [] },
            { name: "title" }
          ]
        };
      }
      case Service: {
        return {
          name: Service,
          properties: [
            { name: "icon" },
            { name: "iconText" },
            { name: "id" },
            { name: "in" },
            { name: "title" }
          ]
        };
      }
      case Direction: {
        return {
          name: Direction,
          properties: [
            { name: "accDescr" },
            { name: "accTitle" },
            { name: "dir" },
            { name: "statements", defaultValue: [] },
            { name: "title" }
          ]
        };
      }
      default: {
        return {
          name: type,
          properties: []
        };
      }
    }
  }
};
var reflection = new MermaidAstReflection();

// src/language/generated/grammar.ts

var loadedInfoGrammar;
var InfoGrammar = /* @__PURE__ */ __name(() => loadedInfoGrammar ?? (loadedInfoGrammar = loadGrammarFromJson('{"$type":"Grammar","isDeclared":true,"name":"Info","imports":[],"rules":[{"$type":"ParserRule","entry":true,"name":"Info","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[],"cardinality":"*"},{"$type":"Keyword","value":"info"},{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[],"cardinality":"*"},{"$type":"Group","elements":[{"$type":"Keyword","value":"showInfo"},{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[],"cardinality":"*"}],"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[],"cardinality":"?"}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false}],"definesHiddenTokens":false,"hiddenTokens":[],"interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]}],"types":[],"usedGrammars":[]}')), "InfoGrammar");
var loadedPacketGrammar;
var PacketGrammar = /* @__PURE__ */ __name(() => loadedPacketGrammar ?? (loadedPacketGrammar = loadGrammarFromJson(`{"$type":"Grammar","isDeclared":true,"name":"Packet","imports":[],"rules":[{"$type":"ParserRule","entry":true,"name":"Packet","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"},{"$type":"Keyword","value":"packet-beta"},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]},{"$type":"Assignment","feature":"blocks","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]},"cardinality":"*"}]},{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"+"},{"$type":"Assignment","feature":"blocks","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]},"cardinality":"+"}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"}]}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"PacketBlock","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"start","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":"-"},{"$type":"Assignment","feature":"end","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}}],"cardinality":"?"},{"$type":"Keyword","value":":"},{"$type":"Assignment","feature":"label","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"INT","type":{"$type":"ReturnType","name":"number"},"definition":{"$type":"RegexToken","regex":"/0|[1-9][0-9]*/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"STRING","definition":{"$type":"RegexToken","regex":"/\\"[^\\"]*\\"|'[^']*'/"},"fragment":false,"hidden":false},{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@7"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@8"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@9"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false}],"definesHiddenTokens":false,"hiddenTokens":[],"interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]}],"types":[],"usedGrammars":[]}`)), "PacketGrammar");
var loadedPieGrammar;
var PieGrammar = /* @__PURE__ */ __name(() => loadedPieGrammar ?? (loadedPieGrammar = loadGrammarFromJson('{"$type":"Grammar","isDeclared":true,"name":"Pie","imports":[],"rules":[{"$type":"ParserRule","entry":true,"name":"Pie","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"},{"$type":"Keyword","value":"pie"},{"$type":"Assignment","feature":"showData","operator":"?=","terminal":{"$type":"Keyword","value":"showData"},"cardinality":"?"},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]},{"$type":"Assignment","feature":"sections","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]},"cardinality":"*"}]},{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"+"},{"$type":"Assignment","feature":"sections","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]},"cardinality":"+"}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"*"}]}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"PieSection","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"label","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}},{"$type":"Keyword","value":":"},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"PIE_SECTION_LABEL","definition":{"$type":"RegexToken","regex":"/\\"[^\\"]+\\"/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"PIE_SECTION_VALUE","type":{"$type":"ReturnType","name":"number"},"definition":{"$type":"RegexToken","regex":"/(0|[1-9][0-9]*)(\\\\.[0-9]+)?/"},"fragment":false,"hidden":false},{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@7"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@8"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@9"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false}],"definesHiddenTokens":false,"hiddenTokens":[],"interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]}],"types":[],"usedGrammars":[]}')), "PieGrammar");
var loadedArchitectureGrammar;
var ArchitectureGrammar = /* @__PURE__ */ __name(() => loadedArchitectureGrammar ?? (loadedArchitectureGrammar = loadGrammarFromJson('{"$type":"Grammar","isDeclared":true,"name":"Architecture","imports":[],"rules":[{"$type":"ParserRule","entry":true,"name":"Architecture","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[],"cardinality":"*"},{"$type":"Keyword","value":"architecture-beta"},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@16"},"arguments":[]}]},{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[],"cardinality":"*"}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[],"cardinality":"*"}]}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"Statement","definition":{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"groups","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}},{"$type":"Assignment","feature":"services","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@6"},"arguments":[]}},{"$type":"Assignment","feature":"junctions","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@7"},"arguments":[]}},{"$type":"Assignment","feature":"edges","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@8"},"arguments":[]}}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"LeftPort","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":":"},{"$type":"Assignment","feature":"lhsDir","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@9"},"arguments":[]}}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"RightPort","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"rhsDir","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@9"},"arguments":[]}},{"$type":"Keyword","value":":"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"Arrow","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]},{"$type":"Assignment","feature":"lhsInto","operator":"?=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@15"},"arguments":[]},"cardinality":"?"},{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"--"},{"$type":"Group","elements":[{"$type":"Keyword","value":"-"},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]}},{"$type":"Keyword","value":"-"}]}]},{"$type":"Assignment","feature":"rhsInto","operator":"?=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@15"},"arguments":[]},"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Group","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"group"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}},{"$type":"Assignment","feature":"icon","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@12"},"arguments":[]},"cardinality":"?"},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]},"cardinality":"?"},{"$type":"Group","elements":[{"$type":"Keyword","value":"in"},{"$type":"Assignment","feature":"in","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}}],"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Service","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"service"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}},{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"iconText","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@11"},"arguments":[]}},{"$type":"Assignment","feature":"icon","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@12"},"arguments":[]}}],"cardinality":"?"},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]},"cardinality":"?"},{"$type":"Group","elements":[{"$type":"Keyword","value":"in"},{"$type":"Assignment","feature":"in","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}}],"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Junction","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"junction"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":"in"},{"$type":"Assignment","feature":"in","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}}],"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Edge","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"lhsId","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}},{"$type":"Assignment","feature":"lhsGroup","operator":"?=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@14"},"arguments":[]},"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]},{"$type":"Assignment","feature":"rhsId","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@10"},"arguments":[]}},{"$type":"Assignment","feature":"rhsGroup","operator":"?=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@14"},"arguments":[]},"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"ARROW_DIRECTION","definition":{"$type":"TerminalAlternatives","elements":[{"$type":"TerminalAlternatives","elements":[{"$type":"TerminalAlternatives","elements":[{"$type":"CharacterRange","left":{"$type":"Keyword","value":"L"}},{"$type":"CharacterRange","left":{"$type":"Keyword","value":"R"}}]},{"$type":"CharacterRange","left":{"$type":"Keyword","value":"T"}}]},{"$type":"CharacterRange","left":{"$type":"Keyword","value":"B"}}]},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARCH_ID","definition":{"$type":"RegexToken","regex":"/[\\\\w]+/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARCH_TEXT_ICON","definition":{"$type":"RegexToken","regex":"/\\\\(\\"[^\\"]+\\"\\\\)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARCH_ICON","definition":{"$type":"RegexToken","regex":"/\\\\([\\\\w-:]+\\\\)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARCH_TITLE","definition":{"$type":"RegexToken","regex":"/\\\\[[\\\\w ]+\\\\]/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARROW_GROUP","definition":{"$type":"RegexToken","regex":"/\\\\{group\\\\}/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ARROW_INTO","definition":{"$type":"RegexToken","regex":"/<|>/"},"fragment":false,"hidden":false},{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@19"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@21"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false}],"definesHiddenTokens":false,"hiddenTokens":[],"interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]}],"types":[],"usedGrammars":[]}')), "ArchitectureGrammar");
var loadedGitGraphGrammar;
var GitGraphGrammar = /* @__PURE__ */ __name(() => loadedGitGraphGrammar ?? (loadedGitGraphGrammar = loadGrammarFromJson(`{"$type":"Grammar","isDeclared":true,"name":"GitGraph","interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]}],"rules":[{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false},{"$type":"ParserRule","entry":true,"name":"GitGraph","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"gitGraph"},{"$type":"Group","elements":[{"$type":"Keyword","value":"gitGraph"},{"$type":"Keyword","value":":"}]},{"$type":"Keyword","value":"gitGraph:"},{"$type":"Group","elements":[{"$type":"Keyword","value":"gitGraph"},{"$type":"RuleCall","rule":{"$ref":"#/rules@12"},"arguments":[]},{"$type":"Keyword","value":":"}]}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@0"},"arguments":[]},{"$type":"Assignment","feature":"statements","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@11"},"arguments":[]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}],"cardinality":"*"}]}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Statement","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@14"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@15"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@16"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Direction","definition":{"$type":"Assignment","feature":"dir","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"LR"},{"$type":"Keyword","value":"TB"},{"$type":"Keyword","value":"BT"}]}},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Commit","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"commit"},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"Keyword","value":"id:"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"msg:","cardinality":"?"},{"$type":"Assignment","feature":"message","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"tag:"},{"$type":"Assignment","feature":"tags","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"type:"},{"$type":"Assignment","feature":"type","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"NORMAL"},{"$type":"Keyword","value":"REVERSE"},{"$type":"Keyword","value":"HIGHLIGHT"}]}}]}],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Branch","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"branch"},{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@19"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}]}},{"$type":"Group","elements":[{"$type":"Keyword","value":"order:"},{"$type":"Assignment","feature":"order","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}}],"cardinality":"?"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Merge","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"merge"},{"$type":"Assignment","feature":"branch","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@19"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}]}},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"Keyword","value":"id:"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"tag:"},{"$type":"Assignment","feature":"tags","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"type:"},{"$type":"Assignment","feature":"type","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"NORMAL"},{"$type":"Keyword","value":"REVERSE"},{"$type":"Keyword","value":"HIGHLIGHT"}]}}]}],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Checkout","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"checkout"},{"$type":"Keyword","value":"switch"}]},{"$type":"Assignment","feature":"branch","operator":"=","terminal":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@19"},"arguments":[]},{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"CherryPicking","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"cherry-pick"},{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"Keyword","value":"id:"},{"$type":"Assignment","feature":"id","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"tag:"},{"$type":"Assignment","feature":"tags","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"parent:"},{"$type":"Assignment","feature":"parent","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]}],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"INT","type":{"$type":"ReturnType","name":"number"},"definition":{"$type":"RegexToken","regex":"/[0-9]+(?=\\\\s)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ID","type":{"$type":"ReturnType","name":"string"},"definition":{"$type":"RegexToken","regex":"/\\\\w([-\\\\./\\\\w]*[-\\\\w])?/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"STRING","definition":{"$type":"RegexToken","regex":"/\\"[^\\"]*\\"|'[^']*'/"},"fragment":false,"hidden":false}],"definesHiddenTokens":false,"hiddenTokens":[],"imports":[],"types":[],"usedGrammars":[]}`)), "GitGraphGrammar");
var loadedRadarGrammar;
var RadarGrammar = /* @__PURE__ */ __name(() => loadedRadarGrammar ?? (loadedRadarGrammar = loadGrammarFromJson(`{"$type":"Grammar","isDeclared":true,"name":"Radar","interfaces":[{"$type":"Interface","name":"Common","attributes":[{"$type":"TypeAttribute","name":"accDescr","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"accTitle","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}},{"$type":"TypeAttribute","name":"title","isOptional":true,"type":{"$type":"SimpleType","primitiveType":"string"}}],"superTypes":[]},{"$type":"Interface","name":"Entry","attributes":[{"$type":"TypeAttribute","name":"axis","isOptional":true,"type":{"$type":"ReferenceType","referenceType":{"$type":"SimpleType","typeRef":{"$ref":"#/rules@12"}}}},{"$type":"TypeAttribute","name":"value","type":{"$type":"SimpleType","primitiveType":"number"},"isOptional":false}],"superTypes":[]}],"rules":[{"$type":"ParserRule","fragment":true,"name":"TitleAndAccessibilities","definition":{"$type":"Group","elements":[{"$type":"Alternatives","elements":[{"$type":"Assignment","feature":"accDescr","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@3"},"arguments":[]}},{"$type":"Assignment","feature":"accTitle","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@4"},"arguments":[]}},{"$type":"Assignment","feature":"title","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@5"},"arguments":[]}}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@1"},"arguments":[]}],"cardinality":"+"},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"EOL","dataType":"string","definition":{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"+"},{"$type":"EndOfFile"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NEWLINE","definition":{"$type":"RegexToken","regex":"/\\\\r?\\\\n/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_DESCR","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accDescr(?:[\\\\t ]*:([^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)|\\\\s*{([^}]*)})/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ACC_TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*accTitle[\\\\t ]*:(?:[^\\\\n\\\\r]*?(?=%%)|[^\\\\n\\\\r]*)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"TITLE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*title(?:[\\\\t ][^\\\\n\\\\r]*?(?=%%)|[\\\\t ][^\\\\n\\\\r]*|)/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","hidden":true,"name":"WHITESPACE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]+/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"YAML","definition":{"$type":"RegexToken","regex":"/---[\\\\t ]*\\\\r?\\\\n(?:[\\\\S\\\\s]*?\\\\r?\\\\n)?---(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"DIRECTIVE","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%{[\\\\S\\\\s]*?}%%(?:\\\\r?\\\\n|(?!\\\\S))/"},"fragment":false},{"$type":"TerminalRule","hidden":true,"name":"SINGLE_LINE_COMMENT","definition":{"$type":"RegexToken","regex":"/[\\\\t ]*%%[^\\\\n\\\\r]*/"},"fragment":false},{"$type":"ParserRule","entry":true,"name":"Radar","definition":{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Alternatives","elements":[{"$type":"Keyword","value":"radar-beta"},{"$type":"Keyword","value":"radar-beta:"},{"$type":"Group","elements":[{"$type":"Keyword","value":"radar-beta"},{"$type":"Keyword","value":":"}]}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Alternatives","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@0"},"arguments":[]},{"$type":"Group","elements":[{"$type":"Keyword","value":"axis"},{"$type":"Assignment","feature":"axes","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@12"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":","},{"$type":"Assignment","feature":"axes","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@12"},"arguments":[]}}],"cardinality":"*"}]},{"$type":"Group","elements":[{"$type":"Keyword","value":"curve"},{"$type":"Assignment","feature":"curves","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":","},{"$type":"Assignment","feature":"curves","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@13"},"arguments":[]}}],"cardinality":"*"}]},{"$type":"Group","elements":[{"$type":"Assignment","feature":"options","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":","},{"$type":"Assignment","feature":"options","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@17"},"arguments":[]}}],"cardinality":"*"}]},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[]}],"cardinality":"*"}]},"definesHiddenTokens":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"Label","definition":{"$type":"Group","elements":[{"$type":"Keyword","value":"["},{"$type":"Assignment","feature":"label","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@22"},"arguments":[]}},{"$type":"Keyword","value":"]"}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Axis","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@21"},"arguments":[]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@11"},"arguments":[],"cardinality":"?"}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Curve","definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@21"},"arguments":[]}},{"$type":"RuleCall","rule":{"$ref":"#/rules@11"},"arguments":[],"cardinality":"?"},{"$type":"Keyword","value":"{"},{"$type":"RuleCall","rule":{"$ref":"#/rules@14"},"arguments":[]},{"$type":"Keyword","value":"}"}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","fragment":true,"name":"Entries","definition":{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Assignment","feature":"entries","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@16"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":","},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Assignment","feature":"entries","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@16"},"arguments":[]}}],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"}]},{"$type":"Group","elements":[{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Assignment","feature":"entries","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@15"},"arguments":[]}},{"$type":"Group","elements":[{"$type":"Keyword","value":","},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"},{"$type":"Assignment","feature":"entries","operator":"+=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@15"},"arguments":[]}}],"cardinality":"*"},{"$type":"RuleCall","rule":{"$ref":"#/rules@2"},"arguments":[],"cardinality":"*"}]}]},"definesHiddenTokens":false,"entry":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"DetailedEntry","returnType":{"$ref":"#/interfaces@1"},"definition":{"$type":"Group","elements":[{"$type":"Assignment","feature":"axis","operator":"=","terminal":{"$type":"CrossReference","type":{"$ref":"#/rules@12"},"terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@21"},"arguments":[]},"deprecatedSyntax":false}},{"$type":"Keyword","value":":","cardinality":"?"},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"NumberEntry","returnType":{"$ref":"#/interfaces@1"},"definition":{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"ParserRule","name":"Option","definition":{"$type":"Alternatives","elements":[{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Keyword","value":"showLegend"}},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@19"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Keyword","value":"ticks"}},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Keyword","value":"max"}},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Keyword","value":"min"}},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@18"},"arguments":[]}}]},{"$type":"Group","elements":[{"$type":"Assignment","feature":"name","operator":"=","terminal":{"$type":"Keyword","value":"graticule"}},{"$type":"Assignment","feature":"value","operator":"=","terminal":{"$type":"RuleCall","rule":{"$ref":"#/rules@20"},"arguments":[]}}]}]},"definesHiddenTokens":false,"entry":false,"fragment":false,"hiddenTokens":[],"parameters":[],"wildcard":false},{"$type":"TerminalRule","name":"NUMBER","type":{"$type":"ReturnType","name":"number"},"definition":{"$type":"RegexToken","regex":"/(0|[1-9][0-9]*)(\\\\.[0-9]+)?/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"BOOLEAN","type":{"$type":"ReturnType","name":"boolean"},"definition":{"$type":"TerminalAlternatives","elements":[{"$type":"CharacterRange","left":{"$type":"Keyword","value":"true"}},{"$type":"CharacterRange","left":{"$type":"Keyword","value":"false"}}]},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"GRATICULE","type":{"$type":"ReturnType","name":"string"},"definition":{"$type":"TerminalAlternatives","elements":[{"$type":"CharacterRange","left":{"$type":"Keyword","value":"circle"}},{"$type":"CharacterRange","left":{"$type":"Keyword","value":"polygon"}}]},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"ID","type":{"$type":"ReturnType","name":"string"},"definition":{"$type":"RegexToken","regex":"/[a-zA-Z_][a-zA-Z0-9\\\\-_]*/"},"fragment":false,"hidden":false},{"$type":"TerminalRule","name":"STRING","definition":{"$type":"RegexToken","regex":"/\\"[^\\"]*\\"|'[^']*'/"},"fragment":false,"hidden":false}],"definesHiddenTokens":false,"hiddenTokens":[],"imports":[],"types":[],"usedGrammars":[]}`)), "RadarGrammar");

// src/language/generated/module.ts
var InfoLanguageMetaData = {
  languageId: "info",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var PacketLanguageMetaData = {
  languageId: "packet",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var PieLanguageMetaData = {
  languageId: "pie",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var ArchitectureLanguageMetaData = {
  languageId: "architecture",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var GitGraphLanguageMetaData = {
  languageId: "gitGraph",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var RadarLanguageMetaData = {
  languageId: "radar",
  fileExtensions: [".mmd", ".mermaid"],
  caseInsensitive: false,
  mode: "production"
};
var MermaidGeneratedSharedModule = {
  AstReflection: /* @__PURE__ */ __name(() => new MermaidAstReflection(), "AstReflection")
};
var InfoGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => InfoGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => InfoLanguageMetaData, "LanguageMetaData"),
  parser: {}
};
var PacketGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => PacketGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => PacketLanguageMetaData, "LanguageMetaData"),
  parser: {}
};
var PieGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => PieGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => PieLanguageMetaData, "LanguageMetaData"),
  parser: {}
};
var ArchitectureGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => ArchitectureGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => ArchitectureLanguageMetaData, "LanguageMetaData"),
  parser: {}
};
var GitGraphGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => GitGraphGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => GitGraphLanguageMetaData, "LanguageMetaData"),
  parser: {}
};
var RadarGeneratedModule = {
  Grammar: /* @__PURE__ */ __name(() => RadarGrammar(), "Grammar"),
  LanguageMetaData: /* @__PURE__ */ __name(() => RadarLanguageMetaData, "LanguageMetaData"),
  parser: {}
};

// src/language/common/valueConverter.ts


// src/language/common/matcher.ts
var accessibilityDescrRegex = /accDescr(?:[\t ]*:([^\n\r]*)|\s*{([^}]*)})/;
var accessibilityTitleRegex = /accTitle[\t ]*:([^\n\r]*)/;
var titleRegex = /title([\t ][^\n\r]*|)/;

// src/language/common/valueConverter.ts
var rulesRegexes = {
  ACC_DESCR: accessibilityDescrRegex,
  ACC_TITLE: accessibilityTitleRegex,
  TITLE: titleRegex
};
var AbstractMermaidValueConverter = class extends value_converter/* DefaultValueConverter */.t {
  static {
    __name(this, "AbstractMermaidValueConverter");
  }
  runConverter(rule, input, cstNode) {
    let value = this.runCommonConverter(rule, input, cstNode);
    if (value === void 0) {
      value = this.runCustomConverter(rule, input, cstNode);
    }
    if (value === void 0) {
      return super.runConverter(rule, input, cstNode);
    }
    return value;
  }
  runCommonConverter(rule, input, _cstNode) {
    const regex = rulesRegexes[rule.name];
    if (regex === void 0) {
      return void 0;
    }
    const match = regex.exec(input);
    if (match === null) {
      return void 0;
    }
    if (match[1] !== void 0) {
      return match[1].trim().replace(/[\t ]{2,}/gm, " ");
    }
    if (match[2] !== void 0) {
      return match[2].replace(/^\s*/gm, "").replace(/\s+$/gm, "").replace(/[\t ]{2,}/gm, " ").replace(/[\n\r]{2,}/gm, "\n");
    }
    return void 0;
  }
};
var CommonValueConverter = class extends AbstractMermaidValueConverter {
  static {
    __name(this, "CommonValueConverter");
  }
  runCustomConverter(_rule, _input, _cstNode) {
    return void 0;
  }
};

// src/language/common/tokenBuilder.ts

var AbstractMermaidTokenBuilder = class extends token_builder/* DefaultTokenBuilder */.P {
  static {
    __name(this, "AbstractMermaidTokenBuilder");
  }
  constructor(keywords) {
    super();
    this.keywords = new Set(keywords);
  }
  buildKeywordTokens(rules, terminalTokens, options) {
    const tokenTypes = super.buildKeywordTokens(rules, terminalTokens, options);
    tokenTypes.forEach((tokenType) => {
      if (this.keywords.has(tokenType.name) && tokenType.PATTERN !== void 0) {
        tokenType.PATTERN = new RegExp(tokenType.PATTERN.toString() + "(?:(?=%%)|(?!\\S))");
      }
    });
    return tokenTypes;
  }
};
var CommonTokenBuilder = class extends AbstractMermaidTokenBuilder {
  static {
    __name(this, "CommonTokenBuilder");
  }
};




/***/ }),

/***/ 95626:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   i: () => (/* binding */ createArchitectureServices)
/* harmony export */ });
/* unused harmony export ArchitectureModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/architecture/module.ts


// src/language/architecture/tokenBuilder.ts
var ArchitectureTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "ArchitectureTokenBuilder");
  }
  constructor() {
    super(["architecture"]);
  }
};

// src/language/architecture/valueConverter.ts
var ArchitectureValueConverter = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidValueConverter */ .kb {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "ArchitectureValueConverter");
  }
  runCustomConverter(rule, input, _cstNode) {
    if (rule.name === "ARCH_ICON") {
      return input.replace(/[()]/g, "").trim();
    } else if (rule.name === "ARCH_TEXT_ICON") {
      return input.replace(/["()]/g, "");
    } else if (rule.name === "ARCH_TITLE") {
      return input.replace(/[[\]]/g, "").trim();
    }
    return void 0;
  }
};

// src/language/architecture/module.ts
var ArchitectureModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new ArchitectureTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new ArchitectureValueConverter(), "ValueConverter")
  }
};
function createArchitectureServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const Architecture = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .ArchitectureGeneratedModule */ .Qr,
    ArchitectureModule
  );
  shared.ServiceRegistry.register(Architecture);
  return { shared, Architecture };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createArchitectureServices, "createArchitectureServices");




/***/ }),

/***/ 63809:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   M: () => (/* binding */ createInfoServices)
/* harmony export */ });
/* unused harmony export InfoModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/info/module.ts


// src/language/info/tokenBuilder.ts
var InfoTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "InfoTokenBuilder");
  }
  constructor() {
    super(["info", "showInfo"]);
  }
};

// src/language/info/module.ts
var InfoModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new InfoTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .CommonValueConverter */ .nr(), "ValueConverter")
  }
};
function createInfoServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const Info = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .InfoGeneratedModule */ .F_,
    InfoModule
  );
  shared.ServiceRegistry.register(Info);
  return { shared, Info };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createInfoServices, "createInfoServices");




/***/ }),

/***/ 76150:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   l: () => (/* binding */ createPieServices)
/* harmony export */ });
/* unused harmony export PieModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/pie/module.ts


// src/language/pie/tokenBuilder.ts
var PieTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "PieTokenBuilder");
  }
  constructor() {
    super(["pie", "showData"]);
  }
};

// src/language/pie/valueConverter.ts
var PieValueConverter = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidValueConverter */ .kb {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "PieValueConverter");
  }
  runCustomConverter(rule, input, _cstNode) {
    if (rule.name !== "PIE_SECTION_LABEL") {
      return void 0;
    }
    return input.replace(/"/g, "").trim();
  }
};

// src/language/pie/module.ts
var PieModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new PieTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new PieValueConverter(), "ValueConverter")
  }
};
function createPieServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const Pie = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .PieGeneratedModule */ .WH,
    PieModule
  );
  shared.ServiceRegistry.register(Pie);
  return { shared, Pie };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createPieServices, "createPieServices");




/***/ }),

/***/ 89053:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   g: () => (/* binding */ createPacketServices)
/* harmony export */ });
/* unused harmony export PacketModule */
/* harmony import */ var _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60960);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(44014);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(81210);
/* harmony import */ var langium__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(73001);


// src/language/packet/module.ts


// src/language/packet/tokenBuilder.ts
var PacketTokenBuilder = class extends _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .AbstractMermaidTokenBuilder */ .T7 {
  static {
    (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "PacketTokenBuilder");
  }
  constructor() {
    super(["packet-beta"]);
  }
};

// src/language/packet/module.ts
var PacketModule = {
  parser: {
    TokenBuilder: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new PacketTokenBuilder(), "TokenBuilder"),
    ValueConverter: /* @__PURE__ */ (0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => new _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .CommonValueConverter */ .nr(), "ValueConverter")
  }
};
function createPacketServices(context = langium__WEBPACK_IMPORTED_MODULE_1__/* .EmptyFileSystem */ .u) {
  const shared = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultSharedCoreModule */ .T)(context),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .MermaidGeneratedSharedModule */ .GS
  );
  const Packet = (0,langium__WEBPACK_IMPORTED_MODULE_2__/* .inject */ .f3)(
    (0,langium__WEBPACK_IMPORTED_MODULE_3__/* .createDefaultCoreModule */ .Q)({ shared }),
    _chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .PacketGeneratedModule */ .bb,
    PacketModule
  );
  shared.ServiceRegistry.register(Packet);
  return { shared, Packet };
}
(0,_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(createPacketServices, "createPacketServices");




/***/ }),

/***/ 13197:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Qc: () => (/* binding */ parse)
/* harmony export */ });
/* unused harmony export MermaidParseError */
/* harmony import */ var _chunks_mermaid_parser_core_chunk_2NYFTIL2_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(20078);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_EXZZNE6F_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(63809);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_V4Q32G6S_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(89053);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_ROXG7S4E_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(76150);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_C4OEIS7N_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(95626);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_2O5ZK7RR_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(51991);
/* harmony import */ var _chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(60960);








// src/parse.ts
var parsers = {};
var initializers = {
  info: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createInfoServices: createInfoServices2 } = await __webpack_require__.e(/* import() */ 4828).then(__webpack_require__.bind(__webpack_require__, 14828));
    const parser = createInfoServices2().Info.parser.LangiumParser;
    parsers.info = parser;
  }, "info"),
  packet: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createPacketServices: createPacketServices2 } = await __webpack_require__.e(/* import() */ 339).then(__webpack_require__.bind(__webpack_require__, 90339));
    const parser = createPacketServices2().Packet.parser.LangiumParser;
    parsers.packet = parser;
  }, "packet"),
  pie: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createPieServices: createPieServices2 } = await __webpack_require__.e(/* import() */ 6561).then(__webpack_require__.bind(__webpack_require__, 96561));
    const parser = createPieServices2().Pie.parser.LangiumParser;
    parsers.pie = parser;
  }, "pie"),
  architecture: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createArchitectureServices: createArchitectureServices2 } = await __webpack_require__.e(/* import() */ 6800).then(__webpack_require__.bind(__webpack_require__, 96800));
    const parser = createArchitectureServices2().Architecture.parser.LangiumParser;
    parsers.architecture = parser;
  }, "architecture"),
  gitGraph: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createGitGraphServices: createGitGraphServices2 } = await __webpack_require__.e(/* import() */ 9322).then(__webpack_require__.bind(__webpack_require__, 59322));
    const parser = createGitGraphServices2().GitGraph.parser.LangiumParser;
    parsers.gitGraph = parser;
  }, "gitGraph"),
  radar: /* @__PURE__ */ (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async () => {
    const { createRadarServices: createRadarServices2 } = await __webpack_require__.e(/* import() */ 3076).then(__webpack_require__.bind(__webpack_require__, 73076));
    const parser = createRadarServices2().Radar.parser.LangiumParser;
    parsers.radar = parser;
  }, "radar")
};
async function parse(diagramType, text) {
  const initializer = initializers[diagramType];
  if (!initializer) {
    throw new Error(`Unknown diagram type: ${diagramType}`);
  }
  if (!parsers[diagramType]) {
    await initializer();
  }
  const parser = parsers[diagramType];
  const result = parser.parse(text);
  if (result.lexerErrors.length > 0 || result.parserErrors.length > 0) {
    throw new MermaidParseError(result);
  }
  return result.value;
}
(0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(parse, "parse");
var MermaidParseError = class extends Error {
  constructor(result) {
    const lexerErrors = result.lexerErrors.map((err) => err.message).join("\n");
    const parserErrors = result.parserErrors.map((err) => err.message).join("\n");
    super(`Parsing failed: ${lexerErrors} ${parserErrors}`);
    this.result = result;
  }
  static {
    (0,_chunks_mermaid_parser_core_chunk_7PKI6E2E_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(this, "MermaidParseError");
  }
};



/***/ }),

/***/ 34326:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  ue: () => (/* reexport */ Alternation),
  _o: () => (/* reexport */ EMPTY_ALT),
  sd: () => (/* reexport */ EOF),
  nu: () => (/* reexport */ EmbeddedActionsParser),
  dV: () => (/* reexport */ LLkLookaheadStrategy),
  hW: () => (/* reexport */ Lexer),
  Sj: () => (/* reexport */ NonTerminal),
  Wx: () => (/* reexport */ Option),
  hI: () => (/* reexport */ Repetition),
  ej: () => (/* reexport */ RepetitionMandatory),
  fK: () => (/* reexport */ RepetitionMandatoryWithSeparator),
  pT: () => (/* reexport */ RepetitionWithSeparator),
  oI: () => (/* reexport */ Terminal),
  ZW: () => (/* reexport */ defaultLexerErrorProvider),
  Hs: () => (/* reexport */ defaultParserErrorProvider),
  oC: () => (/* reexport */ getLookaheadPaths),
  l$: () => (/* reexport */ tokens_public_tokenLabel),
  ol: () => (/* reexport */ tokenMatcher)
});

// UNUSED EXPORTS: Alternative, CstParser, EarlyExitException, GAstVisitor, LexerDefinitionErrorType, MismatchedTokenException, NoViableAltException, NotAllInputParsedException, Parser, ParserDefinitionErrorType, Rule, VERSION, clearCache, createSyntaxDiagramsCode, createToken, createTokenInstance, generateCstDts, isRecognitionException, serializeGrammar, serializeProduction, tokenName

// EXTERNAL MODULE: ../node_modules/lodash-es/forEach.js
var forEach = __webpack_require__(21845);
// EXTERNAL MODULE: ../node_modules/lodash-es/values.js + 1 modules
var values = __webpack_require__(88873);
// EXTERNAL MODULE: ../node_modules/lodash-es/isEmpty.js
var isEmpty = __webpack_require__(66400);
// EXTERNAL MODULE: ../node_modules/lodash-es/map.js
var map = __webpack_require__(12930);
// EXTERNAL MODULE: ../node_modules/lodash-es/has.js + 1 modules
var has = __webpack_require__(49493);
// EXTERNAL MODULE: ../node_modules/lodash-es/clone.js
var lodash_es_clone = __webpack_require__(20190);
;// CONCATENATED MODULE: ../node_modules/@chevrotain/utils/lib/src/to-fast-properties.js
// based on: https://github.com/petkaantonov/bluebird/blob/b97c0d2d487e8c5076e8bd897e0dcd4622d31846/src/util.js#L201-L216
function toFastProperties(toBecomeFast) {
    function FakeConstructor() { }
    // If our object is used as a constructor, it would receive
    FakeConstructor.prototype = toBecomeFast;
    const fakeInstance = new FakeConstructor();
    function fakeAccess() {
        return typeof fakeInstance.bar;
    }
    // help V8 understand this is a "real" prototype by actually using
    // the fake instance.
    fakeAccess();
    fakeAccess();
    // Always true condition to suppress the Firefox warning of unreachable
    // code after a return statement.
    if (true)
        return toBecomeFast;
    // Eval prevents optimization of this method (even though this is dead code)
    // - https://esbuild.github.io/content-types/#direct-eval
    /* istanbul ignore next */
    // tslint:disable-next-line
    (0, eval)(toBecomeFast);
}
//# sourceMappingURL=to-fast-properties.js.map
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseSlice.js
/**
 * The base implementation of `_.slice` without an iteratee call guard.
 *
 * @private
 * @param {Array} array The array to slice.
 * @param {number} [start=0] The start position.
 * @param {number} [end=array.length] The end position.
 * @returns {Array} Returns the slice of `array`.
 */
function baseSlice(array, start, end) {
  var index = -1,
      length = array.length;

  if (start < 0) {
    start = -start > length ? 0 : (length + start);
  }
  end = end > length ? length : end;
  if (end < 0) {
    end += length;
  }
  length = start > end ? 0 : ((end - start) >>> 0);
  start >>>= 0;

  var result = Array(length);
  while (++index < length) {
    result[index] = array[index + start];
  }
  return result;
}

/* harmony default export */ const _baseSlice = (baseSlice);

// EXTERNAL MODULE: ../node_modules/lodash-es/toInteger.js
var toInteger = __webpack_require__(98670);
;// CONCATENATED MODULE: ../node_modules/lodash-es/drop.js



/**
 * Creates a slice of `array` with `n` elements dropped from the beginning.
 *
 * @static
 * @memberOf _
 * @since 0.5.0
 * @category Array
 * @param {Array} array The array to query.
 * @param {number} [n=1] The number of elements to drop.
 * @param- {Object} [guard] Enables use as an iteratee for methods like `_.map`.
 * @returns {Array} Returns the slice of `array`.
 * @example
 *
 * _.drop([1, 2, 3]);
 * // => [2, 3]
 *
 * _.drop([1, 2, 3], 2);
 * // => [3]
 *
 * _.drop([1, 2, 3], 5);
 * // => []
 *
 * _.drop([1, 2, 3], 0);
 * // => [1, 2, 3]
 */
function drop(array, n, guard) {
  var length = array == null ? 0 : array.length;
  if (!length) {
    return [];
  }
  n = (guard || n === undefined) ? 1 : (0,toInteger/* default */.Z)(n);
  return _baseSlice(array, n < 0 ? 0 : n, length);
}

/* harmony default export */ const lodash_es_drop = (drop);

// EXTERNAL MODULE: ../node_modules/lodash-es/isString.js
var isString = __webpack_require__(75732);
// EXTERNAL MODULE: ../node_modules/lodash-es/_assignValue.js
var _assignValue = __webpack_require__(15561);
// EXTERNAL MODULE: ../node_modules/lodash-es/_copyObject.js
var _copyObject = __webpack_require__(47313);
// EXTERNAL MODULE: ../node_modules/lodash-es/_createAssigner.js
var _createAssigner = __webpack_require__(40690);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArrayLike.js
var isArrayLike = __webpack_require__(69959);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isPrototype.js
var _isPrototype = __webpack_require__(89418);
// EXTERNAL MODULE: ../node_modules/lodash-es/keys.js
var keys = __webpack_require__(11723);
;// CONCATENATED MODULE: ../node_modules/lodash-es/assign.js







/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var assign_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Assigns own enumerable string keyed properties of source objects to the
 * destination object. Source objects are applied from left to right.
 * Subsequent sources overwrite property assignments of previous sources.
 *
 * **Note:** This method mutates `object` and is loosely based on
 * [`Object.assign`](https://mdn.io/Object/assign).
 *
 * @static
 * @memberOf _
 * @since 0.10.0
 * @category Object
 * @param {Object} object The destination object.
 * @param {...Object} [sources] The source objects.
 * @returns {Object} Returns `object`.
 * @see _.assignIn
 * @example
 *
 * function Foo() {
 *   this.a = 1;
 * }
 *
 * function Bar() {
 *   this.c = 3;
 * }
 *
 * Foo.prototype.b = 2;
 * Bar.prototype.d = 4;
 *
 * _.assign({ 'a': 0 }, new Foo, new Bar);
 * // => { 'a': 1, 'c': 3 }
 */
var assign_assign = (0,_createAssigner/* default */.Z)(function(object, source) {
  if ((0,_isPrototype/* default */.Z)(source) || (0,isArrayLike/* default */.Z)(source)) {
    (0,_copyObject/* default */.Z)(source, (0,keys/* default */.Z)(source), object);
    return;
  }
  for (var key in source) {
    if (assign_hasOwnProperty.call(source, key)) {
      (0,_assignValue/* default */.Z)(object, key, source[key]);
    }
  }
});

/* harmony default export */ const lodash_es_assign = (assign_assign);

// EXTERNAL MODULE: ../node_modules/lodash-es/_arrayMap.js
var _arrayMap = __webpack_require__(33043);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseIteratee.js + 15 modules
var _baseIteratee = __webpack_require__(86494);
// EXTERNAL MODULE: ../node_modules/lodash-es/_basePickBy.js + 1 modules
var _basePickBy = __webpack_require__(73338);
// EXTERNAL MODULE: ../node_modules/lodash-es/_getAllKeysIn.js
var _getAllKeysIn = __webpack_require__(5206);
;// CONCATENATED MODULE: ../node_modules/lodash-es/pickBy.js





/**
 * Creates an object composed of the `object` properties `predicate` returns
 * truthy for. The predicate is invoked with two arguments: (value, key).
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Object
 * @param {Object} object The source object.
 * @param {Function} [predicate=_.identity] The function invoked per property.
 * @returns {Object} Returns the new object.
 * @example
 *
 * var object = { 'a': 1, 'b': '2', 'c': 3 };
 *
 * _.pickBy(object, _.isNumber);
 * // => { 'a': 1, 'c': 3 }
 */
function pickBy(object, predicate) {
  if (object == null) {
    return {};
  }
  var props = (0,_arrayMap/* default */.Z)((0,_getAllKeysIn/* default */.Z)(object), function(prop) {
    return [prop];
  });
  predicate = (0,_baseIteratee/* default */.Z)(predicate);
  return (0,_basePickBy/* default */.Z)(object, props, function(value, path) {
    return predicate(value, path[0]);
  });
}

/* harmony default export */ const lodash_es_pickBy = (pickBy);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGetTag.js + 2 modules
var _baseGetTag = __webpack_require__(77070);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObjectLike.js
var isObjectLike = __webpack_require__(9615);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseIsRegExp.js



/** `Object#toString` result references. */
var regexpTag = '[object RegExp]';

/**
 * The base implementation of `_.isRegExp` without Node.js optimizations.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a regexp, else `false`.
 */
function baseIsRegExp(value) {
  return (0,isObjectLike/* default */.Z)(value) && (0,_baseGetTag/* default */.Z)(value) == regexpTag;
}

/* harmony default export */ const _baseIsRegExp = (baseIsRegExp);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseUnary.js
var _baseUnary = __webpack_require__(20274);
// EXTERNAL MODULE: ../node_modules/lodash-es/_nodeUtil.js
var _nodeUtil = __webpack_require__(53594);
;// CONCATENATED MODULE: ../node_modules/lodash-es/isRegExp.js




/* Node.js helper references. */
var nodeIsRegExp = _nodeUtil/* default */.Z && _nodeUtil/* default */.Z.isRegExp;

/**
 * Checks if `value` is classified as a `RegExp` object.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a regexp, else `false`.
 * @example
 *
 * _.isRegExp(/abc/);
 * // => true
 *
 * _.isRegExp('/abc/');
 * // => false
 */
var isRegExp = nodeIsRegExp ? (0,_baseUnary/* default */.Z)(nodeIsRegExp) : _baseIsRegExp;

/* harmony default export */ const lodash_es_isRegExp = (isRegExp);

;// CONCATENATED MODULE: ../node_modules/@chevrotain/gast/lib/src/model.js

// TODO: duplicated code to avoid extracting another sub-package -- how to avoid?
function tokenLabel(tokType) {
    if (hasTokenLabel(tokType)) {
        return tokType.LABEL;
    }
    else {
        return tokType.name;
    }
}
// TODO: duplicated code to avoid extracting another sub-package -- how to avoid?
function hasTokenLabel(obj) {
    return (0,isString/* default */.Z)(obj.LABEL) && obj.LABEL !== "";
}
class AbstractProduction {
    get definition() {
        return this._definition;
    }
    set definition(value) {
        this._definition = value;
    }
    constructor(_definition) {
        this._definition = _definition;
    }
    accept(visitor) {
        visitor.visit(this);
        (0,forEach/* default */.Z)(this.definition, (prod) => {
            prod.accept(visitor);
        });
    }
}
class NonTerminal extends AbstractProduction {
    constructor(options) {
        super([]);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
    set definition(definition) {
        // immutable
    }
    get definition() {
        if (this.referencedRule !== undefined) {
            return this.referencedRule.definition;
        }
        return [];
    }
    accept(visitor) {
        visitor.visit(this);
        // don't visit children of a reference, we will get cyclic infinite loops if we do so
    }
}
class Rule extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.orgText = "";
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class Alternative extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.ignoreAmbiguities = false;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class Option extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class RepetitionMandatory extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class RepetitionMandatoryWithSeparator extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class Repetition extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class RepetitionWithSeparator extends AbstractProduction {
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class Alternation extends AbstractProduction {
    get definition() {
        return this._definition;
    }
    set definition(value) {
        this._definition = value;
    }
    constructor(options) {
        super(options.definition);
        this.idx = 1;
        this.ignoreAmbiguities = false;
        this.hasPredicates = false;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
}
class Terminal {
    constructor(options) {
        this.idx = 1;
        lodash_es_assign(this, lodash_es_pickBy(options, (v) => v !== undefined));
    }
    accept(visitor) {
        visitor.visit(this);
    }
}
function serializeGrammar(topRules) {
    return (0,map/* default */.Z)(topRules, serializeProduction);
}
function serializeProduction(node) {
    function convertDefinition(definition) {
        return (0,map/* default */.Z)(definition, serializeProduction);
    }
    /* istanbul ignore else */
    if (node instanceof NonTerminal) {
        const serializedNonTerminal = {
            type: "NonTerminal",
            name: node.nonTerminalName,
            idx: node.idx,
        };
        if ((0,isString/* default */.Z)(node.label)) {
            serializedNonTerminal.label = node.label;
        }
        return serializedNonTerminal;
    }
    else if (node instanceof Alternative) {
        return {
            type: "Alternative",
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof Option) {
        return {
            type: "Option",
            idx: node.idx,
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof RepetitionMandatory) {
        return {
            type: "RepetitionMandatory",
            idx: node.idx,
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof RepetitionMandatoryWithSeparator) {
        return {
            type: "RepetitionMandatoryWithSeparator",
            idx: node.idx,
            separator: (serializeProduction(new Terminal({ terminalType: node.separator }))),
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof RepetitionWithSeparator) {
        return {
            type: "RepetitionWithSeparator",
            idx: node.idx,
            separator: (serializeProduction(new Terminal({ terminalType: node.separator }))),
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof Repetition) {
        return {
            type: "Repetition",
            idx: node.idx,
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof Alternation) {
        return {
            type: "Alternation",
            idx: node.idx,
            definition: convertDefinition(node.definition),
        };
    }
    else if (node instanceof Terminal) {
        const serializedTerminal = {
            type: "Terminal",
            name: node.terminalType.name,
            label: tokenLabel(node.terminalType),
            idx: node.idx,
        };
        if ((0,isString/* default */.Z)(node.label)) {
            serializedTerminal.terminalLabel = node.label;
        }
        const pattern = node.terminalType.PATTERN;
        if (node.terminalType.PATTERN) {
            serializedTerminal.pattern = lodash_es_isRegExp(pattern)
                ? pattern.source
                : pattern;
        }
        return serializedTerminal;
    }
    else if (node instanceof Rule) {
        return {
            type: "Rule",
            name: node.name,
            orgText: node.orgText,
            definition: convertDefinition(node.definition),
        };
        /* c8 ignore next 3 */
    }
    else {
        throw Error("non exhaustive match");
    }
}
//# sourceMappingURL=model.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/rest.js


/**
 *  A Grammar Walker that computes the "remaining" grammar "after" a productions in the grammar.
 */
class RestWalker {
    walk(prod, prevRest = []) {
        (0,forEach/* default */.Z)(prod.definition, (subProd, index) => {
            const currRest = lodash_es_drop(prod.definition, index + 1);
            /* istanbul ignore else */
            if (subProd instanceof NonTerminal) {
                this.walkProdRef(subProd, currRest, prevRest);
            }
            else if (subProd instanceof Terminal) {
                this.walkTerminal(subProd, currRest, prevRest);
            }
            else if (subProd instanceof Alternative) {
                this.walkFlat(subProd, currRest, prevRest);
            }
            else if (subProd instanceof Option) {
                this.walkOption(subProd, currRest, prevRest);
            }
            else if (subProd instanceof RepetitionMandatory) {
                this.walkAtLeastOne(subProd, currRest, prevRest);
            }
            else if (subProd instanceof RepetitionMandatoryWithSeparator) {
                this.walkAtLeastOneSep(subProd, currRest, prevRest);
            }
            else if (subProd instanceof RepetitionWithSeparator) {
                this.walkManySep(subProd, currRest, prevRest);
            }
            else if (subProd instanceof Repetition) {
                this.walkMany(subProd, currRest, prevRest);
            }
            else if (subProd instanceof Alternation) {
                this.walkOr(subProd, currRest, prevRest);
            }
            else {
                throw Error("non exhaustive match");
            }
        });
    }
    walkTerminal(terminal, currRest, prevRest) { }
    walkProdRef(refProd, currRest, prevRest) { }
    walkFlat(flatProd, currRest, prevRest) {
        // ABCDEF => after the D the rest is EF
        const fullOrRest = currRest.concat(prevRest);
        this.walk(flatProd, fullOrRest);
    }
    walkOption(optionProd, currRest, prevRest) {
        // ABC(DE)?F => after the (DE)? the rest is F
        const fullOrRest = currRest.concat(prevRest);
        this.walk(optionProd, fullOrRest);
    }
    walkAtLeastOne(atLeastOneProd, currRest, prevRest) {
        // ABC(DE)+F => after the (DE)+ the rest is (DE)?F
        const fullAtLeastOneRest = [
            new Option({ definition: atLeastOneProd.definition }),
        ].concat(currRest, prevRest);
        this.walk(atLeastOneProd, fullAtLeastOneRest);
    }
    walkAtLeastOneSep(atLeastOneSepProd, currRest, prevRest) {
        // ABC DE(,DE)* F => after the (,DE)+ the rest is (,DE)?F
        const fullAtLeastOneSepRest = restForRepetitionWithSeparator(atLeastOneSepProd, currRest, prevRest);
        this.walk(atLeastOneSepProd, fullAtLeastOneSepRest);
    }
    walkMany(manyProd, currRest, prevRest) {
        // ABC(DE)*F => after the (DE)* the rest is (DE)?F
        const fullManyRest = [
            new Option({ definition: manyProd.definition }),
        ].concat(currRest, prevRest);
        this.walk(manyProd, fullManyRest);
    }
    walkManySep(manySepProd, currRest, prevRest) {
        // ABC (DE(,DE)*)? F => after the (,DE)* the rest is (,DE)?F
        const fullManySepRest = restForRepetitionWithSeparator(manySepProd, currRest, prevRest);
        this.walk(manySepProd, fullManySepRest);
    }
    walkOr(orProd, currRest, prevRest) {
        // ABC(D|E|F)G => when finding the (D|E|F) the rest is G
        const fullOrRest = currRest.concat(prevRest);
        // walk all different alternatives
        (0,forEach/* default */.Z)(orProd.definition, (alt) => {
            // wrapping each alternative in a single definition wrapper
            // to avoid errors in computing the rest of that alternative in the invocation to computeInProdFollows
            // (otherwise for OR([alt1,alt2]) alt2 will be considered in 'rest' of alt1
            const prodWrapper = new Alternative({ definition: [alt] });
            this.walk(prodWrapper, fullOrRest);
        });
    }
}
function restForRepetitionWithSeparator(repSepProd, currRest, prevRest) {
    const repSepRest = [
        new Option({
            definition: [
                new Terminal({ terminalType: repSepProd.separator }),
            ].concat(repSepProd.definition),
        }),
    ];
    const fullRepSepRest = repSepRest.concat(currRest, prevRest);
    return fullRepSepRest;
}
//# sourceMappingURL=rest.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseUniq.js + 1 modules
var _baseUniq = __webpack_require__(99633);
;// CONCATENATED MODULE: ../node_modules/lodash-es/uniq.js


/**
 * Creates a duplicate-free version of an array, using
 * [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * for equality comparisons, in which only the first occurrence of each element
 * is kept. The order of result values is determined by the order they occur
 * in the array.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to inspect.
 * @returns {Array} Returns the new duplicate free array.
 * @example
 *
 * _.uniq([2, 1, 2]);
 * // => [2, 1]
 */
function uniq(array) {
  return (array && array.length) ? (0,_baseUniq/* default */.Z)(array) : [];
}

/* harmony default export */ const lodash_es_uniq = (uniq);

// EXTERNAL MODULE: ../node_modules/lodash-es/flatten.js
var flatten = __webpack_require__(28099);
// EXTERNAL MODULE: ../node_modules/lodash-es/_arraySome.js
var _arraySome = __webpack_require__(82964);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseEach.js + 1 modules
var _baseEach = __webpack_require__(77201);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseSome.js


/**
 * The base implementation of `_.some` without support for iteratee shorthands.
 *
 * @private
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} predicate The function invoked per iteration.
 * @returns {boolean} Returns `true` if any element passes the predicate check,
 *  else `false`.
 */
function baseSome(collection, predicate) {
  var result;

  (0,_baseEach/* default */.Z)(collection, function(value, index, collection) {
    result = predicate(value, index, collection);
    return !result;
  });
  return !!result;
}

/* harmony default export */ const _baseSome = (baseSome);

// EXTERNAL MODULE: ../node_modules/lodash-es/isArray.js
var isArray = __webpack_require__(64058);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isIterateeCall.js
var _isIterateeCall = __webpack_require__(47952);
;// CONCATENATED MODULE: ../node_modules/lodash-es/some.js






/**
 * Checks if `predicate` returns truthy for **any** element of `collection`.
 * Iteration is stopped once `predicate` returns truthy. The predicate is
 * invoked with three arguments: (value, index|key, collection).
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [predicate=_.identity] The function invoked per iteration.
 * @param- {Object} [guard] Enables use as an iteratee for methods like `_.map`.
 * @returns {boolean} Returns `true` if any element passes the predicate check,
 *  else `false`.
 * @example
 *
 * _.some([null, 0, 'yes', false], Boolean);
 * // => true
 *
 * var users = [
 *   { 'user': 'barney', 'active': true },
 *   { 'user': 'fred',   'active': false }
 * ];
 *
 * // The `_.matches` iteratee shorthand.
 * _.some(users, { 'user': 'barney', 'active': false });
 * // => false
 *
 * // The `_.matchesProperty` iteratee shorthand.
 * _.some(users, ['active', false]);
 * // => true
 *
 * // The `_.property` iteratee shorthand.
 * _.some(users, 'active');
 * // => true
 */
function some(collection, predicate, guard) {
  var func = (0,isArray/* default */.Z)(collection) ? _arraySome/* default */.Z : _baseSome;
  if (guard && (0,_isIterateeCall/* default */.Z)(collection, predicate, guard)) {
    predicate = undefined;
  }
  return func(collection, (0,_baseIteratee/* default */.Z)(predicate, 3));
}

/* harmony default export */ const lodash_es_some = (some);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseIndexOf.js + 2 modules
var _baseIndexOf = __webpack_require__(90393);
;// CONCATENATED MODULE: ../node_modules/lodash-es/includes.js






/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeMax = Math.max;

/**
 * Checks if `value` is in `collection`. If `collection` is a string, it's
 * checked for a substring of `value`, otherwise
 * [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * is used for equality comparisons. If `fromIndex` is negative, it's used as
 * the offset from the end of `collection`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object|string} collection The collection to inspect.
 * @param {*} value The value to search for.
 * @param {number} [fromIndex=0] The index to search from.
 * @param- {Object} [guard] Enables use as an iteratee for methods like `_.reduce`.
 * @returns {boolean} Returns `true` if `value` is found, else `false`.
 * @example
 *
 * _.includes([1, 2, 3], 1);
 * // => true
 *
 * _.includes([1, 2, 3], 1, 2);
 * // => false
 *
 * _.includes({ 'a': 1, 'b': 2 }, 1);
 * // => true
 *
 * _.includes('abcd', 'bc');
 * // => true
 */
function includes(collection, value, fromIndex, guard) {
  collection = (0,isArrayLike/* default */.Z)(collection) ? collection : (0,values/* default */.Z)(collection);
  fromIndex = (fromIndex && !guard) ? (0,toInteger/* default */.Z)(fromIndex) : 0;

  var length = collection.length;
  if (fromIndex < 0) {
    fromIndex = nativeMax(length + fromIndex, 0);
  }
  return (0,isString/* default */.Z)(collection)
    ? (fromIndex <= length && collection.indexOf(value, fromIndex) > -1)
    : (!!length && (0,_baseIndexOf/* default */.Z)(collection, value, fromIndex) > -1);
}

/* harmony default export */ const lodash_es_includes = (includes);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_arrayEvery.js
/**
 * A specialized version of `_.every` for arrays without support for
 * iteratee shorthands.
 *
 * @private
 * @param {Array} [array] The array to iterate over.
 * @param {Function} predicate The function invoked per iteration.
 * @returns {boolean} Returns `true` if all elements pass the predicate check,
 *  else `false`.
 */
function arrayEvery(array, predicate) {
  var index = -1,
      length = array == null ? 0 : array.length;

  while (++index < length) {
    if (!predicate(array[index], index, array)) {
      return false;
    }
  }
  return true;
}

/* harmony default export */ const _arrayEvery = (arrayEvery);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseEvery.js


/**
 * The base implementation of `_.every` without support for iteratee shorthands.
 *
 * @private
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} predicate The function invoked per iteration.
 * @returns {boolean} Returns `true` if all elements pass the predicate check,
 *  else `false`
 */
function baseEvery(collection, predicate) {
  var result = true;
  (0,_baseEach/* default */.Z)(collection, function(value, index, collection) {
    result = !!predicate(value, index, collection);
    return result;
  });
  return result;
}

/* harmony default export */ const _baseEvery = (baseEvery);

;// CONCATENATED MODULE: ../node_modules/lodash-es/every.js






/**
 * Checks if `predicate` returns truthy for **all** elements of `collection`.
 * Iteration is stopped once `predicate` returns falsey. The predicate is
 * invoked with three arguments: (value, index|key, collection).
 *
 * **Note:** This method returns `true` for
 * [empty collections](https://en.wikipedia.org/wiki/Empty_set) because
 * [everything is true](https://en.wikipedia.org/wiki/Vacuous_truth) of
 * elements of empty collections.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [predicate=_.identity] The function invoked per iteration.
 * @param- {Object} [guard] Enables use as an iteratee for methods like `_.map`.
 * @returns {boolean} Returns `true` if all elements pass the predicate check,
 *  else `false`.
 * @example
 *
 * _.every([true, 1, null, 'yes'], Boolean);
 * // => false
 *
 * var users = [
 *   { 'user': 'barney', 'age': 36, 'active': false },
 *   { 'user': 'fred',   'age': 40, 'active': false }
 * ];
 *
 * // The `_.matches` iteratee shorthand.
 * _.every(users, { 'user': 'barney', 'active': false });
 * // => false
 *
 * // The `_.matchesProperty` iteratee shorthand.
 * _.every(users, ['active', false]);
 * // => true
 *
 * // The `_.property` iteratee shorthand.
 * _.every(users, 'active');
 * // => false
 */
function every(collection, predicate, guard) {
  var func = (0,isArray/* default */.Z)(collection) ? _arrayEvery : _baseEvery;
  if (guard && (0,_isIterateeCall/* default */.Z)(collection, predicate, guard)) {
    predicate = undefined;
  }
  return func(collection, (0,_baseIteratee/* default */.Z)(predicate, 3));
}

/* harmony default export */ const lodash_es_every = (every);

;// CONCATENATED MODULE: ../node_modules/@chevrotain/gast/lib/src/helpers.js


function isSequenceProd(prod) {
    return (prod instanceof Alternative ||
        prod instanceof Option ||
        prod instanceof Repetition ||
        prod instanceof RepetitionMandatory ||
        prod instanceof RepetitionMandatoryWithSeparator ||
        prod instanceof RepetitionWithSeparator ||
        prod instanceof Terminal ||
        prod instanceof Rule);
}
function isOptionalProd(prod, alreadyVisited = []) {
    const isDirectlyOptional = prod instanceof Option ||
        prod instanceof Repetition ||
        prod instanceof RepetitionWithSeparator;
    if (isDirectlyOptional) {
        return true;
    }
    // note that this can cause infinite loop if one optional empty TOP production has a cyclic dependency with another
    // empty optional top rule
    // may be indirectly optional ((A?B?C?) | (D?E?F?))
    if (prod instanceof Alternation) {
        // for OR its enough for just one of the alternatives to be optional
        return lodash_es_some(prod.definition, (subProd) => {
            return isOptionalProd(subProd, alreadyVisited);
        });
    }
    else if (prod instanceof NonTerminal && lodash_es_includes(alreadyVisited, prod)) {
        // avoiding stack overflow due to infinite recursion
        return false;
    }
    else if (prod instanceof AbstractProduction) {
        if (prod instanceof NonTerminal) {
            alreadyVisited.push(prod);
        }
        return lodash_es_every(prod.definition, (subProd) => {
            return isOptionalProd(subProd, alreadyVisited);
        });
    }
    else {
        return false;
    }
}
function isBranchingProd(prod) {
    return prod instanceof Alternation;
}
function getProductionDslName(prod) {
    /* istanbul ignore else */
    if (prod instanceof NonTerminal) {
        return "SUBRULE";
    }
    else if (prod instanceof Option) {
        return "OPTION";
    }
    else if (prod instanceof Alternation) {
        return "OR";
    }
    else if (prod instanceof RepetitionMandatory) {
        return "AT_LEAST_ONE";
    }
    else if (prod instanceof RepetitionMandatoryWithSeparator) {
        return "AT_LEAST_ONE_SEP";
    }
    else if (prod instanceof RepetitionWithSeparator) {
        return "MANY_SEP";
    }
    else if (prod instanceof Repetition) {
        return "MANY";
    }
    else if (prod instanceof Terminal) {
        return "CONSUME";
        /* c8 ignore next 3 */
    }
    else {
        throw Error("non exhaustive match");
    }
}
//# sourceMappingURL=helpers.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/first.js


function first(prod) {
    /* istanbul ignore else */
    if (prod instanceof NonTerminal) {
        // this could in theory cause infinite loops if
        // (1) prod A refs prod B.
        // (2) prod B refs prod A
        // (3) AB can match the empty set
        // in other words a cycle where everything is optional so the first will keep
        // looking ahead for the next optional part and will never exit
        // currently there is no safeguard for this unique edge case because
        // (1) not sure a grammar in which this can happen is useful for anything (productive)
        return first(prod.referencedRule);
    }
    else if (prod instanceof Terminal) {
        return firstForTerminal(prod);
    }
    else if (isSequenceProd(prod)) {
        return firstForSequence(prod);
    }
    else if (isBranchingProd(prod)) {
        return firstForBranching(prod);
    }
    else {
        throw Error("non exhaustive match");
    }
}
function firstForSequence(prod) {
    let firstSet = [];
    const seq = prod.definition;
    let nextSubProdIdx = 0;
    let hasInnerProdsRemaining = seq.length > nextSubProdIdx;
    let currSubProd;
    // so we enter the loop at least once (if the definition is not empty
    let isLastInnerProdOptional = true;
    // scan a sequence until it's end or until we have found a NONE optional production in it
    while (hasInnerProdsRemaining && isLastInnerProdOptional) {
        currSubProd = seq[nextSubProdIdx];
        isLastInnerProdOptional = isOptionalProd(currSubProd);
        firstSet = firstSet.concat(first(currSubProd));
        nextSubProdIdx = nextSubProdIdx + 1;
        hasInnerProdsRemaining = seq.length > nextSubProdIdx;
    }
    return lodash_es_uniq(firstSet);
}
function firstForBranching(prod) {
    const allAlternativesFirsts = (0,map/* default */.Z)(prod.definition, (innerProd) => {
        return first(innerProd);
    });
    return lodash_es_uniq((0,flatten/* default */.Z)(allAlternativesFirsts));
}
function firstForTerminal(terminal) {
    return [terminal.terminalType];
}
//# sourceMappingURL=first.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/constants.js
// TODO: can this be removed? where is it used?
const constants_IN = "_~IN~_";
//# sourceMappingURL=constants.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/follow.js





// This ResyncFollowsWalker computes all of the follows required for RESYNC
// (skipping reference production).
class ResyncFollowsWalker extends RestWalker {
    constructor(topProd) {
        super();
        this.topProd = topProd;
        this.follows = {};
    }
    startWalking() {
        this.walk(this.topProd);
        return this.follows;
    }
    walkTerminal(terminal, currRest, prevRest) {
        // do nothing! just like in the public sector after 13:00
    }
    walkProdRef(refProd, currRest, prevRest) {
        const followName = buildBetweenProdsFollowPrefix(refProd.referencedRule, refProd.idx) +
            this.topProd.name;
        const fullRest = currRest.concat(prevRest);
        const restProd = new Alternative({ definition: fullRest });
        const t_in_topProd_follows = first(restProd);
        this.follows[followName] = t_in_topProd_follows;
    }
}
function computeAllProdsFollows(topProductions) {
    const reSyncFollows = {};
    (0,forEach/* default */.Z)(topProductions, (topProd) => {
        const currRefsFollow = new ResyncFollowsWalker(topProd).startWalking();
        lodash_es_assign(reSyncFollows, currRefsFollow);
    });
    return reSyncFollows;
}
function buildBetweenProdsFollowPrefix(inner, occurenceInParent) {
    return inner.name + occurenceInParent + constants_IN;
}
function buildInProdFollowPrefix(terminal) {
    const terminalName = terminal.terminalType.name;
    return terminalName + terminal.idx + IN;
}
//# sourceMappingURL=follow.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/isUndefined.js
var isUndefined = __webpack_require__(52307);
// EXTERNAL MODULE: ../node_modules/@chevrotain/regexp-to-ast/lib/src/api.js + 4 modules
var api = __webpack_require__(77647);
// EXTERNAL MODULE: ../node_modules/lodash-es/defaults.js
var defaults = __webpack_require__(65479);
// EXTERNAL MODULE: ../node_modules/lodash-es/_arrayFilter.js
var _arrayFilter = __webpack_require__(11819);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseFilter.js
var _baseFilter = __webpack_require__(45701);
;// CONCATENATED MODULE: ../node_modules/lodash-es/negate.js
/** Error message constants. */
var FUNC_ERROR_TEXT = 'Expected a function';

/**
 * Creates a function that negates the result of the predicate `func`. The
 * `func` predicate is invoked with the `this` binding and arguments of the
 * created function.
 *
 * @static
 * @memberOf _
 * @since 3.0.0
 * @category Function
 * @param {Function} predicate The predicate to negate.
 * @returns {Function} Returns the new negated function.
 * @example
 *
 * function isEven(n) {
 *   return n % 2 == 0;
 * }
 *
 * _.filter([1, 2, 3, 4, 5, 6], _.negate(isEven));
 * // => [1, 3, 5]
 */
function negate(predicate) {
  if (typeof predicate != 'function') {
    throw new TypeError(FUNC_ERROR_TEXT);
  }
  return function() {
    var args = arguments;
    switch (args.length) {
      case 0: return !predicate.call(this);
      case 1: return !predicate.call(this, args[0]);
      case 2: return !predicate.call(this, args[0], args[1]);
      case 3: return !predicate.call(this, args[0], args[1], args[2]);
    }
    return !predicate.apply(this, args);
  };
}

/* harmony default export */ const lodash_es_negate = (negate);

;// CONCATENATED MODULE: ../node_modules/lodash-es/reject.js






/**
 * The opposite of `_.filter`; this method returns the elements of `collection`
 * that `predicate` does **not** return truthy for.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [predicate=_.identity] The function invoked per iteration.
 * @returns {Array} Returns the new filtered array.
 * @see _.filter
 * @example
 *
 * var users = [
 *   { 'user': 'barney', 'age': 36, 'active': false },
 *   { 'user': 'fred',   'age': 40, 'active': true }
 * ];
 *
 * _.reject(users, function(o) { return !o.active; });
 * // => objects for ['fred']
 *
 * // The `_.matches` iteratee shorthand.
 * _.reject(users, { 'age': 40, 'active': true });
 * // => objects for ['barney']
 *
 * // The `_.matchesProperty` iteratee shorthand.
 * _.reject(users, ['active', false]);
 * // => objects for ['fred']
 *
 * // The `_.property` iteratee shorthand.
 * _.reject(users, 'active');
 * // => objects for ['barney']
 */
function reject(collection, predicate) {
  var func = (0,isArray/* default */.Z)(collection) ? _arrayFilter/* default */.Z : _baseFilter/* default */.Z;
  return func(collection, lodash_es_negate((0,_baseIteratee/* default */.Z)(predicate, 3)));
}

/* harmony default export */ const lodash_es_reject = (reject);

// EXTERNAL MODULE: ../node_modules/lodash-es/isFunction.js
var isFunction = __webpack_require__(48489);
;// CONCATENATED MODULE: ../node_modules/lodash-es/indexOf.js



/* Built-in method references for those with the same name as other `lodash` methods. */
var indexOf_nativeMax = Math.max;

/**
 * Gets the index at which the first occurrence of `value` is found in `array`
 * using [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * for equality comparisons. If `fromIndex` is negative, it's used as the
 * offset from the end of `array`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to inspect.
 * @param {*} value The value to search for.
 * @param {number} [fromIndex=0] The index to search from.
 * @returns {number} Returns the index of the matched value, else `-1`.
 * @example
 *
 * _.indexOf([1, 2, 1, 2], 2);
 * // => 1
 *
 * // Search from the `fromIndex`.
 * _.indexOf([1, 2, 1, 2], 2, 2);
 * // => 3
 */
function indexOf(array, value, fromIndex) {
  var length = array == null ? 0 : array.length;
  if (!length) {
    return -1;
  }
  var index = fromIndex == null ? 0 : (0,toInteger/* default */.Z)(fromIndex);
  if (index < 0) {
    index = indexOf_nativeMax(length + index, 0);
  }
  return (0,_baseIndexOf/* default */.Z)(array, value, index);
}

/* harmony default export */ const lodash_es_indexOf = (indexOf);

// EXTERNAL MODULE: ../node_modules/lodash-es/reduce.js + 2 modules
var reduce = __webpack_require__(99413);
// EXTERNAL MODULE: ../node_modules/lodash-es/filter.js
var filter = __webpack_require__(11382);
// EXTERNAL MODULE: ../node_modules/lodash-es/_SetCache.js + 2 modules
var _SetCache = __webpack_require__(40105);
// EXTERNAL MODULE: ../node_modules/lodash-es/_arrayIncludes.js
var _arrayIncludes = __webpack_require__(79036);
// EXTERNAL MODULE: ../node_modules/lodash-es/_arrayIncludesWith.js
var _arrayIncludesWith = __webpack_require__(58084);
// EXTERNAL MODULE: ../node_modules/lodash-es/_cacheHas.js
var _cacheHas = __webpack_require__(8142);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseDifference.js







/** Used as the size to enable large array optimizations. */
var LARGE_ARRAY_SIZE = 200;

/**
 * The base implementation of methods like `_.difference` without support
 * for excluding multiple arrays or iteratee shorthands.
 *
 * @private
 * @param {Array} array The array to inspect.
 * @param {Array} values The values to exclude.
 * @param {Function} [iteratee] The iteratee invoked per element.
 * @param {Function} [comparator] The comparator invoked per element.
 * @returns {Array} Returns the new array of filtered values.
 */
function baseDifference(array, values, iteratee, comparator) {
  var index = -1,
      includes = _arrayIncludes/* default */.Z,
      isCommon = true,
      length = array.length,
      result = [],
      valuesLength = values.length;

  if (!length) {
    return result;
  }
  if (iteratee) {
    values = (0,_arrayMap/* default */.Z)(values, (0,_baseUnary/* default */.Z)(iteratee));
  }
  if (comparator) {
    includes = _arrayIncludesWith/* default */.Z;
    isCommon = false;
  }
  else if (values.length >= LARGE_ARRAY_SIZE) {
    includes = _cacheHas/* default */.Z;
    isCommon = false;
    values = new _SetCache/* default */.Z(values);
  }
  outer:
  while (++index < length) {
    var value = array[index],
        computed = iteratee == null ? value : iteratee(value);

    value = (comparator || value !== 0) ? value : 0;
    if (isCommon && computed === computed) {
      var valuesIndex = valuesLength;
      while (valuesIndex--) {
        if (values[valuesIndex] === computed) {
          continue outer;
        }
      }
      result.push(value);
    }
    else if (!includes(values, computed, comparator)) {
      result.push(value);
    }
  }
  return result;
}

/* harmony default export */ const _baseDifference = (baseDifference);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseFlatten.js + 1 modules
var _baseFlatten = __webpack_require__(65029);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseRest.js
var _baseRest = __webpack_require__(99719);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArrayLikeObject.js
var isArrayLikeObject = __webpack_require__(60492);
;// CONCATENATED MODULE: ../node_modules/lodash-es/difference.js





/**
 * Creates an array of `array` values not included in the other given arrays
 * using [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * for equality comparisons. The order and references of result values are
 * determined by the first array.
 *
 * **Note:** Unlike `_.pullAll`, this method returns a new array.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to inspect.
 * @param {...Array} [values] The values to exclude.
 * @returns {Array} Returns the new array of filtered values.
 * @see _.without, _.xor
 * @example
 *
 * _.difference([2, 1], [2, 3]);
 * // => [1]
 */
var difference = (0,_baseRest/* default */.Z)(function(array, values) {
  return (0,isArrayLikeObject/* default */.Z)(array)
    ? _baseDifference(array, (0,_baseFlatten/* default */.Z)(values, 1, isArrayLikeObject/* default */.Z, true))
    : [];
});

/* harmony default export */ const lodash_es_difference = (difference);

;// CONCATENATED MODULE: ../node_modules/lodash-es/compact.js
/**
 * Creates an array with all falsey values removed. The values `false`, `null`,
 * `0`, `""`, `undefined`, and `NaN` are falsey.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to compact.
 * @returns {Array} Returns the new array of filtered values.
 * @example
 *
 * _.compact([0, 1, false, 2, '', 3]);
 * // => [1, 2, 3]
 */
function compact(array) {
  var index = -1,
      length = array == null ? 0 : array.length,
      resIndex = 0,
      result = [];

  while (++index < length) {
    var value = array[index];
    if (value) {
      result[resIndex++] = value;
    }
  }
  return result;
}

/* harmony default export */ const lodash_es_compact = (compact);

;// CONCATENATED MODULE: ../node_modules/lodash-es/head.js
/**
 * Gets the first element of `array`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @alias first
 * @category Array
 * @param {Array} array The array to query.
 * @returns {*} Returns the first element of `array`.
 * @example
 *
 * _.head([1, 2, 3]);
 * // => 1
 *
 * _.head([]);
 * // => undefined
 */
function head(array) {
  return (array && array.length) ? array[0] : undefined;
}

/* harmony default export */ const lodash_es_head = (head);

// EXTERNAL MODULE: ../node_modules/lodash-es/find.js + 2 modules
var find = __webpack_require__(90970);
;// CONCATENATED MODULE: ../node_modules/@chevrotain/utils/lib/src/print.js
function PRINT_ERROR(msg) {
    /* istanbul ignore else - can't override global.console in node.js */
    if (console && console.error) {
        console.error(`Error: ${msg}`);
    }
}
function PRINT_WARNING(msg) {
    /* istanbul ignore else - can't override global.console in node.js*/
    if (console && console.warn) {
        // TODO: modify docs accordingly
        console.warn(`Warning: ${msg}`);
    }
}
//# sourceMappingURL=print.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/reg_exp_parser.js

let regExpAstCache = {};
const regExpParser = new api/* RegExpParser */.O();
function getRegExpAst(regExp) {
    const regExpStr = regExp.toString();
    if (regExpAstCache.hasOwnProperty(regExpStr)) {
        return regExpAstCache[regExpStr];
    }
    else {
        const regExpAst = regExpParser.pattern(regExpStr);
        regExpAstCache[regExpStr] = regExpAst;
        return regExpAst;
    }
}
function clearRegExpParserCache() {
    regExpAstCache = {};
}
//# sourceMappingURL=reg_exp_parser.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/reg_exp.js





const complementErrorMessage = "Complement Sets are not supported for first char optimization";
const failedOptimizationPrefixMsg = 'Unable to use "first char" lexer optimizations:\n';
function getOptimizedStartCodesIndices(regExp, ensureOptimizations = false) {
    try {
        const ast = getRegExpAst(regExp);
        const firstChars = firstCharOptimizedIndices(ast.value, {}, ast.flags.ignoreCase);
        return firstChars;
    }
    catch (e) {
        /* istanbul ignore next */
        // Testing this relies on the regexp-to-ast library having a bug... */
        // TODO: only the else branch needs to be ignored, try to fix with newer prettier / tsc
        if (e.message === complementErrorMessage) {
            if (ensureOptimizations) {
                PRINT_WARNING(`${failedOptimizationPrefixMsg}` +
                    `\tUnable to optimize: < ${regExp.toString()} >\n` +
                    "\tComplement Sets cannot be automatically optimized.\n" +
                    "\tThis will disable the lexer's first char optimizations.\n" +
                    "\tSee: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#COMPLEMENT for details.");
            }
        }
        else {
            let msgSuffix = "";
            if (ensureOptimizations) {
                msgSuffix =
                    "\n\tThis will disable the lexer's first char optimizations.\n" +
                        "\tSee: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#REGEXP_PARSING for details.";
            }
            PRINT_ERROR(`${failedOptimizationPrefixMsg}\n` +
                `\tFailed parsing: < ${regExp.toString()} >\n` +
                `\tUsing the @chevrotain/regexp-to-ast library\n` +
                "\tPlease open an issue at: https://github.com/chevrotain/chevrotain/issues" +
                msgSuffix);
        }
    }
    return [];
}
function firstCharOptimizedIndices(ast, result, ignoreCase) {
    switch (ast.type) {
        case "Disjunction":
            for (let i = 0; i < ast.value.length; i++) {
                firstCharOptimizedIndices(ast.value[i], result, ignoreCase);
            }
            break;
        case "Alternative":
            const terms = ast.value;
            for (let i = 0; i < terms.length; i++) {
                const term = terms[i];
                // skip terms that cannot effect the first char results
                switch (term.type) {
                    case "EndAnchor":
                    // A group back reference cannot affect potential starting char.
                    // because if a back reference is the first production than automatically
                    // the group being referenced has had to come BEFORE so its codes have already been added
                    case "GroupBackReference":
                    // assertions do not affect potential starting codes
                    case "Lookahead":
                    case "NegativeLookahead":
                    case "StartAnchor":
                    case "WordBoundary":
                    case "NonWordBoundary":
                        continue;
                }
                const atom = term;
                switch (atom.type) {
                    case "Character":
                        addOptimizedIdxToResult(atom.value, result, ignoreCase);
                        break;
                    case "Set":
                        if (atom.complement === true) {
                            throw Error(complementErrorMessage);
                        }
                        (0,forEach/* default */.Z)(atom.value, (code) => {
                            if (typeof code === "number") {
                                addOptimizedIdxToResult(code, result, ignoreCase);
                            }
                            else {
                                // range
                                const range = code;
                                // cannot optimize when ignoreCase is
                                if (ignoreCase === true) {
                                    for (let rangeCode = range.from; rangeCode <= range.to; rangeCode++) {
                                        addOptimizedIdxToResult(rangeCode, result, ignoreCase);
                                    }
                                }
                                // Optimization (2 orders of magnitude less work for very large ranges)
                                else {
                                    // handle unoptimized values
                                    for (let rangeCode = range.from; rangeCode <= range.to && rangeCode < minOptimizationVal; rangeCode++) {
                                        addOptimizedIdxToResult(rangeCode, result, ignoreCase);
                                    }
                                    // Less common charCode where we optimize for faster init time, by using larger "buckets"
                                    if (range.to >= minOptimizationVal) {
                                        const minUnOptVal = range.from >= minOptimizationVal
                                            ? range.from
                                            : minOptimizationVal;
                                        const maxUnOptVal = range.to;
                                        const minOptIdx = charCodeToOptimizedIndex(minUnOptVal);
                                        const maxOptIdx = charCodeToOptimizedIndex(maxUnOptVal);
                                        for (let currOptIdx = minOptIdx; currOptIdx <= maxOptIdx; currOptIdx++) {
                                            result[currOptIdx] = currOptIdx;
                                        }
                                    }
                                }
                            }
                        });
                        break;
                    case "Group":
                        firstCharOptimizedIndices(atom.value, result, ignoreCase);
                        break;
                    /* istanbul ignore next */
                    default:
                        throw Error("Non Exhaustive Match");
                }
                // reached a mandatory production, no more **start** codes can be found on this alternative
                const isOptionalQuantifier = atom.quantifier !== undefined && atom.quantifier.atLeast === 0;
                if (
                // A group may be optional due to empty contents /(?:)/
                // or if everything inside it is optional /((a)?)/
                (atom.type === "Group" && isWholeOptional(atom) === false) ||
                    // If this term is not a group it may only be optional if it has an optional quantifier
                    (atom.type !== "Group" && isOptionalQuantifier === false)) {
                    break;
                }
            }
            break;
        /* istanbul ignore next */
        default:
            throw Error("non exhaustive match!");
    }
    // console.log(Object.keys(result).length)
    return (0,values/* default */.Z)(result);
}
function addOptimizedIdxToResult(code, result, ignoreCase) {
    const optimizedCharIdx = charCodeToOptimizedIndex(code);
    result[optimizedCharIdx] = optimizedCharIdx;
    if (ignoreCase === true) {
        handleIgnoreCase(code, result);
    }
}
function handleIgnoreCase(code, result) {
    const char = String.fromCharCode(code);
    const upperChar = char.toUpperCase();
    /* istanbul ignore else */
    if (upperChar !== char) {
        const optimizedCharIdx = charCodeToOptimizedIndex(upperChar.charCodeAt(0));
        result[optimizedCharIdx] = optimizedCharIdx;
    }
    else {
        const lowerChar = char.toLowerCase();
        if (lowerChar !== char) {
            const optimizedCharIdx = charCodeToOptimizedIndex(lowerChar.charCodeAt(0));
            result[optimizedCharIdx] = optimizedCharIdx;
        }
    }
}
function findCode(setNode, targetCharCodes) {
    return (0,find/* default */.Z)(setNode.value, (codeOrRange) => {
        if (typeof codeOrRange === "number") {
            return lodash_es_includes(targetCharCodes, codeOrRange);
        }
        else {
            // range
            const range = codeOrRange;
            return ((0,find/* default */.Z)(targetCharCodes, (targetCode) => range.from <= targetCode && targetCode <= range.to) !== undefined);
        }
    });
}
function isWholeOptional(ast) {
    const quantifier = ast.quantifier;
    if (quantifier && quantifier.atLeast === 0) {
        return true;
    }
    if (!ast.value) {
        return false;
    }
    return (0,isArray/* default */.Z)(ast.value)
        ? lodash_es_every(ast.value, isWholeOptional)
        : isWholeOptional(ast.value);
}
class CharCodeFinder extends api/* BaseRegExpVisitor */.e {
    constructor(targetCharCodes) {
        super();
        this.targetCharCodes = targetCharCodes;
        this.found = false;
    }
    visitChildren(node) {
        // No need to keep looking...
        if (this.found === true) {
            return;
        }
        // switch lookaheads as they do not actually consume any characters thus
        // finding a charCode at lookahead context does not mean that regexp can actually contain it in a match.
        switch (node.type) {
            case "Lookahead":
                this.visitLookahead(node);
                return;
            case "NegativeLookahead":
                this.visitNegativeLookahead(node);
                return;
        }
        super.visitChildren(node);
    }
    visitCharacter(node) {
        if (lodash_es_includes(this.targetCharCodes, node.value)) {
            this.found = true;
        }
    }
    visitSet(node) {
        if (node.complement) {
            if (findCode(node, this.targetCharCodes) === undefined) {
                this.found = true;
            }
        }
        else {
            if (findCode(node, this.targetCharCodes) !== undefined) {
                this.found = true;
            }
        }
    }
}
function canMatchCharCode(charCodes, pattern) {
    if (pattern instanceof RegExp) {
        const ast = getRegExpAst(pattern);
        const charCodeFinder = new CharCodeFinder(charCodes);
        charCodeFinder.visit(ast);
        return charCodeFinder.found;
    }
    else {
        return ((0,find/* default */.Z)(pattern, (char) => {
            return lodash_es_includes(charCodes, char.charCodeAt(0));
        }) !== undefined);
    }
}
//# sourceMappingURL=reg_exp.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/lexer.js






const PATTERN = "PATTERN";
const DEFAULT_MODE = "defaultMode";
const MODES = "modes";
let SUPPORT_STICKY = typeof new RegExp("(?:)").sticky === "boolean";
function disableSticky() {
    SUPPORT_STICKY = false;
}
function enableSticky() {
    SUPPORT_STICKY = true;
}
function analyzeTokenTypes(tokenTypes, options) {
    options = (0,defaults/* default */.Z)(options, {
        useSticky: SUPPORT_STICKY,
        debug: false,
        safeMode: false,
        positionTracking: "full",
        lineTerminatorCharacters: ["\r", "\n"],
        tracer: (msg, action) => action(),
    });
    const tracer = options.tracer;
    tracer("initCharCodeToOptimizedIndexMap", () => {
        initCharCodeToOptimizedIndexMap();
    });
    let onlyRelevantTypes;
    tracer("Reject Lexer.NA", () => {
        onlyRelevantTypes = lodash_es_reject(tokenTypes, (currType) => {
            return currType[PATTERN] === Lexer.NA;
        });
    });
    let hasCustom = false;
    let allTransformedPatterns;
    tracer("Transform Patterns", () => {
        hasCustom = false;
        allTransformedPatterns = (0,map/* default */.Z)(onlyRelevantTypes, (currType) => {
            const currPattern = currType[PATTERN];
            /* istanbul ignore else */
            if (lodash_es_isRegExp(currPattern)) {
                const regExpSource = currPattern.source;
                if (regExpSource.length === 1 &&
                    // only these regExp meta characters which can appear in a length one regExp
                    regExpSource !== "^" &&
                    regExpSource !== "$" &&
                    regExpSource !== "." &&
                    !currPattern.ignoreCase) {
                    return regExpSource;
                }
                else if (regExpSource.length === 2 &&
                    regExpSource[0] === "\\" &&
                    // not a meta character
                    !lodash_es_includes([
                        "d",
                        "D",
                        "s",
                        "S",
                        "t",
                        "r",
                        "n",
                        "t",
                        "0",
                        "c",
                        "b",
                        "B",
                        "f",
                        "v",
                        "w",
                        "W",
                    ], regExpSource[1])) {
                    // escaped meta Characters: /\+/ /\[/
                    // or redundant escaping: /\a/
                    // without the escaping "\"
                    return regExpSource[1];
                }
                else {
                    return options.useSticky
                        ? addStickyFlag(currPattern)
                        : addStartOfInput(currPattern);
                }
            }
            else if ((0,isFunction/* default */.Z)(currPattern)) {
                hasCustom = true;
                // CustomPatternMatcherFunc - custom patterns do not require any transformations, only wrapping in a RegExp Like object
                return { exec: currPattern };
            }
            else if (typeof currPattern === "object") {
                hasCustom = true;
                // ICustomPattern
                return currPattern;
            }
            else if (typeof currPattern === "string") {
                if (currPattern.length === 1) {
                    return currPattern;
                }
                else {
                    const escapedRegExpString = currPattern.replace(/[\\^$.*+?()[\]{}|]/g, "\\$&");
                    const wrappedRegExp = new RegExp(escapedRegExpString);
                    return options.useSticky
                        ? addStickyFlag(wrappedRegExp)
                        : addStartOfInput(wrappedRegExp);
                }
            }
            else {
                throw Error("non exhaustive match");
            }
        });
    });
    let patternIdxToType;
    let patternIdxToGroup;
    let patternIdxToLongerAltIdxArr;
    let patternIdxToPushMode;
    let patternIdxToPopMode;
    tracer("misc mapping", () => {
        patternIdxToType = (0,map/* default */.Z)(onlyRelevantTypes, (currType) => currType.tokenTypeIdx);
        patternIdxToGroup = (0,map/* default */.Z)(onlyRelevantTypes, (clazz) => {
            const groupName = clazz.GROUP;
            /* istanbul ignore next */
            if (groupName === Lexer.SKIPPED) {
                return undefined;
            }
            else if ((0,isString/* default */.Z)(groupName)) {
                return groupName;
            }
            else if ((0,isUndefined/* default */.Z)(groupName)) {
                return false;
            }
            else {
                throw Error("non exhaustive match");
            }
        });
        patternIdxToLongerAltIdxArr = (0,map/* default */.Z)(onlyRelevantTypes, (clazz) => {
            const longerAltType = clazz.LONGER_ALT;
            if (longerAltType) {
                const longerAltIdxArr = (0,isArray/* default */.Z)(longerAltType)
                    ? (0,map/* default */.Z)(longerAltType, (type) => lodash_es_indexOf(onlyRelevantTypes, type))
                    : [lodash_es_indexOf(onlyRelevantTypes, longerAltType)];
                return longerAltIdxArr;
            }
        });
        patternIdxToPushMode = (0,map/* default */.Z)(onlyRelevantTypes, (clazz) => clazz.PUSH_MODE);
        patternIdxToPopMode = (0,map/* default */.Z)(onlyRelevantTypes, (clazz) => (0,has/* default */.Z)(clazz, "POP_MODE"));
    });
    let patternIdxToCanLineTerminator;
    tracer("Line Terminator Handling", () => {
        const lineTerminatorCharCodes = getCharCodes(options.lineTerminatorCharacters);
        patternIdxToCanLineTerminator = (0,map/* default */.Z)(onlyRelevantTypes, (tokType) => false);
        if (options.positionTracking !== "onlyOffset") {
            patternIdxToCanLineTerminator = (0,map/* default */.Z)(onlyRelevantTypes, (tokType) => {
                if ((0,has/* default */.Z)(tokType, "LINE_BREAKS")) {
                    return !!tokType.LINE_BREAKS;
                }
                else {
                    return (checkLineBreaksIssues(tokType, lineTerminatorCharCodes) === false &&
                        canMatchCharCode(lineTerminatorCharCodes, tokType.PATTERN));
                }
            });
        }
    });
    let patternIdxToIsCustom;
    let patternIdxToShort;
    let emptyGroups;
    let patternIdxToConfig;
    tracer("Misc Mapping #2", () => {
        patternIdxToIsCustom = (0,map/* default */.Z)(onlyRelevantTypes, isCustomPattern);
        patternIdxToShort = (0,map/* default */.Z)(allTransformedPatterns, isShortPattern);
        emptyGroups = (0,reduce/* default */.Z)(onlyRelevantTypes, (acc, clazz) => {
            const groupName = clazz.GROUP;
            if ((0,isString/* default */.Z)(groupName) && !(groupName === Lexer.SKIPPED)) {
                acc[groupName] = [];
            }
            return acc;
        }, {});
        patternIdxToConfig = (0,map/* default */.Z)(allTransformedPatterns, (x, idx) => {
            return {
                pattern: allTransformedPatterns[idx],
                longerAlt: patternIdxToLongerAltIdxArr[idx],
                canLineTerminator: patternIdxToCanLineTerminator[idx],
                isCustom: patternIdxToIsCustom[idx],
                short: patternIdxToShort[idx],
                group: patternIdxToGroup[idx],
                push: patternIdxToPushMode[idx],
                pop: patternIdxToPopMode[idx],
                tokenTypeIdx: patternIdxToType[idx],
                tokenType: onlyRelevantTypes[idx],
            };
        });
    });
    let canBeOptimized = true;
    let charCodeToPatternIdxToConfig = [];
    if (!options.safeMode) {
        tracer("First Char Optimization", () => {
            charCodeToPatternIdxToConfig = (0,reduce/* default */.Z)(onlyRelevantTypes, (result, currTokType, idx) => {
                if (typeof currTokType.PATTERN === "string") {
                    const charCode = currTokType.PATTERN.charCodeAt(0);
                    const optimizedIdx = charCodeToOptimizedIndex(charCode);
                    addToMapOfArrays(result, optimizedIdx, patternIdxToConfig[idx]);
                }
                else if ((0,isArray/* default */.Z)(currTokType.START_CHARS_HINT)) {
                    let lastOptimizedIdx;
                    (0,forEach/* default */.Z)(currTokType.START_CHARS_HINT, (charOrInt) => {
                        const charCode = typeof charOrInt === "string"
                            ? charOrInt.charCodeAt(0)
                            : charOrInt;
                        const currOptimizedIdx = charCodeToOptimizedIndex(charCode);
                        // Avoid adding the config multiple times
                        /* istanbul ignore else */
                        // - Difficult to check this scenario effects as it is only a performance
                        //   optimization that does not change correctness
                        if (lastOptimizedIdx !== currOptimizedIdx) {
                            lastOptimizedIdx = currOptimizedIdx;
                            addToMapOfArrays(result, currOptimizedIdx, patternIdxToConfig[idx]);
                        }
                    });
                }
                else if (lodash_es_isRegExp(currTokType.PATTERN)) {
                    if (currTokType.PATTERN.unicode) {
                        canBeOptimized = false;
                        if (options.ensureOptimizations) {
                            PRINT_ERROR(`${failedOptimizationPrefixMsg}` +
                                `\tUnable to analyze < ${currTokType.PATTERN.toString()} > pattern.\n` +
                                "\tThe regexp unicode flag is not currently supported by the regexp-to-ast library.\n" +
                                "\tThis will disable the lexer's first char optimizations.\n" +
                                "\tFor details See: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#UNICODE_OPTIMIZE");
                        }
                    }
                    else {
                        const optimizedCodes = getOptimizedStartCodesIndices(currTokType.PATTERN, options.ensureOptimizations);
                        /* istanbul ignore if */
                        // start code will only be empty given an empty regExp or failure of regexp-to-ast library
                        // the first should be a different validation and the second cannot be tested.
                        if ((0,isEmpty/* default */.Z)(optimizedCodes)) {
                            // we cannot understand what codes may start possible matches
                            // The optimization correctness requires knowing start codes for ALL patterns.
                            // Not actually sure this is an error, no debug message
                            canBeOptimized = false;
                        }
                        (0,forEach/* default */.Z)(optimizedCodes, (code) => {
                            addToMapOfArrays(result, code, patternIdxToConfig[idx]);
                        });
                    }
                }
                else {
                    if (options.ensureOptimizations) {
                        PRINT_ERROR(`${failedOptimizationPrefixMsg}` +
                            `\tTokenType: <${currTokType.name}> is using a custom token pattern without providing <start_chars_hint> parameter.\n` +
                            "\tThis will disable the lexer's first char optimizations.\n" +
                            "\tFor details See: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#CUSTOM_OPTIMIZE");
                    }
                    canBeOptimized = false;
                }
                return result;
            }, []);
        });
    }
    return {
        emptyGroups: emptyGroups,
        patternIdxToConfig: patternIdxToConfig,
        charCodeToPatternIdxToConfig: charCodeToPatternIdxToConfig,
        hasCustom: hasCustom,
        canBeOptimized: canBeOptimized,
    };
}
function validatePatterns(tokenTypes, validModesNames) {
    let errors = [];
    const missingResult = findMissingPatterns(tokenTypes);
    errors = errors.concat(missingResult.errors);
    const invalidResult = findInvalidPatterns(missingResult.valid);
    const validTokenTypes = invalidResult.valid;
    errors = errors.concat(invalidResult.errors);
    errors = errors.concat(validateRegExpPattern(validTokenTypes));
    errors = errors.concat(findInvalidGroupType(validTokenTypes));
    errors = errors.concat(findModesThatDoNotExist(validTokenTypes, validModesNames));
    errors = errors.concat(findUnreachablePatterns(validTokenTypes));
    return errors;
}
function validateRegExpPattern(tokenTypes) {
    let errors = [];
    const withRegExpPatterns = (0,filter/* default */.Z)(tokenTypes, (currTokType) => lodash_es_isRegExp(currTokType[PATTERN]));
    errors = errors.concat(findEndOfInputAnchor(withRegExpPatterns));
    errors = errors.concat(findStartOfInputAnchor(withRegExpPatterns));
    errors = errors.concat(findUnsupportedFlags(withRegExpPatterns));
    errors = errors.concat(findDuplicatePatterns(withRegExpPatterns));
    errors = errors.concat(findEmptyMatchRegExps(withRegExpPatterns));
    return errors;
}
function findMissingPatterns(tokenTypes) {
    const tokenTypesWithMissingPattern = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        return !(0,has/* default */.Z)(currType, PATTERN);
    });
    const errors = (0,map/* default */.Z)(tokenTypesWithMissingPattern, (currType) => {
        return {
            message: "Token Type: ->" +
                currType.name +
                "<- missing static 'PATTERN' property",
            type: LexerDefinitionErrorType.MISSING_PATTERN,
            tokenTypes: [currType],
        };
    });
    const valid = lodash_es_difference(tokenTypes, tokenTypesWithMissingPattern);
    return { errors, valid };
}
function findInvalidPatterns(tokenTypes) {
    const tokenTypesWithInvalidPattern = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        const pattern = currType[PATTERN];
        return (!lodash_es_isRegExp(pattern) &&
            !(0,isFunction/* default */.Z)(pattern) &&
            !(0,has/* default */.Z)(pattern, "exec") &&
            !(0,isString/* default */.Z)(pattern));
    });
    const errors = (0,map/* default */.Z)(tokenTypesWithInvalidPattern, (currType) => {
        return {
            message: "Token Type: ->" +
                currType.name +
                "<- static 'PATTERN' can only be a RegExp, a" +
                " Function matching the {CustomPatternMatcherFunc} type or an Object matching the {ICustomPattern} interface.",
            type: LexerDefinitionErrorType.INVALID_PATTERN,
            tokenTypes: [currType],
        };
    });
    const valid = lodash_es_difference(tokenTypes, tokenTypesWithInvalidPattern);
    return { errors, valid };
}
const end_of_input = /[^\\][$]/;
function findEndOfInputAnchor(tokenTypes) {
    class EndAnchorFinder extends api/* BaseRegExpVisitor */.e {
        constructor() {
            super(...arguments);
            this.found = false;
        }
        visitEndAnchor(node) {
            this.found = true;
        }
    }
    const invalidRegex = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        const pattern = currType.PATTERN;
        try {
            const regexpAst = getRegExpAst(pattern);
            const endAnchorVisitor = new EndAnchorFinder();
            endAnchorVisitor.visit(regexpAst);
            return endAnchorVisitor.found;
        }
        catch (e) {
            // old behavior in case of runtime exceptions with regexp-to-ast.
            /* istanbul ignore next - cannot ensure an error in regexp-to-ast*/
            return end_of_input.test(pattern.source);
        }
    });
    const errors = (0,map/* default */.Z)(invalidRegex, (currType) => {
        return {
            message: "Unexpected RegExp Anchor Error:\n" +
                "\tToken Type: ->" +
                currType.name +
                "<- static 'PATTERN' cannot contain end of input anchor '$'\n" +
                "\tSee chevrotain.io/docs/guide/resolving_lexer_errors.html#ANCHORS" +
                "\tfor details.",
            type: LexerDefinitionErrorType.EOI_ANCHOR_FOUND,
            tokenTypes: [currType],
        };
    });
    return errors;
}
function findEmptyMatchRegExps(tokenTypes) {
    const matchesEmptyString = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        const pattern = currType.PATTERN;
        return pattern.test("");
    });
    const errors = (0,map/* default */.Z)(matchesEmptyString, (currType) => {
        return {
            message: "Token Type: ->" +
                currType.name +
                "<- static 'PATTERN' must not match an empty string",
            type: LexerDefinitionErrorType.EMPTY_MATCH_PATTERN,
            tokenTypes: [currType],
        };
    });
    return errors;
}
const start_of_input = /[^\\[][\^]|^\^/;
function findStartOfInputAnchor(tokenTypes) {
    class StartAnchorFinder extends api/* BaseRegExpVisitor */.e {
        constructor() {
            super(...arguments);
            this.found = false;
        }
        visitStartAnchor(node) {
            this.found = true;
        }
    }
    const invalidRegex = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        const pattern = currType.PATTERN;
        try {
            const regexpAst = getRegExpAst(pattern);
            const startAnchorVisitor = new StartAnchorFinder();
            startAnchorVisitor.visit(regexpAst);
            return startAnchorVisitor.found;
        }
        catch (e) {
            // old behavior in case of runtime exceptions with regexp-to-ast.
            /* istanbul ignore next - cannot ensure an error in regexp-to-ast*/
            return start_of_input.test(pattern.source);
        }
    });
    const errors = (0,map/* default */.Z)(invalidRegex, (currType) => {
        return {
            message: "Unexpected RegExp Anchor Error:\n" +
                "\tToken Type: ->" +
                currType.name +
                "<- static 'PATTERN' cannot contain start of input anchor '^'\n" +
                "\tSee https://chevrotain.io/docs/guide/resolving_lexer_errors.html#ANCHORS" +
                "\tfor details.",
            type: LexerDefinitionErrorType.SOI_ANCHOR_FOUND,
            tokenTypes: [currType],
        };
    });
    return errors;
}
function findUnsupportedFlags(tokenTypes) {
    const invalidFlags = (0,filter/* default */.Z)(tokenTypes, (currType) => {
        const pattern = currType[PATTERN];
        return pattern instanceof RegExp && (pattern.multiline || pattern.global);
    });
    const errors = (0,map/* default */.Z)(invalidFlags, (currType) => {
        return {
            message: "Token Type: ->" +
                currType.name +
                "<- static 'PATTERN' may NOT contain global('g') or multiline('m')",
            type: LexerDefinitionErrorType.UNSUPPORTED_FLAGS_FOUND,
            tokenTypes: [currType],
        };
    });
    return errors;
}
// This can only test for identical duplicate RegExps, not semantically equivalent ones.
function findDuplicatePatterns(tokenTypes) {
    const found = [];
    let identicalPatterns = (0,map/* default */.Z)(tokenTypes, (outerType) => {
        return (0,reduce/* default */.Z)(tokenTypes, (result, innerType) => {
            if (outerType.PATTERN.source === innerType.PATTERN.source &&
                !lodash_es_includes(found, innerType) &&
                innerType.PATTERN !== Lexer.NA) {
                // this avoids duplicates in the result, each Token Type may only appear in one "set"
                // in essence we are creating Equivalence classes on equality relation.
                found.push(innerType);
                result.push(innerType);
                return result;
            }
            return result;
        }, []);
    });
    identicalPatterns = lodash_es_compact(identicalPatterns);
    const duplicatePatterns = (0,filter/* default */.Z)(identicalPatterns, (currIdenticalSet) => {
        return currIdenticalSet.length > 1;
    });
    const errors = (0,map/* default */.Z)(duplicatePatterns, (setOfIdentical) => {
        const tokenTypeNames = (0,map/* default */.Z)(setOfIdentical, (currType) => {
            return currType.name;
        });
        const dupPatternSrc = lodash_es_head(setOfIdentical).PATTERN;
        return {
            message: `The same RegExp pattern ->${dupPatternSrc}<-` +
                `has been used in all of the following Token Types: ${tokenTypeNames.join(", ")} <-`,
            type: LexerDefinitionErrorType.DUPLICATE_PATTERNS_FOUND,
            tokenTypes: setOfIdentical,
        };
    });
    return errors;
}
function findInvalidGroupType(tokenTypes) {
    const invalidTypes = (0,filter/* default */.Z)(tokenTypes, (clazz) => {
        if (!(0,has/* default */.Z)(clazz, "GROUP")) {
            return false;
        }
        const group = clazz.GROUP;
        return group !== Lexer.SKIPPED && group !== Lexer.NA && !(0,isString/* default */.Z)(group);
    });
    const errors = (0,map/* default */.Z)(invalidTypes, (currType) => {
        return {
            message: "Token Type: ->" +
                currType.name +
                "<- static 'GROUP' can only be Lexer.SKIPPED/Lexer.NA/A String",
            type: LexerDefinitionErrorType.INVALID_GROUP_TYPE_FOUND,
            tokenTypes: [currType],
        };
    });
    return errors;
}
function findModesThatDoNotExist(tokenTypes, validModes) {
    const invalidModes = (0,filter/* default */.Z)(tokenTypes, (clazz) => {
        return (clazz.PUSH_MODE !== undefined && !lodash_es_includes(validModes, clazz.PUSH_MODE));
    });
    const errors = (0,map/* default */.Z)(invalidModes, (tokType) => {
        const msg = `Token Type: ->${tokType.name}<- static 'PUSH_MODE' value cannot refer to a Lexer Mode ->${tokType.PUSH_MODE}<-` +
            `which does not exist`;
        return {
            message: msg,
            type: LexerDefinitionErrorType.PUSH_MODE_DOES_NOT_EXIST,
            tokenTypes: [tokType],
        };
    });
    return errors;
}
function findUnreachablePatterns(tokenTypes) {
    const errors = [];
    const canBeTested = (0,reduce/* default */.Z)(tokenTypes, (result, tokType, idx) => {
        const pattern = tokType.PATTERN;
        if (pattern === Lexer.NA) {
            return result;
        }
        // a more comprehensive validation for all forms of regExps would require
        // deeper regExp analysis capabilities
        if ((0,isString/* default */.Z)(pattern)) {
            result.push({ str: pattern, idx, tokenType: tokType });
        }
        else if (lodash_es_isRegExp(pattern) && noMetaChar(pattern)) {
            result.push({ str: pattern.source, idx, tokenType: tokType });
        }
        return result;
    }, []);
    (0,forEach/* default */.Z)(tokenTypes, (tokType, testIdx) => {
        (0,forEach/* default */.Z)(canBeTested, ({ str, idx, tokenType }) => {
            if (testIdx < idx && testTokenType(str, tokType.PATTERN)) {
                const msg = `Token: ->${tokenType.name}<- can never be matched.\n` +
                    `Because it appears AFTER the Token Type ->${tokType.name}<-` +
                    `in the lexer's definition.\n` +
                    `See https://chevrotain.io/docs/guide/resolving_lexer_errors.html#UNREACHABLE`;
                errors.push({
                    message: msg,
                    type: LexerDefinitionErrorType.UNREACHABLE_PATTERN,
                    tokenTypes: [tokType, tokenType],
                });
            }
        });
    });
    return errors;
}
function testTokenType(str, pattern) {
    /* istanbul ignore else */
    if (lodash_es_isRegExp(pattern)) {
        const regExpArray = pattern.exec(str);
        return regExpArray !== null && regExpArray.index === 0;
    }
    else if ((0,isFunction/* default */.Z)(pattern)) {
        // maintain the API of custom patterns
        return pattern(str, 0, [], {});
    }
    else if ((0,has/* default */.Z)(pattern, "exec")) {
        // maintain the API of custom patterns
        return pattern.exec(str, 0, [], {});
    }
    else if (typeof pattern === "string") {
        return pattern === str;
    }
    else {
        throw Error("non exhaustive match");
    }
}
function noMetaChar(regExp) {
    //https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/RegExp
    const metaChars = [
        ".",
        "\\",
        "[",
        "]",
        "|",
        "^",
        "$",
        "(",
        ")",
        "?",
        "*",
        "+",
        "{",
    ];
    return ((0,find/* default */.Z)(metaChars, (char) => regExp.source.indexOf(char) !== -1) === undefined);
}
function addStartOfInput(pattern) {
    const flags = pattern.ignoreCase ? "i" : "";
    // always wrapping in a none capturing group preceded by '^' to make sure matching can only work on start of input.
    // duplicate/redundant start of input markers have no meaning (/^^^^A/ === /^A/)
    return new RegExp(`^(?:${pattern.source})`, flags);
}
function addStickyFlag(pattern) {
    const flags = pattern.ignoreCase ? "iy" : "y";
    // always wrapping in a none capturing group preceded by '^' to make sure matching can only work on start of input.
    // duplicate/redundant start of input markers have no meaning (/^^^^A/ === /^A/)
    return new RegExp(`${pattern.source}`, flags);
}
function performRuntimeChecks(lexerDefinition, trackLines, lineTerminatorCharacters) {
    const errors = [];
    // some run time checks to help the end users.
    if (!(0,has/* default */.Z)(lexerDefinition, DEFAULT_MODE)) {
        errors.push({
            message: "A MultiMode Lexer cannot be initialized without a <" +
                DEFAULT_MODE +
                "> property in its definition\n",
            type: LexerDefinitionErrorType.MULTI_MODE_LEXER_WITHOUT_DEFAULT_MODE,
        });
    }
    if (!(0,has/* default */.Z)(lexerDefinition, MODES)) {
        errors.push({
            message: "A MultiMode Lexer cannot be initialized without a <" +
                MODES +
                "> property in its definition\n",
            type: LexerDefinitionErrorType.MULTI_MODE_LEXER_WITHOUT_MODES_PROPERTY,
        });
    }
    if ((0,has/* default */.Z)(lexerDefinition, MODES) &&
        (0,has/* default */.Z)(lexerDefinition, DEFAULT_MODE) &&
        !(0,has/* default */.Z)(lexerDefinition.modes, lexerDefinition.defaultMode)) {
        errors.push({
            message: `A MultiMode Lexer cannot be initialized with a ${DEFAULT_MODE}: <${lexerDefinition.defaultMode}>` +
                `which does not exist\n`,
            type: LexerDefinitionErrorType.MULTI_MODE_LEXER_DEFAULT_MODE_VALUE_DOES_NOT_EXIST,
        });
    }
    if ((0,has/* default */.Z)(lexerDefinition, MODES)) {
        (0,forEach/* default */.Z)(lexerDefinition.modes, (currModeValue, currModeName) => {
            (0,forEach/* default */.Z)(currModeValue, (currTokType, currIdx) => {
                if ((0,isUndefined/* default */.Z)(currTokType)) {
                    errors.push({
                        message: `A Lexer cannot be initialized using an undefined Token Type. Mode:` +
                            `<${currModeName}> at index: <${currIdx}>\n`,
                        type: LexerDefinitionErrorType.LEXER_DEFINITION_CANNOT_CONTAIN_UNDEFINED,
                    });
                }
                else if ((0,has/* default */.Z)(currTokType, "LONGER_ALT")) {
                    const longerAlt = (0,isArray/* default */.Z)(currTokType.LONGER_ALT)
                        ? currTokType.LONGER_ALT
                        : [currTokType.LONGER_ALT];
                    (0,forEach/* default */.Z)(longerAlt, (currLongerAlt) => {
                        if (!(0,isUndefined/* default */.Z)(currLongerAlt) &&
                            !lodash_es_includes(currModeValue, currLongerAlt)) {
                            errors.push({
                                message: `A MultiMode Lexer cannot be initialized with a longer_alt <${currLongerAlt.name}> on token <${currTokType.name}> outside of mode <${currModeName}>\n`,
                                type: LexerDefinitionErrorType.MULTI_MODE_LEXER_LONGER_ALT_NOT_IN_CURRENT_MODE,
                            });
                        }
                    });
                }
            });
        });
    }
    return errors;
}
function performWarningRuntimeChecks(lexerDefinition, trackLines, lineTerminatorCharacters) {
    const warnings = [];
    let hasAnyLineBreak = false;
    const allTokenTypes = lodash_es_compact((0,flatten/* default */.Z)((0,values/* default */.Z)(lexerDefinition.modes)));
    const concreteTokenTypes = lodash_es_reject(allTokenTypes, (currType) => currType[PATTERN] === Lexer.NA);
    const terminatorCharCodes = getCharCodes(lineTerminatorCharacters);
    if (trackLines) {
        (0,forEach/* default */.Z)(concreteTokenTypes, (tokType) => {
            const currIssue = checkLineBreaksIssues(tokType, terminatorCharCodes);
            if (currIssue !== false) {
                const message = buildLineBreakIssueMessage(tokType, currIssue);
                const warningDescriptor = {
                    message,
                    type: currIssue.issue,
                    tokenType: tokType,
                };
                warnings.push(warningDescriptor);
            }
            else {
                // we don't want to attempt to scan if the user explicitly specified the line_breaks option.
                if ((0,has/* default */.Z)(tokType, "LINE_BREAKS")) {
                    if (tokType.LINE_BREAKS === true) {
                        hasAnyLineBreak = true;
                    }
                }
                else {
                    if (canMatchCharCode(terminatorCharCodes, tokType.PATTERN)) {
                        hasAnyLineBreak = true;
                    }
                }
            }
        });
    }
    if (trackLines && !hasAnyLineBreak) {
        warnings.push({
            message: "Warning: No LINE_BREAKS Found.\n" +
                "\tThis Lexer has been defined to track line and column information,\n" +
                "\tBut none of the Token Types can be identified as matching a line terminator.\n" +
                "\tSee https://chevrotain.io/docs/guide/resolving_lexer_errors.html#LINE_BREAKS \n" +
                "\tfor details.",
            type: LexerDefinitionErrorType.NO_LINE_BREAKS_FLAGS,
        });
    }
    return warnings;
}
function cloneEmptyGroups(emptyGroups) {
    const clonedResult = {};
    const groupKeys = (0,keys/* default */.Z)(emptyGroups);
    (0,forEach/* default */.Z)(groupKeys, (currKey) => {
        const currGroupValue = emptyGroups[currKey];
        /* istanbul ignore else */
        if ((0,isArray/* default */.Z)(currGroupValue)) {
            clonedResult[currKey] = [];
        }
        else {
            throw Error("non exhaustive match");
        }
    });
    return clonedResult;
}
// TODO: refactor to avoid duplication
function isCustomPattern(tokenType) {
    const pattern = tokenType.PATTERN;
    /* istanbul ignore else */
    if (lodash_es_isRegExp(pattern)) {
        return false;
    }
    else if ((0,isFunction/* default */.Z)(pattern)) {
        // CustomPatternMatcherFunc - custom patterns do not require any transformations, only wrapping in a RegExp Like object
        return true;
    }
    else if ((0,has/* default */.Z)(pattern, "exec")) {
        // ICustomPattern
        return true;
    }
    else if ((0,isString/* default */.Z)(pattern)) {
        return false;
    }
    else {
        throw Error("non exhaustive match");
    }
}
function isShortPattern(pattern) {
    if ((0,isString/* default */.Z)(pattern) && pattern.length === 1) {
        return pattern.charCodeAt(0);
    }
    else {
        return false;
    }
}
/**
 * Faster than using a RegExp for default newline detection during lexing.
 */
const LineTerminatorOptimizedTester = {
    // implements /\n|\r\n?/g.test
    test: function (text) {
        const len = text.length;
        for (let i = this.lastIndex; i < len; i++) {
            const c = text.charCodeAt(i);
            if (c === 10) {
                this.lastIndex = i + 1;
                return true;
            }
            else if (c === 13) {
                if (text.charCodeAt(i + 1) === 10) {
                    this.lastIndex = i + 2;
                }
                else {
                    this.lastIndex = i + 1;
                }
                return true;
            }
        }
        return false;
    },
    lastIndex: 0,
};
function checkLineBreaksIssues(tokType, lineTerminatorCharCodes) {
    if ((0,has/* default */.Z)(tokType, "LINE_BREAKS")) {
        // if the user explicitly declared the line_breaks option we will respect their choice
        // and assume it is correct.
        return false;
    }
    else {
        /* istanbul ignore else */
        if (lodash_es_isRegExp(tokType.PATTERN)) {
            try {
                // TODO: why is the casting suddenly needed?
                canMatchCharCode(lineTerminatorCharCodes, tokType.PATTERN);
            }
            catch (e) {
                /* istanbul ignore next - to test this we would have to mock <canMatchCharCode> to throw an error */
                return {
                    issue: LexerDefinitionErrorType.IDENTIFY_TERMINATOR,
                    errMsg: e.message,
                };
            }
            return false;
        }
        else if ((0,isString/* default */.Z)(tokType.PATTERN)) {
            // string literal patterns can always be analyzed to detect line terminator usage
            return false;
        }
        else if (isCustomPattern(tokType)) {
            // custom token types
            return { issue: LexerDefinitionErrorType.CUSTOM_LINE_BREAK };
        }
        else {
            throw Error("non exhaustive match");
        }
    }
}
function buildLineBreakIssueMessage(tokType, details) {
    /* istanbul ignore else */
    if (details.issue === LexerDefinitionErrorType.IDENTIFY_TERMINATOR) {
        return ("Warning: unable to identify line terminator usage in pattern.\n" +
            `\tThe problem is in the <${tokType.name}> Token Type\n` +
            `\t Root cause: ${details.errMsg}.\n` +
            "\tFor details See: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#IDENTIFY_TERMINATOR");
    }
    else if (details.issue === LexerDefinitionErrorType.CUSTOM_LINE_BREAK) {
        return ("Warning: A Custom Token Pattern should specify the <line_breaks> option.\n" +
            `\tThe problem is in the <${tokType.name}> Token Type\n` +
            "\tFor details See: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#CUSTOM_LINE_BREAK");
    }
    else {
        throw Error("non exhaustive match");
    }
}
function getCharCodes(charsOrCodes) {
    const charCodes = (0,map/* default */.Z)(charsOrCodes, (numOrString) => {
        if ((0,isString/* default */.Z)(numOrString)) {
            return numOrString.charCodeAt(0);
        }
        else {
            return numOrString;
        }
    });
    return charCodes;
}
function addToMapOfArrays(map, key, value) {
    if (map[key] === undefined) {
        map[key] = [value];
    }
    else {
        map[key].push(value);
    }
}
const minOptimizationVal = 256;
/**
 * We are mapping charCode above ASCI (256) into buckets each in the size of 256.
 * This is because ASCI are the most common start chars so each one of those will get its own
 * possible token configs vector.
 *
 * Tokens starting with charCodes "above" ASCI are uncommon, so we can "afford"
 * to place these into buckets of possible token configs, What we gain from
 * this is avoiding the case of creating an optimization 'charCodeToPatternIdxToConfig'
 * which would contain 10,000+ arrays of small size (e.g unicode Identifiers scenario).
 * Our 'charCodeToPatternIdxToConfig' max size will now be:
 * 256 + (2^16 / 2^8) - 1 === 511
 *
 * note the hack for fast division integer part extraction
 * See: https://stackoverflow.com/a/4228528
 */
let charCodeToOptimizedIdxMap = [];
function charCodeToOptimizedIndex(charCode) {
    return charCode < minOptimizationVal
        ? charCode
        : charCodeToOptimizedIdxMap[charCode];
}
/**
 * This is a compromise between cold start / hot running performance
 * Creating this array takes ~3ms on a modern machine,
 * But if we perform the computation at runtime as needed the CSS Lexer benchmark
 * performance degrades by ~10%
 *
 * TODO: Perhaps it should be lazy initialized only if a charCode > 255 is used.
 */
function initCharCodeToOptimizedIndexMap() {
    if ((0,isEmpty/* default */.Z)(charCodeToOptimizedIdxMap)) {
        charCodeToOptimizedIdxMap = new Array(65536);
        for (let i = 0; i < 65536; i++) {
            charCodeToOptimizedIdxMap[i] = i > 255 ? 255 + ~~(i / 255) : i;
        }
    }
}
//# sourceMappingURL=lexer.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/identity.js
var identity = __webpack_require__(64056);
// EXTERNAL MODULE: ../node_modules/lodash-es/noop.js
var noop = __webpack_require__(10152);
// EXTERNAL MODULE: ../node_modules/lodash-es/last.js
var last = __webpack_require__(36411);
;// CONCATENATED MODULE: ../node_modules/@chevrotain/utils/lib/src/timer.js
function timer(func) {
    const start = new Date().getTime();
    const val = func();
    const end = new Date().getTime();
    const total = end - start;
    return { time: total, value: val };
}
//# sourceMappingURL=timer.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/tokens.js

function tokenStructuredMatcher(tokInstance, tokConstructor) {
    const instanceType = tokInstance.tokenTypeIdx;
    if (instanceType === tokConstructor.tokenTypeIdx) {
        return true;
    }
    else {
        return (tokConstructor.isParent === true &&
            tokConstructor.categoryMatchesMap[instanceType] === true);
    }
}
// Optimized tokenMatcher in case our grammar does not use token categories
// Being so tiny it is much more likely to be in-lined and this avoid the function call overhead
function tokenStructuredMatcherNoCategories(token, tokType) {
    return token.tokenTypeIdx === tokType.tokenTypeIdx;
}
let tokenShortNameIdx = 1;
const tokenIdxToClass = {};
function augmentTokenTypes(tokenTypes) {
    // collect the parent Token Types as well.
    const tokenTypesAndParents = expandCategories(tokenTypes);
    // add required tokenType and categoryMatches properties
    assignTokenDefaultProps(tokenTypesAndParents);
    // fill up the categoryMatches
    assignCategoriesMapProp(tokenTypesAndParents);
    assignCategoriesTokensProp(tokenTypesAndParents);
    (0,forEach/* default */.Z)(tokenTypesAndParents, (tokType) => {
        tokType.isParent = tokType.categoryMatches.length > 0;
    });
}
function expandCategories(tokenTypes) {
    let result = (0,lodash_es_clone/* default */.Z)(tokenTypes);
    let categories = tokenTypes;
    let searching = true;
    while (searching) {
        categories = lodash_es_compact((0,flatten/* default */.Z)((0,map/* default */.Z)(categories, (currTokType) => currTokType.CATEGORIES)));
        const newCategories = lodash_es_difference(categories, result);
        result = result.concat(newCategories);
        if ((0,isEmpty/* default */.Z)(newCategories)) {
            searching = false;
        }
        else {
            categories = newCategories;
        }
    }
    return result;
}
function assignTokenDefaultProps(tokenTypes) {
    (0,forEach/* default */.Z)(tokenTypes, (currTokType) => {
        if (!hasShortKeyProperty(currTokType)) {
            tokenIdxToClass[tokenShortNameIdx] = currTokType;
            currTokType.tokenTypeIdx = tokenShortNameIdx++;
        }
        // CATEGORIES? : TokenType | TokenType[]
        if (hasCategoriesProperty(currTokType) &&
            !(0,isArray/* default */.Z)(currTokType.CATEGORIES)
        // &&
        // !isUndefined(currTokType.CATEGORIES.PATTERN)
        ) {
            currTokType.CATEGORIES = [currTokType.CATEGORIES];
        }
        if (!hasCategoriesProperty(currTokType)) {
            currTokType.CATEGORIES = [];
        }
        if (!hasExtendingTokensTypesProperty(currTokType)) {
            currTokType.categoryMatches = [];
        }
        if (!hasExtendingTokensTypesMapProperty(currTokType)) {
            currTokType.categoryMatchesMap = {};
        }
    });
}
function assignCategoriesTokensProp(tokenTypes) {
    (0,forEach/* default */.Z)(tokenTypes, (currTokType) => {
        // avoid duplications
        currTokType.categoryMatches = [];
        (0,forEach/* default */.Z)(currTokType.categoryMatchesMap, (val, key) => {
            currTokType.categoryMatches.push(tokenIdxToClass[key].tokenTypeIdx);
        });
    });
}
function assignCategoriesMapProp(tokenTypes) {
    (0,forEach/* default */.Z)(tokenTypes, (currTokType) => {
        singleAssignCategoriesToksMap([], currTokType);
    });
}
function singleAssignCategoriesToksMap(path, nextNode) {
    (0,forEach/* default */.Z)(path, (pathNode) => {
        nextNode.categoryMatchesMap[pathNode.tokenTypeIdx] = true;
    });
    (0,forEach/* default */.Z)(nextNode.CATEGORIES, (nextCategory) => {
        const newPath = path.concat(nextNode);
        // avoids infinite loops due to cyclic categories.
        if (!lodash_es_includes(newPath, nextCategory)) {
            singleAssignCategoriesToksMap(newPath, nextCategory);
        }
    });
}
function hasShortKeyProperty(tokType) {
    return (0,has/* default */.Z)(tokType, "tokenTypeIdx");
}
function hasCategoriesProperty(tokType) {
    return (0,has/* default */.Z)(tokType, "CATEGORIES");
}
function hasExtendingTokensTypesProperty(tokType) {
    return (0,has/* default */.Z)(tokType, "categoryMatches");
}
function hasExtendingTokensTypesMapProperty(tokType) {
    return (0,has/* default */.Z)(tokType, "categoryMatchesMap");
}
function isTokenType(tokType) {
    return (0,has/* default */.Z)(tokType, "tokenTypeIdx");
}
//# sourceMappingURL=tokens.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/lexer_errors_public.js
const defaultLexerErrorProvider = {
    buildUnableToPopLexerModeMessage(token) {
        return `Unable to pop Lexer Mode after encountering Token ->${token.image}<- The Mode Stack is empty`;
    },
    buildUnexpectedCharactersMessage(fullText, startOffset, length, line, column) {
        return (`unexpected character: ->${fullText.charAt(startOffset)}<- at offset: ${startOffset},` + ` skipped ${length} characters.`);
    },
};
//# sourceMappingURL=lexer_errors_public.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/lexer_public.js






var LexerDefinitionErrorType;
(function (LexerDefinitionErrorType) {
    LexerDefinitionErrorType[LexerDefinitionErrorType["MISSING_PATTERN"] = 0] = "MISSING_PATTERN";
    LexerDefinitionErrorType[LexerDefinitionErrorType["INVALID_PATTERN"] = 1] = "INVALID_PATTERN";
    LexerDefinitionErrorType[LexerDefinitionErrorType["EOI_ANCHOR_FOUND"] = 2] = "EOI_ANCHOR_FOUND";
    LexerDefinitionErrorType[LexerDefinitionErrorType["UNSUPPORTED_FLAGS_FOUND"] = 3] = "UNSUPPORTED_FLAGS_FOUND";
    LexerDefinitionErrorType[LexerDefinitionErrorType["DUPLICATE_PATTERNS_FOUND"] = 4] = "DUPLICATE_PATTERNS_FOUND";
    LexerDefinitionErrorType[LexerDefinitionErrorType["INVALID_GROUP_TYPE_FOUND"] = 5] = "INVALID_GROUP_TYPE_FOUND";
    LexerDefinitionErrorType[LexerDefinitionErrorType["PUSH_MODE_DOES_NOT_EXIST"] = 6] = "PUSH_MODE_DOES_NOT_EXIST";
    LexerDefinitionErrorType[LexerDefinitionErrorType["MULTI_MODE_LEXER_WITHOUT_DEFAULT_MODE"] = 7] = "MULTI_MODE_LEXER_WITHOUT_DEFAULT_MODE";
    LexerDefinitionErrorType[LexerDefinitionErrorType["MULTI_MODE_LEXER_WITHOUT_MODES_PROPERTY"] = 8] = "MULTI_MODE_LEXER_WITHOUT_MODES_PROPERTY";
    LexerDefinitionErrorType[LexerDefinitionErrorType["MULTI_MODE_LEXER_DEFAULT_MODE_VALUE_DOES_NOT_EXIST"] = 9] = "MULTI_MODE_LEXER_DEFAULT_MODE_VALUE_DOES_NOT_EXIST";
    LexerDefinitionErrorType[LexerDefinitionErrorType["LEXER_DEFINITION_CANNOT_CONTAIN_UNDEFINED"] = 10] = "LEXER_DEFINITION_CANNOT_CONTAIN_UNDEFINED";
    LexerDefinitionErrorType[LexerDefinitionErrorType["SOI_ANCHOR_FOUND"] = 11] = "SOI_ANCHOR_FOUND";
    LexerDefinitionErrorType[LexerDefinitionErrorType["EMPTY_MATCH_PATTERN"] = 12] = "EMPTY_MATCH_PATTERN";
    LexerDefinitionErrorType[LexerDefinitionErrorType["NO_LINE_BREAKS_FLAGS"] = 13] = "NO_LINE_BREAKS_FLAGS";
    LexerDefinitionErrorType[LexerDefinitionErrorType["UNREACHABLE_PATTERN"] = 14] = "UNREACHABLE_PATTERN";
    LexerDefinitionErrorType[LexerDefinitionErrorType["IDENTIFY_TERMINATOR"] = 15] = "IDENTIFY_TERMINATOR";
    LexerDefinitionErrorType[LexerDefinitionErrorType["CUSTOM_LINE_BREAK"] = 16] = "CUSTOM_LINE_BREAK";
    LexerDefinitionErrorType[LexerDefinitionErrorType["MULTI_MODE_LEXER_LONGER_ALT_NOT_IN_CURRENT_MODE"] = 17] = "MULTI_MODE_LEXER_LONGER_ALT_NOT_IN_CURRENT_MODE";
})(LexerDefinitionErrorType || (LexerDefinitionErrorType = {}));
const DEFAULT_LEXER_CONFIG = {
    deferDefinitionErrorsHandling: false,
    positionTracking: "full",
    lineTerminatorsPattern: /\n|\r\n?/g,
    lineTerminatorCharacters: ["\n", "\r"],
    ensureOptimizations: false,
    safeMode: false,
    errorMessageProvider: defaultLexerErrorProvider,
    traceInitPerf: false,
    skipValidations: false,
    recoveryEnabled: true,
};
Object.freeze(DEFAULT_LEXER_CONFIG);
class Lexer {
    constructor(lexerDefinition, config = DEFAULT_LEXER_CONFIG) {
        this.lexerDefinition = lexerDefinition;
        this.lexerDefinitionErrors = [];
        this.lexerDefinitionWarning = [];
        this.patternIdxToConfig = {};
        this.charCodeToPatternIdxToConfig = {};
        this.modes = [];
        this.emptyGroups = {};
        this.trackStartLines = true;
        this.trackEndLines = true;
        this.hasCustom = false;
        this.canModeBeOptimized = {};
        // Duplicated from the parser's perf trace trait to allow future extraction
        // of the lexer to a separate package.
        this.TRACE_INIT = (phaseDesc, phaseImpl) => {
            // No need to optimize this using NOOP pattern because
            // It is not called in a hot spot...
            if (this.traceInitPerf === true) {
                this.traceInitIndent++;
                const indent = new Array(this.traceInitIndent + 1).join("\t");
                if (this.traceInitIndent < this.traceInitMaxIdent) {
                    console.log(`${indent}--> <${phaseDesc}>`);
                }
                const { time, value } = timer(phaseImpl);
                /* istanbul ignore next - Difficult to reproduce specific performance behavior (>10ms) in tests */
                const traceMethod = time > 10 ? console.warn : console.log;
                if (this.traceInitIndent < this.traceInitMaxIdent) {
                    traceMethod(`${indent}<-- <${phaseDesc}> time: ${time}ms`);
                }
                this.traceInitIndent--;
                return value;
            }
            else {
                return phaseImpl();
            }
        };
        if (typeof config === "boolean") {
            throw Error("The second argument to the Lexer constructor is now an ILexerConfig Object.\n" +
                "a boolean 2nd argument is no longer supported");
        }
        // todo: defaults func?
        this.config = lodash_es_assign({}, DEFAULT_LEXER_CONFIG, config);
        const traceInitVal = this.config.traceInitPerf;
        if (traceInitVal === true) {
            this.traceInitMaxIdent = Infinity;
            this.traceInitPerf = true;
        }
        else if (typeof traceInitVal === "number") {
            this.traceInitMaxIdent = traceInitVal;
            this.traceInitPerf = true;
        }
        this.traceInitIndent = -1;
        this.TRACE_INIT("Lexer Constructor", () => {
            let actualDefinition;
            let hasOnlySingleMode = true;
            this.TRACE_INIT("Lexer Config handling", () => {
                if (this.config.lineTerminatorsPattern ===
                    DEFAULT_LEXER_CONFIG.lineTerminatorsPattern) {
                    // optimized built-in implementation for the defaults definition of lineTerminators
                    this.config.lineTerminatorsPattern = LineTerminatorOptimizedTester;
                }
                else {
                    if (this.config.lineTerminatorCharacters ===
                        DEFAULT_LEXER_CONFIG.lineTerminatorCharacters) {
                        throw Error("Error: Missing <lineTerminatorCharacters> property on the Lexer config.\n" +
                            "\tFor details See: https://chevrotain.io/docs/guide/resolving_lexer_errors.html#MISSING_LINE_TERM_CHARS");
                    }
                }
                if (config.safeMode && config.ensureOptimizations) {
                    throw Error('"safeMode" and "ensureOptimizations" flags are mutually exclusive.');
                }
                this.trackStartLines = /full|onlyStart/i.test(this.config.positionTracking);
                this.trackEndLines = /full/i.test(this.config.positionTracking);
                // Convert SingleModeLexerDefinition into a IMultiModeLexerDefinition.
                if ((0,isArray/* default */.Z)(lexerDefinition)) {
                    actualDefinition = {
                        modes: { defaultMode: (0,lodash_es_clone/* default */.Z)(lexerDefinition) },
                        defaultMode: DEFAULT_MODE,
                    };
                }
                else {
                    // no conversion needed, input should already be a IMultiModeLexerDefinition
                    hasOnlySingleMode = false;
                    actualDefinition = (0,lodash_es_clone/* default */.Z)(lexerDefinition);
                }
            });
            if (this.config.skipValidations === false) {
                this.TRACE_INIT("performRuntimeChecks", () => {
                    this.lexerDefinitionErrors = this.lexerDefinitionErrors.concat(performRuntimeChecks(actualDefinition, this.trackStartLines, this.config.lineTerminatorCharacters));
                });
                this.TRACE_INIT("performWarningRuntimeChecks", () => {
                    this.lexerDefinitionWarning = this.lexerDefinitionWarning.concat(performWarningRuntimeChecks(actualDefinition, this.trackStartLines, this.config.lineTerminatorCharacters));
                });
            }
            // for extra robustness to avoid throwing an none informative error message
            actualDefinition.modes = actualDefinition.modes
                ? actualDefinition.modes
                : {};
            // an error of undefined TokenTypes will be detected in "performRuntimeChecks" above.
            // this transformation is to increase robustness in the case of partially invalid lexer definition.
            (0,forEach/* default */.Z)(actualDefinition.modes, (currModeValue, currModeName) => {
                actualDefinition.modes[currModeName] = lodash_es_reject(currModeValue, (currTokType) => (0,isUndefined/* default */.Z)(currTokType));
            });
            const allModeNames = (0,keys/* default */.Z)(actualDefinition.modes);
            (0,forEach/* default */.Z)(actualDefinition.modes, (currModDef, currModName) => {
                this.TRACE_INIT(`Mode: <${currModName}> processing`, () => {
                    this.modes.push(currModName);
                    if (this.config.skipValidations === false) {
                        this.TRACE_INIT(`validatePatterns`, () => {
                            this.lexerDefinitionErrors = this.lexerDefinitionErrors.concat(validatePatterns(currModDef, allModeNames));
                        });
                    }
                    // If definition errors were encountered, the analysis phase may fail unexpectedly/
                    // Considering a lexer with definition errors may never be used, there is no point
                    // to performing the analysis anyhow...
                    if ((0,isEmpty/* default */.Z)(this.lexerDefinitionErrors)) {
                        augmentTokenTypes(currModDef);
                        let currAnalyzeResult;
                        this.TRACE_INIT(`analyzeTokenTypes`, () => {
                            currAnalyzeResult = analyzeTokenTypes(currModDef, {
                                lineTerminatorCharacters: this.config.lineTerminatorCharacters,
                                positionTracking: config.positionTracking,
                                ensureOptimizations: config.ensureOptimizations,
                                safeMode: config.safeMode,
                                tracer: this.TRACE_INIT,
                            });
                        });
                        this.patternIdxToConfig[currModName] =
                            currAnalyzeResult.patternIdxToConfig;
                        this.charCodeToPatternIdxToConfig[currModName] =
                            currAnalyzeResult.charCodeToPatternIdxToConfig;
                        this.emptyGroups = lodash_es_assign({}, this.emptyGroups, currAnalyzeResult.emptyGroups);
                        this.hasCustom = currAnalyzeResult.hasCustom || this.hasCustom;
                        this.canModeBeOptimized[currModName] =
                            currAnalyzeResult.canBeOptimized;
                    }
                });
            });
            this.defaultMode = actualDefinition.defaultMode;
            if (!(0,isEmpty/* default */.Z)(this.lexerDefinitionErrors) &&
                !this.config.deferDefinitionErrorsHandling) {
                const allErrMessages = (0,map/* default */.Z)(this.lexerDefinitionErrors, (error) => {
                    return error.message;
                });
                const allErrMessagesString = allErrMessages.join("-----------------------\n");
                throw new Error("Errors detected in definition of Lexer:\n" + allErrMessagesString);
            }
            // Only print warning if there are no errors, This will avoid pl
            (0,forEach/* default */.Z)(this.lexerDefinitionWarning, (warningDescriptor) => {
                PRINT_WARNING(warningDescriptor.message);
            });
            this.TRACE_INIT("Choosing sub-methods implementations", () => {
                // Choose the relevant internal implementations for this specific parser.
                // These implementations should be in-lined by the JavaScript engine
                // to provide optimal performance in each scenario.
                if (SUPPORT_STICKY) {
                    this.chopInput = identity/* default */.Z;
                    this.match = this.matchWithTest;
                }
                else {
                    this.updateLastIndex = noop/* default */.Z;
                    this.match = this.matchWithExec;
                }
                if (hasOnlySingleMode) {
                    this.handleModes = noop/* default */.Z;
                }
                if (this.trackStartLines === false) {
                    this.computeNewColumn = identity/* default */.Z;
                }
                if (this.trackEndLines === false) {
                    this.updateTokenEndLineColumnLocation = noop/* default */.Z;
                }
                if (/full/i.test(this.config.positionTracking)) {
                    this.createTokenInstance = this.createFullToken;
                }
                else if (/onlyStart/i.test(this.config.positionTracking)) {
                    this.createTokenInstance = this.createStartOnlyToken;
                }
                else if (/onlyOffset/i.test(this.config.positionTracking)) {
                    this.createTokenInstance = this.createOffsetOnlyToken;
                }
                else {
                    throw Error(`Invalid <positionTracking> config option: "${this.config.positionTracking}"`);
                }
                if (this.hasCustom) {
                    this.addToken = this.addTokenUsingPush;
                    this.handlePayload = this.handlePayloadWithCustom;
                }
                else {
                    this.addToken = this.addTokenUsingMemberAccess;
                    this.handlePayload = this.handlePayloadNoCustom;
                }
            });
            this.TRACE_INIT("Failed Optimization Warnings", () => {
                const unOptimizedModes = (0,reduce/* default */.Z)(this.canModeBeOptimized, (cannotBeOptimized, canBeOptimized, modeName) => {
                    if (canBeOptimized === false) {
                        cannotBeOptimized.push(modeName);
                    }
                    return cannotBeOptimized;
                }, []);
                if (config.ensureOptimizations && !(0,isEmpty/* default */.Z)(unOptimizedModes)) {
                    throw Error(`Lexer Modes: < ${unOptimizedModes.join(", ")} > cannot be optimized.\n` +
                        '\t Disable the "ensureOptimizations" lexer config flag to silently ignore this and run the lexer in an un-optimized mode.\n' +
                        "\t Or inspect the console log for details on how to resolve these issues.");
                }
            });
            this.TRACE_INIT("clearRegExpParserCache", () => {
                clearRegExpParserCache();
            });
            this.TRACE_INIT("toFastProperties", () => {
                toFastProperties(this);
            });
        });
    }
    tokenize(text, initialMode = this.defaultMode) {
        if (!(0,isEmpty/* default */.Z)(this.lexerDefinitionErrors)) {
            const allErrMessages = (0,map/* default */.Z)(this.lexerDefinitionErrors, (error) => {
                return error.message;
            });
            const allErrMessagesString = allErrMessages.join("-----------------------\n");
            throw new Error("Unable to Tokenize because Errors detected in definition of Lexer:\n" +
                allErrMessagesString);
        }
        return this.tokenizeInternal(text, initialMode);
    }
    // There is quite a bit of duplication between this and "tokenizeInternalLazy"
    // This is intentional due to performance considerations.
    // this method also used quite a bit of `!` none null assertions because it is too optimized
    // for `tsc` to always understand it is "safe"
    tokenizeInternal(text, initialMode) {
        let i, j, k, matchAltImage, longerAlt, matchedImage, payload, altPayload, imageLength, group, tokType, newToken, errLength, droppedChar, msg, match;
        const orgText = text;
        const orgLength = orgText.length;
        let offset = 0;
        let matchedTokensIndex = 0;
        // initializing the tokensArray to the "guessed" size.
        // guessing too little will still reduce the number of array re-sizes on pushes.
        // guessing too large (Tested by guessing x4 too large) may cost a bit more of memory
        // but would still have a faster runtime by avoiding (All but one) array resizing.
        const guessedNumberOfTokens = this.hasCustom
            ? 0 // will break custom token pattern APIs the matchedTokens array will contain undefined elements.
            : Math.floor(text.length / 10);
        const matchedTokens = new Array(guessedNumberOfTokens);
        const errors = [];
        let line = this.trackStartLines ? 1 : undefined;
        let column = this.trackStartLines ? 1 : undefined;
        const groups = cloneEmptyGroups(this.emptyGroups);
        const trackLines = this.trackStartLines;
        const lineTerminatorPattern = this.config.lineTerminatorsPattern;
        let currModePatternsLength = 0;
        let patternIdxToConfig = [];
        let currCharCodeToPatternIdxToConfig = [];
        const modeStack = [];
        const emptyArray = [];
        Object.freeze(emptyArray);
        let getPossiblePatterns;
        function getPossiblePatternsSlow() {
            return patternIdxToConfig;
        }
        function getPossiblePatternsOptimized(charCode) {
            const optimizedCharIdx = charCodeToOptimizedIndex(charCode);
            const possiblePatterns = currCharCodeToPatternIdxToConfig[optimizedCharIdx];
            if (possiblePatterns === undefined) {
                return emptyArray;
            }
            else {
                return possiblePatterns;
            }
        }
        const pop_mode = (popToken) => {
            // TODO: perhaps avoid this error in the edge case there is no more input?
            if (modeStack.length === 1 &&
                // if we have both a POP_MODE and a PUSH_MODE this is in-fact a "transition"
                // So no error should occur.
                popToken.tokenType.PUSH_MODE === undefined) {
                // if we try to pop the last mode there lexer will no longer have ANY mode.
                // thus the pop is ignored, an error will be created and the lexer will continue parsing in the previous mode.
                const msg = this.config.errorMessageProvider.buildUnableToPopLexerModeMessage(popToken);
                errors.push({
                    offset: popToken.startOffset,
                    line: popToken.startLine,
                    column: popToken.startColumn,
                    length: popToken.image.length,
                    message: msg,
                });
            }
            else {
                modeStack.pop();
                const newMode = (0,last/* default */.Z)(modeStack);
                patternIdxToConfig = this.patternIdxToConfig[newMode];
                currCharCodeToPatternIdxToConfig =
                    this.charCodeToPatternIdxToConfig[newMode];
                currModePatternsLength = patternIdxToConfig.length;
                const modeCanBeOptimized = this.canModeBeOptimized[newMode] && this.config.safeMode === false;
                if (currCharCodeToPatternIdxToConfig && modeCanBeOptimized) {
                    getPossiblePatterns = getPossiblePatternsOptimized;
                }
                else {
                    getPossiblePatterns = getPossiblePatternsSlow;
                }
            }
        };
        function push_mode(newMode) {
            modeStack.push(newMode);
            currCharCodeToPatternIdxToConfig =
                this.charCodeToPatternIdxToConfig[newMode];
            patternIdxToConfig = this.patternIdxToConfig[newMode];
            currModePatternsLength = patternIdxToConfig.length;
            currModePatternsLength = patternIdxToConfig.length;
            const modeCanBeOptimized = this.canModeBeOptimized[newMode] && this.config.safeMode === false;
            if (currCharCodeToPatternIdxToConfig && modeCanBeOptimized) {
                getPossiblePatterns = getPossiblePatternsOptimized;
            }
            else {
                getPossiblePatterns = getPossiblePatternsSlow;
            }
        }
        // this pattern seems to avoid a V8 de-optimization, although that de-optimization does not
        // seem to matter performance wise.
        push_mode.call(this, initialMode);
        let currConfig;
        const recoveryEnabled = this.config.recoveryEnabled;
        while (offset < orgLength) {
            matchedImage = null;
            const nextCharCode = orgText.charCodeAt(offset);
            const chosenPatternIdxToConfig = getPossiblePatterns(nextCharCode);
            const chosenPatternsLength = chosenPatternIdxToConfig.length;
            for (i = 0; i < chosenPatternsLength; i++) {
                currConfig = chosenPatternIdxToConfig[i];
                const currPattern = currConfig.pattern;
                payload = null;
                // manually in-lined because > 600 chars won't be in-lined in V8
                const singleCharCode = currConfig.short;
                if (singleCharCode !== false) {
                    if (nextCharCode === singleCharCode) {
                        // single character string
                        matchedImage = currPattern;
                    }
                }
                else if (currConfig.isCustom === true) {
                    match = currPattern.exec(orgText, offset, matchedTokens, groups);
                    if (match !== null) {
                        matchedImage = match[0];
                        if (match.payload !== undefined) {
                            payload = match.payload;
                        }
                    }
                    else {
                        matchedImage = null;
                    }
                }
                else {
                    this.updateLastIndex(currPattern, offset);
                    matchedImage = this.match(currPattern, text, offset);
                }
                if (matchedImage !== null) {
                    // even though this pattern matched we must try a another longer alternative.
                    // this can be used to prioritize keywords over identifiers
                    longerAlt = currConfig.longerAlt;
                    if (longerAlt !== undefined) {
                        // TODO: micro optimize, avoid extra prop access
                        // by saving/linking longerAlt on the original config?
                        const longerAltLength = longerAlt.length;
                        for (k = 0; k < longerAltLength; k++) {
                            const longerAltConfig = patternIdxToConfig[longerAlt[k]];
                            const longerAltPattern = longerAltConfig.pattern;
                            altPayload = null;
                            // single Char can never be a longer alt so no need to test it.
                            // manually in-lined because > 600 chars won't be in-lined in V8
                            if (longerAltConfig.isCustom === true) {
                                match = longerAltPattern.exec(orgText, offset, matchedTokens, groups);
                                if (match !== null) {
                                    matchAltImage = match[0];
                                    if (match.payload !== undefined) {
                                        altPayload = match.payload;
                                    }
                                }
                                else {
                                    matchAltImage = null;
                                }
                            }
                            else {
                                this.updateLastIndex(longerAltPattern, offset);
                                matchAltImage = this.match(longerAltPattern, text, offset);
                            }
                            if (matchAltImage && matchAltImage.length > matchedImage.length) {
                                matchedImage = matchAltImage;
                                payload = altPayload;
                                currConfig = longerAltConfig;
                                // Exit the loop early after matching one of the longer alternatives
                                // The first matched alternative takes precedence
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            // successful match
            if (matchedImage !== null) {
                imageLength = matchedImage.length;
                group = currConfig.group;
                if (group !== undefined) {
                    tokType = currConfig.tokenTypeIdx;
                    // TODO: "offset + imageLength" and the new column may be computed twice in case of "full" location information inside
                    // createFullToken method
                    newToken = this.createTokenInstance(matchedImage, offset, tokType, currConfig.tokenType, line, column, imageLength);
                    this.handlePayload(newToken, payload);
                    // TODO: optimize NOOP in case there are no special groups?
                    if (group === false) {
                        matchedTokensIndex = this.addToken(matchedTokens, matchedTokensIndex, newToken);
                    }
                    else {
                        groups[group].push(newToken);
                    }
                }
                text = this.chopInput(text, imageLength);
                offset = offset + imageLength;
                // TODO: with newlines the column may be assigned twice
                column = this.computeNewColumn(column, imageLength);
                if (trackLines === true && currConfig.canLineTerminator === true) {
                    let numOfLTsInMatch = 0;
                    let foundTerminator;
                    let lastLTEndOffset;
                    lineTerminatorPattern.lastIndex = 0;
                    do {
                        foundTerminator = lineTerminatorPattern.test(matchedImage);
                        if (foundTerminator === true) {
                            lastLTEndOffset = lineTerminatorPattern.lastIndex - 1;
                            numOfLTsInMatch++;
                        }
                    } while (foundTerminator === true);
                    if (numOfLTsInMatch !== 0) {
                        line = line + numOfLTsInMatch;
                        column = imageLength - lastLTEndOffset;
                        this.updateTokenEndLineColumnLocation(newToken, group, lastLTEndOffset, numOfLTsInMatch, line, column, imageLength);
                    }
                }
                // will be NOOP if no modes present
                this.handleModes(currConfig, pop_mode, push_mode, newToken);
            }
            else {
                // error recovery, drop characters until we identify a valid token's start point
                const errorStartOffset = offset;
                const errorLine = line;
                const errorColumn = column;
                let foundResyncPoint = recoveryEnabled === false;
                while (foundResyncPoint === false && offset < orgLength) {
                    // Identity Func (when sticky flag is enabled)
                    text = this.chopInput(text, 1);
                    offset++;
                    for (j = 0; j < currModePatternsLength; j++) {
                        const currConfig = patternIdxToConfig[j];
                        const currPattern = currConfig.pattern;
                        // manually in-lined because > 600 chars won't be in-lined in V8
                        const singleCharCode = currConfig.short;
                        if (singleCharCode !== false) {
                            if (orgText.charCodeAt(offset) === singleCharCode) {
                                // single character string
                                foundResyncPoint = true;
                            }
                        }
                        else if (currConfig.isCustom === true) {
                            foundResyncPoint =
                                currPattern.exec(orgText, offset, matchedTokens, groups) !== null;
                        }
                        else {
                            this.updateLastIndex(currPattern, offset);
                            foundResyncPoint = currPattern.exec(text) !== null;
                        }
                        if (foundResyncPoint === true) {
                            break;
                        }
                    }
                }
                errLength = offset - errorStartOffset;
                column = this.computeNewColumn(column, errLength);
                // at this point we either re-synced or reached the end of the input text
                msg = this.config.errorMessageProvider.buildUnexpectedCharactersMessage(orgText, errorStartOffset, errLength, errorLine, errorColumn);
                errors.push({
                    offset: errorStartOffset,
                    line: errorLine,
                    column: errorColumn,
                    length: errLength,
                    message: msg,
                });
                if (recoveryEnabled === false) {
                    break;
                }
            }
        }
        // if we do have custom patterns which push directly into the
        // TODO: custom tokens should not push directly??
        if (!this.hasCustom) {
            // if we guessed a too large size for the tokens array this will shrink it to the right size.
            matchedTokens.length = matchedTokensIndex;
        }
        return {
            tokens: matchedTokens,
            groups: groups,
            errors: errors,
        };
    }
    handleModes(config, pop_mode, push_mode, newToken) {
        if (config.pop === true) {
            // need to save the PUSH_MODE property as if the mode is popped
            // patternIdxToPopMode is updated to reflect the new mode after popping the stack
            const pushMode = config.push;
            pop_mode(newToken);
            if (pushMode !== undefined) {
                push_mode.call(this, pushMode);
            }
        }
        else if (config.push !== undefined) {
            push_mode.call(this, config.push);
        }
    }
    chopInput(text, length) {
        return text.substring(length);
    }
    updateLastIndex(regExp, newLastIndex) {
        regExp.lastIndex = newLastIndex;
    }
    // TODO: decrease this under 600 characters? inspect stripping comments option in TSC compiler
    updateTokenEndLineColumnLocation(newToken, group, lastLTIdx, numOfLTsInMatch, line, column, imageLength) {
        let lastCharIsLT, fixForEndingInLT;
        if (group !== undefined) {
            // a none skipped multi line Token, need to update endLine/endColumn
            lastCharIsLT = lastLTIdx === imageLength - 1;
            fixForEndingInLT = lastCharIsLT ? -1 : 0;
            if (!(numOfLTsInMatch === 1 && lastCharIsLT === true)) {
                // if a token ends in a LT that last LT only affects the line numbering of following Tokens
                newToken.endLine = line + fixForEndingInLT;
                // the last LT in a token does not affect the endColumn either as the [columnStart ... columnEnd)
                // inclusive to exclusive range.
                newToken.endColumn = column - 1 + -fixForEndingInLT;
            }
            // else single LT in the last character of a token, no need to modify the endLine/EndColumn
        }
    }
    computeNewColumn(oldColumn, imageLength) {
        return oldColumn + imageLength;
    }
    createOffsetOnlyToken(image, startOffset, tokenTypeIdx, tokenType) {
        return {
            image,
            startOffset,
            tokenTypeIdx,
            tokenType,
        };
    }
    createStartOnlyToken(image, startOffset, tokenTypeIdx, tokenType, startLine, startColumn) {
        return {
            image,
            startOffset,
            startLine,
            startColumn,
            tokenTypeIdx,
            tokenType,
        };
    }
    createFullToken(image, startOffset, tokenTypeIdx, tokenType, startLine, startColumn, imageLength) {
        return {
            image,
            startOffset,
            endOffset: startOffset + imageLength - 1,
            startLine,
            endLine: startLine,
            startColumn,
            endColumn: startColumn + imageLength - 1,
            tokenTypeIdx,
            tokenType,
        };
    }
    addTokenUsingPush(tokenVector, index, tokenToAdd) {
        tokenVector.push(tokenToAdd);
        return index;
    }
    addTokenUsingMemberAccess(tokenVector, index, tokenToAdd) {
        tokenVector[index] = tokenToAdd;
        index++;
        return index;
    }
    handlePayloadNoCustom(token, payload) { }
    handlePayloadWithCustom(token, payload) {
        if (payload !== null) {
            token.payload = payload;
        }
    }
    matchWithTest(pattern, text, offset) {
        const found = pattern.test(text);
        if (found === true) {
            return text.substring(offset, pattern.lastIndex);
        }
        return null;
    }
    matchWithExec(pattern, text) {
        const regExpArray = pattern.exec(text);
        return regExpArray !== null ? regExpArray[0] : null;
    }
}
Lexer.SKIPPED = "This marks a skipped Token pattern, this means each token identified by it will" +
    "be consumed and then thrown into oblivion, this can be used to for example to completely ignore whitespace.";
Lexer.NA = /NOT_APPLICABLE/;
//# sourceMappingURL=lexer_public.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/scan/tokens_public.js



function tokens_public_tokenLabel(tokType) {
    if (tokens_public_hasTokenLabel(tokType)) {
        return tokType.LABEL;
    }
    else {
        return tokType.name;
    }
}
function tokenName(tokType) {
    return tokType.name;
}
function tokens_public_hasTokenLabel(obj) {
    return (0,isString/* default */.Z)(obj.LABEL) && obj.LABEL !== "";
}
const PARENT = "parent";
const CATEGORIES = "categories";
const LABEL = "label";
const GROUP = "group";
const PUSH_MODE = "push_mode";
const POP_MODE = "pop_mode";
const LONGER_ALT = "longer_alt";
const LINE_BREAKS = "line_breaks";
const START_CHARS_HINT = "start_chars_hint";
function createToken(config) {
    return createTokenInternal(config);
}
function createTokenInternal(config) {
    const pattern = config.pattern;
    const tokenType = {};
    tokenType.name = config.name;
    if (!(0,isUndefined/* default */.Z)(pattern)) {
        tokenType.PATTERN = pattern;
    }
    if ((0,has/* default */.Z)(config, PARENT)) {
        throw ("The parent property is no longer supported.\n" +
            "See: https://github.com/chevrotain/chevrotain/issues/564#issuecomment-349062346 for details.");
    }
    if ((0,has/* default */.Z)(config, CATEGORIES)) {
        // casting to ANY as this will be fixed inside `augmentTokenTypes``
        tokenType.CATEGORIES = config[CATEGORIES];
    }
    augmentTokenTypes([tokenType]);
    if ((0,has/* default */.Z)(config, LABEL)) {
        tokenType.LABEL = config[LABEL];
    }
    if ((0,has/* default */.Z)(config, GROUP)) {
        tokenType.GROUP = config[GROUP];
    }
    if ((0,has/* default */.Z)(config, POP_MODE)) {
        tokenType.POP_MODE = config[POP_MODE];
    }
    if ((0,has/* default */.Z)(config, PUSH_MODE)) {
        tokenType.PUSH_MODE = config[PUSH_MODE];
    }
    if ((0,has/* default */.Z)(config, LONGER_ALT)) {
        tokenType.LONGER_ALT = config[LONGER_ALT];
    }
    if ((0,has/* default */.Z)(config, LINE_BREAKS)) {
        tokenType.LINE_BREAKS = config[LINE_BREAKS];
    }
    if ((0,has/* default */.Z)(config, START_CHARS_HINT)) {
        tokenType.START_CHARS_HINT = config[START_CHARS_HINT];
    }
    return tokenType;
}
const EOF = createToken({ name: "EOF", pattern: Lexer.NA });
augmentTokenTypes([EOF]);
function createTokenInstance(tokType, image, startOffset, endOffset, startLine, endLine, startColumn, endColumn) {
    return {
        image,
        startOffset,
        endOffset,
        startLine,
        endLine,
        startColumn,
        endColumn,
        tokenTypeIdx: tokType.tokenTypeIdx,
        tokenType: tokType,
    };
}
function tokenMatcher(token, tokType) {
    return tokenStructuredMatcher(token, tokType);
}
//# sourceMappingURL=tokens_public.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/errors_public.js



const defaultParserErrorProvider = {
    buildMismatchTokenMessage({ expected, actual, previous, ruleName }) {
        const hasLabel = tokens_public_hasTokenLabel(expected);
        const expectedMsg = hasLabel
            ? `--> ${tokens_public_tokenLabel(expected)} <--`
            : `token of type --> ${expected.name} <--`;
        const msg = `Expecting ${expectedMsg} but found --> '${actual.image}' <--`;
        return msg;
    },
    buildNotAllInputParsedMessage({ firstRedundant, ruleName }) {
        return "Redundant input, expecting EOF but found: " + firstRedundant.image;
    },
    buildNoViableAltMessage({ expectedPathsPerAlt, actual, previous, customUserDescription, ruleName, }) {
        const errPrefix = "Expecting: ";
        // TODO: issue: No Viable Alternative Error may have incomplete details. #502
        const actualText = lodash_es_head(actual).image;
        const errSuffix = "\nbut found: '" + actualText + "'";
        if (customUserDescription) {
            return errPrefix + customUserDescription + errSuffix;
        }
        else {
            const allLookAheadPaths = (0,reduce/* default */.Z)(expectedPathsPerAlt, (result, currAltPaths) => result.concat(currAltPaths), []);
            const nextValidTokenSequences = (0,map/* default */.Z)(allLookAheadPaths, (currPath) => `[${(0,map/* default */.Z)(currPath, (currTokenType) => tokens_public_tokenLabel(currTokenType)).join(", ")}]`);
            const nextValidSequenceItems = (0,map/* default */.Z)(nextValidTokenSequences, (itemMsg, idx) => `  ${idx + 1}. ${itemMsg}`);
            const calculatedDescription = `one of these possible Token sequences:\n${nextValidSequenceItems.join("\n")}`;
            return errPrefix + calculatedDescription + errSuffix;
        }
    },
    buildEarlyExitMessage({ expectedIterationPaths, actual, customUserDescription, ruleName, }) {
        const errPrefix = "Expecting: ";
        // TODO: issue: No Viable Alternative Error may have incomplete details. #502
        const actualText = lodash_es_head(actual).image;
        const errSuffix = "\nbut found: '" + actualText + "'";
        if (customUserDescription) {
            return errPrefix + customUserDescription + errSuffix;
        }
        else {
            const nextValidTokenSequences = (0,map/* default */.Z)(expectedIterationPaths, (currPath) => `[${(0,map/* default */.Z)(currPath, (currTokenType) => tokens_public_tokenLabel(currTokenType)).join(",")}]`);
            const calculatedDescription = `expecting at least one iteration which starts with one of these possible Token sequences::\n  ` +
                `<${nextValidTokenSequences.join(" ,")}>`;
            return errPrefix + calculatedDescription + errSuffix;
        }
    },
};
Object.freeze(defaultParserErrorProvider);
const defaultGrammarResolverErrorProvider = {
    buildRuleNotFoundError(topLevelRule, undefinedRule) {
        const msg = "Invalid grammar, reference to a rule which is not defined: ->" +
            undefinedRule.nonTerminalName +
            "<-\n" +
            "inside top level rule: ->" +
            topLevelRule.name +
            "<-";
        return msg;
    },
};
const defaultGrammarValidatorErrorProvider = {
    buildDuplicateFoundError(topLevelRule, duplicateProds) {
        function getExtraProductionArgument(prod) {
            if (prod instanceof Terminal) {
                return prod.terminalType.name;
            }
            else if (prod instanceof NonTerminal) {
                return prod.nonTerminalName;
            }
            else {
                return "";
            }
        }
        const topLevelName = topLevelRule.name;
        const duplicateProd = lodash_es_head(duplicateProds);
        const index = duplicateProd.idx;
        const dslName = getProductionDslName(duplicateProd);
        const extraArgument = getExtraProductionArgument(duplicateProd);
        const hasExplicitIndex = index > 0;
        let msg = `->${dslName}${hasExplicitIndex ? index : ""}<- ${extraArgument ? `with argument: ->${extraArgument}<-` : ""}
                  appears more than once (${duplicateProds.length} times) in the top level rule: ->${topLevelName}<-.                  
                  For further details see: https://chevrotain.io/docs/FAQ.html#NUMERICAL_SUFFIXES 
                  `;
        // white space trimming time! better to trim afterwards as it allows to use WELL formatted multi line template strings...
        msg = msg.replace(/[ \t]+/g, " ");
        msg = msg.replace(/\s\s+/g, "\n");
        return msg;
    },
    buildNamespaceConflictError(rule) {
        const errMsg = `Namespace conflict found in grammar.\n` +
            `The grammar has both a Terminal(Token) and a Non-Terminal(Rule) named: <${rule.name}>.\n` +
            `To resolve this make sure each Terminal and Non-Terminal names are unique\n` +
            `This is easy to accomplish by using the convention that Terminal names start with an uppercase letter\n` +
            `and Non-Terminal names start with a lower case letter.`;
        return errMsg;
    },
    buildAlternationPrefixAmbiguityError(options) {
        const pathMsg = (0,map/* default */.Z)(options.prefixPath, (currTok) => tokens_public_tokenLabel(currTok)).join(", ");
        const occurrence = options.alternation.idx === 0 ? "" : options.alternation.idx;
        const errMsg = `Ambiguous alternatives: <${options.ambiguityIndices.join(" ,")}> due to common lookahead prefix\n` +
            `in <OR${occurrence}> inside <${options.topLevelRule.name}> Rule,\n` +
            `<${pathMsg}> may appears as a prefix path in all these alternatives.\n` +
            `See: https://chevrotain.io/docs/guide/resolving_grammar_errors.html#COMMON_PREFIX\n` +
            `For Further details.`;
        return errMsg;
    },
    buildAlternationAmbiguityError(options) {
        const pathMsg = (0,map/* default */.Z)(options.prefixPath, (currtok) => tokens_public_tokenLabel(currtok)).join(", ");
        const occurrence = options.alternation.idx === 0 ? "" : options.alternation.idx;
        let currMessage = `Ambiguous Alternatives Detected: <${options.ambiguityIndices.join(" ,")}> in <OR${occurrence}>` +
            ` inside <${options.topLevelRule.name}> Rule,\n` +
            `<${pathMsg}> may appears as a prefix path in all these alternatives.\n`;
        currMessage =
            currMessage +
                `See: https://chevrotain.io/docs/guide/resolving_grammar_errors.html#AMBIGUOUS_ALTERNATIVES\n` +
                `For Further details.`;
        return currMessage;
    },
    buildEmptyRepetitionError(options) {
        let dslName = getProductionDslName(options.repetition);
        if (options.repetition.idx !== 0) {
            dslName += options.repetition.idx;
        }
        const errMsg = `The repetition <${dslName}> within Rule <${options.topLevelRule.name}> can never consume any tokens.\n` +
            `This could lead to an infinite loop.`;
        return errMsg;
    },
    // TODO: remove - `errors_public` from nyc.config.js exclude
    //       once this method is fully removed from this file
    buildTokenNameError(options) {
        /* istanbul ignore next */
        return "deprecated";
    },
    buildEmptyAlternationError(options) {
        const errMsg = `Ambiguous empty alternative: <${options.emptyChoiceIdx + 1}>` +
            ` in <OR${options.alternation.idx}> inside <${options.topLevelRule.name}> Rule.\n` +
            `Only the last alternative may be an empty alternative.`;
        return errMsg;
    },
    buildTooManyAlternativesError(options) {
        const errMsg = `An Alternation cannot have more than 256 alternatives:\n` +
            `<OR${options.alternation.idx}> inside <${options.topLevelRule.name}> Rule.\n has ${options.alternation.definition.length + 1} alternatives.`;
        return errMsg;
    },
    buildLeftRecursionError(options) {
        const ruleName = options.topLevelRule.name;
        const pathNames = (0,map/* default */.Z)(options.leftRecursionPath, (currRule) => currRule.name);
        const leftRecursivePath = `${ruleName} --> ${pathNames
            .concat([ruleName])
            .join(" --> ")}`;
        const errMsg = `Left Recursion found in grammar.\n` +
            `rule: <${ruleName}> can be invoked from itself (directly or indirectly)\n` +
            `without consuming any Tokens. The grammar path that causes this is: \n ${leftRecursivePath}\n` +
            ` To fix this refactor your grammar to remove the left recursion.\n` +
            `see: https://en.wikipedia.org/wiki/LL_parser#Left_factoring.`;
        return errMsg;
    },
    // TODO: remove - `errors_public` from nyc.config.js exclude
    //       once this method is fully removed from this file
    buildInvalidRuleNameError(options) {
        /* istanbul ignore next */
        return "deprecated";
    },
    buildDuplicateRuleNameError(options) {
        let ruleName;
        if (options.topLevelRule instanceof Rule) {
            ruleName = options.topLevelRule.name;
        }
        else {
            ruleName = options.topLevelRule;
        }
        const errMsg = `Duplicate definition, rule: ->${ruleName}<- is already defined in the grammar: ->${options.grammarName}<-`;
        return errMsg;
    },
};
//# sourceMappingURL=errors_public.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/gast/lib/src/visitor.js

class GAstVisitor {
    visit(node) {
        const nodeAny = node;
        switch (nodeAny.constructor) {
            case NonTerminal:
                return this.visitNonTerminal(nodeAny);
            case Alternative:
                return this.visitAlternative(nodeAny);
            case Option:
                return this.visitOption(nodeAny);
            case RepetitionMandatory:
                return this.visitRepetitionMandatory(nodeAny);
            case RepetitionMandatoryWithSeparator:
                return this.visitRepetitionMandatoryWithSeparator(nodeAny);
            case RepetitionWithSeparator:
                return this.visitRepetitionWithSeparator(nodeAny);
            case Repetition:
                return this.visitRepetition(nodeAny);
            case Alternation:
                return this.visitAlternation(nodeAny);
            case Terminal:
                return this.visitTerminal(nodeAny);
            case Rule:
                return this.visitRule(nodeAny);
            /* c8 ignore next 2 */
            default:
                throw Error("non exhaustive match");
        }
    }
    /* c8 ignore next */
    visitNonTerminal(node) { }
    /* c8 ignore next */
    visitAlternative(node) { }
    /* c8 ignore next */
    visitOption(node) { }
    /* c8 ignore next */
    visitRepetition(node) { }
    /* c8 ignore next */
    visitRepetitionMandatory(node) { }
    /* c8 ignore next 3 */
    visitRepetitionMandatoryWithSeparator(node) { }
    /* c8 ignore next */
    visitRepetitionWithSeparator(node) { }
    /* c8 ignore next */
    visitAlternation(node) { }
    /* c8 ignore next */
    visitTerminal(node) { }
    /* c8 ignore next */
    visitRule(node) { }
}
//# sourceMappingURL=visitor.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/resolver.js



function resolveGrammar(topLevels, errMsgProvider) {
    const refResolver = new GastRefResolverVisitor(topLevels, errMsgProvider);
    refResolver.resolveRefs();
    return refResolver.errors;
}
class GastRefResolverVisitor extends GAstVisitor {
    constructor(nameToTopRule, errMsgProvider) {
        super();
        this.nameToTopRule = nameToTopRule;
        this.errMsgProvider = errMsgProvider;
        this.errors = [];
    }
    resolveRefs() {
        (0,forEach/* default */.Z)((0,values/* default */.Z)(this.nameToTopRule), (prod) => {
            this.currTopLevel = prod;
            prod.accept(this);
        });
    }
    visitNonTerminal(node) {
        const ref = this.nameToTopRule[node.nonTerminalName];
        if (!ref) {
            const msg = this.errMsgProvider.buildRuleNotFoundError(this.currTopLevel, node);
            this.errors.push({
                message: msg,
                type: ParserDefinitionErrorType.UNRESOLVED_SUBRULE_REF,
                ruleName: this.currTopLevel.name,
                unresolvedRefName: node.nonTerminalName,
            });
        }
        else {
            node.referencedRule = ref;
        }
    }
}
//# sourceMappingURL=resolver.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/flatMap.js
var flatMap = __webpack_require__(34134);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseAssignValue.js
var _baseAssignValue = __webpack_require__(93586);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_arrayAggregator.js
/**
 * A specialized version of `baseAggregator` for arrays.
 *
 * @private
 * @param {Array} [array] The array to iterate over.
 * @param {Function} setter The function to set `accumulator` values.
 * @param {Function} iteratee The iteratee to transform keys.
 * @param {Object} accumulator The initial aggregated object.
 * @returns {Function} Returns `accumulator`.
 */
function arrayAggregator(array, setter, iteratee, accumulator) {
  var index = -1,
      length = array == null ? 0 : array.length;

  while (++index < length) {
    var value = array[index];
    setter(accumulator, value, iteratee(value), array);
  }
  return accumulator;
}

/* harmony default export */ const _arrayAggregator = (arrayAggregator);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseAggregator.js


/**
 * Aggregates elements of `collection` on `accumulator` with keys transformed
 * by `iteratee` and values set by `setter`.
 *
 * @private
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} setter The function to set `accumulator` values.
 * @param {Function} iteratee The iteratee to transform keys.
 * @param {Object} accumulator The initial aggregated object.
 * @returns {Function} Returns `accumulator`.
 */
function baseAggregator(collection, setter, iteratee, accumulator) {
  (0,_baseEach/* default */.Z)(collection, function(value, key, collection) {
    setter(accumulator, value, iteratee(value), collection);
  });
  return accumulator;
}

/* harmony default export */ const _baseAggregator = (baseAggregator);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_createAggregator.js





/**
 * Creates a function like `_.groupBy`.
 *
 * @private
 * @param {Function} setter The function to set accumulator values.
 * @param {Function} [initializer] The accumulator object initializer.
 * @returns {Function} Returns the new aggregator function.
 */
function createAggregator(setter, initializer) {
  return function(collection, iteratee) {
    var func = (0,isArray/* default */.Z)(collection) ? _arrayAggregator : _baseAggregator,
        accumulator = initializer ? initializer() : {};

    return func(collection, setter, (0,_baseIteratee/* default */.Z)(iteratee, 2), accumulator);
  };
}

/* harmony default export */ const _createAggregator = (createAggregator);

;// CONCATENATED MODULE: ../node_modules/lodash-es/groupBy.js



/** Used for built-in method references. */
var groupBy_objectProto = Object.prototype;

/** Used to check objects for own properties. */
var groupBy_hasOwnProperty = groupBy_objectProto.hasOwnProperty;

/**
 * Creates an object composed of keys generated from the results of running
 * each element of `collection` thru `iteratee`. The order of grouped values
 * is determined by the order they occur in `collection`. The corresponding
 * value of each key is an array of elements responsible for generating the
 * key. The iteratee is invoked with one argument: (value).
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [iteratee=_.identity] The iteratee to transform keys.
 * @returns {Object} Returns the composed aggregate object.
 * @example
 *
 * _.groupBy([6.1, 4.2, 6.3], Math.floor);
 * // => { '4': [4.2], '6': [6.1, 6.3] }
 *
 * // The `_.property` iteratee shorthand.
 * _.groupBy(['one', 'two', 'three'], 'length');
 * // => { '3': ['one', 'two'], '5': ['three'] }
 */
var groupBy = _createAggregator(function(result, value, key) {
  if (groupBy_hasOwnProperty.call(result, key)) {
    result[key].push(value);
  } else {
    (0,_baseAssignValue/* default */.Z)(result, key, [value]);
  }
});

/* harmony default export */ const lodash_es_groupBy = (groupBy);

;// CONCATENATED MODULE: ../node_modules/lodash-es/dropRight.js



/**
 * Creates a slice of `array` with `n` elements dropped from the end.
 *
 * @static
 * @memberOf _
 * @since 3.0.0
 * @category Array
 * @param {Array} array The array to query.
 * @param {number} [n=1] The number of elements to drop.
 * @param- {Object} [guard] Enables use as an iteratee for methods like `_.map`.
 * @returns {Array} Returns the slice of `array`.
 * @example
 *
 * _.dropRight([1, 2, 3]);
 * // => [1, 2]
 *
 * _.dropRight([1, 2, 3], 2);
 * // => [1]
 *
 * _.dropRight([1, 2, 3], 5);
 * // => []
 *
 * _.dropRight([1, 2, 3], 0);
 * // => [1, 2, 3]
 */
function dropRight(array, n, guard) {
  var length = array == null ? 0 : array.length;
  if (!length) {
    return [];
  }
  n = (guard || n === undefined) ? 1 : (0,toInteger/* default */.Z)(n);
  n = length - n;
  return _baseSlice(array, 0, n < 0 ? 0 : n);
}

/* harmony default export */ const lodash_es_dropRight = (dropRight);

;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/interpreter.js




class AbstractNextPossibleTokensWalker extends RestWalker {
    constructor(topProd, path) {
        super();
        this.topProd = topProd;
        this.path = path;
        this.possibleTokTypes = [];
        this.nextProductionName = "";
        this.nextProductionOccurrence = 0;
        this.found = false;
        this.isAtEndOfPath = false;
    }
    startWalking() {
        this.found = false;
        if (this.path.ruleStack[0] !== this.topProd.name) {
            throw Error("The path does not start with the walker's top Rule!");
        }
        // immutable for the win
        this.ruleStack = (0,lodash_es_clone/* default */.Z)(this.path.ruleStack).reverse(); // intelij bug requires assertion
        this.occurrenceStack = (0,lodash_es_clone/* default */.Z)(this.path.occurrenceStack).reverse(); // intelij bug requires assertion
        // already verified that the first production is valid, we now seek the 2nd production
        this.ruleStack.pop();
        this.occurrenceStack.pop();
        this.updateExpectedNext();
        this.walk(this.topProd);
        return this.possibleTokTypes;
    }
    walk(prod, prevRest = []) {
        // stop scanning once we found the path
        if (!this.found) {
            super.walk(prod, prevRest);
        }
    }
    walkProdRef(refProd, currRest, prevRest) {
        // found the next production, need to keep walking in it
        if (refProd.referencedRule.name === this.nextProductionName &&
            refProd.idx === this.nextProductionOccurrence) {
            const fullRest = currRest.concat(prevRest);
            this.updateExpectedNext();
            this.walk(refProd.referencedRule, fullRest);
        }
    }
    updateExpectedNext() {
        // need to consume the Terminal
        if ((0,isEmpty/* default */.Z)(this.ruleStack)) {
            // must reset nextProductionXXX to avoid walking down another Top Level production while what we are
            // really seeking is the last Terminal...
            this.nextProductionName = "";
            this.nextProductionOccurrence = 0;
            this.isAtEndOfPath = true;
        }
        else {
            this.nextProductionName = this.ruleStack.pop();
            this.nextProductionOccurrence = this.occurrenceStack.pop();
        }
    }
}
class NextAfterTokenWalker extends AbstractNextPossibleTokensWalker {
    constructor(topProd, path) {
        super(topProd, path);
        this.path = path;
        this.nextTerminalName = "";
        this.nextTerminalOccurrence = 0;
        this.nextTerminalName = this.path.lastTok.name;
        this.nextTerminalOccurrence = this.path.lastTokOccurrence;
    }
    walkTerminal(terminal, currRest, prevRest) {
        if (this.isAtEndOfPath &&
            terminal.terminalType.name === this.nextTerminalName &&
            terminal.idx === this.nextTerminalOccurrence &&
            !this.found) {
            const fullRest = currRest.concat(prevRest);
            const restProd = new Alternative({ definition: fullRest });
            this.possibleTokTypes = first(restProd);
            this.found = true;
        }
    }
}
/**
 * This walker only "walks" a single "TOP" level in the Grammar Ast, this means
 * it never "follows" production refs
 */
class AbstractNextTerminalAfterProductionWalker extends RestWalker {
    constructor(topRule, occurrence) {
        super();
        this.topRule = topRule;
        this.occurrence = occurrence;
        this.result = {
            token: undefined,
            occurrence: undefined,
            isEndOfRule: undefined,
        };
    }
    startWalking() {
        this.walk(this.topRule);
        return this.result;
    }
}
class NextTerminalAfterManyWalker extends AbstractNextTerminalAfterProductionWalker {
    walkMany(manyProd, currRest, prevRest) {
        if (manyProd.idx === this.occurrence) {
            const firstAfterMany = lodash_es_head(currRest.concat(prevRest));
            this.result.isEndOfRule = firstAfterMany === undefined;
            if (firstAfterMany instanceof Terminal) {
                this.result.token = firstAfterMany.terminalType;
                this.result.occurrence = firstAfterMany.idx;
            }
        }
        else {
            super.walkMany(manyProd, currRest, prevRest);
        }
    }
}
class NextTerminalAfterManySepWalker extends AbstractNextTerminalAfterProductionWalker {
    walkManySep(manySepProd, currRest, prevRest) {
        if (manySepProd.idx === this.occurrence) {
            const firstAfterManySep = lodash_es_head(currRest.concat(prevRest));
            this.result.isEndOfRule = firstAfterManySep === undefined;
            if (firstAfterManySep instanceof Terminal) {
                this.result.token = firstAfterManySep.terminalType;
                this.result.occurrence = firstAfterManySep.idx;
            }
        }
        else {
            super.walkManySep(manySepProd, currRest, prevRest);
        }
    }
}
class NextTerminalAfterAtLeastOneWalker extends AbstractNextTerminalAfterProductionWalker {
    walkAtLeastOne(atLeastOneProd, currRest, prevRest) {
        if (atLeastOneProd.idx === this.occurrence) {
            const firstAfterAtLeastOne = lodash_es_head(currRest.concat(prevRest));
            this.result.isEndOfRule = firstAfterAtLeastOne === undefined;
            if (firstAfterAtLeastOne instanceof Terminal) {
                this.result.token = firstAfterAtLeastOne.terminalType;
                this.result.occurrence = firstAfterAtLeastOne.idx;
            }
        }
        else {
            super.walkAtLeastOne(atLeastOneProd, currRest, prevRest);
        }
    }
}
// TODO: reduce code duplication in the AfterWalkers
class NextTerminalAfterAtLeastOneSepWalker extends AbstractNextTerminalAfterProductionWalker {
    walkAtLeastOneSep(atleastOneSepProd, currRest, prevRest) {
        if (atleastOneSepProd.idx === this.occurrence) {
            const firstAfterfirstAfterAtLeastOneSep = lodash_es_head(currRest.concat(prevRest));
            this.result.isEndOfRule = firstAfterfirstAfterAtLeastOneSep === undefined;
            if (firstAfterfirstAfterAtLeastOneSep instanceof Terminal) {
                this.result.token = firstAfterfirstAfterAtLeastOneSep.terminalType;
                this.result.occurrence = firstAfterfirstAfterAtLeastOneSep.idx;
            }
        }
        else {
            super.walkAtLeastOneSep(atleastOneSepProd, currRest, prevRest);
        }
    }
}
function possiblePathsFrom(targetDef, maxLength, currPath = []) {
    // avoid side effects
    currPath = (0,lodash_es_clone/* default */.Z)(currPath);
    let result = [];
    let i = 0;
    // TODO: avoid inner funcs
    function remainingPathWith(nextDef) {
        return nextDef.concat(lodash_es_drop(targetDef, i + 1));
    }
    // TODO: avoid inner funcs
    function getAlternativesForProd(definition) {
        const alternatives = possiblePathsFrom(remainingPathWith(definition), maxLength, currPath);
        return result.concat(alternatives);
    }
    /**
     * Mandatory productions will halt the loop as the paths computed from their recursive calls will already contain the
     * following (rest) of the targetDef.
     *
     * For optional productions (Option/Repetition/...) the loop will continue to represent the paths that do not include the
     * the optional production.
     */
    while (currPath.length < maxLength && i < targetDef.length) {
        const prod = targetDef[i];
        /* istanbul ignore else */
        if (prod instanceof Alternative) {
            return getAlternativesForProd(prod.definition);
        }
        else if (prod instanceof NonTerminal) {
            return getAlternativesForProd(prod.definition);
        }
        else if (prod instanceof Option) {
            result = getAlternativesForProd(prod.definition);
        }
        else if (prod instanceof RepetitionMandatory) {
            const newDef = prod.definition.concat([
                new Repetition({
                    definition: prod.definition,
                }),
            ]);
            return getAlternativesForProd(newDef);
        }
        else if (prod instanceof RepetitionMandatoryWithSeparator) {
            const newDef = [
                new Alternative({ definition: prod.definition }),
                new Repetition({
                    definition: [new Terminal({ terminalType: prod.separator })].concat(prod.definition),
                }),
            ];
            return getAlternativesForProd(newDef);
        }
        else if (prod instanceof RepetitionWithSeparator) {
            const newDef = prod.definition.concat([
                new Repetition({
                    definition: [new Terminal({ terminalType: prod.separator })].concat(prod.definition),
                }),
            ]);
            result = getAlternativesForProd(newDef);
        }
        else if (prod instanceof Repetition) {
            const newDef = prod.definition.concat([
                new Repetition({
                    definition: prod.definition,
                }),
            ]);
            result = getAlternativesForProd(newDef);
        }
        else if (prod instanceof Alternation) {
            (0,forEach/* default */.Z)(prod.definition, (currAlt) => {
                // TODO: this is a limited check for empty alternatives
                //   It would prevent a common case of infinite loops during parser initialization.
                //   However **in-directly** empty alternatives may still cause issues.
                if ((0,isEmpty/* default */.Z)(currAlt.definition) === false) {
                    result = getAlternativesForProd(currAlt.definition);
                }
            });
            return result;
        }
        else if (prod instanceof Terminal) {
            currPath.push(prod.terminalType);
        }
        else {
            throw Error("non exhaustive match");
        }
        i++;
    }
    result.push({
        partialPath: currPath,
        suffixDef: lodash_es_drop(targetDef, i),
    });
    return result;
}
function nextPossibleTokensAfter(initialDef, tokenVector, tokMatcher, maxLookAhead) {
    const EXIT_NON_TERMINAL = "EXIT_NONE_TERMINAL";
    // to avoid creating a new Array each time.
    const EXIT_NON_TERMINAL_ARR = [EXIT_NON_TERMINAL];
    const EXIT_ALTERNATIVE = "EXIT_ALTERNATIVE";
    let foundCompletePath = false;
    const tokenVectorLength = tokenVector.length;
    const minimalAlternativesIndex = tokenVectorLength - maxLookAhead - 1;
    const result = [];
    const possiblePaths = [];
    possiblePaths.push({
        idx: -1,
        def: initialDef,
        ruleStack: [],
        occurrenceStack: [],
    });
    while (!(0,isEmpty/* default */.Z)(possiblePaths)) {
        const currPath = possiblePaths.pop();
        // skip alternatives if no more results can be found (assuming deterministic grammar with fixed lookahead)
        if (currPath === EXIT_ALTERNATIVE) {
            if (foundCompletePath &&
                (0,last/* default */.Z)(possiblePaths).idx <= minimalAlternativesIndex) {
                // remove irrelevant alternative
                possiblePaths.pop();
            }
            continue;
        }
        const currDef = currPath.def;
        const currIdx = currPath.idx;
        const currRuleStack = currPath.ruleStack;
        const currOccurrenceStack = currPath.occurrenceStack;
        // For Example: an empty path could exist in a valid grammar in the case of an EMPTY_ALT
        if ((0,isEmpty/* default */.Z)(currDef)) {
            continue;
        }
        const prod = currDef[0];
        /* istanbul ignore else */
        if (prod === EXIT_NON_TERMINAL) {
            const nextPath = {
                idx: currIdx,
                def: lodash_es_drop(currDef),
                ruleStack: lodash_es_dropRight(currRuleStack),
                occurrenceStack: lodash_es_dropRight(currOccurrenceStack),
            };
            possiblePaths.push(nextPath);
        }
        else if (prod instanceof Terminal) {
            /* istanbul ignore else */
            if (currIdx < tokenVectorLength - 1) {
                const nextIdx = currIdx + 1;
                const actualToken = tokenVector[nextIdx];
                if (tokMatcher(actualToken, prod.terminalType)) {
                    const nextPath = {
                        idx: nextIdx,
                        def: lodash_es_drop(currDef),
                        ruleStack: currRuleStack,
                        occurrenceStack: currOccurrenceStack,
                    };
                    possiblePaths.push(nextPath);
                }
                // end of the line
            }
            else if (currIdx === tokenVectorLength - 1) {
                // IGNORE ABOVE ELSE
                result.push({
                    nextTokenType: prod.terminalType,
                    nextTokenOccurrence: prod.idx,
                    ruleStack: currRuleStack,
                    occurrenceStack: currOccurrenceStack,
                });
                foundCompletePath = true;
            }
            else {
                throw Error("non exhaustive match");
            }
        }
        else if (prod instanceof NonTerminal) {
            const newRuleStack = (0,lodash_es_clone/* default */.Z)(currRuleStack);
            newRuleStack.push(prod.nonTerminalName);
            const newOccurrenceStack = (0,lodash_es_clone/* default */.Z)(currOccurrenceStack);
            newOccurrenceStack.push(prod.idx);
            const nextPath = {
                idx: currIdx,
                def: prod.definition.concat(EXIT_NON_TERMINAL_ARR, lodash_es_drop(currDef)),
                ruleStack: newRuleStack,
                occurrenceStack: newOccurrenceStack,
            };
            possiblePaths.push(nextPath);
        }
        else if (prod instanceof Option) {
            // the order of alternatives is meaningful, FILO (Last path will be traversed first).
            const nextPathWithout = {
                idx: currIdx,
                def: lodash_es_drop(currDef),
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWithout);
            // required marker to avoid backtracking paths whose higher priority alternatives already matched
            possiblePaths.push(EXIT_ALTERNATIVE);
            const nextPathWith = {
                idx: currIdx,
                def: prod.definition.concat(lodash_es_drop(currDef)),
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWith);
        }
        else if (prod instanceof RepetitionMandatory) {
            // TODO:(THE NEW operators here take a while...) (convert once?)
            const secondIteration = new Repetition({
                definition: prod.definition,
                idx: prod.idx,
            });
            const nextDef = prod.definition.concat([secondIteration], lodash_es_drop(currDef));
            const nextPath = {
                idx: currIdx,
                def: nextDef,
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPath);
        }
        else if (prod instanceof RepetitionMandatoryWithSeparator) {
            // TODO:(THE NEW operators here take a while...) (convert once?)
            const separatorGast = new Terminal({
                terminalType: prod.separator,
            });
            const secondIteration = new Repetition({
                definition: [separatorGast].concat(prod.definition),
                idx: prod.idx,
            });
            const nextDef = prod.definition.concat([secondIteration], lodash_es_drop(currDef));
            const nextPath = {
                idx: currIdx,
                def: nextDef,
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPath);
        }
        else if (prod instanceof RepetitionWithSeparator) {
            // the order of alternatives is meaningful, FILO (Last path will be traversed first).
            const nextPathWithout = {
                idx: currIdx,
                def: lodash_es_drop(currDef),
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWithout);
            // required marker to avoid backtracking paths whose higher priority alternatives already matched
            possiblePaths.push(EXIT_ALTERNATIVE);
            const separatorGast = new Terminal({
                terminalType: prod.separator,
            });
            const nthRepetition = new Repetition({
                definition: [separatorGast].concat(prod.definition),
                idx: prod.idx,
            });
            const nextDef = prod.definition.concat([nthRepetition], lodash_es_drop(currDef));
            const nextPathWith = {
                idx: currIdx,
                def: nextDef,
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWith);
        }
        else if (prod instanceof Repetition) {
            // the order of alternatives is meaningful, FILO (Last path will be traversed first).
            const nextPathWithout = {
                idx: currIdx,
                def: lodash_es_drop(currDef),
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWithout);
            // required marker to avoid backtracking paths whose higher priority alternatives already matched
            possiblePaths.push(EXIT_ALTERNATIVE);
            // TODO: an empty repetition will cause infinite loops here, will the parser detect this in selfAnalysis?
            const nthRepetition = new Repetition({
                definition: prod.definition,
                idx: prod.idx,
            });
            const nextDef = prod.definition.concat([nthRepetition], lodash_es_drop(currDef));
            const nextPathWith = {
                idx: currIdx,
                def: nextDef,
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            };
            possiblePaths.push(nextPathWith);
        }
        else if (prod instanceof Alternation) {
            // the order of alternatives is meaningful, FILO (Last path will be traversed first).
            for (let i = prod.definition.length - 1; i >= 0; i--) {
                const currAlt = prod.definition[i];
                const currAltPath = {
                    idx: currIdx,
                    def: currAlt.definition.concat(lodash_es_drop(currDef)),
                    ruleStack: currRuleStack,
                    occurrenceStack: currOccurrenceStack,
                };
                possiblePaths.push(currAltPath);
                possiblePaths.push(EXIT_ALTERNATIVE);
            }
        }
        else if (prod instanceof Alternative) {
            possiblePaths.push({
                idx: currIdx,
                def: prod.definition.concat(lodash_es_drop(currDef)),
                ruleStack: currRuleStack,
                occurrenceStack: currOccurrenceStack,
            });
        }
        else if (prod instanceof Rule) {
            // last because we should only encounter at most a single one of these per invocation.
            possiblePaths.push(expandTopLevelRule(prod, currIdx, currRuleStack, currOccurrenceStack));
        }
        else {
            throw Error("non exhaustive match");
        }
    }
    return result;
}
function expandTopLevelRule(topRule, currIdx, currRuleStack, currOccurrenceStack) {
    const newRuleStack = (0,lodash_es_clone/* default */.Z)(currRuleStack);
    newRuleStack.push(topRule.name);
    const newCurrOccurrenceStack = (0,lodash_es_clone/* default */.Z)(currOccurrenceStack);
    // top rule is always assumed to have been called with occurrence index 1
    newCurrOccurrenceStack.push(1);
    return {
        idx: currIdx,
        def: topRule.definition,
        ruleStack: newRuleStack,
        occurrenceStack: newCurrOccurrenceStack,
    };
}
//# sourceMappingURL=interpreter.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/lookahead.js





var PROD_TYPE;
(function (PROD_TYPE) {
    PROD_TYPE[PROD_TYPE["OPTION"] = 0] = "OPTION";
    PROD_TYPE[PROD_TYPE["REPETITION"] = 1] = "REPETITION";
    PROD_TYPE[PROD_TYPE["REPETITION_MANDATORY"] = 2] = "REPETITION_MANDATORY";
    PROD_TYPE[PROD_TYPE["REPETITION_MANDATORY_WITH_SEPARATOR"] = 3] = "REPETITION_MANDATORY_WITH_SEPARATOR";
    PROD_TYPE[PROD_TYPE["REPETITION_WITH_SEPARATOR"] = 4] = "REPETITION_WITH_SEPARATOR";
    PROD_TYPE[PROD_TYPE["ALTERNATION"] = 5] = "ALTERNATION";
})(PROD_TYPE || (PROD_TYPE = {}));
function getProdType(prod) {
    /* istanbul ignore else */
    if (prod instanceof Option || prod === "Option") {
        return PROD_TYPE.OPTION;
    }
    else if (prod instanceof Repetition || prod === "Repetition") {
        return PROD_TYPE.REPETITION;
    }
    else if (prod instanceof RepetitionMandatory ||
        prod === "RepetitionMandatory") {
        return PROD_TYPE.REPETITION_MANDATORY;
    }
    else if (prod instanceof RepetitionMandatoryWithSeparator ||
        prod === "RepetitionMandatoryWithSeparator") {
        return PROD_TYPE.REPETITION_MANDATORY_WITH_SEPARATOR;
    }
    else if (prod instanceof RepetitionWithSeparator ||
        prod === "RepetitionWithSeparator") {
        return PROD_TYPE.REPETITION_WITH_SEPARATOR;
    }
    else if (prod instanceof Alternation || prod === "Alternation") {
        return PROD_TYPE.ALTERNATION;
    }
    else {
        throw Error("non exhaustive match");
    }
}
function getLookaheadPaths(options) {
    const { occurrence, rule, prodType, maxLookahead } = options;
    const type = getProdType(prodType);
    if (type === PROD_TYPE.ALTERNATION) {
        return getLookaheadPathsForOr(occurrence, rule, maxLookahead);
    }
    else {
        return getLookaheadPathsForOptionalProd(occurrence, rule, type, maxLookahead);
    }
}
function buildLookaheadFuncForOr(occurrence, ruleGrammar, maxLookahead, hasPredicates, dynamicTokensEnabled, laFuncBuilder) {
    const lookAheadPaths = getLookaheadPathsForOr(occurrence, ruleGrammar, maxLookahead);
    const tokenMatcher = areTokenCategoriesNotUsed(lookAheadPaths)
        ? tokenStructuredMatcherNoCategories
        : tokenStructuredMatcher;
    return laFuncBuilder(lookAheadPaths, hasPredicates, tokenMatcher, dynamicTokensEnabled);
}
/**
 *  When dealing with an Optional production (OPTION/MANY/2nd iteration of AT_LEAST_ONE/...) we need to compare
 *  the lookahead "inside" the production and the lookahead immediately "after" it in the same top level rule (context free).
 *
 *  Example: given a production:
 *  ABC(DE)?DF
 *
 *  The optional '(DE)?' should only be entered if we see 'DE'. a single Token 'D' is not sufficient to distinguish between the two
 *  alternatives.
 *
 *  @returns A Lookahead function which will return true IFF the parser should parse the Optional production.
 */
function buildLookaheadFuncForOptionalProd(occurrence, ruleGrammar, k, dynamicTokensEnabled, prodType, lookaheadBuilder) {
    const lookAheadPaths = getLookaheadPathsForOptionalProd(occurrence, ruleGrammar, prodType, k);
    const tokenMatcher = areTokenCategoriesNotUsed(lookAheadPaths)
        ? tokenStructuredMatcherNoCategories
        : tokenStructuredMatcher;
    return lookaheadBuilder(lookAheadPaths[0], tokenMatcher, dynamicTokensEnabled);
}
function buildAlternativesLookAheadFunc(alts, hasPredicates, tokenMatcher, dynamicTokensEnabled) {
    const numOfAlts = alts.length;
    const areAllOneTokenLookahead = lodash_es_every(alts, (currAlt) => {
        return lodash_es_every(currAlt, (currPath) => {
            return currPath.length === 1;
        });
    });
    // This version takes into account the predicates as well.
    if (hasPredicates) {
        /**
         * @returns {number} - The chosen alternative index
         */
        return function (orAlts) {
            // unfortunately the predicates must be extracted every single time
            // as they cannot be cached due to references to parameters(vars) which are no longer valid.
            // note that in the common case of no predicates, no cpu time will be wasted on this (see else block)
            const predicates = (0,map/* default */.Z)(orAlts, (currAlt) => currAlt.GATE);
            for (let t = 0; t < numOfAlts; t++) {
                const currAlt = alts[t];
                const currNumOfPaths = currAlt.length;
                const currPredicate = predicates[t];
                if (currPredicate !== undefined && currPredicate.call(this) === false) {
                    // if the predicate does not match there is no point in checking the paths
                    continue;
                }
                nextPath: for (let j = 0; j < currNumOfPaths; j++) {
                    const currPath = currAlt[j];
                    const currPathLength = currPath.length;
                    for (let i = 0; i < currPathLength; i++) {
                        const nextToken = this.LA(i + 1);
                        if (tokenMatcher(nextToken, currPath[i]) === false) {
                            // mismatch in current path
                            // try the next pth
                            continue nextPath;
                        }
                    }
                    // found a full path that matches.
                    // this will also work for an empty ALT as the loop will be skipped
                    return t;
                }
                // none of the paths for the current alternative matched
                // try the next alternative
            }
            // none of the alternatives could be matched
            return undefined;
        };
    }
    else if (areAllOneTokenLookahead && !dynamicTokensEnabled) {
        // optimized (common) case of all the lookaheads paths requiring only
        // a single token lookahead. These Optimizations cannot work if dynamically defined Tokens are used.
        const singleTokenAlts = (0,map/* default */.Z)(alts, (currAlt) => {
            return (0,flatten/* default */.Z)(currAlt);
        });
        const choiceToAlt = (0,reduce/* default */.Z)(singleTokenAlts, (result, currAlt, idx) => {
            (0,forEach/* default */.Z)(currAlt, (currTokType) => {
                if (!(0,has/* default */.Z)(result, currTokType.tokenTypeIdx)) {
                    result[currTokType.tokenTypeIdx] = idx;
                }
                (0,forEach/* default */.Z)(currTokType.categoryMatches, (currExtendingType) => {
                    if (!(0,has/* default */.Z)(result, currExtendingType)) {
                        result[currExtendingType] = idx;
                    }
                });
            });
            return result;
        }, {});
        /**
         * @returns {number} - The chosen alternative index
         */
        return function () {
            const nextToken = this.LA(1);
            return choiceToAlt[nextToken.tokenTypeIdx];
        };
    }
    else {
        // optimized lookahead without needing to check the predicates at all.
        // this causes code duplication which is intentional to improve performance.
        /**
         * @returns {number} - The chosen alternative index
         */
        return function () {
            for (let t = 0; t < numOfAlts; t++) {
                const currAlt = alts[t];
                const currNumOfPaths = currAlt.length;
                nextPath: for (let j = 0; j < currNumOfPaths; j++) {
                    const currPath = currAlt[j];
                    const currPathLength = currPath.length;
                    for (let i = 0; i < currPathLength; i++) {
                        const nextToken = this.LA(i + 1);
                        if (tokenMatcher(nextToken, currPath[i]) === false) {
                            // mismatch in current path
                            // try the next pth
                            continue nextPath;
                        }
                    }
                    // found a full path that matches.
                    // this will also work for an empty ALT as the loop will be skipped
                    return t;
                }
                // none of the paths for the current alternative matched
                // try the next alternative
            }
            // none of the alternatives could be matched
            return undefined;
        };
    }
}
function buildSingleAlternativeLookaheadFunction(alt, tokenMatcher, dynamicTokensEnabled) {
    const areAllOneTokenLookahead = lodash_es_every(alt, (currPath) => {
        return currPath.length === 1;
    });
    const numOfPaths = alt.length;
    // optimized (common) case of all the lookaheads paths requiring only
    // a single token lookahead.
    if (areAllOneTokenLookahead && !dynamicTokensEnabled) {
        const singleTokensTypes = (0,flatten/* default */.Z)(alt);
        if (singleTokensTypes.length === 1 &&
            (0,isEmpty/* default */.Z)(singleTokensTypes[0].categoryMatches)) {
            const expectedTokenType = singleTokensTypes[0];
            const expectedTokenUniqueKey = expectedTokenType.tokenTypeIdx;
            return function () {
                return this.LA(1).tokenTypeIdx === expectedTokenUniqueKey;
            };
        }
        else {
            const choiceToAlt = (0,reduce/* default */.Z)(singleTokensTypes, (result, currTokType, idx) => {
                result[currTokType.tokenTypeIdx] = true;
                (0,forEach/* default */.Z)(currTokType.categoryMatches, (currExtendingType) => {
                    result[currExtendingType] = true;
                });
                return result;
            }, []);
            return function () {
                const nextToken = this.LA(1);
                return choiceToAlt[nextToken.tokenTypeIdx] === true;
            };
        }
    }
    else {
        return function () {
            nextPath: for (let j = 0; j < numOfPaths; j++) {
                const currPath = alt[j];
                const currPathLength = currPath.length;
                for (let i = 0; i < currPathLength; i++) {
                    const nextToken = this.LA(i + 1);
                    if (tokenMatcher(nextToken, currPath[i]) === false) {
                        // mismatch in current path
                        // try the next pth
                        continue nextPath;
                    }
                }
                // found a full path that matches.
                return true;
            }
            // none of the paths matched
            return false;
        };
    }
}
class RestDefinitionFinderWalker extends RestWalker {
    constructor(topProd, targetOccurrence, targetProdType) {
        super();
        this.topProd = topProd;
        this.targetOccurrence = targetOccurrence;
        this.targetProdType = targetProdType;
    }
    startWalking() {
        this.walk(this.topProd);
        return this.restDef;
    }
    checkIsTarget(node, expectedProdType, currRest, prevRest) {
        if (node.idx === this.targetOccurrence &&
            this.targetProdType === expectedProdType) {
            this.restDef = currRest.concat(prevRest);
            return true;
        }
        // performance optimization, do not iterate over the entire Grammar ast after we have found the target
        return false;
    }
    walkOption(optionProd, currRest, prevRest) {
        if (!this.checkIsTarget(optionProd, PROD_TYPE.OPTION, currRest, prevRest)) {
            super.walkOption(optionProd, currRest, prevRest);
        }
    }
    walkAtLeastOne(atLeastOneProd, currRest, prevRest) {
        if (!this.checkIsTarget(atLeastOneProd, PROD_TYPE.REPETITION_MANDATORY, currRest, prevRest)) {
            super.walkOption(atLeastOneProd, currRest, prevRest);
        }
    }
    walkAtLeastOneSep(atLeastOneSepProd, currRest, prevRest) {
        if (!this.checkIsTarget(atLeastOneSepProd, PROD_TYPE.REPETITION_MANDATORY_WITH_SEPARATOR, currRest, prevRest)) {
            super.walkOption(atLeastOneSepProd, currRest, prevRest);
        }
    }
    walkMany(manyProd, currRest, prevRest) {
        if (!this.checkIsTarget(manyProd, PROD_TYPE.REPETITION, currRest, prevRest)) {
            super.walkOption(manyProd, currRest, prevRest);
        }
    }
    walkManySep(manySepProd, currRest, prevRest) {
        if (!this.checkIsTarget(manySepProd, PROD_TYPE.REPETITION_WITH_SEPARATOR, currRest, prevRest)) {
            super.walkOption(manySepProd, currRest, prevRest);
        }
    }
}
/**
 * Returns the definition of a target production in a top level level rule.
 */
class InsideDefinitionFinderVisitor extends GAstVisitor {
    constructor(targetOccurrence, targetProdType, targetRef) {
        super();
        this.targetOccurrence = targetOccurrence;
        this.targetProdType = targetProdType;
        this.targetRef = targetRef;
        this.result = [];
    }
    checkIsTarget(node, expectedProdName) {
        if (node.idx === this.targetOccurrence &&
            this.targetProdType === expectedProdName &&
            (this.targetRef === undefined || node === this.targetRef)) {
            this.result = node.definition;
        }
    }
    visitOption(node) {
        this.checkIsTarget(node, PROD_TYPE.OPTION);
    }
    visitRepetition(node) {
        this.checkIsTarget(node, PROD_TYPE.REPETITION);
    }
    visitRepetitionMandatory(node) {
        this.checkIsTarget(node, PROD_TYPE.REPETITION_MANDATORY);
    }
    visitRepetitionMandatoryWithSeparator(node) {
        this.checkIsTarget(node, PROD_TYPE.REPETITION_MANDATORY_WITH_SEPARATOR);
    }
    visitRepetitionWithSeparator(node) {
        this.checkIsTarget(node, PROD_TYPE.REPETITION_WITH_SEPARATOR);
    }
    visitAlternation(node) {
        this.checkIsTarget(node, PROD_TYPE.ALTERNATION);
    }
}
function initializeArrayOfArrays(size) {
    const result = new Array(size);
    for (let i = 0; i < size; i++) {
        result[i] = [];
    }
    return result;
}
/**
 * A sort of hash function between a Path in the grammar and a string.
 * Note that this returns multiple "hashes" to support the scenario of token categories.
 * -  A single path with categories may match multiple **actual** paths.
 */
function pathToHashKeys(path) {
    let keys = [""];
    for (let i = 0; i < path.length; i++) {
        const tokType = path[i];
        const longerKeys = [];
        for (let j = 0; j < keys.length; j++) {
            const currShorterKey = keys[j];
            longerKeys.push(currShorterKey + "_" + tokType.tokenTypeIdx);
            for (let t = 0; t < tokType.categoryMatches.length; t++) {
                const categoriesKeySuffix = "_" + tokType.categoryMatches[t];
                longerKeys.push(currShorterKey + categoriesKeySuffix);
            }
        }
        keys = longerKeys;
    }
    return keys;
}
/**
 * Imperative style due to being called from a hot spot
 */
function isUniquePrefixHash(altKnownPathsKeys, searchPathKeys, idx) {
    for (let currAltIdx = 0; currAltIdx < altKnownPathsKeys.length; currAltIdx++) {
        // We only want to test vs the other alternatives
        if (currAltIdx === idx) {
            continue;
        }
        const otherAltKnownPathsKeys = altKnownPathsKeys[currAltIdx];
        for (let searchIdx = 0; searchIdx < searchPathKeys.length; searchIdx++) {
            const searchKey = searchPathKeys[searchIdx];
            if (otherAltKnownPathsKeys[searchKey] === true) {
                return false;
            }
        }
    }
    // None of the SearchPathKeys were found in any of the other alternatives
    return true;
}
function lookAheadSequenceFromAlternatives(altsDefs, k) {
    const partialAlts = (0,map/* default */.Z)(altsDefs, (currAlt) => possiblePathsFrom([currAlt], 1));
    const finalResult = initializeArrayOfArrays(partialAlts.length);
    const altsHashes = (0,map/* default */.Z)(partialAlts, (currAltPaths) => {
        const dict = {};
        (0,forEach/* default */.Z)(currAltPaths, (item) => {
            const keys = pathToHashKeys(item.partialPath);
            (0,forEach/* default */.Z)(keys, (currKey) => {
                dict[currKey] = true;
            });
        });
        return dict;
    });
    let newData = partialAlts;
    // maxLookahead loop
    for (let pathLength = 1; pathLength <= k; pathLength++) {
        const currDataset = newData;
        newData = initializeArrayOfArrays(currDataset.length);
        // alternatives loop
        for (let altIdx = 0; altIdx < currDataset.length; altIdx++) {
            const currAltPathsAndSuffixes = currDataset[altIdx];
            // paths in current alternative loop
            for (let currPathIdx = 0; currPathIdx < currAltPathsAndSuffixes.length; currPathIdx++) {
                const currPathPrefix = currAltPathsAndSuffixes[currPathIdx].partialPath;
                const suffixDef = currAltPathsAndSuffixes[currPathIdx].suffixDef;
                const prefixKeys = pathToHashKeys(currPathPrefix);
                const isUnique = isUniquePrefixHash(altsHashes, prefixKeys, altIdx);
                // End of the line for this path.
                if (isUnique || (0,isEmpty/* default */.Z)(suffixDef) || currPathPrefix.length === k) {
                    const currAltResult = finalResult[altIdx];
                    // TODO: Can we implement a containsPath using Maps/Dictionaries?
                    if (containsPath(currAltResult, currPathPrefix) === false) {
                        currAltResult.push(currPathPrefix);
                        // Update all new  keys for the current path.
                        for (let j = 0; j < prefixKeys.length; j++) {
                            const currKey = prefixKeys[j];
                            altsHashes[altIdx][currKey] = true;
                        }
                    }
                }
                // Expand longer paths
                else {
                    const newPartialPathsAndSuffixes = possiblePathsFrom(suffixDef, pathLength + 1, currPathPrefix);
                    newData[altIdx] = newData[altIdx].concat(newPartialPathsAndSuffixes);
                    // Update keys for new known paths
                    (0,forEach/* default */.Z)(newPartialPathsAndSuffixes, (item) => {
                        const prefixKeys = pathToHashKeys(item.partialPath);
                        (0,forEach/* default */.Z)(prefixKeys, (key) => {
                            altsHashes[altIdx][key] = true;
                        });
                    });
                }
            }
        }
    }
    return finalResult;
}
function getLookaheadPathsForOr(occurrence, ruleGrammar, k, orProd) {
    const visitor = new InsideDefinitionFinderVisitor(occurrence, PROD_TYPE.ALTERNATION, orProd);
    ruleGrammar.accept(visitor);
    return lookAheadSequenceFromAlternatives(visitor.result, k);
}
function getLookaheadPathsForOptionalProd(occurrence, ruleGrammar, prodType, k) {
    const insideDefVisitor = new InsideDefinitionFinderVisitor(occurrence, prodType);
    ruleGrammar.accept(insideDefVisitor);
    const insideDef = insideDefVisitor.result;
    const afterDefWalker = new RestDefinitionFinderWalker(ruleGrammar, occurrence, prodType);
    const afterDef = afterDefWalker.startWalking();
    const insideFlat = new Alternative({ definition: insideDef });
    const afterFlat = new Alternative({ definition: afterDef });
    return lookAheadSequenceFromAlternatives([insideFlat, afterFlat], k);
}
function containsPath(alternative, searchPath) {
    compareOtherPath: for (let i = 0; i < alternative.length; i++) {
        const otherPath = alternative[i];
        if (otherPath.length !== searchPath.length) {
            continue;
        }
        for (let j = 0; j < otherPath.length; j++) {
            const searchTok = searchPath[j];
            const otherTok = otherPath[j];
            const matchingTokens = searchTok === otherTok ||
                otherTok.categoryMatchesMap[searchTok.tokenTypeIdx] !== undefined;
            if (matchingTokens === false) {
                continue compareOtherPath;
            }
        }
        return true;
    }
    return false;
}
function isStrictPrefixOfPath(prefix, other) {
    return (prefix.length < other.length &&
        lodash_es_every(prefix, (tokType, idx) => {
            const otherTokType = other[idx];
            return (tokType === otherTokType ||
                otherTokType.categoryMatchesMap[tokType.tokenTypeIdx]);
        }));
}
function areTokenCategoriesNotUsed(lookAheadPaths) {
    return lodash_es_every(lookAheadPaths, (singleAltPaths) => lodash_es_every(singleAltPaths, (singlePath) => lodash_es_every(singlePath, (token) => (0,isEmpty/* default */.Z)(token.categoryMatches))));
}
//# sourceMappingURL=lookahead.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/checks.js






function validateLookahead(options) {
    const lookaheadValidationErrorMessages = options.lookaheadStrategy.validate({
        rules: options.rules,
        tokenTypes: options.tokenTypes,
        grammarName: options.grammarName,
    });
    return (0,map/* default */.Z)(lookaheadValidationErrorMessages, (errorMessage) => (Object.assign({ type: ParserDefinitionErrorType.CUSTOM_LOOKAHEAD_VALIDATION }, errorMessage)));
}
function validateGrammar(topLevels, tokenTypes, errMsgProvider, grammarName) {
    const duplicateErrors = (0,flatMap/* default */.Z)(topLevels, (currTopLevel) => validateDuplicateProductions(currTopLevel, errMsgProvider));
    const termsNamespaceConflictErrors = checkTerminalAndNoneTerminalsNameSpace(topLevels, tokenTypes, errMsgProvider);
    const tooManyAltsErrors = (0,flatMap/* default */.Z)(topLevels, (curRule) => validateTooManyAlts(curRule, errMsgProvider));
    const duplicateRulesError = (0,flatMap/* default */.Z)(topLevels, (curRule) => validateRuleDoesNotAlreadyExist(curRule, topLevels, grammarName, errMsgProvider));
    return duplicateErrors.concat(termsNamespaceConflictErrors, tooManyAltsErrors, duplicateRulesError);
}
function validateDuplicateProductions(topLevelRule, errMsgProvider) {
    const collectorVisitor = new OccurrenceValidationCollector();
    topLevelRule.accept(collectorVisitor);
    const allRuleProductions = collectorVisitor.allProductions;
    const productionGroups = lodash_es_groupBy(allRuleProductions, identifyProductionForDuplicates);
    const duplicates = lodash_es_pickBy(productionGroups, (currGroup) => {
        return currGroup.length > 1;
    });
    const errors = (0,map/* default */.Z)((0,values/* default */.Z)(duplicates), (currDuplicates) => {
        const firstProd = lodash_es_head(currDuplicates);
        const msg = errMsgProvider.buildDuplicateFoundError(topLevelRule, currDuplicates);
        const dslName = getProductionDslName(firstProd);
        const defError = {
            message: msg,
            type: ParserDefinitionErrorType.DUPLICATE_PRODUCTIONS,
            ruleName: topLevelRule.name,
            dslName: dslName,
            occurrence: firstProd.idx,
        };
        const param = getExtraProductionArgument(firstProd);
        if (param) {
            defError.parameter = param;
        }
        return defError;
    });
    return errors;
}
function identifyProductionForDuplicates(prod) {
    return `${getProductionDslName(prod)}_#_${prod.idx}_#_${getExtraProductionArgument(prod)}`;
}
function getExtraProductionArgument(prod) {
    if (prod instanceof Terminal) {
        return prod.terminalType.name;
    }
    else if (prod instanceof NonTerminal) {
        return prod.nonTerminalName;
    }
    else {
        return "";
    }
}
class OccurrenceValidationCollector extends GAstVisitor {
    constructor() {
        super(...arguments);
        this.allProductions = [];
    }
    visitNonTerminal(subrule) {
        this.allProductions.push(subrule);
    }
    visitOption(option) {
        this.allProductions.push(option);
    }
    visitRepetitionWithSeparator(manySep) {
        this.allProductions.push(manySep);
    }
    visitRepetitionMandatory(atLeastOne) {
        this.allProductions.push(atLeastOne);
    }
    visitRepetitionMandatoryWithSeparator(atLeastOneSep) {
        this.allProductions.push(atLeastOneSep);
    }
    visitRepetition(many) {
        this.allProductions.push(many);
    }
    visitAlternation(or) {
        this.allProductions.push(or);
    }
    visitTerminal(terminal) {
        this.allProductions.push(terminal);
    }
}
function validateRuleDoesNotAlreadyExist(rule, allRules, className, errMsgProvider) {
    const errors = [];
    const occurrences = (0,reduce/* default */.Z)(allRules, (result, curRule) => {
        if (curRule.name === rule.name) {
            return result + 1;
        }
        return result;
    }, 0);
    if (occurrences > 1) {
        const errMsg = errMsgProvider.buildDuplicateRuleNameError({
            topLevelRule: rule,
            grammarName: className,
        });
        errors.push({
            message: errMsg,
            type: ParserDefinitionErrorType.DUPLICATE_RULE_NAME,
            ruleName: rule.name,
        });
    }
    return errors;
}
// TODO: is there anyway to get only the rule names of rules inherited from the super grammars?
// This is not part of the IGrammarErrorProvider because the validation cannot be performed on
// The grammar structure, only at runtime.
function validateRuleIsOverridden(ruleName, definedRulesNames, className) {
    const errors = [];
    let errMsg;
    if (!lodash_es_includes(definedRulesNames, ruleName)) {
        errMsg =
            `Invalid rule override, rule: ->${ruleName}<- cannot be overridden in the grammar: ->${className}<-` +
                `as it is not defined in any of the super grammars `;
        errors.push({
            message: errMsg,
            type: ParserDefinitionErrorType.INVALID_RULE_OVERRIDE,
            ruleName: ruleName,
        });
    }
    return errors;
}
function validateNoLeftRecursion(topRule, currRule, errMsgProvider, path = []) {
    const errors = [];
    const nextNonTerminals = getFirstNoneTerminal(currRule.definition);
    if ((0,isEmpty/* default */.Z)(nextNonTerminals)) {
        return [];
    }
    else {
        const ruleName = topRule.name;
        const foundLeftRecursion = lodash_es_includes(nextNonTerminals, topRule);
        if (foundLeftRecursion) {
            errors.push({
                message: errMsgProvider.buildLeftRecursionError({
                    topLevelRule: topRule,
                    leftRecursionPath: path,
                }),
                type: ParserDefinitionErrorType.LEFT_RECURSION,
                ruleName: ruleName,
            });
        }
        // we are only looking for cyclic paths leading back to the specific topRule
        // other cyclic paths are ignored, we still need this difference to avoid infinite loops...
        const validNextSteps = lodash_es_difference(nextNonTerminals, path.concat([topRule]));
        const errorsFromNextSteps = (0,flatMap/* default */.Z)(validNextSteps, (currRefRule) => {
            const newPath = (0,lodash_es_clone/* default */.Z)(path);
            newPath.push(currRefRule);
            return validateNoLeftRecursion(topRule, currRefRule, errMsgProvider, newPath);
        });
        return errors.concat(errorsFromNextSteps);
    }
}
function getFirstNoneTerminal(definition) {
    let result = [];
    if ((0,isEmpty/* default */.Z)(definition)) {
        return result;
    }
    const firstProd = lodash_es_head(definition);
    /* istanbul ignore else */
    if (firstProd instanceof NonTerminal) {
        result.push(firstProd.referencedRule);
    }
    else if (firstProd instanceof Alternative ||
        firstProd instanceof Option ||
        firstProd instanceof RepetitionMandatory ||
        firstProd instanceof RepetitionMandatoryWithSeparator ||
        firstProd instanceof RepetitionWithSeparator ||
        firstProd instanceof Repetition) {
        result = result.concat(getFirstNoneTerminal(firstProd.definition));
    }
    else if (firstProd instanceof Alternation) {
        // each sub definition in alternation is a FLAT
        result = (0,flatten/* default */.Z)((0,map/* default */.Z)(firstProd.definition, (currSubDef) => getFirstNoneTerminal(currSubDef.definition)));
    }
    else if (firstProd instanceof Terminal) {
        // nothing to see, move along
    }
    else {
        throw Error("non exhaustive match");
    }
    const isFirstOptional = isOptionalProd(firstProd);
    const hasMore = definition.length > 1;
    if (isFirstOptional && hasMore) {
        const rest = lodash_es_drop(definition);
        return result.concat(getFirstNoneTerminal(rest));
    }
    else {
        return result;
    }
}
class OrCollector extends GAstVisitor {
    constructor() {
        super(...arguments);
        this.alternations = [];
    }
    visitAlternation(node) {
        this.alternations.push(node);
    }
}
function validateEmptyOrAlternative(topLevelRule, errMsgProvider) {
    const orCollector = new OrCollector();
    topLevelRule.accept(orCollector);
    const ors = orCollector.alternations;
    const errors = (0,flatMap/* default */.Z)(ors, (currOr) => {
        const exceptLast = lodash_es_dropRight(currOr.definition);
        return (0,flatMap/* default */.Z)(exceptLast, (currAlternative, currAltIdx) => {
            const possibleFirstInAlt = nextPossibleTokensAfter([currAlternative], [], tokenStructuredMatcher, 1);
            if ((0,isEmpty/* default */.Z)(possibleFirstInAlt)) {
                return [
                    {
                        message: errMsgProvider.buildEmptyAlternationError({
                            topLevelRule: topLevelRule,
                            alternation: currOr,
                            emptyChoiceIdx: currAltIdx,
                        }),
                        type: ParserDefinitionErrorType.NONE_LAST_EMPTY_ALT,
                        ruleName: topLevelRule.name,
                        occurrence: currOr.idx,
                        alternative: currAltIdx + 1,
                    },
                ];
            }
            else {
                return [];
            }
        });
    });
    return errors;
}
function validateAmbiguousAlternationAlternatives(topLevelRule, globalMaxLookahead, errMsgProvider) {
    const orCollector = new OrCollector();
    topLevelRule.accept(orCollector);
    let ors = orCollector.alternations;
    // New Handling of ignoring ambiguities
    // - https://github.com/chevrotain/chevrotain/issues/869
    ors = lodash_es_reject(ors, (currOr) => currOr.ignoreAmbiguities === true);
    const errors = (0,flatMap/* default */.Z)(ors, (currOr) => {
        const currOccurrence = currOr.idx;
        const actualMaxLookahead = currOr.maxLookahead || globalMaxLookahead;
        const alternatives = getLookaheadPathsForOr(currOccurrence, topLevelRule, actualMaxLookahead, currOr);
        const altsAmbiguityErrors = checkAlternativesAmbiguities(alternatives, currOr, topLevelRule, errMsgProvider);
        const altsPrefixAmbiguityErrors = checkPrefixAlternativesAmbiguities(alternatives, currOr, topLevelRule, errMsgProvider);
        return altsAmbiguityErrors.concat(altsPrefixAmbiguityErrors);
    });
    return errors;
}
class RepetitionCollector extends GAstVisitor {
    constructor() {
        super(...arguments);
        this.allProductions = [];
    }
    visitRepetitionWithSeparator(manySep) {
        this.allProductions.push(manySep);
    }
    visitRepetitionMandatory(atLeastOne) {
        this.allProductions.push(atLeastOne);
    }
    visitRepetitionMandatoryWithSeparator(atLeastOneSep) {
        this.allProductions.push(atLeastOneSep);
    }
    visitRepetition(many) {
        this.allProductions.push(many);
    }
}
function validateTooManyAlts(topLevelRule, errMsgProvider) {
    const orCollector = new OrCollector();
    topLevelRule.accept(orCollector);
    const ors = orCollector.alternations;
    const errors = (0,flatMap/* default */.Z)(ors, (currOr) => {
        if (currOr.definition.length > 255) {
            return [
                {
                    message: errMsgProvider.buildTooManyAlternativesError({
                        topLevelRule: topLevelRule,
                        alternation: currOr,
                    }),
                    type: ParserDefinitionErrorType.TOO_MANY_ALTS,
                    ruleName: topLevelRule.name,
                    occurrence: currOr.idx,
                },
            ];
        }
        else {
            return [];
        }
    });
    return errors;
}
function validateSomeNonEmptyLookaheadPath(topLevelRules, maxLookahead, errMsgProvider) {
    const errors = [];
    (0,forEach/* default */.Z)(topLevelRules, (currTopRule) => {
        const collectorVisitor = new RepetitionCollector();
        currTopRule.accept(collectorVisitor);
        const allRuleProductions = collectorVisitor.allProductions;
        (0,forEach/* default */.Z)(allRuleProductions, (currProd) => {
            const prodType = getProdType(currProd);
            const actualMaxLookahead = currProd.maxLookahead || maxLookahead;
            const currOccurrence = currProd.idx;
            const paths = getLookaheadPathsForOptionalProd(currOccurrence, currTopRule, prodType, actualMaxLookahead);
            const pathsInsideProduction = paths[0];
            if ((0,isEmpty/* default */.Z)((0,flatten/* default */.Z)(pathsInsideProduction))) {
                const errMsg = errMsgProvider.buildEmptyRepetitionError({
                    topLevelRule: currTopRule,
                    repetition: currProd,
                });
                errors.push({
                    message: errMsg,
                    type: ParserDefinitionErrorType.NO_NON_EMPTY_LOOKAHEAD,
                    ruleName: currTopRule.name,
                });
            }
        });
    });
    return errors;
}
function checkAlternativesAmbiguities(alternatives, alternation, rule, errMsgProvider) {
    const foundAmbiguousPaths = [];
    const identicalAmbiguities = (0,reduce/* default */.Z)(alternatives, (result, currAlt, currAltIdx) => {
        // ignore (skip) ambiguities with this alternative
        if (alternation.definition[currAltIdx].ignoreAmbiguities === true) {
            return result;
        }
        (0,forEach/* default */.Z)(currAlt, (currPath) => {
            const altsCurrPathAppearsIn = [currAltIdx];
            (0,forEach/* default */.Z)(alternatives, (currOtherAlt, currOtherAltIdx) => {
                if (currAltIdx !== currOtherAltIdx &&
                    containsPath(currOtherAlt, currPath) &&
                    // ignore (skip) ambiguities with this "other" alternative
                    alternation.definition[currOtherAltIdx].ignoreAmbiguities !== true) {
                    altsCurrPathAppearsIn.push(currOtherAltIdx);
                }
            });
            if (altsCurrPathAppearsIn.length > 1 &&
                !containsPath(foundAmbiguousPaths, currPath)) {
                foundAmbiguousPaths.push(currPath);
                result.push({
                    alts: altsCurrPathAppearsIn,
                    path: currPath,
                });
            }
        });
        return result;
    }, []);
    const currErrors = (0,map/* default */.Z)(identicalAmbiguities, (currAmbDescriptor) => {
        const ambgIndices = (0,map/* default */.Z)(currAmbDescriptor.alts, (currAltIdx) => currAltIdx + 1);
        const currMessage = errMsgProvider.buildAlternationAmbiguityError({
            topLevelRule: rule,
            alternation: alternation,
            ambiguityIndices: ambgIndices,
            prefixPath: currAmbDescriptor.path,
        });
        return {
            message: currMessage,
            type: ParserDefinitionErrorType.AMBIGUOUS_ALTS,
            ruleName: rule.name,
            occurrence: alternation.idx,
            alternatives: currAmbDescriptor.alts,
        };
    });
    return currErrors;
}
function checkPrefixAlternativesAmbiguities(alternatives, alternation, rule, errMsgProvider) {
    // flatten
    const pathsAndIndices = (0,reduce/* default */.Z)(alternatives, (result, currAlt, idx) => {
        const currPathsAndIdx = (0,map/* default */.Z)(currAlt, (currPath) => {
            return { idx: idx, path: currPath };
        });
        return result.concat(currPathsAndIdx);
    }, []);
    const errors = lodash_es_compact((0,flatMap/* default */.Z)(pathsAndIndices, (currPathAndIdx) => {
        const alternativeGast = alternation.definition[currPathAndIdx.idx];
        // ignore (skip) ambiguities with this alternative
        if (alternativeGast.ignoreAmbiguities === true) {
            return [];
        }
        const targetIdx = currPathAndIdx.idx;
        const targetPath = currPathAndIdx.path;
        const prefixAmbiguitiesPathsAndIndices = (0,filter/* default */.Z)(pathsAndIndices, (searchPathAndIdx) => {
            // prefix ambiguity can only be created from lower idx (higher priority) path
            return (
            // ignore (skip) ambiguities with this "other" alternative
            alternation.definition[searchPathAndIdx.idx].ignoreAmbiguities !==
                true &&
                searchPathAndIdx.idx < targetIdx &&
                // checking for strict prefix because identical lookaheads
                // will be be detected using a different validation.
                isStrictPrefixOfPath(searchPathAndIdx.path, targetPath));
        });
        const currPathPrefixErrors = (0,map/* default */.Z)(prefixAmbiguitiesPathsAndIndices, (currAmbPathAndIdx) => {
            const ambgIndices = [currAmbPathAndIdx.idx + 1, targetIdx + 1];
            const occurrence = alternation.idx === 0 ? "" : alternation.idx;
            const message = errMsgProvider.buildAlternationPrefixAmbiguityError({
                topLevelRule: rule,
                alternation: alternation,
                ambiguityIndices: ambgIndices,
                prefixPath: currAmbPathAndIdx.path,
            });
            return {
                message: message,
                type: ParserDefinitionErrorType.AMBIGUOUS_PREFIX_ALTS,
                ruleName: rule.name,
                occurrence: occurrence,
                alternatives: ambgIndices,
            };
        });
        return currPathPrefixErrors;
    }));
    return errors;
}
function checkTerminalAndNoneTerminalsNameSpace(topLevels, tokenTypes, errMsgProvider) {
    const errors = [];
    const tokenNames = (0,map/* default */.Z)(tokenTypes, (currToken) => currToken.name);
    (0,forEach/* default */.Z)(topLevels, (currRule) => {
        const currRuleName = currRule.name;
        if (lodash_es_includes(tokenNames, currRuleName)) {
            const errMsg = errMsgProvider.buildNamespaceConflictError(currRule);
            errors.push({
                message: errMsg,
                type: ParserDefinitionErrorType.CONFLICT_TOKENS_RULES_NAMESPACE,
                ruleName: currRuleName,
            });
        }
    });
    return errors;
}
//# sourceMappingURL=checks.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/gast/gast_resolver_public.js




function gast_resolver_public_resolveGrammar(options) {
    const actualOptions = (0,defaults/* default */.Z)(options, {
        errMsgProvider: defaultGrammarResolverErrorProvider,
    });
    const topRulesTable = {};
    (0,forEach/* default */.Z)(options.rules, (rule) => {
        topRulesTable[rule.name] = rule;
    });
    return resolveGrammar(topRulesTable, actualOptions.errMsgProvider);
}
function gast_resolver_public_validateGrammar(options) {
    options = (0,defaults/* default */.Z)(options, {
        errMsgProvider: defaultGrammarValidatorErrorProvider,
    });
    return validateGrammar(options.rules, options.tokenTypes, options.errMsgProvider, options.grammarName);
}
//# sourceMappingURL=gast_resolver_public.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/exceptions_public.js

const MISMATCHED_TOKEN_EXCEPTION = "MismatchedTokenException";
const NO_VIABLE_ALT_EXCEPTION = "NoViableAltException";
const EARLY_EXIT_EXCEPTION = "EarlyExitException";
const NOT_ALL_INPUT_PARSED_EXCEPTION = "NotAllInputParsedException";
const RECOGNITION_EXCEPTION_NAMES = [
    MISMATCHED_TOKEN_EXCEPTION,
    NO_VIABLE_ALT_EXCEPTION,
    EARLY_EXIT_EXCEPTION,
    NOT_ALL_INPUT_PARSED_EXCEPTION,
];
Object.freeze(RECOGNITION_EXCEPTION_NAMES);
// hacks to bypass no support for custom Errors in javascript/typescript
function isRecognitionException(error) {
    // can't do instanceof on hacked custom js exceptions
    return lodash_es_includes(RECOGNITION_EXCEPTION_NAMES, error.name);
}
class RecognitionException extends Error {
    constructor(message, token) {
        super(message);
        this.token = token;
        this.resyncedTokens = [];
        // fix prototype chain when typescript target is ES5
        Object.setPrototypeOf(this, new.target.prototype);
        /* istanbul ignore next - V8 workaround to remove constructor from stacktrace when typescript target is ES5 */
        if (Error.captureStackTrace) {
            Error.captureStackTrace(this, this.constructor);
        }
    }
}
class MismatchedTokenException extends RecognitionException {
    constructor(message, token, previousToken) {
        super(message, token);
        this.previousToken = previousToken;
        this.name = MISMATCHED_TOKEN_EXCEPTION;
    }
}
class NoViableAltException extends RecognitionException {
    constructor(message, token, previousToken) {
        super(message, token);
        this.previousToken = previousToken;
        this.name = NO_VIABLE_ALT_EXCEPTION;
    }
}
class NotAllInputParsedException extends RecognitionException {
    constructor(message, token) {
        super(message, token);
        this.name = NOT_ALL_INPUT_PARSED_EXCEPTION;
    }
}
class EarlyExitException extends RecognitionException {
    constructor(message, token, previousToken) {
        super(message, token);
        this.previousToken = previousToken;
        this.name = EARLY_EXIT_EXCEPTION;
    }
}
//# sourceMappingURL=exceptions_public.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/recoverable.js





const EOF_FOLLOW_KEY = {};
const IN_RULE_RECOVERY_EXCEPTION = "InRuleRecoveryException";
class InRuleRecoveryException extends Error {
    constructor(message) {
        super(message);
        this.name = IN_RULE_RECOVERY_EXCEPTION;
    }
}
/**
 * This trait is responsible for the error recovery and fault tolerant logic
 */
class Recoverable {
    initRecoverable(config) {
        this.firstAfterRepMap = {};
        this.resyncFollows = {};
        this.recoveryEnabled = (0,has/* default */.Z)(config, "recoveryEnabled")
            ? config.recoveryEnabled // assumes end user provides the correct config value/type
            : DEFAULT_PARSER_CONFIG.recoveryEnabled;
        // performance optimization, NOOP will be inlined which
        // effectively means that this optional feature does not exist
        // when not used.
        if (this.recoveryEnabled) {
            this.attemptInRepetitionRecovery = attemptInRepetitionRecovery;
        }
    }
    getTokenToInsert(tokType) {
        const tokToInsert = createTokenInstance(tokType, "", NaN, NaN, NaN, NaN, NaN, NaN);
        tokToInsert.isInsertedInRecovery = true;
        return tokToInsert;
    }
    canTokenTypeBeInsertedInRecovery(tokType) {
        return true;
    }
    canTokenTypeBeDeletedInRecovery(tokType) {
        return true;
    }
    tryInRepetitionRecovery(grammarRule, grammarRuleArgs, lookAheadFunc, expectedTokType) {
        // TODO: can the resyncTokenType be cached?
        const reSyncTokType = this.findReSyncTokenType();
        const savedLexerState = this.exportLexerState();
        const resyncedTokens = [];
        let passedResyncPoint = false;
        const nextTokenWithoutResync = this.LA(1);
        let currToken = this.LA(1);
        const generateErrorMessage = () => {
            const previousToken = this.LA(0);
            // we are preemptively re-syncing before an error has been detected, therefor we must reproduce
            // the error that would have been thrown
            const msg = this.errorMessageProvider.buildMismatchTokenMessage({
                expected: expectedTokType,
                actual: nextTokenWithoutResync,
                previous: previousToken,
                ruleName: this.getCurrRuleFullName(),
            });
            const error = new MismatchedTokenException(msg, nextTokenWithoutResync, this.LA(0));
            // the first token here will be the original cause of the error, this is not part of the resyncedTokens property.
            error.resyncedTokens = lodash_es_dropRight(resyncedTokens);
            this.SAVE_ERROR(error);
        };
        while (!passedResyncPoint) {
            // re-synced to a point where we can safely exit the repetition/
            if (this.tokenMatcher(currToken, expectedTokType)) {
                generateErrorMessage();
                return; // must return here to avoid reverting the inputIdx
            }
            else if (lookAheadFunc.call(this)) {
                // we skipped enough tokens so we can resync right back into another iteration of the repetition grammar rule
                generateErrorMessage();
                // recursive invocation in other to support multiple re-syncs in the same top level repetition grammar rule
                grammarRule.apply(this, grammarRuleArgs);
                return; // must return here to avoid reverting the inputIdx
            }
            else if (this.tokenMatcher(currToken, reSyncTokType)) {
                passedResyncPoint = true;
            }
            else {
                currToken = this.SKIP_TOKEN();
                this.addToResyncTokens(currToken, resyncedTokens);
            }
        }
        // we were unable to find a CLOSER point to resync inside the Repetition, reset the state.
        // The parsing exception we were trying to prevent will happen in the NEXT parsing step. it may be handled by
        // "between rules" resync recovery later in the flow.
        this.importLexerState(savedLexerState);
    }
    shouldInRepetitionRecoveryBeTried(expectTokAfterLastMatch, nextTokIdx, notStuck) {
        // Edge case of arriving from a MANY repetition which is stuck
        // Attempting recovery in this case could cause an infinite loop
        if (notStuck === false) {
            return false;
        }
        // no need to recover, next token is what we expect...
        if (this.tokenMatcher(this.LA(1), expectTokAfterLastMatch)) {
            return false;
        }
        // error recovery is disabled during backtracking as it can make the parser ignore a valid grammar path
        // and prefer some backtracking path that includes recovered errors.
        if (this.isBackTracking()) {
            return false;
        }
        // if we can perform inRule recovery (single token insertion or deletion) we always prefer that recovery algorithm
        // because if it works, it makes the least amount of changes to the input stream (greedy algorithm)
        //noinspection RedundantIfStatementJS
        if (this.canPerformInRuleRecovery(expectTokAfterLastMatch, this.getFollowsForInRuleRecovery(expectTokAfterLastMatch, nextTokIdx))) {
            return false;
        }
        return true;
    }
    // Error Recovery functionality
    getFollowsForInRuleRecovery(tokType, tokIdxInRule) {
        const grammarPath = this.getCurrentGrammarPath(tokType, tokIdxInRule);
        const follows = this.getNextPossibleTokenTypes(grammarPath);
        return follows;
    }
    tryInRuleRecovery(expectedTokType, follows) {
        if (this.canRecoverWithSingleTokenInsertion(expectedTokType, follows)) {
            const tokToInsert = this.getTokenToInsert(expectedTokType);
            return tokToInsert;
        }
        if (this.canRecoverWithSingleTokenDeletion(expectedTokType)) {
            const nextTok = this.SKIP_TOKEN();
            this.consumeToken();
            return nextTok;
        }
        throw new InRuleRecoveryException("sad sad panda");
    }
    canPerformInRuleRecovery(expectedToken, follows) {
        return (this.canRecoverWithSingleTokenInsertion(expectedToken, follows) ||
            this.canRecoverWithSingleTokenDeletion(expectedToken));
    }
    canRecoverWithSingleTokenInsertion(expectedTokType, follows) {
        if (!this.canTokenTypeBeInsertedInRecovery(expectedTokType)) {
            return false;
        }
        // must know the possible following tokens to perform single token insertion
        if ((0,isEmpty/* default */.Z)(follows)) {
            return false;
        }
        const mismatchedTok = this.LA(1);
        const isMisMatchedTokInFollows = (0,find/* default */.Z)(follows, (possibleFollowsTokType) => {
            return this.tokenMatcher(mismatchedTok, possibleFollowsTokType);
        }) !== undefined;
        return isMisMatchedTokInFollows;
    }
    canRecoverWithSingleTokenDeletion(expectedTokType) {
        if (!this.canTokenTypeBeDeletedInRecovery(expectedTokType)) {
            return false;
        }
        const isNextTokenWhatIsExpected = this.tokenMatcher(this.LA(2), expectedTokType);
        return isNextTokenWhatIsExpected;
    }
    isInCurrentRuleReSyncSet(tokenTypeIdx) {
        const followKey = this.getCurrFollowKey();
        const currentRuleReSyncSet = this.getFollowSetFromFollowKey(followKey);
        return lodash_es_includes(currentRuleReSyncSet, tokenTypeIdx);
    }
    findReSyncTokenType() {
        const allPossibleReSyncTokTypes = this.flattenFollowSet();
        // this loop will always terminate as EOF is always in the follow stack and also always (virtually) in the input
        let nextToken = this.LA(1);
        let k = 2;
        while (true) {
            const foundMatch = (0,find/* default */.Z)(allPossibleReSyncTokTypes, (resyncTokType) => {
                const canMatch = tokenMatcher(nextToken, resyncTokType);
                return canMatch;
            });
            if (foundMatch !== undefined) {
                return foundMatch;
            }
            nextToken = this.LA(k);
            k++;
        }
    }
    getCurrFollowKey() {
        // the length is at least one as we always add the ruleName to the stack before invoking the rule.
        if (this.RULE_STACK.length === 1) {
            return EOF_FOLLOW_KEY;
        }
        const currRuleShortName = this.getLastExplicitRuleShortName();
        const currRuleIdx = this.getLastExplicitRuleOccurrenceIndex();
        const prevRuleShortName = this.getPreviousExplicitRuleShortName();
        return {
            ruleName: this.shortRuleNameToFullName(currRuleShortName),
            idxInCallingRule: currRuleIdx,
            inRule: this.shortRuleNameToFullName(prevRuleShortName),
        };
    }
    buildFullFollowKeyStack() {
        const explicitRuleStack = this.RULE_STACK;
        const explicitOccurrenceStack = this.RULE_OCCURRENCE_STACK;
        return (0,map/* default */.Z)(explicitRuleStack, (ruleName, idx) => {
            if (idx === 0) {
                return EOF_FOLLOW_KEY;
            }
            return {
                ruleName: this.shortRuleNameToFullName(ruleName),
                idxInCallingRule: explicitOccurrenceStack[idx],
                inRule: this.shortRuleNameToFullName(explicitRuleStack[idx - 1]),
            };
        });
    }
    flattenFollowSet() {
        const followStack = (0,map/* default */.Z)(this.buildFullFollowKeyStack(), (currKey) => {
            return this.getFollowSetFromFollowKey(currKey);
        });
        return (0,flatten/* default */.Z)(followStack);
    }
    getFollowSetFromFollowKey(followKey) {
        if (followKey === EOF_FOLLOW_KEY) {
            return [EOF];
        }
        const followName = followKey.ruleName + followKey.idxInCallingRule + constants_IN + followKey.inRule;
        return this.resyncFollows[followName];
    }
    // It does not make any sense to include a virtual EOF token in the list of resynced tokens
    // as EOF does not really exist and thus does not contain any useful information (line/column numbers)
    addToResyncTokens(token, resyncTokens) {
        if (!this.tokenMatcher(token, EOF)) {
            resyncTokens.push(token);
        }
        return resyncTokens;
    }
    reSyncTo(tokType) {
        const resyncedTokens = [];
        let nextTok = this.LA(1);
        while (this.tokenMatcher(nextTok, tokType) === false) {
            nextTok = this.SKIP_TOKEN();
            this.addToResyncTokens(nextTok, resyncedTokens);
        }
        // the last token is not part of the error.
        return lodash_es_dropRight(resyncedTokens);
    }
    attemptInRepetitionRecovery(prodFunc, args, lookaheadFunc, dslMethodIdx, prodOccurrence, nextToksWalker, notStuck) {
        // by default this is a NO-OP
        // The actual implementation is with the function(not method) below
    }
    getCurrentGrammarPath(tokType, tokIdxInRule) {
        const pathRuleStack = this.getHumanReadableRuleStack();
        const pathOccurrenceStack = (0,lodash_es_clone/* default */.Z)(this.RULE_OCCURRENCE_STACK);
        const grammarPath = {
            ruleStack: pathRuleStack,
            occurrenceStack: pathOccurrenceStack,
            lastTok: tokType,
            lastTokOccurrence: tokIdxInRule,
        };
        return grammarPath;
    }
    getHumanReadableRuleStack() {
        return (0,map/* default */.Z)(this.RULE_STACK, (currShortName) => this.shortRuleNameToFullName(currShortName));
    }
}
function attemptInRepetitionRecovery(prodFunc, args, lookaheadFunc, dslMethodIdx, prodOccurrence, nextToksWalker, notStuck) {
    const key = this.getKeyForAutomaticLookahead(dslMethodIdx, prodOccurrence);
    let firstAfterRepInfo = this.firstAfterRepMap[key];
    if (firstAfterRepInfo === undefined) {
        const currRuleName = this.getCurrRuleFullName();
        const ruleGrammar = this.getGAstProductions()[currRuleName];
        const walker = new nextToksWalker(ruleGrammar, prodOccurrence);
        firstAfterRepInfo = walker.startWalking();
        this.firstAfterRepMap[key] = firstAfterRepInfo;
    }
    let expectTokAfterLastMatch = firstAfterRepInfo.token;
    let nextTokIdx = firstAfterRepInfo.occurrence;
    const isEndOfRule = firstAfterRepInfo.isEndOfRule;
    // special edge case of a TOP most repetition after which the input should END.
    // this will force an attempt for inRule recovery in that scenario.
    if (this.RULE_STACK.length === 1 &&
        isEndOfRule &&
        expectTokAfterLastMatch === undefined) {
        expectTokAfterLastMatch = EOF;
        nextTokIdx = 1;
    }
    // We don't have anything to re-sync to...
    // this condition was extracted from `shouldInRepetitionRecoveryBeTried` to act as a type-guard
    if (expectTokAfterLastMatch === undefined || nextTokIdx === undefined) {
        return;
    }
    if (this.shouldInRepetitionRecoveryBeTried(expectTokAfterLastMatch, nextTokIdx, notStuck)) {
        // TODO: performance optimization: instead of passing the original args here, we modify
        // the args param (or create a new one) and make sure the lookahead func is explicitly provided
        // to avoid searching the cache for it once more.
        this.tryInRepetitionRecovery(prodFunc, args, lookaheadFunc, expectTokAfterLastMatch);
    }
}
//# sourceMappingURL=recoverable.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/keys.js
// Lookahead keys are 32Bit integers in the form
// TTTTTTTT-ZZZZZZZZZZZZ-YYYY-XXXXXXXX
// XXXX -> Occurrence Index bitmap.
// YYYY -> DSL Method Type bitmap.
// ZZZZZZZZZZZZZZZ -> Rule short Index bitmap.
// TTTTTTTTT -> alternation alternative index bitmap
const BITS_FOR_METHOD_TYPE = 4;
const BITS_FOR_OCCURRENCE_IDX = 8;
const BITS_FOR_RULE_IDX = 12;
// TODO: validation, this means that there may at most 2^8 --> 256 alternatives for an alternation.
const BITS_FOR_ALT_IDX = 8;
// short string used as part of mapping keys.
// being short improves the performance when composing KEYS for maps out of these
// The 5 - 8 bits (16 possible values, are reserved for the DSL method indices)
const OR_IDX = 1 << BITS_FOR_OCCURRENCE_IDX;
const OPTION_IDX = 2 << BITS_FOR_OCCURRENCE_IDX;
const MANY_IDX = 3 << BITS_FOR_OCCURRENCE_IDX;
const AT_LEAST_ONE_IDX = 4 << BITS_FOR_OCCURRENCE_IDX;
const MANY_SEP_IDX = 5 << BITS_FOR_OCCURRENCE_IDX;
const AT_LEAST_ONE_SEP_IDX = 6 << BITS_FOR_OCCURRENCE_IDX;
// this actually returns a number, but it is always used as a string (object prop key)
function getKeyForAutomaticLookahead(ruleIdx, dslMethodIdx, occurrence) {
    return occurrence | dslMethodIdx | ruleIdx;
}
const BITS_START_FOR_ALT_IDX = 32 - BITS_FOR_ALT_IDX;
//# sourceMappingURL=keys.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/grammar/llk_lookahead.js





class LLkLookaheadStrategy {
    constructor(options) {
        var _a;
        this.maxLookahead =
            (_a = options === null || options === void 0 ? void 0 : options.maxLookahead) !== null && _a !== void 0 ? _a : DEFAULT_PARSER_CONFIG.maxLookahead;
    }
    validate(options) {
        const leftRecursionErrors = this.validateNoLeftRecursion(options.rules);
        if ((0,isEmpty/* default */.Z)(leftRecursionErrors)) {
            const emptyAltErrors = this.validateEmptyOrAlternatives(options.rules);
            const ambiguousAltsErrors = this.validateAmbiguousAlternationAlternatives(options.rules, this.maxLookahead);
            const emptyRepetitionErrors = this.validateSomeNonEmptyLookaheadPath(options.rules, this.maxLookahead);
            const allErrors = [
                ...leftRecursionErrors,
                ...emptyAltErrors,
                ...ambiguousAltsErrors,
                ...emptyRepetitionErrors,
            ];
            return allErrors;
        }
        return leftRecursionErrors;
    }
    validateNoLeftRecursion(rules) {
        return (0,flatMap/* default */.Z)(rules, (currTopRule) => validateNoLeftRecursion(currTopRule, currTopRule, defaultGrammarValidatorErrorProvider));
    }
    validateEmptyOrAlternatives(rules) {
        return (0,flatMap/* default */.Z)(rules, (currTopRule) => validateEmptyOrAlternative(currTopRule, defaultGrammarValidatorErrorProvider));
    }
    validateAmbiguousAlternationAlternatives(rules, maxLookahead) {
        return (0,flatMap/* default */.Z)(rules, (currTopRule) => validateAmbiguousAlternationAlternatives(currTopRule, maxLookahead, defaultGrammarValidatorErrorProvider));
    }
    validateSomeNonEmptyLookaheadPath(rules, maxLookahead) {
        return validateSomeNonEmptyLookaheadPath(rules, maxLookahead, defaultGrammarValidatorErrorProvider);
    }
    buildLookaheadForAlternation(options) {
        return buildLookaheadFuncForOr(options.prodOccurrence, options.rule, options.maxLookahead, options.hasPredicates, options.dynamicTokensEnabled, buildAlternativesLookAheadFunc);
    }
    buildLookaheadForOptional(options) {
        return buildLookaheadFuncForOptionalProd(options.prodOccurrence, options.rule, options.maxLookahead, options.dynamicTokensEnabled, getProdType(options.prodType), buildSingleAlternativeLookaheadFunction);
    }
}
//# sourceMappingURL=llk_lookahead.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/looksahead.js





/**
 * Trait responsible for the lookahead related utilities and optimizations.
 */
class LooksAhead {
    initLooksAhead(config) {
        this.dynamicTokensEnabled = (0,has/* default */.Z)(config, "dynamicTokensEnabled")
            ? config.dynamicTokensEnabled // assumes end user provides the correct config value/type
            : DEFAULT_PARSER_CONFIG.dynamicTokensEnabled;
        this.maxLookahead = (0,has/* default */.Z)(config, "maxLookahead")
            ? config.maxLookahead // assumes end user provides the correct config value/type
            : DEFAULT_PARSER_CONFIG.maxLookahead;
        this.lookaheadStrategy = (0,has/* default */.Z)(config, "lookaheadStrategy")
            ? config.lookaheadStrategy // assumes end user provides the correct config value/type
            : new LLkLookaheadStrategy({ maxLookahead: this.maxLookahead });
        this.lookAheadFuncsCache = new Map();
    }
    preComputeLookaheadFunctions(rules) {
        (0,forEach/* default */.Z)(rules, (currRule) => {
            this.TRACE_INIT(`${currRule.name} Rule Lookahead`, () => {
                const { alternation, repetition, option, repetitionMandatory, repetitionMandatoryWithSeparator, repetitionWithSeparator, } = collectMethods(currRule);
                (0,forEach/* default */.Z)(alternation, (currProd) => {
                    const prodIdx = currProd.idx === 0 ? "" : currProd.idx;
                    this.TRACE_INIT(`${getProductionDslName(currProd)}${prodIdx}`, () => {
                        const laFunc = this.lookaheadStrategy.buildLookaheadForAlternation({
                            prodOccurrence: currProd.idx,
                            rule: currRule,
                            maxLookahead: currProd.maxLookahead || this.maxLookahead,
                            hasPredicates: currProd.hasPredicates,
                            dynamicTokensEnabled: this.dynamicTokensEnabled,
                        });
                        const key = getKeyForAutomaticLookahead(this.fullRuleNameToShort[currRule.name], OR_IDX, currProd.idx);
                        this.setLaFuncCache(key, laFunc);
                    });
                });
                (0,forEach/* default */.Z)(repetition, (currProd) => {
                    this.computeLookaheadFunc(currRule, currProd.idx, MANY_IDX, "Repetition", currProd.maxLookahead, getProductionDslName(currProd));
                });
                (0,forEach/* default */.Z)(option, (currProd) => {
                    this.computeLookaheadFunc(currRule, currProd.idx, OPTION_IDX, "Option", currProd.maxLookahead, getProductionDslName(currProd));
                });
                (0,forEach/* default */.Z)(repetitionMandatory, (currProd) => {
                    this.computeLookaheadFunc(currRule, currProd.idx, AT_LEAST_ONE_IDX, "RepetitionMandatory", currProd.maxLookahead, getProductionDslName(currProd));
                });
                (0,forEach/* default */.Z)(repetitionMandatoryWithSeparator, (currProd) => {
                    this.computeLookaheadFunc(currRule, currProd.idx, AT_LEAST_ONE_SEP_IDX, "RepetitionMandatoryWithSeparator", currProd.maxLookahead, getProductionDslName(currProd));
                });
                (0,forEach/* default */.Z)(repetitionWithSeparator, (currProd) => {
                    this.computeLookaheadFunc(currRule, currProd.idx, MANY_SEP_IDX, "RepetitionWithSeparator", currProd.maxLookahead, getProductionDslName(currProd));
                });
            });
        });
    }
    computeLookaheadFunc(rule, prodOccurrence, prodKey, prodType, prodMaxLookahead, dslMethodName) {
        this.TRACE_INIT(`${dslMethodName}${prodOccurrence === 0 ? "" : prodOccurrence}`, () => {
            const laFunc = this.lookaheadStrategy.buildLookaheadForOptional({
                prodOccurrence,
                rule,
                maxLookahead: prodMaxLookahead || this.maxLookahead,
                dynamicTokensEnabled: this.dynamicTokensEnabled,
                prodType,
            });
            const key = getKeyForAutomaticLookahead(this.fullRuleNameToShort[rule.name], prodKey, prodOccurrence);
            this.setLaFuncCache(key, laFunc);
        });
    }
    // this actually returns a number, but it is always used as a string (object prop key)
    getKeyForAutomaticLookahead(dslMethodIdx, occurrence) {
        const currRuleShortName = this.getLastExplicitRuleShortName();
        return getKeyForAutomaticLookahead(currRuleShortName, dslMethodIdx, occurrence);
    }
    getLaFuncFromCache(key) {
        return this.lookAheadFuncsCache.get(key);
    }
    /* istanbul ignore next */
    setLaFuncCache(key, value) {
        this.lookAheadFuncsCache.set(key, value);
    }
}
class DslMethodsCollectorVisitor extends GAstVisitor {
    constructor() {
        super(...arguments);
        this.dslMethods = {
            option: [],
            alternation: [],
            repetition: [],
            repetitionWithSeparator: [],
            repetitionMandatory: [],
            repetitionMandatoryWithSeparator: [],
        };
    }
    reset() {
        this.dslMethods = {
            option: [],
            alternation: [],
            repetition: [],
            repetitionWithSeparator: [],
            repetitionMandatory: [],
            repetitionMandatoryWithSeparator: [],
        };
    }
    visitOption(option) {
        this.dslMethods.option.push(option);
    }
    visitRepetitionWithSeparator(manySep) {
        this.dslMethods.repetitionWithSeparator.push(manySep);
    }
    visitRepetitionMandatory(atLeastOne) {
        this.dslMethods.repetitionMandatory.push(atLeastOne);
    }
    visitRepetitionMandatoryWithSeparator(atLeastOneSep) {
        this.dslMethods.repetitionMandatoryWithSeparator.push(atLeastOneSep);
    }
    visitRepetition(many) {
        this.dslMethods.repetition.push(many);
    }
    visitAlternation(or) {
        this.dslMethods.alternation.push(or);
    }
}
const collectorVisitor = new DslMethodsCollectorVisitor();
function collectMethods(rule) {
    collectorVisitor.reset();
    rule.accept(collectorVisitor);
    const dslMethods = collectorVisitor.dslMethods;
    // avoid uncleaned references
    collectorVisitor.reset();
    return dslMethods;
}
//# sourceMappingURL=looksahead.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/cst/cst.js
/**
 * This nodeLocation tracking is not efficient and should only be used
 * when error recovery is enabled or the Token Vector contains virtual Tokens
 * (e.g, Python Indent/Outdent)
 * As it executes the calculation for every single terminal/nonTerminal
 * and does not rely on the fact the token vector is **sorted**
 */
function setNodeLocationOnlyOffset(currNodeLocation, newLocationInfo) {
    // First (valid) update for this cst node
    if (isNaN(currNodeLocation.startOffset) === true) {
        // assumption1: Token location information is either NaN or a valid number
        // assumption2: Token location information is fully valid if it exist
        // (both start/end offsets exist and are numbers).
        currNodeLocation.startOffset = newLocationInfo.startOffset;
        currNodeLocation.endOffset = newLocationInfo.endOffset;
    }
    // Once the startOffset has been updated with a valid number it should never receive
    // any farther updates as the Token vector is sorted.
    // We still have to check this this condition for every new possible location info
    // because with error recovery enabled we may encounter invalid tokens (NaN location props)
    else if (currNodeLocation.endOffset < newLocationInfo.endOffset === true) {
        currNodeLocation.endOffset = newLocationInfo.endOffset;
    }
}
/**
 * This nodeLocation tracking is not efficient and should only be used
 * when error recovery is enabled or the Token Vector contains virtual Tokens
 * (e.g, Python Indent/Outdent)
 * As it executes the calculation for every single terminal/nonTerminal
 * and does not rely on the fact the token vector is **sorted**
 */
function setNodeLocationFull(currNodeLocation, newLocationInfo) {
    // First (valid) update for this cst node
    if (isNaN(currNodeLocation.startOffset) === true) {
        // assumption1: Token location information is either NaN or a valid number
        // assumption2: Token location information is fully valid if it exist
        // (all start/end props exist and are numbers).
        currNodeLocation.startOffset = newLocationInfo.startOffset;
        currNodeLocation.startColumn = newLocationInfo.startColumn;
        currNodeLocation.startLine = newLocationInfo.startLine;
        currNodeLocation.endOffset = newLocationInfo.endOffset;
        currNodeLocation.endColumn = newLocationInfo.endColumn;
        currNodeLocation.endLine = newLocationInfo.endLine;
    }
    // Once the start props has been updated with a valid number it should never receive
    // any farther updates as the Token vector is sorted.
    // We still have to check this this condition for every new possible location info
    // because with error recovery enabled we may encounter invalid tokens (NaN location props)
    else if (currNodeLocation.endOffset < newLocationInfo.endOffset === true) {
        currNodeLocation.endOffset = newLocationInfo.endOffset;
        currNodeLocation.endColumn = newLocationInfo.endColumn;
        currNodeLocation.endLine = newLocationInfo.endLine;
    }
}
function addTerminalToCst(node, token, tokenTypeName) {
    if (node.children[tokenTypeName] === undefined) {
        node.children[tokenTypeName] = [token];
    }
    else {
        node.children[tokenTypeName].push(token);
    }
}
function addNoneTerminalToCst(node, ruleName, ruleResult) {
    if (node.children[ruleName] === undefined) {
        node.children[ruleName] = [ruleResult];
    }
    else {
        node.children[ruleName].push(ruleResult);
    }
}
//# sourceMappingURL=cst.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/lang/lang_extensions.js
const NAME = "name";
function defineNameProp(obj, nameValue) {
    Object.defineProperty(obj, NAME, {
        enumerable: false,
        configurable: true,
        writable: false,
        value: nameValue,
    });
}
//# sourceMappingURL=lang_extensions.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/cst/cst_visitor.js


function defaultVisit(ctx, param) {
    const childrenNames = (0,keys/* default */.Z)(ctx);
    const childrenNamesLength = childrenNames.length;
    for (let i = 0; i < childrenNamesLength; i++) {
        const currChildName = childrenNames[i];
        const currChildArray = ctx[currChildName];
        const currChildArrayLength = currChildArray.length;
        for (let j = 0; j < currChildArrayLength; j++) {
            const currChild = currChildArray[j];
            // distinction between Tokens Children and CstNode children
            if (currChild.tokenTypeIdx === undefined) {
                this[currChild.name](currChild.children, param);
            }
        }
    }
    // defaultVisit does not support generic out param
}
function createBaseSemanticVisitorConstructor(grammarName, ruleNames) {
    const derivedConstructor = function () { };
    // can be overwritten according to:
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Function/
    // name?redirectlocale=en-US&redirectslug=JavaScript%2FReference%2FGlobal_Objects%2FFunction%2Fname
    defineNameProp(derivedConstructor, grammarName + "BaseSemantics");
    const semanticProto = {
        visit: function (cstNode, param) {
            // enables writing more concise visitor methods when CstNode has only a single child
            if ((0,isArray/* default */.Z)(cstNode)) {
                // A CST Node's children dictionary can never have empty arrays as values
                // If a key is defined there will be at least one element in the corresponding value array.
                cstNode = cstNode[0];
            }
            // enables passing optional CstNodes concisely.
            if ((0,isUndefined/* default */.Z)(cstNode)) {
                return undefined;
            }
            return this[cstNode.name](cstNode.children, param);
        },
        validateVisitor: function () {
            const semanticDefinitionErrors = validateVisitor(this, ruleNames);
            if (!(0,isEmpty/* default */.Z)(semanticDefinitionErrors)) {
                const errorMessages = (0,map/* default */.Z)(semanticDefinitionErrors, (currDefError) => currDefError.msg);
                throw Error(`Errors Detected in CST Visitor <${this.constructor.name}>:\n\t` +
                    `${errorMessages.join("\n\n").replace(/\n/g, "\n\t")}`);
            }
        },
    };
    derivedConstructor.prototype = semanticProto;
    derivedConstructor.prototype.constructor = derivedConstructor;
    derivedConstructor._RULE_NAMES = ruleNames;
    return derivedConstructor;
}
function createBaseVisitorConstructorWithDefaults(grammarName, ruleNames, baseConstructor) {
    const derivedConstructor = function () { };
    // can be overwritten according to:
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Function/
    // name?redirectlocale=en-US&redirectslug=JavaScript%2FReference%2FGlobal_Objects%2FFunction%2Fname
    defineNameProp(derivedConstructor, grammarName + "BaseSemanticsWithDefaults");
    const withDefaultsProto = Object.create(baseConstructor.prototype);
    (0,forEach/* default */.Z)(ruleNames, (ruleName) => {
        withDefaultsProto[ruleName] = defaultVisit;
    });
    derivedConstructor.prototype = withDefaultsProto;
    derivedConstructor.prototype.constructor = derivedConstructor;
    return derivedConstructor;
}
var CstVisitorDefinitionError;
(function (CstVisitorDefinitionError) {
    CstVisitorDefinitionError[CstVisitorDefinitionError["REDUNDANT_METHOD"] = 0] = "REDUNDANT_METHOD";
    CstVisitorDefinitionError[CstVisitorDefinitionError["MISSING_METHOD"] = 1] = "MISSING_METHOD";
})(CstVisitorDefinitionError || (CstVisitorDefinitionError = {}));
function validateVisitor(visitorInstance, ruleNames) {
    const missingErrors = validateMissingCstMethods(visitorInstance, ruleNames);
    return missingErrors;
}
function validateMissingCstMethods(visitorInstance, ruleNames) {
    const missingRuleNames = (0,filter/* default */.Z)(ruleNames, (currRuleName) => {
        return (0,isFunction/* default */.Z)(visitorInstance[currRuleName]) === false;
    });
    const errors = (0,map/* default */.Z)(missingRuleNames, (currRuleName) => {
        return {
            msg: `Missing visitor method: <${currRuleName}> on ${(visitorInstance.constructor.name)} CST Visitor.`,
            type: CstVisitorDefinitionError.MISSING_METHOD,
            methodName: currRuleName,
        };
    });
    return lodash_es_compact(errors);
}
//# sourceMappingURL=cst_visitor.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/tree_builder.js




/**
 * This trait is responsible for the CST building logic.
 */
class TreeBuilder {
    initTreeBuilder(config) {
        this.CST_STACK = [];
        // outputCst is no longer exposed/defined in the pubic API
        this.outputCst = config.outputCst;
        this.nodeLocationTracking = (0,has/* default */.Z)(config, "nodeLocationTracking")
            ? config.nodeLocationTracking // assumes end user provides the correct config value/type
            : DEFAULT_PARSER_CONFIG.nodeLocationTracking;
        if (!this.outputCst) {
            this.cstInvocationStateUpdate = noop/* default */.Z;
            this.cstFinallyStateUpdate = noop/* default */.Z;
            this.cstPostTerminal = noop/* default */.Z;
            this.cstPostNonTerminal = noop/* default */.Z;
            this.cstPostRule = noop/* default */.Z;
        }
        else {
            if (/full/i.test(this.nodeLocationTracking)) {
                if (this.recoveryEnabled) {
                    this.setNodeLocationFromToken = setNodeLocationFull;
                    this.setNodeLocationFromNode = setNodeLocationFull;
                    this.cstPostRule = noop/* default */.Z;
                    this.setInitialNodeLocation = this.setInitialNodeLocationFullRecovery;
                }
                else {
                    this.setNodeLocationFromToken = noop/* default */.Z;
                    this.setNodeLocationFromNode = noop/* default */.Z;
                    this.cstPostRule = this.cstPostRuleFull;
                    this.setInitialNodeLocation = this.setInitialNodeLocationFullRegular;
                }
            }
            else if (/onlyOffset/i.test(this.nodeLocationTracking)) {
                if (this.recoveryEnabled) {
                    this.setNodeLocationFromToken = setNodeLocationOnlyOffset;
                    this.setNodeLocationFromNode = setNodeLocationOnlyOffset;
                    this.cstPostRule = noop/* default */.Z;
                    this.setInitialNodeLocation =
                        this.setInitialNodeLocationOnlyOffsetRecovery;
                }
                else {
                    this.setNodeLocationFromToken = noop/* default */.Z;
                    this.setNodeLocationFromNode = noop/* default */.Z;
                    this.cstPostRule = this.cstPostRuleOnlyOffset;
                    this.setInitialNodeLocation =
                        this.setInitialNodeLocationOnlyOffsetRegular;
                }
            }
            else if (/none/i.test(this.nodeLocationTracking)) {
                this.setNodeLocationFromToken = noop/* default */.Z;
                this.setNodeLocationFromNode = noop/* default */.Z;
                this.cstPostRule = noop/* default */.Z;
                this.setInitialNodeLocation = noop/* default */.Z;
            }
            else {
                throw Error(`Invalid <nodeLocationTracking> config option: "${config.nodeLocationTracking}"`);
            }
        }
    }
    setInitialNodeLocationOnlyOffsetRecovery(cstNode) {
        cstNode.location = {
            startOffset: NaN,
            endOffset: NaN,
        };
    }
    setInitialNodeLocationOnlyOffsetRegular(cstNode) {
        cstNode.location = {
            // without error recovery the starting Location of a new CstNode is guaranteed
            // To be the next Token's startOffset (for valid inputs).
            // For invalid inputs there won't be any CSTOutput so this potential
            // inaccuracy does not matter
            startOffset: this.LA(1).startOffset,
            endOffset: NaN,
        };
    }
    setInitialNodeLocationFullRecovery(cstNode) {
        cstNode.location = {
            startOffset: NaN,
            startLine: NaN,
            startColumn: NaN,
            endOffset: NaN,
            endLine: NaN,
            endColumn: NaN,
        };
    }
    /**
       *  @see setInitialNodeLocationOnlyOffsetRegular for explanation why this work
  
       * @param cstNode
       */
    setInitialNodeLocationFullRegular(cstNode) {
        const nextToken = this.LA(1);
        cstNode.location = {
            startOffset: nextToken.startOffset,
            startLine: nextToken.startLine,
            startColumn: nextToken.startColumn,
            endOffset: NaN,
            endLine: NaN,
            endColumn: NaN,
        };
    }
    cstInvocationStateUpdate(fullRuleName) {
        const cstNode = {
            name: fullRuleName,
            children: Object.create(null),
        };
        this.setInitialNodeLocation(cstNode);
        this.CST_STACK.push(cstNode);
    }
    cstFinallyStateUpdate() {
        this.CST_STACK.pop();
    }
    cstPostRuleFull(ruleCstNode) {
        // casts to `required<CstNodeLocation>` are safe because `cstPostRuleFull` should only be invoked when full location is enabled
        const prevToken = this.LA(0);
        const loc = ruleCstNode.location;
        // If this condition is true it means we consumed at least one Token
        // In this CstNode.
        if (loc.startOffset <= prevToken.startOffset === true) {
            loc.endOffset = prevToken.endOffset;
            loc.endLine = prevToken.endLine;
            loc.endColumn = prevToken.endColumn;
        }
        // "empty" CstNode edge case
        else {
            loc.startOffset = NaN;
            loc.startLine = NaN;
            loc.startColumn = NaN;
        }
    }
    cstPostRuleOnlyOffset(ruleCstNode) {
        const prevToken = this.LA(0);
        // `location' is not null because `cstPostRuleOnlyOffset` will only be invoked when location tracking is enabled.
        const loc = ruleCstNode.location;
        // If this condition is true it means we consumed at least one Token
        // In this CstNode.
        if (loc.startOffset <= prevToken.startOffset === true) {
            loc.endOffset = prevToken.endOffset;
        }
        // "empty" CstNode edge case
        else {
            loc.startOffset = NaN;
        }
    }
    cstPostTerminal(key, consumedToken) {
        const rootCst = this.CST_STACK[this.CST_STACK.length - 1];
        addTerminalToCst(rootCst, consumedToken, key);
        // This is only used when **both** error recovery and CST Output are enabled.
        this.setNodeLocationFromToken(rootCst.location, consumedToken);
    }
    cstPostNonTerminal(ruleCstResult, ruleName) {
        const preCstNode = this.CST_STACK[this.CST_STACK.length - 1];
        addNoneTerminalToCst(preCstNode, ruleName, ruleCstResult);
        // This is only used when **both** error recovery and CST Output are enabled.
        this.setNodeLocationFromNode(preCstNode.location, ruleCstResult.location);
    }
    getBaseCstVisitorConstructor() {
        if ((0,isUndefined/* default */.Z)(this.baseCstVisitorConstructor)) {
            const newBaseCstVisitorConstructor = createBaseSemanticVisitorConstructor(this.className, (0,keys/* default */.Z)(this.gastProductionsCache));
            this.baseCstVisitorConstructor = newBaseCstVisitorConstructor;
            return newBaseCstVisitorConstructor;
        }
        return this.baseCstVisitorConstructor;
    }
    getBaseCstVisitorConstructorWithDefaults() {
        if ((0,isUndefined/* default */.Z)(this.baseCstVisitorWithDefaultsConstructor)) {
            const newConstructor = createBaseVisitorConstructorWithDefaults(this.className, (0,keys/* default */.Z)(this.gastProductionsCache), this.getBaseCstVisitorConstructor());
            this.baseCstVisitorWithDefaultsConstructor = newConstructor;
            return newConstructor;
        }
        return this.baseCstVisitorWithDefaultsConstructor;
    }
    getLastExplicitRuleShortName() {
        const ruleStack = this.RULE_STACK;
        return ruleStack[ruleStack.length - 1];
    }
    getPreviousExplicitRuleShortName() {
        const ruleStack = this.RULE_STACK;
        return ruleStack[ruleStack.length - 2];
    }
    getLastExplicitRuleOccurrenceIndex() {
        const occurrenceStack = this.RULE_OCCURRENCE_STACK;
        return occurrenceStack[occurrenceStack.length - 1];
    }
}
//# sourceMappingURL=tree_builder.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/lexer_adapter.js

/**
 * Trait responsible abstracting over the interaction with Lexer output (Token vector).
 *
 * This could be generalized to support other kinds of lexers, e.g.
 * - Just in Time Lexing / Lexer-Less parsing.
 * - Streaming Lexer.
 */
class LexerAdapter {
    initLexerAdapter() {
        this.tokVector = [];
        this.tokVectorLength = 0;
        this.currIdx = -1;
    }
    set input(newInput) {
        // @ts-ignore - `this parameter` not supported in setters/getters
        //   - https://www.typescriptlang.org/docs/handbook/functions.html#this-parameters
        if (this.selfAnalysisDone !== true) {
            throw Error(`Missing <performSelfAnalysis> invocation at the end of the Parser's constructor.`);
        }
        // @ts-ignore - `this parameter` not supported in setters/getters
        //   - https://www.typescriptlang.org/docs/handbook/functions.html#this-parameters
        this.reset();
        this.tokVector = newInput;
        this.tokVectorLength = newInput.length;
    }
    get input() {
        return this.tokVector;
    }
    // skips a token and returns the next token
    SKIP_TOKEN() {
        if (this.currIdx <= this.tokVector.length - 2) {
            this.consumeToken();
            return this.LA(1);
        }
        else {
            return END_OF_FILE;
        }
    }
    // Lexer (accessing Token vector) related methods which can be overridden to implement lazy lexers
    // or lexers dependent on parser context.
    LA(howMuch) {
        const soughtIdx = this.currIdx + howMuch;
        if (soughtIdx < 0 || this.tokVectorLength <= soughtIdx) {
            return END_OF_FILE;
        }
        else {
            return this.tokVector[soughtIdx];
        }
    }
    consumeToken() {
        this.currIdx++;
    }
    exportLexerState() {
        return this.currIdx;
    }
    importLexerState(newState) {
        this.currIdx = newState;
    }
    resetLexerState() {
        this.currIdx = -1;
    }
    moveToTerminatedState() {
        this.currIdx = this.tokVector.length - 1;
    }
    getLexerPosition() {
        return this.exportLexerState();
    }
}
//# sourceMappingURL=lexer_adapter.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/recognizer_api.js






/**
 * This trait is responsible for implementing the public API
 * for defining Chevrotain parsers, i.e:
 * - CONSUME
 * - RULE
 * - OPTION
 * - ...
 */
class RecognizerApi {
    ACTION(impl) {
        return impl.call(this);
    }
    consume(idx, tokType, options) {
        return this.consumeInternal(tokType, idx, options);
    }
    subrule(idx, ruleToCall, options) {
        return this.subruleInternal(ruleToCall, idx, options);
    }
    option(idx, actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, idx);
    }
    or(idx, altsOrOpts) {
        return this.orInternal(altsOrOpts, idx);
    }
    many(idx, actionORMethodDef) {
        return this.manyInternal(idx, actionORMethodDef);
    }
    atLeastOne(idx, actionORMethodDef) {
        return this.atLeastOneInternal(idx, actionORMethodDef);
    }
    CONSUME(tokType, options) {
        return this.consumeInternal(tokType, 0, options);
    }
    CONSUME1(tokType, options) {
        return this.consumeInternal(tokType, 1, options);
    }
    CONSUME2(tokType, options) {
        return this.consumeInternal(tokType, 2, options);
    }
    CONSUME3(tokType, options) {
        return this.consumeInternal(tokType, 3, options);
    }
    CONSUME4(tokType, options) {
        return this.consumeInternal(tokType, 4, options);
    }
    CONSUME5(tokType, options) {
        return this.consumeInternal(tokType, 5, options);
    }
    CONSUME6(tokType, options) {
        return this.consumeInternal(tokType, 6, options);
    }
    CONSUME7(tokType, options) {
        return this.consumeInternal(tokType, 7, options);
    }
    CONSUME8(tokType, options) {
        return this.consumeInternal(tokType, 8, options);
    }
    CONSUME9(tokType, options) {
        return this.consumeInternal(tokType, 9, options);
    }
    SUBRULE(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 0, options);
    }
    SUBRULE1(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 1, options);
    }
    SUBRULE2(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 2, options);
    }
    SUBRULE3(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 3, options);
    }
    SUBRULE4(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 4, options);
    }
    SUBRULE5(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 5, options);
    }
    SUBRULE6(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 6, options);
    }
    SUBRULE7(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 7, options);
    }
    SUBRULE8(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 8, options);
    }
    SUBRULE9(ruleToCall, options) {
        return this.subruleInternal(ruleToCall, 9, options);
    }
    OPTION(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 0);
    }
    OPTION1(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 1);
    }
    OPTION2(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 2);
    }
    OPTION3(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 3);
    }
    OPTION4(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 4);
    }
    OPTION5(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 5);
    }
    OPTION6(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 6);
    }
    OPTION7(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 7);
    }
    OPTION8(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 8);
    }
    OPTION9(actionORMethodDef) {
        return this.optionInternal(actionORMethodDef, 9);
    }
    OR(altsOrOpts) {
        return this.orInternal(altsOrOpts, 0);
    }
    OR1(altsOrOpts) {
        return this.orInternal(altsOrOpts, 1);
    }
    OR2(altsOrOpts) {
        return this.orInternal(altsOrOpts, 2);
    }
    OR3(altsOrOpts) {
        return this.orInternal(altsOrOpts, 3);
    }
    OR4(altsOrOpts) {
        return this.orInternal(altsOrOpts, 4);
    }
    OR5(altsOrOpts) {
        return this.orInternal(altsOrOpts, 5);
    }
    OR6(altsOrOpts) {
        return this.orInternal(altsOrOpts, 6);
    }
    OR7(altsOrOpts) {
        return this.orInternal(altsOrOpts, 7);
    }
    OR8(altsOrOpts) {
        return this.orInternal(altsOrOpts, 8);
    }
    OR9(altsOrOpts) {
        return this.orInternal(altsOrOpts, 9);
    }
    MANY(actionORMethodDef) {
        this.manyInternal(0, actionORMethodDef);
    }
    MANY1(actionORMethodDef) {
        this.manyInternal(1, actionORMethodDef);
    }
    MANY2(actionORMethodDef) {
        this.manyInternal(2, actionORMethodDef);
    }
    MANY3(actionORMethodDef) {
        this.manyInternal(3, actionORMethodDef);
    }
    MANY4(actionORMethodDef) {
        this.manyInternal(4, actionORMethodDef);
    }
    MANY5(actionORMethodDef) {
        this.manyInternal(5, actionORMethodDef);
    }
    MANY6(actionORMethodDef) {
        this.manyInternal(6, actionORMethodDef);
    }
    MANY7(actionORMethodDef) {
        this.manyInternal(7, actionORMethodDef);
    }
    MANY8(actionORMethodDef) {
        this.manyInternal(8, actionORMethodDef);
    }
    MANY9(actionORMethodDef) {
        this.manyInternal(9, actionORMethodDef);
    }
    MANY_SEP(options) {
        this.manySepFirstInternal(0, options);
    }
    MANY_SEP1(options) {
        this.manySepFirstInternal(1, options);
    }
    MANY_SEP2(options) {
        this.manySepFirstInternal(2, options);
    }
    MANY_SEP3(options) {
        this.manySepFirstInternal(3, options);
    }
    MANY_SEP4(options) {
        this.manySepFirstInternal(4, options);
    }
    MANY_SEP5(options) {
        this.manySepFirstInternal(5, options);
    }
    MANY_SEP6(options) {
        this.manySepFirstInternal(6, options);
    }
    MANY_SEP7(options) {
        this.manySepFirstInternal(7, options);
    }
    MANY_SEP8(options) {
        this.manySepFirstInternal(8, options);
    }
    MANY_SEP9(options) {
        this.manySepFirstInternal(9, options);
    }
    AT_LEAST_ONE(actionORMethodDef) {
        this.atLeastOneInternal(0, actionORMethodDef);
    }
    AT_LEAST_ONE1(actionORMethodDef) {
        return this.atLeastOneInternal(1, actionORMethodDef);
    }
    AT_LEAST_ONE2(actionORMethodDef) {
        this.atLeastOneInternal(2, actionORMethodDef);
    }
    AT_LEAST_ONE3(actionORMethodDef) {
        this.atLeastOneInternal(3, actionORMethodDef);
    }
    AT_LEAST_ONE4(actionORMethodDef) {
        this.atLeastOneInternal(4, actionORMethodDef);
    }
    AT_LEAST_ONE5(actionORMethodDef) {
        this.atLeastOneInternal(5, actionORMethodDef);
    }
    AT_LEAST_ONE6(actionORMethodDef) {
        this.atLeastOneInternal(6, actionORMethodDef);
    }
    AT_LEAST_ONE7(actionORMethodDef) {
        this.atLeastOneInternal(7, actionORMethodDef);
    }
    AT_LEAST_ONE8(actionORMethodDef) {
        this.atLeastOneInternal(8, actionORMethodDef);
    }
    AT_LEAST_ONE9(actionORMethodDef) {
        this.atLeastOneInternal(9, actionORMethodDef);
    }
    AT_LEAST_ONE_SEP(options) {
        this.atLeastOneSepFirstInternal(0, options);
    }
    AT_LEAST_ONE_SEP1(options) {
        this.atLeastOneSepFirstInternal(1, options);
    }
    AT_LEAST_ONE_SEP2(options) {
        this.atLeastOneSepFirstInternal(2, options);
    }
    AT_LEAST_ONE_SEP3(options) {
        this.atLeastOneSepFirstInternal(3, options);
    }
    AT_LEAST_ONE_SEP4(options) {
        this.atLeastOneSepFirstInternal(4, options);
    }
    AT_LEAST_ONE_SEP5(options) {
        this.atLeastOneSepFirstInternal(5, options);
    }
    AT_LEAST_ONE_SEP6(options) {
        this.atLeastOneSepFirstInternal(6, options);
    }
    AT_LEAST_ONE_SEP7(options) {
        this.atLeastOneSepFirstInternal(7, options);
    }
    AT_LEAST_ONE_SEP8(options) {
        this.atLeastOneSepFirstInternal(8, options);
    }
    AT_LEAST_ONE_SEP9(options) {
        this.atLeastOneSepFirstInternal(9, options);
    }
    RULE(name, implementation, config = DEFAULT_RULE_CONFIG) {
        if (lodash_es_includes(this.definedRulesNames, name)) {
            const errMsg = defaultGrammarValidatorErrorProvider.buildDuplicateRuleNameError({
                topLevelRule: name,
                grammarName: this.className,
            });
            const error = {
                message: errMsg,
                type: ParserDefinitionErrorType.DUPLICATE_RULE_NAME,
                ruleName: name,
            };
            this.definitionErrors.push(error);
        }
        this.definedRulesNames.push(name);
        const ruleImplementation = this.defineRule(name, implementation, config);
        this[name] = ruleImplementation;
        return ruleImplementation;
    }
    OVERRIDE_RULE(name, impl, config = DEFAULT_RULE_CONFIG) {
        const ruleErrors = validateRuleIsOverridden(name, this.definedRulesNames, this.className);
        this.definitionErrors = this.definitionErrors.concat(ruleErrors);
        const ruleImplementation = this.defineRule(name, impl, config);
        this[name] = ruleImplementation;
        return ruleImplementation;
    }
    BACKTRACK(grammarRule, args) {
        return function () {
            // save org state
            this.isBackTrackingStack.push(1);
            const orgState = this.saveRecogState();
            try {
                grammarRule.apply(this, args);
                // if no exception was thrown we have succeed parsing the rule.
                return true;
            }
            catch (e) {
                if (isRecognitionException(e)) {
                    return false;
                }
                else {
                    throw e;
                }
            }
            finally {
                this.reloadRecogState(orgState);
                this.isBackTrackingStack.pop();
            }
        };
    }
    // GAST export APIs
    getGAstProductions() {
        return this.gastProductionsCache;
    }
    getSerializedGastProductions() {
        return serializeGrammar((0,values/* default */.Z)(this.gastProductionsCache));
    }
}
//# sourceMappingURL=recognizer_api.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/recognizer_engine.js









/**
 * This trait is responsible for the runtime parsing engine
 * Used by the official API (recognizer_api.ts)
 */
class RecognizerEngine {
    initRecognizerEngine(tokenVocabulary, config) {
        this.className = this.constructor.name;
        // TODO: would using an ES6 Map or plain object be faster (CST building scenario)
        this.shortRuleNameToFull = {};
        this.fullRuleNameToShort = {};
        this.ruleShortNameIdx = 256;
        this.tokenMatcher = tokenStructuredMatcherNoCategories;
        this.subruleIdx = 0;
        this.definedRulesNames = [];
        this.tokensMap = {};
        this.isBackTrackingStack = [];
        this.RULE_STACK = [];
        this.RULE_OCCURRENCE_STACK = [];
        this.gastProductionsCache = {};
        if ((0,has/* default */.Z)(config, "serializedGrammar")) {
            throw Error("The Parser's configuration can no longer contain a <serializedGrammar> property.\n" +
                "\tSee: https://chevrotain.io/docs/changes/BREAKING_CHANGES.html#_6-0-0\n" +
                "\tFor Further details.");
        }
        if ((0,isArray/* default */.Z)(tokenVocabulary)) {
            // This only checks for Token vocabularies provided as arrays.
            // That is good enough because the main objective is to detect users of pre-V4.0 APIs
            // rather than all edge cases of empty Token vocabularies.
            if ((0,isEmpty/* default */.Z)(tokenVocabulary)) {
                throw Error("A Token Vocabulary cannot be empty.\n" +
                    "\tNote that the first argument for the parser constructor\n" +
                    "\tis no longer a Token vector (since v4.0).");
            }
            if (typeof tokenVocabulary[0].startOffset === "number") {
                throw Error("The Parser constructor no longer accepts a token vector as the first argument.\n" +
                    "\tSee: https://chevrotain.io/docs/changes/BREAKING_CHANGES.html#_4-0-0\n" +
                    "\tFor Further details.");
            }
        }
        if ((0,isArray/* default */.Z)(tokenVocabulary)) {
            this.tokensMap = (0,reduce/* default */.Z)(tokenVocabulary, (acc, tokType) => {
                acc[tokType.name] = tokType;
                return acc;
            }, {});
        }
        else if ((0,has/* default */.Z)(tokenVocabulary, "modes") &&
            lodash_es_every((0,flatten/* default */.Z)((0,values/* default */.Z)(tokenVocabulary.modes)), isTokenType)) {
            const allTokenTypes = (0,flatten/* default */.Z)((0,values/* default */.Z)(tokenVocabulary.modes));
            const uniqueTokens = lodash_es_uniq(allTokenTypes);
            this.tokensMap = (0,reduce/* default */.Z)(uniqueTokens, (acc, tokType) => {
                acc[tokType.name] = tokType;
                return acc;
            }, {});
        }
        else if ((0,isObject/* default */.Z)(tokenVocabulary)) {
            this.tokensMap = (0,lodash_es_clone/* default */.Z)(tokenVocabulary);
        }
        else {
            throw new Error("<tokensDictionary> argument must be An Array of Token constructors," +
                " A dictionary of Token constructors or an IMultiModeLexerDefinition");
        }
        // always add EOF to the tokenNames -> constructors map. it is useful to assure all the input has been
        // parsed with a clear error message ("expecting EOF but found ...")
        this.tokensMap["EOF"] = EOF;
        const allTokenTypes = (0,has/* default */.Z)(tokenVocabulary, "modes")
            ? (0,flatten/* default */.Z)((0,values/* default */.Z)(tokenVocabulary.modes))
            : (0,values/* default */.Z)(tokenVocabulary);
        const noTokenCategoriesUsed = lodash_es_every(allTokenTypes, (tokenConstructor) => (0,isEmpty/* default */.Z)(tokenConstructor.categoryMatches));
        this.tokenMatcher = noTokenCategoriesUsed
            ? tokenStructuredMatcherNoCategories
            : tokenStructuredMatcher;
        // Because ES2015+ syntax should be supported for creating Token classes
        // We cannot assume that the Token classes were created using the "extendToken" utilities
        // Therefore we must augment the Token classes both on Lexer initialization and on Parser initialization
        augmentTokenTypes((0,values/* default */.Z)(this.tokensMap));
    }
    defineRule(ruleName, impl, config) {
        if (this.selfAnalysisDone) {
            throw Error(`Grammar rule <${ruleName}> may not be defined after the 'performSelfAnalysis' method has been called'\n` +
                `Make sure that all grammar rule definitions are done before 'performSelfAnalysis' is called.`);
        }
        const resyncEnabled = (0,has/* default */.Z)(config, "resyncEnabled")
            ? config.resyncEnabled // assumes end user provides the correct config value/type
            : DEFAULT_RULE_CONFIG.resyncEnabled;
        const recoveryValueFunc = (0,has/* default */.Z)(config, "recoveryValueFunc")
            ? config.recoveryValueFunc // assumes end user provides the correct config value/type
            : DEFAULT_RULE_CONFIG.recoveryValueFunc;
        // performance optimization: Use small integers as keys for the longer human readable "full" rule names.
        // this greatly improves Map access time (as much as 8% for some performance benchmarks).
        const shortName = this.ruleShortNameIdx << (BITS_FOR_METHOD_TYPE + BITS_FOR_OCCURRENCE_IDX);
        this.ruleShortNameIdx++;
        this.shortRuleNameToFull[shortName] = ruleName;
        this.fullRuleNameToShort[ruleName] = shortName;
        let invokeRuleWithTry;
        // Micro optimization, only check the condition **once** on rule definition
        // instead of **every single** rule invocation.
        if (this.outputCst === true) {
            invokeRuleWithTry = function invokeRuleWithTry(...args) {
                try {
                    this.ruleInvocationStateUpdate(shortName, ruleName, this.subruleIdx);
                    impl.apply(this, args);
                    const cst = this.CST_STACK[this.CST_STACK.length - 1];
                    this.cstPostRule(cst);
                    return cst;
                }
                catch (e) {
                    return this.invokeRuleCatch(e, resyncEnabled, recoveryValueFunc);
                }
                finally {
                    this.ruleFinallyStateUpdate();
                }
            };
        }
        else {
            invokeRuleWithTry = function invokeRuleWithTryCst(...args) {
                try {
                    this.ruleInvocationStateUpdate(shortName, ruleName, this.subruleIdx);
                    return impl.apply(this, args);
                }
                catch (e) {
                    return this.invokeRuleCatch(e, resyncEnabled, recoveryValueFunc);
                }
                finally {
                    this.ruleFinallyStateUpdate();
                }
            };
        }
        const wrappedGrammarRule = Object.assign(invokeRuleWithTry, { ruleName, originalGrammarAction: impl });
        return wrappedGrammarRule;
    }
    invokeRuleCatch(e, resyncEnabledConfig, recoveryValueFunc) {
        const isFirstInvokedRule = this.RULE_STACK.length === 1;
        // note the reSync is always enabled for the first rule invocation, because we must always be able to
        // reSync with EOF and just output some INVALID ParseTree
        // during backtracking reSync recovery is disabled, otherwise we can't be certain the backtracking
        // path is really the most valid one
        const reSyncEnabled = resyncEnabledConfig && !this.isBackTracking() && this.recoveryEnabled;
        if (isRecognitionException(e)) {
            const recogError = e;
            if (reSyncEnabled) {
                const reSyncTokType = this.findReSyncTokenType();
                if (this.isInCurrentRuleReSyncSet(reSyncTokType)) {
                    recogError.resyncedTokens = this.reSyncTo(reSyncTokType);
                    if (this.outputCst) {
                        const partialCstResult = this.CST_STACK[this.CST_STACK.length - 1];
                        partialCstResult.recoveredNode = true;
                        return partialCstResult;
                    }
                    else {
                        return recoveryValueFunc(e);
                    }
                }
                else {
                    if (this.outputCst) {
                        const partialCstResult = this.CST_STACK[this.CST_STACK.length - 1];
                        partialCstResult.recoveredNode = true;
                        recogError.partialCstResult = partialCstResult;
                    }
                    // to be handled Further up the call stack
                    throw recogError;
                }
            }
            else if (isFirstInvokedRule) {
                // otherwise a Redundant input error will be created as well and we cannot guarantee that this is indeed the case
                this.moveToTerminatedState();
                // the parser should never throw one of its own errors outside its flow.
                // even if error recovery is disabled
                return recoveryValueFunc(e);
            }
            else {
                // to be recovered Further up the call stack
                throw recogError;
            }
        }
        else {
            // some other Error type which we don't know how to handle (for example a built in JavaScript Error)
            throw e;
        }
    }
    // Implementation of parsing DSL
    optionInternal(actionORMethodDef, occurrence) {
        const key = this.getKeyForAutomaticLookahead(OPTION_IDX, occurrence);
        return this.optionInternalLogic(actionORMethodDef, occurrence, key);
    }
    optionInternalLogic(actionORMethodDef, occurrence, key) {
        let lookAheadFunc = this.getLaFuncFromCache(key);
        let action;
        if (typeof actionORMethodDef !== "function") {
            action = actionORMethodDef.DEF;
            const predicate = actionORMethodDef.GATE;
            // predicate present
            if (predicate !== undefined) {
                const orgLookaheadFunction = lookAheadFunc;
                lookAheadFunc = () => {
                    return predicate.call(this) && orgLookaheadFunction.call(this);
                };
            }
        }
        else {
            action = actionORMethodDef;
        }
        if (lookAheadFunc.call(this) === true) {
            return action.call(this);
        }
        return undefined;
    }
    atLeastOneInternal(prodOccurrence, actionORMethodDef) {
        const laKey = this.getKeyForAutomaticLookahead(AT_LEAST_ONE_IDX, prodOccurrence);
        return this.atLeastOneInternalLogic(prodOccurrence, actionORMethodDef, laKey);
    }
    atLeastOneInternalLogic(prodOccurrence, actionORMethodDef, key) {
        let lookAheadFunc = this.getLaFuncFromCache(key);
        let action;
        if (typeof actionORMethodDef !== "function") {
            action = actionORMethodDef.DEF;
            const predicate = actionORMethodDef.GATE;
            // predicate present
            if (predicate !== undefined) {
                const orgLookaheadFunction = lookAheadFunc;
                lookAheadFunc = () => {
                    return predicate.call(this) && orgLookaheadFunction.call(this);
                };
            }
        }
        else {
            action = actionORMethodDef;
        }
        if (lookAheadFunc.call(this) === true) {
            let notStuck = this.doSingleRepetition(action);
            while (lookAheadFunc.call(this) === true &&
                notStuck === true) {
                notStuck = this.doSingleRepetition(action);
            }
        }
        else {
            throw this.raiseEarlyExitException(prodOccurrence, PROD_TYPE.REPETITION_MANDATORY, actionORMethodDef.ERR_MSG);
        }
        // note that while it may seem that this can cause an error because by using a recursive call to
        // AT_LEAST_ONE we change the grammar to AT_LEAST_TWO, AT_LEAST_THREE ... , the possible recursive call
        // from the tryInRepetitionRecovery(...) will only happen IFF there really are TWO/THREE/.... items.
        // Performance optimization: "attemptInRepetitionRecovery" will be defined as NOOP unless recovery is enabled
        this.attemptInRepetitionRecovery(this.atLeastOneInternal, [prodOccurrence, actionORMethodDef], lookAheadFunc, AT_LEAST_ONE_IDX, prodOccurrence, NextTerminalAfterAtLeastOneWalker);
    }
    atLeastOneSepFirstInternal(prodOccurrence, options) {
        const laKey = this.getKeyForAutomaticLookahead(AT_LEAST_ONE_SEP_IDX, prodOccurrence);
        this.atLeastOneSepFirstInternalLogic(prodOccurrence, options, laKey);
    }
    atLeastOneSepFirstInternalLogic(prodOccurrence, options, key) {
        const action = options.DEF;
        const separator = options.SEP;
        const firstIterationLookaheadFunc = this.getLaFuncFromCache(key);
        // 1st iteration
        if (firstIterationLookaheadFunc.call(this) === true) {
            action.call(this);
            //  TODO: Optimization can move this function construction into "attemptInRepetitionRecovery"
            //  because it is only needed in error recovery scenarios.
            const separatorLookAheadFunc = () => {
                return this.tokenMatcher(this.LA(1), separator);
            };
            // 2nd..nth iterations
            while (this.tokenMatcher(this.LA(1), separator) === true) {
                // note that this CONSUME will never enter recovery because
                // the separatorLookAheadFunc checks that the separator really does exist.
                this.CONSUME(separator);
                // No need for checking infinite loop here due to consuming the separator.
                action.call(this);
            }
            // Performance optimization: "attemptInRepetitionRecovery" will be defined as NOOP unless recovery is enabled
            this.attemptInRepetitionRecovery(this.repetitionSepSecondInternal, [
                prodOccurrence,
                separator,
                separatorLookAheadFunc,
                action,
                NextTerminalAfterAtLeastOneSepWalker,
            ], separatorLookAheadFunc, AT_LEAST_ONE_SEP_IDX, prodOccurrence, NextTerminalAfterAtLeastOneSepWalker);
        }
        else {
            throw this.raiseEarlyExitException(prodOccurrence, PROD_TYPE.REPETITION_MANDATORY_WITH_SEPARATOR, options.ERR_MSG);
        }
    }
    manyInternal(prodOccurrence, actionORMethodDef) {
        const laKey = this.getKeyForAutomaticLookahead(MANY_IDX, prodOccurrence);
        return this.manyInternalLogic(prodOccurrence, actionORMethodDef, laKey);
    }
    manyInternalLogic(prodOccurrence, actionORMethodDef, key) {
        let lookaheadFunction = this.getLaFuncFromCache(key);
        let action;
        if (typeof actionORMethodDef !== "function") {
            action = actionORMethodDef.DEF;
            const predicate = actionORMethodDef.GATE;
            // predicate present
            if (predicate !== undefined) {
                const orgLookaheadFunction = lookaheadFunction;
                lookaheadFunction = () => {
                    return predicate.call(this) && orgLookaheadFunction.call(this);
                };
            }
        }
        else {
            action = actionORMethodDef;
        }
        let notStuck = true;
        while (lookaheadFunction.call(this) === true && notStuck === true) {
            notStuck = this.doSingleRepetition(action);
        }
        // Performance optimization: "attemptInRepetitionRecovery" will be defined as NOOP unless recovery is enabled
        this.attemptInRepetitionRecovery(this.manyInternal, [prodOccurrence, actionORMethodDef], lookaheadFunction, MANY_IDX, prodOccurrence, NextTerminalAfterManyWalker, 
        // The notStuck parameter is only relevant when "attemptInRepetitionRecovery"
        // is invoked from manyInternal, in the MANY_SEP case and AT_LEAST_ONE[_SEP]
        // An infinite loop cannot occur as:
        // - Either the lookahead is guaranteed to consume something (Single Token Separator)
        // - AT_LEAST_ONE by definition is guaranteed to consume something (or error out).
        notStuck);
    }
    manySepFirstInternal(prodOccurrence, options) {
        const laKey = this.getKeyForAutomaticLookahead(MANY_SEP_IDX, prodOccurrence);
        this.manySepFirstInternalLogic(prodOccurrence, options, laKey);
    }
    manySepFirstInternalLogic(prodOccurrence, options, key) {
        const action = options.DEF;
        const separator = options.SEP;
        const firstIterationLaFunc = this.getLaFuncFromCache(key);
        // 1st iteration
        if (firstIterationLaFunc.call(this) === true) {
            action.call(this);
            const separatorLookAheadFunc = () => {
                return this.tokenMatcher(this.LA(1), separator);
            };
            // 2nd..nth iterations
            while (this.tokenMatcher(this.LA(1), separator) === true) {
                // note that this CONSUME will never enter recovery because
                // the separatorLookAheadFunc checks that the separator really does exist.
                this.CONSUME(separator);
                // No need for checking infinite loop here due to consuming the separator.
                action.call(this);
            }
            // Performance optimization: "attemptInRepetitionRecovery" will be defined as NOOP unless recovery is enabled
            this.attemptInRepetitionRecovery(this.repetitionSepSecondInternal, [
                prodOccurrence,
                separator,
                separatorLookAheadFunc,
                action,
                NextTerminalAfterManySepWalker,
            ], separatorLookAheadFunc, MANY_SEP_IDX, prodOccurrence, NextTerminalAfterManySepWalker);
        }
    }
    repetitionSepSecondInternal(prodOccurrence, separator, separatorLookAheadFunc, action, nextTerminalAfterWalker) {
        while (separatorLookAheadFunc()) {
            // note that this CONSUME will never enter recovery because
            // the separatorLookAheadFunc checks that the separator really does exist.
            this.CONSUME(separator);
            action.call(this);
        }
        // we can only arrive to this function after an error
        // has occurred (hence the name 'second') so the following
        // IF will always be entered, its possible to remove it...
        // however it is kept to avoid confusion and be consistent.
        // Performance optimization: "attemptInRepetitionRecovery" will be defined as NOOP unless recovery is enabled
        /* istanbul ignore else */
        this.attemptInRepetitionRecovery(this.repetitionSepSecondInternal, [
            prodOccurrence,
            separator,
            separatorLookAheadFunc,
            action,
            nextTerminalAfterWalker,
        ], separatorLookAheadFunc, AT_LEAST_ONE_SEP_IDX, prodOccurrence, nextTerminalAfterWalker);
    }
    doSingleRepetition(action) {
        const beforeIteration = this.getLexerPosition();
        action.call(this);
        const afterIteration = this.getLexerPosition();
        // This boolean will indicate if this repetition progressed
        // or if we are "stuck" (potential infinite loop in the repetition).
        return afterIteration > beforeIteration;
    }
    orInternal(altsOrOpts, occurrence) {
        const laKey = this.getKeyForAutomaticLookahead(OR_IDX, occurrence);
        const alts = (0,isArray/* default */.Z)(altsOrOpts) ? altsOrOpts : altsOrOpts.DEF;
        const laFunc = this.getLaFuncFromCache(laKey);
        const altIdxToTake = laFunc.call(this, alts);
        if (altIdxToTake !== undefined) {
            const chosenAlternative = alts[altIdxToTake];
            return chosenAlternative.ALT.call(this);
        }
        this.raiseNoAltException(occurrence, altsOrOpts.ERR_MSG);
    }
    ruleFinallyStateUpdate() {
        this.RULE_STACK.pop();
        this.RULE_OCCURRENCE_STACK.pop();
        // NOOP when cst is disabled
        this.cstFinallyStateUpdate();
        if (this.RULE_STACK.length === 0 && this.isAtEndOfInput() === false) {
            const firstRedundantTok = this.LA(1);
            const errMsg = this.errorMessageProvider.buildNotAllInputParsedMessage({
                firstRedundant: firstRedundantTok,
                ruleName: this.getCurrRuleFullName(),
            });
            this.SAVE_ERROR(new NotAllInputParsedException(errMsg, firstRedundantTok));
        }
    }
    subruleInternal(ruleToCall, idx, options) {
        let ruleResult;
        try {
            const args = options !== undefined ? options.ARGS : undefined;
            this.subruleIdx = idx;
            ruleResult = ruleToCall.apply(this, args);
            this.cstPostNonTerminal(ruleResult, options !== undefined && options.LABEL !== undefined
                ? options.LABEL
                : ruleToCall.ruleName);
            return ruleResult;
        }
        catch (e) {
            throw this.subruleInternalError(e, options, ruleToCall.ruleName);
        }
    }
    subruleInternalError(e, options, ruleName) {
        if (isRecognitionException(e) && e.partialCstResult !== undefined) {
            this.cstPostNonTerminal(e.partialCstResult, options !== undefined && options.LABEL !== undefined
                ? options.LABEL
                : ruleName);
            delete e.partialCstResult;
        }
        throw e;
    }
    consumeInternal(tokType, idx, options) {
        let consumedToken;
        try {
            const nextToken = this.LA(1);
            if (this.tokenMatcher(nextToken, tokType) === true) {
                this.consumeToken();
                consumedToken = nextToken;
            }
            else {
                this.consumeInternalError(tokType, nextToken, options);
            }
        }
        catch (eFromConsumption) {
            consumedToken = this.consumeInternalRecovery(tokType, idx, eFromConsumption);
        }
        this.cstPostTerminal(options !== undefined && options.LABEL !== undefined
            ? options.LABEL
            : tokType.name, consumedToken);
        return consumedToken;
    }
    consumeInternalError(tokType, nextToken, options) {
        let msg;
        const previousToken = this.LA(0);
        if (options !== undefined && options.ERR_MSG) {
            msg = options.ERR_MSG;
        }
        else {
            msg = this.errorMessageProvider.buildMismatchTokenMessage({
                expected: tokType,
                actual: nextToken,
                previous: previousToken,
                ruleName: this.getCurrRuleFullName(),
            });
        }
        throw this.SAVE_ERROR(new MismatchedTokenException(msg, nextToken, previousToken));
    }
    consumeInternalRecovery(tokType, idx, eFromConsumption) {
        // no recovery allowed during backtracking, otherwise backtracking may recover invalid syntax and accept it
        // but the original syntax could have been parsed successfully without any backtracking + recovery
        if (this.recoveryEnabled &&
            // TODO: more robust checking of the exception type. Perhaps Typescript extending expressions?
            eFromConsumption.name === "MismatchedTokenException" &&
            !this.isBackTracking()) {
            const follows = this.getFollowsForInRuleRecovery(tokType, idx);
            try {
                return this.tryInRuleRecovery(tokType, follows);
            }
            catch (eFromInRuleRecovery) {
                if (eFromInRuleRecovery.name === IN_RULE_RECOVERY_EXCEPTION) {
                    // failed in RuleRecovery.
                    // throw the original error in order to trigger reSync error recovery
                    throw eFromConsumption;
                }
                else {
                    throw eFromInRuleRecovery;
                }
            }
        }
        else {
            throw eFromConsumption;
        }
    }
    saveRecogState() {
        // errors is a getter which will clone the errors array
        const savedErrors = this.errors;
        const savedRuleStack = (0,lodash_es_clone/* default */.Z)(this.RULE_STACK);
        return {
            errors: savedErrors,
            lexerState: this.exportLexerState(),
            RULE_STACK: savedRuleStack,
            CST_STACK: this.CST_STACK,
        };
    }
    reloadRecogState(newState) {
        this.errors = newState.errors;
        this.importLexerState(newState.lexerState);
        this.RULE_STACK = newState.RULE_STACK;
    }
    ruleInvocationStateUpdate(shortName, fullName, idxInCallingRule) {
        this.RULE_OCCURRENCE_STACK.push(idxInCallingRule);
        this.RULE_STACK.push(shortName);
        // NOOP when cst is disabled
        this.cstInvocationStateUpdate(fullName);
    }
    isBackTracking() {
        return this.isBackTrackingStack.length !== 0;
    }
    getCurrRuleFullName() {
        const shortName = this.getLastExplicitRuleShortName();
        return this.shortRuleNameToFull[shortName];
    }
    shortRuleNameToFullName(shortName) {
        return this.shortRuleNameToFull[shortName];
    }
    isAtEndOfInput() {
        return this.tokenMatcher(this.LA(1), EOF);
    }
    reset() {
        this.resetLexerState();
        this.subruleIdx = 0;
        this.isBackTrackingStack = [];
        this.errors = [];
        this.RULE_STACK = [];
        // TODO: extract a specific reset for TreeBuilder trait
        this.CST_STACK = [];
        this.RULE_OCCURRENCE_STACK = [];
    }
}
//# sourceMappingURL=recognizer_engine.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/error_handler.js




/**
 * Trait responsible for runtime parsing errors.
 */
class ErrorHandler {
    initErrorHandler(config) {
        this._errors = [];
        this.errorMessageProvider = (0,has/* default */.Z)(config, "errorMessageProvider")
            ? config.errorMessageProvider // assumes end user provides the correct config value/type
            : DEFAULT_PARSER_CONFIG.errorMessageProvider;
    }
    SAVE_ERROR(error) {
        if (isRecognitionException(error)) {
            error.context = {
                ruleStack: this.getHumanReadableRuleStack(),
                ruleOccurrenceStack: (0,lodash_es_clone/* default */.Z)(this.RULE_OCCURRENCE_STACK),
            };
            this._errors.push(error);
            return error;
        }
        else {
            throw Error("Trying to save an Error which is not a RecognitionException");
        }
    }
    get errors() {
        return (0,lodash_es_clone/* default */.Z)(this._errors);
    }
    set errors(newErrors) {
        this._errors = newErrors;
    }
    // TODO: consider caching the error message computed information
    raiseEarlyExitException(occurrence, prodType, userDefinedErrMsg) {
        const ruleName = this.getCurrRuleFullName();
        const ruleGrammar = this.getGAstProductions()[ruleName];
        const lookAheadPathsPerAlternative = getLookaheadPathsForOptionalProd(occurrence, ruleGrammar, prodType, this.maxLookahead);
        const insideProdPaths = lookAheadPathsPerAlternative[0];
        const actualTokens = [];
        for (let i = 1; i <= this.maxLookahead; i++) {
            actualTokens.push(this.LA(i));
        }
        const msg = this.errorMessageProvider.buildEarlyExitMessage({
            expectedIterationPaths: insideProdPaths,
            actual: actualTokens,
            previous: this.LA(0),
            customUserDescription: userDefinedErrMsg,
            ruleName: ruleName,
        });
        throw this.SAVE_ERROR(new EarlyExitException(msg, this.LA(1), this.LA(0)));
    }
    // TODO: consider caching the error message computed information
    raiseNoAltException(occurrence, errMsgTypes) {
        const ruleName = this.getCurrRuleFullName();
        const ruleGrammar = this.getGAstProductions()[ruleName];
        // TODO: getLookaheadPathsForOr can be slow for large enough maxLookahead and certain grammars, consider caching ?
        const lookAheadPathsPerAlternative = getLookaheadPathsForOr(occurrence, ruleGrammar, this.maxLookahead);
        const actualTokens = [];
        for (let i = 1; i <= this.maxLookahead; i++) {
            actualTokens.push(this.LA(i));
        }
        const previousToken = this.LA(0);
        const errMsg = this.errorMessageProvider.buildNoViableAltMessage({
            expectedPathsPerAlt: lookAheadPathsPerAlternative,
            actual: actualTokens,
            previous: previousToken,
            customUserDescription: errMsgTypes,
            ruleName: this.getCurrRuleFullName(),
        });
        throw this.SAVE_ERROR(new NoViableAltException(errMsg, this.LA(1), previousToken));
    }
}
//# sourceMappingURL=error_handler.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/context_assist.js


class ContentAssist {
    initContentAssist() { }
    computeContentAssist(startRuleName, precedingInput) {
        const startRuleGast = this.gastProductionsCache[startRuleName];
        if ((0,isUndefined/* default */.Z)(startRuleGast)) {
            throw Error(`Rule ->${startRuleName}<- does not exist in this grammar.`);
        }
        return nextPossibleTokensAfter([startRuleGast], precedingInput, this.tokenMatcher, this.maxLookahead);
    }
    // TODO: should this be a member method or a utility? it does not have any state or usage of 'this'...
    // TODO: should this be more explicitly part of the public API?
    getNextPossibleTokenTypes(grammarPath) {
        const topRuleName = lodash_es_head(grammarPath.ruleStack);
        const gastProductions = this.getGAstProductions();
        const topProduction = gastProductions[topRuleName];
        const nextPossibleTokenTypes = new NextAfterTokenWalker(topProduction, grammarPath).startWalking();
        return nextPossibleTokenTypes;
    }
}
//# sourceMappingURL=context_assist.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/gast_recorder.js







const RECORDING_NULL_OBJECT = {
    description: "This Object indicates the Parser is during Recording Phase",
};
Object.freeze(RECORDING_NULL_OBJECT);
const HANDLE_SEPARATOR = true;
const MAX_METHOD_IDX = Math.pow(2, BITS_FOR_OCCURRENCE_IDX) - 1;
const RFT = createToken({ name: "RECORDING_PHASE_TOKEN", pattern: Lexer.NA });
augmentTokenTypes([RFT]);
const RECORDING_PHASE_TOKEN = createTokenInstance(RFT, "This IToken indicates the Parser is in Recording Phase\n\t" +
    "" +
    "See: https://chevrotain.io/docs/guide/internals.html#grammar-recording for details", 
// Using "-1" instead of NaN (as in EOF) because an actual number is less likely to
// cause errors if the output of LA or CONSUME would be (incorrectly) used during the recording phase.
-1, -1, -1, -1, -1, -1);
Object.freeze(RECORDING_PHASE_TOKEN);
const RECORDING_PHASE_CSTNODE = {
    name: "This CSTNode indicates the Parser is in Recording Phase\n\t" +
        "See: https://chevrotain.io/docs/guide/internals.html#grammar-recording for details",
    children: {},
};
/**
 * This trait handles the creation of the GAST structure for Chevrotain Grammars
 */
class GastRecorder {
    initGastRecorder(config) {
        this.recordingProdStack = [];
        this.RECORDING_PHASE = false;
    }
    enableRecording() {
        this.RECORDING_PHASE = true;
        this.TRACE_INIT("Enable Recording", () => {
            /**
             * Warning Dark Voodoo Magic upcoming!
             * We are "replacing" the public parsing DSL methods API
             * With **new** alternative implementations on the Parser **instance**
             *
             * So far this is the only way I've found to avoid performance regressions during parsing time.
             * - Approx 30% performance regression was measured on Chrome 75 Canary when attempting to replace the "internal"
             *   implementations directly instead.
             */
            for (let i = 0; i < 10; i++) {
                const idx = i > 0 ? i : "";
                this[`CONSUME${idx}`] = function (arg1, arg2) {
                    return this.consumeInternalRecord(arg1, i, arg2);
                };
                this[`SUBRULE${idx}`] = function (arg1, arg2) {
                    return this.subruleInternalRecord(arg1, i, arg2);
                };
                this[`OPTION${idx}`] = function (arg1) {
                    return this.optionInternalRecord(arg1, i);
                };
                this[`OR${idx}`] = function (arg1) {
                    return this.orInternalRecord(arg1, i);
                };
                this[`MANY${idx}`] = function (arg1) {
                    this.manyInternalRecord(i, arg1);
                };
                this[`MANY_SEP${idx}`] = function (arg1) {
                    this.manySepFirstInternalRecord(i, arg1);
                };
                this[`AT_LEAST_ONE${idx}`] = function (arg1) {
                    this.atLeastOneInternalRecord(i, arg1);
                };
                this[`AT_LEAST_ONE_SEP${idx}`] = function (arg1) {
                    this.atLeastOneSepFirstInternalRecord(i, arg1);
                };
            }
            // DSL methods with the idx(suffix) as an argument
            this[`consume`] = function (idx, arg1, arg2) {
                return this.consumeInternalRecord(arg1, idx, arg2);
            };
            this[`subrule`] = function (idx, arg1, arg2) {
                return this.subruleInternalRecord(arg1, idx, arg2);
            };
            this[`option`] = function (idx, arg1) {
                return this.optionInternalRecord(arg1, idx);
            };
            this[`or`] = function (idx, arg1) {
                return this.orInternalRecord(arg1, idx);
            };
            this[`many`] = function (idx, arg1) {
                this.manyInternalRecord(idx, arg1);
            };
            this[`atLeastOne`] = function (idx, arg1) {
                this.atLeastOneInternalRecord(idx, arg1);
            };
            this.ACTION = this.ACTION_RECORD;
            this.BACKTRACK = this.BACKTRACK_RECORD;
            this.LA = this.LA_RECORD;
        });
    }
    disableRecording() {
        this.RECORDING_PHASE = false;
        // By deleting these **instance** properties, any future invocation
        // will be deferred to the original methods on the **prototype** object
        // This seems to get rid of any incorrect optimizations that V8 may
        // do during the recording phase.
        this.TRACE_INIT("Deleting Recording methods", () => {
            const that = this;
            for (let i = 0; i < 10; i++) {
                const idx = i > 0 ? i : "";
                delete that[`CONSUME${idx}`];
                delete that[`SUBRULE${idx}`];
                delete that[`OPTION${idx}`];
                delete that[`OR${idx}`];
                delete that[`MANY${idx}`];
                delete that[`MANY_SEP${idx}`];
                delete that[`AT_LEAST_ONE${idx}`];
                delete that[`AT_LEAST_ONE_SEP${idx}`];
            }
            delete that[`consume`];
            delete that[`subrule`];
            delete that[`option`];
            delete that[`or`];
            delete that[`many`];
            delete that[`atLeastOne`];
            delete that.ACTION;
            delete that.BACKTRACK;
            delete that.LA;
        });
    }
    //   Parser methods are called inside an ACTION?
    //   Maybe try/catch/finally on ACTIONS while disabling the recorders state changes?
    // @ts-expect-error -- noop place holder
    ACTION_RECORD(impl) {
        // NO-OP during recording
    }
    // Executing backtracking logic will break our recording logic assumptions
    BACKTRACK_RECORD(grammarRule, args) {
        return () => true;
    }
    // LA is part of the official API and may be used for custom lookahead logic
    // by end users who may forget to wrap it in ACTION or inside a GATE
    LA_RECORD(howMuch) {
        // We cannot use the RECORD_PHASE_TOKEN here because someone may depend
        // On LA return EOF at the end of the input so an infinite loop may occur.
        return END_OF_FILE;
    }
    topLevelRuleRecord(name, def) {
        try {
            const newTopLevelRule = new Rule({ definition: [], name: name });
            newTopLevelRule.name = name;
            this.recordingProdStack.push(newTopLevelRule);
            def.call(this);
            this.recordingProdStack.pop();
            return newTopLevelRule;
        }
        catch (originalError) {
            if (originalError.KNOWN_RECORDER_ERROR !== true) {
                try {
                    originalError.message =
                        originalError.message +
                            '\n\t This error was thrown during the "grammar recording phase" For more info see:\n\t' +
                            "https://chevrotain.io/docs/guide/internals.html#grammar-recording";
                }
                catch (mutabilityError) {
                    // We may not be able to modify the original error object
                    throw originalError;
                }
            }
            throw originalError;
        }
    }
    // Implementation of parsing DSL
    optionInternalRecord(actionORMethodDef, occurrence) {
        return recordProd.call(this, Option, actionORMethodDef, occurrence);
    }
    atLeastOneInternalRecord(occurrence, actionORMethodDef) {
        recordProd.call(this, RepetitionMandatory, actionORMethodDef, occurrence);
    }
    atLeastOneSepFirstInternalRecord(occurrence, options) {
        recordProd.call(this, RepetitionMandatoryWithSeparator, options, occurrence, HANDLE_SEPARATOR);
    }
    manyInternalRecord(occurrence, actionORMethodDef) {
        recordProd.call(this, Repetition, actionORMethodDef, occurrence);
    }
    manySepFirstInternalRecord(occurrence, options) {
        recordProd.call(this, RepetitionWithSeparator, options, occurrence, HANDLE_SEPARATOR);
    }
    orInternalRecord(altsOrOpts, occurrence) {
        return recordOrProd.call(this, altsOrOpts, occurrence);
    }
    subruleInternalRecord(ruleToCall, occurrence, options) {
        assertMethodIdxIsValid(occurrence);
        if (!ruleToCall || (0,has/* default */.Z)(ruleToCall, "ruleName") === false) {
            const error = new Error(`<SUBRULE${getIdxSuffix(occurrence)}> argument is invalid` +
                ` expecting a Parser method reference but got: <${JSON.stringify(ruleToCall)}>` +
                `\n inside top level rule: <${this.recordingProdStack[0].name}>`);
            error.KNOWN_RECORDER_ERROR = true;
            throw error;
        }
        const prevProd = (0,last/* default */.Z)(this.recordingProdStack);
        const ruleName = ruleToCall.ruleName;
        const newNoneTerminal = new NonTerminal({
            idx: occurrence,
            nonTerminalName: ruleName,
            label: options === null || options === void 0 ? void 0 : options.LABEL,
            // The resolving of the `referencedRule` property will be done once all the Rule's GASTs have been created
            referencedRule: undefined,
        });
        prevProd.definition.push(newNoneTerminal);
        return this.outputCst
            ? RECORDING_PHASE_CSTNODE
            : RECORDING_NULL_OBJECT;
    }
    consumeInternalRecord(tokType, occurrence, options) {
        assertMethodIdxIsValid(occurrence);
        if (!hasShortKeyProperty(tokType)) {
            const error = new Error(`<CONSUME${getIdxSuffix(occurrence)}> argument is invalid` +
                ` expecting a TokenType reference but got: <${JSON.stringify(tokType)}>` +
                `\n inside top level rule: <${this.recordingProdStack[0].name}>`);
            error.KNOWN_RECORDER_ERROR = true;
            throw error;
        }
        const prevProd = (0,last/* default */.Z)(this.recordingProdStack);
        const newNoneTerminal = new Terminal({
            idx: occurrence,
            terminalType: tokType,
            label: options === null || options === void 0 ? void 0 : options.LABEL,
        });
        prevProd.definition.push(newNoneTerminal);
        return RECORDING_PHASE_TOKEN;
    }
}
function recordProd(prodConstructor, mainProdArg, occurrence, handleSep = false) {
    assertMethodIdxIsValid(occurrence);
    const prevProd = (0,last/* default */.Z)(this.recordingProdStack);
    const grammarAction = (0,isFunction/* default */.Z)(mainProdArg) ? mainProdArg : mainProdArg.DEF;
    const newProd = new prodConstructor({ definition: [], idx: occurrence });
    if (handleSep) {
        newProd.separator = mainProdArg.SEP;
    }
    if ((0,has/* default */.Z)(mainProdArg, "MAX_LOOKAHEAD")) {
        newProd.maxLookahead = mainProdArg.MAX_LOOKAHEAD;
    }
    this.recordingProdStack.push(newProd);
    grammarAction.call(this);
    prevProd.definition.push(newProd);
    this.recordingProdStack.pop();
    return RECORDING_NULL_OBJECT;
}
function recordOrProd(mainProdArg, occurrence) {
    assertMethodIdxIsValid(occurrence);
    const prevProd = (0,last/* default */.Z)(this.recordingProdStack);
    // Only an array of alternatives
    const hasOptions = (0,isArray/* default */.Z)(mainProdArg) === false;
    const alts = hasOptions === false ? mainProdArg : mainProdArg.DEF;
    const newOrProd = new Alternation({
        definition: [],
        idx: occurrence,
        ignoreAmbiguities: hasOptions && mainProdArg.IGNORE_AMBIGUITIES === true,
    });
    if ((0,has/* default */.Z)(mainProdArg, "MAX_LOOKAHEAD")) {
        newOrProd.maxLookahead = mainProdArg.MAX_LOOKAHEAD;
    }
    const hasPredicates = lodash_es_some(alts, (currAlt) => (0,isFunction/* default */.Z)(currAlt.GATE));
    newOrProd.hasPredicates = hasPredicates;
    prevProd.definition.push(newOrProd);
    (0,forEach/* default */.Z)(alts, (currAlt) => {
        const currAltFlat = new Alternative({ definition: [] });
        newOrProd.definition.push(currAltFlat);
        if ((0,has/* default */.Z)(currAlt, "IGNORE_AMBIGUITIES")) {
            currAltFlat.ignoreAmbiguities = currAlt.IGNORE_AMBIGUITIES; // assumes end user provides the correct config value/type
        }
        // **implicit** ignoreAmbiguities due to usage of gate
        else if ((0,has/* default */.Z)(currAlt, "GATE")) {
            currAltFlat.ignoreAmbiguities = true;
        }
        this.recordingProdStack.push(currAltFlat);
        currAlt.ALT.call(this);
        this.recordingProdStack.pop();
    });
    return RECORDING_NULL_OBJECT;
}
function getIdxSuffix(idx) {
    return idx === 0 ? "" : `${idx}`;
}
function assertMethodIdxIsValid(idx) {
    if (idx < 0 || idx > MAX_METHOD_IDX) {
        const error = new Error(
        // The stack trace will contain all the needed details
        `Invalid DSL Method idx value: <${idx}>\n\t` +
            `Idx value must be a none negative value smaller than ${MAX_METHOD_IDX + 1}`);
        error.KNOWN_RECORDER_ERROR = true;
        throw error;
    }
}
//# sourceMappingURL=gast_recorder.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/traits/perf_tracer.js



/**
 * Trait responsible for runtime parsing errors.
 */
class PerformanceTracer {
    initPerformanceTracer(config) {
        if ((0,has/* default */.Z)(config, "traceInitPerf")) {
            const userTraceInitPerf = config.traceInitPerf;
            const traceIsNumber = typeof userTraceInitPerf === "number";
            this.traceInitMaxIdent = traceIsNumber
                ? userTraceInitPerf
                : Infinity;
            this.traceInitPerf = traceIsNumber
                ? userTraceInitPerf > 0
                : userTraceInitPerf; // assumes end user provides the correct config value/type
        }
        else {
            this.traceInitMaxIdent = 0;
            this.traceInitPerf = DEFAULT_PARSER_CONFIG.traceInitPerf;
        }
        this.traceInitIndent = -1;
    }
    TRACE_INIT(phaseDesc, phaseImpl) {
        // No need to optimize this using NOOP pattern because
        // It is not called in a hot spot...
        if (this.traceInitPerf === true) {
            this.traceInitIndent++;
            const indent = new Array(this.traceInitIndent + 1).join("\t");
            if (this.traceInitIndent < this.traceInitMaxIdent) {
                console.log(`${indent}--> <${phaseDesc}>`);
            }
            const { time, value } = timer(phaseImpl);
            /* istanbul ignore next - Difficult to reproduce specific performance behavior (>10ms) in tests */
            const traceMethod = time > 10 ? console.warn : console.log;
            if (this.traceInitIndent < this.traceInitMaxIdent) {
                traceMethod(`${indent}<-- <${phaseDesc}> time: ${time}ms`);
            }
            this.traceInitIndent--;
            return value;
        }
        else {
            return phaseImpl();
        }
    }
}
//# sourceMappingURL=perf_tracer.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/utils/apply_mixins.js
function applyMixins(derivedCtor, baseCtors) {
    baseCtors.forEach((baseCtor) => {
        const baseProto = baseCtor.prototype;
        Object.getOwnPropertyNames(baseProto).forEach((propName) => {
            if (propName === "constructor") {
                return;
            }
            const basePropDescriptor = Object.getOwnPropertyDescriptor(baseProto, propName);
            // Handle Accessors
            if (basePropDescriptor &&
                (basePropDescriptor.get || basePropDescriptor.set)) {
                Object.defineProperty(derivedCtor.prototype, propName, basePropDescriptor);
            }
            else {
                derivedCtor.prototype[propName] = baseCtor.prototype[propName];
            }
        });
    });
}
//# sourceMappingURL=apply_mixins.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/parse/parser/parser.js


















const END_OF_FILE = createTokenInstance(EOF, "", NaN, NaN, NaN, NaN, NaN, NaN);
Object.freeze(END_OF_FILE);
const DEFAULT_PARSER_CONFIG = Object.freeze({
    recoveryEnabled: false,
    maxLookahead: 3,
    dynamicTokensEnabled: false,
    outputCst: true,
    errorMessageProvider: defaultParserErrorProvider,
    nodeLocationTracking: "none",
    traceInitPerf: false,
    skipValidations: false,
});
const DEFAULT_RULE_CONFIG = Object.freeze({
    recoveryValueFunc: () => undefined,
    resyncEnabled: true,
});
var ParserDefinitionErrorType;
(function (ParserDefinitionErrorType) {
    ParserDefinitionErrorType[ParserDefinitionErrorType["INVALID_RULE_NAME"] = 0] = "INVALID_RULE_NAME";
    ParserDefinitionErrorType[ParserDefinitionErrorType["DUPLICATE_RULE_NAME"] = 1] = "DUPLICATE_RULE_NAME";
    ParserDefinitionErrorType[ParserDefinitionErrorType["INVALID_RULE_OVERRIDE"] = 2] = "INVALID_RULE_OVERRIDE";
    ParserDefinitionErrorType[ParserDefinitionErrorType["DUPLICATE_PRODUCTIONS"] = 3] = "DUPLICATE_PRODUCTIONS";
    ParserDefinitionErrorType[ParserDefinitionErrorType["UNRESOLVED_SUBRULE_REF"] = 4] = "UNRESOLVED_SUBRULE_REF";
    ParserDefinitionErrorType[ParserDefinitionErrorType["LEFT_RECURSION"] = 5] = "LEFT_RECURSION";
    ParserDefinitionErrorType[ParserDefinitionErrorType["NONE_LAST_EMPTY_ALT"] = 6] = "NONE_LAST_EMPTY_ALT";
    ParserDefinitionErrorType[ParserDefinitionErrorType["AMBIGUOUS_ALTS"] = 7] = "AMBIGUOUS_ALTS";
    ParserDefinitionErrorType[ParserDefinitionErrorType["CONFLICT_TOKENS_RULES_NAMESPACE"] = 8] = "CONFLICT_TOKENS_RULES_NAMESPACE";
    ParserDefinitionErrorType[ParserDefinitionErrorType["INVALID_TOKEN_NAME"] = 9] = "INVALID_TOKEN_NAME";
    ParserDefinitionErrorType[ParserDefinitionErrorType["NO_NON_EMPTY_LOOKAHEAD"] = 10] = "NO_NON_EMPTY_LOOKAHEAD";
    ParserDefinitionErrorType[ParserDefinitionErrorType["AMBIGUOUS_PREFIX_ALTS"] = 11] = "AMBIGUOUS_PREFIX_ALTS";
    ParserDefinitionErrorType[ParserDefinitionErrorType["TOO_MANY_ALTS"] = 12] = "TOO_MANY_ALTS";
    ParserDefinitionErrorType[ParserDefinitionErrorType["CUSTOM_LOOKAHEAD_VALIDATION"] = 13] = "CUSTOM_LOOKAHEAD_VALIDATION";
})(ParserDefinitionErrorType || (ParserDefinitionErrorType = {}));
function EMPTY_ALT(value = undefined) {
    return function () {
        return value;
    };
}
class Parser {
    /**
     *  @deprecated use the **instance** method with the same name instead
     */
    static performSelfAnalysis(parserInstance) {
        throw Error("The **static** `performSelfAnalysis` method has been deprecated." +
            "\t\nUse the **instance** method with the same name instead.");
    }
    performSelfAnalysis() {
        this.TRACE_INIT("performSelfAnalysis", () => {
            let defErrorsMsgs;
            this.selfAnalysisDone = true;
            const className = this.className;
            this.TRACE_INIT("toFastProps", () => {
                // Without this voodoo magic the parser would be x3-x4 slower
                // It seems it is better to invoke `toFastProperties` **before**
                // Any manipulations of the `this` object done during the recording phase.
                toFastProperties(this);
            });
            this.TRACE_INIT("Grammar Recording", () => {
                try {
                    this.enableRecording();
                    // Building the GAST
                    (0,forEach/* default */.Z)(this.definedRulesNames, (currRuleName) => {
                        const wrappedRule = this[currRuleName];
                        const originalGrammarAction = wrappedRule["originalGrammarAction"];
                        let recordedRuleGast;
                        this.TRACE_INIT(`${currRuleName} Rule`, () => {
                            recordedRuleGast = this.topLevelRuleRecord(currRuleName, originalGrammarAction);
                        });
                        this.gastProductionsCache[currRuleName] = recordedRuleGast;
                    });
                }
                finally {
                    this.disableRecording();
                }
            });
            let resolverErrors = [];
            this.TRACE_INIT("Grammar Resolving", () => {
                resolverErrors = gast_resolver_public_resolveGrammar({
                    rules: (0,values/* default */.Z)(this.gastProductionsCache),
                });
                this.definitionErrors = this.definitionErrors.concat(resolverErrors);
            });
            this.TRACE_INIT("Grammar Validations", () => {
                // only perform additional grammar validations IFF no resolving errors have occurred.
                // as unresolved grammar may lead to unhandled runtime exceptions in the follow up validations.
                if ((0,isEmpty/* default */.Z)(resolverErrors) && this.skipValidations === false) {
                    const validationErrors = gast_resolver_public_validateGrammar({
                        rules: (0,values/* default */.Z)(this.gastProductionsCache),
                        tokenTypes: (0,values/* default */.Z)(this.tokensMap),
                        errMsgProvider: defaultGrammarValidatorErrorProvider,
                        grammarName: className,
                    });
                    const lookaheadValidationErrors = validateLookahead({
                        lookaheadStrategy: this.lookaheadStrategy,
                        rules: (0,values/* default */.Z)(this.gastProductionsCache),
                        tokenTypes: (0,values/* default */.Z)(this.tokensMap),
                        grammarName: className,
                    });
                    this.definitionErrors = this.definitionErrors.concat(validationErrors, lookaheadValidationErrors);
                }
            });
            // this analysis may fail if the grammar is not perfectly valid
            if ((0,isEmpty/* default */.Z)(this.definitionErrors)) {
                // The results of these computations are not needed unless error recovery is enabled.
                if (this.recoveryEnabled) {
                    this.TRACE_INIT("computeAllProdsFollows", () => {
                        const allFollows = computeAllProdsFollows((0,values/* default */.Z)(this.gastProductionsCache));
                        this.resyncFollows = allFollows;
                    });
                }
                this.TRACE_INIT("ComputeLookaheadFunctions", () => {
                    var _a, _b;
                    (_b = (_a = this.lookaheadStrategy).initialize) === null || _b === void 0 ? void 0 : _b.call(_a, {
                        rules: (0,values/* default */.Z)(this.gastProductionsCache),
                    });
                    this.preComputeLookaheadFunctions((0,values/* default */.Z)(this.gastProductionsCache));
                });
            }
            if (!Parser.DEFER_DEFINITION_ERRORS_HANDLING &&
                !(0,isEmpty/* default */.Z)(this.definitionErrors)) {
                defErrorsMsgs = (0,map/* default */.Z)(this.definitionErrors, (defError) => defError.message);
                throw new Error(`Parser Definition Errors detected:\n ${defErrorsMsgs.join("\n-------------------------------\n")}`);
            }
        });
    }
    constructor(tokenVocabulary, config) {
        this.definitionErrors = [];
        this.selfAnalysisDone = false;
        const that = this;
        that.initErrorHandler(config);
        that.initLexerAdapter();
        that.initLooksAhead(config);
        that.initRecognizerEngine(tokenVocabulary, config);
        that.initRecoverable(config);
        that.initTreeBuilder(config);
        that.initContentAssist();
        that.initGastRecorder(config);
        that.initPerformanceTracer(config);
        if ((0,has/* default */.Z)(config, "ignoredIssues")) {
            throw new Error("The <ignoredIssues> IParserConfig property has been deprecated.\n\t" +
                "Please use the <IGNORE_AMBIGUITIES> flag on the relevant DSL method instead.\n\t" +
                "See: https://chevrotain.io/docs/guide/resolving_grammar_errors.html#IGNORING_AMBIGUITIES\n\t" +
                "For further details.");
        }
        this.skipValidations = (0,has/* default */.Z)(config, "skipValidations")
            ? config.skipValidations // casting assumes the end user passing the correct type
            : DEFAULT_PARSER_CONFIG.skipValidations;
    }
}
// Set this flag to true if you don't want the Parser to throw error when problems in it's definition are detected.
// (normally during the parser's constructor).
// This is a design time flag, it will not affect the runtime error handling of the parser, just design time errors,
// for example: duplicate rule names, referencing an unresolved subrule, ect...
// This flag should not be enabled during normal usage, it is used in special situations, for example when
// needing to display the parser definition errors in some GUI(online playground).
Parser.DEFER_DEFINITION_ERRORS_HANDLING = false;
applyMixins(Parser, [
    Recoverable,
    LooksAhead,
    TreeBuilder,
    LexerAdapter,
    RecognizerEngine,
    RecognizerApi,
    ErrorHandler,
    ContentAssist,
    GastRecorder,
    PerformanceTracer,
]);
class CstParser extends (/* unused pure expression or super */ null && (Parser)) {
    constructor(tokenVocabulary, config = DEFAULT_PARSER_CONFIG) {
        const configClone = clone(config);
        configClone.outputCst = true;
        super(tokenVocabulary, configClone);
    }
}
class EmbeddedActionsParser extends Parser {
    constructor(tokenVocabulary, config = DEFAULT_PARSER_CONFIG) {
        const configClone = (0,lodash_es_clone/* default */.Z)(config);
        configClone.outputCst = false;
        super(tokenVocabulary, configClone);
    }
}
//# sourceMappingURL=parser.js.map
;// CONCATENATED MODULE: ../node_modules/@chevrotain/cst-dts-gen/lib/src/api.js


const defaultOptions = {
    includeVisitorInterface: true,
    visitorInterfaceName: "ICstNodeVisitor",
};
function generateCstDts(productions, options) {
    const effectiveOptions = Object.assign(Object.assign({}, defaultOptions), options);
    const model = buildModel(productions);
    return genDts(model, effectiveOptions);
}
//# sourceMappingURL=api.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain/lib/src/api.js
/* istanbul ignore file - tricky to import some things from this module during testing */
// semantic version



// Tokens utilities

// Lookahead


// Other Utilities



// grammar reflection API

// GAST Utilities


/* istanbul ignore next */
function clearCache() {
    console.warn("The clearCache function was 'soft' removed from the Chevrotain API." +
        "\n\t It performs no action other than printing this message." +
        "\n\t Please avoid using it as it will be completely removed in the future");
}

class api_Parser {
    constructor() {
        throw new Error("The Parser class has been deprecated, use CstParser or EmbeddedActionsParser instead.\t\n" +
            "See: https://chevrotain.io/docs/changes/BREAKING_CHANGES.html#_7-0-0");
    }
}
//# sourceMappingURL=api.js.map

/***/ }),

/***/ 73001:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Q: () => (/* binding */ createDefaultCoreModule),
  T: () => (/* binding */ createDefaultSharedCoreModule)
});

// EXTERNAL MODULE: ../node_modules/langium/lib/utils/cst-utils.js
var cst_utils = __webpack_require__(13871);
// EXTERNAL MODULE: ../node_modules/langium/lib/utils/grammar-utils.js
var grammar_utils = __webpack_require__(30447);
// EXTERNAL MODULE: ../node_modules/langium/lib/utils/regexp-utils.js
var regexp_utils = __webpack_require__(43078);
// EXTERNAL MODULE: ../node_modules/langium/lib/languages/generated/ast.js
var ast = __webpack_require__(34905);
;// CONCATENATED MODULE: ../node_modules/langium/lib/languages/grammar-config.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/




/**
 * Create the default grammar configuration (used by `createDefaultModule`). This can be overridden in a
 * language-specific module.
 */
function createGrammarConfig(services) {
    const rules = [];
    const grammar = services.Grammar;
    for (const rule of grammar.rules) {
        if ((0,ast/* isTerminalRule */.MS)(rule) && (0,grammar_utils/* isCommentTerminal */.md)(rule) && (0,regexp_utils/* isMultilineComment */.Rn)((0,grammar_utils/* terminalRegex */.s1)(rule))) {
            rules.push(rule.name);
        }
    }
    return {
        multilineCommentRules: rules,
        nameRegexp: cst_utils/* DefaultNameRegexp */.uz
    };
}
//# sourceMappingURL=grammar-config.js.map
// EXTERNAL MODULE: ../node_modules/chevrotain/lib/src/api.js + 67 modules
var api = __webpack_require__(34326);
// EXTERNAL MODULE: ../node_modules/lodash-es/map.js
var map = __webpack_require__(12930);
// EXTERNAL MODULE: ../node_modules/lodash-es/filter.js
var filter = __webpack_require__(11382);
;// CONCATENATED MODULE: ../node_modules/chevrotain-allstar/lib/atn.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/



function buildATNKey(rule, type, occurrence) {
    return `${rule.name}_${type}_${occurrence}`;
}
const ATN_INVALID_TYPE = 0;
const ATN_BASIC = 1;
const ATN_RULE_START = 2;
const ATN_PLUS_BLOCK_START = 4;
const ATN_STAR_BLOCK_START = 5;
// Currently unused as the ATN is not used for lexing
const ATN_TOKEN_START = 6;
const ATN_RULE_STOP = 7;
const ATN_BLOCK_END = 8;
const ATN_STAR_LOOP_BACK = 9;
const ATN_STAR_LOOP_ENTRY = 10;
const ATN_PLUS_LOOP_BACK = 11;
const ATN_LOOP_END = 12;
class AbstractTransition {
    constructor(target) {
        this.target = target;
    }
    isEpsilon() {
        return false;
    }
}
class AtomTransition extends AbstractTransition {
    constructor(target, tokenType) {
        super(target);
        this.tokenType = tokenType;
    }
}
class EpsilonTransition extends AbstractTransition {
    constructor(target) {
        super(target);
    }
    isEpsilon() {
        return true;
    }
}
class RuleTransition extends AbstractTransition {
    constructor(ruleStart, rule, followState) {
        super(ruleStart);
        this.rule = rule;
        this.followState = followState;
    }
    isEpsilon() {
        return true;
    }
}
function createATN(rules) {
    const atn = {
        decisionMap: {},
        decisionStates: [],
        ruleToStartState: new Map(),
        ruleToStopState: new Map(),
        states: []
    };
    createRuleStartAndStopATNStates(atn, rules);
    const ruleLength = rules.length;
    for (let i = 0; i < ruleLength; i++) {
        const rule = rules[i];
        const ruleBlock = block(atn, rule, rule);
        if (ruleBlock === undefined) {
            continue;
        }
        buildRuleHandle(atn, rule, ruleBlock);
    }
    return atn;
}
function createRuleStartAndStopATNStates(atn, rules) {
    const ruleLength = rules.length;
    for (let i = 0; i < ruleLength; i++) {
        const rule = rules[i];
        const start = newState(atn, rule, undefined, {
            type: ATN_RULE_START
        });
        const stop = newState(atn, rule, undefined, {
            type: ATN_RULE_STOP
        });
        start.stop = stop;
        atn.ruleToStartState.set(rule, start);
        atn.ruleToStopState.set(rule, stop);
    }
}
function atom(atn, rule, production) {
    if (production instanceof api/* Terminal */.oI) {
        return tokenRef(atn, rule, production.terminalType, production);
    }
    else if (production instanceof api/* NonTerminal */.Sj) {
        return ruleRef(atn, rule, production);
    }
    else if (production instanceof api/* Alternation */.ue) {
        return alternation(atn, rule, production);
    }
    else if (production instanceof api/* Option */.Wx) {
        return atn_option(atn, rule, production);
    }
    else if (production instanceof api/* Repetition */.hI) {
        return repetition(atn, rule, production);
    }
    else if (production instanceof api/* RepetitionWithSeparator */.pT) {
        return repetitionSep(atn, rule, production);
    }
    else if (production instanceof api/* RepetitionMandatory */.ej) {
        return repetitionMandatory(atn, rule, production);
    }
    else if (production instanceof api/* RepetitionMandatoryWithSeparator */.fK) {
        return repetitionMandatorySep(atn, rule, production);
    }
    else {
        return block(atn, rule, production);
    }
}
function repetition(atn, rule, repetition) {
    const starState = newState(atn, rule, repetition, {
        type: ATN_STAR_BLOCK_START
    });
    defineDecisionState(atn, starState);
    const handle = makeAlts(atn, rule, starState, repetition, block(atn, rule, repetition));
    return star(atn, rule, repetition, handle);
}
function repetitionSep(atn, rule, repetition) {
    const starState = newState(atn, rule, repetition, {
        type: ATN_STAR_BLOCK_START
    });
    defineDecisionState(atn, starState);
    const handle = makeAlts(atn, rule, starState, repetition, block(atn, rule, repetition));
    const sep = tokenRef(atn, rule, repetition.separator, repetition);
    return star(atn, rule, repetition, handle, sep);
}
function repetitionMandatory(atn, rule, repetition) {
    const plusState = newState(atn, rule, repetition, {
        type: ATN_PLUS_BLOCK_START
    });
    defineDecisionState(atn, plusState);
    const handle = makeAlts(atn, rule, plusState, repetition, block(atn, rule, repetition));
    return plus(atn, rule, repetition, handle);
}
function repetitionMandatorySep(atn, rule, repetition) {
    const plusState = newState(atn, rule, repetition, {
        type: ATN_PLUS_BLOCK_START
    });
    defineDecisionState(atn, plusState);
    const handle = makeAlts(atn, rule, plusState, repetition, block(atn, rule, repetition));
    const sep = tokenRef(atn, rule, repetition.separator, repetition);
    return plus(atn, rule, repetition, handle, sep);
}
function alternation(atn, rule, alternation) {
    const start = newState(atn, rule, alternation, {
        type: ATN_BASIC
    });
    defineDecisionState(atn, start);
    const alts = (0,map/* default */.Z)(alternation.definition, (e) => atom(atn, rule, e));
    const handle = makeAlts(atn, rule, start, alternation, ...alts);
    return handle;
}
function atn_option(atn, rule, option) {
    const start = newState(atn, rule, option, {
        type: ATN_BASIC
    });
    defineDecisionState(atn, start);
    const handle = makeAlts(atn, rule, start, option, block(atn, rule, option));
    return optional(atn, rule, option, handle);
}
function block(atn, rule, block) {
    const handles = (0,filter/* default */.Z)((0,map/* default */.Z)(block.definition, (e) => atom(atn, rule, e)), (e) => e !== undefined);
    if (handles.length === 1) {
        return handles[0];
    }
    else if (handles.length === 0) {
        return undefined;
    }
    else {
        return makeBlock(atn, handles);
    }
}
function plus(atn, rule, plus, handle, sep) {
    const blkStart = handle.left;
    const blkEnd = handle.right;
    const loop = newState(atn, rule, plus, {
        type: ATN_PLUS_LOOP_BACK
    });
    defineDecisionState(atn, loop);
    const end = newState(atn, rule, plus, {
        type: ATN_LOOP_END
    });
    blkStart.loopback = loop;
    end.loopback = loop;
    atn.decisionMap[buildATNKey(rule, sep ? 'RepetitionMandatoryWithSeparator' : 'RepetitionMandatory', plus.idx)] = loop;
    epsilon(blkEnd, loop); // block can see loop back
    // Depending on whether we have a separator we put the exit transition at index 1 or 0
    // This influences the chosen option in the lookahead DFA
    if (sep === undefined) {
        epsilon(loop, blkStart); // loop back to start
        epsilon(loop, end); // exit
    }
    else {
        epsilon(loop, end); // exit
        // loop back to start with separator
        epsilon(loop, sep.left);
        epsilon(sep.right, blkStart);
    }
    return {
        left: blkStart,
        right: end
    };
}
function star(atn, rule, star, handle, sep) {
    const start = handle.left;
    const end = handle.right;
    const entry = newState(atn, rule, star, {
        type: ATN_STAR_LOOP_ENTRY
    });
    defineDecisionState(atn, entry);
    const loopEnd = newState(atn, rule, star, {
        type: ATN_LOOP_END
    });
    const loop = newState(atn, rule, star, {
        type: ATN_STAR_LOOP_BACK
    });
    entry.loopback = loop;
    loopEnd.loopback = loop;
    epsilon(entry, start); // loop enter edge (alt 2)
    epsilon(entry, loopEnd); // bypass loop edge (alt 1)
    epsilon(end, loop); // block end hits loop back
    if (sep !== undefined) {
        epsilon(loop, loopEnd); // end loop
        // loop back to start of handle using separator
        epsilon(loop, sep.left);
        epsilon(sep.right, start);
    }
    else {
        epsilon(loop, entry); // loop back to entry/exit decision
    }
    atn.decisionMap[buildATNKey(rule, sep ? 'RepetitionWithSeparator' : 'Repetition', star.idx)] = entry;
    return {
        left: entry,
        right: loopEnd
    };
}
function optional(atn, rule, optional, handle) {
    const start = handle.left;
    const end = handle.right;
    epsilon(start, end);
    atn.decisionMap[buildATNKey(rule, 'Option', optional.idx)] = start;
    return handle;
}
function defineDecisionState(atn, state) {
    atn.decisionStates.push(state);
    state.decision = atn.decisionStates.length - 1;
    return state.decision;
}
function makeAlts(atn, rule, start, production, ...alts) {
    const end = newState(atn, rule, production, {
        type: ATN_BLOCK_END,
        start
    });
    start.end = end;
    for (const alt of alts) {
        if (alt !== undefined) {
            // hook alts up to decision block
            epsilon(start, alt.left);
            epsilon(alt.right, end);
        }
        else {
            epsilon(start, end);
        }
    }
    const handle = {
        left: start,
        right: end
    };
    atn.decisionMap[buildATNKey(rule, getProdType(production), production.idx)] = start;
    return handle;
}
function getProdType(production) {
    if (production instanceof api/* Alternation */.ue) {
        return 'Alternation';
    }
    else if (production instanceof api/* Option */.Wx) {
        return 'Option';
    }
    else if (production instanceof api/* Repetition */.hI) {
        return 'Repetition';
    }
    else if (production instanceof api/* RepetitionWithSeparator */.pT) {
        return 'RepetitionWithSeparator';
    }
    else if (production instanceof api/* RepetitionMandatory */.ej) {
        return 'RepetitionMandatory';
    }
    else if (production instanceof api/* RepetitionMandatoryWithSeparator */.fK) {
        return 'RepetitionMandatoryWithSeparator';
    }
    else {
        throw new Error('Invalid production type encountered');
    }
}
function makeBlock(atn, alts) {
    const altsLength = alts.length;
    for (let i = 0; i < altsLength - 1; i++) {
        const handle = alts[i];
        let transition;
        if (handle.left.transitions.length === 1) {
            transition = handle.left.transitions[0];
        }
        const isRuleTransition = transition instanceof RuleTransition;
        const ruleTransition = transition;
        const next = alts[i + 1].left;
        if (handle.left.type === ATN_BASIC &&
            handle.right.type === ATN_BASIC &&
            transition !== undefined &&
            ((isRuleTransition && ruleTransition.followState === handle.right) ||
                transition.target === handle.right)) {
            // we can avoid epsilon edge to next element
            if (isRuleTransition) {
                ruleTransition.followState = next;
            }
            else {
                transition.target = next;
            }
            removeState(atn, handle.right); // we skipped over this state
        }
        else {
            // need epsilon if previous block's right end node is complex
            epsilon(handle.right, next);
        }
    }
    const first = alts[0];
    const last = alts[altsLength - 1];
    return {
        left: first.left,
        right: last.right
    };
}
function tokenRef(atn, rule, tokenType, production) {
    const left = newState(atn, rule, production, {
        type: ATN_BASIC
    });
    const right = newState(atn, rule, production, {
        type: ATN_BASIC
    });
    addTransition(left, new AtomTransition(right, tokenType));
    return {
        left,
        right
    };
}
function ruleRef(atn, currentRule, nonTerminal) {
    const rule = nonTerminal.referencedRule;
    const start = atn.ruleToStartState.get(rule);
    const left = newState(atn, currentRule, nonTerminal, {
        type: ATN_BASIC
    });
    const right = newState(atn, currentRule, nonTerminal, {
        type: ATN_BASIC
    });
    const call = new RuleTransition(start, rule, right);
    addTransition(left, call);
    return {
        left,
        right
    };
}
function buildRuleHandle(atn, rule, block) {
    const start = atn.ruleToStartState.get(rule);
    epsilon(start, block.left);
    const stop = atn.ruleToStopState.get(rule);
    epsilon(block.right, stop);
    const handle = {
        left: start,
        right: stop
    };
    return handle;
}
function epsilon(a, b) {
    const transition = new EpsilonTransition(b);
    addTransition(a, transition);
}
function newState(atn, rule, production, partial) {
    const t = Object.assign({ atn,
        production, epsilonOnlyTransitions: false, rule, transitions: [], nextTokenWithinRule: [], stateNumber: atn.states.length }, partial);
    atn.states.push(t);
    return t;
}
function addTransition(state, transition) {
    // A single ATN state can only contain epsilon transitions or non-epsilon transitions
    // Because they are never mixed, only setting the property for the first transition is fine
    if (state.transitions.length === 0) {
        state.epsilonOnlyTransitions = transition.isEpsilon();
    }
    state.transitions.push(transition);
}
function removeState(atn, state) {
    atn.states.splice(atn.states.indexOf(state), 1);
}
//# sourceMappingURL=atn.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain-allstar/lib/dfa.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

const DFA_ERROR = {};
class ATNConfigSet {
    constructor() {
        this.map = {};
        this.configs = [];
    }
    get size() {
        return this.configs.length;
    }
    finalize() {
        // Empties the map to free up memory
        this.map = {};
    }
    add(config) {
        const key = getATNConfigKey(config);
        // Only add configs which don't exist in our map already
        // While this does not influence the actual algorithm, adding them anyway would massively increase memory consumption
        if (!(key in this.map)) {
            this.map[key] = this.configs.length;
            this.configs.push(config);
        }
    }
    get elements() {
        return this.configs;
    }
    get alts() {
        return (0,map/* default */.Z)(this.configs, (e) => e.alt);
    }
    get key() {
        let value = "";
        for (const k in this.map) {
            value += k + ":";
        }
        return value;
    }
}
function getATNConfigKey(config, alt = true) {
    return `${alt ? `a${config.alt}` : ""}s${config.state.stateNumber}:${config.stack.map((e) => e.stateNumber.toString()).join("_")}`;
}
//# sourceMappingURL=dfa.js.map
// EXTERNAL MODULE: ../node_modules/lodash-es/min.js
var min = __webpack_require__(18519);
// EXTERNAL MODULE: ../node_modules/lodash-es/flatMap.js
var flatMap = __webpack_require__(34134);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseIteratee.js + 15 modules
var _baseIteratee = __webpack_require__(86494);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseUniq.js + 1 modules
var _baseUniq = __webpack_require__(99633);
;// CONCATENATED MODULE: ../node_modules/lodash-es/uniqBy.js



/**
 * This method is like `_.uniq` except that it accepts `iteratee` which is
 * invoked for each element in `array` to generate the criterion by which
 * uniqueness is computed. The order of result values is determined by the
 * order they occur in the array. The iteratee is invoked with one argument:
 * (value).
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Array
 * @param {Array} array The array to inspect.
 * @param {Function} [iteratee=_.identity] The iteratee invoked per element.
 * @returns {Array} Returns the new duplicate free array.
 * @example
 *
 * _.uniqBy([2.1, 1.2, 2.3], Math.floor);
 * // => [2.1, 1.2]
 *
 * // The `_.property` iteratee shorthand.
 * _.uniqBy([{ 'x': 1 }, { 'x': 2 }, { 'x': 1 }], 'x');
 * // => [{ 'x': 1 }, { 'x': 2 }]
 */
function uniqBy(array, iteratee) {
  return (array && array.length) ? (0,_baseUniq/* default */.Z)(array, (0,_baseIteratee/* default */.Z)(iteratee, 2)) : [];
}

/* harmony default export */ const lodash_es_uniqBy = (uniqBy);

// EXTERNAL MODULE: ../node_modules/lodash-es/flatten.js
var flatten = __webpack_require__(28099);
// EXTERNAL MODULE: ../node_modules/lodash-es/forEach.js
var forEach = __webpack_require__(21845);
// EXTERNAL MODULE: ../node_modules/lodash-es/isEmpty.js
var isEmpty = __webpack_require__(66400);
// EXTERNAL MODULE: ../node_modules/lodash-es/reduce.js + 2 modules
var reduce = __webpack_require__(99413);
;// CONCATENATED MODULE: ../node_modules/chevrotain-allstar/lib/all-star-lookahead.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/











function createDFACache(startState, decision) {
    const map = {};
    return (predicateSet) => {
        const key = predicateSet.toString();
        let existing = map[key];
        if (existing !== undefined) {
            return existing;
        }
        else {
            existing = {
                atnStartState: startState,
                decision,
                states: {}
            };
            map[key] = existing;
            return existing;
        }
    };
}
class PredicateSet {
    constructor() {
        this.predicates = [];
    }
    is(index) {
        return index >= this.predicates.length || this.predicates[index];
    }
    set(index, value) {
        this.predicates[index] = value;
    }
    toString() {
        let value = "";
        const size = this.predicates.length;
        for (let i = 0; i < size; i++) {
            value += this.predicates[i] === true ? "1" : "0";
        }
        return value;
    }
}
const EMPTY_PREDICATES = new PredicateSet();
class LLStarLookaheadStrategy extends api/* LLkLookaheadStrategy */.dV {
    constructor(options) {
        var _a;
        super();
        this.logging = (_a = options === null || options === void 0 ? void 0 : options.logging) !== null && _a !== void 0 ? _a : ((message) => console.log(message));
    }
    initialize(options) {
        this.atn = createATN(options.rules);
        this.dfas = initATNSimulator(this.atn);
    }
    validateAmbiguousAlternationAlternatives() {
        return [];
    }
    validateEmptyOrAlternatives() {
        return [];
    }
    buildLookaheadForAlternation(options) {
        const { prodOccurrence, rule, hasPredicates, dynamicTokensEnabled } = options;
        const dfas = this.dfas;
        const logging = this.logging;
        const key = buildATNKey(rule, 'Alternation', prodOccurrence);
        const decisionState = this.atn.decisionMap[key];
        const decisionIndex = decisionState.decision;
        const partialAlts = (0,map/* default */.Z)((0,api/* getLookaheadPaths */.oC)({
            maxLookahead: 1,
            occurrence: prodOccurrence,
            prodType: "Alternation",
            rule: rule
        }), (currAlt) => (0,map/* default */.Z)(currAlt, (path) => path[0]));
        if (isLL1Sequence(partialAlts, false) && !dynamicTokensEnabled) {
            const choiceToAlt = (0,reduce/* default */.Z)(partialAlts, (result, currAlt, idx) => {
                (0,forEach/* default */.Z)(currAlt, (currTokType) => {
                    if (currTokType) {
                        result[currTokType.tokenTypeIdx] = idx;
                        (0,forEach/* default */.Z)(currTokType.categoryMatches, (currExtendingType) => {
                            result[currExtendingType] = idx;
                        });
                    }
                });
                return result;
            }, {});
            if (hasPredicates) {
                return function (orAlts) {
                    var _a;
                    const nextToken = this.LA(1);
                    const prediction = choiceToAlt[nextToken.tokenTypeIdx];
                    if (orAlts !== undefined && prediction !== undefined) {
                        const gate = (_a = orAlts[prediction]) === null || _a === void 0 ? void 0 : _a.GATE;
                        if (gate !== undefined && gate.call(this) === false) {
                            return undefined;
                        }
                    }
                    return prediction;
                };
            }
            else {
                return function () {
                    const nextToken = this.LA(1);
                    return choiceToAlt[nextToken.tokenTypeIdx];
                };
            }
        }
        else if (hasPredicates) {
            return function (orAlts) {
                const predicates = new PredicateSet();
                const length = orAlts === undefined ? 0 : orAlts.length;
                for (let i = 0; i < length; i++) {
                    const gate = orAlts === null || orAlts === void 0 ? void 0 : orAlts[i].GATE;
                    predicates.set(i, gate === undefined || gate.call(this));
                }
                const result = adaptivePredict.call(this, dfas, decisionIndex, predicates, logging);
                return typeof result === 'number' ? result : undefined;
            };
        }
        else {
            return function () {
                const result = adaptivePredict.call(this, dfas, decisionIndex, EMPTY_PREDICATES, logging);
                return typeof result === 'number' ? result : undefined;
            };
        }
    }
    buildLookaheadForOptional(options) {
        const { prodOccurrence, rule, prodType, dynamicTokensEnabled } = options;
        const dfas = this.dfas;
        const logging = this.logging;
        const key = buildATNKey(rule, prodType, prodOccurrence);
        const decisionState = this.atn.decisionMap[key];
        const decisionIndex = decisionState.decision;
        const alts = (0,map/* default */.Z)((0,api/* getLookaheadPaths */.oC)({
            maxLookahead: 1,
            occurrence: prodOccurrence,
            prodType,
            rule
        }), (e) => {
            return (0,map/* default */.Z)(e, (g) => g[0]);
        });
        if (isLL1Sequence(alts) && alts[0][0] && !dynamicTokensEnabled) {
            const alt = alts[0];
            const singleTokensTypes = (0,flatten/* default */.Z)(alt);
            if (singleTokensTypes.length === 1 &&
                (0,isEmpty/* default */.Z)(singleTokensTypes[0].categoryMatches)) {
                const expectedTokenType = singleTokensTypes[0];
                const expectedTokenUniqueKey = expectedTokenType.tokenTypeIdx;
                return function () {
                    return this.LA(1).tokenTypeIdx === expectedTokenUniqueKey;
                };
            }
            else {
                const choiceToAlt = (0,reduce/* default */.Z)(singleTokensTypes, (result, currTokType) => {
                    if (currTokType !== undefined) {
                        result[currTokType.tokenTypeIdx] = true;
                        (0,forEach/* default */.Z)(currTokType.categoryMatches, (currExtendingType) => {
                            result[currExtendingType] = true;
                        });
                    }
                    return result;
                }, {});
                return function () {
                    const nextToken = this.LA(1);
                    return choiceToAlt[nextToken.tokenTypeIdx] === true;
                };
            }
        }
        return function () {
            const result = adaptivePredict.call(this, dfas, decisionIndex, EMPTY_PREDICATES, logging);
            return typeof result === "object" ? false : result === 0;
        };
    }
}
function isLL1Sequence(sequences, allowEmpty = true) {
    const fullSet = new Set();
    for (const alt of sequences) {
        const altSet = new Set();
        for (const tokType of alt) {
            if (tokType === undefined) {
                if (allowEmpty) {
                    // Epsilon production encountered
                    break;
                }
                else {
                    return false;
                }
            }
            const indices = [tokType.tokenTypeIdx].concat(tokType.categoryMatches);
            for (const index of indices) {
                if (fullSet.has(index)) {
                    if (!altSet.has(index)) {
                        return false;
                    }
                }
                else {
                    fullSet.add(index);
                    altSet.add(index);
                }
            }
        }
    }
    return true;
}
function initATNSimulator(atn) {
    const decisionLength = atn.decisionStates.length;
    const decisionToDFA = Array(decisionLength);
    for (let i = 0; i < decisionLength; i++) {
        decisionToDFA[i] = createDFACache(atn.decisionStates[i], i);
    }
    return decisionToDFA;
}
function adaptivePredict(dfaCaches, decision, predicateSet, logging) {
    const dfa = dfaCaches[decision](predicateSet);
    let start = dfa.start;
    if (start === undefined) {
        const closure = computeStartState(dfa.atnStartState);
        start = addDFAState(dfa, newDFAState(closure));
        dfa.start = start;
    }
    const alt = performLookahead.apply(this, [dfa, start, predicateSet, logging]);
    return alt;
}
function performLookahead(dfa, s0, predicateSet, logging) {
    let previousD = s0;
    let i = 1;
    const path = [];
    let t = this.LA(i++);
    while (true) {
        let d = getExistingTargetState(previousD, t);
        if (d === undefined) {
            d = computeLookaheadTarget.apply(this, [dfa, previousD, t, i, predicateSet, logging]);
        }
        if (d === DFA_ERROR) {
            return buildAdaptivePredictError(path, previousD, t);
        }
        if (d.isAcceptState === true) {
            return d.prediction;
        }
        previousD = d;
        path.push(t);
        t = this.LA(i++);
    }
}
function computeLookaheadTarget(dfa, previousD, token, lookahead, predicateSet, logging) {
    const reach = computeReachSet(previousD.configs, token, predicateSet);
    if (reach.size === 0) {
        addDFAEdge(dfa, previousD, token, DFA_ERROR);
        return DFA_ERROR;
    }
    let newState = newDFAState(reach);
    const predictedAlt = getUniqueAlt(reach, predicateSet);
    if (predictedAlt !== undefined) {
        newState.isAcceptState = true;
        newState.prediction = predictedAlt;
        newState.configs.uniqueAlt = predictedAlt;
    }
    else if (hasConflictTerminatingPrediction(reach)) {
        const prediction = (0,min/* default */.Z)(reach.alts);
        newState.isAcceptState = true;
        newState.prediction = prediction;
        newState.configs.uniqueAlt = prediction;
        reportLookaheadAmbiguity.apply(this, [dfa, lookahead, reach.alts, logging]);
    }
    newState = addDFAEdge(dfa, previousD, token, newState);
    return newState;
}
function reportLookaheadAmbiguity(dfa, lookahead, ambiguityIndices, logging) {
    const prefixPath = [];
    for (let i = 1; i <= lookahead; i++) {
        prefixPath.push(this.LA(i).tokenType);
    }
    const atnState = dfa.atnStartState;
    const topLevelRule = atnState.rule;
    const production = atnState.production;
    const message = buildAmbiguityError({
        topLevelRule,
        ambiguityIndices,
        production,
        prefixPath
    });
    logging(message);
}
function buildAmbiguityError(options) {
    const pathMsg = (0,map/* default */.Z)(options.prefixPath, (currtok) => (0,api/* tokenLabel */.l$)(currtok)).join(", ");
    const occurrence = options.production.idx === 0 ? "" : options.production.idx;
    let currMessage = `Ambiguous Alternatives Detected: <${options.ambiguityIndices.join(", ")}> in <${getProductionDslName(options.production)}${occurrence}>` +
        ` inside <${options.topLevelRule.name}> Rule,\n` +
        `<${pathMsg}> may appears as a prefix path in all these alternatives.\n`;
    currMessage =
        currMessage +
            `See: https://chevrotain.io/docs/guide/resolving_grammar_errors.html#AMBIGUOUS_ALTERNATIVES\n` +
            `For Further details.`;
    return currMessage;
}
function getProductionDslName(prod) {
    if (prod instanceof api/* NonTerminal */.Sj) {
        return "SUBRULE";
    }
    else if (prod instanceof api/* Option */.Wx) {
        return "OPTION";
    }
    else if (prod instanceof api/* Alternation */.ue) {
        return "OR";
    }
    else if (prod instanceof api/* RepetitionMandatory */.ej) {
        return "AT_LEAST_ONE";
    }
    else if (prod instanceof api/* RepetitionMandatoryWithSeparator */.fK) {
        return "AT_LEAST_ONE_SEP";
    }
    else if (prod instanceof api/* RepetitionWithSeparator */.pT) {
        return "MANY_SEP";
    }
    else if (prod instanceof api/* Repetition */.hI) {
        return "MANY";
    }
    else if (prod instanceof api/* Terminal */.oI) {
        return "CONSUME";
    }
    else {
        throw Error("non exhaustive match");
    }
}
function buildAdaptivePredictError(path, previous, current) {
    const nextTransitions = (0,flatMap/* default */.Z)(previous.configs.elements, (e) => e.state.transitions);
    const nextTokenTypes = lodash_es_uniqBy(nextTransitions
        .filter((e) => e instanceof AtomTransition)
        .map((e) => e.tokenType), (e) => e.tokenTypeIdx);
    return {
        actualToken: current,
        possibleTokenTypes: nextTokenTypes,
        tokenPath: path
    };
}
function getExistingTargetState(state, token) {
    return state.edges[token.tokenTypeIdx];
}
function computeReachSet(configs, token, predicateSet) {
    const intermediate = new ATNConfigSet();
    const skippedStopStates = [];
    for (const c of configs.elements) {
        if (predicateSet.is(c.alt) === false) {
            continue;
        }
        if (c.state.type === ATN_RULE_STOP) {
            skippedStopStates.push(c);
            continue;
        }
        const transitionLength = c.state.transitions.length;
        for (let i = 0; i < transitionLength; i++) {
            const transition = c.state.transitions[i];
            const target = getReachableTarget(transition, token);
            if (target !== undefined) {
                intermediate.add({
                    state: target,
                    alt: c.alt,
                    stack: c.stack
                });
            }
        }
    }
    let reach;
    if (skippedStopStates.length === 0 && intermediate.size === 1) {
        reach = intermediate;
    }
    if (reach === undefined) {
        reach = new ATNConfigSet();
        for (const c of intermediate.elements) {
            closure(c, reach);
        }
    }
    if (skippedStopStates.length > 0 && !hasConfigInRuleStopState(reach)) {
        for (const c of skippedStopStates) {
            reach.add(c);
        }
    }
    return reach;
}
function getReachableTarget(transition, token) {
    if (transition instanceof AtomTransition &&
        (0,api/* tokenMatcher */.ol)(token, transition.tokenType)) {
        return transition.target;
    }
    return undefined;
}
function getUniqueAlt(configs, predicateSet) {
    let alt;
    for (const c of configs.elements) {
        if (predicateSet.is(c.alt) === true) {
            if (alt === undefined) {
                alt = c.alt;
            }
            else if (alt !== c.alt) {
                return undefined;
            }
        }
    }
    return alt;
}
function newDFAState(closure) {
    return {
        configs: closure,
        edges: {},
        isAcceptState: false,
        prediction: -1
    };
}
function addDFAEdge(dfa, from, token, to) {
    to = addDFAState(dfa, to);
    from.edges[token.tokenTypeIdx] = to;
    return to;
}
function addDFAState(dfa, state) {
    if (state === DFA_ERROR) {
        return state;
    }
    // Repetitions have the same config set
    // Therefore, storing the key of the config in a map allows us to create a loop in our DFA
    const mapKey = state.configs.key;
    const existing = dfa.states[mapKey];
    if (existing !== undefined) {
        return existing;
    }
    state.configs.finalize();
    dfa.states[mapKey] = state;
    return state;
}
function computeStartState(atnState) {
    const configs = new ATNConfigSet();
    const numberOfTransitions = atnState.transitions.length;
    for (let i = 0; i < numberOfTransitions; i++) {
        const target = atnState.transitions[i].target;
        const config = {
            state: target,
            alt: i,
            stack: []
        };
        closure(config, configs);
    }
    return configs;
}
function closure(config, configs) {
    const p = config.state;
    if (p.type === ATN_RULE_STOP) {
        if (config.stack.length > 0) {
            const atnStack = [...config.stack];
            const followState = atnStack.pop();
            const followConfig = {
                state: followState,
                alt: config.alt,
                stack: atnStack
            };
            closure(followConfig, configs);
        }
        else {
            // Dipping into outer context, simply add the config
            // This will stop computation once every config is at the rule stop state
            configs.add(config);
        }
        return;
    }
    if (!p.epsilonOnlyTransitions) {
        configs.add(config);
    }
    const transitionLength = p.transitions.length;
    for (let i = 0; i < transitionLength; i++) {
        const transition = p.transitions[i];
        const c = getEpsilonTarget(config, transition);
        if (c !== undefined) {
            closure(c, configs);
        }
    }
}
function getEpsilonTarget(config, transition) {
    if (transition instanceof EpsilonTransition) {
        return {
            state: transition.target,
            alt: config.alt,
            stack: config.stack
        };
    }
    else if (transition instanceof RuleTransition) {
        const stack = [...config.stack, transition.followState];
        return {
            state: transition.target,
            alt: config.alt,
            stack
        };
    }
    return undefined;
}
function hasConfigInRuleStopState(configs) {
    for (const c of configs.elements) {
        if (c.state.type === ATN_RULE_STOP) {
            return true;
        }
    }
    return false;
}
function allConfigsInRuleStopStates(configs) {
    for (const c of configs.elements) {
        if (c.state.type !== ATN_RULE_STOP) {
            return false;
        }
    }
    return true;
}
function hasConflictTerminatingPrediction(configs) {
    if (allConfigsInRuleStopStates(configs)) {
        return true;
    }
    const altSets = getConflictingAltSets(configs.elements);
    const heuristic = hasConflictingAltSet(altSets) && !hasStateAssociatedWithOneAlt(altSets);
    return heuristic;
}
function getConflictingAltSets(configs) {
    const configToAlts = new Map();
    for (const c of configs) {
        const key = getATNConfigKey(c, false);
        let alts = configToAlts.get(key);
        if (alts === undefined) {
            alts = {};
            configToAlts.set(key, alts);
        }
        alts[c.alt] = true;
    }
    return configToAlts;
}
function hasConflictingAltSet(altSets) {
    for (const value of Array.from(altSets.values())) {
        if (Object.keys(value).length > 1) {
            return true;
        }
    }
    return false;
}
function hasStateAssociatedWithOneAlt(altSets) {
    for (const value of Array.from(altSets.values())) {
        if (Object.keys(value).length === 1) {
            return true;
        }
    }
    return false;
}
//# sourceMappingURL=all-star-lookahead.js.map
;// CONCATENATED MODULE: ../node_modules/chevrotain-allstar/lib/index.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

//# sourceMappingURL=index.js.map
// EXTERNAL MODULE: ../node_modules/langium/lib/utils/ast-utils.js
var ast_utils = __webpack_require__(74857);
;// CONCATENATED MODULE: ../node_modules/vscode-languageserver-types/lib/esm/main.js
/* --------------------------------------------------------------------------------------------
 * Copyright (c) Microsoft Corporation. All rights reserved.
 * Licensed under the MIT License. See License.txt in the project root for license information.
 * ------------------------------------------------------------------------------------------ */

var DocumentUri;
(function (DocumentUri) {
    function is(value) {
        return typeof value === 'string';
    }
    DocumentUri.is = is;
})(DocumentUri || (DocumentUri = {}));
var URI;
(function (URI) {
    function is(value) {
        return typeof value === 'string';
    }
    URI.is = is;
})(URI || (URI = {}));
var integer;
(function (integer) {
    integer.MIN_VALUE = -2147483648;
    integer.MAX_VALUE = 2147483647;
    function is(value) {
        return typeof value === 'number' && integer.MIN_VALUE <= value && value <= integer.MAX_VALUE;
    }
    integer.is = is;
})(integer || (integer = {}));
var uinteger;
(function (uinteger) {
    uinteger.MIN_VALUE = 0;
    uinteger.MAX_VALUE = 2147483647;
    function is(value) {
        return typeof value === 'number' && uinteger.MIN_VALUE <= value && value <= uinteger.MAX_VALUE;
    }
    uinteger.is = is;
})(uinteger || (uinteger = {}));
/**
 * The Position namespace provides helper functions to work with
 * {@link Position} literals.
 */
var Position;
(function (Position) {
    /**
     * Creates a new Position literal from the given line and character.
     * @param line The position's line.
     * @param character The position's character.
     */
    function create(line, character) {
        if (line === Number.MAX_VALUE) {
            line = uinteger.MAX_VALUE;
        }
        if (character === Number.MAX_VALUE) {
            character = uinteger.MAX_VALUE;
        }
        return { line, character };
    }
    Position.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Position} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.objectLiteral(candidate) && Is.uinteger(candidate.line) && Is.uinteger(candidate.character);
    }
    Position.is = is;
})(Position || (Position = {}));
/**
 * The Range namespace provides helper functions to work with
 * {@link Range} literals.
 */
var Range;
(function (Range) {
    function create(one, two, three, four) {
        if (Is.uinteger(one) && Is.uinteger(two) && Is.uinteger(three) && Is.uinteger(four)) {
            return { start: Position.create(one, two), end: Position.create(three, four) };
        }
        else if (Position.is(one) && Position.is(two)) {
            return { start: one, end: two };
        }
        else {
            throw new Error(`Range#create called with invalid arguments[${one}, ${two}, ${three}, ${four}]`);
        }
    }
    Range.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Range} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.objectLiteral(candidate) && Position.is(candidate.start) && Position.is(candidate.end);
    }
    Range.is = is;
})(Range || (Range = {}));
/**
 * The Location namespace provides helper functions to work with
 * {@link Location} literals.
 */
var Location;
(function (Location) {
    /**
     * Creates a Location literal.
     * @param uri The location's uri.
     * @param range The location's range.
     */
    function create(uri, range) {
        return { uri, range };
    }
    Location.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Location} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.objectLiteral(candidate) && Range.is(candidate.range) && (Is.string(candidate.uri) || Is.undefined(candidate.uri));
    }
    Location.is = is;
})(Location || (Location = {}));
/**
 * The LocationLink namespace provides helper functions to work with
 * {@link LocationLink} literals.
 */
var LocationLink;
(function (LocationLink) {
    /**
     * Creates a LocationLink literal.
     * @param targetUri The definition's uri.
     * @param targetRange The full range of the definition.
     * @param targetSelectionRange The span of the symbol definition at the target.
     * @param originSelectionRange The span of the symbol being defined in the originating source file.
     */
    function create(targetUri, targetRange, targetSelectionRange, originSelectionRange) {
        return { targetUri, targetRange, targetSelectionRange, originSelectionRange };
    }
    LocationLink.create = create;
    /**
     * Checks whether the given literal conforms to the {@link LocationLink} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.objectLiteral(candidate) && Range.is(candidate.targetRange) && Is.string(candidate.targetUri)
            && Range.is(candidate.targetSelectionRange)
            && (Range.is(candidate.originSelectionRange) || Is.undefined(candidate.originSelectionRange));
    }
    LocationLink.is = is;
})(LocationLink || (LocationLink = {}));
/**
 * The Color namespace provides helper functions to work with
 * {@link Color} literals.
 */
var Color;
(function (Color) {
    /**
     * Creates a new Color literal.
     */
    function create(red, green, blue, alpha) {
        return {
            red,
            green,
            blue,
            alpha,
        };
    }
    Color.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Color} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Is.numberRange(candidate.red, 0, 1)
            && Is.numberRange(candidate.green, 0, 1)
            && Is.numberRange(candidate.blue, 0, 1)
            && Is.numberRange(candidate.alpha, 0, 1);
    }
    Color.is = is;
})(Color || (Color = {}));
/**
 * The ColorInformation namespace provides helper functions to work with
 * {@link ColorInformation} literals.
 */
var ColorInformation;
(function (ColorInformation) {
    /**
     * Creates a new ColorInformation literal.
     */
    function create(range, color) {
        return {
            range,
            color,
        };
    }
    ColorInformation.create = create;
    /**
     * Checks whether the given literal conforms to the {@link ColorInformation} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Range.is(candidate.range) && Color.is(candidate.color);
    }
    ColorInformation.is = is;
})(ColorInformation || (ColorInformation = {}));
/**
 * The Color namespace provides helper functions to work with
 * {@link ColorPresentation} literals.
 */
var ColorPresentation;
(function (ColorPresentation) {
    /**
     * Creates a new ColorInformation literal.
     */
    function create(label, textEdit, additionalTextEdits) {
        return {
            label,
            textEdit,
            additionalTextEdits,
        };
    }
    ColorPresentation.create = create;
    /**
     * Checks whether the given literal conforms to the {@link ColorInformation} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Is.string(candidate.label)
            && (Is.undefined(candidate.textEdit) || TextEdit.is(candidate))
            && (Is.undefined(candidate.additionalTextEdits) || Is.typedArray(candidate.additionalTextEdits, TextEdit.is));
    }
    ColorPresentation.is = is;
})(ColorPresentation || (ColorPresentation = {}));
/**
 * A set of predefined range kinds.
 */
var FoldingRangeKind;
(function (FoldingRangeKind) {
    /**
     * Folding range for a comment
     */
    FoldingRangeKind.Comment = 'comment';
    /**
     * Folding range for an import or include
     */
    FoldingRangeKind.Imports = 'imports';
    /**
     * Folding range for a region (e.g. `#region`)
     */
    FoldingRangeKind.Region = 'region';
})(FoldingRangeKind || (FoldingRangeKind = {}));
/**
 * The folding range namespace provides helper functions to work with
 * {@link FoldingRange} literals.
 */
var FoldingRange;
(function (FoldingRange) {
    /**
     * Creates a new FoldingRange literal.
     */
    function create(startLine, endLine, startCharacter, endCharacter, kind, collapsedText) {
        const result = {
            startLine,
            endLine
        };
        if (Is.defined(startCharacter)) {
            result.startCharacter = startCharacter;
        }
        if (Is.defined(endCharacter)) {
            result.endCharacter = endCharacter;
        }
        if (Is.defined(kind)) {
            result.kind = kind;
        }
        if (Is.defined(collapsedText)) {
            result.collapsedText = collapsedText;
        }
        return result;
    }
    FoldingRange.create = create;
    /**
     * Checks whether the given literal conforms to the {@link FoldingRange} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Is.uinteger(candidate.startLine) && Is.uinteger(candidate.startLine)
            && (Is.undefined(candidate.startCharacter) || Is.uinteger(candidate.startCharacter))
            && (Is.undefined(candidate.endCharacter) || Is.uinteger(candidate.endCharacter))
            && (Is.undefined(candidate.kind) || Is.string(candidate.kind));
    }
    FoldingRange.is = is;
})(FoldingRange || (FoldingRange = {}));
/**
 * The DiagnosticRelatedInformation namespace provides helper functions to work with
 * {@link DiagnosticRelatedInformation} literals.
 */
var DiagnosticRelatedInformation;
(function (DiagnosticRelatedInformation) {
    /**
     * Creates a new DiagnosticRelatedInformation literal.
     */
    function create(location, message) {
        return {
            location,
            message
        };
    }
    DiagnosticRelatedInformation.create = create;
    /**
     * Checks whether the given literal conforms to the {@link DiagnosticRelatedInformation} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Location.is(candidate.location) && Is.string(candidate.message);
    }
    DiagnosticRelatedInformation.is = is;
})(DiagnosticRelatedInformation || (DiagnosticRelatedInformation = {}));
/**
 * The diagnostic's severity.
 */
var DiagnosticSeverity;
(function (DiagnosticSeverity) {
    /**
     * Reports an error.
     */
    DiagnosticSeverity.Error = 1;
    /**
     * Reports a warning.
     */
    DiagnosticSeverity.Warning = 2;
    /**
     * Reports an information.
     */
    DiagnosticSeverity.Information = 3;
    /**
     * Reports a hint.
     */
    DiagnosticSeverity.Hint = 4;
})(DiagnosticSeverity || (DiagnosticSeverity = {}));
/**
 * The diagnostic tags.
 *
 * @since 3.15.0
 */
var DiagnosticTag;
(function (DiagnosticTag) {
    /**
     * Unused or unnecessary code.
     *
     * Clients are allowed to render diagnostics with this tag faded out instead of having
     * an error squiggle.
     */
    DiagnosticTag.Unnecessary = 1;
    /**
     * Deprecated or obsolete code.
     *
     * Clients are allowed to rendered diagnostics with this tag strike through.
     */
    DiagnosticTag.Deprecated = 2;
})(DiagnosticTag || (DiagnosticTag = {}));
/**
 * The CodeDescription namespace provides functions to deal with descriptions for diagnostic codes.
 *
 * @since 3.16.0
 */
var CodeDescription;
(function (CodeDescription) {
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Is.string(candidate.href);
    }
    CodeDescription.is = is;
})(CodeDescription || (CodeDescription = {}));
/**
 * The Diagnostic namespace provides helper functions to work with
 * {@link Diagnostic} literals.
 */
var Diagnostic;
(function (Diagnostic) {
    /**
     * Creates a new Diagnostic literal.
     */
    function create(range, message, severity, code, source, relatedInformation) {
        let result = { range, message };
        if (Is.defined(severity)) {
            result.severity = severity;
        }
        if (Is.defined(code)) {
            result.code = code;
        }
        if (Is.defined(source)) {
            result.source = source;
        }
        if (Is.defined(relatedInformation)) {
            result.relatedInformation = relatedInformation;
        }
        return result;
    }
    Diagnostic.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Diagnostic} interface.
     */
    function is(value) {
        var _a;
        let candidate = value;
        return Is.defined(candidate)
            && Range.is(candidate.range)
            && Is.string(candidate.message)
            && (Is.number(candidate.severity) || Is.undefined(candidate.severity))
            && (Is.integer(candidate.code) || Is.string(candidate.code) || Is.undefined(candidate.code))
            && (Is.undefined(candidate.codeDescription) || (Is.string((_a = candidate.codeDescription) === null || _a === void 0 ? void 0 : _a.href)))
            && (Is.string(candidate.source) || Is.undefined(candidate.source))
            && (Is.undefined(candidate.relatedInformation) || Is.typedArray(candidate.relatedInformation, DiagnosticRelatedInformation.is));
    }
    Diagnostic.is = is;
})(Diagnostic || (Diagnostic = {}));
/**
 * The Command namespace provides helper functions to work with
 * {@link Command} literals.
 */
var Command;
(function (Command) {
    /**
     * Creates a new Command literal.
     */
    function create(title, command, ...args) {
        let result = { title, command };
        if (Is.defined(args) && args.length > 0) {
            result.arguments = args;
        }
        return result;
    }
    Command.create = create;
    /**
     * Checks whether the given literal conforms to the {@link Command} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.title) && Is.string(candidate.command);
    }
    Command.is = is;
})(Command || (Command = {}));
/**
 * The TextEdit namespace provides helper function to create replace,
 * insert and delete edits more easily.
 */
var TextEdit;
(function (TextEdit) {
    /**
     * Creates a replace text edit.
     * @param range The range of text to be replaced.
     * @param newText The new text.
     */
    function replace(range, newText) {
        return { range, newText };
    }
    TextEdit.replace = replace;
    /**
     * Creates an insert text edit.
     * @param position The position to insert the text at.
     * @param newText The text to be inserted.
     */
    function insert(position, newText) {
        return { range: { start: position, end: position }, newText };
    }
    TextEdit.insert = insert;
    /**
     * Creates a delete text edit.
     * @param range The range of text to be deleted.
     */
    function del(range) {
        return { range, newText: '' };
    }
    TextEdit.del = del;
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate)
            && Is.string(candidate.newText)
            && Range.is(candidate.range);
    }
    TextEdit.is = is;
})(TextEdit || (TextEdit = {}));
var ChangeAnnotation;
(function (ChangeAnnotation) {
    function create(label, needsConfirmation, description) {
        const result = { label };
        if (needsConfirmation !== undefined) {
            result.needsConfirmation = needsConfirmation;
        }
        if (description !== undefined) {
            result.description = description;
        }
        return result;
    }
    ChangeAnnotation.create = create;
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Is.string(candidate.label) &&
            (Is.boolean(candidate.needsConfirmation) || candidate.needsConfirmation === undefined) &&
            (Is.string(candidate.description) || candidate.description === undefined);
    }
    ChangeAnnotation.is = is;
})(ChangeAnnotation || (ChangeAnnotation = {}));
var ChangeAnnotationIdentifier;
(function (ChangeAnnotationIdentifier) {
    function is(value) {
        const candidate = value;
        return Is.string(candidate);
    }
    ChangeAnnotationIdentifier.is = is;
})(ChangeAnnotationIdentifier || (ChangeAnnotationIdentifier = {}));
var AnnotatedTextEdit;
(function (AnnotatedTextEdit) {
    /**
     * Creates an annotated replace text edit.
     *
     * @param range The range of text to be replaced.
     * @param newText The new text.
     * @param annotation The annotation.
     */
    function replace(range, newText, annotation) {
        return { range, newText, annotationId: annotation };
    }
    AnnotatedTextEdit.replace = replace;
    /**
     * Creates an annotated insert text edit.
     *
     * @param position The position to insert the text at.
     * @param newText The text to be inserted.
     * @param annotation The annotation.
     */
    function insert(position, newText, annotation) {
        return { range: { start: position, end: position }, newText, annotationId: annotation };
    }
    AnnotatedTextEdit.insert = insert;
    /**
     * Creates an annotated delete text edit.
     *
     * @param range The range of text to be deleted.
     * @param annotation The annotation.
     */
    function del(range, annotation) {
        return { range, newText: '', annotationId: annotation };
    }
    AnnotatedTextEdit.del = del;
    function is(value) {
        const candidate = value;
        return TextEdit.is(candidate) && (ChangeAnnotation.is(candidate.annotationId) || ChangeAnnotationIdentifier.is(candidate.annotationId));
    }
    AnnotatedTextEdit.is = is;
})(AnnotatedTextEdit || (AnnotatedTextEdit = {}));
/**
 * The TextDocumentEdit namespace provides helper function to create
 * an edit that manipulates a text document.
 */
var TextDocumentEdit;
(function (TextDocumentEdit) {
    /**
     * Creates a new `TextDocumentEdit`
     */
    function create(textDocument, edits) {
        return { textDocument, edits };
    }
    TextDocumentEdit.create = create;
    function is(value) {
        let candidate = value;
        return Is.defined(candidate)
            && OptionalVersionedTextDocumentIdentifier.is(candidate.textDocument)
            && Array.isArray(candidate.edits);
    }
    TextDocumentEdit.is = is;
})(TextDocumentEdit || (TextDocumentEdit = {}));
var CreateFile;
(function (CreateFile) {
    function create(uri, options, annotation) {
        let result = {
            kind: 'create',
            uri
        };
        if (options !== undefined && (options.overwrite !== undefined || options.ignoreIfExists !== undefined)) {
            result.options = options;
        }
        if (annotation !== undefined) {
            result.annotationId = annotation;
        }
        return result;
    }
    CreateFile.create = create;
    function is(value) {
        let candidate = value;
        return candidate && candidate.kind === 'create' && Is.string(candidate.uri) && (candidate.options === undefined ||
            ((candidate.options.overwrite === undefined || Is.boolean(candidate.options.overwrite)) && (candidate.options.ignoreIfExists === undefined || Is.boolean(candidate.options.ignoreIfExists)))) && (candidate.annotationId === undefined || ChangeAnnotationIdentifier.is(candidate.annotationId));
    }
    CreateFile.is = is;
})(CreateFile || (CreateFile = {}));
var RenameFile;
(function (RenameFile) {
    function create(oldUri, newUri, options, annotation) {
        let result = {
            kind: 'rename',
            oldUri,
            newUri
        };
        if (options !== undefined && (options.overwrite !== undefined || options.ignoreIfExists !== undefined)) {
            result.options = options;
        }
        if (annotation !== undefined) {
            result.annotationId = annotation;
        }
        return result;
    }
    RenameFile.create = create;
    function is(value) {
        let candidate = value;
        return candidate && candidate.kind === 'rename' && Is.string(candidate.oldUri) && Is.string(candidate.newUri) && (candidate.options === undefined ||
            ((candidate.options.overwrite === undefined || Is.boolean(candidate.options.overwrite)) && (candidate.options.ignoreIfExists === undefined || Is.boolean(candidate.options.ignoreIfExists)))) && (candidate.annotationId === undefined || ChangeAnnotationIdentifier.is(candidate.annotationId));
    }
    RenameFile.is = is;
})(RenameFile || (RenameFile = {}));
var DeleteFile;
(function (DeleteFile) {
    function create(uri, options, annotation) {
        let result = {
            kind: 'delete',
            uri
        };
        if (options !== undefined && (options.recursive !== undefined || options.ignoreIfNotExists !== undefined)) {
            result.options = options;
        }
        if (annotation !== undefined) {
            result.annotationId = annotation;
        }
        return result;
    }
    DeleteFile.create = create;
    function is(value) {
        let candidate = value;
        return candidate && candidate.kind === 'delete' && Is.string(candidate.uri) && (candidate.options === undefined ||
            ((candidate.options.recursive === undefined || Is.boolean(candidate.options.recursive)) && (candidate.options.ignoreIfNotExists === undefined || Is.boolean(candidate.options.ignoreIfNotExists)))) && (candidate.annotationId === undefined || ChangeAnnotationIdentifier.is(candidate.annotationId));
    }
    DeleteFile.is = is;
})(DeleteFile || (DeleteFile = {}));
var WorkspaceEdit;
(function (WorkspaceEdit) {
    function is(value) {
        let candidate = value;
        return candidate &&
            (candidate.changes !== undefined || candidate.documentChanges !== undefined) &&
            (candidate.documentChanges === undefined || candidate.documentChanges.every((change) => {
                if (Is.string(change.kind)) {
                    return CreateFile.is(change) || RenameFile.is(change) || DeleteFile.is(change);
                }
                else {
                    return TextDocumentEdit.is(change);
                }
            }));
    }
    WorkspaceEdit.is = is;
})(WorkspaceEdit || (WorkspaceEdit = {}));
class TextEditChangeImpl {
    constructor(edits, changeAnnotations) {
        this.edits = edits;
        this.changeAnnotations = changeAnnotations;
    }
    insert(position, newText, annotation) {
        let edit;
        let id;
        if (annotation === undefined) {
            edit = TextEdit.insert(position, newText);
        }
        else if (ChangeAnnotationIdentifier.is(annotation)) {
            id = annotation;
            edit = AnnotatedTextEdit.insert(position, newText, annotation);
        }
        else {
            this.assertChangeAnnotations(this.changeAnnotations);
            id = this.changeAnnotations.manage(annotation);
            edit = AnnotatedTextEdit.insert(position, newText, id);
        }
        this.edits.push(edit);
        if (id !== undefined) {
            return id;
        }
    }
    replace(range, newText, annotation) {
        let edit;
        let id;
        if (annotation === undefined) {
            edit = TextEdit.replace(range, newText);
        }
        else if (ChangeAnnotationIdentifier.is(annotation)) {
            id = annotation;
            edit = AnnotatedTextEdit.replace(range, newText, annotation);
        }
        else {
            this.assertChangeAnnotations(this.changeAnnotations);
            id = this.changeAnnotations.manage(annotation);
            edit = AnnotatedTextEdit.replace(range, newText, id);
        }
        this.edits.push(edit);
        if (id !== undefined) {
            return id;
        }
    }
    delete(range, annotation) {
        let edit;
        let id;
        if (annotation === undefined) {
            edit = TextEdit.del(range);
        }
        else if (ChangeAnnotationIdentifier.is(annotation)) {
            id = annotation;
            edit = AnnotatedTextEdit.del(range, annotation);
        }
        else {
            this.assertChangeAnnotations(this.changeAnnotations);
            id = this.changeAnnotations.manage(annotation);
            edit = AnnotatedTextEdit.del(range, id);
        }
        this.edits.push(edit);
        if (id !== undefined) {
            return id;
        }
    }
    add(edit) {
        this.edits.push(edit);
    }
    all() {
        return this.edits;
    }
    clear() {
        this.edits.splice(0, this.edits.length);
    }
    assertChangeAnnotations(value) {
        if (value === undefined) {
            throw new Error(`Text edit change is not configured to manage change annotations.`);
        }
    }
}
/**
 * A helper class
 */
class ChangeAnnotations {
    constructor(annotations) {
        this._annotations = annotations === undefined ? Object.create(null) : annotations;
        this._counter = 0;
        this._size = 0;
    }
    all() {
        return this._annotations;
    }
    get size() {
        return this._size;
    }
    manage(idOrAnnotation, annotation) {
        let id;
        if (ChangeAnnotationIdentifier.is(idOrAnnotation)) {
            id = idOrAnnotation;
        }
        else {
            id = this.nextId();
            annotation = idOrAnnotation;
        }
        if (this._annotations[id] !== undefined) {
            throw new Error(`Id ${id} is already in use.`);
        }
        if (annotation === undefined) {
            throw new Error(`No annotation provided for id ${id}`);
        }
        this._annotations[id] = annotation;
        this._size++;
        return id;
    }
    nextId() {
        this._counter++;
        return this._counter.toString();
    }
}
/**
 * A workspace change helps constructing changes to a workspace.
 */
class WorkspaceChange {
    constructor(workspaceEdit) {
        this._textEditChanges = Object.create(null);
        if (workspaceEdit !== undefined) {
            this._workspaceEdit = workspaceEdit;
            if (workspaceEdit.documentChanges) {
                this._changeAnnotations = new ChangeAnnotations(workspaceEdit.changeAnnotations);
                workspaceEdit.changeAnnotations = this._changeAnnotations.all();
                workspaceEdit.documentChanges.forEach((change) => {
                    if (TextDocumentEdit.is(change)) {
                        const textEditChange = new TextEditChangeImpl(change.edits, this._changeAnnotations);
                        this._textEditChanges[change.textDocument.uri] = textEditChange;
                    }
                });
            }
            else if (workspaceEdit.changes) {
                Object.keys(workspaceEdit.changes).forEach((key) => {
                    const textEditChange = new TextEditChangeImpl(workspaceEdit.changes[key]);
                    this._textEditChanges[key] = textEditChange;
                });
            }
        }
        else {
            this._workspaceEdit = {};
        }
    }
    /**
     * Returns the underlying {@link WorkspaceEdit} literal
     * use to be returned from a workspace edit operation like rename.
     */
    get edit() {
        this.initDocumentChanges();
        if (this._changeAnnotations !== undefined) {
            if (this._changeAnnotations.size === 0) {
                this._workspaceEdit.changeAnnotations = undefined;
            }
            else {
                this._workspaceEdit.changeAnnotations = this._changeAnnotations.all();
            }
        }
        return this._workspaceEdit;
    }
    getTextEditChange(key) {
        if (OptionalVersionedTextDocumentIdentifier.is(key)) {
            this.initDocumentChanges();
            if (this._workspaceEdit.documentChanges === undefined) {
                throw new Error('Workspace edit is not configured for document changes.');
            }
            const textDocument = { uri: key.uri, version: key.version };
            let result = this._textEditChanges[textDocument.uri];
            if (!result) {
                const edits = [];
                const textDocumentEdit = {
                    textDocument,
                    edits
                };
                this._workspaceEdit.documentChanges.push(textDocumentEdit);
                result = new TextEditChangeImpl(edits, this._changeAnnotations);
                this._textEditChanges[textDocument.uri] = result;
            }
            return result;
        }
        else {
            this.initChanges();
            if (this._workspaceEdit.changes === undefined) {
                throw new Error('Workspace edit is not configured for normal text edit changes.');
            }
            let result = this._textEditChanges[key];
            if (!result) {
                let edits = [];
                this._workspaceEdit.changes[key] = edits;
                result = new TextEditChangeImpl(edits);
                this._textEditChanges[key] = result;
            }
            return result;
        }
    }
    initDocumentChanges() {
        if (this._workspaceEdit.documentChanges === undefined && this._workspaceEdit.changes === undefined) {
            this._changeAnnotations = new ChangeAnnotations();
            this._workspaceEdit.documentChanges = [];
            this._workspaceEdit.changeAnnotations = this._changeAnnotations.all();
        }
    }
    initChanges() {
        if (this._workspaceEdit.documentChanges === undefined && this._workspaceEdit.changes === undefined) {
            this._workspaceEdit.changes = Object.create(null);
        }
    }
    createFile(uri, optionsOrAnnotation, options) {
        this.initDocumentChanges();
        if (this._workspaceEdit.documentChanges === undefined) {
            throw new Error('Workspace edit is not configured for document changes.');
        }
        let annotation;
        if (ChangeAnnotation.is(optionsOrAnnotation) || ChangeAnnotationIdentifier.is(optionsOrAnnotation)) {
            annotation = optionsOrAnnotation;
        }
        else {
            options = optionsOrAnnotation;
        }
        let operation;
        let id;
        if (annotation === undefined) {
            operation = CreateFile.create(uri, options);
        }
        else {
            id = ChangeAnnotationIdentifier.is(annotation) ? annotation : this._changeAnnotations.manage(annotation);
            operation = CreateFile.create(uri, options, id);
        }
        this._workspaceEdit.documentChanges.push(operation);
        if (id !== undefined) {
            return id;
        }
    }
    renameFile(oldUri, newUri, optionsOrAnnotation, options) {
        this.initDocumentChanges();
        if (this._workspaceEdit.documentChanges === undefined) {
            throw new Error('Workspace edit is not configured for document changes.');
        }
        let annotation;
        if (ChangeAnnotation.is(optionsOrAnnotation) || ChangeAnnotationIdentifier.is(optionsOrAnnotation)) {
            annotation = optionsOrAnnotation;
        }
        else {
            options = optionsOrAnnotation;
        }
        let operation;
        let id;
        if (annotation === undefined) {
            operation = RenameFile.create(oldUri, newUri, options);
        }
        else {
            id = ChangeAnnotationIdentifier.is(annotation) ? annotation : this._changeAnnotations.manage(annotation);
            operation = RenameFile.create(oldUri, newUri, options, id);
        }
        this._workspaceEdit.documentChanges.push(operation);
        if (id !== undefined) {
            return id;
        }
    }
    deleteFile(uri, optionsOrAnnotation, options) {
        this.initDocumentChanges();
        if (this._workspaceEdit.documentChanges === undefined) {
            throw new Error('Workspace edit is not configured for document changes.');
        }
        let annotation;
        if (ChangeAnnotation.is(optionsOrAnnotation) || ChangeAnnotationIdentifier.is(optionsOrAnnotation)) {
            annotation = optionsOrAnnotation;
        }
        else {
            options = optionsOrAnnotation;
        }
        let operation;
        let id;
        if (annotation === undefined) {
            operation = DeleteFile.create(uri, options);
        }
        else {
            id = ChangeAnnotationIdentifier.is(annotation) ? annotation : this._changeAnnotations.manage(annotation);
            operation = DeleteFile.create(uri, options, id);
        }
        this._workspaceEdit.documentChanges.push(operation);
        if (id !== undefined) {
            return id;
        }
    }
}
/**
 * The TextDocumentIdentifier namespace provides helper functions to work with
 * {@link TextDocumentIdentifier} literals.
 */
var TextDocumentIdentifier;
(function (TextDocumentIdentifier) {
    /**
     * Creates a new TextDocumentIdentifier literal.
     * @param uri The document's uri.
     */
    function create(uri) {
        return { uri };
    }
    TextDocumentIdentifier.create = create;
    /**
     * Checks whether the given literal conforms to the {@link TextDocumentIdentifier} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.uri);
    }
    TextDocumentIdentifier.is = is;
})(TextDocumentIdentifier || (TextDocumentIdentifier = {}));
/**
 * The VersionedTextDocumentIdentifier namespace provides helper functions to work with
 * {@link VersionedTextDocumentIdentifier} literals.
 */
var VersionedTextDocumentIdentifier;
(function (VersionedTextDocumentIdentifier) {
    /**
     * Creates a new VersionedTextDocumentIdentifier literal.
     * @param uri The document's uri.
     * @param version The document's version.
     */
    function create(uri, version) {
        return { uri, version };
    }
    VersionedTextDocumentIdentifier.create = create;
    /**
     * Checks whether the given literal conforms to the {@link VersionedTextDocumentIdentifier} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.uri) && Is.integer(candidate.version);
    }
    VersionedTextDocumentIdentifier.is = is;
})(VersionedTextDocumentIdentifier || (VersionedTextDocumentIdentifier = {}));
/**
 * The OptionalVersionedTextDocumentIdentifier namespace provides helper functions to work with
 * {@link OptionalVersionedTextDocumentIdentifier} literals.
 */
var OptionalVersionedTextDocumentIdentifier;
(function (OptionalVersionedTextDocumentIdentifier) {
    /**
     * Creates a new OptionalVersionedTextDocumentIdentifier literal.
     * @param uri The document's uri.
     * @param version The document's version.
     */
    function create(uri, version) {
        return { uri, version };
    }
    OptionalVersionedTextDocumentIdentifier.create = create;
    /**
     * Checks whether the given literal conforms to the {@link OptionalVersionedTextDocumentIdentifier} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.uri) && (candidate.version === null || Is.integer(candidate.version));
    }
    OptionalVersionedTextDocumentIdentifier.is = is;
})(OptionalVersionedTextDocumentIdentifier || (OptionalVersionedTextDocumentIdentifier = {}));
/**
 * The TextDocumentItem namespace provides helper functions to work with
 * {@link TextDocumentItem} literals.
 */
var TextDocumentItem;
(function (TextDocumentItem) {
    /**
     * Creates a new TextDocumentItem literal.
     * @param uri The document's uri.
     * @param languageId The document's language identifier.
     * @param version The document's version number.
     * @param text The document's text.
     */
    function create(uri, languageId, version, text) {
        return { uri, languageId, version, text };
    }
    TextDocumentItem.create = create;
    /**
     * Checks whether the given literal conforms to the {@link TextDocumentItem} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.uri) && Is.string(candidate.languageId) && Is.integer(candidate.version) && Is.string(candidate.text);
    }
    TextDocumentItem.is = is;
})(TextDocumentItem || (TextDocumentItem = {}));
/**
 * Describes the content type that a client supports in various
 * result literals like `Hover`, `ParameterInfo` or `CompletionItem`.
 *
 * Please note that `MarkupKinds` must not start with a `$`. This kinds
 * are reserved for internal usage.
 */
var MarkupKind;
(function (MarkupKind) {
    /**
     * Plain text is supported as a content format
     */
    MarkupKind.PlainText = 'plaintext';
    /**
     * Markdown is supported as a content format
     */
    MarkupKind.Markdown = 'markdown';
    /**
     * Checks whether the given value is a value of the {@link MarkupKind} type.
     */
    function is(value) {
        const candidate = value;
        return candidate === MarkupKind.PlainText || candidate === MarkupKind.Markdown;
    }
    MarkupKind.is = is;
})(MarkupKind || (MarkupKind = {}));
var MarkupContent;
(function (MarkupContent) {
    /**
     * Checks whether the given value conforms to the {@link MarkupContent} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(value) && MarkupKind.is(candidate.kind) && Is.string(candidate.value);
    }
    MarkupContent.is = is;
})(MarkupContent || (MarkupContent = {}));
/**
 * The kind of a completion entry.
 */
var CompletionItemKind;
(function (CompletionItemKind) {
    CompletionItemKind.Text = 1;
    CompletionItemKind.Method = 2;
    CompletionItemKind.Function = 3;
    CompletionItemKind.Constructor = 4;
    CompletionItemKind.Field = 5;
    CompletionItemKind.Variable = 6;
    CompletionItemKind.Class = 7;
    CompletionItemKind.Interface = 8;
    CompletionItemKind.Module = 9;
    CompletionItemKind.Property = 10;
    CompletionItemKind.Unit = 11;
    CompletionItemKind.Value = 12;
    CompletionItemKind.Enum = 13;
    CompletionItemKind.Keyword = 14;
    CompletionItemKind.Snippet = 15;
    CompletionItemKind.Color = 16;
    CompletionItemKind.File = 17;
    CompletionItemKind.Reference = 18;
    CompletionItemKind.Folder = 19;
    CompletionItemKind.EnumMember = 20;
    CompletionItemKind.Constant = 21;
    CompletionItemKind.Struct = 22;
    CompletionItemKind.Event = 23;
    CompletionItemKind.Operator = 24;
    CompletionItemKind.TypeParameter = 25;
})(CompletionItemKind || (CompletionItemKind = {}));
/**
 * Defines whether the insert text in a completion item should be interpreted as
 * plain text or a snippet.
 */
var InsertTextFormat;
(function (InsertTextFormat) {
    /**
     * The primary text to be inserted is treated as a plain string.
     */
    InsertTextFormat.PlainText = 1;
    /**
     * The primary text to be inserted is treated as a snippet.
     *
     * A snippet can define tab stops and placeholders with `$1`, `$2`
     * and `${3:foo}`. `$0` defines the final tab stop, it defaults to
     * the end of the snippet. Placeholders with equal identifiers are linked,
     * that is typing in one will update others too.
     *
     * See also: https://microsoft.github.io/language-server-protocol/specifications/specification-current/#snippet_syntax
     */
    InsertTextFormat.Snippet = 2;
})(InsertTextFormat || (InsertTextFormat = {}));
/**
 * Completion item tags are extra annotations that tweak the rendering of a completion
 * item.
 *
 * @since 3.15.0
 */
var CompletionItemTag;
(function (CompletionItemTag) {
    /**
     * Render a completion as obsolete, usually using a strike-out.
     */
    CompletionItemTag.Deprecated = 1;
})(CompletionItemTag || (CompletionItemTag = {}));
/**
 * The InsertReplaceEdit namespace provides functions to deal with insert / replace edits.
 *
 * @since 3.16.0
 */
var InsertReplaceEdit;
(function (InsertReplaceEdit) {
    /**
     * Creates a new insert / replace edit
     */
    function create(newText, insert, replace) {
        return { newText, insert, replace };
    }
    InsertReplaceEdit.create = create;
    /**
     * Checks whether the given literal conforms to the {@link InsertReplaceEdit} interface.
     */
    function is(value) {
        const candidate = value;
        return candidate && Is.string(candidate.newText) && Range.is(candidate.insert) && Range.is(candidate.replace);
    }
    InsertReplaceEdit.is = is;
})(InsertReplaceEdit || (InsertReplaceEdit = {}));
/**
 * How whitespace and indentation is handled during completion
 * item insertion.
 *
 * @since 3.16.0
 */
var InsertTextMode;
(function (InsertTextMode) {
    /**
     * The insertion or replace strings is taken as it is. If the
     * value is multi line the lines below the cursor will be
     * inserted using the indentation defined in the string value.
     * The client will not apply any kind of adjustments to the
     * string.
     */
    InsertTextMode.asIs = 1;
    /**
     * The editor adjusts leading whitespace of new lines so that
     * they match the indentation up to the cursor of the line for
     * which the item is accepted.
     *
     * Consider a line like this: <2tabs><cursor><3tabs>foo. Accepting a
     * multi line completion item is indented using 2 tabs and all
     * following lines inserted will be indented using 2 tabs as well.
     */
    InsertTextMode.adjustIndentation = 2;
})(InsertTextMode || (InsertTextMode = {}));
var CompletionItemLabelDetails;
(function (CompletionItemLabelDetails) {
    function is(value) {
        const candidate = value;
        return candidate && (Is.string(candidate.detail) || candidate.detail === undefined) &&
            (Is.string(candidate.description) || candidate.description === undefined);
    }
    CompletionItemLabelDetails.is = is;
})(CompletionItemLabelDetails || (CompletionItemLabelDetails = {}));
/**
 * The CompletionItem namespace provides functions to deal with
 * completion items.
 */
var CompletionItem;
(function (CompletionItem) {
    /**
     * Create a completion item and seed it with a label.
     * @param label The completion item's label
     */
    function create(label) {
        return { label };
    }
    CompletionItem.create = create;
})(CompletionItem || (CompletionItem = {}));
/**
 * The CompletionList namespace provides functions to deal with
 * completion lists.
 */
var CompletionList;
(function (CompletionList) {
    /**
     * Creates a new completion list.
     *
     * @param items The completion items.
     * @param isIncomplete The list is not complete.
     */
    function create(items, isIncomplete) {
        return { items: items ? items : [], isIncomplete: !!isIncomplete };
    }
    CompletionList.create = create;
})(CompletionList || (CompletionList = {}));
var MarkedString;
(function (MarkedString) {
    /**
     * Creates a marked string from plain text.
     *
     * @param plainText The plain text.
     */
    function fromPlainText(plainText) {
        return plainText.replace(/[\\`*_{}[\]()#+\-.!]/g, '\\$&'); // escape markdown syntax tokens: http://daringfireball.net/projects/markdown/syntax#backslash
    }
    MarkedString.fromPlainText = fromPlainText;
    /**
     * Checks whether the given value conforms to the {@link MarkedString} type.
     */
    function is(value) {
        const candidate = value;
        return Is.string(candidate) || (Is.objectLiteral(candidate) && Is.string(candidate.language) && Is.string(candidate.value));
    }
    MarkedString.is = is;
})(MarkedString || (MarkedString = {}));
var Hover;
(function (Hover) {
    /**
     * Checks whether the given value conforms to the {@link Hover} interface.
     */
    function is(value) {
        let candidate = value;
        return !!candidate && Is.objectLiteral(candidate) && (MarkupContent.is(candidate.contents) ||
            MarkedString.is(candidate.contents) ||
            Is.typedArray(candidate.contents, MarkedString.is)) && (value.range === undefined || Range.is(value.range));
    }
    Hover.is = is;
})(Hover || (Hover = {}));
/**
 * The ParameterInformation namespace provides helper functions to work with
 * {@link ParameterInformation} literals.
 */
var ParameterInformation;
(function (ParameterInformation) {
    /**
     * Creates a new parameter information literal.
     *
     * @param label A label string.
     * @param documentation A doc string.
     */
    function create(label, documentation) {
        return documentation ? { label, documentation } : { label };
    }
    ParameterInformation.create = create;
})(ParameterInformation || (ParameterInformation = {}));
/**
 * The SignatureInformation namespace provides helper functions to work with
 * {@link SignatureInformation} literals.
 */
var SignatureInformation;
(function (SignatureInformation) {
    function create(label, documentation, ...parameters) {
        let result = { label };
        if (Is.defined(documentation)) {
            result.documentation = documentation;
        }
        if (Is.defined(parameters)) {
            result.parameters = parameters;
        }
        else {
            result.parameters = [];
        }
        return result;
    }
    SignatureInformation.create = create;
})(SignatureInformation || (SignatureInformation = {}));
/**
 * A document highlight kind.
 */
var DocumentHighlightKind;
(function (DocumentHighlightKind) {
    /**
     * A textual occurrence.
     */
    DocumentHighlightKind.Text = 1;
    /**
     * Read-access of a symbol, like reading a variable.
     */
    DocumentHighlightKind.Read = 2;
    /**
     * Write-access of a symbol, like writing to a variable.
     */
    DocumentHighlightKind.Write = 3;
})(DocumentHighlightKind || (DocumentHighlightKind = {}));
/**
 * DocumentHighlight namespace to provide helper functions to work with
 * {@link DocumentHighlight} literals.
 */
var DocumentHighlight;
(function (DocumentHighlight) {
    /**
     * Create a DocumentHighlight object.
     * @param range The range the highlight applies to.
     * @param kind The highlight kind
     */
    function create(range, kind) {
        let result = { range };
        if (Is.number(kind)) {
            result.kind = kind;
        }
        return result;
    }
    DocumentHighlight.create = create;
})(DocumentHighlight || (DocumentHighlight = {}));
/**
 * A symbol kind.
 */
var SymbolKind;
(function (SymbolKind) {
    SymbolKind.File = 1;
    SymbolKind.Module = 2;
    SymbolKind.Namespace = 3;
    SymbolKind.Package = 4;
    SymbolKind.Class = 5;
    SymbolKind.Method = 6;
    SymbolKind.Property = 7;
    SymbolKind.Field = 8;
    SymbolKind.Constructor = 9;
    SymbolKind.Enum = 10;
    SymbolKind.Interface = 11;
    SymbolKind.Function = 12;
    SymbolKind.Variable = 13;
    SymbolKind.Constant = 14;
    SymbolKind.String = 15;
    SymbolKind.Number = 16;
    SymbolKind.Boolean = 17;
    SymbolKind.Array = 18;
    SymbolKind.Object = 19;
    SymbolKind.Key = 20;
    SymbolKind.Null = 21;
    SymbolKind.EnumMember = 22;
    SymbolKind.Struct = 23;
    SymbolKind.Event = 24;
    SymbolKind.Operator = 25;
    SymbolKind.TypeParameter = 26;
})(SymbolKind || (SymbolKind = {}));
/**
 * Symbol tags are extra annotations that tweak the rendering of a symbol.
 *
 * @since 3.16
 */
var SymbolTag;
(function (SymbolTag) {
    /**
     * Render a symbol as obsolete, usually using a strike-out.
     */
    SymbolTag.Deprecated = 1;
})(SymbolTag || (SymbolTag = {}));
var SymbolInformation;
(function (SymbolInformation) {
    /**
     * Creates a new symbol information literal.
     *
     * @param name The name of the symbol.
     * @param kind The kind of the symbol.
     * @param range The range of the location of the symbol.
     * @param uri The resource of the location of symbol.
     * @param containerName The name of the symbol containing the symbol.
     */
    function create(name, kind, range, uri, containerName) {
        let result = {
            name,
            kind,
            location: { uri, range }
        };
        if (containerName) {
            result.containerName = containerName;
        }
        return result;
    }
    SymbolInformation.create = create;
})(SymbolInformation || (SymbolInformation = {}));
var WorkspaceSymbol;
(function (WorkspaceSymbol) {
    /**
     * Create a new workspace symbol.
     *
     * @param name The name of the symbol.
     * @param kind The kind of the symbol.
     * @param uri The resource of the location of the symbol.
     * @param range An options range of the location.
     * @returns A WorkspaceSymbol.
     */
    function create(name, kind, uri, range) {
        return range !== undefined
            ? { name, kind, location: { uri, range } }
            : { name, kind, location: { uri } };
    }
    WorkspaceSymbol.create = create;
})(WorkspaceSymbol || (WorkspaceSymbol = {}));
var DocumentSymbol;
(function (DocumentSymbol) {
    /**
     * Creates a new symbol information literal.
     *
     * @param name The name of the symbol.
     * @param detail The detail of the symbol.
     * @param kind The kind of the symbol.
     * @param range The range of the symbol.
     * @param selectionRange The selectionRange of the symbol.
     * @param children Children of the symbol.
     */
    function create(name, detail, kind, range, selectionRange, children) {
        let result = {
            name,
            detail,
            kind,
            range,
            selectionRange
        };
        if (children !== undefined) {
            result.children = children;
        }
        return result;
    }
    DocumentSymbol.create = create;
    /**
     * Checks whether the given literal conforms to the {@link DocumentSymbol} interface.
     */
    function is(value) {
        let candidate = value;
        return candidate &&
            Is.string(candidate.name) && Is.number(candidate.kind) &&
            Range.is(candidate.range) && Range.is(candidate.selectionRange) &&
            (candidate.detail === undefined || Is.string(candidate.detail)) &&
            (candidate.deprecated === undefined || Is.boolean(candidate.deprecated)) &&
            (candidate.children === undefined || Array.isArray(candidate.children)) &&
            (candidate.tags === undefined || Array.isArray(candidate.tags));
    }
    DocumentSymbol.is = is;
})(DocumentSymbol || (DocumentSymbol = {}));
/**
 * A set of predefined code action kinds
 */
var CodeActionKind;
(function (CodeActionKind) {
    /**
     * Empty kind.
     */
    CodeActionKind.Empty = '';
    /**
     * Base kind for quickfix actions: 'quickfix'
     */
    CodeActionKind.QuickFix = 'quickfix';
    /**
     * Base kind for refactoring actions: 'refactor'
     */
    CodeActionKind.Refactor = 'refactor';
    /**
     * Base kind for refactoring extraction actions: 'refactor.extract'
     *
     * Example extract actions:
     *
     * - Extract method
     * - Extract function
     * - Extract variable
     * - Extract interface from class
     * - ...
     */
    CodeActionKind.RefactorExtract = 'refactor.extract';
    /**
     * Base kind for refactoring inline actions: 'refactor.inline'
     *
     * Example inline actions:
     *
     * - Inline function
     * - Inline variable
     * - Inline constant
     * - ...
     */
    CodeActionKind.RefactorInline = 'refactor.inline';
    /**
     * Base kind for refactoring rewrite actions: 'refactor.rewrite'
     *
     * Example rewrite actions:
     *
     * - Convert JavaScript function to class
     * - Add or remove parameter
     * - Encapsulate field
     * - Make method static
     * - Move method to base class
     * - ...
     */
    CodeActionKind.RefactorRewrite = 'refactor.rewrite';
    /**
     * Base kind for source actions: `source`
     *
     * Source code actions apply to the entire file.
     */
    CodeActionKind.Source = 'source';
    /**
     * Base kind for an organize imports source action: `source.organizeImports`
     */
    CodeActionKind.SourceOrganizeImports = 'source.organizeImports';
    /**
     * Base kind for auto-fix source actions: `source.fixAll`.
     *
     * Fix all actions automatically fix errors that have a clear fix that do not require user input.
     * They should not suppress errors or perform unsafe fixes such as generating new types or classes.
     *
     * @since 3.15.0
     */
    CodeActionKind.SourceFixAll = 'source.fixAll';
})(CodeActionKind || (CodeActionKind = {}));
/**
 * The reason why code actions were requested.
 *
 * @since 3.17.0
 */
var CodeActionTriggerKind;
(function (CodeActionTriggerKind) {
    /**
     * Code actions were explicitly requested by the user or by an extension.
     */
    CodeActionTriggerKind.Invoked = 1;
    /**
     * Code actions were requested automatically.
     *
     * This typically happens when current selection in a file changes, but can
     * also be triggered when file content changes.
     */
    CodeActionTriggerKind.Automatic = 2;
})(CodeActionTriggerKind || (CodeActionTriggerKind = {}));
/**
 * The CodeActionContext namespace provides helper functions to work with
 * {@link CodeActionContext} literals.
 */
var CodeActionContext;
(function (CodeActionContext) {
    /**
     * Creates a new CodeActionContext literal.
     */
    function create(diagnostics, only, triggerKind) {
        let result = { diagnostics };
        if (only !== undefined && only !== null) {
            result.only = only;
        }
        if (triggerKind !== undefined && triggerKind !== null) {
            result.triggerKind = triggerKind;
        }
        return result;
    }
    CodeActionContext.create = create;
    /**
     * Checks whether the given literal conforms to the {@link CodeActionContext} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.typedArray(candidate.diagnostics, Diagnostic.is)
            && (candidate.only === undefined || Is.typedArray(candidate.only, Is.string))
            && (candidate.triggerKind === undefined || candidate.triggerKind === CodeActionTriggerKind.Invoked || candidate.triggerKind === CodeActionTriggerKind.Automatic);
    }
    CodeActionContext.is = is;
})(CodeActionContext || (CodeActionContext = {}));
var CodeAction;
(function (CodeAction) {
    function create(title, kindOrCommandOrEdit, kind) {
        let result = { title };
        let checkKind = true;
        if (typeof kindOrCommandOrEdit === 'string') {
            checkKind = false;
            result.kind = kindOrCommandOrEdit;
        }
        else if (Command.is(kindOrCommandOrEdit)) {
            result.command = kindOrCommandOrEdit;
        }
        else {
            result.edit = kindOrCommandOrEdit;
        }
        if (checkKind && kind !== undefined) {
            result.kind = kind;
        }
        return result;
    }
    CodeAction.create = create;
    function is(value) {
        let candidate = value;
        return candidate && Is.string(candidate.title) &&
            (candidate.diagnostics === undefined || Is.typedArray(candidate.diagnostics, Diagnostic.is)) &&
            (candidate.kind === undefined || Is.string(candidate.kind)) &&
            (candidate.edit !== undefined || candidate.command !== undefined) &&
            (candidate.command === undefined || Command.is(candidate.command)) &&
            (candidate.isPreferred === undefined || Is.boolean(candidate.isPreferred)) &&
            (candidate.edit === undefined || WorkspaceEdit.is(candidate.edit));
    }
    CodeAction.is = is;
})(CodeAction || (CodeAction = {}));
/**
 * The CodeLens namespace provides helper functions to work with
 * {@link CodeLens} literals.
 */
var CodeLens;
(function (CodeLens) {
    /**
     * Creates a new CodeLens literal.
     */
    function create(range, data) {
        let result = { range };
        if (Is.defined(data)) {
            result.data = data;
        }
        return result;
    }
    CodeLens.create = create;
    /**
     * Checks whether the given literal conforms to the {@link CodeLens} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Range.is(candidate.range) && (Is.undefined(candidate.command) || Command.is(candidate.command));
    }
    CodeLens.is = is;
})(CodeLens || (CodeLens = {}));
/**
 * The FormattingOptions namespace provides helper functions to work with
 * {@link FormattingOptions} literals.
 */
var FormattingOptions;
(function (FormattingOptions) {
    /**
     * Creates a new FormattingOptions literal.
     */
    function create(tabSize, insertSpaces) {
        return { tabSize, insertSpaces };
    }
    FormattingOptions.create = create;
    /**
     * Checks whether the given literal conforms to the {@link FormattingOptions} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.uinteger(candidate.tabSize) && Is.boolean(candidate.insertSpaces);
    }
    FormattingOptions.is = is;
})(FormattingOptions || (FormattingOptions = {}));
/**
 * The DocumentLink namespace provides helper functions to work with
 * {@link DocumentLink} literals.
 */
var DocumentLink;
(function (DocumentLink) {
    /**
     * Creates a new DocumentLink literal.
     */
    function create(range, target, data) {
        return { range, target, data };
    }
    DocumentLink.create = create;
    /**
     * Checks whether the given literal conforms to the {@link DocumentLink} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Range.is(candidate.range) && (Is.undefined(candidate.target) || Is.string(candidate.target));
    }
    DocumentLink.is = is;
})(DocumentLink || (DocumentLink = {}));
/**
 * The SelectionRange namespace provides helper function to work with
 * SelectionRange literals.
 */
var SelectionRange;
(function (SelectionRange) {
    /**
     * Creates a new SelectionRange
     * @param range the range.
     * @param parent an optional parent.
     */
    function create(range, parent) {
        return { range, parent };
    }
    SelectionRange.create = create;
    function is(value) {
        let candidate = value;
        return Is.objectLiteral(candidate) && Range.is(candidate.range) && (candidate.parent === undefined || SelectionRange.is(candidate.parent));
    }
    SelectionRange.is = is;
})(SelectionRange || (SelectionRange = {}));
/**
 * A set of predefined token types. This set is not fixed
 * an clients can specify additional token types via the
 * corresponding client capabilities.
 *
 * @since 3.16.0
 */
var SemanticTokenTypes;
(function (SemanticTokenTypes) {
    SemanticTokenTypes["namespace"] = "namespace";
    /**
     * Represents a generic type. Acts as a fallback for types which can't be mapped to
     * a specific type like class or enum.
     */
    SemanticTokenTypes["type"] = "type";
    SemanticTokenTypes["class"] = "class";
    SemanticTokenTypes["enum"] = "enum";
    SemanticTokenTypes["interface"] = "interface";
    SemanticTokenTypes["struct"] = "struct";
    SemanticTokenTypes["typeParameter"] = "typeParameter";
    SemanticTokenTypes["parameter"] = "parameter";
    SemanticTokenTypes["variable"] = "variable";
    SemanticTokenTypes["property"] = "property";
    SemanticTokenTypes["enumMember"] = "enumMember";
    SemanticTokenTypes["event"] = "event";
    SemanticTokenTypes["function"] = "function";
    SemanticTokenTypes["method"] = "method";
    SemanticTokenTypes["macro"] = "macro";
    SemanticTokenTypes["keyword"] = "keyword";
    SemanticTokenTypes["modifier"] = "modifier";
    SemanticTokenTypes["comment"] = "comment";
    SemanticTokenTypes["string"] = "string";
    SemanticTokenTypes["number"] = "number";
    SemanticTokenTypes["regexp"] = "regexp";
    SemanticTokenTypes["operator"] = "operator";
    /**
     * @since 3.17.0
     */
    SemanticTokenTypes["decorator"] = "decorator";
})(SemanticTokenTypes || (SemanticTokenTypes = {}));
/**
 * A set of predefined token modifiers. This set is not fixed
 * an clients can specify additional token types via the
 * corresponding client capabilities.
 *
 * @since 3.16.0
 */
var SemanticTokenModifiers;
(function (SemanticTokenModifiers) {
    SemanticTokenModifiers["declaration"] = "declaration";
    SemanticTokenModifiers["definition"] = "definition";
    SemanticTokenModifiers["readonly"] = "readonly";
    SemanticTokenModifiers["static"] = "static";
    SemanticTokenModifiers["deprecated"] = "deprecated";
    SemanticTokenModifiers["abstract"] = "abstract";
    SemanticTokenModifiers["async"] = "async";
    SemanticTokenModifiers["modification"] = "modification";
    SemanticTokenModifiers["documentation"] = "documentation";
    SemanticTokenModifiers["defaultLibrary"] = "defaultLibrary";
})(SemanticTokenModifiers || (SemanticTokenModifiers = {}));
/**
 * @since 3.16.0
 */
var SemanticTokens;
(function (SemanticTokens) {
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && (candidate.resultId === undefined || typeof candidate.resultId === 'string') &&
            Array.isArray(candidate.data) && (candidate.data.length === 0 || typeof candidate.data[0] === 'number');
    }
    SemanticTokens.is = is;
})(SemanticTokens || (SemanticTokens = {}));
/**
 * The InlineValueText namespace provides functions to deal with InlineValueTexts.
 *
 * @since 3.17.0
 */
var InlineValueText;
(function (InlineValueText) {
    /**
     * Creates a new InlineValueText literal.
     */
    function create(range, text) {
        return { range, text };
    }
    InlineValueText.create = create;
    function is(value) {
        const candidate = value;
        return candidate !== undefined && candidate !== null && Range.is(candidate.range) && Is.string(candidate.text);
    }
    InlineValueText.is = is;
})(InlineValueText || (InlineValueText = {}));
/**
 * The InlineValueVariableLookup namespace provides functions to deal with InlineValueVariableLookups.
 *
 * @since 3.17.0
 */
var InlineValueVariableLookup;
(function (InlineValueVariableLookup) {
    /**
     * Creates a new InlineValueText literal.
     */
    function create(range, variableName, caseSensitiveLookup) {
        return { range, variableName, caseSensitiveLookup };
    }
    InlineValueVariableLookup.create = create;
    function is(value) {
        const candidate = value;
        return candidate !== undefined && candidate !== null && Range.is(candidate.range) && Is.boolean(candidate.caseSensitiveLookup)
            && (Is.string(candidate.variableName) || candidate.variableName === undefined);
    }
    InlineValueVariableLookup.is = is;
})(InlineValueVariableLookup || (InlineValueVariableLookup = {}));
/**
 * The InlineValueEvaluatableExpression namespace provides functions to deal with InlineValueEvaluatableExpression.
 *
 * @since 3.17.0
 */
var InlineValueEvaluatableExpression;
(function (InlineValueEvaluatableExpression) {
    /**
     * Creates a new InlineValueEvaluatableExpression literal.
     */
    function create(range, expression) {
        return { range, expression };
    }
    InlineValueEvaluatableExpression.create = create;
    function is(value) {
        const candidate = value;
        return candidate !== undefined && candidate !== null && Range.is(candidate.range)
            && (Is.string(candidate.expression) || candidate.expression === undefined);
    }
    InlineValueEvaluatableExpression.is = is;
})(InlineValueEvaluatableExpression || (InlineValueEvaluatableExpression = {}));
/**
 * The InlineValueContext namespace provides helper functions to work with
 * {@link InlineValueContext} literals.
 *
 * @since 3.17.0
 */
var InlineValueContext;
(function (InlineValueContext) {
    /**
     * Creates a new InlineValueContext literal.
     */
    function create(frameId, stoppedLocation) {
        return { frameId, stoppedLocation };
    }
    InlineValueContext.create = create;
    /**
     * Checks whether the given literal conforms to the {@link InlineValueContext} interface.
     */
    function is(value) {
        const candidate = value;
        return Is.defined(candidate) && Range.is(value.stoppedLocation);
    }
    InlineValueContext.is = is;
})(InlineValueContext || (InlineValueContext = {}));
/**
 * Inlay hint kinds.
 *
 * @since 3.17.0
 */
var InlayHintKind;
(function (InlayHintKind) {
    /**
     * An inlay hint that for a type annotation.
     */
    InlayHintKind.Type = 1;
    /**
     * An inlay hint that is for a parameter.
     */
    InlayHintKind.Parameter = 2;
    function is(value) {
        return value === 1 || value === 2;
    }
    InlayHintKind.is = is;
})(InlayHintKind || (InlayHintKind = {}));
var InlayHintLabelPart;
(function (InlayHintLabelPart) {
    function create(value) {
        return { value };
    }
    InlayHintLabelPart.create = create;
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate)
            && (candidate.tooltip === undefined || Is.string(candidate.tooltip) || MarkupContent.is(candidate.tooltip))
            && (candidate.location === undefined || Location.is(candidate.location))
            && (candidate.command === undefined || Command.is(candidate.command));
    }
    InlayHintLabelPart.is = is;
})(InlayHintLabelPart || (InlayHintLabelPart = {}));
var InlayHint;
(function (InlayHint) {
    function create(position, label, kind) {
        const result = { position, label };
        if (kind !== undefined) {
            result.kind = kind;
        }
        return result;
    }
    InlayHint.create = create;
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && Position.is(candidate.position)
            && (Is.string(candidate.label) || Is.typedArray(candidate.label, InlayHintLabelPart.is))
            && (candidate.kind === undefined || InlayHintKind.is(candidate.kind))
            && (candidate.textEdits === undefined) || Is.typedArray(candidate.textEdits, TextEdit.is)
            && (candidate.tooltip === undefined || Is.string(candidate.tooltip) || MarkupContent.is(candidate.tooltip))
            && (candidate.paddingLeft === undefined || Is.boolean(candidate.paddingLeft))
            && (candidate.paddingRight === undefined || Is.boolean(candidate.paddingRight));
    }
    InlayHint.is = is;
})(InlayHint || (InlayHint = {}));
var StringValue;
(function (StringValue) {
    function createSnippet(value) {
        return { kind: 'snippet', value };
    }
    StringValue.createSnippet = createSnippet;
})(StringValue || (StringValue = {}));
var InlineCompletionItem;
(function (InlineCompletionItem) {
    function create(insertText, filterText, range, command) {
        return { insertText, filterText, range, command };
    }
    InlineCompletionItem.create = create;
})(InlineCompletionItem || (InlineCompletionItem = {}));
var InlineCompletionList;
(function (InlineCompletionList) {
    function create(items) {
        return { items };
    }
    InlineCompletionList.create = create;
})(InlineCompletionList || (InlineCompletionList = {}));
/**
 * Describes how an {@link InlineCompletionItemProvider inline completion provider} was triggered.
 *
 * @since 3.18.0
 * @proposed
 */
var InlineCompletionTriggerKind;
(function (InlineCompletionTriggerKind) {
    /**
     * Completion was triggered explicitly by a user gesture.
     */
    InlineCompletionTriggerKind.Invoked = 0;
    /**
     * Completion was triggered automatically while editing.
     */
    InlineCompletionTriggerKind.Automatic = 1;
})(InlineCompletionTriggerKind || (InlineCompletionTriggerKind = {}));
var SelectedCompletionInfo;
(function (SelectedCompletionInfo) {
    function create(range, text) {
        return { range, text };
    }
    SelectedCompletionInfo.create = create;
})(SelectedCompletionInfo || (SelectedCompletionInfo = {}));
var InlineCompletionContext;
(function (InlineCompletionContext) {
    function create(triggerKind, selectedCompletionInfo) {
        return { triggerKind, selectedCompletionInfo };
    }
    InlineCompletionContext.create = create;
})(InlineCompletionContext || (InlineCompletionContext = {}));
var WorkspaceFolder;
(function (WorkspaceFolder) {
    function is(value) {
        const candidate = value;
        return Is.objectLiteral(candidate) && URI.is(candidate.uri) && Is.string(candidate.name);
    }
    WorkspaceFolder.is = is;
})(WorkspaceFolder || (WorkspaceFolder = {}));
const EOL = (/* unused pure expression or super */ null && (['\n', '\r\n', '\r']));
/**
 * @deprecated Use the text document from the new vscode-languageserver-textdocument package.
 */
var TextDocument;
(function (TextDocument) {
    /**
     * Creates a new ITextDocument literal from the given uri and content.
     * @param uri The document's uri.
     * @param languageId The document's language Id.
     * @param version The document's version.
     * @param content The document's content.
     */
    function create(uri, languageId, version, content) {
        return new FullTextDocument(uri, languageId, version, content);
    }
    TextDocument.create = create;
    /**
     * Checks whether the given literal conforms to the {@link ITextDocument} interface.
     */
    function is(value) {
        let candidate = value;
        return Is.defined(candidate) && Is.string(candidate.uri) && (Is.undefined(candidate.languageId) || Is.string(candidate.languageId)) && Is.uinteger(candidate.lineCount)
            && Is.func(candidate.getText) && Is.func(candidate.positionAt) && Is.func(candidate.offsetAt) ? true : false;
    }
    TextDocument.is = is;
    function applyEdits(document, edits) {
        let text = document.getText();
        let sortedEdits = mergeSort(edits, (a, b) => {
            let diff = a.range.start.line - b.range.start.line;
            if (diff === 0) {
                return a.range.start.character - b.range.start.character;
            }
            return diff;
        });
        let lastModifiedOffset = text.length;
        for (let i = sortedEdits.length - 1; i >= 0; i--) {
            let e = sortedEdits[i];
            let startOffset = document.offsetAt(e.range.start);
            let endOffset = document.offsetAt(e.range.end);
            if (endOffset <= lastModifiedOffset) {
                text = text.substring(0, startOffset) + e.newText + text.substring(endOffset, text.length);
            }
            else {
                throw new Error('Overlapping edit');
            }
            lastModifiedOffset = startOffset;
        }
        return text;
    }
    TextDocument.applyEdits = applyEdits;
    function mergeSort(data, compare) {
        if (data.length <= 1) {
            // sorted
            return data;
        }
        const p = (data.length / 2) | 0;
        const left = data.slice(0, p);
        const right = data.slice(p);
        mergeSort(left, compare);
        mergeSort(right, compare);
        let leftIdx = 0;
        let rightIdx = 0;
        let i = 0;
        while (leftIdx < left.length && rightIdx < right.length) {
            let ret = compare(left[leftIdx], right[rightIdx]);
            if (ret <= 0) {
                // smaller_equal -> take left to preserve order
                data[i++] = left[leftIdx++];
            }
            else {
                // greater -> take right
                data[i++] = right[rightIdx++];
            }
        }
        while (leftIdx < left.length) {
            data[i++] = left[leftIdx++];
        }
        while (rightIdx < right.length) {
            data[i++] = right[rightIdx++];
        }
        return data;
    }
})(TextDocument || (TextDocument = {}));
/**
 * @deprecated Use the text document from the new vscode-languageserver-textdocument package.
 */
class FullTextDocument {
    constructor(uri, languageId, version, content) {
        this._uri = uri;
        this._languageId = languageId;
        this._version = version;
        this._content = content;
        this._lineOffsets = undefined;
    }
    get uri() {
        return this._uri;
    }
    get languageId() {
        return this._languageId;
    }
    get version() {
        return this._version;
    }
    getText(range) {
        if (range) {
            let start = this.offsetAt(range.start);
            let end = this.offsetAt(range.end);
            return this._content.substring(start, end);
        }
        return this._content;
    }
    update(event, version) {
        this._content = event.text;
        this._version = version;
        this._lineOffsets = undefined;
    }
    getLineOffsets() {
        if (this._lineOffsets === undefined) {
            let lineOffsets = [];
            let text = this._content;
            let isLineStart = true;
            for (let i = 0; i < text.length; i++) {
                if (isLineStart) {
                    lineOffsets.push(i);
                    isLineStart = false;
                }
                let ch = text.charAt(i);
                isLineStart = (ch === '\r' || ch === '\n');
                if (ch === '\r' && i + 1 < text.length && text.charAt(i + 1) === '\n') {
                    i++;
                }
            }
            if (isLineStart && text.length > 0) {
                lineOffsets.push(text.length);
            }
            this._lineOffsets = lineOffsets;
        }
        return this._lineOffsets;
    }
    positionAt(offset) {
        offset = Math.max(Math.min(offset, this._content.length), 0);
        let lineOffsets = this.getLineOffsets();
        let low = 0, high = lineOffsets.length;
        if (high === 0) {
            return Position.create(0, offset);
        }
        while (low < high) {
            let mid = Math.floor((low + high) / 2);
            if (lineOffsets[mid] > offset) {
                high = mid;
            }
            else {
                low = mid + 1;
            }
        }
        // low is the least x for which the line offset is larger than the current offset
        // or array.length if no line offset is larger than the current offset
        let line = low - 1;
        return Position.create(line, offset - lineOffsets[line]);
    }
    offsetAt(position) {
        let lineOffsets = this.getLineOffsets();
        if (position.line >= lineOffsets.length) {
            return this._content.length;
        }
        else if (position.line < 0) {
            return 0;
        }
        let lineOffset = lineOffsets[position.line];
        let nextLineOffset = (position.line + 1 < lineOffsets.length) ? lineOffsets[position.line + 1] : this._content.length;
        return Math.max(Math.min(lineOffset + position.character, nextLineOffset), lineOffset);
    }
    get lineCount() {
        return this.getLineOffsets().length;
    }
}
var Is;
(function (Is) {
    const toString = Object.prototype.toString;
    function defined(value) {
        return typeof value !== 'undefined';
    }
    Is.defined = defined;
    function undefined(value) {
        return typeof value === 'undefined';
    }
    Is.undefined = undefined;
    function boolean(value) {
        return value === true || value === false;
    }
    Is.boolean = boolean;
    function string(value) {
        return toString.call(value) === '[object String]';
    }
    Is.string = string;
    function number(value) {
        return toString.call(value) === '[object Number]';
    }
    Is.number = number;
    function numberRange(value, min, max) {
        return toString.call(value) === '[object Number]' && min <= value && value <= max;
    }
    Is.numberRange = numberRange;
    function integer(value) {
        return toString.call(value) === '[object Number]' && -2147483648 <= value && value <= 2147483647;
    }
    Is.integer = integer;
    function uinteger(value) {
        return toString.call(value) === '[object Number]' && 0 <= value && value <= 2147483647;
    }
    Is.uinteger = uinteger;
    function func(value) {
        return toString.call(value) === '[object Function]';
    }
    Is.func = func;
    function objectLiteral(value) {
        // Strictly speaking class instances pass this check as well. Since the LSP
        // doesn't use classes we ignore this for now. If we do we need to add something
        // like this: `Object.getPrototypeOf(Object.getPrototypeOf(x)) === null`
        return value !== null && typeof value === 'object';
    }
    Is.objectLiteral = objectLiteral;
    function typedArray(value, check) {
        return Array.isArray(value) && value.every(check);
    }
    Is.typedArray = typedArray;
})(Is || (Is = {}));

;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/cst-node-builder.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


class CstNodeBuilder {
    constructor() {
        this.nodeStack = [];
    }
    get current() {
        var _a;
        return (_a = this.nodeStack[this.nodeStack.length - 1]) !== null && _a !== void 0 ? _a : this.rootNode;
    }
    buildRootNode(input) {
        this.rootNode = new RootCstNodeImpl(input);
        this.rootNode.root = this.rootNode;
        this.nodeStack = [this.rootNode];
        return this.rootNode;
    }
    buildCompositeNode(feature) {
        const compositeNode = new CompositeCstNodeImpl();
        compositeNode.grammarSource = feature;
        compositeNode.root = this.rootNode;
        this.current.content.push(compositeNode);
        this.nodeStack.push(compositeNode);
        return compositeNode;
    }
    buildLeafNode(token, feature) {
        const leafNode = new LeafCstNodeImpl(token.startOffset, token.image.length, (0,cst_utils/* tokenToRange */.sp)(token), token.tokenType, !feature);
        leafNode.grammarSource = feature;
        leafNode.root = this.rootNode;
        this.current.content.push(leafNode);
        return leafNode;
    }
    removeNode(node) {
        const parent = node.container;
        if (parent) {
            const index = parent.content.indexOf(node);
            if (index >= 0) {
                parent.content.splice(index, 1);
            }
        }
    }
    addHiddenNodes(tokens) {
        const nodes = [];
        for (const token of tokens) {
            const leafNode = new LeafCstNodeImpl(token.startOffset, token.image.length, (0,cst_utils/* tokenToRange */.sp)(token), token.tokenType, true);
            leafNode.root = this.rootNode;
            nodes.push(leafNode);
        }
        let current = this.current;
        let added = false;
        // If we are within a composite node, we add the hidden nodes to the content
        if (current.content.length > 0) {
            current.content.push(...nodes);
            return;
        }
        // Otherwise we are at a newly created node
        // Instead of adding the hidden nodes here, we search for the first parent node with content
        while (current.container) {
            const index = current.container.content.indexOf(current);
            if (index > 0) {
                // Add the hidden nodes before the current node
                current.container.content.splice(index, 0, ...nodes);
                added = true;
                break;
            }
            current = current.container;
        }
        // If we arrive at the root node, we add the hidden nodes at the beginning
        // This is the case if the hidden nodes are the first nodes in the tree
        if (!added) {
            this.rootNode.content.unshift(...nodes);
        }
    }
    construct(item) {
        const current = this.current;
        // The specified item could be a datatype ($type is symbol) or a fragment ($type is undefined)
        // Only if the $type is a string, we actually assign the element
        if (typeof item.$type === 'string') {
            this.current.astNode = item;
        }
        item.$cstNode = current;
        const node = this.nodeStack.pop();
        // Empty composite nodes are not valid
        // Simply remove the node from the tree
        if ((node === null || node === void 0 ? void 0 : node.content.length) === 0) {
            this.removeNode(node);
        }
    }
}
class AbstractCstNode {
    /** @deprecated use `container` instead. */
    get parent() {
        return this.container;
    }
    /** @deprecated use `grammarSource` instead. */
    get feature() {
        return this.grammarSource;
    }
    get hidden() {
        return false;
    }
    get astNode() {
        var _a, _b;
        const node = typeof ((_a = this._astNode) === null || _a === void 0 ? void 0 : _a.$type) === 'string' ? this._astNode : (_b = this.container) === null || _b === void 0 ? void 0 : _b.astNode;
        if (!node) {
            throw new Error('This node has no associated AST element');
        }
        return node;
    }
    set astNode(value) {
        this._astNode = value;
    }
    /** @deprecated use `astNode` instead. */
    get element() {
        return this.astNode;
    }
    get text() {
        return this.root.fullText.substring(this.offset, this.end);
    }
}
class LeafCstNodeImpl extends AbstractCstNode {
    get offset() {
        return this._offset;
    }
    get length() {
        return this._length;
    }
    get end() {
        return this._offset + this._length;
    }
    get hidden() {
        return this._hidden;
    }
    get tokenType() {
        return this._tokenType;
    }
    get range() {
        return this._range;
    }
    constructor(offset, length, range, tokenType, hidden = false) {
        super();
        this._hidden = hidden;
        this._offset = offset;
        this._tokenType = tokenType;
        this._length = length;
        this._range = range;
    }
}
class CompositeCstNodeImpl extends AbstractCstNode {
    constructor() {
        super(...arguments);
        this.content = new CstNodeContainer(this);
    }
    /** @deprecated use `content` instead. */
    get children() {
        return this.content;
    }
    get offset() {
        var _a, _b;
        return (_b = (_a = this.firstNonHiddenNode) === null || _a === void 0 ? void 0 : _a.offset) !== null && _b !== void 0 ? _b : 0;
    }
    get length() {
        return this.end - this.offset;
    }
    get end() {
        var _a, _b;
        return (_b = (_a = this.lastNonHiddenNode) === null || _a === void 0 ? void 0 : _a.end) !== null && _b !== void 0 ? _b : 0;
    }
    get range() {
        const firstNode = this.firstNonHiddenNode;
        const lastNode = this.lastNonHiddenNode;
        if (firstNode && lastNode) {
            if (this._rangeCache === undefined) {
                const { range: firstRange } = firstNode;
                const { range: lastRange } = lastNode;
                this._rangeCache = { start: firstRange.start, end: lastRange.end.line < firstRange.start.line ? firstRange.start : lastRange.end };
            }
            return this._rangeCache;
        }
        else {
            return { start: Position.create(0, 0), end: Position.create(0, 0) };
        }
    }
    get firstNonHiddenNode() {
        for (const child of this.content) {
            if (!child.hidden) {
                return child;
            }
        }
        return this.content[0];
    }
    get lastNonHiddenNode() {
        for (let i = this.content.length - 1; i >= 0; i--) {
            const child = this.content[i];
            if (!child.hidden) {
                return child;
            }
        }
        return this.content[this.content.length - 1];
    }
}
class CstNodeContainer extends Array {
    constructor(parent) {
        super();
        this.parent = parent;
        Object.setPrototypeOf(this, CstNodeContainer.prototype);
    }
    push(...items) {
        this.addParents(items);
        return super.push(...items);
    }
    unshift(...items) {
        this.addParents(items);
        return super.unshift(...items);
    }
    splice(start, count, ...items) {
        this.addParents(items);
        return super.splice(start, count, ...items);
    }
    addParents(items) {
        for (const item of items) {
            item.container = this.parent;
        }
    }
}
class RootCstNodeImpl extends CompositeCstNodeImpl {
    get text() {
        return this._text.substring(this.offset, this.end);
    }
    get fullText() {
        return this._text;
    }
    constructor(input) {
        super();
        this._text = '';
        this._text = input !== null && input !== void 0 ? input : '';
    }
}
//# sourceMappingURL=cst-node-builder.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/langium-parser.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






const DatatypeSymbol = Symbol('Datatype');
function isDataTypeNode(node) {
    return node.$type === DatatypeSymbol;
}
const ruleSuffix = '\u200B';
const withRuleSuffix = (name) => name.endsWith(ruleSuffix) ? name : name + ruleSuffix;
class AbstractLangiumParser {
    constructor(services) {
        this._unorderedGroups = new Map();
        this.allRules = new Map();
        this.lexer = services.parser.Lexer;
        const tokens = this.lexer.definition;
        const production = services.LanguageMetaData.mode === 'production';
        this.wrapper = new ChevrotainWrapper(tokens, Object.assign(Object.assign({}, services.parser.ParserConfig), { skipValidations: production, errorMessageProvider: services.parser.ParserErrorMessageProvider }));
    }
    alternatives(idx, choices) {
        this.wrapper.wrapOr(idx, choices);
    }
    optional(idx, callback) {
        this.wrapper.wrapOption(idx, callback);
    }
    many(idx, callback) {
        this.wrapper.wrapMany(idx, callback);
    }
    atLeastOne(idx, callback) {
        this.wrapper.wrapAtLeastOne(idx, callback);
    }
    getRule(name) {
        return this.allRules.get(name);
    }
    isRecording() {
        return this.wrapper.IS_RECORDING;
    }
    get unorderedGroups() {
        return this._unorderedGroups;
    }
    getRuleStack() {
        return this.wrapper.RULE_STACK;
    }
    finalize() {
        this.wrapper.wrapSelfAnalysis();
    }
}
class LangiumParser extends AbstractLangiumParser {
    get current() {
        return this.stack[this.stack.length - 1];
    }
    constructor(services) {
        super(services);
        this.nodeBuilder = new CstNodeBuilder();
        this.stack = [];
        this.assignmentMap = new Map();
        this.linker = services.references.Linker;
        this.converter = services.parser.ValueConverter;
        this.astReflection = services.shared.AstReflection;
    }
    rule(rule, impl) {
        const type = this.computeRuleType(rule);
        const ruleMethod = this.wrapper.DEFINE_RULE(withRuleSuffix(rule.name), this.startImplementation(type, impl).bind(this));
        this.allRules.set(rule.name, ruleMethod);
        if (rule.entry) {
            this.mainRule = ruleMethod;
        }
        return ruleMethod;
    }
    computeRuleType(rule) {
        if (rule.fragment) {
            return undefined;
        }
        else if ((0,grammar_utils/* isDataTypeRule */.UP)(rule)) {
            return DatatypeSymbol;
        }
        else {
            const explicit = (0,grammar_utils/* getExplicitRuleType */.$G)(rule);
            return explicit !== null && explicit !== void 0 ? explicit : rule.name;
        }
    }
    parse(input, options = {}) {
        this.nodeBuilder.buildRootNode(input);
        const lexerResult = this.lexerResult = this.lexer.tokenize(input);
        this.wrapper.input = lexerResult.tokens;
        const ruleMethod = options.rule ? this.allRules.get(options.rule) : this.mainRule;
        if (!ruleMethod) {
            throw new Error(options.rule ? `No rule found with name '${options.rule}'` : 'No main rule available.');
        }
        const result = ruleMethod.call(this.wrapper, {});
        this.nodeBuilder.addHiddenNodes(lexerResult.hidden);
        this.unorderedGroups.clear();
        this.lexerResult = undefined;
        return {
            value: result,
            lexerErrors: lexerResult.errors,
            lexerReport: lexerResult.report,
            parserErrors: this.wrapper.errors
        };
    }
    startImplementation($type, implementation) {
        return (args) => {
            // Only create a new AST node in case the calling rule is not a fragment rule
            const createNode = !this.isRecording() && $type !== undefined;
            if (createNode) {
                const node = { $type };
                this.stack.push(node);
                if ($type === DatatypeSymbol) {
                    node.value = '';
                }
            }
            let result;
            try {
                result = implementation(args);
            }
            catch (err) {
                result = undefined;
            }
            if (result === undefined && createNode) {
                result = this.construct();
            }
            return result;
        };
    }
    extractHiddenTokens(token) {
        const hiddenTokens = this.lexerResult.hidden;
        if (!hiddenTokens.length) {
            return [];
        }
        const offset = token.startOffset;
        for (let i = 0; i < hiddenTokens.length; i++) {
            const token = hiddenTokens[i];
            if (token.startOffset > offset) {
                return hiddenTokens.splice(0, i);
            }
        }
        return hiddenTokens.splice(0, hiddenTokens.length);
    }
    consume(idx, tokenType, feature) {
        const token = this.wrapper.wrapConsume(idx, tokenType);
        if (!this.isRecording() && this.isValidToken(token)) {
            const hiddenTokens = this.extractHiddenTokens(token);
            this.nodeBuilder.addHiddenNodes(hiddenTokens);
            const leafNode = this.nodeBuilder.buildLeafNode(token, feature);
            const { assignment, isCrossRef } = this.getAssignment(feature);
            const current = this.current;
            if (assignment) {
                const convertedValue = (0,ast/* isKeyword */.p1)(feature) ? token.image : this.converter.convert(token.image, leafNode);
                this.assign(assignment.operator, assignment.feature, convertedValue, leafNode, isCrossRef);
            }
            else if (isDataTypeNode(current)) {
                let text = token.image;
                if (!(0,ast/* isKeyword */.p1)(feature)) {
                    text = this.converter.convert(text, leafNode).toString();
                }
                current.value += text;
            }
        }
    }
    /**
     * Most consumed parser tokens are valid. However there are two cases in which they are not valid:
     *
     * 1. They were inserted during error recovery by the parser. These tokens don't really exist and should not be further processed
     * 2. They contain invalid token ranges. This might include the special EOF token, or other tokens produced by invalid token builders.
     */
    isValidToken(token) {
        return !token.isInsertedInRecovery && !isNaN(token.startOffset) && typeof token.endOffset === 'number' && !isNaN(token.endOffset);
    }
    subrule(idx, rule, fragment, feature, args) {
        let cstNode;
        if (!this.isRecording() && !fragment) {
            // We only want to create a new CST node if the subrule actually creates a new AST node.
            // In other cases like calls of fragment rules the current CST/AST is populated further.
            // Note that skipping this initialization and leaving cstNode unassigned also skips the subrule assignment later on.
            // This is intended, as fragment rules only enrich the current AST node
            cstNode = this.nodeBuilder.buildCompositeNode(feature);
        }
        const subruleResult = this.wrapper.wrapSubrule(idx, rule, args);
        if (!this.isRecording() && cstNode && cstNode.length > 0) {
            this.performSubruleAssignment(subruleResult, feature, cstNode);
        }
    }
    performSubruleAssignment(result, feature, cstNode) {
        const { assignment, isCrossRef } = this.getAssignment(feature);
        if (assignment) {
            this.assign(assignment.operator, assignment.feature, result, cstNode, isCrossRef);
        }
        else if (!assignment) {
            // If we call a subrule without an assignment we either:
            // 1. append the result of the subrule (data type rule)
            // 2. override the current object with the newly parsed object
            // If the current element is an AST node and the result of the subrule
            // is a data type rule, we can safely discard the results.
            const current = this.current;
            if (isDataTypeNode(current)) {
                current.value += result.toString();
            }
            else if (typeof result === 'object' && result) {
                const object = this.assignWithoutOverride(result, current);
                const newItem = object;
                this.stack.pop();
                this.stack.push(newItem);
            }
        }
    }
    action($type, action) {
        if (!this.isRecording()) {
            let last = this.current;
            if (action.feature && action.operator) {
                last = this.construct();
                this.nodeBuilder.removeNode(last.$cstNode);
                const node = this.nodeBuilder.buildCompositeNode(action);
                node.content.push(last.$cstNode);
                const newItem = { $type };
                this.stack.push(newItem);
                this.assign(action.operator, action.feature, last, last.$cstNode, false);
            }
            else {
                last.$type = $type;
            }
        }
    }
    construct() {
        if (this.isRecording()) {
            return undefined;
        }
        const obj = this.current;
        (0,ast_utils/* linkContentToContainer */.b2)(obj);
        this.nodeBuilder.construct(obj);
        this.stack.pop();
        if (isDataTypeNode(obj)) {
            return this.converter.convert(obj.value, obj.$cstNode);
        }
        else {
            (0,ast_utils/* assignMandatoryProperties */.a1)(this.astReflection, obj);
        }
        return obj;
    }
    getAssignment(feature) {
        if (!this.assignmentMap.has(feature)) {
            const assignment = (0,ast_utils/* getContainerOfType */.V_)(feature, ast/* isAssignment */.B7);
            this.assignmentMap.set(feature, {
                assignment: assignment,
                isCrossRef: assignment ? (0,ast/* isCrossReference */.Ki)(assignment.terminal) : false
            });
        }
        return this.assignmentMap.get(feature);
    }
    assign(operator, feature, value, cstNode, isCrossRef) {
        const obj = this.current;
        let item;
        if (isCrossRef && typeof value === 'string') {
            item = this.linker.buildReference(obj, feature, cstNode, value);
        }
        else {
            item = value;
        }
        switch (operator) {
            case '=': {
                obj[feature] = item;
                break;
            }
            case '?=': {
                obj[feature] = true;
                break;
            }
            case '+=': {
                if (!Array.isArray(obj[feature])) {
                    obj[feature] = [];
                }
                obj[feature].push(item);
            }
        }
    }
    assignWithoutOverride(target, source) {
        for (const [name, existingValue] of Object.entries(source)) {
            const newValue = target[name];
            if (newValue === undefined) {
                target[name] = existingValue;
            }
            else if (Array.isArray(newValue) && Array.isArray(existingValue)) {
                existingValue.push(...newValue);
                target[name] = existingValue;
            }
        }
        // The target was parsed from a unassigned subrule
        // After the subrule construction, it received a cst node
        // This CST node will later be overriden by the cst node builder
        // To prevent references to stale AST nodes in the CST,
        // we need to remove the reference here
        const targetCstNode = target.$cstNode;
        if (targetCstNode) {
            targetCstNode.astNode = undefined;
            target.$cstNode = undefined;
        }
        return target;
    }
    get definitionErrors() {
        return this.wrapper.definitionErrors;
    }
}
class AbstractParserErrorMessageProvider {
    buildMismatchTokenMessage(options) {
        return api/* defaultParserErrorProvider */.Hs.buildMismatchTokenMessage(options);
    }
    buildNotAllInputParsedMessage(options) {
        return api/* defaultParserErrorProvider */.Hs.buildNotAllInputParsedMessage(options);
    }
    buildNoViableAltMessage(options) {
        return api/* defaultParserErrorProvider */.Hs.buildNoViableAltMessage(options);
    }
    buildEarlyExitMessage(options) {
        return api/* defaultParserErrorProvider */.Hs.buildEarlyExitMessage(options);
    }
}
class LangiumParserErrorMessageProvider extends AbstractParserErrorMessageProvider {
    buildMismatchTokenMessage({ expected, actual }) {
        const expectedMsg = expected.LABEL
            ? '`' + expected.LABEL + '`'
            : expected.name.endsWith(':KW')
                ? `keyword '${expected.name.substring(0, expected.name.length - 3)}'`
                : `token of type '${expected.name}'`;
        return `Expecting ${expectedMsg} but found \`${actual.image}\`.`;
    }
    buildNotAllInputParsedMessage({ firstRedundant }) {
        return `Expecting end of file but found \`${firstRedundant.image}\`.`;
    }
}
class LangiumCompletionParser extends AbstractLangiumParser {
    constructor() {
        super(...arguments);
        this.tokens = [];
        this.elementStack = [];
        this.lastElementStack = [];
        this.nextTokenIndex = 0;
        this.stackSize = 0;
    }
    action() {
        // NOOP
    }
    construct() {
        // NOOP
        return undefined;
    }
    parse(input) {
        this.resetState();
        const tokens = this.lexer.tokenize(input, { mode: 'partial' });
        this.tokens = tokens.tokens;
        this.wrapper.input = [...this.tokens];
        this.mainRule.call(this.wrapper, {});
        this.unorderedGroups.clear();
        return {
            tokens: this.tokens,
            elementStack: [...this.lastElementStack],
            tokenIndex: this.nextTokenIndex
        };
    }
    rule(rule, impl) {
        const ruleMethod = this.wrapper.DEFINE_RULE(withRuleSuffix(rule.name), this.startImplementation(impl).bind(this));
        this.allRules.set(rule.name, ruleMethod);
        if (rule.entry) {
            this.mainRule = ruleMethod;
        }
        return ruleMethod;
    }
    resetState() {
        this.elementStack = [];
        this.lastElementStack = [];
        this.nextTokenIndex = 0;
        this.stackSize = 0;
    }
    startImplementation(implementation) {
        return (args) => {
            const size = this.keepStackSize();
            try {
                implementation(args);
            }
            finally {
                this.resetStackSize(size);
            }
        };
    }
    removeUnexpectedElements() {
        this.elementStack.splice(this.stackSize);
    }
    keepStackSize() {
        const size = this.elementStack.length;
        this.stackSize = size;
        return size;
    }
    resetStackSize(size) {
        this.removeUnexpectedElements();
        this.stackSize = size;
    }
    consume(idx, tokenType, feature) {
        this.wrapper.wrapConsume(idx, tokenType);
        if (!this.isRecording()) {
            this.lastElementStack = [...this.elementStack, feature];
            this.nextTokenIndex = this.currIdx + 1;
        }
    }
    subrule(idx, rule, fragment, feature, args) {
        this.before(feature);
        this.wrapper.wrapSubrule(idx, rule, args);
        this.after(feature);
    }
    before(element) {
        if (!this.isRecording()) {
            this.elementStack.push(element);
        }
    }
    after(element) {
        if (!this.isRecording()) {
            const index = this.elementStack.lastIndexOf(element);
            if (index >= 0) {
                this.elementStack.splice(index);
            }
        }
    }
    get currIdx() {
        return this.wrapper.currIdx;
    }
}
const defaultConfig = {
    recoveryEnabled: true,
    nodeLocationTracking: 'full',
    skipValidations: true,
    errorMessageProvider: new LangiumParserErrorMessageProvider()
};
/**
 * This class wraps the embedded actions parser of chevrotain and exposes protected methods.
 * This way, we can build the `LangiumParser` as a composition.
 */
class ChevrotainWrapper extends api/* EmbeddedActionsParser */.nu {
    constructor(tokens, config) {
        const useDefaultLookahead = config && 'maxLookahead' in config;
        super(tokens, Object.assign(Object.assign(Object.assign({}, defaultConfig), { lookaheadStrategy: useDefaultLookahead
                ? new api/* LLkLookaheadStrategy */.dV({ maxLookahead: config.maxLookahead })
                : new LLStarLookaheadStrategy({
                    // If validations are skipped, don't log the lookahead warnings
                    logging: config.skipValidations ? () => { } : undefined
                }) }), config));
    }
    get IS_RECORDING() {
        return this.RECORDING_PHASE;
    }
    DEFINE_RULE(name, impl) {
        return this.RULE(name, impl);
    }
    wrapSelfAnalysis() {
        this.performSelfAnalysis();
    }
    wrapConsume(idx, tokenType) {
        return this.consume(idx, tokenType);
    }
    wrapSubrule(idx, rule, args) {
        return this.subrule(idx, rule, {
            ARGS: [args]
        });
    }
    wrapOr(idx, choices) {
        this.or(idx, choices);
    }
    wrapOption(idx, callback) {
        this.option(idx, callback);
    }
    wrapMany(idx, callback) {
        this.many(idx, callback);
    }
    wrapAtLeastOne(idx, callback) {
        this.atLeastOne(idx, callback);
    }
}
//# sourceMappingURL=langium-parser.js.map
// EXTERNAL MODULE: ../node_modules/langium/lib/utils/errors.js
var errors = __webpack_require__(45209);
// EXTERNAL MODULE: ../node_modules/langium/lib/utils/stream.js
var stream = __webpack_require__(99293);
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/parser-builder-base.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/





function createParser(grammar, parser, tokens) {
    const parserContext = {
        parser,
        tokens,
        ruleNames: new Map()
    };
    buildRules(parserContext, grammar);
    return parser;
}
function buildRules(parserContext, grammar) {
    const reachable = (0,grammar_utils/* getAllReachableRules */.VD)(grammar, false);
    const parserRules = (0,stream/* stream */.Vw)(grammar.rules).filter(ast/* isParserRule */.F9).filter(rule => reachable.has(rule));
    for (const rule of parserRules) {
        const ctx = Object.assign(Object.assign({}, parserContext), { consume: 1, optional: 1, subrule: 1, many: 1, or: 1 });
        parserContext.parser.rule(rule, buildElement(ctx, rule.definition));
    }
}
function buildElement(ctx, element, ignoreGuard = false) {
    let method;
    if ((0,ast/* isKeyword */.p1)(element)) {
        method = buildKeyword(ctx, element);
    }
    else if ((0,ast/* isAction */.LG)(element)) {
        method = buildAction(ctx, element);
    }
    else if ((0,ast/* isAssignment */.B7)(element)) {
        method = buildElement(ctx, element.terminal);
    }
    else if ((0,ast/* isCrossReference */.Ki)(element)) {
        method = buildCrossReference(ctx, element);
    }
    else if ((0,ast/* isRuleCall */.t3)(element)) {
        method = buildRuleCall(ctx, element);
    }
    else if ((0,ast/* isAlternatives */.MZ)(element)) {
        method = buildAlternatives(ctx, element);
    }
    else if ((0,ast/* isUnorderedGroup */.W1)(element)) {
        method = buildUnorderedGroup(ctx, element);
    }
    else if ((0,ast/* isGroup */.ty)(element)) {
        method = buildGroup(ctx, element);
    }
    else if ((0,ast/* isEndOfFile */.rT)(element)) {
        const idx = ctx.consume++;
        method = () => ctx.parser.consume(idx, api/* EOF */.sd, element);
    }
    else {
        throw new errors/* ErrorWithLocation */.h(element.$cstNode, `Unexpected element type: ${element.$type}`);
    }
    return wrap(ctx, ignoreGuard ? undefined : getGuardCondition(element), method, element.cardinality);
}
function buildAction(ctx, action) {
    const actionType = (0,grammar_utils/* getTypeName */.z$)(action);
    return () => ctx.parser.action(actionType, action);
}
function buildRuleCall(ctx, ruleCall) {
    const rule = ruleCall.rule.ref;
    if ((0,ast/* isParserRule */.F9)(rule)) {
        const idx = ctx.subrule++;
        const fragment = rule.fragment;
        const predicate = ruleCall.arguments.length > 0 ? buildRuleCallPredicate(rule, ruleCall.arguments) : () => ({});
        return (args) => ctx.parser.subrule(idx, getRule(ctx, rule), fragment, ruleCall, predicate(args));
    }
    else if ((0,ast/* isTerminalRule */.MS)(rule)) {
        const idx = ctx.consume++;
        const method = getToken(ctx, rule.name);
        return () => ctx.parser.consume(idx, method, ruleCall);
    }
    else if (!rule) {
        throw new errors/* ErrorWithLocation */.h(ruleCall.$cstNode, `Undefined rule: ${ruleCall.rule.$refText}`);
    }
    else {
        (0,errors/* assertUnreachable */.U)(rule);
    }
}
function buildRuleCallPredicate(rule, namedArgs) {
    const predicates = namedArgs.map(e => buildPredicate(e.value));
    return (args) => {
        const ruleArgs = {};
        for (let i = 0; i < predicates.length; i++) {
            const ruleTarget = rule.parameters[i];
            const predicate = predicates[i];
            ruleArgs[ruleTarget.name] = predicate(args);
        }
        return ruleArgs;
    };
}
function buildPredicate(condition) {
    if ((0,ast/* isDisjunction */.F8)(condition)) {
        const left = buildPredicate(condition.left);
        const right = buildPredicate(condition.right);
        return (args) => (left(args) || right(args));
    }
    else if ((0,ast/* isConjunction */.TB)(condition)) {
        const left = buildPredicate(condition.left);
        const right = buildPredicate(condition.right);
        return (args) => (left(args) && right(args));
    }
    else if ((0,ast/* isNegation */.Ii)(condition)) {
        const value = buildPredicate(condition.value);
        return (args) => !value(args);
    }
    else if ((0,ast/* isParameterReference */.yW)(condition)) {
        const name = condition.parameter.ref.name;
        return (args) => args !== undefined && args[name] === true;
    }
    else if ((0,ast/* isBooleanLiteral */.L)(condition)) {
        const value = Boolean(condition.true);
        return () => value;
    }
    (0,errors/* assertUnreachable */.U)(condition);
}
function buildAlternatives(ctx, alternatives) {
    if (alternatives.elements.length === 1) {
        return buildElement(ctx, alternatives.elements[0]);
    }
    else {
        const methods = [];
        for (const element of alternatives.elements) {
            const predicatedMethod = {
                // Since we handle the guard condition in the alternative already
                // We can ignore the group guard condition inside
                ALT: buildElement(ctx, element, true)
            };
            const guard = getGuardCondition(element);
            if (guard) {
                predicatedMethod.GATE = buildPredicate(guard);
            }
            methods.push(predicatedMethod);
        }
        const idx = ctx.or++;
        return (args) => ctx.parser.alternatives(idx, methods.map(method => {
            const alt = {
                ALT: () => method.ALT(args)
            };
            const gate = method.GATE;
            if (gate) {
                alt.GATE = () => gate(args);
            }
            return alt;
        }));
    }
}
function buildUnorderedGroup(ctx, group) {
    if (group.elements.length === 1) {
        return buildElement(ctx, group.elements[0]);
    }
    const methods = [];
    for (const element of group.elements) {
        const predicatedMethod = {
            // Since we handle the guard condition in the alternative already
            // We can ignore the group guard condition inside
            ALT: buildElement(ctx, element, true)
        };
        const guard = getGuardCondition(element);
        if (guard) {
            predicatedMethod.GATE = buildPredicate(guard);
        }
        methods.push(predicatedMethod);
    }
    const orIdx = ctx.or++;
    const idFunc = (groupIdx, lParser) => {
        const stackId = lParser.getRuleStack().join('-');
        return `uGroup_${groupIdx}_${stackId}`;
    };
    const alternatives = (args) => ctx.parser.alternatives(orIdx, methods.map((method, idx) => {
        const alt = { ALT: () => true };
        const parser = ctx.parser;
        alt.ALT = () => {
            method.ALT(args);
            if (!parser.isRecording()) {
                const key = idFunc(orIdx, parser);
                if (!parser.unorderedGroups.get(key)) {
                    // init after clear state
                    parser.unorderedGroups.set(key, []);
                }
                const groupState = parser.unorderedGroups.get(key);
                if (typeof (groupState === null || groupState === void 0 ? void 0 : groupState[idx]) === 'undefined') {
                    // Not accessed yet
                    groupState[idx] = true;
                }
            }
        };
        const gate = method.GATE;
        if (gate) {
            alt.GATE = () => gate(args);
        }
        else {
            alt.GATE = () => {
                const trackedAlternatives = parser.unorderedGroups.get(idFunc(orIdx, parser));
                const allow = !(trackedAlternatives === null || trackedAlternatives === void 0 ? void 0 : trackedAlternatives[idx]);
                return allow;
            };
        }
        return alt;
    }));
    const wrapped = wrap(ctx, getGuardCondition(group), alternatives, '*');
    return (args) => {
        wrapped(args);
        if (!ctx.parser.isRecording()) {
            ctx.parser.unorderedGroups.delete(idFunc(orIdx, ctx.parser));
        }
    };
}
function buildGroup(ctx, group) {
    const methods = group.elements.map(e => buildElement(ctx, e));
    return (args) => methods.forEach(method => method(args));
}
function getGuardCondition(element) {
    if ((0,ast/* isGroup */.ty)(element)) {
        return element.guardCondition;
    }
    return undefined;
}
function buildCrossReference(ctx, crossRef, terminal = crossRef.terminal) {
    if (!terminal) {
        if (!crossRef.type.ref) {
            throw new Error('Could not resolve reference to type: ' + crossRef.type.$refText);
        }
        const assignment = (0,grammar_utils/* findNameAssignment */.ib)(crossRef.type.ref);
        const assignTerminal = assignment === null || assignment === void 0 ? void 0 : assignment.terminal;
        if (!assignTerminal) {
            throw new Error('Could not find name assignment for type: ' + (0,grammar_utils/* getTypeName */.z$)(crossRef.type.ref));
        }
        return buildCrossReference(ctx, crossRef, assignTerminal);
    }
    else if ((0,ast/* isRuleCall */.t3)(terminal) && (0,ast/* isParserRule */.F9)(terminal.rule.ref)) {
        // The terminal is a data type rule here. Everything else will result in a validation error.
        const rule = terminal.rule.ref;
        const idx = ctx.subrule++;
        return (args) => ctx.parser.subrule(idx, getRule(ctx, rule), false, crossRef, args);
    }
    else if ((0,ast/* isRuleCall */.t3)(terminal) && (0,ast/* isTerminalRule */.MS)(terminal.rule.ref)) {
        const idx = ctx.consume++;
        const terminalRule = getToken(ctx, terminal.rule.ref.name);
        return () => ctx.parser.consume(idx, terminalRule, crossRef);
    }
    else if ((0,ast/* isKeyword */.p1)(terminal)) {
        const idx = ctx.consume++;
        const keyword = getToken(ctx, terminal.value);
        return () => ctx.parser.consume(idx, keyword, crossRef);
    }
    else {
        throw new Error('Could not build cross reference parser');
    }
}
function buildKeyword(ctx, keyword) {
    const idx = ctx.consume++;
    const token = ctx.tokens[keyword.value];
    if (!token) {
        throw new Error('Could not find token for keyword: ' + keyword.value);
    }
    return () => ctx.parser.consume(idx, token, keyword);
}
function wrap(ctx, guard, method, cardinality) {
    const gate = guard && buildPredicate(guard);
    if (!cardinality) {
        if (gate) {
            const idx = ctx.or++;
            return (args) => ctx.parser.alternatives(idx, [
                {
                    ALT: () => method(args),
                    GATE: () => gate(args)
                },
                {
                    ALT: (0,api/* EMPTY_ALT */._o)(),
                    GATE: () => !gate(args)
                }
            ]);
        }
        else {
            return method;
        }
    }
    if (cardinality === '*') {
        const idx = ctx.many++;
        return (args) => ctx.parser.many(idx, {
            DEF: () => method(args),
            GATE: gate ? () => gate(args) : undefined
        });
    }
    else if (cardinality === '+') {
        const idx = ctx.many++;
        if (gate) {
            const orIdx = ctx.or++;
            // In the case of a guard condition for the `+` group
            // We combine it with an empty alternative
            // If the condition returns true, it needs to parse at least a single iteration
            // If its false, it is not allowed to parse anything
            return (args) => ctx.parser.alternatives(orIdx, [
                {
                    ALT: () => ctx.parser.atLeastOne(idx, {
                        DEF: () => method(args)
                    }),
                    GATE: () => gate(args)
                },
                {
                    ALT: (0,api/* EMPTY_ALT */._o)(),
                    GATE: () => !gate(args)
                }
            ]);
        }
        else {
            return (args) => ctx.parser.atLeastOne(idx, {
                DEF: () => method(args),
            });
        }
    }
    else if (cardinality === '?') {
        const idx = ctx.optional++;
        return (args) => ctx.parser.optional(idx, {
            DEF: () => method(args),
            GATE: gate ? () => gate(args) : undefined
        });
    }
    else {
        (0,errors/* assertUnreachable */.U)(cardinality);
    }
}
function getRule(ctx, element) {
    const name = getRuleName(ctx, element);
    const rule = ctx.parser.getRule(name);
    if (!rule)
        throw new Error(`Rule "${name}" not found."`);
    return rule;
}
function getRuleName(ctx, element) {
    if ((0,ast/* isParserRule */.F9)(element)) {
        return element.name;
    }
    else if (ctx.ruleNames.has(element)) {
        return ctx.ruleNames.get(element);
    }
    else {
        let item = element;
        let parent = item.$container;
        let ruleName = element.$type;
        while (!(0,ast/* isParserRule */.F9)(parent)) {
            if ((0,ast/* isGroup */.ty)(parent) || (0,ast/* isAlternatives */.MZ)(parent) || (0,ast/* isUnorderedGroup */.W1)(parent)) {
                const index = parent.elements.indexOf(item);
                ruleName = index.toString() + ':' + ruleName;
            }
            item = parent;
            parent = parent.$container;
        }
        const rule = parent;
        ruleName = rule.name + ':' + ruleName;
        ctx.ruleNames.set(element, ruleName);
        return ruleName;
    }
}
function getToken(ctx, name) {
    const token = ctx.tokens[name];
    if (!token)
        throw new Error(`Token "${name}" not found."`);
    return token;
}
//# sourceMappingURL=parser-builder-base.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/completion-parser-builder.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


function createCompletionParser(services) {
    const grammar = services.Grammar;
    const lexer = services.parser.Lexer;
    const parser = new LangiumCompletionParser(services);
    createParser(grammar, parser, lexer.definition);
    parser.finalize();
    return parser;
}
//# sourceMappingURL=completion-parser-builder.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/langium-parser-builder.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


/**
 * Create and finalize a Langium parser. The parser rules are derived from the grammar, which is
 * available at `services.Grammar`.
 */
function createLangiumParser(services) {
    const parser = prepareLangiumParser(services);
    parser.finalize();
    return parser;
}
/**
 * Create a Langium parser without finalizing it. This is used to extract more detailed error
 * information when the parser is initially validated.
 */
function prepareLangiumParser(services) {
    const grammar = services.Grammar;
    const lexer = services.parser.Lexer;
    const parser = new LangiumParser(services);
    return createParser(grammar, parser, lexer.definition);
}
//# sourceMappingURL=langium-parser-builder.js.map
// EXTERNAL MODULE: ../node_modules/langium/lib/parser/token-builder.js
var token_builder = __webpack_require__(35481);
// EXTERNAL MODULE: ../node_modules/langium/lib/parser/value-converter.js
var value_converter = __webpack_require__(46174);
// EXTERNAL MODULE: ../node_modules/vscode-jsonrpc/lib/common/cancellation.js
var cancellation = __webpack_require__(97770);
// EXTERNAL MODULE: ../node_modules/langium/lib/syntax-tree.js
var syntax_tree = __webpack_require__(91303);
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/promise-utils.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

/**
 * Delays the execution of the current code to the next tick of the event loop.
 * Don't call this method directly in a tight loop to prevent too many promises from being created.
 */
function delayNextTick() {
    return new Promise(resolve => {
        // In case we are running in a non-node environment, `setImmediate` isn't available.
        // Using `setTimeout` of the browser API accomplishes the same result.
        if (typeof setImmediate === 'undefined') {
            setTimeout(resolve, 0);
        }
        else {
            setImmediate(resolve);
        }
    });
}
let lastTick = 0;
let globalInterruptionPeriod = 10;
/**
 * Reset the global interruption period and create a cancellation token source.
 */
function startCancelableOperation() {
    lastTick = performance.now();
    return new cancellation.CancellationTokenSource();
}
/**
 * Change the period duration for `interruptAndCheck` to the given number of milliseconds.
 * The default value is 10ms.
 */
function setInterruptionPeriod(period) {
    globalInterruptionPeriod = period;
}
/**
 * This symbol may be thrown in an asynchronous context by any Langium service that receives
 * a `CancellationToken`. This means that the promise returned by such a service is rejected with
 * this symbol as rejection reason.
 */
const promise_utils_OperationCancelled = Symbol('OperationCancelled');
/**
 * Use this in a `catch` block to check whether the thrown object indicates that the operation
 * has been cancelled.
 */
function isOperationCancelled(err) {
    return err === promise_utils_OperationCancelled;
}
/**
 * This function does two things:
 *  1. Check the elapsed time since the last call to this function or to `startCancelableOperation`. If the predefined
 *     period (configured with `setInterruptionPeriod`) is exceeded, execution is delayed with `delayNextTick`.
 *  2. If the predefined period is not met yet or execution is resumed after an interruption, the given cancellation
 *     token is checked, and if cancellation is requested, `OperationCanceled` is thrown.
 *
 * All services in Langium that receive a `CancellationToken` may potentially call this function, so the
 * `CancellationToken` must be caught (with an `async` try-catch block or a `catch` callback attached to
 * the promise) to avoid that event being exposed as an error.
 */
async function interruptAndCheck(token) {
    if (token === cancellation.CancellationToken.None) {
        // Early exit in case cancellation was disabled by the caller
        return;
    }
    const current = performance.now();
    if (current - lastTick >= globalInterruptionPeriod) {
        lastTick = current;
        await delayNextTick();
        // prevent calling delayNextTick every iteration of loop
        // where delayNextTick takes up the majority or all of the
        // globalInterruptionPeriod itself
        lastTick = performance.now();
    }
    if (token.isCancellationRequested) {
        throw promise_utils_OperationCancelled;
    }
}
/**
 * Simple implementation of the deferred pattern.
 * An object that exposes a promise and functions to resolve and reject it.
 */
class promise_utils_Deferred {
    constructor() {
        this.promise = new Promise((resolve, reject) => {
            this.resolve = (arg) => {
                resolve(arg);
                return this;
            };
            this.reject = (err) => {
                reject(err);
                return this;
            };
        });
    }
}
//# sourceMappingURL=promise-utils.js.map
;// CONCATENATED MODULE: ../node_modules/vscode-languageserver-textdocument/lib/esm/main.js
/* --------------------------------------------------------------------------------------------
 * Copyright (c) Microsoft Corporation. All rights reserved.
 * Licensed under the MIT License. See License.txt in the project root for license information.
 * ------------------------------------------------------------------------------------------ */

class main_FullTextDocument {
    constructor(uri, languageId, version, content) {
        this._uri = uri;
        this._languageId = languageId;
        this._version = version;
        this._content = content;
        this._lineOffsets = undefined;
    }
    get uri() {
        return this._uri;
    }
    get languageId() {
        return this._languageId;
    }
    get version() {
        return this._version;
    }
    getText(range) {
        if (range) {
            const start = this.offsetAt(range.start);
            const end = this.offsetAt(range.end);
            return this._content.substring(start, end);
        }
        return this._content;
    }
    update(changes, version) {
        for (const change of changes) {
            if (main_FullTextDocument.isIncremental(change)) {
                // makes sure start is before end
                const range = getWellformedRange(change.range);
                // update content
                const startOffset = this.offsetAt(range.start);
                const endOffset = this.offsetAt(range.end);
                this._content = this._content.substring(0, startOffset) + change.text + this._content.substring(endOffset, this._content.length);
                // update the offsets
                const startLine = Math.max(range.start.line, 0);
                const endLine = Math.max(range.end.line, 0);
                let lineOffsets = this._lineOffsets;
                const addedLineOffsets = computeLineOffsets(change.text, false, startOffset);
                if (endLine - startLine === addedLineOffsets.length) {
                    for (let i = 0, len = addedLineOffsets.length; i < len; i++) {
                        lineOffsets[i + startLine + 1] = addedLineOffsets[i];
                    }
                }
                else {
                    if (addedLineOffsets.length < 10000) {
                        lineOffsets.splice(startLine + 1, endLine - startLine, ...addedLineOffsets);
                    }
                    else { // avoid too many arguments for splice
                        this._lineOffsets = lineOffsets = lineOffsets.slice(0, startLine + 1).concat(addedLineOffsets, lineOffsets.slice(endLine + 1));
                    }
                }
                const diff = change.text.length - (endOffset - startOffset);
                if (diff !== 0) {
                    for (let i = startLine + 1 + addedLineOffsets.length, len = lineOffsets.length; i < len; i++) {
                        lineOffsets[i] = lineOffsets[i] + diff;
                    }
                }
            }
            else if (main_FullTextDocument.isFull(change)) {
                this._content = change.text;
                this._lineOffsets = undefined;
            }
            else {
                throw new Error('Unknown change event received');
            }
        }
        this._version = version;
    }
    getLineOffsets() {
        if (this._lineOffsets === undefined) {
            this._lineOffsets = computeLineOffsets(this._content, true);
        }
        return this._lineOffsets;
    }
    positionAt(offset) {
        offset = Math.max(Math.min(offset, this._content.length), 0);
        const lineOffsets = this.getLineOffsets();
        let low = 0, high = lineOffsets.length;
        if (high === 0) {
            return { line: 0, character: offset };
        }
        while (low < high) {
            const mid = Math.floor((low + high) / 2);
            if (lineOffsets[mid] > offset) {
                high = mid;
            }
            else {
                low = mid + 1;
            }
        }
        // low is the least x for which the line offset is larger than the current offset
        // or array.length if no line offset is larger than the current offset
        const line = low - 1;
        offset = this.ensureBeforeEOL(offset, lineOffsets[line]);
        return { line, character: offset - lineOffsets[line] };
    }
    offsetAt(position) {
        const lineOffsets = this.getLineOffsets();
        if (position.line >= lineOffsets.length) {
            return this._content.length;
        }
        else if (position.line < 0) {
            return 0;
        }
        const lineOffset = lineOffsets[position.line];
        if (position.character <= 0) {
            return lineOffset;
        }
        const nextLineOffset = (position.line + 1 < lineOffsets.length) ? lineOffsets[position.line + 1] : this._content.length;
        const offset = Math.min(lineOffset + position.character, nextLineOffset);
        return this.ensureBeforeEOL(offset, lineOffset);
    }
    ensureBeforeEOL(offset, lineOffset) {
        while (offset > lineOffset && isEOL(this._content.charCodeAt(offset - 1))) {
            offset--;
        }
        return offset;
    }
    get lineCount() {
        return this.getLineOffsets().length;
    }
    static isIncremental(event) {
        const candidate = event;
        return candidate !== undefined && candidate !== null &&
            typeof candidate.text === 'string' && candidate.range !== undefined &&
            (candidate.rangeLength === undefined || typeof candidate.rangeLength === 'number');
    }
    static isFull(event) {
        const candidate = event;
        return candidate !== undefined && candidate !== null &&
            typeof candidate.text === 'string' && candidate.range === undefined && candidate.rangeLength === undefined;
    }
}
var main_TextDocument;
(function (TextDocument) {
    /**
     * Creates a new text document.
     *
     * @param uri The document's uri.
     * @param languageId  The document's language Id.
     * @param version The document's initial version number.
     * @param content The document's content.
     */
    function create(uri, languageId, version, content) {
        return new main_FullTextDocument(uri, languageId, version, content);
    }
    TextDocument.create = create;
    /**
     * Updates a TextDocument by modifying its content.
     *
     * @param document the document to update. Only documents created by TextDocument.create are valid inputs.
     * @param changes the changes to apply to the document.
     * @param version the changes version for the document.
     * @returns The updated TextDocument. Note: That's the same document instance passed in as first parameter.
     *
     */
    function update(document, changes, version) {
        if (document instanceof main_FullTextDocument) {
            document.update(changes, version);
            return document;
        }
        else {
            throw new Error('TextDocument.update: document must be created by TextDocument.create');
        }
    }
    TextDocument.update = update;
    function applyEdits(document, edits) {
        const text = document.getText();
        const sortedEdits = mergeSort(edits.map(getWellformedEdit), (a, b) => {
            const diff = a.range.start.line - b.range.start.line;
            if (diff === 0) {
                return a.range.start.character - b.range.start.character;
            }
            return diff;
        });
        let lastModifiedOffset = 0;
        const spans = [];
        for (const e of sortedEdits) {
            const startOffset = document.offsetAt(e.range.start);
            if (startOffset < lastModifiedOffset) {
                throw new Error('Overlapping edit');
            }
            else if (startOffset > lastModifiedOffset) {
                spans.push(text.substring(lastModifiedOffset, startOffset));
            }
            if (e.newText.length) {
                spans.push(e.newText);
            }
            lastModifiedOffset = document.offsetAt(e.range.end);
        }
        spans.push(text.substr(lastModifiedOffset));
        return spans.join('');
    }
    TextDocument.applyEdits = applyEdits;
})(main_TextDocument || (main_TextDocument = {}));
function mergeSort(data, compare) {
    if (data.length <= 1) {
        // sorted
        return data;
    }
    const p = (data.length / 2) | 0;
    const left = data.slice(0, p);
    const right = data.slice(p);
    mergeSort(left, compare);
    mergeSort(right, compare);
    let leftIdx = 0;
    let rightIdx = 0;
    let i = 0;
    while (leftIdx < left.length && rightIdx < right.length) {
        const ret = compare(left[leftIdx], right[rightIdx]);
        if (ret <= 0) {
            // smaller_equal -> take left to preserve order
            data[i++] = left[leftIdx++];
        }
        else {
            // greater -> take right
            data[i++] = right[rightIdx++];
        }
    }
    while (leftIdx < left.length) {
        data[i++] = left[leftIdx++];
    }
    while (rightIdx < right.length) {
        data[i++] = right[rightIdx++];
    }
    return data;
}
function computeLineOffsets(text, isAtLineStart, textOffset = 0) {
    const result = isAtLineStart ? [textOffset] : [];
    for (let i = 0; i < text.length; i++) {
        const ch = text.charCodeAt(i);
        if (isEOL(ch)) {
            if (ch === 13 /* CharCode.CarriageReturn */ && i + 1 < text.length && text.charCodeAt(i + 1) === 10 /* CharCode.LineFeed */) {
                i++;
            }
            result.push(textOffset + i + 1);
        }
    }
    return result;
}
function isEOL(char) {
    return char === 13 /* CharCode.CarriageReturn */ || char === 10 /* CharCode.LineFeed */;
}
function getWellformedRange(range) {
    const start = range.start;
    const end = range.end;
    if (start.line > end.line || (start.line === end.line && start.character > end.character)) {
        return { start: end, end: start };
    }
    return range;
}
function getWellformedEdit(textEdit) {
    const range = getWellformedRange(textEdit.range);
    if (range !== textEdit.range) {
        return { newText: textEdit.newText, range };
    }
    return textEdit;
}

// EXTERNAL MODULE: ../node_modules/vscode-uri/lib/esm/index.mjs
var esm = __webpack_require__(37943);
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/documents.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
/**
 * Re-export 'TextDocument' from 'vscode-languageserver-textdocument' for convenience,
 *  including both type _and_ symbol (namespace), as we here and there also refer to the symbol,
 *  the overhead is very small, just a few kilobytes.
 * Everything else of that package (at the time contributing) is also defined
 *  in 'vscode-languageserver-protocol' or 'vscode-languageserver-types'.
 */





/**
 * A document is subject to several phases that are run in predefined order. Any state value implies that
 * smaller state values are finished as well.
 */
var DocumentState;
(function (DocumentState) {
    /**
     * The text content has changed and needs to be parsed again. The AST held by this outdated
     * document instance is no longer valid.
     */
    DocumentState[DocumentState["Changed"] = 0] = "Changed";
    /**
     * An AST has been created from the text content. The document structure can be traversed,
     * but cross-references cannot be resolved yet. If necessary, the structure can be manipulated
     * at this stage as a preprocessing step.
     */
    DocumentState[DocumentState["Parsed"] = 1] = "Parsed";
    /**
     * The `IndexManager` service has processed AST nodes of this document. This means the
     * exported symbols are available in the global scope and can be resolved from other documents.
     */
    DocumentState[DocumentState["IndexedContent"] = 2] = "IndexedContent";
    /**
     * The `ScopeComputation` service has processed this document. This means the local symbols
     * are stored in a MultiMap so they can be looked up by the `ScopeProvider` service.
     * Once a document has reached this state, you may follow every reference - it will lazily
     * resolve its `ref` property and yield either the target AST node or `undefined` in case
     * the target is not in scope.
     */
    DocumentState[DocumentState["ComputedScopes"] = 3] = "ComputedScopes";
    /**
     * The `Linker` service has processed this document. All outgoing references have been
     * resolved or marked as erroneous.
     */
    DocumentState[DocumentState["Linked"] = 4] = "Linked";
    /**
     * The `IndexManager` service has processed AST node references of this document. This is
     * necessary to determine which documents are affected by a change in one of the workspace
     * documents.
     */
    DocumentState[DocumentState["IndexedReferences"] = 5] = "IndexedReferences";
    /**
     * The `DocumentValidator` service has processed this document. The language server listens
     * to the results of this phase and sends diagnostics to the client.
     */
    DocumentState[DocumentState["Validated"] = 6] = "Validated";
})(DocumentState || (DocumentState = {}));
class DefaultLangiumDocumentFactory {
    constructor(services) {
        this.serviceRegistry = services.ServiceRegistry;
        this.textDocuments = services.workspace.TextDocuments;
        this.fileSystemProvider = services.workspace.FileSystemProvider;
    }
    async fromUri(uri, cancellationToken = cancellation.CancellationToken.None) {
        const content = await this.fileSystemProvider.readFile(uri);
        return this.createAsync(uri, content, cancellationToken);
    }
    fromTextDocument(textDocument, uri, token) {
        uri = uri !== null && uri !== void 0 ? uri : esm/* URI */.o.parse(textDocument.uri);
        if (cancellation.CancellationToken.is(token)) {
            return this.createAsync(uri, textDocument, token);
        }
        else {
            return this.create(uri, textDocument, token);
        }
    }
    fromString(text, uri, token) {
        if (cancellation.CancellationToken.is(token)) {
            return this.createAsync(uri, text, token);
        }
        else {
            return this.create(uri, text, token);
        }
    }
    fromModel(model, uri) {
        return this.create(uri, { $model: model });
    }
    create(uri, content, options) {
        if (typeof content === 'string') {
            const parseResult = this.parse(uri, content, options);
            return this.createLangiumDocument(parseResult, uri, undefined, content);
        }
        else if ('$model' in content) {
            const parseResult = { value: content.$model, parserErrors: [], lexerErrors: [] };
            return this.createLangiumDocument(parseResult, uri);
        }
        else {
            const parseResult = this.parse(uri, content.getText(), options);
            return this.createLangiumDocument(parseResult, uri, content);
        }
    }
    async createAsync(uri, content, cancelToken) {
        if (typeof content === 'string') {
            const parseResult = await this.parseAsync(uri, content, cancelToken);
            return this.createLangiumDocument(parseResult, uri, undefined, content);
        }
        else {
            const parseResult = await this.parseAsync(uri, content.getText(), cancelToken);
            return this.createLangiumDocument(parseResult, uri, content);
        }
    }
    /**
     * Create a LangiumDocument from a given parse result.
     *
     * A TextDocument is created on demand if it is not provided as argument here. Usually this
     * should not be necessary because the main purpose of the TextDocument is to convert between
     * text ranges and offsets, which is done solely in LSP request handling.
     *
     * With the introduction of {@link update} below this method is supposed to be mainly called
     * during workspace initialization and on addition/recognition of new files, while changes in
     * existing documents are processed via {@link update}.
     */
    createLangiumDocument(parseResult, uri, textDocument, text) {
        let document;
        if (textDocument) {
            document = {
                parseResult,
                uri,
                state: DocumentState.Parsed,
                references: [],
                textDocument
            };
        }
        else {
            const textDocumentGetter = this.createTextDocumentGetter(uri, text);
            document = {
                parseResult,
                uri,
                state: DocumentState.Parsed,
                references: [],
                get textDocument() {
                    return textDocumentGetter();
                }
            };
        }
        parseResult.value.$document = document;
        return document;
    }
    async update(document, cancellationToken) {
        var _a, _b;
        // The CST full text property contains the original text that was used to create the AST.
        const oldText = (_a = document.parseResult.value.$cstNode) === null || _a === void 0 ? void 0 : _a.root.fullText;
        const textDocument = (_b = this.textDocuments) === null || _b === void 0 ? void 0 : _b.get(document.uri.toString());
        const text = textDocument ? textDocument.getText() : await this.fileSystemProvider.readFile(document.uri);
        if (textDocument) {
            Object.defineProperty(document, 'textDocument', {
                value: textDocument
            });
        }
        else {
            const textDocumentGetter = this.createTextDocumentGetter(document.uri, text);
            Object.defineProperty(document, 'textDocument', {
                get: textDocumentGetter
            });
        }
        // Some of these documents can be pretty large, so parsing them again can be quite expensive.
        // Therefore, we only parse if the text has actually changed.
        if (oldText !== text) {
            document.parseResult = await this.parseAsync(document.uri, text, cancellationToken);
            document.parseResult.value.$document = document;
        }
        document.state = DocumentState.Parsed;
        return document;
    }
    parse(uri, text, options) {
        const services = this.serviceRegistry.getServices(uri);
        return services.parser.LangiumParser.parse(text, options);
    }
    parseAsync(uri, text, cancellationToken) {
        const services = this.serviceRegistry.getServices(uri);
        return services.parser.AsyncParser.parse(text, cancellationToken);
    }
    createTextDocumentGetter(uri, text) {
        const serviceRegistry = this.serviceRegistry;
        let textDoc = undefined;
        return () => {
            return textDoc !== null && textDoc !== void 0 ? textDoc : (textDoc = main_TextDocument.create(uri.toString(), serviceRegistry.getServices(uri).LanguageMetaData.languageId, 0, text !== null && text !== void 0 ? text : ''));
        };
    }
}
class DefaultLangiumDocuments {
    constructor(services) {
        this.documentMap = new Map();
        this.langiumDocumentFactory = services.workspace.LangiumDocumentFactory;
        this.serviceRegistry = services.ServiceRegistry;
    }
    get all() {
        return (0,stream/* stream */.Vw)(this.documentMap.values());
    }
    addDocument(document) {
        const uriString = document.uri.toString();
        if (this.documentMap.has(uriString)) {
            throw new Error(`A document with the URI '${uriString}' is already present.`);
        }
        this.documentMap.set(uriString, document);
    }
    getDocument(uri) {
        const uriString = uri.toString();
        return this.documentMap.get(uriString);
    }
    async getOrCreateDocument(uri, cancellationToken) {
        let document = this.getDocument(uri);
        if (document) {
            return document;
        }
        document = await this.langiumDocumentFactory.fromUri(uri, cancellationToken);
        this.addDocument(document);
        return document;
    }
    createDocument(uri, text, cancellationToken) {
        if (cancellationToken) {
            return this.langiumDocumentFactory.fromString(text, uri, cancellationToken).then(document => {
                this.addDocument(document);
                return document;
            });
        }
        else {
            const document = this.langiumDocumentFactory.fromString(text, uri);
            this.addDocument(document);
            return document;
        }
    }
    hasDocument(uri) {
        return this.documentMap.has(uri.toString());
    }
    invalidateDocument(uri) {
        const uriString = uri.toString();
        const langiumDoc = this.documentMap.get(uriString);
        if (langiumDoc) {
            const linker = this.serviceRegistry.getServices(uri).references.Linker;
            linker.unlink(langiumDoc);
            langiumDoc.state = DocumentState.Changed;
            langiumDoc.precomputedScopes = undefined;
            langiumDoc.diagnostics = undefined;
        }
        return langiumDoc;
    }
    deleteDocument(uri) {
        const uriString = uri.toString();
        const langiumDoc = this.documentMap.get(uriString);
        if (langiumDoc) {
            langiumDoc.state = DocumentState.Changed;
            this.documentMap.delete(uriString);
        }
        return langiumDoc;
    }
}
//# sourceMappingURL=documents.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/linker.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/





const ref_resolving = Symbol('ref_resolving');
class DefaultLinker {
    constructor(services) {
        this.reflection = services.shared.AstReflection;
        this.langiumDocuments = () => services.shared.workspace.LangiumDocuments;
        this.scopeProvider = services.references.ScopeProvider;
        this.astNodeLocator = services.workspace.AstNodeLocator;
    }
    async link(document, cancelToken = cancellation.CancellationToken.None) {
        for (const node of (0,ast_utils/* streamAst */.Zc)(document.parseResult.value)) {
            await interruptAndCheck(cancelToken);
            (0,ast_utils/* streamReferences */.fy)(node).forEach(ref => this.doLink(ref, document));
        }
    }
    doLink(refInfo, document) {
        var _a;
        const ref = refInfo.reference;
        // The reference may already have been resolved lazily by accessing its `ref` property.
        if (ref._ref === undefined) {
            ref._ref = ref_resolving;
            try {
                const description = this.getCandidate(refInfo);
                if ((0,syntax_tree/* isLinkingError */.et)(description)) {
                    ref._ref = description;
                }
                else {
                    ref._nodeDescription = description;
                    if (this.langiumDocuments().hasDocument(description.documentUri)) {
                        // The target document is already loaded
                        const linkedNode = this.loadAstNode(description);
                        ref._ref = linkedNode !== null && linkedNode !== void 0 ? linkedNode : this.createLinkingError(refInfo, description);
                    }
                    else {
                        // Try to load the target AST node later using the already provided description
                        ref._ref = undefined;
                    }
                }
            }
            catch (err) {
                console.error(`An error occurred while resolving reference to '${ref.$refText}':`, err);
                const errorMessage = (_a = err.message) !== null && _a !== void 0 ? _a : String(err);
                ref._ref = Object.assign(Object.assign({}, refInfo), { message: `An error occurred while resolving reference to '${ref.$refText}': ${errorMessage}` });
            }
            // Add the reference to the document's array of references
            // Only add if the reference has been not been resolved earlier
            // Otherwise we end up with duplicates
            // See also implementation of `buildReference`
            document.references.push(ref);
        }
    }
    unlink(document) {
        for (const ref of document.references) {
            delete ref._ref;
            delete ref._nodeDescription;
        }
        document.references = [];
    }
    getCandidate(refInfo) {
        const scope = this.scopeProvider.getScope(refInfo);
        const description = scope.getElement(refInfo.reference.$refText);
        return description !== null && description !== void 0 ? description : this.createLinkingError(refInfo);
    }
    buildReference(node, property, refNode, refText) {
        // See behavior description in doc of Linker, update that on changes in here.
        // eslint-disable-next-line @typescript-eslint/no-this-alias
        const linker = this;
        const reference = {
            $refNode: refNode,
            $refText: refText,
            get ref() {
                var _a;
                if ((0,syntax_tree/* isAstNode */.xA)(this._ref)) {
                    // Most frequent case: the target is already resolved.
                    return this._ref;
                }
                else if ((0,syntax_tree/* isAstNodeDescription */.SI)(this._nodeDescription)) {
                    // A candidate has been found before, but it is not loaded yet.
                    const linkedNode = linker.loadAstNode(this._nodeDescription);
                    this._ref = linkedNode !== null && linkedNode !== void 0 ? linkedNode : linker.createLinkingError({ reference, container: node, property }, this._nodeDescription);
                }
                else if (this._ref === undefined) {
                    // The reference has not been linked yet, so do that now.
                    this._ref = ref_resolving;
                    const document = (0,ast_utils/* findRootNode */.E$)(node).$document;
                    const refData = linker.getLinkedNode({ reference, container: node, property });
                    if (refData.error && document && document.state < DocumentState.ComputedScopes) {
                        // Document scope is not ready, don't set `this._ref` so linker can retry later.
                        return this._ref = undefined;
                    }
                    this._ref = (_a = refData.node) !== null && _a !== void 0 ? _a : refData.error;
                    this._nodeDescription = refData.descr;
                    document === null || document === void 0 ? void 0 : document.references.push(this);
                }
                else if (this._ref === ref_resolving) {
                    throw new Error(`Cyclic reference resolution detected: ${linker.astNodeLocator.getAstNodePath(node)}/${property} (symbol '${refText}')`);
                }
                return (0,syntax_tree/* isAstNode */.xA)(this._ref) ? this._ref : undefined;
            },
            get $nodeDescription() {
                return this._nodeDescription;
            },
            get error() {
                return (0,syntax_tree/* isLinkingError */.et)(this._ref) ? this._ref : undefined;
            }
        };
        return reference;
    }
    getLinkedNode(refInfo) {
        var _a;
        try {
            const description = this.getCandidate(refInfo);
            if ((0,syntax_tree/* isLinkingError */.et)(description)) {
                return { error: description };
            }
            const linkedNode = this.loadAstNode(description);
            if (linkedNode) {
                return { node: linkedNode, descr: description };
            }
            else {
                return {
                    descr: description,
                    error: this.createLinkingError(refInfo, description)
                };
            }
        }
        catch (err) {
            console.error(`An error occurred while resolving reference to '${refInfo.reference.$refText}':`, err);
            const errorMessage = (_a = err.message) !== null && _a !== void 0 ? _a : String(err);
            return {
                error: Object.assign(Object.assign({}, refInfo), { message: `An error occurred while resolving reference to '${refInfo.reference.$refText}': ${errorMessage}` })
            };
        }
    }
    loadAstNode(nodeDescription) {
        if (nodeDescription.node) {
            return nodeDescription.node;
        }
        const doc = this.langiumDocuments().getDocument(nodeDescription.documentUri);
        if (!doc) {
            return undefined;
        }
        return this.astNodeLocator.getAstNode(doc.parseResult.value, nodeDescription.path);
    }
    createLinkingError(refInfo, targetDescription) {
        // Check whether the document is sufficiently processed by the DocumentBuilder. If not, this is a hint for a bug
        // in the language implementation.
        const document = (0,ast_utils/* findRootNode */.E$)(refInfo.container).$document;
        if (document && document.state < DocumentState.ComputedScopes) {
            console.warn(`Attempted reference resolution before document reached ComputedScopes state (${document.uri}).`);
        }
        const referenceType = this.reflection.getReferenceType(refInfo);
        return Object.assign(Object.assign({}, refInfo), { message: `Could not resolve reference to ${referenceType} named '${refInfo.reference.$refText}'.`, targetDescription });
    }
}
//# sourceMappingURL=linker.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/name-provider.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

function isNamed(node) {
    return typeof node.name === 'string';
}
class DefaultNameProvider {
    getName(node) {
        if (isNamed(node)) {
            return node.name;
        }
        return undefined;
    }
    getNameNode(node) {
        return (0,grammar_utils/* findNodeForProperty */.vb)(node.$cstNode, 'name');
    }
}
//# sourceMappingURL=name-provider.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/uri-utils.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


var UriUtils;
(function (UriUtils) {
    UriUtils.basename = esm/* Utils */.c.basename;
    UriUtils.dirname = esm/* Utils */.c.dirname;
    UriUtils.extname = esm/* Utils */.c.extname;
    UriUtils.joinPath = esm/* Utils */.c.joinPath;
    UriUtils.resolvePath = esm/* Utils */.c.resolvePath;
    function equals(a, b) {
        return (a === null || a === void 0 ? void 0 : a.toString()) === (b === null || b === void 0 ? void 0 : b.toString());
    }
    UriUtils.equals = equals;
    function relative(from, to) {
        const fromPath = typeof from === 'string' ? from : from.path;
        const toPath = typeof to === 'string' ? to : to.path;
        const fromParts = fromPath.split('/').filter(e => e.length > 0);
        const toParts = toPath.split('/').filter(e => e.length > 0);
        let i = 0;
        for (; i < fromParts.length; i++) {
            if (fromParts[i] !== toParts[i]) {
                break;
            }
        }
        const backPart = '../'.repeat(fromParts.length - i);
        const toPart = toParts.slice(i).join('/');
        return backPart + toPart;
    }
    UriUtils.relative = relative;
    function normalize(uri) {
        return esm/* URI */.o.parse(uri.toString()).toString();
    }
    UriUtils.normalize = normalize;
})(UriUtils || (UriUtils = {}));
//# sourceMappingURL=uri-utils.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/references.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






class DefaultReferences {
    constructor(services) {
        this.nameProvider = services.references.NameProvider;
        this.index = services.shared.workspace.IndexManager;
        this.nodeLocator = services.workspace.AstNodeLocator;
    }
    findDeclaration(sourceCstNode) {
        if (sourceCstNode) {
            const assignment = (0,grammar_utils/* findAssignment */.h7)(sourceCstNode);
            const nodeElem = sourceCstNode.astNode;
            if (assignment && nodeElem) {
                const reference = nodeElem[assignment.feature];
                if ((0,syntax_tree/* isReference */.Yk)(reference)) {
                    return reference.ref;
                }
                else if (Array.isArray(reference)) {
                    for (const ref of reference) {
                        if ((0,syntax_tree/* isReference */.Yk)(ref) && ref.$refNode
                            && ref.$refNode.offset <= sourceCstNode.offset
                            && ref.$refNode.end >= sourceCstNode.end) {
                            return ref.ref;
                        }
                    }
                }
            }
            if (nodeElem) {
                const nameNode = this.nameProvider.getNameNode(nodeElem);
                // Only return the targeted node in case the targeted cst node is the name node or part of it
                if (nameNode && (nameNode === sourceCstNode || (0,cst_utils/* isChildNode */.OB)(sourceCstNode, nameNode))) {
                    return nodeElem;
                }
            }
        }
        return undefined;
    }
    findDeclarationNode(sourceCstNode) {
        const astNode = this.findDeclaration(sourceCstNode);
        if (astNode === null || astNode === void 0 ? void 0 : astNode.$cstNode) {
            const targetNode = this.nameProvider.getNameNode(astNode);
            return targetNode !== null && targetNode !== void 0 ? targetNode : astNode.$cstNode;
        }
        return undefined;
    }
    findReferences(targetNode, options) {
        const refs = [];
        if (options.includeDeclaration) {
            const ref = this.getReferenceToSelf(targetNode);
            if (ref) {
                refs.push(ref);
            }
        }
        let indexReferences = this.index.findAllReferences(targetNode, this.nodeLocator.getAstNodePath(targetNode));
        if (options.documentUri) {
            indexReferences = indexReferences.filter(ref => UriUtils.equals(ref.sourceUri, options.documentUri));
        }
        refs.push(...indexReferences);
        return (0,stream/* stream */.Vw)(refs);
    }
    getReferenceToSelf(targetNode) {
        const nameNode = this.nameProvider.getNameNode(targetNode);
        if (nameNode) {
            const doc = (0,ast_utils/* getDocument */.Me)(targetNode);
            const path = this.nodeLocator.getAstNodePath(targetNode);
            return {
                sourceUri: doc.uri,
                sourcePath: path,
                targetUri: doc.uri,
                targetPath: path,
                segment: (0,cst_utils/* toDocumentSegment */.yn)(nameNode),
                local: true
            };
        }
        return undefined;
    }
}
//# sourceMappingURL=references.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/collections.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

/**
 * A multimap is a variation of a Map that has potentially multiple values for every key.
 */
class MultiMap {
    constructor(elements) {
        this.map = new Map();
        if (elements) {
            for (const [key, value] of elements) {
                this.add(key, value);
            }
        }
    }
    /**
     * The total number of values in the multimap.
     */
    get size() {
        return stream/* Reduction */.IH.sum((0,stream/* stream */.Vw)(this.map.values()).map(a => a.length));
    }
    /**
     * Clear all entries in the multimap.
     */
    clear() {
        this.map.clear();
    }
    /**
     * Operates differently depending on whether a `value` is given:
     *  * With a value, this method deletes the specific key / value pair from the multimap.
     *  * Without a value, all values associated with the given key are deleted.
     *
     * @returns `true` if a value existed and has been removed, or `false` if the specified
     *     key / value does not exist.
     */
    delete(key, value) {
        if (value === undefined) {
            return this.map.delete(key);
        }
        else {
            const values = this.map.get(key);
            if (values) {
                const index = values.indexOf(value);
                if (index >= 0) {
                    if (values.length === 1) {
                        this.map.delete(key);
                    }
                    else {
                        values.splice(index, 1);
                    }
                    return true;
                }
            }
            return false;
        }
    }
    /**
     * Returns an array of all values associated with the given key. If no value exists,
     * an empty array is returned.
     *
     * _Note:_ The returned array is assumed not to be modified. Use the `set` method to add a
     * value and `delete` to remove a value from the multimap.
     */
    get(key) {
        var _a;
        return (_a = this.map.get(key)) !== null && _a !== void 0 ? _a : [];
    }
    /**
     * Operates differently depending on whether a `value` is given:
     *  * With a value, this method returns `true` if the specific key / value pair is present in the multimap.
     *  * Without a value, this method returns `true` if the given key is present in the multimap.
     */
    has(key, value) {
        if (value === undefined) {
            return this.map.has(key);
        }
        else {
            const values = this.map.get(key);
            if (values) {
                return values.indexOf(value) >= 0;
            }
            return false;
        }
    }
    /**
     * Add the given key / value pair to the multimap.
     */
    add(key, value) {
        if (this.map.has(key)) {
            this.map.get(key).push(value);
        }
        else {
            this.map.set(key, [value]);
        }
        return this;
    }
    /**
     * Add the given set of key / value pairs to the multimap.
     */
    addAll(key, values) {
        if (this.map.has(key)) {
            this.map.get(key).push(...values);
        }
        else {
            this.map.set(key, Array.from(values));
        }
        return this;
    }
    /**
     * Invokes the given callback function for every key / value pair in the multimap.
     */
    forEach(callbackfn) {
        this.map.forEach((array, key) => array.forEach(value => callbackfn(value, key, this)));
    }
    /**
     * Returns an iterator of key, value pairs for every entry in the map.
     */
    [Symbol.iterator]() {
        return this.entries().iterator();
    }
    /**
     * Returns a stream of key, value pairs for every entry in the map.
     */
    entries() {
        return (0,stream/* stream */.Vw)(this.map.entries())
            .flatMap(([key, array]) => array.map(value => [key, value]));
    }
    /**
     * Returns a stream of keys in the map.
     */
    keys() {
        return (0,stream/* stream */.Vw)(this.map.keys());
    }
    /**
     * Returns a stream of values in the map.
     */
    values() {
        return (0,stream/* stream */.Vw)(this.map.values()).flat();
    }
    /**
     * Returns a stream of key, value set pairs for every key in the map.
     */
    entriesGroupedByKey() {
        return (0,stream/* stream */.Vw)(this.map.entries());
    }
}
class BiMap {
    get size() {
        return this.map.size;
    }
    constructor(elements) {
        this.map = new Map();
        this.inverse = new Map();
        if (elements) {
            for (const [key, value] of elements) {
                this.set(key, value);
            }
        }
    }
    clear() {
        this.map.clear();
        this.inverse.clear();
    }
    set(key, value) {
        this.map.set(key, value);
        this.inverse.set(value, key);
        return this;
    }
    get(key) {
        return this.map.get(key);
    }
    getKey(value) {
        return this.inverse.get(value);
    }
    delete(key) {
        const value = this.map.get(key);
        if (value !== undefined) {
            this.map.delete(key);
            this.inverse.delete(value);
            return true;
        }
        return false;
    }
}
//# sourceMappingURL=collections.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/scope-computation.js
/******************************************************************************
 * Copyright 2021-2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/




/**
 * The default scope computation creates and collectes descriptions of the AST nodes to be exported into the
 * _global_ scope from the given document. By default those are the document's root AST node and its directly
 * contained child nodes.
 *
 * Besides, it gathers all AST nodes that have a name (according to the `NameProvider` service) and includes them
 * in the local scope of their particular container nodes. As a result, for every cross-reference in the AST,
 * target elements from the same level (siblings) and further up towards the root (parents and siblings of parents)
 * are visible. Elements being nested inside lower levels (children, children of siblings and parents' siblings)
 * are _invisible_ by default, but that can be changed by customizing this service.
 */
class DefaultScopeComputation {
    constructor(services) {
        this.nameProvider = services.references.NameProvider;
        this.descriptions = services.workspace.AstNodeDescriptionProvider;
    }
    async computeExports(document, cancelToken = cancellation.CancellationToken.None) {
        return this.computeExportsForNode(document.parseResult.value, document, undefined, cancelToken);
    }
    /**
     * Creates {@link AstNodeDescription AstNodeDescriptions} for the given {@link AstNode parentNode} and its children.
     * The list of children to be considered is determined by the function parameter {@link children}.
     * By default only the direct children of {@link parentNode} are visited, nested nodes are not exported.
     *
     * @param parentNode AST node to be exported, i.e., of which an {@link AstNodeDescription} shall be added to the returned list.
     * @param document The document containing the AST node to be exported.
     * @param children A function called with {@link parentNode} as single argument and returning an {@link Iterable} supplying the children to be visited, which must be directly or transitively contained in {@link parentNode}.
     * @param cancelToken Indicates when to cancel the current operation.
     * @throws `OperationCancelled` if a user action occurs during execution.
     * @returns A list of {@link AstNodeDescription AstNodeDescriptions} to be published to index.
     */
    async computeExportsForNode(parentNode, document, children = ast_utils/* streamContents */.sx, cancelToken = cancellation.CancellationToken.None) {
        const exports = [];
        this.exportNode(parentNode, exports, document);
        for (const node of children(parentNode)) {
            await interruptAndCheck(cancelToken);
            this.exportNode(node, exports, document);
        }
        return exports;
    }
    /**
     * Add a single node to the list of exports if it has a name. Override this method to change how
     * symbols are exported, e.g. by modifying their exported name.
     */
    exportNode(node, exports, document) {
        const name = this.nameProvider.getName(node);
        if (name) {
            exports.push(this.descriptions.createDescription(node, name, document));
        }
    }
    async computeLocalScopes(document, cancelToken = cancellation.CancellationToken.None) {
        const rootNode = document.parseResult.value;
        const scopes = new MultiMap();
        // Here we navigate the full AST - local scopes shall be available in the whole document
        for (const node of (0,ast_utils/* streamAllContents */.VY)(rootNode)) {
            await interruptAndCheck(cancelToken);
            this.processNode(node, document, scopes);
        }
        return scopes;
    }
    /**
     * Process a single node during scopes computation. The default implementation makes the node visible
     * in the subtree of its container (if the node has a name). Override this method to change this,
     * e.g. by increasing the visibility to a higher level in the AST.
     */
    processNode(node, document, scopes) {
        const container = node.$container;
        if (container) {
            const name = this.nameProvider.getName(node);
            if (name) {
                scopes.add(container, this.descriptions.createDescription(node, name, document));
            }
        }
    }
}
//# sourceMappingURL=scope-computation.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/scope.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

/**
 * The default scope implementation is based on a `Stream`. It has an optional _outer scope_ describing
 * the next level of elements, which are queried when a target element is not found in the stream provided
 * to this scope.
 */
class StreamScope {
    constructor(elements, outerScope, options) {
        var _a;
        this.elements = elements;
        this.outerScope = outerScope;
        this.caseInsensitive = (_a = options === null || options === void 0 ? void 0 : options.caseInsensitive) !== null && _a !== void 0 ? _a : false;
    }
    getAllElements() {
        if (this.outerScope) {
            return this.elements.concat(this.outerScope.getAllElements());
        }
        else {
            return this.elements;
        }
    }
    getElement(name) {
        const local = this.caseInsensitive
            ? this.elements.find(e => e.name.toLowerCase() === name.toLowerCase())
            : this.elements.find(e => e.name === name);
        if (local) {
            return local;
        }
        if (this.outerScope) {
            return this.outerScope.getElement(name);
        }
        return undefined;
    }
}
class MapScope {
    constructor(elements, outerScope, options) {
        var _a;
        this.elements = new Map();
        this.caseInsensitive = (_a = options === null || options === void 0 ? void 0 : options.caseInsensitive) !== null && _a !== void 0 ? _a : false;
        for (const element of elements) {
            const name = this.caseInsensitive
                ? element.name.toLowerCase()
                : element.name;
            this.elements.set(name, element);
        }
        this.outerScope = outerScope;
    }
    getElement(name) {
        const localName = this.caseInsensitive ? name.toLowerCase() : name;
        const local = this.elements.get(localName);
        if (local) {
            return local;
        }
        if (this.outerScope) {
            return this.outerScope.getElement(name);
        }
        return undefined;
    }
    getAllElements() {
        let elementStream = (0,stream/* stream */.Vw)(this.elements.values());
        if (this.outerScope) {
            elementStream = elementStream.concat(this.outerScope.getAllElements());
        }
        return elementStream;
    }
}
const EMPTY_SCOPE = {
    getElement() {
        return undefined;
    },
    getAllElements() {
        return stream/* EMPTY_STREAM */.Cl;
    }
};
//# sourceMappingURL=scope.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/caching.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
class DisposableCache {
    constructor() {
        this.toDispose = [];
        this.isDisposed = false;
    }
    onDispose(disposable) {
        this.toDispose.push(disposable);
    }
    dispose() {
        this.throwIfDisposed();
        this.clear();
        this.isDisposed = true;
        this.toDispose.forEach(disposable => disposable.dispose());
    }
    throwIfDisposed() {
        if (this.isDisposed) {
            throw new Error('This cache has already been disposed');
        }
    }
}
class SimpleCache extends DisposableCache {
    constructor() {
        super(...arguments);
        this.cache = new Map();
    }
    has(key) {
        this.throwIfDisposed();
        return this.cache.has(key);
    }
    set(key, value) {
        this.throwIfDisposed();
        this.cache.set(key, value);
    }
    get(key, provider) {
        this.throwIfDisposed();
        if (this.cache.has(key)) {
            return this.cache.get(key);
        }
        else if (provider) {
            const value = provider();
            this.cache.set(key, value);
            return value;
        }
        else {
            return undefined;
        }
    }
    delete(key) {
        this.throwIfDisposed();
        return this.cache.delete(key);
    }
    clear() {
        this.throwIfDisposed();
        this.cache.clear();
    }
}
class ContextCache extends DisposableCache {
    constructor(converter) {
        super();
        this.cache = new Map();
        this.converter = converter !== null && converter !== void 0 ? converter : (value => value);
    }
    has(contextKey, key) {
        this.throwIfDisposed();
        return this.cacheForContext(contextKey).has(key);
    }
    set(contextKey, key, value) {
        this.throwIfDisposed();
        this.cacheForContext(contextKey).set(key, value);
    }
    get(contextKey, key, provider) {
        this.throwIfDisposed();
        const contextCache = this.cacheForContext(contextKey);
        if (contextCache.has(key)) {
            return contextCache.get(key);
        }
        else if (provider) {
            const value = provider();
            contextCache.set(key, value);
            return value;
        }
        else {
            return undefined;
        }
    }
    delete(contextKey, key) {
        this.throwIfDisposed();
        return this.cacheForContext(contextKey).delete(key);
    }
    clear(contextKey) {
        this.throwIfDisposed();
        if (contextKey) {
            const mapKey = this.converter(contextKey);
            this.cache.delete(mapKey);
        }
        else {
            this.cache.clear();
        }
    }
    cacheForContext(contextKey) {
        const mapKey = this.converter(contextKey);
        let documentCache = this.cache.get(mapKey);
        if (!documentCache) {
            documentCache = new Map();
            this.cache.set(mapKey, documentCache);
        }
        return documentCache;
    }
}
/**
 * Every key/value pair in this cache is scoped to a document.
 * If this document is changed or deleted, all associated key/value pairs are deleted.
 */
class DocumentCache extends (/* unused pure expression or super */ null && (ContextCache)) {
    /**
     * Creates a new document cache.
     *
     * @param sharedServices Service container instance to hook into document lifecycle events.
     * @param state Optional document state on which the cache should evict.
     * If not provided, the cache will evict on `DocumentBuilder#onUpdate`.
     * *Deleted* documents are considered in both cases.
     *
     * Providing a state here will use `DocumentBuilder#onDocumentPhase` instead,
     * which triggers on all documents that have been affected by this change, assuming that the
     * state is `DocumentState.Linked` or a later state.
     */
    constructor(sharedServices, state) {
        super(uri => uri.toString());
        if (state) {
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onDocumentPhase(state, document => {
                this.clear(document.uri.toString());
            }));
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onUpdate((_changed, deleted) => {
                for (const uri of deleted) { // react only on deleted documents
                    this.clear(uri);
                }
            }));
        }
        else {
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onUpdate((changed, deleted) => {
                const allUris = changed.concat(deleted); // react on both changed and deleted documents
                for (const uri of allUris) {
                    this.clear(uri);
                }
            }));
        }
    }
}
/**
 * Every key/value pair in this cache is scoped to the whole workspace.
 * If any document in the workspace is added, changed or deleted, the whole cache is evicted.
 */
class WorkspaceCache extends SimpleCache {
    /**
     * Creates a new workspace cache.
     *
     * @param sharedServices Service container instance to hook into document lifecycle events.
     * @param state Optional document state on which the cache should evict.
     * If not provided, the cache will evict on `DocumentBuilder#onUpdate`.
     * *Deleted* documents are considered in both cases.
     */
    constructor(sharedServices, state) {
        super();
        if (state) {
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onBuildPhase(state, () => {
                this.clear();
            }));
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onUpdate((_changed, deleted) => {
                if (deleted.length > 0) { // react only on deleted documents
                    this.clear();
                }
            }));
        }
        else {
            this.toDispose.push(sharedServices.workspace.DocumentBuilder.onUpdate(() => {
                this.clear();
            }));
        }
    }
}
//# sourceMappingURL=caching.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/references/scope-provider.js
/******************************************************************************
 * Copyright 2021-2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/




class DefaultScopeProvider {
    constructor(services) {
        this.reflection = services.shared.AstReflection;
        this.nameProvider = services.references.NameProvider;
        this.descriptions = services.workspace.AstNodeDescriptionProvider;
        this.indexManager = services.shared.workspace.IndexManager;
        this.globalScopeCache = new WorkspaceCache(services.shared);
    }
    getScope(context) {
        const scopes = [];
        const referenceType = this.reflection.getReferenceType(context);
        const precomputed = (0,ast_utils/* getDocument */.Me)(context.container).precomputedScopes;
        if (precomputed) {
            let currentNode = context.container;
            do {
                const allDescriptions = precomputed.get(currentNode);
                if (allDescriptions.length > 0) {
                    scopes.push((0,stream/* stream */.Vw)(allDescriptions).filter(desc => this.reflection.isSubtype(desc.type, referenceType)));
                }
                currentNode = currentNode.$container;
            } while (currentNode);
        }
        let result = this.getGlobalScope(referenceType, context);
        for (let i = scopes.length - 1; i >= 0; i--) {
            result = this.createScope(scopes[i], result);
        }
        return result;
    }
    /**
     * Create a scope for the given collection of AST node descriptions.
     */
    createScope(elements, outerScope, options) {
        return new StreamScope((0,stream/* stream */.Vw)(elements), outerScope, options);
    }
    /**
     * Create a scope for the given collection of AST nodes, which need to be transformed into respective
     * descriptions first. This is done using the `NameProvider` and `AstNodeDescriptionProvider` services.
     */
    createScopeForNodes(elements, outerScope, options) {
        const s = (0,stream/* stream */.Vw)(elements).map(e => {
            const name = this.nameProvider.getName(e);
            if (name) {
                return this.descriptions.createDescription(e, name);
            }
            return undefined;
        }).nonNullable();
        return new StreamScope(s, outerScope, options);
    }
    /**
     * Create a global scope filtered for the given reference type.
     */
    getGlobalScope(referenceType, _context) {
        return this.globalScopeCache.get(referenceType, () => new MapScope(this.indexManager.allElements(referenceType)));
    }
}
//# sourceMappingURL=scope-provider.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/serializer/json-serializer.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/




function isAstNodeWithComment(node) {
    return typeof node.$comment === 'string';
}
function isIntermediateReference(obj) {
    return typeof obj === 'object' && !!obj && ('$ref' in obj || '$error' in obj);
}
class DefaultJsonSerializer {
    constructor(services) {
        /** The set of AstNode properties to be ignored by the serializer. */
        this.ignoreProperties = new Set(['$container', '$containerProperty', '$containerIndex', '$document', '$cstNode']);
        this.langiumDocuments = services.shared.workspace.LangiumDocuments;
        this.astNodeLocator = services.workspace.AstNodeLocator;
        this.nameProvider = services.references.NameProvider;
        this.commentProvider = services.documentation.CommentProvider;
    }
    serialize(node, options) {
        const serializeOptions = options !== null && options !== void 0 ? options : {};
        const specificReplacer = options === null || options === void 0 ? void 0 : options.replacer;
        const defaultReplacer = (key, value) => this.replacer(key, value, serializeOptions);
        const replacer = specificReplacer ? (key, value) => specificReplacer(key, value, defaultReplacer) : defaultReplacer;
        try {
            this.currentDocument = (0,ast_utils/* getDocument */.Me)(node);
            return JSON.stringify(node, replacer, options === null || options === void 0 ? void 0 : options.space);
        }
        finally {
            this.currentDocument = undefined;
        }
    }
    deserialize(content, options) {
        const deserializeOptions = options !== null && options !== void 0 ? options : {};
        const root = JSON.parse(content);
        this.linkNode(root, root, deserializeOptions);
        return root;
    }
    replacer(key, value, { refText, sourceText, textRegions, comments, uriConverter }) {
        var _a, _b, _c, _d;
        if (this.ignoreProperties.has(key)) {
            return undefined;
        }
        else if ((0,syntax_tree/* isReference */.Yk)(value)) {
            const refValue = value.ref;
            const $refText = refText ? value.$refText : undefined;
            if (refValue) {
                const targetDocument = (0,ast_utils/* getDocument */.Me)(refValue);
                let targetUri = '';
                if (this.currentDocument && this.currentDocument !== targetDocument) {
                    if (uriConverter) {
                        targetUri = uriConverter(targetDocument.uri, value);
                    }
                    else {
                        targetUri = targetDocument.uri.toString();
                    }
                }
                const targetPath = this.astNodeLocator.getAstNodePath(refValue);
                return {
                    $ref: `${targetUri}#${targetPath}`,
                    $refText
                };
            }
            else {
                return {
                    $error: (_b = (_a = value.error) === null || _a === void 0 ? void 0 : _a.message) !== null && _b !== void 0 ? _b : 'Could not resolve reference',
                    $refText
                };
            }
        }
        else if ((0,syntax_tree/* isAstNode */.xA)(value)) {
            let astNode = undefined;
            if (textRegions) {
                astNode = this.addAstNodeRegionWithAssignmentsTo(Object.assign({}, value));
                if ((!key || value.$document) && (astNode === null || astNode === void 0 ? void 0 : astNode.$textRegion)) {
                    // The document URI is added to the root node of the resulting JSON tree
                    astNode.$textRegion.documentURI = (_c = this.currentDocument) === null || _c === void 0 ? void 0 : _c.uri.toString();
                }
            }
            if (sourceText && !key) {
                astNode !== null && astNode !== void 0 ? astNode : (astNode = Object.assign({}, value));
                astNode.$sourceText = (_d = value.$cstNode) === null || _d === void 0 ? void 0 : _d.text;
            }
            if (comments) {
                astNode !== null && astNode !== void 0 ? astNode : (astNode = Object.assign({}, value));
                const comment = this.commentProvider.getComment(value);
                if (comment) {
                    astNode.$comment = comment.replace(/\r/g, '');
                }
            }
            return astNode !== null && astNode !== void 0 ? astNode : value;
        }
        else {
            return value;
        }
    }
    addAstNodeRegionWithAssignmentsTo(node) {
        const createDocumentSegment = cstNode => ({
            offset: cstNode.offset,
            end: cstNode.end,
            length: cstNode.length,
            range: cstNode.range,
        });
        if (node.$cstNode) {
            const textRegion = node.$textRegion = createDocumentSegment(node.$cstNode);
            const assignments = textRegion.assignments = {};
            Object.keys(node).filter(key => !key.startsWith('$')).forEach(key => {
                const propertyAssignments = (0,grammar_utils/* findNodesForProperty */.EL)(node.$cstNode, key).map(createDocumentSegment);
                if (propertyAssignments.length !== 0) {
                    assignments[key] = propertyAssignments;
                }
            });
            return node;
        }
        return undefined;
    }
    linkNode(node, root, options, container, containerProperty, containerIndex) {
        for (const [propertyName, item] of Object.entries(node)) {
            if (Array.isArray(item)) {
                for (let index = 0; index < item.length; index++) {
                    const element = item[index];
                    if (isIntermediateReference(element)) {
                        item[index] = this.reviveReference(node, propertyName, root, element, options);
                    }
                    else if ((0,syntax_tree/* isAstNode */.xA)(element)) {
                        this.linkNode(element, root, options, node, propertyName, index);
                    }
                }
            }
            else if (isIntermediateReference(item)) {
                node[propertyName] = this.reviveReference(node, propertyName, root, item, options);
            }
            else if ((0,syntax_tree/* isAstNode */.xA)(item)) {
                this.linkNode(item, root, options, node, propertyName);
            }
        }
        const mutable = node;
        mutable.$container = container;
        mutable.$containerProperty = containerProperty;
        mutable.$containerIndex = containerIndex;
    }
    reviveReference(container, property, root, reference, options) {
        let refText = reference.$refText;
        let error = reference.$error;
        if (reference.$ref) {
            const ref = this.getRefNode(root, reference.$ref, options.uriConverter);
            if ((0,syntax_tree/* isAstNode */.xA)(ref)) {
                if (!refText) {
                    refText = this.nameProvider.getName(ref);
                }
                return {
                    $refText: refText !== null && refText !== void 0 ? refText : '',
                    ref
                };
            }
            else {
                error = ref;
            }
        }
        if (error) {
            const ref = {
                $refText: refText !== null && refText !== void 0 ? refText : ''
            };
            ref.error = {
                container,
                property,
                message: error,
                reference: ref
            };
            return ref;
        }
        else {
            return undefined;
        }
    }
    getRefNode(root, uri, uriConverter) {
        try {
            const fragmentIndex = uri.indexOf('#');
            if (fragmentIndex === 0) {
                const node = this.astNodeLocator.getAstNode(root, uri.substring(1));
                if (!node) {
                    return 'Could not resolve path: ' + uri;
                }
                return node;
            }
            if (fragmentIndex < 0) {
                const documentUri = uriConverter ? uriConverter(uri) : esm/* URI */.o.parse(uri);
                const document = this.langiumDocuments.getDocument(documentUri);
                if (!document) {
                    return 'Could not find document for URI: ' + uri;
                }
                return document.parseResult.value;
            }
            const documentUri = uriConverter ? uriConverter(uri.substring(0, fragmentIndex)) : esm/* URI */.o.parse(uri.substring(0, fragmentIndex));
            const document = this.langiumDocuments.getDocument(documentUri);
            if (!document) {
                return 'Could not find document for URI: ' + uri;
            }
            if (fragmentIndex === uri.length - 1) {
                return document.parseResult.value;
            }
            const node = this.astNodeLocator.getAstNode(document.parseResult.value, uri.substring(fragmentIndex + 1));
            if (!node) {
                return 'Could not resolve URI: ' + uri;
            }
            return node;
        }
        catch (err) {
            return String(err);
        }
    }
}
//# sourceMappingURL=json-serializer.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/service-registry.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

/**
 * Generic registry for Langium services, but capable of being used with extending service sets as well (such as the lsp-complete LangiumCoreServices set)
 */
class DefaultServiceRegistry {
    /**
     * @deprecated Use the new `fileExtensionMap` (or `languageIdMap`) property instead.
     */
    get map() {
        return this.fileExtensionMap;
    }
    constructor(services) {
        this.languageIdMap = new Map();
        this.fileExtensionMap = new Map();
        this.textDocuments = services === null || services === void 0 ? void 0 : services.workspace.TextDocuments;
    }
    register(language) {
        const data = language.LanguageMetaData;
        for (const ext of data.fileExtensions) {
            if (this.fileExtensionMap.has(ext)) {
                console.warn(`The file extension ${ext} is used by multiple languages. It is now assigned to '${data.languageId}'.`);
            }
            this.fileExtensionMap.set(ext, language);
        }
        this.languageIdMap.set(data.languageId, language);
        if (this.languageIdMap.size === 1) {
            this.singleton = language;
        }
        else {
            this.singleton = undefined;
        }
    }
    getServices(uri) {
        var _a, _b;
        if (this.singleton !== undefined) {
            return this.singleton;
        }
        if (this.languageIdMap.size === 0) {
            throw new Error('The service registry is empty. Use `register` to register the services of a language.');
        }
        const languageId = (_b = (_a = this.textDocuments) === null || _a === void 0 ? void 0 : _a.get(uri)) === null || _b === void 0 ? void 0 : _b.languageId;
        if (languageId !== undefined) {
            const services = this.languageIdMap.get(languageId);
            if (services) {
                return services;
            }
        }
        const ext = UriUtils.extname(uri);
        const services = this.fileExtensionMap.get(ext);
        if (!services) {
            if (languageId) {
                throw new Error(`The service registry contains no services for the extension '${ext}' for language '${languageId}'.`);
            }
            else {
                throw new Error(`The service registry contains no services for the extension '${ext}'.`);
            }
        }
        return services;
    }
    hasServices(uri) {
        try {
            this.getServices(uri);
            return true;
        }
        catch (_a) {
            return false;
        }
    }
    get all() {
        return Array.from(this.languageIdMap.values());
    }
}
//# sourceMappingURL=service-registry.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/validation/validation-registry.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/




/**
 * Create DiagnosticData for a given diagnostic code. The result can be put into the `data` field of a DiagnosticInfo.
 */
function diagnosticData(code) {
    return { code };
}
var ValidationCategory;
(function (ValidationCategory) {
    ValidationCategory.all = ['fast', 'slow', 'built-in'];
})(ValidationCategory || (ValidationCategory = {}));
/**
 * Manages a set of `ValidationCheck`s to be applied when documents are validated.
 */
class ValidationRegistry {
    constructor(services) {
        this.entries = new MultiMap();
        this.entriesBefore = [];
        this.entriesAfter = [];
        this.reflection = services.shared.AstReflection;
    }
    /**
     * Register a set of validation checks. Each value in the record can be either a single validation check (i.e. a function)
     * or an array of validation checks.
     *
     * @param checksRecord Set of validation checks to register.
     * @param category Optional category for the validation checks (defaults to `'fast'`).
     * @param thisObj Optional object to be used as `this` when calling the validation check functions.
     */
    register(checksRecord, thisObj = this, category = 'fast') {
        if (category === 'built-in') {
            throw new Error("The 'built-in' category is reserved for lexer, parser, and linker errors.");
        }
        for (const [type, ch] of Object.entries(checksRecord)) {
            const callbacks = ch;
            if (Array.isArray(callbacks)) {
                for (const check of callbacks) {
                    const entry = {
                        check: this.wrapValidationException(check, thisObj),
                        category
                    };
                    this.addEntry(type, entry);
                }
            }
            else if (typeof callbacks === 'function') {
                const entry = {
                    check: this.wrapValidationException(callbacks, thisObj),
                    category
                };
                this.addEntry(type, entry);
            }
            else {
                (0,errors/* assertUnreachable */.U)(callbacks);
            }
        }
    }
    wrapValidationException(check, thisObj) {
        return async (node, accept, cancelToken) => {
            await this.handleException(() => check.call(thisObj, node, accept, cancelToken), 'An error occurred during validation', accept, node);
        };
    }
    async handleException(functionality, messageContext, accept, node) {
        try {
            await functionality();
        }
        catch (err) {
            if (isOperationCancelled(err)) {
                throw err;
            }
            console.error(`${messageContext}:`, err);
            if (err instanceof Error && err.stack) {
                console.error(err.stack);
            }
            const messageDetails = err instanceof Error ? err.message : String(err);
            accept('error', `${messageContext}: ${messageDetails}`, { node });
        }
    }
    addEntry(type, entry) {
        if (type === 'AstNode') {
            this.entries.add('AstNode', entry);
            return;
        }
        for (const subtype of this.reflection.getAllSubTypes(type)) {
            this.entries.add(subtype, entry);
        }
    }
    getChecks(type, categories) {
        let checks = (0,stream/* stream */.Vw)(this.entries.get(type))
            .concat(this.entries.get('AstNode'));
        if (categories) {
            checks = checks.filter(entry => categories.includes(entry.category));
        }
        return checks.map(entry => entry.check);
    }
    /**
     * Register logic which will be executed once before validating all the nodes of an AST/Langium document.
     * This helps to prepare or initialize some information which are required or reusable for the following checks on the AstNodes.
     *
     * As an example, for validating unique fully-qualified names of nodes in the AST,
     * here the map for mapping names to nodes could be established.
     * During the usual checks on the nodes, they are put into this map with their name.
     *
     * Note that this approach makes validations stateful, which is relevant e.g. when cancelling the validation.
     * Therefore it is recommended to clear stored information
     * _before_ validating an AST to validate each AST unaffected from other ASTs
     * AND _after_ validating the AST to free memory by information which are no longer used.
     *
     * @param checkBefore a set-up function which will be called once before actually validating an AST
     * @param thisObj Optional object to be used as `this` when calling the validation check functions.
     */
    registerBeforeDocument(checkBefore, thisObj = this) {
        this.entriesBefore.push(this.wrapPreparationException(checkBefore, 'An error occurred during set-up of the validation', thisObj));
    }
    /**
     * Register logic which will be executed once after validating all the nodes of an AST/Langium document.
     * This helps to finally evaluate information which are collected during the checks on the AstNodes.
     *
     * As an example, for validating unique fully-qualified names of nodes in the AST,
     * here the map with all the collected nodes and their names is checked
     * and validation hints are created for all nodes with the same name.
     *
     * Note that this approach makes validations stateful, which is relevant e.g. when cancelling the validation.
     * Therefore it is recommended to clear stored information
     * _before_ validating an AST to validate each AST unaffected from other ASTs
     * AND _after_ validating the AST to free memory by information which are no longer used.
     *
     * @param checkBefore a set-up function which will be called once before actually validating an AST
     * @param thisObj Optional object to be used as `this` when calling the validation check functions.
     */
    registerAfterDocument(checkAfter, thisObj = this) {
        this.entriesAfter.push(this.wrapPreparationException(checkAfter, 'An error occurred during tear-down of the validation', thisObj));
    }
    wrapPreparationException(check, messageContext, thisObj) {
        return async (rootNode, accept, categories, cancelToken) => {
            await this.handleException(() => check.call(thisObj, rootNode, accept, categories, cancelToken), messageContext, accept, rootNode);
        };
    }
    get checksBefore() {
        return this.entriesBefore;
    }
    get checksAfter() {
        return this.entriesAfter;
    }
}
//# sourceMappingURL=validation-registry.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/validation/document-validator.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






class DefaultDocumentValidator {
    constructor(services) {
        this.validationRegistry = services.validation.ValidationRegistry;
        this.metadata = services.LanguageMetaData;
    }
    async validateDocument(document, options = {}, cancelToken = cancellation.CancellationToken.None) {
        const parseResult = document.parseResult;
        const diagnostics = [];
        await interruptAndCheck(cancelToken);
        if (!options.categories || options.categories.includes('built-in')) {
            this.processLexingErrors(parseResult, diagnostics, options);
            if (options.stopAfterLexingErrors && diagnostics.some(d => { var _a; return ((_a = d.data) === null || _a === void 0 ? void 0 : _a.code) === DocumentValidator.LexingError; })) {
                return diagnostics;
            }
            this.processParsingErrors(parseResult, diagnostics, options);
            if (options.stopAfterParsingErrors && diagnostics.some(d => { var _a; return ((_a = d.data) === null || _a === void 0 ? void 0 : _a.code) === DocumentValidator.ParsingError; })) {
                return diagnostics;
            }
            this.processLinkingErrors(document, diagnostics, options);
            if (options.stopAfterLinkingErrors && diagnostics.some(d => { var _a; return ((_a = d.data) === null || _a === void 0 ? void 0 : _a.code) === DocumentValidator.LinkingError; })) {
                return diagnostics;
            }
        }
        // Process custom validations
        try {
            diagnostics.push(...await this.validateAst(parseResult.value, options, cancelToken));
        }
        catch (err) {
            if (isOperationCancelled(err)) {
                throw err;
            }
            console.error('An error occurred during validation:', err);
        }
        await interruptAndCheck(cancelToken);
        return diagnostics;
    }
    processLexingErrors(parseResult, diagnostics, _options) {
        var _a, _b, _c;
        const lexerDiagnostics = [...parseResult.lexerErrors, ...(_b = (_a = parseResult.lexerReport) === null || _a === void 0 ? void 0 : _a.diagnostics) !== null && _b !== void 0 ? _b : []];
        for (const lexerDiagnostic of lexerDiagnostics) {
            const severity = (_c = lexerDiagnostic.severity) !== null && _c !== void 0 ? _c : 'error';
            const diagnostic = {
                severity: toDiagnosticSeverity(severity),
                range: {
                    start: {
                        line: lexerDiagnostic.line - 1,
                        character: lexerDiagnostic.column - 1
                    },
                    end: {
                        line: lexerDiagnostic.line - 1,
                        character: lexerDiagnostic.column + lexerDiagnostic.length - 1
                    }
                },
                message: lexerDiagnostic.message,
                data: toDiagnosticData(severity),
                source: this.getSource()
            };
            diagnostics.push(diagnostic);
        }
    }
    processParsingErrors(parseResult, diagnostics, _options) {
        for (const parserError of parseResult.parserErrors) {
            let range = undefined;
            // We can run into the chevrotain error recovery here
            // The token contained in the parser error might be automatically inserted
            // In this case every position value will be `NaN`
            if (isNaN(parserError.token.startOffset)) {
                // Some special parser error types contain a `previousToken`
                // We can simply append our diagnostic to that token
                if ('previousToken' in parserError) {
                    const token = parserError.previousToken;
                    if (!isNaN(token.startOffset)) {
                        const position = { line: token.endLine - 1, character: token.endColumn };
                        range = { start: position, end: position };
                    }
                    else {
                        // No valid prev token. Might be empty document or containing only hidden tokens.
                        // Point to document start
                        const position = { line: 0, character: 0 };
                        range = { start: position, end: position };
                    }
                }
            }
            else {
                range = (0,cst_utils/* tokenToRange */.sp)(parserError.token);
            }
            if (range) {
                const diagnostic = {
                    severity: toDiagnosticSeverity('error'),
                    range,
                    message: parserError.message,
                    data: diagnosticData(DocumentValidator.ParsingError),
                    source: this.getSource()
                };
                diagnostics.push(diagnostic);
            }
        }
    }
    processLinkingErrors(document, diagnostics, _options) {
        for (const reference of document.references) {
            const linkingError = reference.error;
            if (linkingError) {
                const info = {
                    node: linkingError.container,
                    property: linkingError.property,
                    index: linkingError.index,
                    data: {
                        code: DocumentValidator.LinkingError,
                        containerType: linkingError.container.$type,
                        property: linkingError.property,
                        refText: linkingError.reference.$refText
                    }
                };
                diagnostics.push(this.toDiagnostic('error', linkingError.message, info));
            }
        }
    }
    async validateAst(rootNode, options, cancelToken = cancellation.CancellationToken.None) {
        const validationItems = [];
        const acceptor = (severity, message, info) => {
            validationItems.push(this.toDiagnostic(severity, message, info));
        };
        await this.validateAstBefore(rootNode, options, acceptor, cancelToken);
        await this.validateAstNodes(rootNode, options, acceptor, cancelToken);
        await this.validateAstAfter(rootNode, options, acceptor, cancelToken);
        return validationItems;
    }
    async validateAstBefore(rootNode, options, acceptor, cancelToken = cancellation.CancellationToken.None) {
        var _a;
        const checksBefore = this.validationRegistry.checksBefore;
        for (const checkBefore of checksBefore) {
            await interruptAndCheck(cancelToken);
            await checkBefore(rootNode, acceptor, (_a = options.categories) !== null && _a !== void 0 ? _a : [], cancelToken);
        }
    }
    async validateAstNodes(rootNode, options, acceptor, cancelToken = cancellation.CancellationToken.None) {
        await Promise.all((0,ast_utils/* streamAst */.Zc)(rootNode).map(async (node) => {
            await interruptAndCheck(cancelToken);
            const checks = this.validationRegistry.getChecks(node.$type, options.categories);
            for (const check of checks) {
                await check(node, acceptor, cancelToken);
            }
        }));
    }
    async validateAstAfter(rootNode, options, acceptor, cancelToken = cancellation.CancellationToken.None) {
        var _a;
        const checksAfter = this.validationRegistry.checksAfter;
        for (const checkAfter of checksAfter) {
            await interruptAndCheck(cancelToken);
            await checkAfter(rootNode, acceptor, (_a = options.categories) !== null && _a !== void 0 ? _a : [], cancelToken);
        }
    }
    toDiagnostic(severity, message, info) {
        return {
            message,
            range: getDiagnosticRange(info),
            severity: toDiagnosticSeverity(severity),
            code: info.code,
            codeDescription: info.codeDescription,
            tags: info.tags,
            relatedInformation: info.relatedInformation,
            data: info.data,
            source: this.getSource()
        };
    }
    getSource() {
        return this.metadata.languageId;
    }
}
function getDiagnosticRange(info) {
    if (info.range) {
        return info.range;
    }
    let cstNode;
    if (typeof info.property === 'string') {
        cstNode = (0,grammar_utils/* findNodeForProperty */.vb)(info.node.$cstNode, info.property, info.index);
    }
    else if (typeof info.keyword === 'string') {
        cstNode = (0,grammar_utils/* findNodeForKeyword */.lA)(info.node.$cstNode, info.keyword, info.index);
    }
    cstNode !== null && cstNode !== void 0 ? cstNode : (cstNode = info.node.$cstNode);
    if (!cstNode) {
        return {
            start: { line: 0, character: 0 },
            end: { line: 0, character: 0 }
        };
    }
    return cstNode.range;
}
/**
 * Transforms the diagnostic severity from the {@link LexingDiagnosticSeverity} format to LSP's `DiagnosticSeverity` format.
 *
 * @param severity The lexing diagnostic severity
 * @returns Diagnostic severity according to `vscode-languageserver-types/lib/esm/main.js#DiagnosticSeverity`
 */
function toDiagnosticSeverity(severity) {
    switch (severity) {
        case 'error':
            return 1;
        case 'warning':
            return 2;
        case 'info':
            return 3;
        case 'hint':
            return 4;
        default:
            throw new Error('Invalid diagnostic severity: ' + severity);
    }
}
function toDiagnosticData(severity) {
    switch (severity) {
        case 'error':
            return diagnosticData(DocumentValidator.LexingError);
        case 'warning':
            return diagnosticData(DocumentValidator.LexingWarning);
        case 'info':
            return diagnosticData(DocumentValidator.LexingInfo);
        case 'hint':
            return diagnosticData(DocumentValidator.LexingHint);
        default:
            throw new Error('Invalid diagnostic severity: ' + severity);
    }
}
var DocumentValidator;
(function (DocumentValidator) {
    DocumentValidator.LexingError = 'lexing-error';
    DocumentValidator.LexingWarning = 'lexing-warning';
    DocumentValidator.LexingInfo = 'lexing-info';
    DocumentValidator.LexingHint = 'lexing-hint';
    DocumentValidator.ParsingError = 'parsing-error';
    DocumentValidator.LinkingError = 'linking-error';
})(DocumentValidator || (DocumentValidator = {}));
//# sourceMappingURL=document-validator.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/ast-descriptions.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






class DefaultAstNodeDescriptionProvider {
    constructor(services) {
        this.astNodeLocator = services.workspace.AstNodeLocator;
        this.nameProvider = services.references.NameProvider;
    }
    createDescription(node, name, document) {
        const doc = document !== null && document !== void 0 ? document : (0,ast_utils/* getDocument */.Me)(node);
        name !== null && name !== void 0 ? name : (name = this.nameProvider.getName(node));
        const path = this.astNodeLocator.getAstNodePath(node);
        if (!name) {
            throw new Error(`Node at path ${path} has no name.`);
        }
        let nameNodeSegment;
        const nameSegmentGetter = () => { var _a; return nameNodeSegment !== null && nameNodeSegment !== void 0 ? nameNodeSegment : (nameNodeSegment = (0,cst_utils/* toDocumentSegment */.yn)((_a = this.nameProvider.getNameNode(node)) !== null && _a !== void 0 ? _a : node.$cstNode)); };
        return {
            node,
            name,
            get nameSegment() {
                return nameSegmentGetter();
            },
            selectionSegment: (0,cst_utils/* toDocumentSegment */.yn)(node.$cstNode),
            type: node.$type,
            documentUri: doc.uri,
            path
        };
    }
}
class DefaultReferenceDescriptionProvider {
    constructor(services) {
        this.nodeLocator = services.workspace.AstNodeLocator;
    }
    async createDescriptions(document, cancelToken = cancellation.CancellationToken.None) {
        const descr = [];
        const rootNode = document.parseResult.value;
        for (const astNode of (0,ast_utils/* streamAst */.Zc)(rootNode)) {
            await interruptAndCheck(cancelToken);
            (0,ast_utils/* streamReferences */.fy)(astNode).filter(refInfo => !(0,syntax_tree/* isLinkingError */.et)(refInfo)).forEach(refInfo => {
                // TODO: Consider logging a warning or throw an exception when DocumentState is < than Linked
                const description = this.createDescription(refInfo);
                if (description) {
                    descr.push(description);
                }
            });
        }
        return descr;
    }
    createDescription(refInfo) {
        const targetNodeDescr = refInfo.reference.$nodeDescription;
        const refCstNode = refInfo.reference.$refNode;
        if (!targetNodeDescr || !refCstNode) {
            return undefined;
        }
        const docUri = (0,ast_utils/* getDocument */.Me)(refInfo.container).uri;
        return {
            sourceUri: docUri,
            sourcePath: this.nodeLocator.getAstNodePath(refInfo.container),
            targetUri: targetNodeDescr.documentUri,
            targetPath: targetNodeDescr.path,
            segment: (0,cst_utils/* toDocumentSegment */.yn)(refCstNode),
            local: UriUtils.equals(targetNodeDescr.documentUri, docUri)
        };
    }
}
//# sourceMappingURL=ast-descriptions.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/ast-node-locator.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
class DefaultAstNodeLocator {
    constructor() {
        this.segmentSeparator = '/';
        this.indexSeparator = '@';
    }
    getAstNodePath(node) {
        if (node.$container) {
            const containerPath = this.getAstNodePath(node.$container);
            const newSegment = this.getPathSegment(node);
            const nodePath = containerPath + this.segmentSeparator + newSegment;
            return nodePath;
        }
        return '';
    }
    getPathSegment({ $containerProperty, $containerIndex }) {
        if (!$containerProperty) {
            throw new Error("Missing '$containerProperty' in AST node.");
        }
        if ($containerIndex !== undefined) {
            return $containerProperty + this.indexSeparator + $containerIndex;
        }
        return $containerProperty;
    }
    getAstNode(node, path) {
        const segments = path.split(this.segmentSeparator);
        return segments.reduce((previousValue, currentValue) => {
            if (!previousValue || currentValue.length === 0) {
                return previousValue;
            }
            const propertyIndex = currentValue.indexOf(this.indexSeparator);
            if (propertyIndex > 0) {
                const property = currentValue.substring(0, propertyIndex);
                const arrayIndex = parseInt(currentValue.substring(propertyIndex + 1));
                const array = previousValue[property];
                return array === null || array === void 0 ? void 0 : array[arrayIndex];
            }
            return previousValue[currentValue];
        }, node);
    }
}
//# sourceMappingURL=ast-node-locator.js.map
// EXTERNAL MODULE: ../node_modules/vscode-jsonrpc/lib/common/events.js
var events = __webpack_require__(345);
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/configuration.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


/**
 * Base configuration provider for building up other configuration providers
 */
class DefaultConfigurationProvider {
    constructor(services) {
        this._ready = new promise_utils_Deferred();
        this.settings = {};
        this.workspaceConfig = false;
        this.onConfigurationSectionUpdateEmitter = new events.Emitter();
        this.serviceRegistry = services.ServiceRegistry;
    }
    get ready() {
        return this._ready.promise;
    }
    initialize(params) {
        var _a, _b;
        this.workspaceConfig = (_b = (_a = params.capabilities.workspace) === null || _a === void 0 ? void 0 : _a.configuration) !== null && _b !== void 0 ? _b : false;
    }
    async initialized(params) {
        if (this.workspaceConfig) {
            if (params.register) {
                // params.register(...) is a function to be provided by the calling language server for the sake of
                //  decoupling this implementation from the concrete LSP implementations, specifically the LSP Connection
                const languages = this.serviceRegistry.all;
                params.register({
                    // Listen to configuration changes for all languages
                    section: languages.map(lang => this.toSectionName(lang.LanguageMetaData.languageId))
                });
            }
            if (params.fetchConfiguration) {
                // params.fetchConfiguration(...) is a function to be provided by the calling language server for the sake of
                //  decoupling this implementation from the concrete LSP implementations, specifically the LSP Connection
                const configToUpdate = this.serviceRegistry.all.map(lang => ({
                    // Fetch the configuration changes for all languages
                    section: this.toSectionName(lang.LanguageMetaData.languageId)
                }));
                // get workspace configurations (default scope URI)
                const configs = await params.fetchConfiguration(configToUpdate);
                configToUpdate.forEach((conf, idx) => {
                    this.updateSectionConfiguration(conf.section, configs[idx]);
                });
            }
        }
        this._ready.resolve();
    }
    /**
     *  Updates the cached configurations using the `change` notification parameters.
     *
     * @param change The parameters of a change configuration notification.
     * `settings` property of the change object could be expressed as `Record<string, Record<string, any>>`
     */
    updateConfiguration(change) {
        if (!change.settings) {
            return;
        }
        Object.keys(change.settings).forEach(section => {
            const configuration = change.settings[section];
            this.updateSectionConfiguration(section, configuration);
            this.onConfigurationSectionUpdateEmitter.fire({ section, configuration });
        });
    }
    updateSectionConfiguration(section, configuration) {
        this.settings[section] = configuration;
    }
    /**
    * Returns a configuration value stored for the given language.
    *
    * @param language The language id
    * @param configuration Configuration name
    */
    async getConfiguration(language, configuration) {
        await this.ready;
        const sectionName = this.toSectionName(language);
        if (this.settings[sectionName]) {
            return this.settings[sectionName][configuration];
        }
    }
    toSectionName(languageId) {
        return `${languageId}`;
    }
    get onConfigurationSectionUpdate() {
        return this.onConfigurationSectionUpdateEmitter.event;
    }
}
//# sourceMappingURL=configuration.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/utils/disposable.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
var Disposable;
(function (Disposable) {
    function create(callback) {
        return {
            dispose: async () => await callback()
        };
    }
    Disposable.create = create;
})(Disposable || (Disposable = {}));
//# sourceMappingURL=disposable.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/document-builder.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/







class DefaultDocumentBuilder {
    constructor(services) {
        this.updateBuildOptions = {
            // Default: run only the built-in validation checks and those in the _fast_ category (includes those without category)
            validation: {
                categories: ['built-in', 'fast']
            }
        };
        this.updateListeners = [];
        this.buildPhaseListeners = new MultiMap();
        this.documentPhaseListeners = new MultiMap();
        this.buildState = new Map();
        this.documentBuildWaiters = new Map();
        this.currentState = DocumentState.Changed;
        this.langiumDocuments = services.workspace.LangiumDocuments;
        this.langiumDocumentFactory = services.workspace.LangiumDocumentFactory;
        this.textDocuments = services.workspace.TextDocuments;
        this.indexManager = services.workspace.IndexManager;
        this.serviceRegistry = services.ServiceRegistry;
    }
    async build(documents, options = {}, cancelToken = cancellation.CancellationToken.None) {
        var _a, _b;
        for (const document of documents) {
            const key = document.uri.toString();
            if (document.state === DocumentState.Validated) {
                if (typeof options.validation === 'boolean' && options.validation) {
                    // Force re-running all validation checks
                    document.state = DocumentState.IndexedReferences;
                    document.diagnostics = undefined;
                    this.buildState.delete(key);
                }
                else if (typeof options.validation === 'object') {
                    const buildState = this.buildState.get(key);
                    const previousCategories = (_a = buildState === null || buildState === void 0 ? void 0 : buildState.result) === null || _a === void 0 ? void 0 : _a.validationChecks;
                    if (previousCategories) {
                        // Validation with explicit options was requested for a document that has already been partly validated.
                        // In this case, we need to merge the previous validation categories with the new ones.
                        const newCategories = (_b = options.validation.categories) !== null && _b !== void 0 ? _b : ValidationCategory.all;
                        const categories = newCategories.filter(c => !previousCategories.includes(c));
                        if (categories.length > 0) {
                            this.buildState.set(key, {
                                completed: false,
                                options: {
                                    validation: Object.assign(Object.assign({}, options.validation), { categories })
                                },
                                result: buildState.result
                            });
                            document.state = DocumentState.IndexedReferences;
                        }
                    }
                }
            }
            else {
                // Default: forget any previous build options
                this.buildState.delete(key);
            }
        }
        this.currentState = DocumentState.Changed;
        await this.emitUpdate(documents.map(e => e.uri), []);
        await this.buildDocuments(documents, options, cancelToken);
    }
    async update(changed, deleted, cancelToken = cancellation.CancellationToken.None) {
        this.currentState = DocumentState.Changed;
        // Remove all metadata of documents that are reported as deleted
        for (const deletedUri of deleted) {
            this.langiumDocuments.deleteDocument(deletedUri);
            this.buildState.delete(deletedUri.toString());
            this.indexManager.remove(deletedUri);
        }
        // Set the state of all changed documents to `Changed` so they are completely rebuilt
        for (const changedUri of changed) {
            const invalidated = this.langiumDocuments.invalidateDocument(changedUri);
            if (!invalidated) {
                // We create an unparsed, invalid document.
                // This will be parsed as soon as we reach the first document builder phase.
                // This allows to cancel the parsing process later in case we need it.
                const newDocument = this.langiumDocumentFactory.fromModel({ $type: 'INVALID' }, changedUri);
                newDocument.state = DocumentState.Changed;
                this.langiumDocuments.addDocument(newDocument);
            }
            this.buildState.delete(changedUri.toString());
        }
        // Set the state of all documents that should be relinked to `ComputedScopes` (if not already lower)
        const allChangedUris = (0,stream/* stream */.Vw)(changed).concat(deleted).map(uri => uri.toString()).toSet();
        this.langiumDocuments.all
            .filter(doc => !allChangedUris.has(doc.uri.toString()) && this.shouldRelink(doc, allChangedUris))
            .forEach(doc => {
            const linker = this.serviceRegistry.getServices(doc.uri).references.Linker;
            linker.unlink(doc);
            doc.state = Math.min(doc.state, DocumentState.ComputedScopes);
            doc.diagnostics = undefined;
        });
        // Notify listeners of the update
        await this.emitUpdate(changed, deleted);
        // Only allow interrupting the execution after all state changes are done
        await interruptAndCheck(cancelToken);
        // Collect and sort all documents that we should rebuild
        const rebuildDocuments = this.sortDocuments(this.langiumDocuments.all
            .filter(doc => {
            var _a;
            // This includes those that were reported as changed and those that we selected for relinking
            return doc.state < DocumentState.Linked
                // This includes those for which a previous build has been cancelled
                || !((_a = this.buildState.get(doc.uri.toString())) === null || _a === void 0 ? void 0 : _a.completed);
        })
            .toArray());
        await this.buildDocuments(rebuildDocuments, this.updateBuildOptions, cancelToken);
    }
    async emitUpdate(changed, deleted) {
        await Promise.all(this.updateListeners.map(listener => listener(changed, deleted)));
    }
    /**
     * Sort the given documents by priority. By default, documents with an open text document are prioritized.
     * This is useful to ensure that visible documents show their diagnostics before all other documents.
     *
     * This improves the responsiveness in large workspaces as users usually don't care about diagnostics
     * in files that are currently not opened in the editor.
     */
    sortDocuments(documents) {
        let left = 0;
        let right = documents.length - 1;
        while (left < right) {
            while (left < documents.length && this.hasTextDocument(documents[left])) {
                left++;
            }
            while (right >= 0 && !this.hasTextDocument(documents[right])) {
                right--;
            }
            if (left < right) {
                [documents[left], documents[right]] = [documents[right], documents[left]];
            }
        }
        return documents;
    }
    hasTextDocument(doc) {
        var _a;
        return Boolean((_a = this.textDocuments) === null || _a === void 0 ? void 0 : _a.get(doc.uri));
    }
    /**
     * Check whether the given document should be relinked after changes were found in the given URIs.
     */
    shouldRelink(document, changedUris) {
        // Relink documents with linking errors -- maybe those references can be resolved now
        if (document.references.some(ref => ref.error !== undefined)) {
            return true;
        }
        // Check whether the document is affected by any of the changed URIs
        return this.indexManager.isAffected(document, changedUris);
    }
    onUpdate(callback) {
        this.updateListeners.push(callback);
        return Disposable.create(() => {
            const index = this.updateListeners.indexOf(callback);
            if (index >= 0) {
                this.updateListeners.splice(index, 1);
            }
        });
    }
    /**
     * Build the given documents by stepping through all build phases. If a document's state indicates
     * that a certain build phase is already done, the phase is skipped for that document.
     *
     * @param documents The documents to build.
     * @param options the {@link BuildOptions} to use.
     * @param cancelToken A cancellation token that can be used to cancel the build.
     * @returns A promise that resolves when the build is done.
     */
    async buildDocuments(documents, options, cancelToken) {
        this.prepareBuild(documents, options);
        // 0. Parse content
        await this.runCancelable(documents, DocumentState.Parsed, cancelToken, doc => this.langiumDocumentFactory.update(doc, cancelToken));
        // 1. Index content
        await this.runCancelable(documents, DocumentState.IndexedContent, cancelToken, doc => this.indexManager.updateContent(doc, cancelToken));
        // 2. Compute scopes
        await this.runCancelable(documents, DocumentState.ComputedScopes, cancelToken, async (doc) => {
            const scopeComputation = this.serviceRegistry.getServices(doc.uri).references.ScopeComputation;
            doc.precomputedScopes = await scopeComputation.computeLocalScopes(doc, cancelToken);
        });
        // 3. Linking
        await this.runCancelable(documents, DocumentState.Linked, cancelToken, doc => {
            const linker = this.serviceRegistry.getServices(doc.uri).references.Linker;
            return linker.link(doc, cancelToken);
        });
        // 4. Index references
        await this.runCancelable(documents, DocumentState.IndexedReferences, cancelToken, doc => this.indexManager.updateReferences(doc, cancelToken));
        // 5. Validation
        const toBeValidated = documents.filter(doc => this.shouldValidate(doc));
        await this.runCancelable(toBeValidated, DocumentState.Validated, cancelToken, doc => this.validate(doc, cancelToken));
        // If we've made it to this point without being cancelled, we can mark the build state as completed.
        for (const doc of documents) {
            const state = this.buildState.get(doc.uri.toString());
            if (state) {
                state.completed = true;
            }
        }
    }
    /**
     * Runs prior to beginning the build process to update the {@link DocumentBuildState} for each document
     *
     * @param documents collection of documents to be built
     * @param options the {@link BuildOptions} to use
     */
    prepareBuild(documents, options) {
        for (const doc of documents) {
            const key = doc.uri.toString();
            const state = this.buildState.get(key);
            // If the document has no previous build state, we set it. If it has one, but it's already marked
            // as completed, we overwrite it. If the previous build was not completed, we keep its state
            // and continue where it was cancelled.
            if (!state || state.completed) {
                this.buildState.set(key, {
                    completed: false,
                    options,
                    result: state === null || state === void 0 ? void 0 : state.result
                });
            }
        }
    }
    /**
     * Runs a cancelable operation on a set of documents to bring them to a specified {@link DocumentState}.
     *
     * @param documents The array of documents to process.
     * @param targetState The target {@link DocumentState} to bring the documents to.
     * @param cancelToken A token that can be used to cancel the operation.
     * @param callback A function to be called for each document.
     * @returns A promise that resolves when all documents have been processed or the operation is canceled.
     * @throws Will throw `OperationCancelled` if the operation is canceled via a `CancellationToken`.
     */
    async runCancelable(documents, targetState, cancelToken, callback) {
        const filtered = documents.filter(doc => doc.state < targetState);
        for (const document of filtered) {
            await interruptAndCheck(cancelToken);
            await callback(document);
            document.state = targetState;
            await this.notifyDocumentPhase(document, targetState, cancelToken);
        }
        // Do not use `filtered` here, as that will miss documents that have previously reached the current target state
        // For example, this happens in case the cancellation triggers between the processing of two documents
        // Or files that were picked up during the workspace initialization
        const targetStateDocs = documents.filter(doc => doc.state === targetState);
        await this.notifyBuildPhase(targetStateDocs, targetState, cancelToken);
        this.currentState = targetState;
    }
    onBuildPhase(targetState, callback) {
        this.buildPhaseListeners.add(targetState, callback);
        return Disposable.create(() => {
            this.buildPhaseListeners.delete(targetState, callback);
        });
    }
    onDocumentPhase(targetState, callback) {
        this.documentPhaseListeners.add(targetState, callback);
        return Disposable.create(() => {
            this.documentPhaseListeners.delete(targetState, callback);
        });
    }
    waitUntil(state, uriOrToken, cancelToken) {
        let uri = undefined;
        if (uriOrToken && 'path' in uriOrToken) {
            uri = uriOrToken;
        }
        else {
            cancelToken = uriOrToken;
        }
        cancelToken !== null && cancelToken !== void 0 ? cancelToken : (cancelToken = cancellation.CancellationToken.None);
        if (uri) {
            const document = this.langiumDocuments.getDocument(uri);
            if (document && document.state > state) {
                return Promise.resolve(uri);
            }
        }
        if (this.currentState >= state) {
            return Promise.resolve(undefined);
        }
        else if (cancelToken.isCancellationRequested) {
            return Promise.reject(promise_utils_OperationCancelled);
        }
        return new Promise((resolve, reject) => {
            const buildDisposable = this.onBuildPhase(state, () => {
                buildDisposable.dispose();
                cancelDisposable.dispose();
                if (uri) {
                    const document = this.langiumDocuments.getDocument(uri);
                    resolve(document === null || document === void 0 ? void 0 : document.uri);
                }
                else {
                    resolve(undefined);
                }
            });
            const cancelDisposable = cancelToken.onCancellationRequested(() => {
                buildDisposable.dispose();
                cancelDisposable.dispose();
                reject(promise_utils_OperationCancelled);
            });
        });
    }
    async notifyDocumentPhase(document, state, cancelToken) {
        const listeners = this.documentPhaseListeners.get(state);
        const listenersCopy = listeners.slice();
        for (const listener of listenersCopy) {
            try {
                await listener(document, cancelToken);
            }
            catch (err) {
                // Ignore cancellation errors
                // We want to finish the listeners before throwing
                if (!isOperationCancelled(err)) {
                    throw err;
                }
            }
        }
    }
    async notifyBuildPhase(documents, state, cancelToken) {
        if (documents.length === 0) {
            // Don't notify when no document has been processed
            return;
        }
        const listeners = this.buildPhaseListeners.get(state);
        const listenersCopy = listeners.slice();
        for (const listener of listenersCopy) {
            await interruptAndCheck(cancelToken);
            await listener(documents, cancelToken);
        }
    }
    /**
     * Determine whether the given document should be validated during a build. The default
     * implementation checks the `validation` property of the build options. If it's set to `true`
     * or a `ValidationOptions` object, the document is included in the validation phase.
     */
    shouldValidate(document) {
        return Boolean(this.getBuildOptions(document).validation);
    }
    /**
     * Run validation checks on the given document and store the resulting diagnostics in the document.
     * If the document already contains diagnostics, the new ones are added to the list.
     */
    async validate(document, cancelToken) {
        var _a, _b;
        const validator = this.serviceRegistry.getServices(document.uri).validation.DocumentValidator;
        const validationSetting = this.getBuildOptions(document).validation;
        const options = typeof validationSetting === 'object' ? validationSetting : undefined;
        const diagnostics = await validator.validateDocument(document, options, cancelToken);
        if (document.diagnostics) {
            document.diagnostics.push(...diagnostics);
        }
        else {
            document.diagnostics = diagnostics;
        }
        // Store information about the executed validation in the build state
        const state = this.buildState.get(document.uri.toString());
        if (state) {
            (_a = state.result) !== null && _a !== void 0 ? _a : (state.result = {});
            const newCategories = (_b = options === null || options === void 0 ? void 0 : options.categories) !== null && _b !== void 0 ? _b : ValidationCategory.all;
            if (state.result.validationChecks) {
                state.result.validationChecks.push(...newCategories);
            }
            else {
                state.result.validationChecks = [...newCategories];
            }
        }
    }
    getBuildOptions(document) {
        var _a, _b;
        return (_b = (_a = this.buildState.get(document.uri.toString())) === null || _a === void 0 ? void 0 : _a.options) !== null && _b !== void 0 ? _b : {};
    }
}
//# sourceMappingURL=document-builder.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/index-manager.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/





class DefaultIndexManager {
    constructor(services) {
        /**
         * The symbol index stores all `AstNodeDescription` items exported by a document.
         * The key used in this map is the string representation of the specific document URI.
         */
        this.symbolIndex = new Map();
        /**
         * This is a cache for the `allElements()` method.
         * It caches the descriptions from `symbolIndex` grouped by types.
         */
        this.symbolByTypeIndex = new ContextCache();
        /**
         * This index keeps track of all `ReferenceDescription` items exported by a document.
         * This is used to compute which elements are affected by a document change
         * and for finding references to an AST node.
         */
        this.referenceIndex = new Map();
        this.documents = services.workspace.LangiumDocuments;
        this.serviceRegistry = services.ServiceRegistry;
        this.astReflection = services.AstReflection;
    }
    findAllReferences(targetNode, astNodePath) {
        const targetDocUri = (0,ast_utils/* getDocument */.Me)(targetNode).uri;
        const result = [];
        this.referenceIndex.forEach(docRefs => {
            docRefs.forEach(refDescr => {
                if (UriUtils.equals(refDescr.targetUri, targetDocUri) && refDescr.targetPath === astNodePath) {
                    result.push(refDescr);
                }
            });
        });
        return (0,stream/* stream */.Vw)(result);
    }
    allElements(nodeType, uris) {
        let documentUris = (0,stream/* stream */.Vw)(this.symbolIndex.keys());
        if (uris) {
            documentUris = documentUris.filter(uri => !uris || uris.has(uri));
        }
        return documentUris
            .map(uri => this.getFileDescriptions(uri, nodeType))
            .flat();
    }
    getFileDescriptions(uri, nodeType) {
        var _a;
        if (!nodeType) {
            return (_a = this.symbolIndex.get(uri)) !== null && _a !== void 0 ? _a : [];
        }
        const descriptions = this.symbolByTypeIndex.get(uri, nodeType, () => {
            var _a;
            const allFileDescriptions = (_a = this.symbolIndex.get(uri)) !== null && _a !== void 0 ? _a : [];
            return allFileDescriptions.filter(e => this.astReflection.isSubtype(e.type, nodeType));
        });
        return descriptions;
    }
    remove(uri) {
        const uriString = uri.toString();
        this.symbolIndex.delete(uriString);
        this.symbolByTypeIndex.clear(uriString);
        this.referenceIndex.delete(uriString);
    }
    async updateContent(document, cancelToken = cancellation.CancellationToken.None) {
        const services = this.serviceRegistry.getServices(document.uri);
        const exports = await services.references.ScopeComputation.computeExports(document, cancelToken);
        const uri = document.uri.toString();
        this.symbolIndex.set(uri, exports);
        this.symbolByTypeIndex.clear(uri);
    }
    async updateReferences(document, cancelToken = cancellation.CancellationToken.None) {
        const services = this.serviceRegistry.getServices(document.uri);
        const indexData = await services.workspace.ReferenceDescriptionProvider.createDescriptions(document, cancelToken);
        this.referenceIndex.set(document.uri.toString(), indexData);
    }
    isAffected(document, changedUris) {
        const references = this.referenceIndex.get(document.uri.toString());
        if (!references) {
            return false;
        }
        return references.some(ref => !ref.local && changedUris.has(ref.targetUri.toString()));
    }
}
//# sourceMappingURL=index-manager.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/workspace-manager.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/



class DefaultWorkspaceManager {
    constructor(services) {
        this.initialBuildOptions = {};
        this._ready = new promise_utils_Deferred();
        this.serviceRegistry = services.ServiceRegistry;
        this.langiumDocuments = services.workspace.LangiumDocuments;
        this.documentBuilder = services.workspace.DocumentBuilder;
        this.fileSystemProvider = services.workspace.FileSystemProvider;
        this.mutex = services.workspace.WorkspaceLock;
    }
    get ready() {
        return this._ready.promise;
    }
    get workspaceFolders() {
        return this.folders;
    }
    initialize(params) {
        var _a;
        this.folders = (_a = params.workspaceFolders) !== null && _a !== void 0 ? _a : undefined;
    }
    initialized(_params) {
        // Initialize the workspace even if there are no workspace folders
        // We still want to load additional documents (language library or similar) during initialization
        return this.mutex.write(token => { var _a; return this.initializeWorkspace((_a = this.folders) !== null && _a !== void 0 ? _a : [], token); });
    }
    async initializeWorkspace(folders, cancelToken = cancellation.CancellationToken.None) {
        const documents = await this.performStartup(folders);
        // Only after creating all documents do we check whether we need to cancel the initialization
        // The document builder will later pick up on all unprocessed documents
        await interruptAndCheck(cancelToken);
        await this.documentBuilder.build(documents, this.initialBuildOptions, cancelToken);
    }
    /**
     * Performs the uninterruptable startup sequence of the workspace manager.
     * This methods loads all documents in the workspace and other documents and returns them.
     */
    async performStartup(folders) {
        const fileExtensions = this.serviceRegistry.all.flatMap(e => e.LanguageMetaData.fileExtensions);
        const documents = [];
        const collector = (document) => {
            documents.push(document);
            if (!this.langiumDocuments.hasDocument(document.uri)) {
                this.langiumDocuments.addDocument(document);
            }
        };
        // Even though we don't await the initialization of the workspace manager,
        // we can still assume that all library documents and file documents are loaded by the time we start building documents.
        // The mutex prevents anything from performing a workspace build until we check the cancellation token
        await this.loadAdditionalDocuments(folders, collector);
        await Promise.all(folders.map(wf => [wf, this.getRootFolder(wf)])
            .map(async (entry) => this.traverseFolder(...entry, fileExtensions, collector)));
        this._ready.resolve();
        return documents;
    }
    /**
     * Load all additional documents that shall be visible in the context of the given workspace
     * folders and add them to the collector. This can be used to include built-in libraries of
     * your language, which can be either loaded from provided files or constructed in memory.
     */
    loadAdditionalDocuments(_folders, _collector) {
        return Promise.resolve();
    }
    /**
     * Determine the root folder of the source documents in the given workspace folder.
     * The default implementation returns the URI of the workspace folder, but you can override
     * this to return a subfolder like `src` instead.
     */
    getRootFolder(workspaceFolder) {
        return esm/* URI */.o.parse(workspaceFolder.uri);
    }
    /**
     * Traverse the file system folder identified by the given URI and its subfolders. All
     * contained files that match the file extensions are added to the collector.
     */
    async traverseFolder(workspaceFolder, folderPath, fileExtensions, collector) {
        const content = await this.fileSystemProvider.readDirectory(folderPath);
        await Promise.all(content.map(async (entry) => {
            if (this.includeEntry(workspaceFolder, entry, fileExtensions)) {
                if (entry.isDirectory) {
                    await this.traverseFolder(workspaceFolder, entry.uri, fileExtensions, collector);
                }
                else if (entry.isFile) {
                    const document = await this.langiumDocuments.getOrCreateDocument(entry.uri);
                    collector(document);
                }
            }
        }));
    }
    /**
     * Determine whether the given folder entry shall be included while indexing the workspace.
     */
    includeEntry(_workspaceFolder, entry, fileExtensions) {
        const name = UriUtils.basename(entry.uri);
        if (name.startsWith('.')) {
            return false;
        }
        if (entry.isDirectory) {
            return name !== 'node_modules' && name !== 'out';
        }
        else if (entry.isFile) {
            const extname = UriUtils.extname(entry.uri);
            return fileExtensions.includes(extname);
        }
        return false;
    }
}
//# sourceMappingURL=workspace-manager.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/lexer.js
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

class DefaultLexerErrorMessageProvider {
    buildUnexpectedCharactersMessage(fullText, startOffset, length, line, column) {
        return api/* defaultLexerErrorProvider */.ZW.buildUnexpectedCharactersMessage(fullText, startOffset, length, line, column);
    }
    buildUnableToPopLexerModeMessage(token) {
        return api/* defaultLexerErrorProvider */.ZW.buildUnableToPopLexerModeMessage(token);
    }
}
const DEFAULT_TOKENIZE_OPTIONS = { mode: 'full' };
class DefaultLexer {
    constructor(services) {
        this.errorMessageProvider = services.parser.LexerErrorMessageProvider;
        this.tokenBuilder = services.parser.TokenBuilder;
        const tokens = this.tokenBuilder.buildTokens(services.Grammar, {
            caseInsensitive: services.LanguageMetaData.caseInsensitive
        });
        this.tokenTypes = this.toTokenTypeDictionary(tokens);
        const lexerTokens = isTokenTypeDictionary(tokens) ? Object.values(tokens) : tokens;
        const production = services.LanguageMetaData.mode === 'production';
        this.chevrotainLexer = new api/* Lexer */.hW(lexerTokens, {
            positionTracking: 'full',
            skipValidations: production,
            errorMessageProvider: this.errorMessageProvider
        });
    }
    get definition() {
        return this.tokenTypes;
    }
    tokenize(text, _options = DEFAULT_TOKENIZE_OPTIONS) {
        var _a, _b, _c;
        const chevrotainResult = this.chevrotainLexer.tokenize(text);
        return {
            tokens: chevrotainResult.tokens,
            errors: chevrotainResult.errors,
            hidden: (_a = chevrotainResult.groups.hidden) !== null && _a !== void 0 ? _a : [],
            report: (_c = (_b = this.tokenBuilder).flushLexingReport) === null || _c === void 0 ? void 0 : _c.call(_b, text)
        };
    }
    toTokenTypeDictionary(buildTokens) {
        if (isTokenTypeDictionary(buildTokens))
            return buildTokens;
        const tokens = isIMultiModeLexerDefinition(buildTokens) ? Object.values(buildTokens.modes).flat() : buildTokens;
        const res = {};
        tokens.forEach(token => res[token.name] = token);
        return res;
    }
}
/**
 * Returns a check whether the given TokenVocabulary is TokenType array
 */
function isTokenTypeArray(tokenVocabulary) {
    return Array.isArray(tokenVocabulary) && (tokenVocabulary.length === 0 || 'name' in tokenVocabulary[0]);
}
/**
 * Returns a check whether the given TokenVocabulary is IMultiModeLexerDefinition
 */
function isIMultiModeLexerDefinition(tokenVocabulary) {
    return tokenVocabulary && 'modes' in tokenVocabulary && 'defaultMode' in tokenVocabulary;
}
/**
 * Returns a check whether the given TokenVocabulary is TokenTypeDictionary
 */
function isTokenTypeDictionary(tokenVocabulary) {
    return !isTokenTypeArray(tokenVocabulary) && !isIMultiModeLexerDefinition(tokenVocabulary);
}
//# sourceMappingURL=lexer.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/documentation/jsdoc.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/



function parseJSDoc(node, start, options) {
    let opts;
    let position;
    if (typeof node === 'string') {
        position = start;
        opts = options;
    }
    else {
        position = node.range.start;
        opts = start;
    }
    if (!position) {
        position = Position.create(0, 0);
    }
    const lines = getLines(node);
    const normalizedOptions = normalizeOptions(opts);
    const tokens = tokenize({
        lines,
        position,
        options: normalizedOptions
    });
    return parseJSDocComment({
        index: 0,
        tokens,
        position
    });
}
function isJSDoc(node, options) {
    const normalizedOptions = normalizeOptions(options);
    const lines = getLines(node);
    if (lines.length === 0) {
        return false;
    }
    const first = lines[0];
    const last = lines[lines.length - 1];
    const firstRegex = normalizedOptions.start;
    const lastRegex = normalizedOptions.end;
    return Boolean(firstRegex === null || firstRegex === void 0 ? void 0 : firstRegex.exec(first)) && Boolean(lastRegex === null || lastRegex === void 0 ? void 0 : lastRegex.exec(last));
}
function getLines(node) {
    let content = '';
    if (typeof node === 'string') {
        content = node;
    }
    else {
        content = node.text;
    }
    const lines = content.split(regexp_utils/* NEWLINE_REGEXP */.K0);
    return lines;
}
const tagRegex = /\s*(@([\p{L}][\p{L}\p{N}]*)?)/uy;
const inlineTagRegex = /\{(@[\p{L}][\p{L}\p{N}]*)(\s*)([^\r\n}]+)?\}/gu;
function tokenize(context) {
    var _a, _b, _c;
    const tokens = [];
    let currentLine = context.position.line;
    let currentCharacter = context.position.character;
    for (let i = 0; i < context.lines.length; i++) {
        const first = i === 0;
        const last = i === context.lines.length - 1;
        let line = context.lines[i];
        let index = 0;
        if (first && context.options.start) {
            const match = (_a = context.options.start) === null || _a === void 0 ? void 0 : _a.exec(line);
            if (match) {
                index = match.index + match[0].length;
            }
        }
        else {
            const match = (_b = context.options.line) === null || _b === void 0 ? void 0 : _b.exec(line);
            if (match) {
                index = match.index + match[0].length;
            }
        }
        if (last) {
            const match = (_c = context.options.end) === null || _c === void 0 ? void 0 : _c.exec(line);
            if (match) {
                line = line.substring(0, match.index);
            }
        }
        line = line.substring(0, lastCharacter(line));
        const whitespaceEnd = skipWhitespace(line, index);
        if (whitespaceEnd >= line.length) {
            // Only create a break token when we already have previous tokens
            if (tokens.length > 0) {
                const position = Position.create(currentLine, currentCharacter);
                tokens.push({
                    type: 'break',
                    content: '',
                    range: Range.create(position, position)
                });
            }
        }
        else {
            tagRegex.lastIndex = index;
            const tagMatch = tagRegex.exec(line);
            if (tagMatch) {
                const fullMatch = tagMatch[0];
                const value = tagMatch[1];
                const start = Position.create(currentLine, currentCharacter + index);
                const end = Position.create(currentLine, currentCharacter + index + fullMatch.length);
                tokens.push({
                    type: 'tag',
                    content: value,
                    range: Range.create(start, end)
                });
                index += fullMatch.length;
                index = skipWhitespace(line, index);
            }
            if (index < line.length) {
                const rest = line.substring(index);
                const inlineTagMatches = Array.from(rest.matchAll(inlineTagRegex));
                tokens.push(...buildInlineTokens(inlineTagMatches, rest, currentLine, currentCharacter + index));
            }
        }
        currentLine++;
        currentCharacter = 0;
    }
    // Remove last break token if there is one
    if (tokens.length > 0 && tokens[tokens.length - 1].type === 'break') {
        return tokens.slice(0, -1);
    }
    return tokens;
}
function buildInlineTokens(tags, line, lineIndex, characterIndex) {
    const tokens = [];
    if (tags.length === 0) {
        const start = Position.create(lineIndex, characterIndex);
        const end = Position.create(lineIndex, characterIndex + line.length);
        tokens.push({
            type: 'text',
            content: line,
            range: Range.create(start, end)
        });
    }
    else {
        let lastIndex = 0;
        for (const match of tags) {
            const matchIndex = match.index;
            const startContent = line.substring(lastIndex, matchIndex);
            if (startContent.length > 0) {
                tokens.push({
                    type: 'text',
                    content: line.substring(lastIndex, matchIndex),
                    range: Range.create(Position.create(lineIndex, lastIndex + characterIndex), Position.create(lineIndex, matchIndex + characterIndex))
                });
            }
            let offset = startContent.length + 1;
            const tagName = match[1];
            tokens.push({
                type: 'inline-tag',
                content: tagName,
                range: Range.create(Position.create(lineIndex, lastIndex + offset + characterIndex), Position.create(lineIndex, lastIndex + offset + tagName.length + characterIndex))
            });
            offset += tagName.length;
            if (match.length === 4) {
                offset += match[2].length;
                const value = match[3];
                tokens.push({
                    type: 'text',
                    content: value,
                    range: Range.create(Position.create(lineIndex, lastIndex + offset + characterIndex), Position.create(lineIndex, lastIndex + offset + value.length + characterIndex))
                });
            }
            else {
                tokens.push({
                    type: 'text',
                    content: '',
                    range: Range.create(Position.create(lineIndex, lastIndex + offset + characterIndex), Position.create(lineIndex, lastIndex + offset + characterIndex))
                });
            }
            lastIndex = matchIndex + match[0].length;
        }
        const endContent = line.substring(lastIndex);
        if (endContent.length > 0) {
            tokens.push({
                type: 'text',
                content: endContent,
                range: Range.create(Position.create(lineIndex, lastIndex + characterIndex), Position.create(lineIndex, lastIndex + characterIndex + endContent.length))
            });
        }
    }
    return tokens;
}
const nonWhitespaceRegex = /\S/;
const whitespaceEndRegex = /\s*$/;
function skipWhitespace(line, index) {
    const match = line.substring(index).match(nonWhitespaceRegex);
    if (match) {
        return index + match.index;
    }
    else {
        return line.length;
    }
}
function lastCharacter(line) {
    const match = line.match(whitespaceEndRegex);
    if (match && typeof match.index === 'number') {
        return match.index;
    }
    return undefined;
}
// Parsing
function parseJSDocComment(context) {
    var _a, _b, _c, _d;
    const startPosition = Position.create(context.position.line, context.position.character);
    if (context.tokens.length === 0) {
        return new JSDocCommentImpl([], Range.create(startPosition, startPosition));
    }
    const elements = [];
    while (context.index < context.tokens.length) {
        const element = parseJSDocElement(context, elements[elements.length - 1]);
        if (element) {
            elements.push(element);
        }
    }
    const start = (_b = (_a = elements[0]) === null || _a === void 0 ? void 0 : _a.range.start) !== null && _b !== void 0 ? _b : startPosition;
    const end = (_d = (_c = elements[elements.length - 1]) === null || _c === void 0 ? void 0 : _c.range.end) !== null && _d !== void 0 ? _d : startPosition;
    return new JSDocCommentImpl(elements, Range.create(start, end));
}
function parseJSDocElement(context, last) {
    const next = context.tokens[context.index];
    if (next.type === 'tag') {
        return parseJSDocTag(context, false);
    }
    else if (next.type === 'text' || next.type === 'inline-tag') {
        return parseJSDocText(context);
    }
    else {
        appendEmptyLine(next, last);
        context.index++;
        return undefined;
    }
}
function appendEmptyLine(token, element) {
    if (element) {
        const line = new JSDocLineImpl('', token.range);
        if ('inlines' in element) {
            element.inlines.push(line);
        }
        else {
            element.content.inlines.push(line);
        }
    }
}
function parseJSDocText(context) {
    let token = context.tokens[context.index];
    const firstToken = token;
    let lastToken = token;
    const lines = [];
    while (token && token.type !== 'break' && token.type !== 'tag') {
        lines.push(parseJSDocInline(context));
        lastToken = token;
        token = context.tokens[context.index];
    }
    return new JSDocTextImpl(lines, Range.create(firstToken.range.start, lastToken.range.end));
}
function parseJSDocInline(context) {
    const token = context.tokens[context.index];
    if (token.type === 'inline-tag') {
        return parseJSDocTag(context, true);
    }
    else {
        return parseJSDocLine(context);
    }
}
function parseJSDocTag(context, inline) {
    const tagToken = context.tokens[context.index++];
    const name = tagToken.content.substring(1);
    const nextToken = context.tokens[context.index];
    if ((nextToken === null || nextToken === void 0 ? void 0 : nextToken.type) === 'text') {
        if (inline) {
            const docLine = parseJSDocLine(context);
            return new JSDocTagImpl(name, new JSDocTextImpl([docLine], docLine.range), inline, Range.create(tagToken.range.start, docLine.range.end));
        }
        else {
            const textDoc = parseJSDocText(context);
            return new JSDocTagImpl(name, textDoc, inline, Range.create(tagToken.range.start, textDoc.range.end));
        }
    }
    else {
        const range = tagToken.range;
        return new JSDocTagImpl(name, new JSDocTextImpl([], range), inline, range);
    }
}
function parseJSDocLine(context) {
    const token = context.tokens[context.index++];
    return new JSDocLineImpl(token.content, token.range);
}
function normalizeOptions(options) {
    if (!options) {
        return normalizeOptions({
            start: '/**',
            end: '*/',
            line: '*'
        });
    }
    const { start, end, line } = options;
    return {
        start: normalizeOption(start, true),
        end: normalizeOption(end, false),
        line: normalizeOption(line, true)
    };
}
function normalizeOption(option, start) {
    if (typeof option === 'string' || typeof option === 'object') {
        const escaped = typeof option === 'string' ? (0,regexp_utils/* escapeRegExp */.hr)(option) : option.source;
        if (start) {
            return new RegExp(`^\\s*${escaped}`);
        }
        else {
            return new RegExp(`\\s*${escaped}\\s*$`);
        }
    }
    else {
        return option;
    }
}
class JSDocCommentImpl {
    constructor(elements, range) {
        this.elements = elements;
        this.range = range;
    }
    getTag(name) {
        return this.getAllTags().find(e => e.name === name);
    }
    getTags(name) {
        return this.getAllTags().filter(e => e.name === name);
    }
    getAllTags() {
        return this.elements.filter((e) => 'name' in e);
    }
    toString() {
        let value = '';
        for (const element of this.elements) {
            if (value.length === 0) {
                value = element.toString();
            }
            else {
                const text = element.toString();
                value += fillNewlines(value) + text;
            }
        }
        return value.trim();
    }
    toMarkdown(options) {
        let value = '';
        for (const element of this.elements) {
            if (value.length === 0) {
                value = element.toMarkdown(options);
            }
            else {
                const text = element.toMarkdown(options);
                value += fillNewlines(value) + text;
            }
        }
        return value.trim();
    }
}
class JSDocTagImpl {
    constructor(name, content, inline, range) {
        this.name = name;
        this.content = content;
        this.inline = inline;
        this.range = range;
    }
    toString() {
        let text = `@${this.name}`;
        const content = this.content.toString();
        if (this.content.inlines.length === 1) {
            text = `${text} ${content}`;
        }
        else if (this.content.inlines.length > 1) {
            text = `${text}\n${content}`;
        }
        if (this.inline) {
            // Inline tags are surrounded by curly braces
            return `{${text}}`;
        }
        else {
            return text;
        }
    }
    toMarkdown(options) {
        var _a, _b;
        return (_b = (_a = options === null || options === void 0 ? void 0 : options.renderTag) === null || _a === void 0 ? void 0 : _a.call(options, this)) !== null && _b !== void 0 ? _b : this.toMarkdownDefault(options);
    }
    toMarkdownDefault(options) {
        const content = this.content.toMarkdown(options);
        if (this.inline) {
            const rendered = renderInlineTag(this.name, content, options !== null && options !== void 0 ? options : {});
            if (typeof rendered === 'string') {
                return rendered;
            }
        }
        let marker = '';
        if ((options === null || options === void 0 ? void 0 : options.tag) === 'italic' || (options === null || options === void 0 ? void 0 : options.tag) === undefined) {
            marker = '*';
        }
        else if ((options === null || options === void 0 ? void 0 : options.tag) === 'bold') {
            marker = '**';
        }
        else if ((options === null || options === void 0 ? void 0 : options.tag) === 'bold-italic') {
            marker = '***';
        }
        let text = `${marker}@${this.name}${marker}`;
        if (this.content.inlines.length === 1) {
            text = `${text} — ${content}`;
        }
        else if (this.content.inlines.length > 1) {
            text = `${text}\n${content}`;
        }
        if (this.inline) {
            // Inline tags are surrounded by curly braces
            return `{${text}}`;
        }
        else {
            return text;
        }
    }
}
function renderInlineTag(tag, content, options) {
    var _a, _b;
    if (tag === 'linkplain' || tag === 'linkcode' || tag === 'link') {
        const index = content.indexOf(' ');
        let display = content;
        if (index > 0) {
            const displayStart = skipWhitespace(content, index);
            display = content.substring(displayStart);
            content = content.substring(0, index);
        }
        if (tag === 'linkcode' || (tag === 'link' && options.link === 'code')) {
            // Surround the display value in a markdown inline code block
            display = `\`${display}\``;
        }
        const renderedLink = (_b = (_a = options.renderLink) === null || _a === void 0 ? void 0 : _a.call(options, content, display)) !== null && _b !== void 0 ? _b : renderLinkDefault(content, display);
        return renderedLink;
    }
    return undefined;
}
function renderLinkDefault(content, display) {
    try {
        esm/* URI */.o.parse(content, true);
        return `[${display}](${content})`;
    }
    catch (_a) {
        return content;
    }
}
class JSDocTextImpl {
    constructor(lines, range) {
        this.inlines = lines;
        this.range = range;
    }
    toString() {
        let text = '';
        for (let i = 0; i < this.inlines.length; i++) {
            const inline = this.inlines[i];
            const next = this.inlines[i + 1];
            text += inline.toString();
            if (next && next.range.start.line > inline.range.start.line) {
                text += '\n';
            }
        }
        return text;
    }
    toMarkdown(options) {
        let text = '';
        for (let i = 0; i < this.inlines.length; i++) {
            const inline = this.inlines[i];
            const next = this.inlines[i + 1];
            text += inline.toMarkdown(options);
            if (next && next.range.start.line > inline.range.start.line) {
                text += '\n';
            }
        }
        return text;
    }
}
class JSDocLineImpl {
    constructor(text, range) {
        this.text = text;
        this.range = range;
    }
    toString() {
        return this.text;
    }
    toMarkdown() {
        return this.text;
    }
}
function fillNewlines(text) {
    if (text.endsWith('\n')) {
        return '\n';
    }
    else {
        return '\n\n';
    }
}
//# sourceMappingURL=jsdoc.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/documentation/documentation-provider.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


class JSDocDocumentationProvider {
    constructor(services) {
        this.indexManager = services.shared.workspace.IndexManager;
        this.commentProvider = services.documentation.CommentProvider;
    }
    getDocumentation(node) {
        const comment = this.commentProvider.getComment(node);
        if (comment && isJSDoc(comment)) {
            const parsedJSDoc = parseJSDoc(comment);
            return parsedJSDoc.toMarkdown({
                renderLink: (link, display) => {
                    return this.documentationLinkRenderer(node, link, display);
                },
                renderTag: (tag) => {
                    return this.documentationTagRenderer(node, tag);
                }
            });
        }
        return undefined;
    }
    documentationLinkRenderer(node, name, display) {
        var _a;
        const description = (_a = this.findNameInPrecomputedScopes(node, name)) !== null && _a !== void 0 ? _a : this.findNameInGlobalScope(node, name);
        if (description && description.nameSegment) {
            const line = description.nameSegment.range.start.line + 1;
            const character = description.nameSegment.range.start.character + 1;
            const uri = description.documentUri.with({ fragment: `L${line},${character}` });
            return `[${display}](${uri.toString()})`;
        }
        else {
            return undefined;
        }
    }
    documentationTagRenderer(_node, _tag) {
        // Fall back to the default tag rendering
        return undefined;
    }
    findNameInPrecomputedScopes(node, name) {
        const document = (0,ast_utils/* getDocument */.Me)(node);
        const precomputed = document.precomputedScopes;
        if (!precomputed) {
            return undefined;
        }
        let currentNode = node;
        do {
            const allDescriptions = precomputed.get(currentNode);
            const description = allDescriptions.find(e => e.name === name);
            if (description) {
                return description;
            }
            currentNode = currentNode.$container;
        } while (currentNode);
        return undefined;
    }
    findNameInGlobalScope(node, name) {
        const description = this.indexManager.allElements().find(e => e.name === name);
        return description;
    }
}
//# sourceMappingURL=documentation-provider.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/documentation/comment-provider.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


class DefaultCommentProvider {
    constructor(services) {
        this.grammarConfig = () => services.parser.GrammarConfig;
    }
    getComment(node) {
        var _a;
        if (isAstNodeWithComment(node)) {
            return node.$comment;
        }
        return (_a = (0,cst_utils/* findCommentNode */.LK)(node.$cstNode, this.grammarConfig().multilineCommentRules)) === null || _a === void 0 ? void 0 : _a.text;
    }
}
//# sourceMappingURL=comment-provider.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/parser/async-parser.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


/**
 * Default implementation of the async parser which simply wraps the sync parser in a promise.
 *
 * @remarks
 * A real implementation would create worker threads or web workers to offload the parsing work.
 */
class DefaultAsyncParser {
    constructor(services) {
        this.syncParser = services.parser.LangiumParser;
    }
    parse(text, _cancelToken) {
        return Promise.resolve(this.syncParser.parse(text));
    }
}
class AbstractThreadedAsyncParser {
    constructor(services) {
        /**
         * The thread count determines how many threads are used to parse files in parallel.
         * The default value is 8. Decreasing this value increases startup performance, but decreases parallel parsing performance.
         */
        this.threadCount = 8;
        /**
         * The termination delay determines how long the parser waits for a thread to finish after a cancellation request.
         * The default value is 200(ms).
         */
        this.terminationDelay = 200;
        this.workerPool = [];
        this.queue = [];
        this.hydrator = services.serializer.Hydrator;
    }
    initializeWorkers() {
        while (this.workerPool.length < this.threadCount) {
            const worker = this.createWorker();
            worker.onReady(() => {
                if (this.queue.length > 0) {
                    const deferred = this.queue.shift();
                    if (deferred) {
                        worker.lock();
                        deferred.resolve(worker);
                    }
                }
            });
            this.workerPool.push(worker);
        }
    }
    async parse(text, cancelToken) {
        const worker = await this.acquireParserWorker(cancelToken);
        const deferred = new Deferred();
        let timeout;
        // If the cancellation token is requested, we wait for a certain time before terminating the worker.
        // Since the cancellation token lives longer than the parsing process, we need to dispose the event listener.
        // Otherwise, we might accidentally terminate the worker after the parsing process has finished.
        const cancellation = cancelToken.onCancellationRequested(() => {
            timeout = setTimeout(() => {
                this.terminateWorker(worker);
            }, this.terminationDelay);
        });
        worker.parse(text).then(result => {
            const hydrated = this.hydrator.hydrate(result);
            deferred.resolve(hydrated);
        }).catch(err => {
            deferred.reject(err);
        }).finally(() => {
            cancellation.dispose();
            clearTimeout(timeout);
        });
        return deferred.promise;
    }
    terminateWorker(worker) {
        worker.terminate();
        const index = this.workerPool.indexOf(worker);
        if (index >= 0) {
            this.workerPool.splice(index, 1);
        }
    }
    async acquireParserWorker(cancelToken) {
        this.initializeWorkers();
        for (const worker of this.workerPool) {
            if (worker.ready) {
                worker.lock();
                return worker;
            }
        }
        const deferred = new Deferred();
        cancelToken.onCancellationRequested(() => {
            const index = this.queue.indexOf(deferred);
            if (index >= 0) {
                this.queue.splice(index, 1);
            }
            deferred.reject(OperationCancelled);
        });
        this.queue.push(deferred);
        return deferred.promise;
    }
}
class ParserWorker {
    get ready() {
        return this._ready;
    }
    get onReady() {
        return this.onReadyEmitter.event;
    }
    constructor(sendMessage, onMessage, onError, terminate) {
        this.onReadyEmitter = new Emitter();
        this.deferred = new Deferred();
        this._ready = true;
        this._parsing = false;
        this.sendMessage = sendMessage;
        this._terminate = terminate;
        onMessage(result => {
            const parseResult = result;
            this.deferred.resolve(parseResult);
            this.unlock();
        });
        onError(error => {
            this.deferred.reject(error);
            this.unlock();
        });
    }
    terminate() {
        this.deferred.reject(OperationCancelled);
        this._terminate();
    }
    lock() {
        this._ready = false;
    }
    unlock() {
        this._parsing = false;
        this._ready = true;
        this.onReadyEmitter.fire();
    }
    parse(text) {
        if (this._parsing) {
            throw new Error('Parser worker is busy');
        }
        this._parsing = true;
        this.deferred = new Deferred();
        this.sendMessage(text);
        return this.deferred.promise;
    }
}
//# sourceMappingURL=async-parser.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/workspace/workspace-lock.js
/******************************************************************************
 * Copyright 2023 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


class DefaultWorkspaceLock {
    constructor() {
        this.previousTokenSource = new cancellation.CancellationTokenSource();
        this.writeQueue = [];
        this.readQueue = [];
        this.done = true;
    }
    write(action) {
        this.cancelWrite();
        const tokenSource = startCancelableOperation();
        this.previousTokenSource = tokenSource;
        return this.enqueue(this.writeQueue, action, tokenSource.token);
    }
    read(action) {
        return this.enqueue(this.readQueue, action);
    }
    enqueue(queue, action, cancellationToken = cancellation.CancellationToken.None) {
        const deferred = new promise_utils_Deferred();
        const entry = {
            action,
            deferred,
            cancellationToken
        };
        queue.push(entry);
        this.performNextOperation();
        return deferred.promise;
    }
    async performNextOperation() {
        if (!this.done) {
            return;
        }
        const entries = [];
        if (this.writeQueue.length > 0) {
            // Just perform the next write action
            entries.push(this.writeQueue.shift());
        }
        else if (this.readQueue.length > 0) {
            // Empty the read queue and perform all actions in parallel
            entries.push(...this.readQueue.splice(0, this.readQueue.length));
        }
        else {
            return;
        }
        this.done = false;
        await Promise.all(entries.map(async ({ action, deferred, cancellationToken }) => {
            try {
                // Move the execution of the action to the next event loop tick via `Promise.resolve()`
                const result = await Promise.resolve().then(() => action(cancellationToken));
                deferred.resolve(result);
            }
            catch (err) {
                if (isOperationCancelled(err)) {
                    // If the operation was cancelled, we don't want to reject the promise
                    deferred.resolve(undefined);
                }
                else {
                    deferred.reject(err);
                }
            }
        }));
        this.done = true;
        this.performNextOperation();
    }
    cancelWrite() {
        this.previousTokenSource.cancel();
    }
}
//# sourceMappingURL=workspace-lock.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/serializer/hydrator.js
/******************************************************************************
 * Copyright 2024 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






class DefaultHydrator {
    constructor(services) {
        this.grammarElementIdMap = new BiMap();
        this.tokenTypeIdMap = new BiMap();
        this.grammar = services.Grammar;
        this.lexer = services.parser.Lexer;
        this.linker = services.references.Linker;
    }
    dehydrate(result) {
        return {
            lexerErrors: result.lexerErrors,
            lexerReport: result.lexerReport ? this.dehydrateLexerReport(result.lexerReport) : undefined,
            // We need to create shallow copies of the errors
            // The original errors inherit from the `Error` class, which is not transferable across worker threads
            parserErrors: result.parserErrors.map(e => (Object.assign(Object.assign({}, e), { message: e.message }))),
            value: this.dehydrateAstNode(result.value, this.createDehyrationContext(result.value))
        };
    }
    dehydrateLexerReport(lexerReport) {
        // By default, lexer reports are serializable
        return lexerReport;
    }
    createDehyrationContext(node) {
        const astNodes = new Map();
        const cstNodes = new Map();
        for (const astNode of (0,ast_utils/* streamAst */.Zc)(node)) {
            astNodes.set(astNode, {});
        }
        if (node.$cstNode) {
            for (const cstNode of (0,cst_utils/* streamCst */._t)(node.$cstNode)) {
                cstNodes.set(cstNode, {});
            }
        }
        return {
            astNodes,
            cstNodes
        };
    }
    dehydrateAstNode(node, context) {
        const obj = context.astNodes.get(node);
        obj.$type = node.$type;
        obj.$containerIndex = node.$containerIndex;
        obj.$containerProperty = node.$containerProperty;
        if (node.$cstNode !== undefined) {
            obj.$cstNode = this.dehydrateCstNode(node.$cstNode, context);
        }
        for (const [name, value] of Object.entries(node)) {
            if (name.startsWith('$')) {
                continue;
            }
            if (Array.isArray(value)) {
                const arr = [];
                obj[name] = arr;
                for (const item of value) {
                    if ((0,syntax_tree/* isAstNode */.xA)(item)) {
                        arr.push(this.dehydrateAstNode(item, context));
                    }
                    else if ((0,syntax_tree/* isReference */.Yk)(item)) {
                        arr.push(this.dehydrateReference(item, context));
                    }
                    else {
                        arr.push(item);
                    }
                }
            }
            else if ((0,syntax_tree/* isAstNode */.xA)(value)) {
                obj[name] = this.dehydrateAstNode(value, context);
            }
            else if ((0,syntax_tree/* isReference */.Yk)(value)) {
                obj[name] = this.dehydrateReference(value, context);
            }
            else if (value !== undefined) {
                obj[name] = value;
            }
        }
        return obj;
    }
    dehydrateReference(reference, context) {
        const obj = {};
        obj.$refText = reference.$refText;
        if (reference.$refNode) {
            obj.$refNode = context.cstNodes.get(reference.$refNode);
        }
        return obj;
    }
    dehydrateCstNode(node, context) {
        const cstNode = context.cstNodes.get(node);
        if ((0,syntax_tree/* isRootCstNode */.U8)(node)) {
            cstNode.fullText = node.fullText;
        }
        else {
            // Note: This returns undefined for hidden nodes (i.e. comments)
            cstNode.grammarSource = this.getGrammarElementId(node.grammarSource);
        }
        cstNode.hidden = node.hidden;
        cstNode.astNode = context.astNodes.get(node.astNode);
        if ((0,syntax_tree/* isCompositeCstNode */.al)(node)) {
            cstNode.content = node.content.map(child => this.dehydrateCstNode(child, context));
        }
        else if ((0,syntax_tree/* isLeafCstNode */.dm)(node)) {
            cstNode.tokenType = node.tokenType.name;
            cstNode.offset = node.offset;
            cstNode.length = node.length;
            cstNode.startLine = node.range.start.line;
            cstNode.startColumn = node.range.start.character;
            cstNode.endLine = node.range.end.line;
            cstNode.endColumn = node.range.end.character;
        }
        return cstNode;
    }
    hydrate(result) {
        const node = result.value;
        const context = this.createHydrationContext(node);
        if ('$cstNode' in node) {
            this.hydrateCstNode(node.$cstNode, context);
        }
        return {
            lexerErrors: result.lexerErrors,
            lexerReport: result.lexerReport,
            parserErrors: result.parserErrors,
            value: this.hydrateAstNode(node, context)
        };
    }
    createHydrationContext(node) {
        const astNodes = new Map();
        const cstNodes = new Map();
        for (const astNode of (0,ast_utils/* streamAst */.Zc)(node)) {
            astNodes.set(astNode, {});
        }
        let root;
        if (node.$cstNode) {
            for (const cstNode of (0,cst_utils/* streamCst */._t)(node.$cstNode)) {
                let cst;
                if ('fullText' in cstNode) {
                    cst = new RootCstNodeImpl(cstNode.fullText);
                    root = cst;
                }
                else if ('content' in cstNode) {
                    cst = new CompositeCstNodeImpl();
                }
                else if ('tokenType' in cstNode) {
                    cst = this.hydrateCstLeafNode(cstNode);
                }
                if (cst) {
                    cstNodes.set(cstNode, cst);
                    cst.root = root;
                }
            }
        }
        return {
            astNodes,
            cstNodes
        };
    }
    hydrateAstNode(node, context) {
        const astNode = context.astNodes.get(node);
        astNode.$type = node.$type;
        astNode.$containerIndex = node.$containerIndex;
        astNode.$containerProperty = node.$containerProperty;
        if (node.$cstNode) {
            astNode.$cstNode = context.cstNodes.get(node.$cstNode);
        }
        for (const [name, value] of Object.entries(node)) {
            if (name.startsWith('$')) {
                continue;
            }
            if (Array.isArray(value)) {
                const arr = [];
                astNode[name] = arr;
                for (const item of value) {
                    if ((0,syntax_tree/* isAstNode */.xA)(item)) {
                        arr.push(this.setParent(this.hydrateAstNode(item, context), astNode));
                    }
                    else if ((0,syntax_tree/* isReference */.Yk)(item)) {
                        arr.push(this.hydrateReference(item, astNode, name, context));
                    }
                    else {
                        arr.push(item);
                    }
                }
            }
            else if ((0,syntax_tree/* isAstNode */.xA)(value)) {
                astNode[name] = this.setParent(this.hydrateAstNode(value, context), astNode);
            }
            else if ((0,syntax_tree/* isReference */.Yk)(value)) {
                astNode[name] = this.hydrateReference(value, astNode, name, context);
            }
            else if (value !== undefined) {
                astNode[name] = value;
            }
        }
        return astNode;
    }
    setParent(node, parent) {
        node.$container = parent;
        return node;
    }
    hydrateReference(reference, node, name, context) {
        return this.linker.buildReference(node, name, context.cstNodes.get(reference.$refNode), reference.$refText);
    }
    hydrateCstNode(cstNode, context, num = 0) {
        const cstNodeObj = context.cstNodes.get(cstNode);
        if (typeof cstNode.grammarSource === 'number') {
            cstNodeObj.grammarSource = this.getGrammarElement(cstNode.grammarSource);
        }
        cstNodeObj.astNode = context.astNodes.get(cstNode.astNode);
        if ((0,syntax_tree/* isCompositeCstNode */.al)(cstNodeObj)) {
            for (const child of cstNode.content) {
                const hydrated = this.hydrateCstNode(child, context, num++);
                cstNodeObj.content.push(hydrated);
            }
        }
        return cstNodeObj;
    }
    hydrateCstLeafNode(cstNode) {
        const tokenType = this.getTokenType(cstNode.tokenType);
        const offset = cstNode.offset;
        const length = cstNode.length;
        const startLine = cstNode.startLine;
        const startColumn = cstNode.startColumn;
        const endLine = cstNode.endLine;
        const endColumn = cstNode.endColumn;
        const hidden = cstNode.hidden;
        const node = new LeafCstNodeImpl(offset, length, {
            start: {
                line: startLine,
                character: startColumn
            },
            end: {
                line: endLine,
                character: endColumn
            }
        }, tokenType, hidden);
        return node;
    }
    getTokenType(name) {
        return this.lexer.definition[name];
    }
    getGrammarElementId(node) {
        if (!node) {
            return undefined;
        }
        if (this.grammarElementIdMap.size === 0) {
            this.createGrammarElementIdMap();
        }
        return this.grammarElementIdMap.get(node);
    }
    getGrammarElement(id) {
        if (this.grammarElementIdMap.size === 0) {
            this.createGrammarElementIdMap();
        }
        const element = this.grammarElementIdMap.getKey(id);
        return element;
    }
    createGrammarElementIdMap() {
        let id = 0;
        for (const element of (0,ast_utils/* streamAst */.Zc)(this.grammar)) {
            if ((0,ast/* isAbstractElement */.zJ)(element)) {
                this.grammarElementIdMap.set(element, id++);
            }
        }
    }
}
//# sourceMappingURL=hydrator.js.map
;// CONCATENATED MODULE: ../node_modules/langium/lib/default-module.js
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
******************************************************************************/




























/**
 * Creates a dependency injection module configuring the default core services.
 * This is a set of services that are dedicated to a specific language.
 */
function createDefaultCoreModule(context) {
    return {
        documentation: {
            CommentProvider: (services) => new DefaultCommentProvider(services),
            DocumentationProvider: (services) => new JSDocDocumentationProvider(services)
        },
        parser: {
            AsyncParser: (services) => new DefaultAsyncParser(services),
            GrammarConfig: (services) => createGrammarConfig(services),
            LangiumParser: (services) => createLangiumParser(services),
            CompletionParser: (services) => createCompletionParser(services),
            ValueConverter: () => new value_converter/* DefaultValueConverter */.t(),
            TokenBuilder: () => new token_builder/* DefaultTokenBuilder */.P(),
            Lexer: (services) => new DefaultLexer(services),
            ParserErrorMessageProvider: () => new LangiumParserErrorMessageProvider(),
            LexerErrorMessageProvider: () => new DefaultLexerErrorMessageProvider()
        },
        workspace: {
            AstNodeLocator: () => new DefaultAstNodeLocator(),
            AstNodeDescriptionProvider: (services) => new DefaultAstNodeDescriptionProvider(services),
            ReferenceDescriptionProvider: (services) => new DefaultReferenceDescriptionProvider(services)
        },
        references: {
            Linker: (services) => new DefaultLinker(services),
            NameProvider: () => new DefaultNameProvider(),
            ScopeProvider: (services) => new DefaultScopeProvider(services),
            ScopeComputation: (services) => new DefaultScopeComputation(services),
            References: (services) => new DefaultReferences(services)
        },
        serializer: {
            Hydrator: (services) => new DefaultHydrator(services),
            JsonSerializer: (services) => new DefaultJsonSerializer(services)
        },
        validation: {
            DocumentValidator: (services) => new DefaultDocumentValidator(services),
            ValidationRegistry: (services) => new ValidationRegistry(services)
        },
        shared: () => context.shared
    };
}
/**
 * Creates a dependency injection module configuring the default shared core services.
 * This is the set of services that are shared between multiple languages.
 */
function createDefaultSharedCoreModule(context) {
    return {
        ServiceRegistry: (services) => new DefaultServiceRegistry(services),
        workspace: {
            LangiumDocuments: (services) => new DefaultLangiumDocuments(services),
            LangiumDocumentFactory: (services) => new DefaultLangiumDocumentFactory(services),
            DocumentBuilder: (services) => new DefaultDocumentBuilder(services),
            IndexManager: (services) => new DefaultIndexManager(services),
            WorkspaceManager: (services) => new DefaultWorkspaceManager(services),
            FileSystemProvider: (services) => context.fileSystemProvider(services),
            WorkspaceLock: () => new DefaultWorkspaceLock(),
            ConfigurationProvider: (services) => new DefaultConfigurationProvider(services)
        }
    };
}
//# sourceMappingURL=default-module.js.map

/***/ }),

/***/ 81210:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   f3: () => (/* binding */ inject)
/* harmony export */ });
/* unused harmony exports Module, eagerLoad */
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
var Module;
(function (Module) {
    Module.merge = (m1, m2) => _merge(_merge({}, m1), m2);
})(Module || (Module = {}));
/**
 * Given a set of modules, the inject function returns a lazily evaluated injector
 * that injects dependencies into the requested service when it is requested the
 * first time. Subsequent requests will return the same service.
 *
 * In the case of cyclic dependencies, an Error will be thrown. This can be fixed
 * by injecting a provider `() => T` instead of a `T`.
 *
 * Please note that the arguments may be objects or arrays. However, the result will
 * be an object. Using it with for..of will have no effect.
 *
 * @param module1 first Module
 * @param module2 (optional) second Module
 * @param module3 (optional) third Module
 * @param module4 (optional) fourth Module
 * @param module5 (optional) fifth Module
 * @param module6 (optional) sixth Module
 * @param module7 (optional) seventh Module
 * @param module8 (optional) eighth Module
 * @param module9 (optional) ninth Module
 * @returns a new object of type I
 */
function inject(module1, module2, module3, module4, module5, module6, module7, module8, module9) {
    const module = [module1, module2, module3, module4, module5, module6, module7, module8, module9].reduce(_merge, {});
    return _inject(module);
}
const isProxy = Symbol('isProxy');
/**
 * Eagerly load all services in the given dependency injection container. This is sometimes
 * necessary because services can register event listeners in their constructors.
 */
function eagerLoad(item) {
    if (item && item[isProxy]) {
        for (const value of Object.values(item)) {
            eagerLoad(value);
        }
    }
    return item;
}
/**
 * Helper function that returns an injector by creating a proxy.
 * Invariant: injector is of type I. If injector is undefined, then T = I.
 */
function _inject(module, injector) {
    const proxy = new Proxy({}, {
        deleteProperty: () => false,
        set: () => {
            throw new Error('Cannot set property on injected service container');
        },
        get: (obj, prop) => {
            if (prop === isProxy) {
                return true;
            }
            else {
                return _resolve(obj, prop, module, injector || proxy);
            }
        },
        getOwnPropertyDescriptor: (obj, prop) => (_resolve(obj, prop, module, injector || proxy), Object.getOwnPropertyDescriptor(obj, prop)), // used by for..in
        has: (_, prop) => prop in module, // used by ..in..
        ownKeys: () => [...Object.getOwnPropertyNames(module)] // used by for..in
    });
    return proxy;
}
/**
 * Internally used to tag a requested dependency, directly before calling the factory.
 * This allows us to find cycles during instance creation.
 */
const __requested__ = Symbol();
/**
 * Returns the value `obj[prop]`. If the value does not exist, yet, it is resolved from
 * the module description. The result of service factories is cached. Groups are
 * recursively proxied.
 *
 * @param obj an object holding all group proxies and services
 * @param prop the key of a value within obj
 * @param module an object containing groups and service factories
 * @param injector the first level proxy that provides access to all values
 * @returns the requested value `obj[prop]`
 * @throws Error if a dependency cycle is detected
 */
function _resolve(obj, prop, module, injector) {
    if (prop in obj) {
        if (obj[prop] instanceof Error) {
            throw new Error('Construction failure. Please make sure that your dependencies are constructable.', { cause: obj[prop] });
        }
        if (obj[prop] === __requested__) {
            throw new Error('Cycle detected. Please make "' + String(prop) + '" lazy. Visit https://langium.org/docs/reference/configuration-services/#resolving-cyclic-dependencies');
        }
        return obj[prop];
    }
    else if (prop in module) {
        const value = module[prop];
        obj[prop] = __requested__;
        try {
            obj[prop] = (typeof value === 'function') ? value(injector) : _inject(value, injector);
        }
        catch (error) {
            obj[prop] = error instanceof Error ? error : undefined;
            throw error;
        }
        return obj[prop];
    }
    else {
        return undefined;
    }
}
/**
 * Performs a deep-merge of two modules by writing source entries into the target module.
 *
 * @param target the module which is written
 * @param source the module which is read
 * @returns the target module
 */
function _merge(target, source) {
    if (source) {
        for (const [key, value2] of Object.entries(source)) {
            if (value2 !== undefined) {
                const value1 = target[key];
                if (value1 !== null && value2 !== null && typeof value1 === 'object' && typeof value2 === 'object') {
                    target[key] = _merge(value1, value2);
                }
                else {
                    target[key] = value2;
                }
            }
        }
    }
    return target;
}
//# sourceMappingURL=dependency-injection.js.map

/***/ }),

/***/ 34905:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   B7: () => (/* binding */ isAssignment),
/* harmony export */   Bf: () => (/* binding */ isCharacterRange),
/* harmony export */   Bi: () => (/* binding */ isNegatedToken),
/* harmony export */   F8: () => (/* binding */ isDisjunction),
/* harmony export */   F9: () => (/* binding */ isParserRule),
/* harmony export */   Ii: () => (/* binding */ isNegation),
/* harmony export */   Iy: () => (/* binding */ isSimpleType),
/* harmony export */   Ki: () => (/* binding */ isCrossReference),
/* harmony export */   L: () => (/* binding */ isBooleanLiteral),
/* harmony export */   LG: () => (/* binding */ isAction),
/* harmony export */   MS: () => (/* binding */ isTerminalRule),
/* harmony export */   MZ: () => (/* binding */ isAlternatives),
/* harmony export */   Mp: () => (/* binding */ isReturnType),
/* harmony export */   OG: () => (/* binding */ isUntilToken),
/* harmony export */   P9: () => (/* binding */ isType),
/* harmony export */   QV: () => (/* binding */ isInterface),
/* harmony export */   SV: () => (/* binding */ LangiumGrammarAstReflection),
/* harmony export */   S_: () => (/* binding */ isInferredType),
/* harmony export */   Sg: () => (/* binding */ isRegexToken),
/* harmony export */   TB: () => (/* binding */ isConjunction),
/* harmony export */   V7: () => (/* binding */ isTerminalAlternatives),
/* harmony export */   W1: () => (/* binding */ isUnorderedGroup),
/* harmony export */   X9: () => (/* binding */ isTerminalGroup),
/* harmony export */   gf: () => (/* binding */ isTerminalRuleCall),
/* harmony export */   p1: () => (/* binding */ isKeyword),
/* harmony export */   qm: () => (/* binding */ isWildcard),
/* harmony export */   rT: () => (/* binding */ isEndOfFile),
/* harmony export */   t3: () => (/* binding */ isRuleCall),
/* harmony export */   ty: () => (/* binding */ isGroup),
/* harmony export */   yW: () => (/* binding */ isParameterReference),
/* harmony export */   zJ: () => (/* binding */ isAbstractElement)
/* harmony export */ });
/* unused harmony exports LangiumGrammarTerminals, AbstractRule, isAbstractRule, AbstractType, isAbstractType, Condition, isCondition, isFeatureName, isPrimitiveType, TypeDefinition, isTypeDefinition, ValueLiteral, isValueLiteral, AbstractElement, ArrayLiteral, isArrayLiteral, ArrayType, isArrayType, BooleanLiteral, Conjunction, Disjunction, Grammar, isGrammar, GrammarImport, isGrammarImport, InferredType, Interface, NamedArgument, isNamedArgument, Negation, NumberLiteral, isNumberLiteral, Parameter, isParameter, ParameterReference, ParserRule, ReferenceType, isReferenceType, ReturnType, SimpleType, StringLiteral, isStringLiteral, TerminalRule, Type, TypeAttribute, isTypeAttribute, UnionType, isUnionType, Action, Alternatives, Assignment, CharacterRange, CrossReference, EndOfFile, Group, Keyword, NegatedToken, RegexToken, RuleCall, TerminalAlternatives, TerminalGroup, TerminalRuleCall, UnorderedGroup, UntilToken, Wildcard, reflection */
/* harmony import */ var _syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(91303);
/******************************************************************************
 * This file was generated by langium-cli 3.3.0.
 * DO NOT EDIT MANUALLY!
 ******************************************************************************/

const LangiumGrammarTerminals = {
    ID: /\^?[_a-zA-Z][\w_]*/,
    STRING: /"(\\.|[^"\\])*"|'(\\.|[^'\\])*'/,
    NUMBER: /NaN|-?((\d*\.\d+|\d+)([Ee][+-]?\d+)?|Infinity)/,
    RegexLiteral: /\/(?![*+?])(?:[^\r\n\[/\\]|\\.|\[(?:[^\r\n\]\\]|\\.)*\])+\/[a-z]*/,
    WS: /\s+/,
    ML_COMMENT: /\/\*[\s\S]*?\*\//,
    SL_COMMENT: /\/\/[^\n\r]*/,
};
const AbstractRule = 'AbstractRule';
function isAbstractRule(item) {
    return reflection.isInstance(item, AbstractRule);
}
const AbstractType = 'AbstractType';
function isAbstractType(item) {
    return reflection.isInstance(item, AbstractType);
}
const Condition = 'Condition';
function isCondition(item) {
    return reflection.isInstance(item, Condition);
}
function isFeatureName(item) {
    return isPrimitiveType(item) || item === 'current' || item === 'entry' || item === 'extends' || item === 'false' || item === 'fragment' || item === 'grammar' || item === 'hidden' || item === 'import' || item === 'interface' || item === 'returns' || item === 'terminal' || item === 'true' || item === 'type' || item === 'infer' || item === 'infers' || item === 'with' || (typeof item === 'string' && (/\^?[_a-zA-Z][\w_]*/.test(item)));
}
function isPrimitiveType(item) {
    return item === 'string' || item === 'number' || item === 'boolean' || item === 'Date' || item === 'bigint';
}
const TypeDefinition = 'TypeDefinition';
function isTypeDefinition(item) {
    return reflection.isInstance(item, TypeDefinition);
}
const ValueLiteral = 'ValueLiteral';
function isValueLiteral(item) {
    return reflection.isInstance(item, ValueLiteral);
}
const AbstractElement = 'AbstractElement';
function isAbstractElement(item) {
    return reflection.isInstance(item, AbstractElement);
}
const ArrayLiteral = 'ArrayLiteral';
function isArrayLiteral(item) {
    return reflection.isInstance(item, ArrayLiteral);
}
const ArrayType = 'ArrayType';
function isArrayType(item) {
    return reflection.isInstance(item, ArrayType);
}
const BooleanLiteral = 'BooleanLiteral';
function isBooleanLiteral(item) {
    return reflection.isInstance(item, BooleanLiteral);
}
const Conjunction = 'Conjunction';
function isConjunction(item) {
    return reflection.isInstance(item, Conjunction);
}
const Disjunction = 'Disjunction';
function isDisjunction(item) {
    return reflection.isInstance(item, Disjunction);
}
const Grammar = 'Grammar';
function isGrammar(item) {
    return reflection.isInstance(item, Grammar);
}
const GrammarImport = 'GrammarImport';
function isGrammarImport(item) {
    return reflection.isInstance(item, GrammarImport);
}
const InferredType = 'InferredType';
function isInferredType(item) {
    return reflection.isInstance(item, InferredType);
}
const Interface = 'Interface';
function isInterface(item) {
    return reflection.isInstance(item, Interface);
}
const NamedArgument = 'NamedArgument';
function isNamedArgument(item) {
    return reflection.isInstance(item, NamedArgument);
}
const Negation = 'Negation';
function isNegation(item) {
    return reflection.isInstance(item, Negation);
}
const NumberLiteral = 'NumberLiteral';
function isNumberLiteral(item) {
    return reflection.isInstance(item, NumberLiteral);
}
const Parameter = 'Parameter';
function isParameter(item) {
    return reflection.isInstance(item, Parameter);
}
const ParameterReference = 'ParameterReference';
function isParameterReference(item) {
    return reflection.isInstance(item, ParameterReference);
}
const ParserRule = 'ParserRule';
function isParserRule(item) {
    return reflection.isInstance(item, ParserRule);
}
const ReferenceType = 'ReferenceType';
function isReferenceType(item) {
    return reflection.isInstance(item, ReferenceType);
}
const ReturnType = 'ReturnType';
function isReturnType(item) {
    return reflection.isInstance(item, ReturnType);
}
const SimpleType = 'SimpleType';
function isSimpleType(item) {
    return reflection.isInstance(item, SimpleType);
}
const StringLiteral = 'StringLiteral';
function isStringLiteral(item) {
    return reflection.isInstance(item, StringLiteral);
}
const TerminalRule = 'TerminalRule';
function isTerminalRule(item) {
    return reflection.isInstance(item, TerminalRule);
}
const Type = 'Type';
function isType(item) {
    return reflection.isInstance(item, Type);
}
const TypeAttribute = 'TypeAttribute';
function isTypeAttribute(item) {
    return reflection.isInstance(item, TypeAttribute);
}
const UnionType = 'UnionType';
function isUnionType(item) {
    return reflection.isInstance(item, UnionType);
}
const Action = 'Action';
function isAction(item) {
    return reflection.isInstance(item, Action);
}
const Alternatives = 'Alternatives';
function isAlternatives(item) {
    return reflection.isInstance(item, Alternatives);
}
const Assignment = 'Assignment';
function isAssignment(item) {
    return reflection.isInstance(item, Assignment);
}
const CharacterRange = 'CharacterRange';
function isCharacterRange(item) {
    return reflection.isInstance(item, CharacterRange);
}
const CrossReference = 'CrossReference';
function isCrossReference(item) {
    return reflection.isInstance(item, CrossReference);
}
const EndOfFile = 'EndOfFile';
function isEndOfFile(item) {
    return reflection.isInstance(item, EndOfFile);
}
const Group = 'Group';
function isGroup(item) {
    return reflection.isInstance(item, Group);
}
const Keyword = 'Keyword';
function isKeyword(item) {
    return reflection.isInstance(item, Keyword);
}
const NegatedToken = 'NegatedToken';
function isNegatedToken(item) {
    return reflection.isInstance(item, NegatedToken);
}
const RegexToken = 'RegexToken';
function isRegexToken(item) {
    return reflection.isInstance(item, RegexToken);
}
const RuleCall = 'RuleCall';
function isRuleCall(item) {
    return reflection.isInstance(item, RuleCall);
}
const TerminalAlternatives = 'TerminalAlternatives';
function isTerminalAlternatives(item) {
    return reflection.isInstance(item, TerminalAlternatives);
}
const TerminalGroup = 'TerminalGroup';
function isTerminalGroup(item) {
    return reflection.isInstance(item, TerminalGroup);
}
const TerminalRuleCall = 'TerminalRuleCall';
function isTerminalRuleCall(item) {
    return reflection.isInstance(item, TerminalRuleCall);
}
const UnorderedGroup = 'UnorderedGroup';
function isUnorderedGroup(item) {
    return reflection.isInstance(item, UnorderedGroup);
}
const UntilToken = 'UntilToken';
function isUntilToken(item) {
    return reflection.isInstance(item, UntilToken);
}
const Wildcard = 'Wildcard';
function isWildcard(item) {
    return reflection.isInstance(item, Wildcard);
}
class LangiumGrammarAstReflection extends _syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .AbstractAstReflection */ .$v {
    getAllTypes() {
        return [AbstractElement, AbstractRule, AbstractType, Action, Alternatives, ArrayLiteral, ArrayType, Assignment, BooleanLiteral, CharacterRange, Condition, Conjunction, CrossReference, Disjunction, EndOfFile, Grammar, GrammarImport, Group, InferredType, Interface, Keyword, NamedArgument, NegatedToken, Negation, NumberLiteral, Parameter, ParameterReference, ParserRule, ReferenceType, RegexToken, ReturnType, RuleCall, SimpleType, StringLiteral, TerminalAlternatives, TerminalGroup, TerminalRule, TerminalRuleCall, Type, TypeAttribute, TypeDefinition, UnionType, UnorderedGroup, UntilToken, ValueLiteral, Wildcard];
    }
    computeIsSubtype(subtype, supertype) {
        switch (subtype) {
            case Action:
            case Alternatives:
            case Assignment:
            case CharacterRange:
            case CrossReference:
            case EndOfFile:
            case Group:
            case Keyword:
            case NegatedToken:
            case RegexToken:
            case RuleCall:
            case TerminalAlternatives:
            case TerminalGroup:
            case TerminalRuleCall:
            case UnorderedGroup:
            case UntilToken:
            case Wildcard: {
                return this.isSubtype(AbstractElement, supertype);
            }
            case ArrayLiteral:
            case NumberLiteral:
            case StringLiteral: {
                return this.isSubtype(ValueLiteral, supertype);
            }
            case ArrayType:
            case ReferenceType:
            case SimpleType:
            case UnionType: {
                return this.isSubtype(TypeDefinition, supertype);
            }
            case BooleanLiteral: {
                return this.isSubtype(Condition, supertype) || this.isSubtype(ValueLiteral, supertype);
            }
            case Conjunction:
            case Disjunction:
            case Negation:
            case ParameterReference: {
                return this.isSubtype(Condition, supertype);
            }
            case InferredType:
            case Interface:
            case Type: {
                return this.isSubtype(AbstractType, supertype);
            }
            case ParserRule: {
                return this.isSubtype(AbstractRule, supertype) || this.isSubtype(AbstractType, supertype);
            }
            case TerminalRule: {
                return this.isSubtype(AbstractRule, supertype);
            }
            default: {
                return false;
            }
        }
    }
    getReferenceType(refInfo) {
        const referenceId = `${refInfo.container.$type}:${refInfo.property}`;
        switch (referenceId) {
            case 'Action:type':
            case 'CrossReference:type':
            case 'Interface:superTypes':
            case 'ParserRule:returnType':
            case 'SimpleType:typeRef': {
                return AbstractType;
            }
            case 'Grammar:hiddenTokens':
            case 'ParserRule:hiddenTokens':
            case 'RuleCall:rule': {
                return AbstractRule;
            }
            case 'Grammar:usedGrammars': {
                return Grammar;
            }
            case 'NamedArgument:parameter':
            case 'ParameterReference:parameter': {
                return Parameter;
            }
            case 'TerminalRuleCall:rule': {
                return TerminalRule;
            }
            default: {
                throw new Error(`${referenceId} is not a valid reference id.`);
            }
        }
    }
    getTypeMetaData(type) {
        switch (type) {
            case AbstractElement: {
                return {
                    name: AbstractElement,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' }
                    ]
                };
            }
            case ArrayLiteral: {
                return {
                    name: ArrayLiteral,
                    properties: [
                        { name: 'elements', defaultValue: [] }
                    ]
                };
            }
            case ArrayType: {
                return {
                    name: ArrayType,
                    properties: [
                        { name: 'elementType' }
                    ]
                };
            }
            case BooleanLiteral: {
                return {
                    name: BooleanLiteral,
                    properties: [
                        { name: 'true', defaultValue: false }
                    ]
                };
            }
            case Conjunction: {
                return {
                    name: Conjunction,
                    properties: [
                        { name: 'left' },
                        { name: 'right' }
                    ]
                };
            }
            case Disjunction: {
                return {
                    name: Disjunction,
                    properties: [
                        { name: 'left' },
                        { name: 'right' }
                    ]
                };
            }
            case Grammar: {
                return {
                    name: Grammar,
                    properties: [
                        { name: 'definesHiddenTokens', defaultValue: false },
                        { name: 'hiddenTokens', defaultValue: [] },
                        { name: 'imports', defaultValue: [] },
                        { name: 'interfaces', defaultValue: [] },
                        { name: 'isDeclared', defaultValue: false },
                        { name: 'name' },
                        { name: 'rules', defaultValue: [] },
                        { name: 'types', defaultValue: [] },
                        { name: 'usedGrammars', defaultValue: [] }
                    ]
                };
            }
            case GrammarImport: {
                return {
                    name: GrammarImport,
                    properties: [
                        { name: 'path' }
                    ]
                };
            }
            case InferredType: {
                return {
                    name: InferredType,
                    properties: [
                        { name: 'name' }
                    ]
                };
            }
            case Interface: {
                return {
                    name: Interface,
                    properties: [
                        { name: 'attributes', defaultValue: [] },
                        { name: 'name' },
                        { name: 'superTypes', defaultValue: [] }
                    ]
                };
            }
            case NamedArgument: {
                return {
                    name: NamedArgument,
                    properties: [
                        { name: 'calledByName', defaultValue: false },
                        { name: 'parameter' },
                        { name: 'value' }
                    ]
                };
            }
            case Negation: {
                return {
                    name: Negation,
                    properties: [
                        { name: 'value' }
                    ]
                };
            }
            case NumberLiteral: {
                return {
                    name: NumberLiteral,
                    properties: [
                        { name: 'value' }
                    ]
                };
            }
            case Parameter: {
                return {
                    name: Parameter,
                    properties: [
                        { name: 'name' }
                    ]
                };
            }
            case ParameterReference: {
                return {
                    name: ParameterReference,
                    properties: [
                        { name: 'parameter' }
                    ]
                };
            }
            case ParserRule: {
                return {
                    name: ParserRule,
                    properties: [
                        { name: 'dataType' },
                        { name: 'definesHiddenTokens', defaultValue: false },
                        { name: 'definition' },
                        { name: 'entry', defaultValue: false },
                        { name: 'fragment', defaultValue: false },
                        { name: 'hiddenTokens', defaultValue: [] },
                        { name: 'inferredType' },
                        { name: 'name' },
                        { name: 'parameters', defaultValue: [] },
                        { name: 'returnType' },
                        { name: 'wildcard', defaultValue: false }
                    ]
                };
            }
            case ReferenceType: {
                return {
                    name: ReferenceType,
                    properties: [
                        { name: 'referenceType' }
                    ]
                };
            }
            case ReturnType: {
                return {
                    name: ReturnType,
                    properties: [
                        { name: 'name' }
                    ]
                };
            }
            case SimpleType: {
                return {
                    name: SimpleType,
                    properties: [
                        { name: 'primitiveType' },
                        { name: 'stringType' },
                        { name: 'typeRef' }
                    ]
                };
            }
            case StringLiteral: {
                return {
                    name: StringLiteral,
                    properties: [
                        { name: 'value' }
                    ]
                };
            }
            case TerminalRule: {
                return {
                    name: TerminalRule,
                    properties: [
                        { name: 'definition' },
                        { name: 'fragment', defaultValue: false },
                        { name: 'hidden', defaultValue: false },
                        { name: 'name' },
                        { name: 'type' }
                    ]
                };
            }
            case Type: {
                return {
                    name: Type,
                    properties: [
                        { name: 'name' },
                        { name: 'type' }
                    ]
                };
            }
            case TypeAttribute: {
                return {
                    name: TypeAttribute,
                    properties: [
                        { name: 'defaultValue' },
                        { name: 'isOptional', defaultValue: false },
                        { name: 'name' },
                        { name: 'type' }
                    ]
                };
            }
            case UnionType: {
                return {
                    name: UnionType,
                    properties: [
                        { name: 'types', defaultValue: [] }
                    ]
                };
            }
            case Action: {
                return {
                    name: Action,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'feature' },
                        { name: 'inferredType' },
                        { name: 'lookahead' },
                        { name: 'operator' },
                        { name: 'type' }
                    ]
                };
            }
            case Alternatives: {
                return {
                    name: Alternatives,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'elements', defaultValue: [] },
                        { name: 'lookahead' }
                    ]
                };
            }
            case Assignment: {
                return {
                    name: Assignment,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'feature' },
                        { name: 'lookahead' },
                        { name: 'operator' },
                        { name: 'terminal' }
                    ]
                };
            }
            case CharacterRange: {
                return {
                    name: CharacterRange,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'left' },
                        { name: 'lookahead' },
                        { name: 'right' }
                    ]
                };
            }
            case CrossReference: {
                return {
                    name: CrossReference,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'deprecatedSyntax', defaultValue: false },
                        { name: 'lookahead' },
                        { name: 'terminal' },
                        { name: 'type' }
                    ]
                };
            }
            case EndOfFile: {
                return {
                    name: EndOfFile,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' }
                    ]
                };
            }
            case Group: {
                return {
                    name: Group,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'elements', defaultValue: [] },
                        { name: 'guardCondition' },
                        { name: 'lookahead' }
                    ]
                };
            }
            case Keyword: {
                return {
                    name: Keyword,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'value' }
                    ]
                };
            }
            case NegatedToken: {
                return {
                    name: NegatedToken,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'terminal' }
                    ]
                };
            }
            case RegexToken: {
                return {
                    name: RegexToken,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'regex' }
                    ]
                };
            }
            case RuleCall: {
                return {
                    name: RuleCall,
                    properties: [
                        { name: 'arguments', defaultValue: [] },
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'rule' }
                    ]
                };
            }
            case TerminalAlternatives: {
                return {
                    name: TerminalAlternatives,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'elements', defaultValue: [] },
                        { name: 'lookahead' }
                    ]
                };
            }
            case TerminalGroup: {
                return {
                    name: TerminalGroup,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'elements', defaultValue: [] },
                        { name: 'lookahead' }
                    ]
                };
            }
            case TerminalRuleCall: {
                return {
                    name: TerminalRuleCall,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'rule' }
                    ]
                };
            }
            case UnorderedGroup: {
                return {
                    name: UnorderedGroup,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'elements', defaultValue: [] },
                        { name: 'lookahead' }
                    ]
                };
            }
            case UntilToken: {
                return {
                    name: UntilToken,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' },
                        { name: 'terminal' }
                    ]
                };
            }
            case Wildcard: {
                return {
                    name: Wildcard,
                    properties: [
                        { name: 'cardinality' },
                        { name: 'lookahead' }
                    ]
                };
            }
            default: {
                return {
                    name: type,
                    properties: []
                };
            }
        }
    }
}
const reflection = new LangiumGrammarAstReflection();
//# sourceMappingURL=ast.js.map

/***/ }),

/***/ 35481:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   P: () => (/* binding */ DefaultTokenBuilder)
/* harmony export */ });
/* harmony import */ var chevrotain__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(34326);
/* harmony import */ var _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(34905);
/* harmony import */ var _utils_ast_utils_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(74857);
/* harmony import */ var _utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(30447);
/* harmony import */ var _utils_regexp_utils_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(43078);
/* harmony import */ var _utils_stream_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(99293);
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






class DefaultTokenBuilder {
    constructor() {
        /**
         * The list of diagnostics stored during the lexing process of a single text.
         */
        this.diagnostics = [];
    }
    buildTokens(grammar, options) {
        const reachableRules = (0,_utils_stream_js__WEBPACK_IMPORTED_MODULE_1__/* .stream */ .Vw)((0,_utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_2__/* .getAllReachableRules */ .VD)(grammar, false));
        const terminalTokens = this.buildTerminalTokens(reachableRules);
        const tokens = this.buildKeywordTokens(reachableRules, terminalTokens, options);
        terminalTokens.forEach(terminalToken => {
            const pattern = terminalToken.PATTERN;
            if (typeof pattern === 'object' && pattern && 'test' in pattern && (0,_utils_regexp_utils_js__WEBPACK_IMPORTED_MODULE_3__/* .isWhitespace */ .cb)(pattern)) {
                tokens.unshift(terminalToken);
            }
            else {
                tokens.push(terminalToken);
            }
        });
        // We don't need to add the EOF token explicitly.
        // It is automatically available at the end of the token stream.
        return tokens;
    }
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    flushLexingReport(text) {
        return { diagnostics: this.popDiagnostics() };
    }
    popDiagnostics() {
        const diagnostics = [...this.diagnostics];
        this.diagnostics = [];
        return diagnostics;
    }
    buildTerminalTokens(rules) {
        return rules.filter(_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_4__/* .isTerminalRule */ .MS).filter(e => !e.fragment)
            .map(terminal => this.buildTerminalToken(terminal)).toArray();
    }
    buildTerminalToken(terminal) {
        const regex = (0,_utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_2__/* .terminalRegex */ .s1)(terminal);
        const pattern = this.requiresCustomPattern(regex) ? this.regexPatternFunction(regex) : regex;
        const tokenType = {
            name: terminal.name,
            PATTERN: pattern,
        };
        if (typeof pattern === 'function') {
            tokenType.LINE_BREAKS = true;
        }
        if (terminal.hidden) {
            // Only skip tokens that are able to accept whitespace
            tokenType.GROUP = (0,_utils_regexp_utils_js__WEBPACK_IMPORTED_MODULE_3__/* .isWhitespace */ .cb)(regex) ? chevrotain__WEBPACK_IMPORTED_MODULE_0__/* .Lexer */ .hW.SKIPPED : 'hidden';
        }
        return tokenType;
    }
    requiresCustomPattern(regex) {
        if (regex.flags.includes('u') || regex.flags.includes('s')) {
            // Unicode and dotall regexes are not supported by Chevrotain.
            return true;
        }
        else if (regex.source.includes('?<=') || regex.source.includes('?<!')) {
            // Negative and positive lookbehind are not supported by Chevrotain yet.
            return true;
        }
        else {
            return false;
        }
    }
    regexPatternFunction(regex) {
        const stickyRegex = new RegExp(regex, regex.flags + 'y');
        return (text, offset) => {
            stickyRegex.lastIndex = offset;
            const execResult = stickyRegex.exec(text);
            return execResult;
        };
    }
    buildKeywordTokens(rules, terminalTokens, options) {
        return rules
            // We filter by parser rules, since keywords in terminal rules get transformed into regex and are not actual tokens
            .filter(_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_4__/* .isParserRule */ .F9)
            .flatMap(rule => (0,_utils_ast_utils_js__WEBPACK_IMPORTED_MODULE_5__/* .streamAllContents */ .VY)(rule).filter(_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_4__/* .isKeyword */ .p1))
            .distinct(e => e.value).toArray()
            // Sort keywords by descending length
            .sort((a, b) => b.value.length - a.value.length)
            .map(keyword => this.buildKeywordToken(keyword, terminalTokens, Boolean(options === null || options === void 0 ? void 0 : options.caseInsensitive)));
    }
    buildKeywordToken(keyword, terminalTokens, caseInsensitive) {
        const keywordPattern = this.buildKeywordPattern(keyword, caseInsensitive);
        const tokenType = {
            name: keyword.value,
            PATTERN: keywordPattern,
            LONGER_ALT: this.findLongerAlt(keyword, terminalTokens)
        };
        if (typeof keywordPattern === 'function') {
            tokenType.LINE_BREAKS = true;
        }
        return tokenType;
    }
    buildKeywordPattern(keyword, caseInsensitive) {
        return caseInsensitive ?
            new RegExp((0,_utils_regexp_utils_js__WEBPACK_IMPORTED_MODULE_3__/* .getCaseInsensitivePattern */ .cp)(keyword.value)) :
            keyword.value;
    }
    findLongerAlt(keyword, terminalTokens) {
        return terminalTokens.reduce((longerAlts, token) => {
            const pattern = token === null || token === void 0 ? void 0 : token.PATTERN;
            if ((pattern === null || pattern === void 0 ? void 0 : pattern.source) && (0,_utils_regexp_utils_js__WEBPACK_IMPORTED_MODULE_3__/* .partialMatches */ .XC)('^' + pattern.source + '$', keyword.value)) {
                longerAlts.push(token);
            }
            return longerAlts;
        }, []);
    }
}
//# sourceMappingURL=token-builder.js.map

/***/ }),

/***/ 46174:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   t: () => (/* binding */ DefaultValueConverter)
/* harmony export */ });
/* unused harmony export ValueConverter */
/* harmony import */ var _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(34905);
/* harmony import */ var _utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(30447);
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


class DefaultValueConverter {
    convert(input, cstNode) {
        let feature = cstNode.grammarSource;
        if ((0,_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isCrossReference */ .Ki)(feature)) {
            feature = (0,_utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .getCrossReferenceTerminal */ .eN)(feature);
        }
        if ((0,_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isRuleCall */ .t3)(feature)) {
            const rule = feature.rule.ref;
            if (!rule) {
                throw new Error('This cst node was not parsed by a rule.');
            }
            return this.runConverter(rule, input, cstNode);
        }
        return input;
    }
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    runConverter(rule, input, cstNode) {
        var _a;
        switch (rule.name.toUpperCase()) {
            case 'INT': return ValueConverter.convertInt(input);
            case 'STRING': return ValueConverter.convertString(input);
            case 'ID': return ValueConverter.convertID(input);
        }
        switch ((_a = (0,_utils_grammar_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .getRuleType */ .mJ)(rule)) === null || _a === void 0 ? void 0 : _a.toLowerCase()) {
            case 'number': return ValueConverter.convertNumber(input);
            case 'boolean': return ValueConverter.convertBoolean(input);
            case 'bigint': return ValueConverter.convertBigint(input);
            case 'date': return ValueConverter.convertDate(input);
            default: return input;
        }
    }
}
var ValueConverter;
(function (ValueConverter) {
    function convertString(input) {
        let result = '';
        for (let i = 1; i < input.length - 1; i++) {
            const c = input.charAt(i);
            if (c === '\\') {
                const c1 = input.charAt(++i);
                result += convertEscapeCharacter(c1);
            }
            else {
                result += c;
            }
        }
        return result;
    }
    ValueConverter.convertString = convertString;
    function convertEscapeCharacter(char) {
        switch (char) {
            case 'b': return '\b';
            case 'f': return '\f';
            case 'n': return '\n';
            case 'r': return '\r';
            case 't': return '\t';
            case 'v': return '\v';
            case '0': return '\0';
            default: return char;
        }
    }
    function convertID(input) {
        if (input.charAt(0) === '^') {
            return input.substring(1);
        }
        else {
            return input;
        }
    }
    ValueConverter.convertID = convertID;
    function convertInt(input) {
        return parseInt(input);
    }
    ValueConverter.convertInt = convertInt;
    function convertBigint(input) {
        return BigInt(input);
    }
    ValueConverter.convertBigint = convertBigint;
    function convertDate(input) {
        return new Date(input);
    }
    ValueConverter.convertDate = convertDate;
    function convertNumber(input) {
        return Number(input);
    }
    ValueConverter.convertNumber = convertNumber;
    function convertBoolean(input) {
        return input.toLowerCase() === 'true';
    }
    ValueConverter.convertBoolean = convertBoolean;
})(ValueConverter || (ValueConverter = {}));
//# sourceMappingURL=value-converter.js.map

/***/ }),

/***/ 91303:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   $v: () => (/* binding */ AbstractAstReflection),
/* harmony export */   SI: () => (/* binding */ isAstNodeDescription),
/* harmony export */   U8: () => (/* binding */ isRootCstNode),
/* harmony export */   Yk: () => (/* binding */ isReference),
/* harmony export */   al: () => (/* binding */ isCompositeCstNode),
/* harmony export */   dm: () => (/* binding */ isLeafCstNode),
/* harmony export */   et: () => (/* binding */ isLinkingError),
/* harmony export */   xA: () => (/* binding */ isAstNode)
/* harmony export */ });
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
function isAstNode(obj) {
    return typeof obj === 'object' && obj !== null && typeof obj.$type === 'string';
}
function isReference(obj) {
    return typeof obj === 'object' && obj !== null && typeof obj.$refText === 'string';
}
function isAstNodeDescription(obj) {
    return typeof obj === 'object' && obj !== null
        && typeof obj.name === 'string'
        && typeof obj.type === 'string'
        && typeof obj.path === 'string';
}
function isLinkingError(obj) {
    return typeof obj === 'object' && obj !== null
        && isAstNode(obj.container)
        && isReference(obj.reference)
        && typeof obj.message === 'string';
}
/**
 * An abstract implementation of the {@link AstReflection} interface.
 * Serves to cache subtype computation results to improve performance throughout different parts of Langium.
 */
class AbstractAstReflection {
    constructor() {
        this.subtypes = {};
        this.allSubtypes = {};
    }
    isInstance(node, type) {
        return isAstNode(node) && this.isSubtype(node.$type, type);
    }
    isSubtype(subtype, supertype) {
        if (subtype === supertype) {
            return true;
        }
        let nested = this.subtypes[subtype];
        if (!nested) {
            nested = this.subtypes[subtype] = {};
        }
        const existing = nested[supertype];
        if (existing !== undefined) {
            return existing;
        }
        else {
            const result = this.computeIsSubtype(subtype, supertype);
            nested[supertype] = result;
            return result;
        }
    }
    getAllSubTypes(type) {
        const existing = this.allSubtypes[type];
        if (existing) {
            return existing;
        }
        else {
            const allTypes = this.getAllTypes();
            const types = [];
            for (const possibleSubType of allTypes) {
                if (this.isSubtype(possibleSubType, type)) {
                    types.push(possibleSubType);
                }
            }
            this.allSubtypes[type] = types;
            return types;
        }
    }
}
function isCompositeCstNode(node) {
    return typeof node === 'object' && node !== null && Array.isArray(node.content);
}
function isLeafCstNode(node) {
    return typeof node === 'object' && node !== null && typeof node.tokenType === 'object';
}
function isRootCstNode(node) {
    return isCompositeCstNode(node) && typeof node.fullText === 'string';
}
//# sourceMappingURL=syntax-tree.js.map

/***/ }),

/***/ 74857:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   E$: () => (/* binding */ findRootNode),
/* harmony export */   Me: () => (/* binding */ getDocument),
/* harmony export */   VY: () => (/* binding */ streamAllContents),
/* harmony export */   V_: () => (/* binding */ getContainerOfType),
/* harmony export */   Zc: () => (/* binding */ streamAst),
/* harmony export */   a1: () => (/* binding */ assignMandatoryProperties),
/* harmony export */   b2: () => (/* binding */ linkContentToContainer),
/* harmony export */   fy: () => (/* binding */ streamReferences),
/* harmony export */   sx: () => (/* binding */ streamContents)
/* harmony export */ });
/* unused harmony exports hasContainerOfType, findLocalReferences, copyAstNode */
/* harmony import */ var _syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(91303);
/* harmony import */ var _stream_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(99293);
/* harmony import */ var _cst_utils_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(13871);
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/



/**
 * Link the `$container` and other related properties of every AST node that is directly contained
 * in the given `node`.
 */
function linkContentToContainer(node) {
    for (const [name, value] of Object.entries(node)) {
        if (!name.startsWith('$')) {
            if (Array.isArray(value)) {
                value.forEach((item, index) => {
                    if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isAstNode */ .xA)(item)) {
                        item.$container = node;
                        item.$containerProperty = name;
                        item.$containerIndex = index;
                    }
                });
            }
            else if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isAstNode */ .xA)(value)) {
                value.$container = node;
                value.$containerProperty = name;
            }
        }
    }
}
/**
 * Walk along the hierarchy of containers from the given AST node to the root and return the first
 * node that matches the type predicate. If the start node itself matches, it is returned.
 * If no container matches, `undefined` is returned.
 */
function getContainerOfType(node, typePredicate) {
    let item = node;
    while (item) {
        if (typePredicate(item)) {
            return item;
        }
        item = item.$container;
    }
    return undefined;
}
/**
 * Walk along the hierarchy of containers from the given AST node to the root and check for existence
 * of a container that matches the given predicate. The start node is included in the checks.
 */
function hasContainerOfType(node, predicate) {
    let item = node;
    while (item) {
        if (predicate(item)) {
            return true;
        }
        item = item.$container;
    }
    return false;
}
/**
 * Retrieve the document in which the given AST node is contained. A reference to the document is
 * usually held by the root node of the AST.
 *
 * @throws an error if the node is not contained in a document.
 */
function getDocument(node) {
    const rootNode = findRootNode(node);
    const result = rootNode.$document;
    if (!result) {
        throw new Error('AST node has no document.');
    }
    return result;
}
/**
 * Returns the root node of the given AST node by following the `$container` references.
 */
function findRootNode(node) {
    while (node.$container) {
        node = node.$container;
    }
    return node;
}
/**
 * Create a stream of all AST nodes that are directly contained in the given node. This includes
 * single-valued as well as multi-valued (array) properties.
 */
function streamContents(node, options) {
    if (!node) {
        throw new Error('Node must be an AstNode.');
    }
    const range = options === null || options === void 0 ? void 0 : options.range;
    return new _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .StreamImpl */ .i(() => ({
        keys: Object.keys(node),
        keyIndex: 0,
        arrayIndex: 0
    }), state => {
        while (state.keyIndex < state.keys.length) {
            const property = state.keys[state.keyIndex];
            if (!property.startsWith('$')) {
                const value = node[property];
                if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isAstNode */ .xA)(value)) {
                    state.keyIndex++;
                    if (isAstNodeInRange(value, range)) {
                        return { done: false, value };
                    }
                }
                else if (Array.isArray(value)) {
                    while (state.arrayIndex < value.length) {
                        const index = state.arrayIndex++;
                        const element = value[index];
                        if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isAstNode */ .xA)(element) && isAstNodeInRange(element, range)) {
                            return { done: false, value: element };
                        }
                    }
                    state.arrayIndex = 0;
                }
            }
            state.keyIndex++;
        }
        return _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .DONE_RESULT */ .Ry;
    });
}
/**
 * Create a stream of all AST nodes that are directly and indirectly contained in the given root node.
 * This does not include the root node itself.
 */
function streamAllContents(root, options) {
    if (!root) {
        throw new Error('Root node must be an AstNode.');
    }
    return new _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .TreeStreamImpl */ .i8(root, node => streamContents(node, options));
}
/**
 * Create a stream of all AST nodes that are directly and indirectly contained in the given root node,
 * including the root node itself.
 */
function streamAst(root, options) {
    if (!root) {
        throw new Error('Root node must be an AstNode.');
    }
    else if ((options === null || options === void 0 ? void 0 : options.range) && !isAstNodeInRange(root, options.range)) {
        // Return an empty stream if the root node isn't in range
        return new _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .TreeStreamImpl */ .i8(root, () => []);
    }
    return new _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .TreeStreamImpl */ .i8(root, node => streamContents(node, options), { includeRoot: true });
}
function isAstNodeInRange(astNode, range) {
    var _a;
    if (!range) {
        return true;
    }
    const nodeRange = (_a = astNode.$cstNode) === null || _a === void 0 ? void 0 : _a.range;
    if (!nodeRange) {
        return false;
    }
    return (0,_cst_utils_js__WEBPACK_IMPORTED_MODULE_2__/* .inRange */ .Z2)(nodeRange, range);
}
/**
 * Create a stream of all cross-references that are held by the given AST node. This includes
 * single-valued as well as multi-valued (array) properties.
 */
function streamReferences(node) {
    return new _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .StreamImpl */ .i(() => ({
        keys: Object.keys(node),
        keyIndex: 0,
        arrayIndex: 0
    }), state => {
        while (state.keyIndex < state.keys.length) {
            const property = state.keys[state.keyIndex];
            if (!property.startsWith('$')) {
                const value = node[property];
                if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isReference */ .Yk)(value)) {
                    state.keyIndex++;
                    return { done: false, value: { reference: value, container: node, property } };
                }
                else if (Array.isArray(value)) {
                    while (state.arrayIndex < value.length) {
                        const index = state.arrayIndex++;
                        const element = value[index];
                        if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_0__/* .isReference */ .Yk)(element)) {
                            return { done: false, value: { reference: element, container: node, property, index } };
                        }
                    }
                    state.arrayIndex = 0;
                }
            }
            state.keyIndex++;
        }
        return _stream_js__WEBPACK_IMPORTED_MODULE_1__/* .DONE_RESULT */ .Ry;
    });
}
/**
 * Returns a Stream of references to the target node from the AstNode tree
 *
 * @param targetNode AstNode we are looking for
 * @param lookup AstNode where we search for references. If not provided, the root node of the document is used as the default value
 */
function findLocalReferences(targetNode, lookup = getDocument(targetNode).parseResult.value) {
    const refs = [];
    streamAst(lookup).forEach(node => {
        streamReferences(node).forEach(refInfo => {
            if (refInfo.reference.ref === targetNode) {
                refs.push(refInfo.reference);
            }
        });
    });
    return stream(refs);
}
/**
 * Assigns all mandatory AST properties to the specified node.
 *
 * @param reflection Reflection object used to gather mandatory properties for the node.
 * @param node Specified node is modified in place and properties are directly assigned.
 */
function assignMandatoryProperties(reflection, node) {
    const typeMetaData = reflection.getTypeMetaData(node.$type);
    const genericNode = node;
    for (const property of typeMetaData.properties) {
        // Only set the value if the property is not already set and if it has a default value
        if (property.defaultValue !== undefined && genericNode[property.name] === undefined) {
            genericNode[property.name] = copyDefaultValue(property.defaultValue);
        }
    }
}
function copyDefaultValue(propertyType) {
    if (Array.isArray(propertyType)) {
        return [...propertyType.map(copyDefaultValue)];
    }
    else {
        return propertyType;
    }
}
/**
 * Creates a deep copy of the specified AST node.
 * The resulting copy will only contain semantically relevant information, such as the `$type` property and AST properties.
 *
 * References are copied without resolved cross reference. The specified function is used to rebuild them.
 */
function copyAstNode(node, buildReference) {
    const copy = { $type: node.$type };
    for (const [name, value] of Object.entries(node)) {
        if (!name.startsWith('$')) {
            if (isAstNode(value)) {
                copy[name] = copyAstNode(value, buildReference);
            }
            else if (isReference(value)) {
                copy[name] = buildReference(copy, name, value.$refNode, value.$refText);
            }
            else if (Array.isArray(value)) {
                const copiedArray = [];
                for (const element of value) {
                    if (isAstNode(element)) {
                        copiedArray.push(copyAstNode(element, buildReference));
                    }
                    else if (isReference(element)) {
                        copiedArray.push(buildReference(copy, name, element.$refNode, element.$refText));
                    }
                    else {
                        copiedArray.push(element);
                    }
                }
                copy[name] = copiedArray;
            }
            else {
                copy[name] = value;
            }
        }
    }
    linkContentToContainer(copy);
    return copy;
}
//# sourceMappingURL=ast-utils.js.map

/***/ }),

/***/ 13871:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   LK: () => (/* binding */ findCommentNode),
/* harmony export */   OB: () => (/* binding */ isChildNode),
/* harmony export */   Z2: () => (/* binding */ inRange),
/* harmony export */   _t: () => (/* binding */ streamCst),
/* harmony export */   sp: () => (/* binding */ tokenToRange),
/* harmony export */   uz: () => (/* binding */ DefaultNameRegexp),
/* harmony export */   yn: () => (/* binding */ toDocumentSegment)
/* harmony export */ });
/* unused harmony exports flattenCst, RangeComparison, compareRange, findDeclarationNodeAtOffset, isCommentNode, findLeafNodeAtOffset, findLeafNodeBeforeOffset, getPreviousNode, getNextNode, getStartlineNode, getInteriorNodes */
/* harmony import */ var _syntax_tree_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(91303);
/* harmony import */ var _stream_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(99293);
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/


/**
 * Create a stream of all CST nodes that are directly and indirectly contained in the given root node,
 * including the root node itself.
 */
function streamCst(node) {
    return new _stream_js__WEBPACK_IMPORTED_MODULE_0__/* .TreeStreamImpl */ .i8(node, element => {
        if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_1__/* .isCompositeCstNode */ .al)(element)) {
            return element.content;
        }
        else {
            return [];
        }
    }, { includeRoot: true });
}
/**
 * Create a stream of all leaf nodes that are directly and indirectly contained in the given root node.
 */
function flattenCst(node) {
    return streamCst(node).filter(isLeafCstNode);
}
/**
 * Determines whether the specified cst node is a child of the specified parent node.
 */
function isChildNode(child, parent) {
    while (child.container) {
        child = child.container;
        if (child === parent) {
            return true;
        }
    }
    return false;
}
function tokenToRange(token) {
    // Chevrotain uses 1-based indices everywhere
    // So we subtract 1 from every value to align with the LSP
    return {
        start: {
            character: token.startColumn - 1,
            line: token.startLine - 1
        },
        end: {
            character: token.endColumn, // endColumn uses the correct index
            line: token.endLine - 1
        }
    };
}
function toDocumentSegment(node) {
    if (!node) {
        return undefined;
    }
    const { offset, end, range } = node;
    return {
        range,
        offset,
        end,
        length: end - offset
    };
}
var RangeComparison;
(function (RangeComparison) {
    RangeComparison[RangeComparison["Before"] = 0] = "Before";
    RangeComparison[RangeComparison["After"] = 1] = "After";
    RangeComparison[RangeComparison["OverlapFront"] = 2] = "OverlapFront";
    RangeComparison[RangeComparison["OverlapBack"] = 3] = "OverlapBack";
    RangeComparison[RangeComparison["Inside"] = 4] = "Inside";
    RangeComparison[RangeComparison["Outside"] = 5] = "Outside";
})(RangeComparison || (RangeComparison = {}));
function compareRange(range, to) {
    if (range.end.line < to.start.line || (range.end.line === to.start.line && range.end.character <= to.start.character)) {
        return RangeComparison.Before;
    }
    else if (range.start.line > to.end.line || (range.start.line === to.end.line && range.start.character >= to.end.character)) {
        return RangeComparison.After;
    }
    const startInside = range.start.line > to.start.line || (range.start.line === to.start.line && range.start.character >= to.start.character);
    const endInside = range.end.line < to.end.line || (range.end.line === to.end.line && range.end.character <= to.end.character);
    if (startInside && endInside) {
        return RangeComparison.Inside;
    }
    else if (startInside) {
        return RangeComparison.OverlapBack;
    }
    else if (endInside) {
        return RangeComparison.OverlapFront;
    }
    else {
        return RangeComparison.Outside;
    }
}
function inRange(range, to) {
    const comparison = compareRange(range, to);
    return comparison > RangeComparison.After;
}
// The \p{L} regex matches any unicode letter character, i.e. characters from non-english alphabets
// Together with \w it matches any kind of character which can commonly appear in IDs
const DefaultNameRegexp = /^[\w\p{L}]$/u;
/**
 * Performs `findLeafNodeAtOffset` with a minor difference: When encountering a character that matches the `nameRegexp` argument,
 * it will instead return the leaf node at the `offset - 1` position.
 *
 * For LSP services, users expect that the declaration of an element is available if the cursor is directly after the element.
 */
function findDeclarationNodeAtOffset(cstNode, offset, nameRegexp = DefaultNameRegexp) {
    if (cstNode) {
        if (offset > 0) {
            const localOffset = offset - cstNode.offset;
            const textAtOffset = cstNode.text.charAt(localOffset);
            if (!nameRegexp.test(textAtOffset)) {
                offset--;
            }
        }
        return findLeafNodeAtOffset(cstNode, offset);
    }
    return undefined;
}
function findCommentNode(cstNode, commentNames) {
    if (cstNode) {
        const previous = getPreviousNode(cstNode, true);
        if (previous && isCommentNode(previous, commentNames)) {
            return previous;
        }
        if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_1__/* .isRootCstNode */ .U8)(cstNode)) {
            // Go from the first non-hidden node through all nodes in reverse order
            // We do this to find the comment node which directly precedes the root node
            const endIndex = cstNode.content.findIndex(e => !e.hidden);
            for (let i = endIndex - 1; i >= 0; i--) {
                const child = cstNode.content[i];
                if (isCommentNode(child, commentNames)) {
                    return child;
                }
            }
        }
    }
    return undefined;
}
function isCommentNode(cstNode, commentNames) {
    return (0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_1__/* .isLeafCstNode */ .dm)(cstNode) && commentNames.includes(cstNode.tokenType.name);
}
/**
 * Finds the leaf CST node at the specified 0-based string offset.
 * Note that the given offset will be within the range of the returned leaf node.
 *
 * If the offset does not point to a CST node (but just white space), this method will return `undefined`.
 *
 * @param node The CST node to search through.
 * @param offset The specified offset.
 * @returns The CST node at the specified offset.
 */
function findLeafNodeAtOffset(node, offset) {
    if (isLeafCstNode(node)) {
        return node;
    }
    else if (isCompositeCstNode(node)) {
        const searchResult = binarySearch(node, offset, false);
        if (searchResult) {
            return findLeafNodeAtOffset(searchResult, offset);
        }
    }
    return undefined;
}
/**
 * Finds the leaf CST node at the specified 0-based string offset.
 * If no CST node exists at the specified position, it will return the leaf node before it.
 *
 * If there is no leaf node before the specified offset, this method will return `undefined`.
 *
 * @param node The CST node to search through.
 * @param offset The specified offset.
 * @returns The CST node closest to the specified offset.
 */
function findLeafNodeBeforeOffset(node, offset) {
    if (isLeafCstNode(node)) {
        return node;
    }
    else if (isCompositeCstNode(node)) {
        const searchResult = binarySearch(node, offset, true);
        if (searchResult) {
            return findLeafNodeBeforeOffset(searchResult, offset);
        }
    }
    return undefined;
}
function binarySearch(node, offset, closest) {
    let left = 0;
    let right = node.content.length - 1;
    let closestNode = undefined;
    while (left <= right) {
        const middle = Math.floor((left + right) / 2);
        const middleNode = node.content[middle];
        if (middleNode.offset <= offset && middleNode.end > offset) {
            // Found an exact match
            return middleNode;
        }
        if (middleNode.end <= offset) {
            // Update the closest node (less than offset) and move to the right half
            closestNode = closest ? middleNode : undefined;
            left = middle + 1;
        }
        else {
            // Move to the left half
            right = middle - 1;
        }
    }
    return closestNode;
}
function getPreviousNode(node, hidden = true) {
    while (node.container) {
        const parent = node.container;
        let index = parent.content.indexOf(node);
        while (index > 0) {
            index--;
            const previous = parent.content[index];
            if (hidden || !previous.hidden) {
                return previous;
            }
        }
        node = parent;
    }
    return undefined;
}
function getNextNode(node, hidden = true) {
    while (node.container) {
        const parent = node.container;
        let index = parent.content.indexOf(node);
        const last = parent.content.length - 1;
        while (index < last) {
            index++;
            const next = parent.content[index];
            if (hidden || !next.hidden) {
                return next;
            }
        }
        node = parent;
    }
    return undefined;
}
function getStartlineNode(node) {
    if (node.range.start.character === 0) {
        return node;
    }
    const line = node.range.start.line;
    let last = node;
    let index;
    while (node.container) {
        const parent = node.container;
        const selfIndex = index !== null && index !== void 0 ? index : parent.content.indexOf(node);
        if (selfIndex === 0) {
            node = parent;
            index = undefined;
        }
        else {
            index = selfIndex - 1;
            node = parent.content[index];
        }
        if (node.range.start.line !== line) {
            break;
        }
        last = node;
    }
    return last;
}
function getInteriorNodes(start, end) {
    const commonParent = getCommonParent(start, end);
    if (!commonParent) {
        return [];
    }
    return commonParent.parent.content.slice(commonParent.a + 1, commonParent.b);
}
function getCommonParent(a, b) {
    const aParents = getParentChain(a);
    const bParents = getParentChain(b);
    let current;
    for (let i = 0; i < aParents.length && i < bParents.length; i++) {
        const aParent = aParents[i];
        const bParent = bParents[i];
        if (aParent.parent === bParent.parent) {
            current = {
                parent: aParent.parent,
                a: aParent.index,
                b: bParent.index
            };
        }
        else {
            break;
        }
    }
    return current;
}
function getParentChain(node) {
    const chain = [];
    while (node.container) {
        const parent = node.container;
        const index = parent.content.indexOf(node);
        chain.push({
            parent,
            index
        });
        node = parent;
    }
    return chain.reverse();
}
//# sourceMappingURL=cst-utils.js.map

/***/ }),

/***/ 45209:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   U: () => (/* binding */ assertUnreachable),
/* harmony export */   h: () => (/* binding */ ErrorWithLocation)
/* harmony export */ });
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
class ErrorWithLocation extends Error {
    constructor(node, message) {
        super(node ? `${message} at ${node.range.start.line}:${node.range.start.character}` : message);
    }
}
function assertUnreachable(_) {
    throw new Error('Error! The input value was not handled.');
}
//# sourceMappingURL=errors.js.map

/***/ }),

/***/ 30447:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   $G: () => (/* binding */ getExplicitRuleType),
/* harmony export */   EL: () => (/* binding */ findNodesForProperty),
/* harmony export */   UP: () => (/* binding */ isDataTypeRule),
/* harmony export */   VD: () => (/* binding */ getAllReachableRules),
/* harmony export */   eN: () => (/* binding */ getCrossReferenceTerminal),
/* harmony export */   h7: () => (/* binding */ findAssignment),
/* harmony export */   ib: () => (/* binding */ findNameAssignment),
/* harmony export */   lA: () => (/* binding */ findNodeForKeyword),
/* harmony export */   mJ: () => (/* binding */ getRuleType),
/* harmony export */   md: () => (/* binding */ isCommentTerminal),
/* harmony export */   s1: () => (/* binding */ terminalRegex),
/* harmony export */   vb: () => (/* binding */ findNodeForProperty),
/* harmony export */   z$: () => (/* binding */ getTypeName)
/* harmony export */ });
/* unused harmony exports getEntryRule, getHiddenRules, findNodesForKeyword, findNodesForKeywordInternal, getActionAtElement, isOptionalCardinality, isArrayCardinality, isArrayOperator, isDataType, getActionType, getRuleTypeName */
/* harmony import */ var _utils_errors_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(45209);
/* harmony import */ var _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(34905);
/* harmony import */ var _syntax_tree_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(91303);
/* harmony import */ var _ast_utils_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74857);
/* harmony import */ var _cst_utils_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(13871);
/* harmony import */ var _regexp_utils_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(43078);
/******************************************************************************
 * Copyright 2021-2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/






/**
 * Returns the entry rule of the given grammar, if any. If the grammar file does not contain an entry rule,
 * the result is `undefined`.
 */
function getEntryRule(grammar) {
    return grammar.rules.find(e => _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(e) && e.entry);
}
/**
 * Returns all hidden terminal rules of the given grammar, if any.
 */
function getHiddenRules(grammar) {
    return grammar.rules.filter((e) => _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalRule */ .MS(e) && e.hidden);
}
/**
 * Returns all rules that can be reached from the topmost rules of the specified grammar (entry and hidden terminal rules).
 *
 * @param grammar The grammar that contains all rules
 * @param allTerminals Whether or not to include terminals that are referenced only by other terminals
 * @returns A list of referenced parser and terminal rules. If the grammar contains no entry rule,
 *      this function returns all rules of the specified grammar.
 */
function getAllReachableRules(grammar, allTerminals) {
    const ruleNames = new Set();
    const entryRule = getEntryRule(grammar);
    if (!entryRule) {
        return new Set(grammar.rules);
    }
    const topMostRules = [entryRule].concat(getHiddenRules(grammar));
    for (const rule of topMostRules) {
        ruleDfs(rule, ruleNames, allTerminals);
    }
    const rules = new Set();
    for (const rule of grammar.rules) {
        if (ruleNames.has(rule.name) || (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalRule */ .MS(rule) && rule.hidden)) {
            rules.add(rule);
        }
    }
    return rules;
}
function ruleDfs(rule, visitedSet, allTerminals) {
    visitedSet.add(rule.name);
    (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .streamAllContents */ .VY)(rule).forEach(node => {
        if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isRuleCall */ .t3(node) || (allTerminals && _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalRuleCall */ .gf(node))) {
            const refRule = node.rule.ref;
            if (refRule && !visitedSet.has(refRule.name)) {
                ruleDfs(refRule, visitedSet, allTerminals);
            }
        }
    });
}
/**
 * Determines the grammar expression used to parse a cross-reference (usually a reference to a terminal rule).
 * A cross-reference can declare this expression explicitly in the form `[Type : Terminal]`, but if `Terminal`
 * is omitted, this function attempts to infer it from the name of the referenced `Type` (using `findNameAssignment`).
 *
 * Returns the grammar expression used to parse the given cross-reference, or `undefined` if it is not declared
 * and cannot be inferred.
 */
function getCrossReferenceTerminal(crossRef) {
    if (crossRef.terminal) {
        return crossRef.terminal;
    }
    else if (crossRef.type.ref) {
        const nameAssigment = findNameAssignment(crossRef.type.ref);
        return nameAssigment === null || nameAssigment === void 0 ? void 0 : nameAssigment.terminal;
    }
    return undefined;
}
/**
 * Determines whether the given terminal rule represents a comment. This is true if the rule is marked
 * as `hidden` and it does not match white space. This means every hidden token (i.e. excluded from the AST)
 * that contains visible characters is considered a comment.
 */
function isCommentTerminal(terminalRule) {
    return terminalRule.hidden && !(0,_regexp_utils_js__WEBPACK_IMPORTED_MODULE_2__/* .isWhitespace */ .cb)(terminalRegex(terminalRule));
}
/**
 * Find all CST nodes within the given node that contribute to the specified property.
 *
 * @param node A CST node in which to look for property assignments. If this is undefined, the result is an empty array.
 * @param property A property name of the constructed AST node. If this is undefined, the result is an empty array.
 */
function findNodesForProperty(node, property) {
    if (!node || !property) {
        return [];
    }
    return findNodesForPropertyInternal(node, property, node.astNode, true);
}
/**
 * Find a single CST node within the given node that contributes to the specified property.
 *
 * @param node A CST node in which to look for property assignments. If this is undefined, the result is `undefined`.
 * @param property A property name of the constructed AST node. If this is undefined, the result is `undefined`.
 * @param index If no index is specified or the index is less than zero, the first found node is returned. If the
 *        specified index exceeds the number of assignments to the property, the last found node is returned. Otherwise,
 *        the node with the specified index is returned.
 */
function findNodeForProperty(node, property, index) {
    if (!node || !property) {
        return undefined;
    }
    const nodes = findNodesForPropertyInternal(node, property, node.astNode, true);
    if (nodes.length === 0) {
        return undefined;
    }
    if (index !== undefined) {
        index = Math.max(0, Math.min(index, nodes.length - 1));
    }
    else {
        index = 0;
    }
    return nodes[index];
}
function findNodesForPropertyInternal(node, property, element, first) {
    if (!first) {
        const nodeFeature = (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .getContainerOfType */ .V_)(node.grammarSource, _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAssignment */ .B7);
        if (nodeFeature && nodeFeature.feature === property) {
            return [node];
        }
    }
    if ((0,_syntax_tree_js__WEBPACK_IMPORTED_MODULE_3__/* .isCompositeCstNode */ .al)(node) && node.astNode === element) {
        return node.content.flatMap(e => findNodesForPropertyInternal(e, property, element, false));
    }
    return [];
}
/**
 * Find all CST nodes within the given node that correspond to the specified keyword.
 *
 * @param node A CST node in which to look for keywords. If this is undefined, the result is an empty array.
 * @param keyword A keyword as specified in the grammar.
 */
function findNodesForKeyword(node, keyword) {
    if (!node) {
        return [];
    }
    return findNodesForKeywordInternal(node, keyword, node === null || node === void 0 ? void 0 : node.astNode);
}
/**
 * Find a single CST node within the given node that corresponds to the specified keyword.
 *
 * @param node A CST node in which to look for keywords. If this is undefined, the result is `undefined`.
 * @param keyword A keyword as specified in the grammar.
 * @param index If no index is specified or the index is less than zero, the first found node is returned. If the
 *        specified index exceeds the number of keyword occurrences, the last found node is returned. Otherwise,
 *        the node with the specified index is returned.
 */
function findNodeForKeyword(node, keyword, index) {
    if (!node) {
        return undefined;
    }
    const nodes = findNodesForKeywordInternal(node, keyword, node === null || node === void 0 ? void 0 : node.astNode);
    if (nodes.length === 0) {
        return undefined;
    }
    if (index !== undefined) {
        index = Math.max(0, Math.min(index, nodes.length - 1));
    }
    else {
        index = 0;
    }
    return nodes[index];
}
function findNodesForKeywordInternal(node, keyword, element) {
    if (node.astNode !== element) {
        return [];
    }
    if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isKeyword */ .p1(node.grammarSource) && node.grammarSource.value === keyword) {
        return [node];
    }
    const treeIterator = (0,_cst_utils_js__WEBPACK_IMPORTED_MODULE_4__/* .streamCst */ ._t)(node).iterator();
    let result;
    const keywordNodes = [];
    do {
        result = treeIterator.next();
        if (!result.done) {
            const childNode = result.value;
            if (childNode.astNode === element) {
                if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isKeyword */ .p1(childNode.grammarSource) && childNode.grammarSource.value === keyword) {
                    keywordNodes.push(childNode);
                }
            }
            else {
                treeIterator.prune();
            }
        }
    } while (!result.done);
    return keywordNodes;
}
/**
 * If the given CST node was parsed in the context of a property assignment, the respective `Assignment` grammar
 * node is returned. If no assignment is found, the result is `undefined`.
 *
 * @param cstNode A CST node for which to find a property assignment.
 */
function findAssignment(cstNode) {
    var _a;
    const astNode = cstNode.astNode;
    // Only search until the ast node of the parent cst node is no longer the original ast node
    // This would make us jump to a preceding rule call, which contains only unrelated assignments
    while (astNode === ((_a = cstNode.container) === null || _a === void 0 ? void 0 : _a.astNode)) {
        const assignment = (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .getContainerOfType */ .V_)(cstNode.grammarSource, _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAssignment */ .B7);
        if (assignment) {
            return assignment;
        }
        cstNode = cstNode.container;
    }
    return undefined;
}
/**
 * Find an assignment to the `name` property for the given grammar type. This requires the `type` to be inferred
 * from a parser rule, and that rule must contain an assignment to the `name` property. In all other cases,
 * this function returns `undefined`.
 */
function findNameAssignment(type) {
    let startNode = type;
    if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isInferredType */ .S_(startNode)) {
        // for inferred types, the location to start searching for the name-assignment is different
        if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAction */ .LG(startNode.$container)) {
            // a type which is explicitly inferred by an action: investigate the sibbling of the Action node, i.e. start searching at the Action's parent
            startNode = startNode.$container.$container;
        }
        else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(startNode.$container)) {
            // investigate the parser rule with the explicitly inferred type
            startNode = startNode.$container;
        }
        else {
            (0,_utils_errors_js__WEBPACK_IMPORTED_MODULE_5__/* .assertUnreachable */ .U)(startNode.$container);
        }
    }
    return findNameAssignmentInternal(type, startNode, new Map());
}
function findNameAssignmentInternal(type, startNode, cache) {
    var _a;
    // the cache is only required to prevent infinite loops
    function go(node, refType) {
        let childAssignment = undefined;
        const parentAssignment = (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .getContainerOfType */ .V_)(node, _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAssignment */ .B7);
        // No parent assignment implies unassigned rule call
        if (!parentAssignment) {
            childAssignment = findNameAssignmentInternal(refType, refType, cache);
        }
        cache.set(type, childAssignment);
        return childAssignment;
    }
    if (cache.has(type)) {
        return cache.get(type);
    }
    cache.set(type, undefined);
    for (const node of (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .streamAllContents */ .VY)(startNode)) {
        if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAssignment */ .B7(node) && node.feature.toLowerCase() === 'name') {
            cache.set(type, node);
            return node;
        }
        else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isRuleCall */ .t3(node) && _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(node.rule.ref)) {
            return go(node, node.rule.ref);
        }
        else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isSimpleType */ .Iy(node) && ((_a = node.typeRef) === null || _a === void 0 ? void 0 : _a.ref)) {
            return go(node, node.typeRef.ref);
        }
    }
    return undefined;
}
function getActionAtElement(element) {
    const parent = element.$container;
    if (ast.isGroup(parent)) {
        const elements = parent.elements;
        const index = elements.indexOf(element);
        for (let i = index - 1; i >= 0; i--) {
            const item = elements[i];
            if (ast.isAction(item)) {
                return item;
            }
            else {
                const action = streamAllContents(elements[i]).find(ast.isAction);
                if (action) {
                    return action;
                }
            }
        }
    }
    if (ast.isAbstractElement(parent)) {
        return getActionAtElement(parent);
    }
    else {
        return undefined;
    }
}
function isOptionalCardinality(cardinality, element) {
    return cardinality === '?' || cardinality === '*' || (ast.isGroup(element) && Boolean(element.guardCondition));
}
function isArrayCardinality(cardinality) {
    return cardinality === '*' || cardinality === '+';
}
function isArrayOperator(operator) {
    return operator === '+=';
}
/**
 * Determines whether the given parser rule is a _data type rule_, meaning that it has a
 * primitive return type like `number`, `boolean`, etc.
 */
function isDataTypeRule(rule) {
    return isDataTypeRuleInternal(rule, new Set());
}
function isDataTypeRuleInternal(rule, visited) {
    if (visited.has(rule)) {
        return true;
    }
    else {
        visited.add(rule);
    }
    for (const node of (0,_ast_utils_js__WEBPACK_IMPORTED_MODULE_1__/* .streamAllContents */ .VY)(rule)) {
        if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isRuleCall */ .t3(node)) {
            if (!node.rule.ref) {
                // RuleCall to unresolved rule. Don't assume `rule` is a DataType rule.
                return false;
            }
            if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(node.rule.ref) && !isDataTypeRuleInternal(node.rule.ref, visited)) {
                return false;
            }
        }
        else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAssignment */ .B7(node)) {
            return false;
        }
        else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAction */ .LG(node)) {
            return false;
        }
    }
    return Boolean(rule.definition);
}
function isDataType(type) {
    return isDataTypeInternal(type.type, new Set());
}
function isDataTypeInternal(type, visited) {
    if (visited.has(type)) {
        return true;
    }
    else {
        visited.add(type);
    }
    if (ast.isArrayType(type)) {
        return false;
    }
    else if (ast.isReferenceType(type)) {
        return false;
    }
    else if (ast.isUnionType(type)) {
        return type.types.every(e => isDataTypeInternal(e, visited));
    }
    else if (ast.isSimpleType(type)) {
        if (type.primitiveType !== undefined) {
            return true;
        }
        else if (type.stringType !== undefined) {
            return true;
        }
        else if (type.typeRef !== undefined) {
            const ref = type.typeRef.ref;
            if (ast.isType(ref)) {
                return isDataTypeInternal(ref.type, visited);
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}
function getExplicitRuleType(rule) {
    if (rule.inferredType) {
        return rule.inferredType.name;
    }
    else if (rule.dataType) {
        return rule.dataType;
    }
    else if (rule.returnType) {
        const refType = rule.returnType.ref;
        if (refType) {
            // check if we need to check Action as return type
            if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(refType)) {
                return refType.name;
            }
            else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isInterface */ .QV(refType) || _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isType */ .P9(refType)) {
                return refType.name;
            }
        }
    }
    return undefined;
}
function getTypeName(type) {
    var _a;
    if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isParserRule */ .F9(type)) {
        return isDataTypeRule(type) ? type.name : (_a = getExplicitRuleType(type)) !== null && _a !== void 0 ? _a : type.name;
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isInterface */ .QV(type) || _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isType */ .P9(type) || _languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isReturnType */ .Mp(type)) {
        return type.name;
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isAction */ .LG(type)) {
        const actionType = getActionType(type);
        if (actionType) {
            return actionType;
        }
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isInferredType */ .S_(type)) {
        return type.name;
    }
    throw new Error('Cannot get name of Unknown Type');
}
function getActionType(action) {
    var _a;
    if (action.inferredType) {
        return action.inferredType.name;
    }
    else if ((_a = action.type) === null || _a === void 0 ? void 0 : _a.ref) {
        return getTypeName(action.type.ref);
    }
    return undefined; // not inferring and not referencing a valid type
}
/**
 * This function is used at development time (for code generation and the internal type system) to get the type of the AST node produced by the given rule.
 * For data type rules, the name of the rule is returned,
 * e.g. "INT_value returns number: MY_INT;" returns "INT_value".
 * @param rule the given rule
 * @returns the name of the AST node type of the rule
 */
function getRuleTypeName(rule) {
    var _a, _b, _c;
    if (ast.isTerminalRule(rule)) {
        return (_b = (_a = rule.type) === null || _a === void 0 ? void 0 : _a.name) !== null && _b !== void 0 ? _b : 'string';
    }
    else {
        return isDataTypeRule(rule) ? rule.name : (_c = getExplicitRuleType(rule)) !== null && _c !== void 0 ? _c : rule.name;
    }
}
/**
 * This function is used at runtime to get the actual type of the values produced by the given rule at runtime.
 * For data type rules, the name of the declared return type of the rule is returned (if any),
 * e.g. "INT_value returns number: MY_INT;" returns "number".
 * @param rule the given rule
 * @returns the name of the type of the produced values of the rule at runtime
 */
function getRuleType(rule) {
    var _a, _b, _c;
    if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalRule */ .MS(rule)) {
        return (_b = (_a = rule.type) === null || _a === void 0 ? void 0 : _a.name) !== null && _b !== void 0 ? _b : 'string';
    }
    else {
        return (_c = getExplicitRuleType(rule)) !== null && _c !== void 0 ? _c : rule.name;
    }
}
function terminalRegex(terminalRule) {
    const flags = {
        s: false,
        i: false,
        u: false
    };
    const source = abstractElementToRegex(terminalRule.definition, flags);
    const flagText = Object.entries(flags).filter(([, value]) => value).map(([name]) => name).join('');
    return new RegExp(source, flagText);
}
// Using [\s\S]* allows to match everything, compared to . which doesn't match line terminators
const WILDCARD = /[\s\S]/.source;
function abstractElementToRegex(element, flags) {
    if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalAlternatives */ .V7(element)) {
        return terminalAlternativesToRegex(element);
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalGroup */ .X9(element)) {
        return terminalGroupToRegex(element);
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isCharacterRange */ .Bf(element)) {
        return characterRangeToRegex(element);
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isTerminalRuleCall */ .gf(element)) {
        const rule = element.rule.ref;
        if (!rule) {
            throw new Error('Missing rule reference.');
        }
        return withCardinality(abstractElementToRegex(rule.definition), {
            cardinality: element.cardinality,
            lookahead: element.lookahead
        });
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isNegatedToken */ .Bi(element)) {
        return negateTokenToRegex(element);
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isUntilToken */ .OG(element)) {
        return untilTokenToRegex(element);
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isRegexToken */ .Sg(element)) {
        const lastSlash = element.regex.lastIndexOf('/');
        const source = element.regex.substring(1, lastSlash);
        const regexFlags = element.regex.substring(lastSlash + 1);
        if (flags) {
            flags.i = regexFlags.includes('i');
            flags.s = regexFlags.includes('s');
            flags.u = regexFlags.includes('u');
        }
        return withCardinality(source, {
            cardinality: element.cardinality,
            lookahead: element.lookahead,
            wrap: false
        });
    }
    else if (_languages_generated_ast_js__WEBPACK_IMPORTED_MODULE_0__/* .isWildcard */ .qm(element)) {
        return withCardinality(WILDCARD, {
            cardinality: element.cardinality,
            lookahead: element.lookahead
        });
    }
    else {
        throw new Error(`Invalid terminal element: ${element === null || element === void 0 ? void 0 : element.$type}`);
    }
}
function terminalAlternativesToRegex(alternatives) {
    return withCardinality(alternatives.elements.map(e => abstractElementToRegex(e)).join('|'), {
        cardinality: alternatives.cardinality,
        lookahead: alternatives.lookahead
    });
}
function terminalGroupToRegex(group) {
    return withCardinality(group.elements.map(e => abstractElementToRegex(e)).join(''), {
        cardinality: group.cardinality,
        lookahead: group.lookahead
    });
}
function untilTokenToRegex(until) {
    return withCardinality(`${WILDCARD}*?${abstractElementToRegex(until.terminal)}`, {
        cardinality: until.cardinality,
        lookahead: until.lookahead
    });
}
function negateTokenToRegex(negate) {
    return withCardinality(`(?!${abstractElementToRegex(negate.terminal)})${WILDCARD}*?`, {
        cardinality: negate.cardinality,
        lookahead: negate.lookahead
    });
}
function characterRangeToRegex(range) {
    if (range.right) {
        return withCardinality(`[${keywordToRegex(range.left)}-${keywordToRegex(range.right)}]`, {
            cardinality: range.cardinality,
            lookahead: range.lookahead,
            wrap: false
        });
    }
    return withCardinality(keywordToRegex(range.left), {
        cardinality: range.cardinality,
        lookahead: range.lookahead,
        wrap: false
    });
}
function keywordToRegex(keyword) {
    return (0,_regexp_utils_js__WEBPACK_IMPORTED_MODULE_2__/* .escapeRegExp */ .hr)(keyword.value);
}
function withCardinality(regex, options) {
    var _a;
    if (options.wrap !== false || options.lookahead) {
        regex = `(${(_a = options.lookahead) !== null && _a !== void 0 ? _a : ''}${regex})`;
    }
    if (options.cardinality) {
        return `${regex}${options.cardinality}`;
    }
    return regex;
}
//# sourceMappingURL=grammar-utils.js.map

/***/ }),

/***/ 43078:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   K0: () => (/* binding */ NEWLINE_REGEXP),
/* harmony export */   Rn: () => (/* binding */ isMultilineComment),
/* harmony export */   XC: () => (/* binding */ partialMatches),
/* harmony export */   cb: () => (/* binding */ isWhitespace),
/* harmony export */   cp: () => (/* binding */ getCaseInsensitivePattern),
/* harmony export */   hr: () => (/* binding */ escapeRegExp)
/* harmony export */ });
/* unused harmony exports getTerminalParts, whitespaceCharacters, partialRegExp */
/* harmony import */ var _chevrotain_regexp_to_ast__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(77647);
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/

const NEWLINE_REGEXP = /\r?\n/gm;
const regexpParser = new _chevrotain_regexp_to_ast__WEBPACK_IMPORTED_MODULE_0__/* .RegExpParser */ .O();
/**
 * This class is in charge of heuristically identifying start/end tokens of terminals.
 *
 * The way this works is by doing the following:
 * 1. Traverse the regular expression in the "start state"
 * 2. Add any encountered sets/single characters to the "start regexp"
 * 3. Once we encounter any variable-length content (i.e. with quantifiers such as +/?/*), we enter the "end state"
 * 4. In the end state, any sets/single characters are added to an "end stack".
 * 5. If we re-encounter any variable-length content we reset the end stack
 * 6. We continue visiting the regex until the end, reseting the end stack and rebuilding it as necessary
 *
 * After traversing a regular expression the `startRegexp/endRegexp` properties allow access to the stored start/end of the terminal
 */
class TerminalRegExpVisitor extends _chevrotain_regexp_to_ast__WEBPACK_IMPORTED_MODULE_0__/* .BaseRegExpVisitor */ .e {
    constructor() {
        super(...arguments);
        this.isStarting = true;
        this.endRegexpStack = [];
        this.multiline = false;
    }
    get endRegex() {
        return this.endRegexpStack.join('');
    }
    reset(regex) {
        this.multiline = false;
        this.regex = regex;
        this.startRegexp = '';
        this.isStarting = true;
        this.endRegexpStack = [];
    }
    visitGroup(node) {
        if (node.quantifier) {
            this.isStarting = false;
            this.endRegexpStack = [];
        }
    }
    visitCharacter(node) {
        const char = String.fromCharCode(node.value);
        if (!this.multiline && char === '\n') {
            this.multiline = true;
        }
        if (node.quantifier) {
            this.isStarting = false;
            this.endRegexpStack = [];
        }
        else {
            const escapedChar = escapeRegExp(char);
            this.endRegexpStack.push(escapedChar);
            if (this.isStarting) {
                this.startRegexp += escapedChar;
            }
        }
    }
    visitSet(node) {
        if (!this.multiline) {
            const set = this.regex.substring(node.loc.begin, node.loc.end);
            const regex = new RegExp(set);
            this.multiline = Boolean('\n'.match(regex));
        }
        if (node.quantifier) {
            this.isStarting = false;
            this.endRegexpStack = [];
        }
        else {
            const set = this.regex.substring(node.loc.begin, node.loc.end);
            this.endRegexpStack.push(set);
            if (this.isStarting) {
                this.startRegexp += set;
            }
        }
    }
    visitChildren(node) {
        if (node.type === 'Group') {
            // Ignore children of groups with quantifier (+/*/?)
            // These groups are unrelated to start/end tokens of terminals
            const group = node;
            if (group.quantifier) {
                return;
            }
        }
        super.visitChildren(node);
    }
}
const visitor = new TerminalRegExpVisitor();
function getTerminalParts(regexp) {
    try {
        if (typeof regexp !== 'string') {
            regexp = regexp.source;
        }
        regexp = `/${regexp}/`;
        const pattern = regexpParser.pattern(regexp);
        const parts = [];
        for (const alternative of pattern.value.value) {
            visitor.reset(regexp);
            visitor.visit(alternative);
            parts.push({
                start: visitor.startRegexp,
                end: visitor.endRegex
            });
        }
        return parts;
    }
    catch (_a) {
        return [];
    }
}
function isMultilineComment(regexp) {
    try {
        if (typeof regexp === 'string') {
            regexp = new RegExp(regexp);
        }
        regexp = regexp.toString();
        visitor.reset(regexp);
        // Parsing the pattern might fail (since it's user code)
        visitor.visit(regexpParser.pattern(regexp));
        return visitor.multiline;
    }
    catch (_a) {
        return false;
    }
}
/**
 * A set of all characters that are considered whitespace by the '\s' RegExp character class.
 * Taken from [MDN](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Regular_expressions/Character_classes).
 */
const whitespaceCharacters = ('\f\n\r\t\v\u0020\u00a0\u1680\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007' +
    '\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000\ufeff').split('');
function isWhitespace(value) {
    const regexp = typeof value === 'string' ? new RegExp(value) : value;
    return whitespaceCharacters.some((ws) => regexp.test(ws));
}
function escapeRegExp(value) {
    return value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
function getCaseInsensitivePattern(keyword) {
    return Array.prototype.map.call(keyword, letter => /\w/.test(letter) ? `[${letter.toLowerCase()}${letter.toUpperCase()}]` : escapeRegExp(letter)).join('');
}
/**
 * Determines whether the given input has a partial match with the specified regex.
 * @param regex The regex to partially match against
 * @param input The input string
 * @returns Whether any match exists.
 */
function partialMatches(regex, input) {
    const partial = partialRegExp(regex);
    const match = input.match(partial);
    return !!match && match[0].length > 0;
}
/**
 * Builds a partial regex from the input regex. A partial regex is able to match incomplete input strings. E.g.
 * a partial regex constructed from `/ab/` is able to match the string `a` without needing a following `b` character. However it won't match `b` alone.
 * @param regex The input regex to be converted.
 * @returns A partial regex constructed from the input regex.
 */
function partialRegExp(regex) {
    if (typeof regex === 'string') {
        regex = new RegExp(regex);
    }
    const re = regex, source = regex.source;
    let i = 0;
    function process() {
        let result = '', tmp;
        function appendRaw(nbChars) {
            result += source.substr(i, nbChars);
            i += nbChars;
        }
        function appendOptional(nbChars) {
            result += '(?:' + source.substr(i, nbChars) + '|$)';
            i += nbChars;
        }
        while (i < source.length) {
            switch (source[i]) {
                case '\\':
                    switch (source[i + 1]) {
                        case 'c':
                            appendOptional(3);
                            break;
                        case 'x':
                            appendOptional(4);
                            break;
                        case 'u':
                            if (re.unicode) {
                                if (source[i + 2] === '{') {
                                    appendOptional(source.indexOf('}', i) - i + 1);
                                }
                                else {
                                    appendOptional(6);
                                }
                            }
                            else {
                                appendOptional(2);
                            }
                            break;
                        case 'p':
                        case 'P':
                            if (re.unicode) {
                                appendOptional(source.indexOf('}', i) - i + 1);
                            }
                            else {
                                appendOptional(2);
                            }
                            break;
                        case 'k':
                            appendOptional(source.indexOf('>', i) - i + 1);
                            break;
                        default:
                            appendOptional(2);
                            break;
                    }
                    break;
                case '[':
                    tmp = /\[(?:\\.|.)*?\]/g;
                    tmp.lastIndex = i;
                    tmp = tmp.exec(source) || [];
                    appendOptional(tmp[0].length);
                    break;
                case '|':
                case '^':
                case '$':
                case '*':
                case '+':
                case '?':
                    appendRaw(1);
                    break;
                case '{':
                    tmp = /\{\d+,?\d*\}/g;
                    tmp.lastIndex = i;
                    tmp = tmp.exec(source);
                    if (tmp) {
                        appendRaw(tmp[0].length);
                    }
                    else {
                        appendOptional(1);
                    }
                    break;
                case '(':
                    if (source[i + 1] === '?') {
                        switch (source[i + 2]) {
                            case ':':
                                result += '(?:';
                                i += 3;
                                result += process() + '|$)';
                                break;
                            case '=':
                                result += '(?=';
                                i += 3;
                                result += process() + ')';
                                break;
                            case '!':
                                tmp = i;
                                i += 3;
                                process();
                                result += source.substr(tmp, i - tmp);
                                break;
                            case '<':
                                switch (source[i + 3]) {
                                    case '=':
                                    case '!':
                                        tmp = i;
                                        i += 4;
                                        process();
                                        result += source.substr(tmp, i - tmp);
                                        break;
                                    default:
                                        appendRaw(source.indexOf('>', i) - i + 1);
                                        result += process() + '|$)';
                                        break;
                                }
                                break;
                        }
                    }
                    else {
                        appendRaw(1);
                        result += process() + '|$)';
                    }
                    break;
                case ')':
                    ++i;
                    return result;
                default:
                    appendOptional(1);
                    break;
            }
        }
        return result;
    }
    return new RegExp(process(), regex.flags);
}
//# sourceMappingURL=regexp-utils.js.map

/***/ }),

/***/ 99293:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Cl: () => (/* binding */ EMPTY_STREAM),
/* harmony export */   IH: () => (/* binding */ Reduction),
/* harmony export */   Ry: () => (/* binding */ DONE_RESULT),
/* harmony export */   Vw: () => (/* binding */ stream),
/* harmony export */   i: () => (/* binding */ StreamImpl),
/* harmony export */   i8: () => (/* binding */ TreeStreamImpl)
/* harmony export */ });
/******************************************************************************
 * Copyright 2021 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
/**
 * The default implementation of `Stream` works with two input functions:
 *  - The first function creates the initial state of an iteration.
 *  - The second function gets the current state as argument and returns an `IteratorResult`.
 */
class StreamImpl {
    constructor(startFn, nextFn) {
        this.startFn = startFn;
        this.nextFn = nextFn;
    }
    iterator() {
        const iterator = {
            state: this.startFn(),
            next: () => this.nextFn(iterator.state),
            [Symbol.iterator]: () => iterator
        };
        return iterator;
    }
    [Symbol.iterator]() {
        return this.iterator();
    }
    isEmpty() {
        const iterator = this.iterator();
        return Boolean(iterator.next().done);
    }
    count() {
        const iterator = this.iterator();
        let count = 0;
        let next = iterator.next();
        while (!next.done) {
            count++;
            next = iterator.next();
        }
        return count;
    }
    toArray() {
        const result = [];
        const iterator = this.iterator();
        let next;
        do {
            next = iterator.next();
            if (next.value !== undefined) {
                result.push(next.value);
            }
        } while (!next.done);
        return result;
    }
    toSet() {
        return new Set(this);
    }
    toMap(keyFn, valueFn) {
        const entryStream = this.map(element => [
            keyFn ? keyFn(element) : element,
            valueFn ? valueFn(element) : element
        ]);
        return new Map(entryStream);
    }
    toString() {
        return this.join();
    }
    concat(other) {
        return new StreamImpl(() => ({ first: this.startFn(), firstDone: false, iterator: other[Symbol.iterator]() }), state => {
            let result;
            if (!state.firstDone) {
                do {
                    result = this.nextFn(state.first);
                    if (!result.done) {
                        return result;
                    }
                } while (!result.done);
                state.firstDone = true;
            }
            do {
                result = state.iterator.next();
                if (!result.done) {
                    return result;
                }
            } while (!result.done);
            return DONE_RESULT;
        });
    }
    join(separator = ',') {
        const iterator = this.iterator();
        let value = '';
        let result;
        let addSeparator = false;
        do {
            result = iterator.next();
            if (!result.done) {
                if (addSeparator) {
                    value += separator;
                }
                value += toString(result.value);
            }
            addSeparator = true;
        } while (!result.done);
        return value;
    }
    indexOf(searchElement, fromIndex = 0) {
        const iterator = this.iterator();
        let index = 0;
        let next = iterator.next();
        while (!next.done) {
            if (index >= fromIndex && next.value === searchElement) {
                return index;
            }
            next = iterator.next();
            index++;
        }
        return -1;
    }
    every(predicate) {
        const iterator = this.iterator();
        let next = iterator.next();
        while (!next.done) {
            if (!predicate(next.value)) {
                return false;
            }
            next = iterator.next();
        }
        return true;
    }
    some(predicate) {
        const iterator = this.iterator();
        let next = iterator.next();
        while (!next.done) {
            if (predicate(next.value)) {
                return true;
            }
            next = iterator.next();
        }
        return false;
    }
    forEach(callbackfn) {
        const iterator = this.iterator();
        let index = 0;
        let next = iterator.next();
        while (!next.done) {
            callbackfn(next.value, index);
            next = iterator.next();
            index++;
        }
    }
    map(callbackfn) {
        return new StreamImpl(this.startFn, (state) => {
            const { done, value } = this.nextFn(state);
            if (done) {
                return DONE_RESULT;
            }
            else {
                return { done: false, value: callbackfn(value) };
            }
        });
    }
    filter(predicate) {
        return new StreamImpl(this.startFn, state => {
            let result;
            do {
                result = this.nextFn(state);
                if (!result.done && predicate(result.value)) {
                    return result;
                }
            } while (!result.done);
            return DONE_RESULT;
        });
    }
    nonNullable() {
        return this.filter(e => e !== undefined && e !== null);
    }
    reduce(callbackfn, initialValue) {
        const iterator = this.iterator();
        let previousValue = initialValue;
        let next = iterator.next();
        while (!next.done) {
            if (previousValue === undefined) {
                previousValue = next.value;
            }
            else {
                previousValue = callbackfn(previousValue, next.value);
            }
            next = iterator.next();
        }
        return previousValue;
    }
    reduceRight(callbackfn, initialValue) {
        return this.recursiveReduce(this.iterator(), callbackfn, initialValue);
    }
    recursiveReduce(iterator, callbackfn, initialValue) {
        const next = iterator.next();
        if (next.done) {
            return initialValue;
        }
        const previousValue = this.recursiveReduce(iterator, callbackfn, initialValue);
        if (previousValue === undefined) {
            return next.value;
        }
        return callbackfn(previousValue, next.value);
    }
    find(predicate) {
        const iterator = this.iterator();
        let next = iterator.next();
        while (!next.done) {
            if (predicate(next.value)) {
                return next.value;
            }
            next = iterator.next();
        }
        return undefined;
    }
    findIndex(predicate) {
        const iterator = this.iterator();
        let index = 0;
        let next = iterator.next();
        while (!next.done) {
            if (predicate(next.value)) {
                return index;
            }
            next = iterator.next();
            index++;
        }
        return -1;
    }
    includes(searchElement) {
        const iterator = this.iterator();
        let next = iterator.next();
        while (!next.done) {
            if (next.value === searchElement) {
                return true;
            }
            next = iterator.next();
        }
        return false;
    }
    flatMap(callbackfn) {
        return new StreamImpl(() => ({ this: this.startFn() }), (state) => {
            do {
                if (state.iterator) {
                    const next = state.iterator.next();
                    if (next.done) {
                        state.iterator = undefined;
                    }
                    else {
                        return next;
                    }
                }
                const { done, value } = this.nextFn(state.this);
                if (!done) {
                    const mapped = callbackfn(value);
                    if (isIterable(mapped)) {
                        state.iterator = mapped[Symbol.iterator]();
                    }
                    else {
                        return { done: false, value: mapped };
                    }
                }
            } while (state.iterator);
            return DONE_RESULT;
        });
    }
    flat(depth) {
        if (depth === undefined) {
            depth = 1;
        }
        if (depth <= 0) {
            return this;
        }
        const stream = depth > 1 ? this.flat(depth - 1) : this;
        return new StreamImpl(() => ({ this: stream.startFn() }), (state) => {
            do {
                if (state.iterator) {
                    const next = state.iterator.next();
                    if (next.done) {
                        state.iterator = undefined;
                    }
                    else {
                        return next;
                    }
                }
                const { done, value } = stream.nextFn(state.this);
                if (!done) {
                    if (isIterable(value)) {
                        state.iterator = value[Symbol.iterator]();
                    }
                    else {
                        return { done: false, value: value };
                    }
                }
            } while (state.iterator);
            return DONE_RESULT;
        });
    }
    head() {
        const iterator = this.iterator();
        const result = iterator.next();
        if (result.done) {
            return undefined;
        }
        return result.value;
    }
    tail(skipCount = 1) {
        return new StreamImpl(() => {
            const state = this.startFn();
            for (let i = 0; i < skipCount; i++) {
                const next = this.nextFn(state);
                if (next.done) {
                    return state;
                }
            }
            return state;
        }, this.nextFn);
    }
    limit(maxSize) {
        return new StreamImpl(() => ({ size: 0, state: this.startFn() }), state => {
            state.size++;
            if (state.size > maxSize) {
                return DONE_RESULT;
            }
            return this.nextFn(state.state);
        });
    }
    distinct(by) {
        return new StreamImpl(() => ({ set: new Set(), internalState: this.startFn() }), state => {
            let result;
            do {
                result = this.nextFn(state.internalState);
                if (!result.done) {
                    const value = by ? by(result.value) : result.value;
                    if (!state.set.has(value)) {
                        state.set.add(value);
                        return result;
                    }
                }
            } while (!result.done);
            return DONE_RESULT;
        });
    }
    exclude(other, key) {
        const otherKeySet = new Set();
        for (const item of other) {
            const value = key ? key(item) : item;
            otherKeySet.add(value);
        }
        return this.filter(e => {
            const ownKey = key ? key(e) : e;
            return !otherKeySet.has(ownKey);
        });
    }
}
function toString(item) {
    if (typeof item === 'string') {
        return item;
    }
    if (typeof item === 'undefined') {
        return 'undefined';
    }
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    if (typeof item.toString === 'function') {
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        return item.toString();
    }
    return Object.prototype.toString.call(item);
}
function isIterable(obj) {
    return !!obj && typeof obj[Symbol.iterator] === 'function';
}
/**
 * An empty stream of any type.
 */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
const EMPTY_STREAM = new StreamImpl(() => undefined, () => DONE_RESULT);
/**
 * Use this `IteratorResult` when implementing a `StreamImpl` to indicate that there are no more elements in the stream.
 */
const DONE_RESULT = Object.freeze({ done: true, value: undefined });
/**
 * Create a stream from one or more iterables or array-likes.
 */
function stream(...collections) {
    if (collections.length === 1) {
        const collection = collections[0];
        if (collection instanceof StreamImpl) {
            return collection;
        }
        if (isIterable(collection)) {
            return new StreamImpl(() => collection[Symbol.iterator](), (iterator) => iterator.next());
        }
        if (typeof collection.length === 'number') {
            return new StreamImpl(() => ({ index: 0 }), (state) => {
                if (state.index < collection.length) {
                    return { done: false, value: collection[state.index++] };
                }
                else {
                    return DONE_RESULT;
                }
            });
        }
    }
    if (collections.length > 1) {
        return new StreamImpl(() => ({ collIndex: 0, arrIndex: 0 }), (state) => {
            do {
                if (state.iterator) {
                    const next = state.iterator.next();
                    if (!next.done) {
                        return next;
                    }
                    state.iterator = undefined;
                }
                if (state.array) {
                    if (state.arrIndex < state.array.length) {
                        return { done: false, value: state.array[state.arrIndex++] };
                    }
                    state.array = undefined;
                    state.arrIndex = 0;
                }
                if (state.collIndex < collections.length) {
                    const collection = collections[state.collIndex++];
                    if (isIterable(collection)) {
                        state.iterator = collection[Symbol.iterator]();
                    }
                    else if (collection && typeof collection.length === 'number') {
                        state.array = collection;
                    }
                }
            } while (state.iterator || state.array || state.collIndex < collections.length);
            return DONE_RESULT;
        });
    }
    return EMPTY_STREAM;
}
/**
 * The default implementation of `TreeStream` takes a root element and a function that computes the
 * children of its argument. Whether the root node included in the stream is controlled with the
 * `includeRoot` option, which defaults to `false`.
 */
class TreeStreamImpl extends StreamImpl {
    constructor(root, children, options) {
        super(() => ({
            iterators: (options === null || options === void 0 ? void 0 : options.includeRoot) ? [[root][Symbol.iterator]()] : [children(root)[Symbol.iterator]()],
            pruned: false
        }), state => {
            if (state.pruned) {
                state.iterators.pop();
                state.pruned = false;
            }
            while (state.iterators.length > 0) {
                const iterator = state.iterators[state.iterators.length - 1];
                const next = iterator.next();
                if (next.done) {
                    state.iterators.pop();
                }
                else {
                    state.iterators.push(children(next.value)[Symbol.iterator]());
                    return next;
                }
            }
            return DONE_RESULT;
        });
    }
    iterator() {
        const iterator = {
            state: this.startFn(),
            next: () => this.nextFn(iterator.state),
            prune: () => {
                iterator.state.pruned = true;
            },
            [Symbol.iterator]: () => iterator
        };
        return iterator;
    }
}
/**
 * A set of utility functions that reduce a stream to a single value.
 */
var Reduction;
(function (Reduction) {
    /**
     * Compute the sum of a number stream.
     */
    function sum(stream) {
        return stream.reduce((a, b) => a + b, 0);
    }
    Reduction.sum = sum;
    /**
     * Compute the product of a number stream.
     */
    function product(stream) {
        return stream.reduce((a, b) => a * b, 0);
    }
    Reduction.product = product;
    /**
     * Compute the minimum of a number stream. Returns `undefined` if the stream is empty.
     */
    function min(stream) {
        return stream.reduce((a, b) => Math.min(a, b));
    }
    Reduction.min = min;
    /**
     * Compute the maximum of a number stream. Returns `undefined` if the stream is empty.
     */
    function max(stream) {
        return stream.reduce((a, b) => Math.max(a, b));
    }
    Reduction.max = max;
})(Reduction || (Reduction = {}));
//# sourceMappingURL=stream.js.map

/***/ }),

/***/ 44014:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   u: () => (/* binding */ EmptyFileSystem)
/* harmony export */ });
/* unused harmony export EmptyFileSystemProvider */
/******************************************************************************
 * Copyright 2022 TypeFox GmbH
 * This program and the accompanying materials are made available under the
 * terms of the MIT License, which is available in the project root.
 ******************************************************************************/
class EmptyFileSystemProvider {
    readFile() {
        throw new Error('No file system is available.');
    }
    async readDirectory() {
        return [];
    }
}
const EmptyFileSystem = {
    fileSystemProvider: () => new EmptyFileSystemProvider()
};
//# sourceMappingURL=file-system-provider.js.map

/***/ }),

/***/ 41589:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _isSymbol_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(59660);


/**
 * The base implementation of methods like `_.max` and `_.min` which accepts a
 * `comparator` to determine the extremum value.
 *
 * @private
 * @param {Array} array The array to iterate over.
 * @param {Function} iteratee The iteratee invoked per iteration.
 * @param {Function} comparator The comparator used to compare values.
 * @returns {*} Returns the extremum value.
 */
function baseExtremum(array, iteratee, comparator) {
  var index = -1,
      length = array.length;

  while (++index < length) {
    var value = array[index],
        current = iteratee(value);

    if (current != null && (computed === undefined
          ? (current === current && !(0,_isSymbol_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(current))
          : comparator(current, computed)
        )) {
      var computed = current,
          result = value;
    }
  }
  return result;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseExtremum);


/***/ }),

/***/ 79520:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * The base implementation of `_.lt` which doesn't coerce arguments.
 *
 * @private
 * @param {*} value The value to compare.
 * @param {*} other The other value to compare.
 * @returns {boolean} Returns `true` if `value` is less than `other`,
 *  else `false`.
 */
function baseLt(value, other) {
  return value < other;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseLt);


/***/ }),

/***/ 15521:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseEach_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(77201);
/* harmony import */ var _isArrayLike_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(69959);



/**
 * The base implementation of `_.map` without support for iteratee shorthands.
 *
 * @private
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} iteratee The function invoked per iteration.
 * @returns {Array} Returns the new mapped array.
 */
function baseMap(collection, iteratee) {
  var index = -1,
      result = (0,_isArrayLike_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(collection) ? Array(collection.length) : [];

  (0,_baseEach_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(collection, function(value, key, collection) {
    result[++index] = iteratee(value, key, collection);
  });
  return result;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseMap);


/***/ }),

/***/ 73338:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _basePickBy)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGet.js
var _baseGet = __webpack_require__(78402);
// EXTERNAL MODULE: ../node_modules/lodash-es/_assignValue.js
var _assignValue = __webpack_require__(15561);
// EXTERNAL MODULE: ../node_modules/lodash-es/_castPath.js + 2 modules
var _castPath = __webpack_require__(94022);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isIndex.js
var _isIndex = __webpack_require__(8616);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
// EXTERNAL MODULE: ../node_modules/lodash-es/_toKey.js
var _toKey = __webpack_require__(13550);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseSet.js






/**
 * The base implementation of `_.set`.
 *
 * @private
 * @param {Object} object The object to modify.
 * @param {Array|string} path The path of the property to set.
 * @param {*} value The value to set.
 * @param {Function} [customizer] The function to customize path creation.
 * @returns {Object} Returns `object`.
 */
function baseSet(object, path, value, customizer) {
  if (!(0,isObject/* default */.Z)(object)) {
    return object;
  }
  path = (0,_castPath/* default */.Z)(path, object);

  var index = -1,
      length = path.length,
      lastIndex = length - 1,
      nested = object;

  while (nested != null && ++index < length) {
    var key = (0,_toKey/* default */.Z)(path[index]),
        newValue = value;

    if (key === '__proto__' || key === 'constructor' || key === 'prototype') {
      return object;
    }

    if (index != lastIndex) {
      var objValue = nested[key];
      newValue = customizer ? customizer(objValue, key, nested) : undefined;
      if (newValue === undefined) {
        newValue = (0,isObject/* default */.Z)(objValue)
          ? objValue
          : ((0,_isIndex/* default */.Z)(path[index + 1]) ? [] : {});
      }
    }
    (0,_assignValue/* default */.Z)(nested, key, newValue);
    nested = nested[key];
  }
  return object;
}

/* harmony default export */ const _baseSet = (baseSet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_basePickBy.js




/**
 * The base implementation of  `_.pickBy` without support for iteratee shorthands.
 *
 * @private
 * @param {Object} object The source object.
 * @param {string[]} paths The property paths to pick.
 * @param {Function} predicate The function invoked per property.
 * @returns {Object} Returns the new object.
 */
function basePickBy(object, paths, predicate) {
  var index = -1,
      length = paths.length,
      result = {};

  while (++index < length) {
    var path = paths[index],
        value = (0,_baseGet/* default */.Z)(object, path);

    if (predicate(value, path)) {
      _baseSet(result, (0,_castPath/* default */.Z)(path, object), value);
    }
  }
  return result;
}

/* harmony default export */ const _basePickBy = (basePickBy);


/***/ }),

/***/ 20190:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseClone_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(52390);


/** Used to compose bitmasks for cloning. */
var CLONE_SYMBOLS_FLAG = 4;

/**
 * Creates a shallow clone of `value`.
 *
 * **Note:** This method is loosely based on the
 * [structured clone algorithm](https://mdn.io/Structured_clone_algorithm)
 * and supports cloning arrays, array buffers, booleans, date objects, maps,
 * numbers, `Object` objects, regexes, sets, strings, symbols, and typed
 * arrays. The own enumerable properties of `arguments` objects are cloned
 * as plain objects. An empty object is returned for uncloneable values such
 * as error objects, functions, DOM nodes, and WeakMaps.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to clone.
 * @returns {*} Returns the cloned value.
 * @see _.cloneDeep
 * @example
 *
 * var objects = [{ 'a': 1 }, { 'b': 2 }];
 *
 * var shallow = _.clone(objects);
 * console.log(shallow[0] === objects[0]);
 * // => true
 */
function clone(value) {
  return (0,_baseClone_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value, CLONE_SYMBOLS_FLAG);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (clone);


/***/ }),

/***/ 65479:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseRest_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(99719);
/* harmony import */ var _eq_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(35050);
/* harmony import */ var _isIterateeCall_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(47952);
/* harmony import */ var _keysIn_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(48441);





/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Assigns own and inherited enumerable string keyed properties of source
 * objects to the destination object for all destination properties that
 * resolve to `undefined`. Source objects are applied from left to right.
 * Once a property is set, additional values of the same property are ignored.
 *
 * **Note:** This method mutates `object`.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Object
 * @param {Object} object The destination object.
 * @param {...Object} [sources] The source objects.
 * @returns {Object} Returns `object`.
 * @see _.defaultsDeep
 * @example
 *
 * _.defaults({ 'a': 1 }, { 'b': 2 }, { 'a': 3 });
 * // => { 'a': 1, 'b': 2 }
 */
var defaults = (0,_baseRest_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(function(object, sources) {
  object = Object(object);

  var index = -1;
  var length = sources.length;
  var guard = length > 2 ? sources[2] : undefined;

  if (guard && (0,_isIterateeCall_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(sources[0], sources[1], guard)) {
    length = 1;
  }

  while (++index < length) {
    var source = sources[index];
    var props = (0,_keysIn_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z)(source);
    var propsIndex = -1;
    var propsLength = props.length;

    while (++propsIndex < propsLength) {
      var key = props[propsIndex];
      var value = object[key];

      if (value === undefined ||
          ((0,_eq_js__WEBPACK_IMPORTED_MODULE_3__/* ["default"] */ .Z)(value, objectProto[key]) && !hasOwnProperty.call(object, key))) {
        object[key] = source[key];
      }
    }
  }

  return object;
});

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (defaults);


/***/ }),

/***/ 90970:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_find)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseIteratee.js + 15 modules
var _baseIteratee = __webpack_require__(86494);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArrayLike.js
var isArrayLike = __webpack_require__(69959);
// EXTERNAL MODULE: ../node_modules/lodash-es/keys.js
var keys = __webpack_require__(11723);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_createFind.js




/**
 * Creates a `_.find` or `_.findLast` function.
 *
 * @private
 * @param {Function} findIndexFunc The function to find the collection index.
 * @returns {Function} Returns the new find function.
 */
function createFind(findIndexFunc) {
  return function(collection, predicate, fromIndex) {
    var iterable = Object(collection);
    if (!(0,isArrayLike/* default */.Z)(collection)) {
      var iteratee = (0,_baseIteratee/* default */.Z)(predicate, 3);
      collection = (0,keys/* default */.Z)(collection);
      predicate = function(key) { return iteratee(iterable[key], key, iterable); };
    }
    var index = findIndexFunc(collection, predicate, fromIndex);
    return index > -1 ? iterable[iteratee ? collection[index] : index] : undefined;
  };
}

/* harmony default export */ const _createFind = (createFind);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseFindIndex.js
var _baseFindIndex = __webpack_require__(9872);
// EXTERNAL MODULE: ../node_modules/lodash-es/toInteger.js
var toInteger = __webpack_require__(98670);
;// CONCATENATED MODULE: ../node_modules/lodash-es/findIndex.js




/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeMax = Math.max;

/**
 * This method is like `_.find` except that it returns the index of the first
 * element `predicate` returns truthy for instead of the element itself.
 *
 * @static
 * @memberOf _
 * @since 1.1.0
 * @category Array
 * @param {Array} array The array to inspect.
 * @param {Function} [predicate=_.identity] The function invoked per iteration.
 * @param {number} [fromIndex=0] The index to search from.
 * @returns {number} Returns the index of the found element, else `-1`.
 * @example
 *
 * var users = [
 *   { 'user': 'barney',  'active': false },
 *   { 'user': 'fred',    'active': false },
 *   { 'user': 'pebbles', 'active': true }
 * ];
 *
 * _.findIndex(users, function(o) { return o.user == 'barney'; });
 * // => 0
 *
 * // The `_.matches` iteratee shorthand.
 * _.findIndex(users, { 'user': 'fred', 'active': false });
 * // => 1
 *
 * // The `_.matchesProperty` iteratee shorthand.
 * _.findIndex(users, ['active', false]);
 * // => 0
 *
 * // The `_.property` iteratee shorthand.
 * _.findIndex(users, 'active');
 * // => 2
 */
function findIndex(array, predicate, fromIndex) {
  var length = array == null ? 0 : array.length;
  if (!length) {
    return -1;
  }
  var index = fromIndex == null ? 0 : (0,toInteger/* default */.Z)(fromIndex);
  if (index < 0) {
    index = nativeMax(length + index, 0);
  }
  return (0,_baseFindIndex/* default */.Z)(array, (0,_baseIteratee/* default */.Z)(predicate, 3), index);
}

/* harmony default export */ const lodash_es_findIndex = (findIndex);

;// CONCATENATED MODULE: ../node_modules/lodash-es/find.js



/**
 * Iterates over elements of `collection`, returning the first element
 * `predicate` returns truthy for. The predicate is invoked with three
 * arguments: (value, index|key, collection).
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to inspect.
 * @param {Function} [predicate=_.identity] The function invoked per iteration.
 * @param {number} [fromIndex=0] The index to search from.
 * @returns {*} Returns the matched element, else `undefined`.
 * @example
 *
 * var users = [
 *   { 'user': 'barney',  'age': 36, 'active': true },
 *   { 'user': 'fred',    'age': 40, 'active': false },
 *   { 'user': 'pebbles', 'age': 1,  'active': true }
 * ];
 *
 * _.find(users, function(o) { return o.age < 40; });
 * // => object for 'barney'
 *
 * // The `_.matches` iteratee shorthand.
 * _.find(users, { 'age': 1, 'active': true });
 * // => object for 'pebbles'
 *
 * // The `_.matchesProperty` iteratee shorthand.
 * _.find(users, ['active', false]);
 * // => object for 'fred'
 *
 * // The `_.property` iteratee shorthand.
 * _.find(users, 'active');
 * // => object for 'barney'
 */
var find = _createFind(lodash_es_findIndex);

/* harmony default export */ const lodash_es_find = (find);


/***/ }),

/***/ 34134:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseFlatten_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(65029);
/* harmony import */ var _map_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(12930);



/**
 * Creates a flattened array of values by running each element in `collection`
 * thru `iteratee` and flattening the mapped results. The iteratee is invoked
 * with three arguments: (value, index|key, collection).
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [iteratee=_.identity] The function invoked per iteration.
 * @returns {Array} Returns the new flattened array.
 * @example
 *
 * function duplicate(n) {
 *   return [n, n];
 * }
 *
 * _.flatMap([1, 2], duplicate);
 * // => [1, 1, 2, 2]
 */
function flatMap(collection, iteratee) {
  return (0,_baseFlatten_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)((0,_map_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(collection, iteratee), 1);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (flatMap);


/***/ }),

/***/ 28099:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseFlatten_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(65029);


/**
 * Flattens `array` a single level deep.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to flatten.
 * @returns {Array} Returns the new flattened array.
 * @example
 *
 * _.flatten([1, [2, [3, [4]], 5]]);
 * // => [1, 2, [3, [4]], 5]
 */
function flatten(array) {
  var length = array == null ? 0 : array.length;
  return length ? (0,_baseFlatten_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(array, 1) : [];
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (flatten);


/***/ }),

/***/ 49493:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_has)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseHas.js
/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _baseHas_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * The base implementation of `_.has` without support for deep paths.
 *
 * @private
 * @param {Object} [object] The object to query.
 * @param {Array|string} key The key to check.
 * @returns {boolean} Returns `true` if `key` exists, else `false`.
 */
function baseHas(object, key) {
  return object != null && _baseHas_hasOwnProperty.call(object, key);
}

/* harmony default export */ const _baseHas = (baseHas);

// EXTERNAL MODULE: ../node_modules/lodash-es/_hasPath.js
var _hasPath = __webpack_require__(18625);
;// CONCATENATED MODULE: ../node_modules/lodash-es/has.js



/**
 * Checks if `path` is a direct property of `object`.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Object
 * @param {Object} object The object to query.
 * @param {Array|string} path The path to check.
 * @returns {boolean} Returns `true` if `path` exists, else `false`.
 * @example
 *
 * var object = { 'a': { 'b': 2 } };
 * var other = _.create({ 'a': _.create({ 'b': 2 }) });
 *
 * _.has(object, 'a');
 * // => true
 *
 * _.has(object, 'a.b');
 * // => true
 *
 * _.has(object, ['a', 'b']);
 * // => true
 *
 * _.has(other, 'a');
 * // => false
 */
function has(object, path) {
  return object != null && (0,_hasPath/* default */.Z)(object, path, _baseHas);
}

/* harmony default export */ const lodash_es_has = (has);


/***/ }),

/***/ 75732:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseGetTag_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(77070);
/* harmony import */ var _isArray_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(64058);
/* harmony import */ var _isObjectLike_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(9615);




/** `Object#toString` result references. */
var stringTag = '[object String]';

/**
 * Checks if `value` is classified as a `String` primitive or object.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a string, else `false`.
 * @example
 *
 * _.isString('abc');
 * // => true
 *
 * _.isString(1);
 * // => false
 */
function isString(value) {
  return typeof value == 'string' ||
    (!(0,_isArray_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value) && (0,_isObjectLike_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(value) && (0,_baseGetTag_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z)(value) == stringTag);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isString);


/***/ }),

/***/ 36411:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Gets the last element of `array`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Array
 * @param {Array} array The array to query.
 * @returns {*} Returns the last element of `array`.
 * @example
 *
 * _.last([1, 2, 3]);
 * // => 3
 */
function last(array) {
  var length = array == null ? 0 : array.length;
  return length ? array[length - 1] : undefined;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (last);


/***/ }),

/***/ 12930:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _arrayMap_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(33043);
/* harmony import */ var _baseIteratee_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86494);
/* harmony import */ var _baseMap_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(15521);
/* harmony import */ var _isArray_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(64058);





/**
 * Creates an array of values by running each element in `collection` thru
 * `iteratee`. The iteratee is invoked with three arguments:
 * (value, index|key, collection).
 *
 * Many lodash methods are guarded to work as iteratees for methods like
 * `_.every`, `_.filter`, `_.map`, `_.mapValues`, `_.reject`, and `_.some`.
 *
 * The guarded methods are:
 * `ary`, `chunk`, `curry`, `curryRight`, `drop`, `dropRight`, `every`,
 * `fill`, `invert`, `parseInt`, `random`, `range`, `rangeRight`, `repeat`,
 * `sampleSize`, `slice`, `some`, `sortBy`, `split`, `take`, `takeRight`,
 * `template`, `trim`, `trimEnd`, `trimStart`, and `words`
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Collection
 * @param {Array|Object} collection The collection to iterate over.
 * @param {Function} [iteratee=_.identity] The function invoked per iteration.
 * @returns {Array} Returns the new mapped array.
 * @example
 *
 * function square(n) {
 *   return n * n;
 * }
 *
 * _.map([4, 8], square);
 * // => [16, 64]
 *
 * _.map({ 'a': 4, 'b': 8 }, square);
 * // => [16, 64] (iteration order is not guaranteed)
 *
 * var users = [
 *   { 'user': 'barney' },
 *   { 'user': 'fred' }
 * ];
 *
 * // The `_.property` iteratee shorthand.
 * _.map(users, 'user');
 * // => ['barney', 'fred']
 */
function map(collection, iteratee) {
  var func = (0,_isArray_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(collection) ? _arrayMap_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z : _baseMap_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z;
  return func(collection, (0,_baseIteratee_js__WEBPACK_IMPORTED_MODULE_3__/* ["default"] */ .Z)(iteratee, 3));
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (map);


/***/ }),

/***/ 18519:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseExtremum_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(41589);
/* harmony import */ var _baseLt_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(79520);
/* harmony import */ var _identity_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(64056);




/**
 * Computes the minimum value of `array`. If `array` is empty or falsey,
 * `undefined` is returned.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Math
 * @param {Array} array The array to iterate over.
 * @returns {*} Returns the minimum value.
 * @example
 *
 * _.min([4, 2, 8, 6]);
 * // => 2
 *
 * _.min([]);
 * // => undefined
 */
function min(array) {
  return (array && array.length)
    ? (0,_baseExtremum_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(array, _identity_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z, _baseLt_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z)
    : undefined;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (min);


/***/ }),

/***/ 41291:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_toFinite)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_trimmedEndIndex.js
/** Used to match a single whitespace character. */
var reWhitespace = /\s/;

/**
 * Used by `_.trim` and `_.trimEnd` to get the index of the last non-whitespace
 * character of `string`.
 *
 * @private
 * @param {string} string The string to inspect.
 * @returns {number} Returns the index of the last non-whitespace character.
 */
function trimmedEndIndex(string) {
  var index = string.length;

  while (index-- && reWhitespace.test(string.charAt(index))) {}
  return index;
}

/* harmony default export */ const _trimmedEndIndex = (trimmedEndIndex);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseTrim.js


/** Used to match leading whitespace. */
var reTrimStart = /^\s+/;

/**
 * The base implementation of `_.trim`.
 *
 * @private
 * @param {string} string The string to trim.
 * @returns {string} Returns the trimmed string.
 */
function baseTrim(string) {
  return string
    ? string.slice(0, _trimmedEndIndex(string) + 1).replace(reTrimStart, '')
    : string;
}

/* harmony default export */ const _baseTrim = (baseTrim);

// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
// EXTERNAL MODULE: ../node_modules/lodash-es/isSymbol.js
var isSymbol = __webpack_require__(59660);
;// CONCATENATED MODULE: ../node_modules/lodash-es/toNumber.js




/** Used as references for various `Number` constants. */
var NAN = 0 / 0;

/** Used to detect bad signed hexadecimal string values. */
var reIsBadHex = /^[-+]0x[0-9a-f]+$/i;

/** Used to detect binary string values. */
var reIsBinary = /^0b[01]+$/i;

/** Used to detect octal string values. */
var reIsOctal = /^0o[0-7]+$/i;

/** Built-in method references without a dependency on `root`. */
var freeParseInt = parseInt;

/**
 * Converts `value` to a number.
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to process.
 * @returns {number} Returns the number.
 * @example
 *
 * _.toNumber(3.2);
 * // => 3.2
 *
 * _.toNumber(Number.MIN_VALUE);
 * // => 5e-324
 *
 * _.toNumber(Infinity);
 * // => Infinity
 *
 * _.toNumber('3.2');
 * // => 3.2
 */
function toNumber(value) {
  if (typeof value == 'number') {
    return value;
  }
  if ((0,isSymbol/* default */.Z)(value)) {
    return NAN;
  }
  if ((0,isObject/* default */.Z)(value)) {
    var other = typeof value.valueOf == 'function' ? value.valueOf() : value;
    value = (0,isObject/* default */.Z)(other) ? (other + '') : other;
  }
  if (typeof value != 'string') {
    return value === 0 ? value : +value;
  }
  value = _baseTrim(value);
  var isBinary = reIsBinary.test(value);
  return (isBinary || reIsOctal.test(value))
    ? freeParseInt(value.slice(2), isBinary ? 2 : 8)
    : (reIsBadHex.test(value) ? NAN : +value);
}

/* harmony default export */ const lodash_es_toNumber = (toNumber);

;// CONCATENATED MODULE: ../node_modules/lodash-es/toFinite.js


/** Used as references for various `Number` constants. */
var INFINITY = 1 / 0,
    MAX_INTEGER = 1.7976931348623157e+308;

/**
 * Converts `value` to a finite number.
 *
 * @static
 * @memberOf _
 * @since 4.12.0
 * @category Lang
 * @param {*} value The value to convert.
 * @returns {number} Returns the converted number.
 * @example
 *
 * _.toFinite(3.2);
 * // => 3.2
 *
 * _.toFinite(Number.MIN_VALUE);
 * // => 5e-324
 *
 * _.toFinite(Infinity);
 * // => 1.7976931348623157e+308
 *
 * _.toFinite('3.2');
 * // => 3.2
 */
function toFinite(value) {
  if (!value) {
    return value === 0 ? value : 0;
  }
  value = lodash_es_toNumber(value);
  if (value === INFINITY || value === -INFINITY) {
    var sign = (value < 0 ? -1 : 1);
    return sign * MAX_INTEGER;
  }
  return value === value ? value : 0;
}

/* harmony default export */ const lodash_es_toFinite = (toFinite);


/***/ }),

/***/ 98670:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _toFinite_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(41291);


/**
 * Converts `value` to an integer.
 *
 * **Note:** This method is loosely based on
 * [`ToInteger`](http://www.ecma-international.org/ecma-262/7.0/#sec-tointeger).
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to convert.
 * @returns {number} Returns the converted integer.
 * @example
 *
 * _.toInteger(3.2);
 * // => 3
 *
 * _.toInteger(Number.MIN_VALUE);
 * // => 0
 *
 * _.toInteger(Infinity);
 * // => 1.7976931348623157e+308
 *
 * _.toInteger('3.2');
 * // => 3
 */
function toInteger(value) {
  var result = (0,_toFinite_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value),
      remainder = result % 1;

  return result === result ? (remainder ? result - remainder : result) : 0;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (toInteger);


/***/ }),

/***/ 27061:
/***/ ((module) => {

// shim for using process in browser
var process = module.exports = {};

// cached from whatever global is present so that test runners that stub it
// don't break things.  But we need to wrap it in a try catch in case it is
// wrapped in strict mode code which doesn't define any globals.  It's inside a
// function because try/catches deoptimize in certain engines.

var cachedSetTimeout;
var cachedClearTimeout;

function defaultSetTimout() {
    throw new Error('setTimeout has not been defined');
}
function defaultClearTimeout () {
    throw new Error('clearTimeout has not been defined');
}
(function () {
    try {
        if (typeof setTimeout === 'function') {
            cachedSetTimeout = setTimeout;
        } else {
            cachedSetTimeout = defaultSetTimout;
        }
    } catch (e) {
        cachedSetTimeout = defaultSetTimout;
    }
    try {
        if (typeof clearTimeout === 'function') {
            cachedClearTimeout = clearTimeout;
        } else {
            cachedClearTimeout = defaultClearTimeout;
        }
    } catch (e) {
        cachedClearTimeout = defaultClearTimeout;
    }
} ())
function runTimeout(fun) {
    if (cachedSetTimeout === setTimeout) {
        //normal enviroments in sane situations
        return setTimeout(fun, 0);
    }
    // if setTimeout wasn't available but was latter defined
    if ((cachedSetTimeout === defaultSetTimout || !cachedSetTimeout) && setTimeout) {
        cachedSetTimeout = setTimeout;
        return setTimeout(fun, 0);
    }
    try {
        // when when somebody has screwed with setTimeout but no I.E. maddness
        return cachedSetTimeout(fun, 0);
    } catch(e){
        try {
            // When we are in I.E. but the script has been evaled so I.E. doesn't trust the global object when called normally
            return cachedSetTimeout.call(null, fun, 0);
        } catch(e){
            // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error
            return cachedSetTimeout.call(this, fun, 0);
        }
    }


}
function runClearTimeout(marker) {
    if (cachedClearTimeout === clearTimeout) {
        //normal enviroments in sane situations
        return clearTimeout(marker);
    }
    // if clearTimeout wasn't available but was latter defined
    if ((cachedClearTimeout === defaultClearTimeout || !cachedClearTimeout) && clearTimeout) {
        cachedClearTimeout = clearTimeout;
        return clearTimeout(marker);
    }
    try {
        // when when somebody has screwed with setTimeout but no I.E. maddness
        return cachedClearTimeout(marker);
    } catch (e){
        try {
            // When we are in I.E. but the script has been evaled so I.E. doesn't  trust the global object when called normally
            return cachedClearTimeout.call(null, marker);
        } catch (e){
            // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error.
            // Some versions of I.E. have different rules for clearTimeout vs setTimeout
            return cachedClearTimeout.call(this, marker);
        }
    }



}
var queue = [];
var draining = false;
var currentQueue;
var queueIndex = -1;

function cleanUpNextTick() {
    if (!draining || !currentQueue) {
        return;
    }
    draining = false;
    if (currentQueue.length) {
        queue = currentQueue.concat(queue);
    } else {
        queueIndex = -1;
    }
    if (queue.length) {
        drainQueue();
    }
}

function drainQueue() {
    if (draining) {
        return;
    }
    var timeout = runTimeout(cleanUpNextTick);
    draining = true;

    var len = queue.length;
    while(len) {
        currentQueue = queue;
        queue = [];
        while (++queueIndex < len) {
            if (currentQueue) {
                currentQueue[queueIndex].run();
            }
        }
        queueIndex = -1;
        len = queue.length;
    }
    currentQueue = null;
    draining = false;
    runClearTimeout(timeout);
}

process.nextTick = function (fun) {
    var args = new Array(arguments.length - 1);
    if (arguments.length > 1) {
        for (var i = 1; i < arguments.length; i++) {
            args[i - 1] = arguments[i];
        }
    }
    queue.push(new Item(fun, args));
    if (queue.length === 1 && !draining) {
        runTimeout(drainQueue);
    }
};

// v8 likes predictible objects
function Item(fun, array) {
    this.fun = fun;
    this.array = array;
}
Item.prototype.run = function () {
    this.fun.apply(null, this.array);
};
process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];
process.version = ''; // empty string to avoid regexp issues
process.versions = {};

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;
process.prependListener = noop;
process.prependOnceListener = noop;

process.listeners = function (name) { return [] }

process.binding = function (name) {
    throw new Error('process.binding is not supported');
};

process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};
process.umask = function() { return 0; };


/***/ }),

/***/ 97770:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

/*---------------------------------------------------------------------------------------------
 *  Copyright (c) Microsoft Corporation. All rights reserved.
 *  Licensed under the MIT License. See License.txt in the project root for license information.
 *--------------------------------------------------------------------------------------------*/
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.CancellationTokenSource = exports.CancellationToken = void 0;
const ral_1 = __webpack_require__(48094);
const Is = __webpack_require__(78472);
const events_1 = __webpack_require__(345);
var CancellationToken;
(function (CancellationToken) {
    CancellationToken.None = Object.freeze({
        isCancellationRequested: false,
        onCancellationRequested: events_1.Event.None
    });
    CancellationToken.Cancelled = Object.freeze({
        isCancellationRequested: true,
        onCancellationRequested: events_1.Event.None
    });
    function is(value) {
        const candidate = value;
        return candidate && (candidate === CancellationToken.None
            || candidate === CancellationToken.Cancelled
            || (Is.boolean(candidate.isCancellationRequested) && !!candidate.onCancellationRequested));
    }
    CancellationToken.is = is;
})(CancellationToken || (exports.CancellationToken = CancellationToken = {}));
const shortcutEvent = Object.freeze(function (callback, context) {
    const handle = (0, ral_1.default)().timer.setTimeout(callback.bind(context), 0);
    return { dispose() { handle.dispose(); } };
});
class MutableToken {
    constructor() {
        this._isCancelled = false;
    }
    cancel() {
        if (!this._isCancelled) {
            this._isCancelled = true;
            if (this._emitter) {
                this._emitter.fire(undefined);
                this.dispose();
            }
        }
    }
    get isCancellationRequested() {
        return this._isCancelled;
    }
    get onCancellationRequested() {
        if (this._isCancelled) {
            return shortcutEvent;
        }
        if (!this._emitter) {
            this._emitter = new events_1.Emitter();
        }
        return this._emitter.event;
    }
    dispose() {
        if (this._emitter) {
            this._emitter.dispose();
            this._emitter = undefined;
        }
    }
}
class CancellationTokenSource {
    get token() {
        if (!this._token) {
            // be lazy and create the token only when
            // actually needed
            this._token = new MutableToken();
        }
        return this._token;
    }
    cancel() {
        if (!this._token) {
            // save an object by returning the default
            // cancelled token when cancellation happens
            // before someone asks for the token
            this._token = CancellationToken.Cancelled;
        }
        else {
            this._token.cancel();
        }
    }
    dispose() {
        if (!this._token) {
            // ensure to initialize with an empty token if we had none
            this._token = CancellationToken.None;
        }
        else if (this._token instanceof MutableToken) {
            // actually dispose
            this._token.dispose();
        }
    }
}
exports.CancellationTokenSource = CancellationTokenSource;


/***/ }),

/***/ 345:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

/* --------------------------------------------------------------------------------------------
 * Copyright (c) Microsoft Corporation. All rights reserved.
 * Licensed under the MIT License. See License.txt in the project root for license information.
 * ------------------------------------------------------------------------------------------ */
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.Emitter = exports.Event = void 0;
const ral_1 = __webpack_require__(48094);
var Event;
(function (Event) {
    const _disposable = { dispose() { } };
    Event.None = function () { return _disposable; };
})(Event || (exports.Event = Event = {}));
class CallbackList {
    add(callback, context = null, bucket) {
        if (!this._callbacks) {
            this._callbacks = [];
            this._contexts = [];
        }
        this._callbacks.push(callback);
        this._contexts.push(context);
        if (Array.isArray(bucket)) {
            bucket.push({ dispose: () => this.remove(callback, context) });
        }
    }
    remove(callback, context = null) {
        if (!this._callbacks) {
            return;
        }
        let foundCallbackWithDifferentContext = false;
        for (let i = 0, len = this._callbacks.length; i < len; i++) {
            if (this._callbacks[i] === callback) {
                if (this._contexts[i] === context) {
                    // callback & context match => remove it
                    this._callbacks.splice(i, 1);
                    this._contexts.splice(i, 1);
                    return;
                }
                else {
                    foundCallbackWithDifferentContext = true;
                }
            }
        }
        if (foundCallbackWithDifferentContext) {
            throw new Error('When adding a listener with a context, you should remove it with the same context');
        }
    }
    invoke(...args) {
        if (!this._callbacks) {
            return [];
        }
        const ret = [], callbacks = this._callbacks.slice(0), contexts = this._contexts.slice(0);
        for (let i = 0, len = callbacks.length; i < len; i++) {
            try {
                ret.push(callbacks[i].apply(contexts[i], args));
            }
            catch (e) {
                // eslint-disable-next-line no-console
                (0, ral_1.default)().console.error(e);
            }
        }
        return ret;
    }
    isEmpty() {
        return !this._callbacks || this._callbacks.length === 0;
    }
    dispose() {
        this._callbacks = undefined;
        this._contexts = undefined;
    }
}
class Emitter {
    constructor(_options) {
        this._options = _options;
    }
    /**
     * For the public to allow to subscribe
     * to events from this Emitter
     */
    get event() {
        if (!this._event) {
            this._event = (listener, thisArgs, disposables) => {
                if (!this._callbacks) {
                    this._callbacks = new CallbackList();
                }
                if (this._options && this._options.onFirstListenerAdd && this._callbacks.isEmpty()) {
                    this._options.onFirstListenerAdd(this);
                }
                this._callbacks.add(listener, thisArgs);
                const result = {
                    dispose: () => {
                        if (!this._callbacks) {
                            // disposable is disposed after emitter is disposed.
                            return;
                        }
                        this._callbacks.remove(listener, thisArgs);
                        result.dispose = Emitter._noop;
                        if (this._options && this._options.onLastListenerRemove && this._callbacks.isEmpty()) {
                            this._options.onLastListenerRemove(this);
                        }
                    }
                };
                if (Array.isArray(disposables)) {
                    disposables.push(result);
                }
                return result;
            };
        }
        return this._event;
    }
    /**
     * To be kept private to fire an event to
     * subscribers
     */
    fire(event) {
        if (this._callbacks) {
            this._callbacks.invoke.call(this._callbacks, event);
        }
    }
    dispose() {
        if (this._callbacks) {
            this._callbacks.dispose();
            this._callbacks = undefined;
        }
    }
}
exports.Emitter = Emitter;
Emitter._noop = function () { };


/***/ }),

/***/ 78472:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

/* --------------------------------------------------------------------------------------------
 * Copyright (c) Microsoft Corporation. All rights reserved.
 * Licensed under the MIT License. See License.txt in the project root for license information.
 * ------------------------------------------------------------------------------------------ */
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.stringArray = exports.array = exports.func = exports.error = exports.number = exports.string = exports.boolean = void 0;
function boolean(value) {
    return value === true || value === false;
}
exports.boolean = boolean;
function string(value) {
    return typeof value === 'string' || value instanceof String;
}
exports.string = string;
function number(value) {
    return typeof value === 'number' || value instanceof Number;
}
exports.number = number;
function error(value) {
    return value instanceof Error;
}
exports.error = error;
function func(value) {
    return typeof value === 'function';
}
exports.func = func;
function array(value) {
    return Array.isArray(value);
}
exports.array = array;
function stringArray(value) {
    return array(value) && value.every(elem => string(elem));
}
exports.stringArray = stringArray;


/***/ }),

/***/ 48094:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

/* --------------------------------------------------------------------------------------------
 * Copyright (c) Microsoft Corporation. All rights reserved.
 * Licensed under the MIT License. See License.txt in the project root for license information.
 * ------------------------------------------------------------------------------------------ */
Object.defineProperty(exports, "__esModule", ({ value: true }));
let _ral;
function RAL() {
    if (_ral === undefined) {
        throw new Error(`No runtime abstraction layer installed`);
    }
    return _ral;
}
(function (RAL) {
    function install(ral) {
        if (ral === undefined) {
            throw new Error(`No runtime abstraction layer provided`);
        }
        _ral = ral;
    }
    RAL.install = install;
})(RAL || (RAL = {}));
exports["default"] = RAL;


/***/ }),

/***/ 37943:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   c: () => (/* binding */ Utils),
/* harmony export */   o: () => (/* binding */ URI)
/* harmony export */ });
/* provided dependency */ var process = __webpack_require__(27061);
var LIB;(()=>{"use strict";var t={470:t=>{function e(t){if("string"!=typeof t)throw new TypeError("Path must be a string. Received "+JSON.stringify(t))}function r(t,e){for(var r,n="",i=0,o=-1,s=0,h=0;h<=t.length;++h){if(h<t.length)r=t.charCodeAt(h);else{if(47===r)break;r=47}if(47===r){if(o===h-1||1===s);else if(o!==h-1&&2===s){if(n.length<2||2!==i||46!==n.charCodeAt(n.length-1)||46!==n.charCodeAt(n.length-2))if(n.length>2){var a=n.lastIndexOf("/");if(a!==n.length-1){-1===a?(n="",i=0):i=(n=n.slice(0,a)).length-1-n.lastIndexOf("/"),o=h,s=0;continue}}else if(2===n.length||1===n.length){n="",i=0,o=h,s=0;continue}e&&(n.length>0?n+="/..":n="..",i=2)}else n.length>0?n+="/"+t.slice(o+1,h):n=t.slice(o+1,h),i=h-o-1;o=h,s=0}else 46===r&&-1!==s?++s:s=-1}return n}var n={resolve:function(){for(var t,n="",i=!1,o=arguments.length-1;o>=-1&&!i;o--){var s;o>=0?s=arguments[o]:(void 0===t&&(t=process.cwd()),s=t),e(s),0!==s.length&&(n=s+"/"+n,i=47===s.charCodeAt(0))}return n=r(n,!i),i?n.length>0?"/"+n:"/":n.length>0?n:"."},normalize:function(t){if(e(t),0===t.length)return".";var n=47===t.charCodeAt(0),i=47===t.charCodeAt(t.length-1);return 0!==(t=r(t,!n)).length||n||(t="."),t.length>0&&i&&(t+="/"),n?"/"+t:t},isAbsolute:function(t){return e(t),t.length>0&&47===t.charCodeAt(0)},join:function(){if(0===arguments.length)return".";for(var t,r=0;r<arguments.length;++r){var i=arguments[r];e(i),i.length>0&&(void 0===t?t=i:t+="/"+i)}return void 0===t?".":n.normalize(t)},relative:function(t,r){if(e(t),e(r),t===r)return"";if((t=n.resolve(t))===(r=n.resolve(r)))return"";for(var i=1;i<t.length&&47===t.charCodeAt(i);++i);for(var o=t.length,s=o-i,h=1;h<r.length&&47===r.charCodeAt(h);++h);for(var a=r.length-h,c=s<a?s:a,f=-1,u=0;u<=c;++u){if(u===c){if(a>c){if(47===r.charCodeAt(h+u))return r.slice(h+u+1);if(0===u)return r.slice(h+u)}else s>c&&(47===t.charCodeAt(i+u)?f=u:0===u&&(f=0));break}var l=t.charCodeAt(i+u);if(l!==r.charCodeAt(h+u))break;47===l&&(f=u)}var g="";for(u=i+f+1;u<=o;++u)u!==o&&47!==t.charCodeAt(u)||(0===g.length?g+="..":g+="/..");return g.length>0?g+r.slice(h+f):(h+=f,47===r.charCodeAt(h)&&++h,r.slice(h))},_makeLong:function(t){return t},dirname:function(t){if(e(t),0===t.length)return".";for(var r=t.charCodeAt(0),n=47===r,i=-1,o=!0,s=t.length-1;s>=1;--s)if(47===(r=t.charCodeAt(s))){if(!o){i=s;break}}else o=!1;return-1===i?n?"/":".":n&&1===i?"//":t.slice(0,i)},basename:function(t,r){if(void 0!==r&&"string"!=typeof r)throw new TypeError('"ext" argument must be a string');e(t);var n,i=0,o=-1,s=!0;if(void 0!==r&&r.length>0&&r.length<=t.length){if(r.length===t.length&&r===t)return"";var h=r.length-1,a=-1;for(n=t.length-1;n>=0;--n){var c=t.charCodeAt(n);if(47===c){if(!s){i=n+1;break}}else-1===a&&(s=!1,a=n+1),h>=0&&(c===r.charCodeAt(h)?-1==--h&&(o=n):(h=-1,o=a))}return i===o?o=a:-1===o&&(o=t.length),t.slice(i,o)}for(n=t.length-1;n>=0;--n)if(47===t.charCodeAt(n)){if(!s){i=n+1;break}}else-1===o&&(s=!1,o=n+1);return-1===o?"":t.slice(i,o)},extname:function(t){e(t);for(var r=-1,n=0,i=-1,o=!0,s=0,h=t.length-1;h>=0;--h){var a=t.charCodeAt(h);if(47!==a)-1===i&&(o=!1,i=h+1),46===a?-1===r?r=h:1!==s&&(s=1):-1!==r&&(s=-1);else if(!o){n=h+1;break}}return-1===r||-1===i||0===s||1===s&&r===i-1&&r===n+1?"":t.slice(r,i)},format:function(t){if(null===t||"object"!=typeof t)throw new TypeError('The "pathObject" argument must be of type Object. Received type '+typeof t);return function(t,e){var r=e.dir||e.root,n=e.base||(e.name||"")+(e.ext||"");return r?r===e.root?r+n:r+"/"+n:n}(0,t)},parse:function(t){e(t);var r={root:"",dir:"",base:"",ext:"",name:""};if(0===t.length)return r;var n,i=t.charCodeAt(0),o=47===i;o?(r.root="/",n=1):n=0;for(var s=-1,h=0,a=-1,c=!0,f=t.length-1,u=0;f>=n;--f)if(47!==(i=t.charCodeAt(f)))-1===a&&(c=!1,a=f+1),46===i?-1===s?s=f:1!==u&&(u=1):-1!==s&&(u=-1);else if(!c){h=f+1;break}return-1===s||-1===a||0===u||1===u&&s===a-1&&s===h+1?-1!==a&&(r.base=r.name=0===h&&o?t.slice(1,a):t.slice(h,a)):(0===h&&o?(r.name=t.slice(1,s),r.base=t.slice(1,a)):(r.name=t.slice(h,s),r.base=t.slice(h,a)),r.ext=t.slice(s,a)),h>0?r.dir=t.slice(0,h-1):o&&(r.dir="/"),r},sep:"/",delimiter:":",win32:null,posix:null};n.posix=n,t.exports=n}},e={};function r(n){var i=e[n];if(void 0!==i)return i.exports;var o=e[n]={exports:{}};return t[n](o,o.exports,r),o.exports}r.d=(t,e)=>{for(var n in e)r.o(e,n)&&!r.o(t,n)&&Object.defineProperty(t,n,{enumerable:!0,get:e[n]})},r.o=(t,e)=>Object.prototype.hasOwnProperty.call(t,e),r.r=t=>{"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(t,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(t,"__esModule",{value:!0})};var n={};(()=>{let t;if(r.r(n),r.d(n,{URI:()=>f,Utils:()=>P}),"object"==typeof process)t="win32"===process.platform;else if("object"==typeof navigator){let e=navigator.userAgent;t=e.indexOf("Windows")>=0}const e=/^\w[\w\d+.-]*$/,i=/^\//,o=/^\/\//;function s(t,r){if(!t.scheme&&r)throw new Error(`[UriError]: Scheme is missing: {scheme: "", authority: "${t.authority}", path: "${t.path}", query: "${t.query}", fragment: "${t.fragment}"}`);if(t.scheme&&!e.test(t.scheme))throw new Error("[UriError]: Scheme contains illegal characters.");if(t.path)if(t.authority){if(!i.test(t.path))throw new Error('[UriError]: If a URI contains an authority component, then the path component must either be empty or begin with a slash ("/") character')}else if(o.test(t.path))throw new Error('[UriError]: If a URI does not contain an authority component, then the path cannot begin with two slash characters ("//")')}const h="",a="/",c=/^(([^:/?#]+?):)?(\/\/([^/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?/;class f{static isUri(t){return t instanceof f||!!t&&"string"==typeof t.authority&&"string"==typeof t.fragment&&"string"==typeof t.path&&"string"==typeof t.query&&"string"==typeof t.scheme&&"string"==typeof t.fsPath&&"function"==typeof t.with&&"function"==typeof t.toString}scheme;authority;path;query;fragment;constructor(t,e,r,n,i,o=!1){"object"==typeof t?(this.scheme=t.scheme||h,this.authority=t.authority||h,this.path=t.path||h,this.query=t.query||h,this.fragment=t.fragment||h):(this.scheme=function(t,e){return t||e?t:"file"}(t,o),this.authority=e||h,this.path=function(t,e){switch(t){case"https":case"http":case"file":e?e[0]!==a&&(e=a+e):e=a}return e}(this.scheme,r||h),this.query=n||h,this.fragment=i||h,s(this,o))}get fsPath(){return m(this,!1)}with(t){if(!t)return this;let{scheme:e,authority:r,path:n,query:i,fragment:o}=t;return void 0===e?e=this.scheme:null===e&&(e=h),void 0===r?r=this.authority:null===r&&(r=h),void 0===n?n=this.path:null===n&&(n=h),void 0===i?i=this.query:null===i&&(i=h),void 0===o?o=this.fragment:null===o&&(o=h),e===this.scheme&&r===this.authority&&n===this.path&&i===this.query&&o===this.fragment?this:new l(e,r,n,i,o)}static parse(t,e=!1){const r=c.exec(t);return r?new l(r[2]||h,C(r[4]||h),C(r[5]||h),C(r[7]||h),C(r[9]||h),e):new l(h,h,h,h,h)}static file(e){let r=h;if(t&&(e=e.replace(/\\/g,a)),e[0]===a&&e[1]===a){const t=e.indexOf(a,2);-1===t?(r=e.substring(2),e=a):(r=e.substring(2,t),e=e.substring(t)||a)}return new l("file",r,e,h,h)}static from(t){const e=new l(t.scheme,t.authority,t.path,t.query,t.fragment);return s(e,!0),e}toString(t=!1){return y(this,t)}toJSON(){return this}static revive(t){if(t){if(t instanceof f)return t;{const e=new l(t);return e._formatted=t.external,e._fsPath=t._sep===u?t.fsPath:null,e}}return t}}const u=t?1:void 0;class l extends f{_formatted=null;_fsPath=null;get fsPath(){return this._fsPath||(this._fsPath=m(this,!1)),this._fsPath}toString(t=!1){return t?y(this,!0):(this._formatted||(this._formatted=y(this,!1)),this._formatted)}toJSON(){const t={$mid:1};return this._fsPath&&(t.fsPath=this._fsPath,t._sep=u),this._formatted&&(t.external=this._formatted),this.path&&(t.path=this.path),this.scheme&&(t.scheme=this.scheme),this.authority&&(t.authority=this.authority),this.query&&(t.query=this.query),this.fragment&&(t.fragment=this.fragment),t}}const g={58:"%3A",47:"%2F",63:"%3F",35:"%23",91:"%5B",93:"%5D",64:"%40",33:"%21",36:"%24",38:"%26",39:"%27",40:"%28",41:"%29",42:"%2A",43:"%2B",44:"%2C",59:"%3B",61:"%3D",32:"%20"};function d(t,e,r){let n,i=-1;for(let o=0;o<t.length;o++){const s=t.charCodeAt(o);if(s>=97&&s<=122||s>=65&&s<=90||s>=48&&s<=57||45===s||46===s||95===s||126===s||e&&47===s||r&&91===s||r&&93===s||r&&58===s)-1!==i&&(n+=encodeURIComponent(t.substring(i,o)),i=-1),void 0!==n&&(n+=t.charAt(o));else{void 0===n&&(n=t.substr(0,o));const e=g[s];void 0!==e?(-1!==i&&(n+=encodeURIComponent(t.substring(i,o)),i=-1),n+=e):-1===i&&(i=o)}}return-1!==i&&(n+=encodeURIComponent(t.substring(i))),void 0!==n?n:t}function p(t){let e;for(let r=0;r<t.length;r++){const n=t.charCodeAt(r);35===n||63===n?(void 0===e&&(e=t.substr(0,r)),e+=g[n]):void 0!==e&&(e+=t[r])}return void 0!==e?e:t}function m(e,r){let n;return n=e.authority&&e.path.length>1&&"file"===e.scheme?`//${e.authority}${e.path}`:47===e.path.charCodeAt(0)&&(e.path.charCodeAt(1)>=65&&e.path.charCodeAt(1)<=90||e.path.charCodeAt(1)>=97&&e.path.charCodeAt(1)<=122)&&58===e.path.charCodeAt(2)?r?e.path.substr(1):e.path[1].toLowerCase()+e.path.substr(2):e.path,t&&(n=n.replace(/\//g,"\\")),n}function y(t,e){const r=e?p:d;let n="",{scheme:i,authority:o,path:s,query:h,fragment:c}=t;if(i&&(n+=i,n+=":"),(o||"file"===i)&&(n+=a,n+=a),o){let t=o.indexOf("@");if(-1!==t){const e=o.substr(0,t);o=o.substr(t+1),t=e.lastIndexOf(":"),-1===t?n+=r(e,!1,!1):(n+=r(e.substr(0,t),!1,!1),n+=":",n+=r(e.substr(t+1),!1,!0)),n+="@"}o=o.toLowerCase(),t=o.lastIndexOf(":"),-1===t?n+=r(o,!1,!0):(n+=r(o.substr(0,t),!1,!0),n+=o.substr(t))}if(s){if(s.length>=3&&47===s.charCodeAt(0)&&58===s.charCodeAt(2)){const t=s.charCodeAt(1);t>=65&&t<=90&&(s=`/${String.fromCharCode(t+32)}:${s.substr(3)}`)}else if(s.length>=2&&58===s.charCodeAt(1)){const t=s.charCodeAt(0);t>=65&&t<=90&&(s=`${String.fromCharCode(t+32)}:${s.substr(2)}`)}n+=r(s,!0,!1)}return h&&(n+="?",n+=r(h,!1,!1)),c&&(n+="#",n+=e?c:d(c,!1,!1)),n}function v(t){try{return decodeURIComponent(t)}catch{return t.length>3?t.substr(0,3)+v(t.substr(3)):t}}const b=/(%[0-9A-Za-z][0-9A-Za-z])+/g;function C(t){return t.match(b)?t.replace(b,(t=>v(t))):t}var A=r(470);const w=A.posix||A,x="/";var P;!function(t){t.joinPath=function(t,...e){return t.with({path:w.join(t.path,...e)})},t.resolvePath=function(t,...e){let r=t.path,n=!1;r[0]!==x&&(r=x+r,n=!0);let i=w.resolve(r,...e);return n&&i[0]===x&&!t.authority&&(i=i.substring(1)),t.with({path:i})},t.dirname=function(t){if(0===t.path.length||t.path===x)return t;let e=w.dirname(t.path);return 1===e.length&&46===e.charCodeAt(0)&&(e=""),t.with({path:e})},t.basename=function(t){return w.basename(t.path)},t.extname=function(t){return w.extname(t.path)}}(P||(P={}))})(),LIB=n})();const{URI,Utils}=LIB;
//# sourceMappingURL=index.mjs.map

/***/ })

}]);
//# sourceMappingURL=3197.132cf892d4ef38649b32.js.map?v=132cf892d4ef38649b32