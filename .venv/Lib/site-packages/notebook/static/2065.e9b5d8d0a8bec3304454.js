"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2065],{

/***/ 32778:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __exportStar = (this && this.__exportStar) || function(m, exports) {
    for (var p in m) if (p !== "default" && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
var __read = (this && this.__read) || function (o, n) {
    var m = typeof Symbol === "function" && o[Symbol.iterator];
    if (!m) return o;
    var i = m.call(o), r, ar = [], e;
    try {
        while ((n === void 0 || n-- > 0) && !(r = i.next()).done) ar.push(r.value);
    }
    catch (error) { e = { error: error }; }
    finally {
        try {
            if (r && !r.done && (m = i["return"])) m.call(i);
        }
        finally { if (e) throw e.error; }
    }
    return ar;
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AddCSS = exports.CHTMLFontData = void 0;
var FontData_js_1 = __webpack_require__(5549);
var Usage_js_1 = __webpack_require__(53250);
var lengths_js_1 = __webpack_require__(56780);
__exportStar(__webpack_require__(5549), exports);
var CHTMLFontData = (function (_super) {
    __extends(CHTMLFontData, _super);
    function CHTMLFontData() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.charUsage = new Usage_js_1.Usage();
        _this.delimUsage = new Usage_js_1.Usage();
        return _this;
    }
    CHTMLFontData.charOptions = function (font, n) {
        return _super.charOptions.call(this, font, n);
    };
    CHTMLFontData.prototype.adaptiveCSS = function (adapt) {
        this.options.adaptiveCSS = adapt;
    };
    CHTMLFontData.prototype.clearCache = function () {
        if (this.options.adaptiveCSS) {
            this.charUsage.clear();
            this.delimUsage.clear();
        }
    };
    CHTMLFontData.prototype.createVariant = function (name, inherit, link) {
        if (inherit === void 0) { inherit = null; }
        if (link === void 0) { link = null; }
        _super.prototype.createVariant.call(this, name, inherit, link);
        var CLASS = this.constructor;
        this.variant[name].classes = CLASS.defaultVariantClasses[name];
        this.variant[name].letter = CLASS.defaultVariantLetters[name];
    };
    CHTMLFontData.prototype.defineChars = function (name, chars) {
        var e_1, _a;
        _super.prototype.defineChars.call(this, name, chars);
        var letter = this.variant[name].letter;
        try {
            for (var _b = __values(Object.keys(chars)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var n = _c.value;
                var options = CHTMLFontData.charOptions(chars, parseInt(n));
                if (options.f === undefined) {
                    options.f = letter;
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
    };
    Object.defineProperty(CHTMLFontData.prototype, "styles", {
        get: function () {
            var CLASS = this.constructor;
            var styles = __assign({}, CLASS.defaultStyles);
            this.addFontURLs(styles, CLASS.defaultFonts, this.options.fontURL);
            if (this.options.adaptiveCSS) {
                this.updateStyles(styles);
            }
            else {
                this.allStyles(styles);
            }
            return styles;
        },
        enumerable: false,
        configurable: true
    });
    CHTMLFontData.prototype.updateStyles = function (styles) {
        var e_2, _a, e_3, _b;
        try {
            for (var _c = __values(this.delimUsage.update()), _d = _c.next(); !_d.done; _d = _c.next()) {
                var N = _d.value;
                this.addDelimiterStyles(styles, N, this.delimiters[N]);
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_2) throw e_2.error; }
        }
        try {
            for (var _e = __values(this.charUsage.update()), _f = _e.next(); !_f.done; _f = _e.next()) {
                var _g = __read(_f.value, 2), name_1 = _g[0], N = _g[1];
                var variant = this.variant[name_1];
                this.addCharStyles(styles, variant.letter, N, variant.chars[N]);
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
            }
            finally { if (e_3) throw e_3.error; }
        }
        return styles;
    };
    CHTMLFontData.prototype.allStyles = function (styles) {
        var e_4, _a, e_5, _b, e_6, _c;
        try {
            for (var _d = __values(Object.keys(this.delimiters)), _e = _d.next(); !_e.done; _e = _d.next()) {
                var n = _e.value;
                var N = parseInt(n);
                this.addDelimiterStyles(styles, N, this.delimiters[N]);
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_e && !_e.done && (_a = _d.return)) _a.call(_d);
            }
            finally { if (e_4) throw e_4.error; }
        }
        try {
            for (var _f = __values(Object.keys(this.variant)), _g = _f.next(); !_g.done; _g = _f.next()) {
                var name_2 = _g.value;
                var variant = this.variant[name_2];
                var vletter = variant.letter;
                try {
                    for (var _h = (e_6 = void 0, __values(Object.keys(variant.chars))), _j = _h.next(); !_j.done; _j = _h.next()) {
                        var n = _j.value;
                        var N = parseInt(n);
                        var char = variant.chars[N];
                        if ((char[3] || {}).smp)
                            continue;
                        if (char.length < 4) {
                            char[3] = {};
                        }
                        this.addCharStyles(styles, vletter, N, char);
                    }
                }
                catch (e_6_1) { e_6 = { error: e_6_1 }; }
                finally {
                    try {
                        if (_j && !_j.done && (_c = _h.return)) _c.call(_h);
                    }
                    finally { if (e_6) throw e_6.error; }
                }
            }
        }
        catch (e_5_1) { e_5 = { error: e_5_1 }; }
        finally {
            try {
                if (_g && !_g.done && (_b = _f.return)) _b.call(_f);
            }
            finally { if (e_5) throw e_5.error; }
        }
    };
    CHTMLFontData.prototype.addFontURLs = function (styles, fonts, url) {
        var e_7, _a;
        try {
            for (var _b = __values(Object.keys(fonts)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var name_3 = _c.value;
                var font = __assign({}, fonts[name_3]);
                font.src = font.src.replace(/%%URL%%/, url);
                styles[name_3] = font;
            }
        }
        catch (e_7_1) { e_7 = { error: e_7_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_7) throw e_7.error; }
        }
    };
    CHTMLFontData.prototype.addDelimiterStyles = function (styles, n, data) {
        var c = this.charSelector(n);
        if (data.c && data.c !== n) {
            c = this.charSelector(data.c);
            styles['.mjx-stretched mjx-c' + c + '::before'] = {
                content: this.charContent(data.c)
            };
        }
        if (!data.stretch)
            return;
        if (data.dir === 1) {
            this.addDelimiterVStyles(styles, c, data);
        }
        else {
            this.addDelimiterHStyles(styles, c, data);
        }
    };
    CHTMLFontData.prototype.addDelimiterVStyles = function (styles, c, data) {
        var HDW = data.HDW;
        var _a = __read(data.stretch, 4), beg = _a[0], ext = _a[1], end = _a[2], mid = _a[3];
        var Hb = this.addDelimiterVPart(styles, c, 'beg', beg, HDW);
        this.addDelimiterVPart(styles, c, 'ext', ext, HDW);
        var He = this.addDelimiterVPart(styles, c, 'end', end, HDW);
        var css = {};
        if (mid) {
            var Hm = this.addDelimiterVPart(styles, c, 'mid', mid, HDW);
            css.height = '50%';
            styles['mjx-stretchy-v' + c + ' > mjx-mid'] = {
                'margin-top': this.em(-Hm / 2),
                'margin-bottom': this.em(-Hm / 2)
            };
        }
        if (Hb) {
            css['border-top-width'] = this.em0(Hb - .03);
        }
        if (He) {
            css['border-bottom-width'] = this.em0(He - .03);
            styles['mjx-stretchy-v' + c + ' > mjx-end'] = { 'margin-top': this.em(-He) };
        }
        if (Object.keys(css).length) {
            styles['mjx-stretchy-v' + c + ' > mjx-ext'] = css;
        }
    };
    CHTMLFontData.prototype.addDelimiterVPart = function (styles, c, part, n, HDW) {
        if (!n)
            return 0;
        var data = this.getDelimiterData(n);
        var dw = (HDW[2] - data[2]) / 2;
        var css = { content: this.charContent(n) };
        if (part !== 'ext') {
            css.padding = this.padding(data, dw);
        }
        else {
            css.width = this.em0(HDW[2]);
            if (dw) {
                css['padding-left'] = this.em0(dw);
            }
        }
        styles['mjx-stretchy-v' + c + ' mjx-' + part + ' mjx-c::before'] = css;
        return data[0] + data[1];
    };
    CHTMLFontData.prototype.addDelimiterHStyles = function (styles, c, data) {
        var _a = __read(data.stretch, 4), beg = _a[0], ext = _a[1], end = _a[2], mid = _a[3];
        var HDW = data.HDW;
        this.addDelimiterHPart(styles, c, 'beg', beg, HDW);
        this.addDelimiterHPart(styles, c, 'ext', ext, HDW);
        this.addDelimiterHPart(styles, c, 'end', end, HDW);
        if (mid) {
            this.addDelimiterHPart(styles, c, 'mid', mid, HDW);
            styles['mjx-stretchy-h' + c + ' > mjx-ext'] = { width: '50%' };
        }
    };
    CHTMLFontData.prototype.addDelimiterHPart = function (styles, c, part, n, HDW) {
        if (!n)
            return;
        var data = this.getDelimiterData(n);
        var options = data[3];
        var css = { content: (options && options.c ? '"' + options.c + '"' : this.charContent(n)) };
        css.padding = this.padding(HDW, 0, -HDW[2]);
        styles['mjx-stretchy-h' + c + ' mjx-' + part + ' mjx-c::before'] = css;
    };
    CHTMLFontData.prototype.addCharStyles = function (styles, vletter, n, data) {
        var options = data[3];
        var letter = (options.f !== undefined ? options.f : vletter);
        var selector = 'mjx-c' + this.charSelector(n) + (letter ? '.TEX-' + letter : '');
        styles[selector + '::before'] = {
            padding: this.padding(data, 0, options.ic || 0),
            content: (options.c != null ? '"' + options.c + '"' : this.charContent(n))
        };
    };
    CHTMLFontData.prototype.getDelimiterData = function (n) {
        return this.getChar('-smallop', n);
    };
    CHTMLFontData.prototype.em = function (n) {
        return (0, lengths_js_1.em)(n);
    };
    CHTMLFontData.prototype.em0 = function (n) {
        return (0, lengths_js_1.em)(Math.max(0, n));
    };
    CHTMLFontData.prototype.padding = function (_a, dw, ic) {
        var _b = __read(_a, 3), h = _b[0], d = _b[1], w = _b[2];
        if (dw === void 0) { dw = 0; }
        if (ic === void 0) { ic = 0; }
        return [h, w + ic, d, dw].map(this.em0).join(' ');
    };
    CHTMLFontData.prototype.charContent = function (n) {
        return '"' + (n >= 0x20 && n <= 0x7E && n !== 0x22 && n !== 0x27 && n !== 0x5C ?
            String.fromCharCode(n) : '\\' + n.toString(16).toUpperCase()) + '"';
    };
    CHTMLFontData.prototype.charSelector = function (n) {
        return '.mjx-c' + n.toString(16).toUpperCase();
    };
    CHTMLFontData.OPTIONS = __assign(__assign({}, FontData_js_1.FontData.OPTIONS), { fontURL: 'js/output/chtml/fonts/tex-woff-v2' });
    CHTMLFontData.JAX = 'CHTML';
    CHTMLFontData.defaultVariantClasses = {};
    CHTMLFontData.defaultVariantLetters = {};
    CHTMLFontData.defaultStyles = {
        'mjx-c::before': {
            display: 'block',
            width: 0
        }
    };
    CHTMLFontData.defaultFonts = {
        '@font-face /* 0 */': {
            'font-family': 'MJXZERO',
            src: 'url("%%URL%%/MathJax_Zero.woff") format("woff")'
        }
    };
    return CHTMLFontData;
}(FontData_js_1.FontData));
exports.CHTMLFontData = CHTMLFontData;
function AddCSS(font, options) {
    var e_8, _a;
    try {
        for (var _b = __values(Object.keys(options)), _c = _b.next(); !_c.done; _c = _b.next()) {
            var c = _c.value;
            var n = parseInt(c);
            Object.assign(FontData_js_1.FontData.charOptions(font, n), options[n]);
        }
    }
    catch (e_8_1) { e_8 = { error: e_8_1 }; }
    finally {
        try {
            if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
        }
        finally { if (e_8) throw e_8.error; }
    }
    return font;
}
exports.AddCSS = AddCSS;
//# sourceMappingURL=FontData.js.map

/***/ }),

/***/ 53250:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.Usage = void 0;
var Usage = (function () {
    function Usage() {
        this.used = new Set();
        this.needsUpdate = [];
    }
    Usage.prototype.add = function (item) {
        var name = JSON.stringify(item);
        if (!this.used.has(name)) {
            this.needsUpdate.push(item);
        }
        this.used.add(name);
    };
    Usage.prototype.has = function (item) {
        return this.used.has(JSON.stringify(item));
    };
    Usage.prototype.clear = function () {
        this.used.clear();
        this.needsUpdate = [];
    };
    Usage.prototype.update = function () {
        var update = this.needsUpdate;
        this.needsUpdate = [];
        return update;
    };
    return Usage;
}());
exports.Usage = Usage;
//# sourceMappingURL=Usage.js.map

/***/ }),

/***/ 42065:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.TeXFont = void 0;
var FontData_js_1 = __webpack_require__(32778);
var tex_js_1 = __webpack_require__(49308);
var bold_italic_js_1 = __webpack_require__(13678);
var bold_js_1 = __webpack_require__(59350);
var double_struck_js_1 = __webpack_require__(23298);
var fraktur_bold_js_1 = __webpack_require__(82399);
var fraktur_js_1 = __webpack_require__(40372);
var italic_js_1 = __webpack_require__(8649);
var largeop_js_1 = __webpack_require__(9485);
var monospace_js_1 = __webpack_require__(50102);
var normal_js_1 = __webpack_require__(43516);
var sans_serif_bold_italic_js_1 = __webpack_require__(92644);
var sans_serif_bold_js_1 = __webpack_require__(7727);
var sans_serif_italic_js_1 = __webpack_require__(62099);
var sans_serif_js_1 = __webpack_require__(53185);
var script_bold_js_1 = __webpack_require__(95662);
var script_js_1 = __webpack_require__(72569);
var smallop_js_1 = __webpack_require__(57366);
var tex_calligraphic_bold_js_1 = __webpack_require__(99983);
var tex_calligraphic_js_1 = __webpack_require__(90769);
var tex_mathit_js_1 = __webpack_require__(71370);
var tex_oldstyle_bold_js_1 = __webpack_require__(17425);
var tex_oldstyle_js_1 = __webpack_require__(41771);
var tex_size3_js_1 = __webpack_require__(21055);
var tex_size4_js_1 = __webpack_require__(39534);
var tex_variant_js_1 = __webpack_require__(27424);
var delimiters_js_1 = __webpack_require__(46922);
var TeXFont = (function (_super) {
    __extends(TeXFont, _super);
    function TeXFont() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    TeXFont.defaultCssFamilyPrefix = 'MJXZERO';
    TeXFont.defaultVariantClasses = {
        'normal': 'mjx-n',
        'bold': 'mjx-b',
        'italic': 'mjx-i',
        'bold-italic': 'mjx-b mjx-i',
        'double-struck': 'mjx-ds mjx-b',
        'fraktur': 'mjx-fr',
        'bold-fraktur': 'mjx-fr mjx-b',
        'script': 'mjx-sc mjx-i',
        'bold-script': 'mjx-sc mjx-b mjx-i',
        'sans-serif': 'mjx-ss',
        'bold-sans-serif': 'mjx-ss mjx-b',
        'sans-serif-italic': 'mjx-ss mjx-i',
        'sans-serif-bold-italic': 'mjx-ss mjx-b mjx-i',
        'monospace': 'mjx-ty',
        '-smallop': 'mjx-sop',
        '-largeop': 'mjx-lop',
        '-size3': 'mjx-s3',
        '-size4': 'mjx-s4',
        '-tex-calligraphic': 'mjx-cal mjx-i',
        '-tex-bold-calligraphic': 'mjx-cal mjx-b',
        '-tex-mathit': 'mjx-mit mjx-i',
        '-tex-oldstyle': 'mjx-os',
        '-tex-bold-oldstyle': 'mjx-os mjx-b',
        '-tex-variant': 'mjx-var'
    };
    TeXFont.defaultVariantLetters = {
        'normal': '',
        'bold': 'B',
        'italic': 'MI',
        'bold-italic': 'BI',
        'double-struck': 'A',
        'fraktur': 'FR',
        'bold-fraktur': 'FRB',
        'script': 'SC',
        'bold-script': 'SCB',
        'sans-serif': 'SS',
        'bold-sans-serif': 'SSB',
        'sans-serif-italic': 'SSI',
        'sans-serif-bold-italic': 'SSBI',
        'monospace': 'T',
        '-smallop': 'S1',
        '-largeop': 'S2',
        '-size3': 'S3',
        '-size4': 'S4',
        '-tex-calligraphic': 'C',
        '-tex-bold-calligraphic': 'CB',
        '-tex-mathit': 'MI',
        '-tex-oldstyle': 'C',
        '-tex-bold-oldstyle': 'CB',
        '-tex-variant': 'A'
    };
    TeXFont.defaultDelimiters = delimiters_js_1.delimiters;
    TeXFont.defaultChars = {
        'normal': normal_js_1.normal,
        'bold': bold_js_1.bold,
        'italic': italic_js_1.italic,
        'bold-italic': bold_italic_js_1.boldItalic,
        'double-struck': double_struck_js_1.doubleStruck,
        'fraktur': fraktur_js_1.fraktur,
        'bold-fraktur': fraktur_bold_js_1.frakturBold,
        'script': script_js_1.script,
        'bold-script': script_bold_js_1.scriptBold,
        'sans-serif': sans_serif_js_1.sansSerif,
        'bold-sans-serif': sans_serif_bold_js_1.sansSerifBold,
        'sans-serif-italic': sans_serif_italic_js_1.sansSerifItalic,
        'sans-serif-bold-italic': sans_serif_bold_italic_js_1.sansSerifBoldItalic,
        'monospace': monospace_js_1.monospace,
        '-smallop': smallop_js_1.smallop,
        '-largeop': largeop_js_1.largeop,
        '-size3': tex_size3_js_1.texSize3,
        '-size4': tex_size4_js_1.texSize4,
        '-tex-calligraphic': tex_calligraphic_js_1.texCalligraphic,
        '-tex-bold-calligraphic': tex_calligraphic_bold_js_1.texCalligraphicBold,
        '-tex-mathit': tex_mathit_js_1.texMathit,
        '-tex-oldstyle': tex_oldstyle_js_1.texOldstyle,
        '-tex-bold-oldstyle': tex_oldstyle_bold_js_1.texOldstyleBold,
        '-tex-variant': tex_variant_js_1.texVariant
    };
    TeXFont.defaultStyles = __assign(__assign({}, FontData_js_1.CHTMLFontData.defaultStyles), { '.MJX-TEX': {
            'font-family': 'MJXZERO, MJXTEX'
        }, '.TEX-B': {
            'font-family': 'MJXZERO, MJXTEX-B'
        }, '.TEX-I': {
            'font-family': 'MJXZERO, MJXTEX-I'
        }, '.TEX-MI': {
            'font-family': 'MJXZERO, MJXTEX-MI'
        }, '.TEX-BI': {
            'font-family': 'MJXZERO, MJXTEX-BI'
        }, '.TEX-S1': {
            'font-family': 'MJXZERO, MJXTEX-S1'
        }, '.TEX-S2': {
            'font-family': 'MJXZERO, MJXTEX-S2'
        }, '.TEX-S3': {
            'font-family': 'MJXZERO, MJXTEX-S3'
        }, '.TEX-S4': {
            'font-family': 'MJXZERO, MJXTEX-S4'
        }, '.TEX-A': {
            'font-family': 'MJXZERO, MJXTEX-A'
        }, '.TEX-C': {
            'font-family': 'MJXZERO, MJXTEX-C'
        }, '.TEX-CB': {
            'font-family': 'MJXZERO, MJXTEX-CB'
        }, '.TEX-FR': {
            'font-family': 'MJXZERO, MJXTEX-FR'
        }, '.TEX-FRB': {
            'font-family': 'MJXZERO, MJXTEX-FRB'
        }, '.TEX-SS': {
            'font-family': 'MJXZERO, MJXTEX-SS'
        }, '.TEX-SSB': {
            'font-family': 'MJXZERO, MJXTEX-SSB'
        }, '.TEX-SSI': {
            'font-family': 'MJXZERO, MJXTEX-SSI'
        }, '.TEX-SC': {
            'font-family': 'MJXZERO, MJXTEX-SC'
        }, '.TEX-T': {
            'font-family': 'MJXZERO, MJXTEX-T'
        }, '.TEX-V': {
            'font-family': 'MJXZERO, MJXTEX-V'
        }, '.TEX-VB': {
            'font-family': 'MJXZERO, MJXTEX-VB'
        }, 'mjx-stretchy-v mjx-c, mjx-stretchy-h mjx-c': {
            'font-family': 'MJXZERO, MJXTEX-S1, MJXTEX-S4, MJXTEX, MJXTEX-A ! important'
        } });
    TeXFont.defaultFonts = __assign(__assign({}, FontData_js_1.CHTMLFontData.defaultFonts), { '@font-face /* 1 */': {
            'font-family': 'MJXTEX',
            src: 'url("%%URL%%/MathJax_Main-Regular.woff") format("woff")'
        }, '@font-face /* 2 */': {
            'font-family': 'MJXTEX-B',
            src: 'url("%%URL%%/MathJax_Main-Bold.woff") format("woff")'
        }, '@font-face /* 3 */': {
            'font-family': 'MJXTEX-I',
            src: 'url("%%URL%%/MathJax_Math-Italic.woff") format("woff")'
        }, '@font-face /* 4 */': {
            'font-family': 'MJXTEX-MI',
            src: 'url("%%URL%%/MathJax_Main-Italic.woff") format("woff")'
        }, '@font-face /* 5 */': {
            'font-family': 'MJXTEX-BI',
            src: 'url("%%URL%%/MathJax_Math-BoldItalic.woff") format("woff")'
        }, '@font-face /* 6 */': {
            'font-family': 'MJXTEX-S1',
            src: 'url("%%URL%%/MathJax_Size1-Regular.woff") format("woff")'
        }, '@font-face /* 7 */': {
            'font-family': 'MJXTEX-S2',
            src: 'url("%%URL%%/MathJax_Size2-Regular.woff") format("woff")'
        }, '@font-face /* 8 */': {
            'font-family': 'MJXTEX-S3',
            src: 'url("%%URL%%/MathJax_Size3-Regular.woff") format("woff")'
        }, '@font-face /* 9 */': {
            'font-family': 'MJXTEX-S4',
            src: 'url("%%URL%%/MathJax_Size4-Regular.woff") format("woff")'
        }, '@font-face /* 10 */': {
            'font-family': 'MJXTEX-A',
            src: 'url("%%URL%%/MathJax_AMS-Regular.woff") format("woff")'
        }, '@font-face /* 11 */': {
            'font-family': 'MJXTEX-C',
            src: 'url("%%URL%%/MathJax_Calligraphic-Regular.woff") format("woff")'
        }, '@font-face /* 12 */': {
            'font-family': 'MJXTEX-CB',
            src: 'url("%%URL%%/MathJax_Calligraphic-Bold.woff") format("woff")'
        }, '@font-face /* 13 */': {
            'font-family': 'MJXTEX-FR',
            src: 'url("%%URL%%/MathJax_Fraktur-Regular.woff") format("woff")'
        }, '@font-face /* 14 */': {
            'font-family': 'MJXTEX-FRB',
            src: 'url("%%URL%%/MathJax_Fraktur-Bold.woff") format("woff")'
        }, '@font-face /* 15 */': {
            'font-family': 'MJXTEX-SS',
            src: 'url("%%URL%%/MathJax_SansSerif-Regular.woff") format("woff")'
        }, '@font-face /* 16 */': {
            'font-family': 'MJXTEX-SSB',
            src: 'url("%%URL%%/MathJax_SansSerif-Bold.woff") format("woff")'
        }, '@font-face /* 17 */': {
            'font-family': 'MJXTEX-SSI',
            src: 'url("%%URL%%/MathJax_SansSerif-Italic.woff") format("woff")'
        }, '@font-face /* 18 */': {
            'font-family': 'MJXTEX-SC',
            src: 'url("%%URL%%/MathJax_Script-Regular.woff") format("woff")'
        }, '@font-face /* 19 */': {
            'font-family': 'MJXTEX-T',
            src: 'url("%%URL%%/MathJax_Typewriter-Regular.woff") format("woff")'
        }, '@font-face /* 20 */': {
            'font-family': 'MJXTEX-V',
            src: 'url("%%URL%%/MathJax_Vector-Regular.woff") format("woff")'
        }, '@font-face /* 21 */': {
            'font-family': 'MJXTEX-VB',
            src: 'url("%%URL%%/MathJax_Vector-Bold.woff") format("woff")'
        } });
    return TeXFont;
}((0, tex_js_1.CommonTeXFontMixin)(FontData_js_1.CHTMLFontData)));
exports.TeXFont = TeXFont;
//# sourceMappingURL=tex.js.map

/***/ }),

/***/ 13678:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.boldItalic = void 0;
var FontData_js_1 = __webpack_require__(32778);
var bold_italic_js_1 = __webpack_require__(95884);
exports.boldItalic = (0, FontData_js_1.AddCSS)(bold_italic_js_1.boldItalic, {
    0x131: { f: 'B' },
    0x237: { f: 'B' },
    0x2044: { c: '/' },
    0x2206: { c: '\\394' },
    0x29F8: { c: '/' },
});
//# sourceMappingURL=bold-italic.js.map

/***/ }),

/***/ 59350:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.bold = void 0;
var FontData_js_1 = __webpack_require__(32778);
var bold_js_1 = __webpack_require__(75083);
exports.bold = (0, FontData_js_1.AddCSS)(bold_js_1.bold, {
    0xB7: { c: '\\22C5' },
    0x131: { f: '' },
    0x237: { f: '' },
    0x2B9: { c: '\\2032' },
    0x2002: { c: '' },
    0x2003: { c: '' },
    0x2004: { c: '' },
    0x2005: { c: '' },
    0x2006: { c: '' },
    0x2009: { c: '' },
    0x200A: { c: '' },
    0x2015: { c: '\\2014' },
    0x2016: { c: '\\2225' },
    0x2017: { c: '_' },
    0x2022: { c: '\\2219' },
    0x2033: { c: '\\2032\\2032' },
    0x2034: { c: '\\2032\\2032\\2032' },
    0x203E: { c: '\\2C9' },
    0x2044: { c: '/' },
    0x2057: { c: '\\2032\\2032\\2032\\2032' },
    0x20D7: { c: '\\2192', f: 'VB' },
    0x219A: { c: '\\2190\\338' },
    0x219B: { c: '\\2192\\338' },
    0x21AE: { c: '\\2194\\338' },
    0x21CD: { c: '\\21D0\\338' },
    0x21CE: { c: '\\21D4\\338' },
    0x21CF: { c: '\\21D2\\338' },
    0x2204: { c: '\\2203\\338' },
    0x2206: { c: '\\394' },
    0x220C: { c: '\\220B\\338' },
    0x2224: { c: '\\2223\\338' },
    0x2226: { c: '\\2225\\338' },
    0x2241: { c: '\\223C\\338' },
    0x2244: { c: '\\2243\\338' },
    0x2247: { c: '\\2245\\338' },
    0x2249: { c: '\\2248\\338' },
    0x2262: { c: '\\2261\\338' },
    0x226D: { c: '\\224D\\338' },
    0x226E: { c: '<\\338' },
    0x226F: { c: '>\\338' },
    0x2270: { c: '\\2264\\338' },
    0x2271: { c: '\\2265\\338' },
    0x2280: { c: '\\227A\\338' },
    0x2281: { c: '\\227B\\338' },
    0x2284: { c: '\\2282\\338' },
    0x2285: { c: '\\2283\\338' },
    0x2288: { c: '\\2286\\338' },
    0x2289: { c: '\\2287\\338' },
    0x22AC: { c: '\\22A2\\338' },
    0x22AD: { c: '\\22A8\\338' },
    0x22E2: { c: '\\2291\\338' },
    0x22E3: { c: '\\2292\\338' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x25B5: { c: '\\25B3' },
    0x25BF: { c: '\\25BD' },
    0x2758: { c: '\\2223' },
    0x29F8: { c: '/', f: 'BI' },
    0x2A2F: { c: '\\D7' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
});
//# sourceMappingURL=bold.js.map

/***/ }),

/***/ 23298:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.doubleStruck = void 0;
var double_struck_js_1 = __webpack_require__(94841);
Object.defineProperty(exports, "doubleStruck", ({ enumerable: true, get: function () { return double_struck_js_1.doubleStruck; } }));
//# sourceMappingURL=double-struck.js.map

/***/ }),

/***/ 82399:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.frakturBold = void 0;
var FontData_js_1 = __webpack_require__(32778);
var fraktur_bold_js_1 = __webpack_require__(19173);
exports.frakturBold = (0, FontData_js_1.AddCSS)(fraktur_bold_js_1.frakturBold, {
    0x2044: { c: '/' },
});
//# sourceMappingURL=fraktur-bold.js.map

/***/ }),

/***/ 40372:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.fraktur = void 0;
var FontData_js_1 = __webpack_require__(32778);
var fraktur_js_1 = __webpack_require__(57246);
exports.fraktur = (0, FontData_js_1.AddCSS)(fraktur_js_1.fraktur, {
    0x2044: { c: '/' },
});
//# sourceMappingURL=fraktur.js.map

/***/ }),

/***/ 8649:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.italic = void 0;
var FontData_js_1 = __webpack_require__(32778);
var italic_js_1 = __webpack_require__(74629);
exports.italic = (0, FontData_js_1.AddCSS)(italic_js_1.italic, {
    0x2F: { f: 'I' },
    0x3DD: { c: '\\E008', f: 'A' },
    0x2015: { c: '\\2014' },
    0x2017: { c: '_' },
    0x2044: { c: '/', f: 'I' },
    0x2206: { c: '\\394', f: 'I' },
    0x29F8: { c: '/', f: 'I' },
});
//# sourceMappingURL=italic.js.map

/***/ }),

/***/ 9485:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.largeop = void 0;
var FontData_js_1 = __webpack_require__(32778);
var largeop_js_1 = __webpack_require__(8709);
exports.largeop = (0, FontData_js_1.AddCSS)(largeop_js_1.largeop, {
    0x2016: { f: 'S1' },
    0x2044: { c: '/' },
    0x2191: { f: 'S1' },
    0x2193: { f: 'S1' },
    0x21D1: { f: 'S1' },
    0x21D3: { f: 'S1' },
    0x2223: { f: 'S1' },
    0x2225: { f: 'S1' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x23D0: { f: 'S1' },
    0x2758: { c: '\\2223', f: 'S1' },
    0x2A0C: { c: '\\222C\\222C' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
});
//# sourceMappingURL=largeop.js.map

/***/ }),

/***/ 50102:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.monospace = void 0;
var FontData_js_1 = __webpack_require__(32778);
var monospace_js_1 = __webpack_require__(60984);
exports.monospace = (0, FontData_js_1.AddCSS)(monospace_js_1.monospace, {
    0x2B9: { c: '\\2032' },
    0x391: { c: 'A' },
    0x392: { c: 'B' },
    0x395: { c: 'E' },
    0x396: { c: 'Z' },
    0x397: { c: 'H' },
    0x399: { c: 'I' },
    0x39A: { c: 'K' },
    0x39C: { c: 'M' },
    0x39D: { c: 'N' },
    0x39F: { c: 'O' },
    0x3A1: { c: 'P' },
    0x3A4: { c: 'T' },
    0x3A7: { c: 'X' },
    0x2017: { c: '_' },
    0x2033: { c: '\\2032\\2032' },
    0x2034: { c: '\\2032\\2032\\2032' },
    0x2044: { c: '/' },
    0x2057: { c: '\\2032\\2032\\2032\\2032' },
    0x2206: { c: '\\394' },
});
//# sourceMappingURL=monospace.js.map

/***/ }),

/***/ 43516:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.normal = void 0;
var FontData_js_1 = __webpack_require__(32778);
var normal_js_1 = __webpack_require__(65037);
exports.normal = (0, FontData_js_1.AddCSS)(normal_js_1.normal, {
    0xA3: { f: 'MI' },
    0xA5: { f: 'A' },
    0xAE: { f: 'A' },
    0xB7: { c: '\\22C5' },
    0xF0: { f: 'A' },
    0x2B9: { c: '\\2032' },
    0x391: { c: 'A' },
    0x392: { c: 'B' },
    0x395: { c: 'E' },
    0x396: { c: 'Z' },
    0x397: { c: 'H' },
    0x399: { c: 'I' },
    0x39A: { c: 'K' },
    0x39C: { c: 'M' },
    0x39D: { c: 'N' },
    0x39F: { c: 'O' },
    0x3A1: { c: 'P' },
    0x3A4: { c: 'T' },
    0x3A7: { c: 'X' },
    0x2000: { c: '' },
    0x2001: { c: '' },
    0x2002: { c: '' },
    0x2003: { c: '' },
    0x2004: { c: '' },
    0x2005: { c: '' },
    0x2006: { c: '' },
    0x2009: { c: '' },
    0x200A: { c: '' },
    0x200B: { c: '' },
    0x200C: { c: '' },
    0x2015: { c: '\\2014' },
    0x2016: { c: '\\2225' },
    0x2017: { c: '_' },
    0x2022: { c: '\\2219' },
    0x2033: { c: '\\2032\\2032' },
    0x2034: { c: '\\2032\\2032\\2032' },
    0x2035: { f: 'A' },
    0x2036: { c: '\\2035\\2035', f: 'A' },
    0x2037: { c: '\\2035\\2035\\2035', f: 'A' },
    0x203E: { c: '\\2C9' },
    0x2044: { c: '/' },
    0x2057: { c: '\\2032\\2032\\2032\\2032' },
    0x2060: { c: '' },
    0x2061: { c: '' },
    0x2062: { c: '' },
    0x2063: { c: '' },
    0x2064: { c: '' },
    0x20D7: { c: '\\2192', f: 'V' },
    0x2102: { c: 'C', f: 'A' },
    0x210B: { c: 'H', f: 'SC' },
    0x210C: { c: 'H', f: 'FR' },
    0x210D: { c: 'H', f: 'A' },
    0x210E: { c: 'h', f: 'I' },
    0x210F: { f: 'A' },
    0x2110: { c: 'I', f: 'SC' },
    0x2111: { c: 'I', f: 'FR' },
    0x2112: { c: 'L', f: 'SC' },
    0x2115: { c: 'N', f: 'A' },
    0x2119: { c: 'P', f: 'A' },
    0x211A: { c: 'Q', f: 'A' },
    0x211B: { c: 'R', f: 'SC' },
    0x211C: { c: 'R', f: 'FR' },
    0x211D: { c: 'R', f: 'A' },
    0x2124: { c: 'Z', f: 'A' },
    0x2126: { c: '\\3A9' },
    0x2127: { f: 'A' },
    0x2128: { c: 'Z', f: 'FR' },
    0x212C: { c: 'B', f: 'SC' },
    0x212D: { c: 'C', f: 'FR' },
    0x2130: { c: 'E', f: 'SC' },
    0x2131: { c: 'F', f: 'SC' },
    0x2132: { f: 'A' },
    0x2133: { c: 'M', f: 'SC' },
    0x2136: { f: 'A' },
    0x2137: { f: 'A' },
    0x2138: { f: 'A' },
    0x2141: { f: 'A' },
    0x219A: { f: 'A' },
    0x219B: { f: 'A' },
    0x219E: { f: 'A' },
    0x21A0: { f: 'A' },
    0x21A2: { f: 'A' },
    0x21A3: { f: 'A' },
    0x21AB: { f: 'A' },
    0x21AC: { f: 'A' },
    0x21AD: { f: 'A' },
    0x21AE: { f: 'A' },
    0x21B0: { f: 'A' },
    0x21B1: { f: 'A' },
    0x21B6: { f: 'A' },
    0x21B7: { f: 'A' },
    0x21BA: { f: 'A' },
    0x21BB: { f: 'A' },
    0x21BE: { f: 'A' },
    0x21BF: { f: 'A' },
    0x21C2: { f: 'A' },
    0x21C3: { f: 'A' },
    0x21C4: { f: 'A' },
    0x21C6: { f: 'A' },
    0x21C7: { f: 'A' },
    0x21C8: { f: 'A' },
    0x21C9: { f: 'A' },
    0x21CA: { f: 'A' },
    0x21CB: { f: 'A' },
    0x21CD: { f: 'A' },
    0x21CE: { f: 'A' },
    0x21CF: { f: 'A' },
    0x21DA: { f: 'A' },
    0x21DB: { f: 'A' },
    0x21DD: { f: 'A' },
    0x21E0: { f: 'A' },
    0x21E2: { f: 'A' },
    0x2201: { f: 'A' },
    0x2204: { c: '\\2203\\338' },
    0x2206: { c: '\\394' },
    0x220C: { c: '\\220B\\338' },
    0x220D: { f: 'A' },
    0x220F: { f: 'S1' },
    0x2210: { f: 'S1' },
    0x2211: { f: 'S1' },
    0x2214: { f: 'A' },
    0x2221: { f: 'A' },
    0x2222: { f: 'A' },
    0x2224: { f: 'A' },
    0x2226: { f: 'A' },
    0x222C: { f: 'S1' },
    0x222D: { f: 'S1' },
    0x222E: { f: 'S1' },
    0x2234: { f: 'A' },
    0x2235: { f: 'A' },
    0x223D: { f: 'A' },
    0x2241: { f: 'A' },
    0x2242: { f: 'A' },
    0x2244: { c: '\\2243\\338' },
    0x2247: { c: '\\2246', f: 'A' },
    0x2249: { c: '\\2248\\338' },
    0x224A: { f: 'A' },
    0x224E: { f: 'A' },
    0x224F: { f: 'A' },
    0x2251: { f: 'A' },
    0x2252: { f: 'A' },
    0x2253: { f: 'A' },
    0x2256: { f: 'A' },
    0x2257: { f: 'A' },
    0x225C: { f: 'A' },
    0x2262: { c: '\\2261\\338' },
    0x2266: { f: 'A' },
    0x2267: { f: 'A' },
    0x2268: { f: 'A' },
    0x2269: { f: 'A' },
    0x226C: { f: 'A' },
    0x226D: { c: '\\224D\\338' },
    0x226E: { f: 'A' },
    0x226F: { f: 'A' },
    0x2270: { f: 'A' },
    0x2271: { f: 'A' },
    0x2272: { f: 'A' },
    0x2273: { f: 'A' },
    0x2274: { c: '\\2272\\338' },
    0x2275: { c: '\\2273\\338' },
    0x2276: { f: 'A' },
    0x2277: { f: 'A' },
    0x2278: { c: '\\2276\\338' },
    0x2279: { c: '\\2277\\338' },
    0x227C: { f: 'A' },
    0x227D: { f: 'A' },
    0x227E: { f: 'A' },
    0x227F: { f: 'A' },
    0x2280: { f: 'A' },
    0x2281: { f: 'A' },
    0x2284: { c: '\\2282\\338' },
    0x2285: { c: '\\2283\\338' },
    0x2288: { f: 'A' },
    0x2289: { f: 'A' },
    0x228A: { f: 'A' },
    0x228B: { f: 'A' },
    0x228F: { f: 'A' },
    0x2290: { f: 'A' },
    0x229A: { f: 'A' },
    0x229B: { f: 'A' },
    0x229D: { f: 'A' },
    0x229E: { f: 'A' },
    0x229F: { f: 'A' },
    0x22A0: { f: 'A' },
    0x22A1: { f: 'A' },
    0x22A9: { f: 'A' },
    0x22AA: { f: 'A' },
    0x22AC: { f: 'A' },
    0x22AD: { f: 'A' },
    0x22AE: { f: 'A' },
    0x22AF: { f: 'A' },
    0x22B2: { f: 'A' },
    0x22B3: { f: 'A' },
    0x22B4: { f: 'A' },
    0x22B5: { f: 'A' },
    0x22B8: { f: 'A' },
    0x22BA: { f: 'A' },
    0x22BB: { f: 'A' },
    0x22BC: { f: 'A' },
    0x22C0: { f: 'S1' },
    0x22C1: { f: 'S1' },
    0x22C2: { f: 'S1' },
    0x22C3: { f: 'S1' },
    0x22C7: { f: 'A' },
    0x22C9: { f: 'A' },
    0x22CA: { f: 'A' },
    0x22CB: { f: 'A' },
    0x22CC: { f: 'A' },
    0x22CD: { f: 'A' },
    0x22CE: { f: 'A' },
    0x22CF: { f: 'A' },
    0x22D0: { f: 'A' },
    0x22D1: { f: 'A' },
    0x22D2: { f: 'A' },
    0x22D3: { f: 'A' },
    0x22D4: { f: 'A' },
    0x22D6: { f: 'A' },
    0x22D7: { f: 'A' },
    0x22D8: { f: 'A' },
    0x22D9: { f: 'A' },
    0x22DA: { f: 'A' },
    0x22DB: { f: 'A' },
    0x22DE: { f: 'A' },
    0x22DF: { f: 'A' },
    0x22E0: { f: 'A' },
    0x22E1: { f: 'A' },
    0x22E2: { c: '\\2291\\338' },
    0x22E3: { c: '\\2292\\338' },
    0x22E6: { f: 'A' },
    0x22E7: { f: 'A' },
    0x22E8: { f: 'A' },
    0x22E9: { f: 'A' },
    0x22EA: { f: 'A' },
    0x22EB: { f: 'A' },
    0x22EC: { f: 'A' },
    0x22ED: { f: 'A' },
    0x2305: { c: '\\22BC', f: 'A' },
    0x2306: { c: '\\2A5E', f: 'A' },
    0x231C: { c: '\\250C', f: 'A' },
    0x231D: { c: '\\2510', f: 'A' },
    0x231E: { c: '\\2514', f: 'A' },
    0x231F: { c: '\\2518', f: 'A' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x23D0: { f: 'S1' },
    0x24C8: { f: 'A' },
    0x250C: { f: 'A' },
    0x2510: { f: 'A' },
    0x2514: { f: 'A' },
    0x2518: { f: 'A' },
    0x2571: { f: 'A' },
    0x2572: { f: 'A' },
    0x25A0: { f: 'A' },
    0x25A1: { f: 'A' },
    0x25AA: { c: '\\25A0', f: 'A' },
    0x25B2: { f: 'A' },
    0x25B4: { c: '\\25B2', f: 'A' },
    0x25B5: { c: '\\25B3' },
    0x25B6: { f: 'A' },
    0x25B8: { c: '\\25B6', f: 'A' },
    0x25BC: { f: 'A' },
    0x25BE: { c: '\\25BC', f: 'A' },
    0x25BF: { c: '\\25BD' },
    0x25C0: { f: 'A' },
    0x25C2: { c: '\\25C0', f: 'A' },
    0x25CA: { f: 'A' },
    0x25FB: { c: '\\25A1', f: 'A' },
    0x25FC: { c: '\\25A0', f: 'A' },
    0x2605: { f: 'A' },
    0x2713: { f: 'A' },
    0x2720: { f: 'A' },
    0x2758: { c: '\\2223' },
    0x29EB: { f: 'A' },
    0x29F8: { c: '/', f: 'I' },
    0x2A00: { f: 'S1' },
    0x2A01: { f: 'S1' },
    0x2A02: { f: 'S1' },
    0x2A04: { f: 'S1' },
    0x2A06: { f: 'S1' },
    0x2A0C: { c: '\\222C\\222C', f: 'S1' },
    0x2A2F: { c: '\\D7' },
    0x2A5E: { f: 'A' },
    0x2A7D: { f: 'A' },
    0x2A7E: { f: 'A' },
    0x2A85: { f: 'A' },
    0x2A86: { f: 'A' },
    0x2A87: { f: 'A' },
    0x2A88: { f: 'A' },
    0x2A89: { f: 'A' },
    0x2A8A: { f: 'A' },
    0x2A8B: { f: 'A' },
    0x2A8C: { f: 'A' },
    0x2A95: { f: 'A' },
    0x2A96: { f: 'A' },
    0x2AB5: { f: 'A' },
    0x2AB6: { f: 'A' },
    0x2AB7: { f: 'A' },
    0x2AB8: { f: 'A' },
    0x2AB9: { f: 'A' },
    0x2ABA: { f: 'A' },
    0x2AC5: { f: 'A' },
    0x2AC6: { f: 'A' },
    0x2ACB: { f: 'A' },
    0x2ACC: { f: 'A' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
    0xE006: { f: 'A' },
    0xE007: { f: 'A' },
    0xE008: { f: 'A' },
    0xE009: { f: 'A' },
    0xE00C: { f: 'A' },
    0xE00D: { f: 'A' },
    0xE00E: { f: 'A' },
    0xE00F: { f: 'A' },
    0xE010: { f: 'A' },
    0xE011: { f: 'A' },
    0xE016: { f: 'A' },
    0xE017: { f: 'A' },
    0xE018: { f: 'A' },
    0xE019: { f: 'A' },
    0xE01A: { f: 'A' },
    0xE01B: { f: 'A' },
    0x1D400: { c: 'A', f: 'B' },
    0x1D401: { c: 'B', f: 'B' },
    0x1D402: { c: 'C', f: 'B' },
    0x1D403: { c: 'D', f: 'B' },
    0x1D404: { c: 'E', f: 'B' },
    0x1D405: { c: 'F', f: 'B' },
    0x1D406: { c: 'G', f: 'B' },
    0x1D407: { c: 'H', f: 'B' },
    0x1D408: { c: 'I', f: 'B' },
    0x1D409: { c: 'J', f: 'B' },
    0x1D40A: { c: 'K', f: 'B' },
    0x1D40B: { c: 'L', f: 'B' },
    0x1D40C: { c: 'M', f: 'B' },
    0x1D40D: { c: 'N', f: 'B' },
    0x1D40E: { c: 'O', f: 'B' },
    0x1D40F: { c: 'P', f: 'B' },
    0x1D410: { c: 'Q', f: 'B' },
    0x1D411: { c: 'R', f: 'B' },
    0x1D412: { c: 'S', f: 'B' },
    0x1D413: { c: 'T', f: 'B' },
    0x1D414: { c: 'U', f: 'B' },
    0x1D415: { c: 'V', f: 'B' },
    0x1D416: { c: 'W', f: 'B' },
    0x1D417: { c: 'X', f: 'B' },
    0x1D418: { c: 'Y', f: 'B' },
    0x1D419: { c: 'Z', f: 'B' },
    0x1D41A: { c: 'a', f: 'B' },
    0x1D41B: { c: 'b', f: 'B' },
    0x1D41C: { c: 'c', f: 'B' },
    0x1D41D: { c: 'd', f: 'B' },
    0x1D41E: { c: 'e', f: 'B' },
    0x1D41F: { c: 'f', f: 'B' },
    0x1D420: { c: 'g', f: 'B' },
    0x1D421: { c: 'h', f: 'B' },
    0x1D422: { c: 'i', f: 'B' },
    0x1D423: { c: 'j', f: 'B' },
    0x1D424: { c: 'k', f: 'B' },
    0x1D425: { c: 'l', f: 'B' },
    0x1D426: { c: 'm', f: 'B' },
    0x1D427: { c: 'n', f: 'B' },
    0x1D428: { c: 'o', f: 'B' },
    0x1D429: { c: 'p', f: 'B' },
    0x1D42A: { c: 'q', f: 'B' },
    0x1D42B: { c: 'r', f: 'B' },
    0x1D42C: { c: 's', f: 'B' },
    0x1D42D: { c: 't', f: 'B' },
    0x1D42E: { c: 'u', f: 'B' },
    0x1D42F: { c: 'v', f: 'B' },
    0x1D430: { c: 'w', f: 'B' },
    0x1D431: { c: 'x', f: 'B' },
    0x1D432: { c: 'y', f: 'B' },
    0x1D433: { c: 'z', f: 'B' },
    0x1D434: { c: 'A', f: 'I' },
    0x1D435: { c: 'B', f: 'I' },
    0x1D436: { c: 'C', f: 'I' },
    0x1D437: { c: 'D', f: 'I' },
    0x1D438: { c: 'E', f: 'I' },
    0x1D439: { c: 'F', f: 'I' },
    0x1D43A: { c: 'G', f: 'I' },
    0x1D43B: { c: 'H', f: 'I' },
    0x1D43C: { c: 'I', f: 'I' },
    0x1D43D: { c: 'J', f: 'I' },
    0x1D43E: { c: 'K', f: 'I' },
    0x1D43F: { c: 'L', f: 'I' },
    0x1D440: { c: 'M', f: 'I' },
    0x1D441: { c: 'N', f: 'I' },
    0x1D442: { c: 'O', f: 'I' },
    0x1D443: { c: 'P', f: 'I' },
    0x1D444: { c: 'Q', f: 'I' },
    0x1D445: { c: 'R', f: 'I' },
    0x1D446: { c: 'S', f: 'I' },
    0x1D447: { c: 'T', f: 'I' },
    0x1D448: { c: 'U', f: 'I' },
    0x1D449: { c: 'V', f: 'I' },
    0x1D44A: { c: 'W', f: 'I' },
    0x1D44B: { c: 'X', f: 'I' },
    0x1D44C: { c: 'Y', f: 'I' },
    0x1D44D: { c: 'Z', f: 'I' },
    0x1D44E: { c: 'a', f: 'I' },
    0x1D44F: { c: 'b', f: 'I' },
    0x1D450: { c: 'c', f: 'I' },
    0x1D451: { c: 'd', f: 'I' },
    0x1D452: { c: 'e', f: 'I' },
    0x1D453: { c: 'f', f: 'I' },
    0x1D454: { c: 'g', f: 'I' },
    0x1D456: { c: 'i', f: 'I' },
    0x1D457: { c: 'j', f: 'I' },
    0x1D458: { c: 'k', f: 'I' },
    0x1D459: { c: 'l', f: 'I' },
    0x1D45A: { c: 'm', f: 'I' },
    0x1D45B: { c: 'n', f: 'I' },
    0x1D45C: { c: 'o', f: 'I' },
    0x1D45D: { c: 'p', f: 'I' },
    0x1D45E: { c: 'q', f: 'I' },
    0x1D45F: { c: 'r', f: 'I' },
    0x1D460: { c: 's', f: 'I' },
    0x1D461: { c: 't', f: 'I' },
    0x1D462: { c: 'u', f: 'I' },
    0x1D463: { c: 'v', f: 'I' },
    0x1D464: { c: 'w', f: 'I' },
    0x1D465: { c: 'x', f: 'I' },
    0x1D466: { c: 'y', f: 'I' },
    0x1D467: { c: 'z', f: 'I' },
    0x1D468: { c: 'A', f: 'BI' },
    0x1D469: { c: 'B', f: 'BI' },
    0x1D46A: { c: 'C', f: 'BI' },
    0x1D46B: { c: 'D', f: 'BI' },
    0x1D46C: { c: 'E', f: 'BI' },
    0x1D46D: { c: 'F', f: 'BI' },
    0x1D46E: { c: 'G', f: 'BI' },
    0x1D46F: { c: 'H', f: 'BI' },
    0x1D470: { c: 'I', f: 'BI' },
    0x1D471: { c: 'J', f: 'BI' },
    0x1D472: { c: 'K', f: 'BI' },
    0x1D473: { c: 'L', f: 'BI' },
    0x1D474: { c: 'M', f: 'BI' },
    0x1D475: { c: 'N', f: 'BI' },
    0x1D476: { c: 'O', f: 'BI' },
    0x1D477: { c: 'P', f: 'BI' },
    0x1D478: { c: 'Q', f: 'BI' },
    0x1D479: { c: 'R', f: 'BI' },
    0x1D47A: { c: 'S', f: 'BI' },
    0x1D47B: { c: 'T', f: 'BI' },
    0x1D47C: { c: 'U', f: 'BI' },
    0x1D47D: { c: 'V', f: 'BI' },
    0x1D47E: { c: 'W', f: 'BI' },
    0x1D47F: { c: 'X', f: 'BI' },
    0x1D480: { c: 'Y', f: 'BI' },
    0x1D481: { c: 'Z', f: 'BI' },
    0x1D482: { c: 'a', f: 'BI' },
    0x1D483: { c: 'b', f: 'BI' },
    0x1D484: { c: 'c', f: 'BI' },
    0x1D485: { c: 'd', f: 'BI' },
    0x1D486: { c: 'e', f: 'BI' },
    0x1D487: { c: 'f', f: 'BI' },
    0x1D488: { c: 'g', f: 'BI' },
    0x1D489: { c: 'h', f: 'BI' },
    0x1D48A: { c: 'i', f: 'BI' },
    0x1D48B: { c: 'j', f: 'BI' },
    0x1D48C: { c: 'k', f: 'BI' },
    0x1D48D: { c: 'l', f: 'BI' },
    0x1D48E: { c: 'm', f: 'BI' },
    0x1D48F: { c: 'n', f: 'BI' },
    0x1D490: { c: 'o', f: 'BI' },
    0x1D491: { c: 'p', f: 'BI' },
    0x1D492: { c: 'q', f: 'BI' },
    0x1D493: { c: 'r', f: 'BI' },
    0x1D494: { c: 's', f: 'BI' },
    0x1D495: { c: 't', f: 'BI' },
    0x1D496: { c: 'u', f: 'BI' },
    0x1D497: { c: 'v', f: 'BI' },
    0x1D498: { c: 'w', f: 'BI' },
    0x1D499: { c: 'x', f: 'BI' },
    0x1D49A: { c: 'y', f: 'BI' },
    0x1D49B: { c: 'z', f: 'BI' },
    0x1D49C: { c: 'A', f: 'SC' },
    0x1D49E: { c: 'C', f: 'SC' },
    0x1D49F: { c: 'D', f: 'SC' },
    0x1D4A2: { c: 'G', f: 'SC' },
    0x1D4A5: { c: 'J', f: 'SC' },
    0x1D4A6: { c: 'K', f: 'SC' },
    0x1D4A9: { c: 'N', f: 'SC' },
    0x1D4AA: { c: 'O', f: 'SC' },
    0x1D4AB: { c: 'P', f: 'SC' },
    0x1D4AC: { c: 'Q', f: 'SC' },
    0x1D4AE: { c: 'S', f: 'SC' },
    0x1D4AF: { c: 'T', f: 'SC' },
    0x1D4B0: { c: 'U', f: 'SC' },
    0x1D4B1: { c: 'V', f: 'SC' },
    0x1D4B2: { c: 'W', f: 'SC' },
    0x1D4B3: { c: 'X', f: 'SC' },
    0x1D4B4: { c: 'Y', f: 'SC' },
    0x1D4B5: { c: 'Z', f: 'SC' },
    0x1D504: { c: 'A', f: 'FR' },
    0x1D505: { c: 'B', f: 'FR' },
    0x1D507: { c: 'D', f: 'FR' },
    0x1D508: { c: 'E', f: 'FR' },
    0x1D509: { c: 'F', f: 'FR' },
    0x1D50A: { c: 'G', f: 'FR' },
    0x1D50D: { c: 'J', f: 'FR' },
    0x1D50E: { c: 'K', f: 'FR' },
    0x1D50F: { c: 'L', f: 'FR' },
    0x1D510: { c: 'M', f: 'FR' },
    0x1D511: { c: 'N', f: 'FR' },
    0x1D512: { c: 'O', f: 'FR' },
    0x1D513: { c: 'P', f: 'FR' },
    0x1D514: { c: 'Q', f: 'FR' },
    0x1D516: { c: 'S', f: 'FR' },
    0x1D517: { c: 'T', f: 'FR' },
    0x1D518: { c: 'U', f: 'FR' },
    0x1D519: { c: 'V', f: 'FR' },
    0x1D51A: { c: 'W', f: 'FR' },
    0x1D51B: { c: 'X', f: 'FR' },
    0x1D51C: { c: 'Y', f: 'FR' },
    0x1D51E: { c: 'a', f: 'FR' },
    0x1D51F: { c: 'b', f: 'FR' },
    0x1D520: { c: 'c', f: 'FR' },
    0x1D521: { c: 'd', f: 'FR' },
    0x1D522: { c: 'e', f: 'FR' },
    0x1D523: { c: 'f', f: 'FR' },
    0x1D524: { c: 'g', f: 'FR' },
    0x1D525: { c: 'h', f: 'FR' },
    0x1D526: { c: 'i', f: 'FR' },
    0x1D527: { c: 'j', f: 'FR' },
    0x1D528: { c: 'k', f: 'FR' },
    0x1D529: { c: 'l', f: 'FR' },
    0x1D52A: { c: 'm', f: 'FR' },
    0x1D52B: { c: 'n', f: 'FR' },
    0x1D52C: { c: 'o', f: 'FR' },
    0x1D52D: { c: 'p', f: 'FR' },
    0x1D52E: { c: 'q', f: 'FR' },
    0x1D52F: { c: 'r', f: 'FR' },
    0x1D530: { c: 's', f: 'FR' },
    0x1D531: { c: 't', f: 'FR' },
    0x1D532: { c: 'u', f: 'FR' },
    0x1D533: { c: 'v', f: 'FR' },
    0x1D534: { c: 'w', f: 'FR' },
    0x1D535: { c: 'x', f: 'FR' },
    0x1D536: { c: 'y', f: 'FR' },
    0x1D537: { c: 'z', f: 'FR' },
    0x1D538: { c: 'A', f: 'A' },
    0x1D539: { c: 'B', f: 'A' },
    0x1D53B: { c: 'D', f: 'A' },
    0x1D53C: { c: 'E', f: 'A' },
    0x1D53D: { c: 'F', f: 'A' },
    0x1D53E: { c: 'G', f: 'A' },
    0x1D540: { c: 'I', f: 'A' },
    0x1D541: { c: 'J', f: 'A' },
    0x1D542: { c: 'K', f: 'A' },
    0x1D543: { c: 'L', f: 'A' },
    0x1D544: { c: 'M', f: 'A' },
    0x1D546: { c: 'O', f: 'A' },
    0x1D54A: { c: 'S', f: 'A' },
    0x1D54B: { c: 'T', f: 'A' },
    0x1D54C: { c: 'U', f: 'A' },
    0x1D54D: { c: 'V', f: 'A' },
    0x1D54E: { c: 'W', f: 'A' },
    0x1D54F: { c: 'X', f: 'A' },
    0x1D550: { c: 'Y', f: 'A' },
    0x1D56C: { c: 'A', f: 'FRB' },
    0x1D56D: { c: 'B', f: 'FRB' },
    0x1D56E: { c: 'C', f: 'FRB' },
    0x1D56F: { c: 'D', f: 'FRB' },
    0x1D570: { c: 'E', f: 'FRB' },
    0x1D571: { c: 'F', f: 'FRB' },
    0x1D572: { c: 'G', f: 'FRB' },
    0x1D573: { c: 'H', f: 'FRB' },
    0x1D574: { c: 'I', f: 'FRB' },
    0x1D575: { c: 'J', f: 'FRB' },
    0x1D576: { c: 'K', f: 'FRB' },
    0x1D577: { c: 'L', f: 'FRB' },
    0x1D578: { c: 'M', f: 'FRB' },
    0x1D579: { c: 'N', f: 'FRB' },
    0x1D57A: { c: 'O', f: 'FRB' },
    0x1D57B: { c: 'P', f: 'FRB' },
    0x1D57C: { c: 'Q', f: 'FRB' },
    0x1D57D: { c: 'R', f: 'FRB' },
    0x1D57E: { c: 'S', f: 'FRB' },
    0x1D57F: { c: 'T', f: 'FRB' },
    0x1D580: { c: 'U', f: 'FRB' },
    0x1D581: { c: 'V', f: 'FRB' },
    0x1D582: { c: 'W', f: 'FRB' },
    0x1D583: { c: 'X', f: 'FRB' },
    0x1D584: { c: 'Y', f: 'FRB' },
    0x1D585: { c: 'Z', f: 'FRB' },
    0x1D586: { c: 'a', f: 'FRB' },
    0x1D587: { c: 'b', f: 'FRB' },
    0x1D588: { c: 'c', f: 'FRB' },
    0x1D589: { c: 'd', f: 'FRB' },
    0x1D58A: { c: 'e', f: 'FRB' },
    0x1D58B: { c: 'f', f: 'FRB' },
    0x1D58C: { c: 'g', f: 'FRB' },
    0x1D58D: { c: 'h', f: 'FRB' },
    0x1D58E: { c: 'i', f: 'FRB' },
    0x1D58F: { c: 'j', f: 'FRB' },
    0x1D590: { c: 'k', f: 'FRB' },
    0x1D591: { c: 'l', f: 'FRB' },
    0x1D592: { c: 'm', f: 'FRB' },
    0x1D593: { c: 'n', f: 'FRB' },
    0x1D594: { c: 'o', f: 'FRB' },
    0x1D595: { c: 'p', f: 'FRB' },
    0x1D596: { c: 'q', f: 'FRB' },
    0x1D597: { c: 'r', f: 'FRB' },
    0x1D598: { c: 's', f: 'FRB' },
    0x1D599: { c: 't', f: 'FRB' },
    0x1D59A: { c: 'u', f: 'FRB' },
    0x1D59B: { c: 'v', f: 'FRB' },
    0x1D59C: { c: 'w', f: 'FRB' },
    0x1D59D: { c: 'x', f: 'FRB' },
    0x1D59E: { c: 'y', f: 'FRB' },
    0x1D59F: { c: 'z', f: 'FRB' },
    0x1D5A0: { c: 'A', f: 'SS' },
    0x1D5A1: { c: 'B', f: 'SS' },
    0x1D5A2: { c: 'C', f: 'SS' },
    0x1D5A3: { c: 'D', f: 'SS' },
    0x1D5A4: { c: 'E', f: 'SS' },
    0x1D5A5: { c: 'F', f: 'SS' },
    0x1D5A6: { c: 'G', f: 'SS' },
    0x1D5A7: { c: 'H', f: 'SS' },
    0x1D5A8: { c: 'I', f: 'SS' },
    0x1D5A9: { c: 'J', f: 'SS' },
    0x1D5AA: { c: 'K', f: 'SS' },
    0x1D5AB: { c: 'L', f: 'SS' },
    0x1D5AC: { c: 'M', f: 'SS' },
    0x1D5AD: { c: 'N', f: 'SS' },
    0x1D5AE: { c: 'O', f: 'SS' },
    0x1D5AF: { c: 'P', f: 'SS' },
    0x1D5B0: { c: 'Q', f: 'SS' },
    0x1D5B1: { c: 'R', f: 'SS' },
    0x1D5B2: { c: 'S', f: 'SS' },
    0x1D5B3: { c: 'T', f: 'SS' },
    0x1D5B4: { c: 'U', f: 'SS' },
    0x1D5B5: { c: 'V', f: 'SS' },
    0x1D5B6: { c: 'W', f: 'SS' },
    0x1D5B7: { c: 'X', f: 'SS' },
    0x1D5B8: { c: 'Y', f: 'SS' },
    0x1D5B9: { c: 'Z', f: 'SS' },
    0x1D5BA: { c: 'a', f: 'SS' },
    0x1D5BB: { c: 'b', f: 'SS' },
    0x1D5BC: { c: 'c', f: 'SS' },
    0x1D5BD: { c: 'd', f: 'SS' },
    0x1D5BE: { c: 'e', f: 'SS' },
    0x1D5BF: { c: 'f', f: 'SS' },
    0x1D5C0: { c: 'g', f: 'SS' },
    0x1D5C1: { c: 'h', f: 'SS' },
    0x1D5C2: { c: 'i', f: 'SS' },
    0x1D5C3: { c: 'j', f: 'SS' },
    0x1D5C4: { c: 'k', f: 'SS' },
    0x1D5C5: { c: 'l', f: 'SS' },
    0x1D5C6: { c: 'm', f: 'SS' },
    0x1D5C7: { c: 'n', f: 'SS' },
    0x1D5C8: { c: 'o', f: 'SS' },
    0x1D5C9: { c: 'p', f: 'SS' },
    0x1D5CA: { c: 'q', f: 'SS' },
    0x1D5CB: { c: 'r', f: 'SS' },
    0x1D5CC: { c: 's', f: 'SS' },
    0x1D5CD: { c: 't', f: 'SS' },
    0x1D5CE: { c: 'u', f: 'SS' },
    0x1D5CF: { c: 'v', f: 'SS' },
    0x1D5D0: { c: 'w', f: 'SS' },
    0x1D5D1: { c: 'x', f: 'SS' },
    0x1D5D2: { c: 'y', f: 'SS' },
    0x1D5D3: { c: 'z', f: 'SS' },
    0x1D5D4: { c: 'A', f: 'SSB' },
    0x1D5D5: { c: 'B', f: 'SSB' },
    0x1D5D6: { c: 'C', f: 'SSB' },
    0x1D5D7: { c: 'D', f: 'SSB' },
    0x1D5D8: { c: 'E', f: 'SSB' },
    0x1D5D9: { c: 'F', f: 'SSB' },
    0x1D5DA: { c: 'G', f: 'SSB' },
    0x1D5DB: { c: 'H', f: 'SSB' },
    0x1D5DC: { c: 'I', f: 'SSB' },
    0x1D5DD: { c: 'J', f: 'SSB' },
    0x1D5DE: { c: 'K', f: 'SSB' },
    0x1D5DF: { c: 'L', f: 'SSB' },
    0x1D5E0: { c: 'M', f: 'SSB' },
    0x1D5E1: { c: 'N', f: 'SSB' },
    0x1D5E2: { c: 'O', f: 'SSB' },
    0x1D5E3: { c: 'P', f: 'SSB' },
    0x1D5E4: { c: 'Q', f: 'SSB' },
    0x1D5E5: { c: 'R', f: 'SSB' },
    0x1D5E6: { c: 'S', f: 'SSB' },
    0x1D5E7: { c: 'T', f: 'SSB' },
    0x1D5E8: { c: 'U', f: 'SSB' },
    0x1D5E9: { c: 'V', f: 'SSB' },
    0x1D5EA: { c: 'W', f: 'SSB' },
    0x1D5EB: { c: 'X', f: 'SSB' },
    0x1D5EC: { c: 'Y', f: 'SSB' },
    0x1D5ED: { c: 'Z', f: 'SSB' },
    0x1D5EE: { c: 'a', f: 'SSB' },
    0x1D5EF: { c: 'b', f: 'SSB' },
    0x1D5F0: { c: 'c', f: 'SSB' },
    0x1D5F1: { c: 'd', f: 'SSB' },
    0x1D5F2: { c: 'e', f: 'SSB' },
    0x1D5F3: { c: 'f', f: 'SSB' },
    0x1D5F4: { c: 'g', f: 'SSB' },
    0x1D5F5: { c: 'h', f: 'SSB' },
    0x1D5F6: { c: 'i', f: 'SSB' },
    0x1D5F7: { c: 'j', f: 'SSB' },
    0x1D5F8: { c: 'k', f: 'SSB' },
    0x1D5F9: { c: 'l', f: 'SSB' },
    0x1D5FA: { c: 'm', f: 'SSB' },
    0x1D5FB: { c: 'n', f: 'SSB' },
    0x1D5FC: { c: 'o', f: 'SSB' },
    0x1D5FD: { c: 'p', f: 'SSB' },
    0x1D5FE: { c: 'q', f: 'SSB' },
    0x1D5FF: { c: 'r', f: 'SSB' },
    0x1D600: { c: 's', f: 'SSB' },
    0x1D601: { c: 't', f: 'SSB' },
    0x1D602: { c: 'u', f: 'SSB' },
    0x1D603: { c: 'v', f: 'SSB' },
    0x1D604: { c: 'w', f: 'SSB' },
    0x1D605: { c: 'x', f: 'SSB' },
    0x1D606: { c: 'y', f: 'SSB' },
    0x1D607: { c: 'z', f: 'SSB' },
    0x1D608: { c: 'A', f: 'SSI' },
    0x1D609: { c: 'B', f: 'SSI' },
    0x1D60A: { c: 'C', f: 'SSI' },
    0x1D60B: { c: 'D', f: 'SSI' },
    0x1D60C: { c: 'E', f: 'SSI' },
    0x1D60D: { c: 'F', f: 'SSI' },
    0x1D60E: { c: 'G', f: 'SSI' },
    0x1D60F: { c: 'H', f: 'SSI' },
    0x1D610: { c: 'I', f: 'SSI' },
    0x1D611: { c: 'J', f: 'SSI' },
    0x1D612: { c: 'K', f: 'SSI' },
    0x1D613: { c: 'L', f: 'SSI' },
    0x1D614: { c: 'M', f: 'SSI' },
    0x1D615: { c: 'N', f: 'SSI' },
    0x1D616: { c: 'O', f: 'SSI' },
    0x1D617: { c: 'P', f: 'SSI' },
    0x1D618: { c: 'Q', f: 'SSI' },
    0x1D619: { c: 'R', f: 'SSI' },
    0x1D61A: { c: 'S', f: 'SSI' },
    0x1D61B: { c: 'T', f: 'SSI' },
    0x1D61C: { c: 'U', f: 'SSI' },
    0x1D61D: { c: 'V', f: 'SSI' },
    0x1D61E: { c: 'W', f: 'SSI' },
    0x1D61F: { c: 'X', f: 'SSI' },
    0x1D620: { c: 'Y', f: 'SSI' },
    0x1D621: { c: 'Z', f: 'SSI' },
    0x1D622: { c: 'a', f: 'SSI' },
    0x1D623: { c: 'b', f: 'SSI' },
    0x1D624: { c: 'c', f: 'SSI' },
    0x1D625: { c: 'd', f: 'SSI' },
    0x1D626: { c: 'e', f: 'SSI' },
    0x1D627: { c: 'f', f: 'SSI' },
    0x1D628: { c: 'g', f: 'SSI' },
    0x1D629: { c: 'h', f: 'SSI' },
    0x1D62A: { c: 'i', f: 'SSI' },
    0x1D62B: { c: 'j', f: 'SSI' },
    0x1D62C: { c: 'k', f: 'SSI' },
    0x1D62D: { c: 'l', f: 'SSI' },
    0x1D62E: { c: 'm', f: 'SSI' },
    0x1D62F: { c: 'n', f: 'SSI' },
    0x1D630: { c: 'o', f: 'SSI' },
    0x1D631: { c: 'p', f: 'SSI' },
    0x1D632: { c: 'q', f: 'SSI' },
    0x1D633: { c: 'r', f: 'SSI' },
    0x1D634: { c: 's', f: 'SSI' },
    0x1D635: { c: 't', f: 'SSI' },
    0x1D636: { c: 'u', f: 'SSI' },
    0x1D637: { c: 'v', f: 'SSI' },
    0x1D638: { c: 'w', f: 'SSI' },
    0x1D639: { c: 'x', f: 'SSI' },
    0x1D63A: { c: 'y', f: 'SSI' },
    0x1D63B: { c: 'z', f: 'SSI' },
    0x1D670: { c: 'A', f: 'T' },
    0x1D671: { c: 'B', f: 'T' },
    0x1D672: { c: 'C', f: 'T' },
    0x1D673: { c: 'D', f: 'T' },
    0x1D674: { c: 'E', f: 'T' },
    0x1D675: { c: 'F', f: 'T' },
    0x1D676: { c: 'G', f: 'T' },
    0x1D677: { c: 'H', f: 'T' },
    0x1D678: { c: 'I', f: 'T' },
    0x1D679: { c: 'J', f: 'T' },
    0x1D67A: { c: 'K', f: 'T' },
    0x1D67B: { c: 'L', f: 'T' },
    0x1D67C: { c: 'M', f: 'T' },
    0x1D67D: { c: 'N', f: 'T' },
    0x1D67E: { c: 'O', f: 'T' },
    0x1D67F: { c: 'P', f: 'T' },
    0x1D680: { c: 'Q', f: 'T' },
    0x1D681: { c: 'R', f: 'T' },
    0x1D682: { c: 'S', f: 'T' },
    0x1D683: { c: 'T', f: 'T' },
    0x1D684: { c: 'U', f: 'T' },
    0x1D685: { c: 'V', f: 'T' },
    0x1D686: { c: 'W', f: 'T' },
    0x1D687: { c: 'X', f: 'T' },
    0x1D688: { c: 'Y', f: 'T' },
    0x1D689: { c: 'Z', f: 'T' },
    0x1D68A: { c: 'a', f: 'T' },
    0x1D68B: { c: 'b', f: 'T' },
    0x1D68C: { c: 'c', f: 'T' },
    0x1D68D: { c: 'd', f: 'T' },
    0x1D68E: { c: 'e', f: 'T' },
    0x1D68F: { c: 'f', f: 'T' },
    0x1D690: { c: 'g', f: 'T' },
    0x1D691: { c: 'h', f: 'T' },
    0x1D692: { c: 'i', f: 'T' },
    0x1D693: { c: 'j', f: 'T' },
    0x1D694: { c: 'k', f: 'T' },
    0x1D695: { c: 'l', f: 'T' },
    0x1D696: { c: 'm', f: 'T' },
    0x1D697: { c: 'n', f: 'T' },
    0x1D698: { c: 'o', f: 'T' },
    0x1D699: { c: 'p', f: 'T' },
    0x1D69A: { c: 'q', f: 'T' },
    0x1D69B: { c: 'r', f: 'T' },
    0x1D69C: { c: 's', f: 'T' },
    0x1D69D: { c: 't', f: 'T' },
    0x1D69E: { c: 'u', f: 'T' },
    0x1D69F: { c: 'v', f: 'T' },
    0x1D6A0: { c: 'w', f: 'T' },
    0x1D6A1: { c: 'x', f: 'T' },
    0x1D6A2: { c: 'y', f: 'T' },
    0x1D6A3: { c: 'z', f: 'T' },
    0x1D6A8: { c: 'A', f: 'B' },
    0x1D6A9: { c: 'B', f: 'B' },
    0x1D6AA: { c: '\\393', f: 'B' },
    0x1D6AB: { c: '\\394', f: 'B' },
    0x1D6AC: { c: 'E', f: 'B' },
    0x1D6AD: { c: 'Z', f: 'B' },
    0x1D6AE: { c: 'H', f: 'B' },
    0x1D6AF: { c: '\\398', f: 'B' },
    0x1D6B0: { c: 'I', f: 'B' },
    0x1D6B1: { c: 'K', f: 'B' },
    0x1D6B2: { c: '\\39B', f: 'B' },
    0x1D6B3: { c: 'M', f: 'B' },
    0x1D6B4: { c: 'N', f: 'B' },
    0x1D6B5: { c: '\\39E', f: 'B' },
    0x1D6B6: { c: 'O', f: 'B' },
    0x1D6B7: { c: '\\3A0', f: 'B' },
    0x1D6B8: { c: 'P', f: 'B' },
    0x1D6BA: { c: '\\3A3', f: 'B' },
    0x1D6BB: { c: 'T', f: 'B' },
    0x1D6BC: { c: '\\3A5', f: 'B' },
    0x1D6BD: { c: '\\3A6', f: 'B' },
    0x1D6BE: { c: 'X', f: 'B' },
    0x1D6BF: { c: '\\3A8', f: 'B' },
    0x1D6C0: { c: '\\3A9', f: 'B' },
    0x1D6C1: { c: '\\2207', f: 'B' },
    0x1D6E2: { c: 'A', f: 'I' },
    0x1D6E3: { c: 'B', f: 'I' },
    0x1D6E4: { c: '\\393', f: 'I' },
    0x1D6E5: { c: '\\394', f: 'I' },
    0x1D6E6: { c: 'E', f: 'I' },
    0x1D6E7: { c: 'Z', f: 'I' },
    0x1D6E8: { c: 'H', f: 'I' },
    0x1D6E9: { c: '\\398', f: 'I' },
    0x1D6EA: { c: 'I', f: 'I' },
    0x1D6EB: { c: 'K', f: 'I' },
    0x1D6EC: { c: '\\39B', f: 'I' },
    0x1D6ED: { c: 'M', f: 'I' },
    0x1D6EE: { c: 'N', f: 'I' },
    0x1D6EF: { c: '\\39E', f: 'I' },
    0x1D6F0: { c: 'O', f: 'I' },
    0x1D6F1: { c: '\\3A0', f: 'I' },
    0x1D6F2: { c: 'P', f: 'I' },
    0x1D6F4: { c: '\\3A3', f: 'I' },
    0x1D6F5: { c: 'T', f: 'I' },
    0x1D6F6: { c: '\\3A5', f: 'I' },
    0x1D6F7: { c: '\\3A6', f: 'I' },
    0x1D6F8: { c: 'X', f: 'I' },
    0x1D6F9: { c: '\\3A8', f: 'I' },
    0x1D6FA: { c: '\\3A9', f: 'I' },
    0x1D6FC: { c: '\\3B1', f: 'I' },
    0x1D6FD: { c: '\\3B2', f: 'I' },
    0x1D6FE: { c: '\\3B3', f: 'I' },
    0x1D6FF: { c: '\\3B4', f: 'I' },
    0x1D700: { c: '\\3B5', f: 'I' },
    0x1D701: { c: '\\3B6', f: 'I' },
    0x1D702: { c: '\\3B7', f: 'I' },
    0x1D703: { c: '\\3B8', f: 'I' },
    0x1D704: { c: '\\3B9', f: 'I' },
    0x1D705: { c: '\\3BA', f: 'I' },
    0x1D706: { c: '\\3BB', f: 'I' },
    0x1D707: { c: '\\3BC', f: 'I' },
    0x1D708: { c: '\\3BD', f: 'I' },
    0x1D709: { c: '\\3BE', f: 'I' },
    0x1D70A: { c: '\\3BF', f: 'I' },
    0x1D70B: { c: '\\3C0', f: 'I' },
    0x1D70C: { c: '\\3C1', f: 'I' },
    0x1D70D: { c: '\\3C2', f: 'I' },
    0x1D70E: { c: '\\3C3', f: 'I' },
    0x1D70F: { c: '\\3C4', f: 'I' },
    0x1D710: { c: '\\3C5', f: 'I' },
    0x1D711: { c: '\\3C6', f: 'I' },
    0x1D712: { c: '\\3C7', f: 'I' },
    0x1D713: { c: '\\3C8', f: 'I' },
    0x1D714: { c: '\\3C9', f: 'I' },
    0x1D715: { c: '\\2202' },
    0x1D716: { c: '\\3F5', f: 'I' },
    0x1D717: { c: '\\3D1', f: 'I' },
    0x1D718: { c: '\\E009', f: 'A' },
    0x1D719: { c: '\\3D5', f: 'I' },
    0x1D71A: { c: '\\3F1', f: 'I' },
    0x1D71B: { c: '\\3D6', f: 'I' },
    0x1D71C: { c: 'A', f: 'BI' },
    0x1D71D: { c: 'B', f: 'BI' },
    0x1D71E: { c: '\\393', f: 'BI' },
    0x1D71F: { c: '\\394', f: 'BI' },
    0x1D720: { c: 'E', f: 'BI' },
    0x1D721: { c: 'Z', f: 'BI' },
    0x1D722: { c: 'H', f: 'BI' },
    0x1D723: { c: '\\398', f: 'BI' },
    0x1D724: { c: 'I', f: 'BI' },
    0x1D725: { c: 'K', f: 'BI' },
    0x1D726: { c: '\\39B', f: 'BI' },
    0x1D727: { c: 'M', f: 'BI' },
    0x1D728: { c: 'N', f: 'BI' },
    0x1D729: { c: '\\39E', f: 'BI' },
    0x1D72A: { c: 'O', f: 'BI' },
    0x1D72B: { c: '\\3A0', f: 'BI' },
    0x1D72C: { c: 'P', f: 'BI' },
    0x1D72E: { c: '\\3A3', f: 'BI' },
    0x1D72F: { c: 'T', f: 'BI' },
    0x1D730: { c: '\\3A5', f: 'BI' },
    0x1D731: { c: '\\3A6', f: 'BI' },
    0x1D732: { c: 'X', f: 'BI' },
    0x1D733: { c: '\\3A8', f: 'BI' },
    0x1D734: { c: '\\3A9', f: 'BI' },
    0x1D736: { c: '\\3B1', f: 'BI' },
    0x1D737: { c: '\\3B2', f: 'BI' },
    0x1D738: { c: '\\3B3', f: 'BI' },
    0x1D739: { c: '\\3B4', f: 'BI' },
    0x1D73A: { c: '\\3B5', f: 'BI' },
    0x1D73B: { c: '\\3B6', f: 'BI' },
    0x1D73C: { c: '\\3B7', f: 'BI' },
    0x1D73D: { c: '\\3B8', f: 'BI' },
    0x1D73E: { c: '\\3B9', f: 'BI' },
    0x1D73F: { c: '\\3BA', f: 'BI' },
    0x1D740: { c: '\\3BB', f: 'BI' },
    0x1D741: { c: '\\3BC', f: 'BI' },
    0x1D742: { c: '\\3BD', f: 'BI' },
    0x1D743: { c: '\\3BE', f: 'BI' },
    0x1D744: { c: '\\3BF', f: 'BI' },
    0x1D745: { c: '\\3C0', f: 'BI' },
    0x1D746: { c: '\\3C1', f: 'BI' },
    0x1D747: { c: '\\3C2', f: 'BI' },
    0x1D748: { c: '\\3C3', f: 'BI' },
    0x1D749: { c: '\\3C4', f: 'BI' },
    0x1D74A: { c: '\\3C5', f: 'BI' },
    0x1D74B: { c: '\\3C6', f: 'BI' },
    0x1D74C: { c: '\\3C7', f: 'BI' },
    0x1D74D: { c: '\\3C8', f: 'BI' },
    0x1D74E: { c: '\\3C9', f: 'BI' },
    0x1D74F: { c: '\\2202', f: 'B' },
    0x1D750: { c: '\\3F5', f: 'BI' },
    0x1D751: { c: '\\3D1', f: 'BI' },
    0x1D752: { c: '\\E009', f: 'A' },
    0x1D753: { c: '\\3D5', f: 'BI' },
    0x1D754: { c: '\\3F1', f: 'BI' },
    0x1D755: { c: '\\3D6', f: 'BI' },
    0x1D756: { c: 'A', f: 'SSB' },
    0x1D757: { c: 'B', f: 'SSB' },
    0x1D758: { c: '\\393', f: 'SSB' },
    0x1D759: { c: '\\394', f: 'SSB' },
    0x1D75A: { c: 'E', f: 'SSB' },
    0x1D75B: { c: 'Z', f: 'SSB' },
    0x1D75C: { c: 'H', f: 'SSB' },
    0x1D75D: { c: '\\398', f: 'SSB' },
    0x1D75E: { c: 'I', f: 'SSB' },
    0x1D75F: { c: 'K', f: 'SSB' },
    0x1D760: { c: '\\39B', f: 'SSB' },
    0x1D761: { c: 'M', f: 'SSB' },
    0x1D762: { c: 'N', f: 'SSB' },
    0x1D763: { c: '\\39E', f: 'SSB' },
    0x1D764: { c: 'O', f: 'SSB' },
    0x1D765: { c: '\\3A0', f: 'SSB' },
    0x1D766: { c: 'P', f: 'SSB' },
    0x1D768: { c: '\\3A3', f: 'SSB' },
    0x1D769: { c: 'T', f: 'SSB' },
    0x1D76A: { c: '\\3A5', f: 'SSB' },
    0x1D76B: { c: '\\3A6', f: 'SSB' },
    0x1D76C: { c: 'X', f: 'SSB' },
    0x1D76D: { c: '\\3A8', f: 'SSB' },
    0x1D76E: { c: '\\3A9', f: 'SSB' },
    0x1D7CE: { c: '0', f: 'B' },
    0x1D7CF: { c: '1', f: 'B' },
    0x1D7D0: { c: '2', f: 'B' },
    0x1D7D1: { c: '3', f: 'B' },
    0x1D7D2: { c: '4', f: 'B' },
    0x1D7D3: { c: '5', f: 'B' },
    0x1D7D4: { c: '6', f: 'B' },
    0x1D7D5: { c: '7', f: 'B' },
    0x1D7D6: { c: '8', f: 'B' },
    0x1D7D7: { c: '9', f: 'B' },
    0x1D7E2: { c: '0', f: 'SS' },
    0x1D7E3: { c: '1', f: 'SS' },
    0x1D7E4: { c: '2', f: 'SS' },
    0x1D7E5: { c: '3', f: 'SS' },
    0x1D7E6: { c: '4', f: 'SS' },
    0x1D7E7: { c: '5', f: 'SS' },
    0x1D7E8: { c: '6', f: 'SS' },
    0x1D7E9: { c: '7', f: 'SS' },
    0x1D7EA: { c: '8', f: 'SS' },
    0x1D7EB: { c: '9', f: 'SS' },
    0x1D7EC: { c: '0', f: 'SSB' },
    0x1D7ED: { c: '1', f: 'SSB' },
    0x1D7EE: { c: '2', f: 'SSB' },
    0x1D7EF: { c: '3', f: 'SSB' },
    0x1D7F0: { c: '4', f: 'SSB' },
    0x1D7F1: { c: '5', f: 'SSB' },
    0x1D7F2: { c: '6', f: 'SSB' },
    0x1D7F3: { c: '7', f: 'SSB' },
    0x1D7F4: { c: '8', f: 'SSB' },
    0x1D7F5: { c: '9', f: 'SSB' },
    0x1D7F6: { c: '0', f: 'T' },
    0x1D7F7: { c: '1', f: 'T' },
    0x1D7F8: { c: '2', f: 'T' },
    0x1D7F9: { c: '3', f: 'T' },
    0x1D7FA: { c: '4', f: 'T' },
    0x1D7FB: { c: '5', f: 'T' },
    0x1D7FC: { c: '6', f: 'T' },
    0x1D7FD: { c: '7', f: 'T' },
    0x1D7FE: { c: '8', f: 'T' },
    0x1D7FF: { c: '9', f: 'T' },
});
//# sourceMappingURL=normal.js.map

/***/ }),

/***/ 92644:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifBoldItalic = void 0;
var FontData_js_1 = __webpack_require__(32778);
var sans_serif_bold_italic_js_1 = __webpack_require__(75656);
exports.sansSerifBoldItalic = (0, FontData_js_1.AddCSS)(sans_serif_bold_italic_js_1.sansSerifBoldItalic, {
    0x131: { f: 'SSB' },
    0x237: { f: 'SSB' },
});
//# sourceMappingURL=sans-serif-bold-italic.js.map

/***/ }),

/***/ 7727:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifBold = void 0;
var FontData_js_1 = __webpack_require__(32778);
var sans_serif_bold_js_1 = __webpack_require__(86650);
exports.sansSerifBold = (0, FontData_js_1.AddCSS)(sans_serif_bold_js_1.sansSerifBold, {
    0x2015: { c: '\\2014' },
    0x2017: { c: '_' },
    0x2044: { c: '/' },
    0x2206: { c: '\\394' },
});
//# sourceMappingURL=sans-serif-bold.js.map

/***/ }),

/***/ 62099:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifItalic = void 0;
var FontData_js_1 = __webpack_require__(32778);
var sans_serif_italic_js_1 = __webpack_require__(51518);
exports.sansSerifItalic = (0, FontData_js_1.AddCSS)(sans_serif_italic_js_1.sansSerifItalic, {
    0x391: { c: 'A' },
    0x392: { c: 'B' },
    0x395: { c: 'E' },
    0x396: { c: 'Z' },
    0x397: { c: 'H' },
    0x399: { c: 'I' },
    0x39A: { c: 'K' },
    0x39C: { c: 'M' },
    0x39D: { c: 'N' },
    0x39F: { c: 'O' },
    0x3A1: { c: 'P' },
    0x3A4: { c: 'T' },
    0x3A7: { c: 'X' },
    0x2015: { c: '\\2014' },
    0x2017: { c: '_' },
    0x2044: { c: '/' },
    0x2206: { c: '\\394' },
});
//# sourceMappingURL=sans-serif-italic.js.map

/***/ }),

/***/ 53185:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerif = void 0;
var FontData_js_1 = __webpack_require__(32778);
var sans_serif_js_1 = __webpack_require__(71290);
exports.sansSerif = (0, FontData_js_1.AddCSS)(sans_serif_js_1.sansSerif, {
    0x391: { c: 'A' },
    0x392: { c: 'B' },
    0x395: { c: 'E' },
    0x396: { c: 'Z' },
    0x397: { c: 'H' },
    0x399: { c: 'I' },
    0x39A: { c: 'K' },
    0x39C: { c: 'M' },
    0x39D: { c: 'N' },
    0x39F: { c: 'O' },
    0x3A1: { c: 'P' },
    0x3A4: { c: 'T' },
    0x3A7: { c: 'X' },
    0x2015: { c: '\\2014' },
    0x2017: { c: '_' },
    0x2044: { c: '/' },
    0x2206: { c: '\\394' },
});
//# sourceMappingURL=sans-serif.js.map

/***/ }),

/***/ 95662:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.scriptBold = void 0;
var script_bold_js_1 = __webpack_require__(89520);
Object.defineProperty(exports, "scriptBold", ({ enumerable: true, get: function () { return script_bold_js_1.scriptBold; } }));
//# sourceMappingURL=script-bold.js.map

/***/ }),

/***/ 72569:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.script = void 0;
var script_js_1 = __webpack_require__(10970);
Object.defineProperty(exports, "script", ({ enumerable: true, get: function () { return script_js_1.script; } }));
//# sourceMappingURL=script.js.map

/***/ }),

/***/ 57366:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.smallop = void 0;
var FontData_js_1 = __webpack_require__(32778);
var smallop_js_1 = __webpack_require__(34664);
exports.smallop = (0, FontData_js_1.AddCSS)(smallop_js_1.smallop, {
    0x2044: { c: '/' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x2758: { c: '\\2223' },
    0x2A0C: { c: '\\222C\\222C' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
});
//# sourceMappingURL=smallop.js.map

/***/ }),

/***/ 99983:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texCalligraphicBold = void 0;
var FontData_js_1 = __webpack_require__(32778);
var tex_calligraphic_bold_js_1 = __webpack_require__(43545);
exports.texCalligraphicBold = (0, FontData_js_1.AddCSS)(tex_calligraphic_bold_js_1.texCalligraphicBold, {
    0x131: { f: 'B' },
    0x237: { f: 'B' },
});
//# sourceMappingURL=tex-calligraphic-bold.js.map

/***/ }),

/***/ 90769:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texCalligraphic = void 0;
var tex_calligraphic_js_1 = __webpack_require__(55389);
Object.defineProperty(exports, "texCalligraphic", ({ enumerable: true, get: function () { return tex_calligraphic_js_1.texCalligraphic; } }));
//# sourceMappingURL=tex-calligraphic.js.map

/***/ }),

/***/ 71370:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texMathit = void 0;
var tex_mathit_js_1 = __webpack_require__(1641);
Object.defineProperty(exports, "texMathit", ({ enumerable: true, get: function () { return tex_mathit_js_1.texMathit; } }));
//# sourceMappingURL=tex-mathit.js.map

/***/ }),

/***/ 17425:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texOldstyleBold = void 0;
var tex_oldstyle_bold_js_1 = __webpack_require__(33944);
Object.defineProperty(exports, "texOldstyleBold", ({ enumerable: true, get: function () { return tex_oldstyle_bold_js_1.texOldstyleBold; } }));
//# sourceMappingURL=tex-oldstyle-bold.js.map

/***/ }),

/***/ 41771:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texOldstyle = void 0;
var tex_oldstyle_js_1 = __webpack_require__(98690);
Object.defineProperty(exports, "texOldstyle", ({ enumerable: true, get: function () { return tex_oldstyle_js_1.texOldstyle; } }));
//# sourceMappingURL=tex-oldstyle.js.map

/***/ }),

/***/ 21055:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texSize3 = void 0;
var FontData_js_1 = __webpack_require__(32778);
var tex_size3_js_1 = __webpack_require__(96464);
exports.texSize3 = (0, FontData_js_1.AddCSS)(tex_size3_js_1.texSize3, {
    0x2044: { c: '/' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
});
//# sourceMappingURL=tex-size3.js.map

/***/ }),

/***/ 39534:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texSize4 = void 0;
var FontData_js_1 = __webpack_require__(32778);
var tex_size4_js_1 = __webpack_require__(58905);
exports.texSize4 = (0, FontData_js_1.AddCSS)(tex_size4_js_1.texSize4, {
    0x2044: { c: '/' },
    0x2329: { c: '\\27E8' },
    0x232A: { c: '\\27E9' },
    0x3008: { c: '\\27E8' },
    0x3009: { c: '\\27E9' },
    0xE155: { c: '\\E153\\E152' },
    0xE156: { c: '\\E151\\E150' },
});
//# sourceMappingURL=tex-size4.js.map

/***/ }),

/***/ 27424:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texVariant = void 0;
var FontData_js_1 = __webpack_require__(32778);
var tex_variant_js_1 = __webpack_require__(46527);
exports.texVariant = (0, FontData_js_1.AddCSS)(tex_variant_js_1.texVariant, {
    0x3F0: { c: '\\E009' },
    0x210F: { f: '' },
    0x2224: { c: '\\E006' },
    0x2226: { c: '\\E007' },
    0x2268: { c: '\\E00C' },
    0x2269: { c: '\\E00D' },
    0x2270: { c: '\\E011' },
    0x2271: { c: '\\E00E' },
    0x2288: { c: '\\E016' },
    0x2289: { c: '\\E018' },
    0x228A: { c: '\\E01A' },
    0x228B: { c: '\\E01B' },
    0x2A87: { c: '\\E010' },
    0x2A88: { c: '\\E00F' },
    0x2ACB: { c: '\\E017' },
    0x2ACC: { c: '\\E019' },
});
//# sourceMappingURL=tex-variant.js.map

/***/ }),

/***/ 5549:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __read = (this && this.__read) || function (o, n) {
    var m = typeof Symbol === "function" && o[Symbol.iterator];
    if (!m) return o;
    var i = m.call(o), r, ar = [], e;
    try {
        while ((n === void 0 || n-- > 0) && !(r = i.next()).done) ar.push(r.value);
    }
    catch (error) { e = { error: error }; }
    finally {
        try {
            if (r && !r.done && (m = i["return"])) m.call(i);
        }
        finally { if (e) throw e.error; }
    }
    return ar;
};
var __spreadArray = (this && this.__spreadArray) || function (to, from, pack) {
    if (pack || arguments.length === 2) for (var i = 0, l = from.length, ar; i < l; i++) {
        if (ar || !(i in from)) {
            if (!ar) ar = Array.prototype.slice.call(from, 0, i);
            ar[i] = from[i];
        }
    }
    return to.concat(ar || Array.prototype.slice.call(from));
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.FontData = exports.NOSTRETCH = exports.H = exports.V = void 0;
var Options_js_1 = __webpack_require__(4498);
exports.V = 1;
exports.H = 2;
exports.NOSTRETCH = { dir: 0 };
var FontData = (function () {
    function FontData(options) {
        var e_1, _a, e_2, _b;
        if (options === void 0) { options = null; }
        this.variant = {};
        this.delimiters = {};
        this.cssFontMap = {};
        this.remapChars = {};
        this.skewIcFactor = .75;
        var CLASS = this.constructor;
        this.options = (0, Options_js_1.userOptions)((0, Options_js_1.defaultOptions)({}, CLASS.OPTIONS), options);
        this.params = __assign({}, CLASS.defaultParams);
        this.sizeVariants = __spreadArray([], __read(CLASS.defaultSizeVariants), false);
        this.stretchVariants = __spreadArray([], __read(CLASS.defaultStretchVariants), false);
        this.cssFontMap = __assign({}, CLASS.defaultCssFonts);
        try {
            for (var _c = __values(Object.keys(this.cssFontMap)), _d = _c.next(); !_d.done; _d = _c.next()) {
                var name_1 = _d.value;
                if (this.cssFontMap[name_1][0] === 'unknown') {
                    this.cssFontMap[name_1][0] = this.options.unknownFamily;
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_1) throw e_1.error; }
        }
        this.cssFamilyPrefix = CLASS.defaultCssFamilyPrefix;
        this.createVariants(CLASS.defaultVariants);
        this.defineDelimiters(CLASS.defaultDelimiters);
        try {
            for (var _e = __values(Object.keys(CLASS.defaultChars)), _f = _e.next(); !_f.done; _f = _e.next()) {
                var name_2 = _f.value;
                this.defineChars(name_2, CLASS.defaultChars[name_2]);
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
            }
            finally { if (e_2) throw e_2.error; }
        }
        this.defineRemap('accent', CLASS.defaultAccentMap);
        this.defineRemap('mo', CLASS.defaultMoMap);
        this.defineRemap('mn', CLASS.defaultMnMap);
    }
    FontData.charOptions = function (font, n) {
        var char = font[n];
        if (char.length === 3) {
            char[3] = {};
        }
        return char[3];
    };
    Object.defineProperty(FontData.prototype, "styles", {
        get: function () {
            return this._styles;
        },
        set: function (style) {
            this._styles = style;
        },
        enumerable: false,
        configurable: true
    });
    FontData.prototype.createVariant = function (name, inherit, link) {
        if (inherit === void 0) { inherit = null; }
        if (link === void 0) { link = null; }
        var variant = {
            linked: [],
            chars: (inherit ? Object.create(this.variant[inherit].chars) : {})
        };
        if (link && this.variant[link]) {
            Object.assign(variant.chars, this.variant[link].chars);
            this.variant[link].linked.push(variant.chars);
            variant.chars = Object.create(variant.chars);
        }
        this.remapSmpChars(variant.chars, name);
        this.variant[name] = variant;
    };
    FontData.prototype.remapSmpChars = function (chars, name) {
        var e_3, _a, e_4, _b;
        var CLASS = this.constructor;
        if (CLASS.VariantSmp[name]) {
            var SmpRemap = CLASS.SmpRemap;
            var SmpGreek = [null, null, CLASS.SmpRemapGreekU, CLASS.SmpRemapGreekL];
            try {
                for (var _c = __values(CLASS.SmpRanges), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var _e = __read(_d.value, 3), i = _e[0], lo = _e[1], hi = _e[2];
                    var base = CLASS.VariantSmp[name][i];
                    if (!base)
                        continue;
                    for (var n = lo; n <= hi; n++) {
                        if (n === 0x3A2)
                            continue;
                        var smp = base + n - lo;
                        chars[n] = this.smpChar(SmpRemap[smp] || smp);
                    }
                    if (SmpGreek[i]) {
                        try {
                            for (var _f = (e_4 = void 0, __values(Object.keys(SmpGreek[i]).map(function (x) { return parseInt(x); }))), _g = _f.next(); !_g.done; _g = _f.next()) {
                                var n = _g.value;
                                chars[n] = this.smpChar(base + SmpGreek[i][n]);
                            }
                        }
                        catch (e_4_1) { e_4 = { error: e_4_1 }; }
                        finally {
                            try {
                                if (_g && !_g.done && (_b = _f.return)) _b.call(_f);
                            }
                            finally { if (e_4) throw e_4.error; }
                        }
                    }
                }
            }
            catch (e_3_1) { e_3 = { error: e_3_1 }; }
            finally {
                try {
                    if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
                }
                finally { if (e_3) throw e_3.error; }
            }
        }
        if (name === 'bold') {
            chars[0x3DC] = this.smpChar(0x1D7CA);
            chars[0x3DD] = this.smpChar(0x1D7CB);
        }
    };
    FontData.prototype.smpChar = function (n) {
        return [, , , { smp: n }];
    };
    FontData.prototype.createVariants = function (variants) {
        var e_5, _a;
        try {
            for (var variants_1 = __values(variants), variants_1_1 = variants_1.next(); !variants_1_1.done; variants_1_1 = variants_1.next()) {
                var variant = variants_1_1.value;
                this.createVariant(variant[0], variant[1], variant[2]);
            }
        }
        catch (e_5_1) { e_5 = { error: e_5_1 }; }
        finally {
            try {
                if (variants_1_1 && !variants_1_1.done && (_a = variants_1.return)) _a.call(variants_1);
            }
            finally { if (e_5) throw e_5.error; }
        }
    };
    FontData.prototype.defineChars = function (name, chars) {
        var e_6, _a;
        var variant = this.variant[name];
        Object.assign(variant.chars, chars);
        try {
            for (var _b = __values(variant.linked), _c = _b.next(); !_c.done; _c = _b.next()) {
                var link = _c.value;
                Object.assign(link, chars);
            }
        }
        catch (e_6_1) { e_6 = { error: e_6_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_6) throw e_6.error; }
        }
    };
    FontData.prototype.defineDelimiters = function (delims) {
        Object.assign(this.delimiters, delims);
    };
    FontData.prototype.defineRemap = function (name, remap) {
        if (!this.remapChars.hasOwnProperty(name)) {
            this.remapChars[name] = {};
        }
        Object.assign(this.remapChars[name], remap);
    };
    FontData.prototype.getDelimiter = function (n) {
        return this.delimiters[n];
    };
    FontData.prototype.getSizeVariant = function (n, i) {
        if (this.delimiters[n].variants) {
            i = this.delimiters[n].variants[i];
        }
        return this.sizeVariants[i];
    };
    FontData.prototype.getStretchVariant = function (n, i) {
        return this.stretchVariants[this.delimiters[n].stretchv ? this.delimiters[n].stretchv[i] : 0];
    };
    FontData.prototype.getChar = function (name, n) {
        return this.variant[name].chars[n];
    };
    FontData.prototype.getVariant = function (name) {
        return this.variant[name];
    };
    FontData.prototype.getCssFont = function (variant) {
        return this.cssFontMap[variant] || ['serif', false, false];
    };
    FontData.prototype.getFamily = function (family) {
        return (this.cssFamilyPrefix ? this.cssFamilyPrefix + ', ' + family : family);
    };
    FontData.prototype.getRemappedChar = function (name, c) {
        var map = this.remapChars[name] || {};
        return map[c];
    };
    FontData.OPTIONS = {
        unknownFamily: 'serif'
    };
    FontData.JAX = 'common';
    FontData.NAME = '';
    FontData.defaultVariants = [
        ['normal'],
        ['bold', 'normal'],
        ['italic', 'normal'],
        ['bold-italic', 'italic', 'bold'],
        ['double-struck', 'bold'],
        ['fraktur', 'normal'],
        ['bold-fraktur', 'bold', 'fraktur'],
        ['script', 'italic'],
        ['bold-script', 'bold-italic', 'script'],
        ['sans-serif', 'normal'],
        ['bold-sans-serif', 'bold', 'sans-serif'],
        ['sans-serif-italic', 'italic', 'sans-serif'],
        ['sans-serif-bold-italic', 'bold-italic', 'bold-sans-serif'],
        ['monospace', 'normal']
    ];
    FontData.defaultCssFonts = {
        normal: ['unknown', false, false],
        bold: ['unknown', false, true],
        italic: ['unknown', true, false],
        'bold-italic': ['unknown', true, true],
        'double-struck': ['unknown', false, true],
        fraktur: ['unknown', false, false],
        'bold-fraktur': ['unknown', false, true],
        script: ['cursive', false, false],
        'bold-script': ['cursive', false, true],
        'sans-serif': ['sans-serif', false, false],
        'bold-sans-serif': ['sans-serif', false, true],
        'sans-serif-italic': ['sans-serif', true, false],
        'sans-serif-bold-italic': ['sans-serif', true, true],
        monospace: ['monospace', false, false]
    };
    FontData.defaultCssFamilyPrefix = '';
    FontData.VariantSmp = {
        bold: [0x1D400, 0x1D41A, 0x1D6A8, 0x1D6C2, 0x1D7CE],
        italic: [0x1D434, 0x1D44E, 0x1D6E2, 0x1D6FC],
        'bold-italic': [0x1D468, 0x1D482, 0x1D71C, 0x1D736],
        script: [0x1D49C, 0x1D4B6],
        'bold-script': [0x1D4D0, 0x1D4EA],
        fraktur: [0x1D504, 0x1D51E],
        'double-struck': [0x1D538, 0x1D552, , , 0x1D7D8],
        'bold-fraktur': [0x1D56C, 0x1D586],
        'sans-serif': [0x1D5A0, 0x1D5BA, , , 0x1D7E2],
        'bold-sans-serif': [0x1D5D4, 0x1D5EE, 0x1D756, 0x1D770, 0x1D7EC],
        'sans-serif-italic': [0x1D608, 0x1D622],
        'sans-serif-bold-italic': [0x1D63C, 0x1D656, 0x1D790, 0x1D7AA],
        'monospace': [0x1D670, 0x1D68A, , , 0x1D7F6]
    };
    FontData.SmpRanges = [
        [0, 0x41, 0x5A],
        [1, 0x61, 0x7A],
        [2, 0x391, 0x3A9],
        [3, 0x3B1, 0x3C9],
        [4, 0x30, 0x39]
    ];
    FontData.SmpRemap = {
        0x1D455: 0x210E,
        0x1D49D: 0x212C,
        0x1D4A0: 0x2130,
        0x1D4A1: 0x2131,
        0x1D4A3: 0x210B,
        0x1D4A4: 0x2110,
        0x1D4A7: 0x2112,
        0x1D4A8: 0x2133,
        0x1D4AD: 0x211B,
        0x1D4BA: 0x212F,
        0x1D4BC: 0x210A,
        0x1D4C4: 0x2134,
        0x1D506: 0x212D,
        0x1D50B: 0x210C,
        0x1D50C: 0x2111,
        0x1D515: 0x211C,
        0x1D51D: 0x2128,
        0x1D53A: 0x2102,
        0x1D53F: 0x210D,
        0x1D545: 0x2115,
        0x1D547: 0x2119,
        0x1D548: 0x211A,
        0x1D549: 0x211D,
        0x1D551: 0x2124,
    };
    FontData.SmpRemapGreekU = {
        0x2207: 0x19,
        0x03F4: 0x11
    };
    FontData.SmpRemapGreekL = {
        0x3D1: 0x1B,
        0x3D5: 0x1D,
        0x3D6: 0x1F,
        0x3F0: 0x1C,
        0x3F1: 0x1E,
        0x3F5: 0x1A,
        0x2202: 0x19
    };
    FontData.defaultAccentMap = {
        0x0300: '\u02CB',
        0x0301: '\u02CA',
        0x0302: '\u02C6',
        0x0303: '\u02DC',
        0x0304: '\u02C9',
        0x0306: '\u02D8',
        0x0307: '\u02D9',
        0x0308: '\u00A8',
        0x030A: '\u02DA',
        0x030C: '\u02C7',
        0x2192: '\u20D7',
        0x2032: '\'',
        0x2033: '\'\'',
        0x2034: '\'\'\'',
        0x2035: '`',
        0x2036: '``',
        0x2037: '```',
        0x2057: '\'\'\'\'',
        0x20D0: '\u21BC',
        0x20D1: '\u21C0',
        0x20D6: '\u2190',
        0x20E1: '\u2194',
        0x20F0: '*',
        0x20DB: '...',
        0x20DC: '....',
        0x20EC: '\u21C1',
        0x20ED: '\u21BD',
        0x20EE: '\u2190',
        0x20EF: '\u2192'
    };
    FontData.defaultMoMap = {
        0x002D: '\u2212'
    };
    FontData.defaultMnMap = {
        0x002D: '\u2212'
    };
    FontData.defaultParams = {
        x_height: .442,
        quad: 1,
        num1: .676,
        num2: .394,
        num3: .444,
        denom1: .686,
        denom2: .345,
        sup1: .413,
        sup2: .363,
        sup3: .289,
        sub1: .15,
        sub2: .247,
        sup_drop: .386,
        sub_drop: .05,
        delim1: 2.39,
        delim2: 1.0,
        axis_height: .25,
        rule_thickness: .06,
        big_op_spacing1: .111,
        big_op_spacing2: .167,
        big_op_spacing3: .2,
        big_op_spacing4: .6,
        big_op_spacing5: .1,
        surd_height: .075,
        scriptspace: .05,
        nulldelimiterspace: .12,
        delimiterfactor: 901,
        delimitershortfall: .3,
        min_rule_thickness: 1.25,
        separation_factor: 1.75,
        extra_ic: .033
    };
    FontData.defaultDelimiters = {};
    FontData.defaultChars = {};
    FontData.defaultSizeVariants = [];
    FontData.defaultStretchVariants = [];
    return FontData;
}());
exports.FontData = FontData;
//# sourceMappingURL=FontData.js.map

/***/ }),

/***/ 49308:
/***/ (function(__unused_webpack_module, exports) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __read = (this && this.__read) || function (o, n) {
    var m = typeof Symbol === "function" && o[Symbol.iterator];
    if (!m) return o;
    var i = m.call(o), r, ar = [], e;
    try {
        while ((n === void 0 || n-- > 0) && !(r = i.next()).done) ar.push(r.value);
    }
    catch (error) { e = { error: error }; }
    finally {
        try {
            if (r && !r.done && (m = i["return"])) m.call(i);
        }
        finally { if (e) throw e.error; }
    }
    return ar;
};
var __spreadArray = (this && this.__spreadArray) || function (to, from, pack) {
    if (pack || arguments.length === 2) for (var i = 0, l = from.length, ar; i < l; i++) {
        if (ar || !(i in from)) {
            if (!ar) ar = Array.prototype.slice.call(from, 0, i);
            ar[i] = from[i];
        }
    }
    return to.concat(ar || Array.prototype.slice.call(from));
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.CommonTeXFontMixin = void 0;
function CommonTeXFontMixin(Base) {
    var _a;
    return _a = (function (_super) {
            __extends(class_1, _super);
            function class_1() {
                return _super !== null && _super.apply(this, arguments) || this;
            }
            class_1.prototype.getDelimiterData = function (n) {
                return this.getChar('-smallop', n) || this.getChar('-size4', n);
            };
            return class_1;
        }(Base)),
        _a.NAME = 'TeX',
        _a.defaultVariants = __spreadArray(__spreadArray([], __read(Base.defaultVariants), false), [
            ['-smallop', 'normal'],
            ['-largeop', 'normal'],
            ['-size3', 'normal'],
            ['-size4', 'normal'],
            ['-tex-calligraphic', 'italic'],
            ['-tex-bold-calligraphic', 'bold-italic'],
            ['-tex-oldstyle', 'normal'],
            ['-tex-bold-oldstyle', 'bold'],
            ['-tex-mathit', 'italic'],
            ['-tex-variant', 'normal']
        ], false),
        _a.defaultCssFonts = __assign(__assign({}, Base.defaultCssFonts), { '-smallop': ['serif', false, false], '-largeop': ['serif', false, false], '-size3': ['serif', false, false], '-size4': ['serif', false, false], '-tex-calligraphic': ['cursive', true, false], '-tex-bold-calligraphic': ['cursive', true, true], '-tex-oldstyle': ['serif', false, false], '-tex-bold-oldstyle': ['serif', false, true], '-tex-mathit': ['serif', true, false] }),
        _a.defaultSizeVariants = ['normal', '-smallop', '-largeop', '-size3', '-size4', '-tex-variant'],
        _a.defaultStretchVariants = ['-size4'],
        _a;
}
exports.CommonTeXFontMixin = CommonTeXFontMixin;
//# sourceMappingURL=tex.js.map

/***/ }),

/***/ 95884:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.boldItalic = void 0;
exports.boldItalic = {
    0x2F: [.711, .21, .894],
    0x131: [.452, .008, .394, { sk: .0319 }],
    0x237: [.451, .201, .439, { sk: .0958 }],
    0x2044: [.711, .21, .894],
    0x2206: [.711, 0, .958, { sk: .192 }],
    0x29F8: [.711, .21, .894],
};
//# sourceMappingURL=bold-italic.js.map

/***/ }),

/***/ 75083:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.bold = void 0;
exports.bold = {
    0x21: [.705, 0, .35],
    0x22: [.694, -0.329, .603],
    0x23: [.694, .193, .958],
    0x24: [.75, .056, .575],
    0x25: [.75, .056, .958],
    0x26: [.705, .011, .894],
    0x27: [.694, -0.329, .319],
    0x28: [.75, .249, .447],
    0x29: [.75, .249, .447],
    0x2A: [.75, -0.306, .575],
    0x2B: [.633, .131, .894],
    0x2C: [.171, .194, .319],
    0x2D: [.278, -0.166, .383],
    0x2E: [.171, 0, .319],
    0x2F: [.75, .25, .575],
    0x3A: [.444, 0, .319],
    0x3B: [.444, .194, .319],
    0x3C: [.587, .085, .894],
    0x3D: [.393, -0.109, .894],
    0x3E: [.587, .085, .894],
    0x3F: [.7, 0, .543],
    0x40: [.699, .006, .894],
    0x5B: [.75, .25, .319],
    0x5C: [.75, .25, .575],
    0x5D: [.75, .25, .319],
    0x5E: [.694, -0.52, .575],
    0x5F: [-0.01, .061, .575],
    0x60: [.706, -0.503, .575],
    0x7B: [.75, .25, .575],
    0x7C: [.75, .249, .319],
    0x7D: [.75, .25, .575],
    0x7E: [.344, -0.202, .575],
    0xA8: [.695, -0.535, .575],
    0xAC: [.371, -0.061, .767],
    0xAF: [.607, -0.54, .575],
    0xB0: [.702, -0.536, .575],
    0xB1: [.728, .035, .894],
    0xB4: [.706, -0.503, .575],
    0xB7: [.336, -0.166, .319],
    0xD7: [.53, .028, .894],
    0xF7: [.597, .096, .894],
    0x131: [.442, 0, .278, { sk: .0278 }],
    0x237: [.442, .205, .306, { sk: .0833 }],
    0x2B9: [.563, -0.033, .344],
    0x2C6: [.694, -0.52, .575],
    0x2C7: [.66, -0.515, .575],
    0x2C9: [.607, -0.54, .575],
    0x2CA: [.706, -0.503, .575],
    0x2CB: [.706, -0.503, .575],
    0x2D8: [.694, -0.5, .575],
    0x2D9: [.695, -0.525, .575],
    0x2DA: [.702, -0.536, .575],
    0x2DC: [.694, -0.552, .575],
    0x300: [.706, -0.503, 0],
    0x301: [.706, -0.503, 0],
    0x302: [.694, -0.52, 0],
    0x303: [.694, -0.552, 0],
    0x304: [.607, -0.54, 0],
    0x306: [.694, -0.5, 0],
    0x307: [.695, -0.525, 0],
    0x308: [.695, -0.535, 0],
    0x30A: [.702, -0.536, 0],
    0x30B: [.714, -0.511, 0],
    0x30C: [.66, -0.515, 0],
    0x338: [.711, .21, 0],
    0x2002: [0, 0, .5],
    0x2003: [0, 0, .999],
    0x2004: [0, 0, .333],
    0x2005: [0, 0, .25],
    0x2006: [0, 0, .167],
    0x2009: [0, 0, .167],
    0x200A: [0, 0, .083],
    0x2013: [.3, -0.249, .575],
    0x2014: [.3, -0.249, 1.15],
    0x2015: [.3, -0.249, 1.15],
    0x2016: [.75, .248, .575],
    0x2017: [-0.01, .061, .575],
    0x2018: [.694, -0.329, .319],
    0x2019: [.694, -0.329, .319],
    0x201C: [.694, -0.329, .603],
    0x201D: [.694, -0.329, .603],
    0x2020: [.702, .211, .511],
    0x2021: [.702, .202, .511],
    0x2022: [.474, -0.028, .575],
    0x2026: [.171, 0, 1.295],
    0x2032: [.563, -0.033, .344],
    0x2033: [.563, 0, .688],
    0x2034: [.563, 0, 1.032],
    0x203E: [.607, -0.54, .575],
    0x2044: [.75, .25, .575],
    0x2057: [.563, 0, 1.376],
    0x20D7: [.723, -0.513, .575],
    0x210F: [.694, .008, .668, { sk: -0.0319 }],
    0x2113: [.702, .019, .474, { sk: .128 }],
    0x2118: [.461, .21, .74],
    0x2135: [.694, 0, .703],
    0x2190: [.518, .017, 1.15],
    0x2191: [.694, .193, .575],
    0x2192: [.518, .017, 1.15],
    0x2193: [.694, .194, .575],
    0x2194: [.518, .017, 1.15],
    0x2195: [.767, .267, .575],
    0x2196: [.724, .194, 1.15],
    0x2197: [.724, .193, 1.15],
    0x2198: [.694, .224, 1.15],
    0x2199: [.694, .224, 1.15],
    0x219A: [.711, .21, 1.15],
    0x219B: [.711, .21, 1.15],
    0x21A6: [.518, .017, 1.15],
    0x21A9: [.518, .017, 1.282],
    0x21AA: [.518, .017, 1.282],
    0x21AE: [.711, .21, 1.15],
    0x21BC: [.518, -0.22, 1.15],
    0x21BD: [.281, .017, 1.15],
    0x21C0: [.518, -0.22, 1.15],
    0x21C1: [.281, .017, 1.15],
    0x21CC: [.718, .017, 1.15],
    0x21CD: [.711, .21, 1.15],
    0x21CE: [.711, .21, 1.15],
    0x21CF: [.711, .21, 1.15],
    0x21D0: [.547, .046, 1.15],
    0x21D1: [.694, .193, .703],
    0x21D2: [.547, .046, 1.15],
    0x21D3: [.694, .194, .703],
    0x21D4: [.547, .046, 1.15],
    0x21D5: [.767, .267, .703],
    0x2200: [.694, .016, .639],
    0x2203: [.694, 0, .639],
    0x2204: [.711, .21, .639],
    0x2205: [.767, .073, .575],
    0x2206: [.698, 0, .958],
    0x2208: [.587, .086, .767],
    0x2209: [.711, .21, .767],
    0x220B: [.587, .086, .767],
    0x220C: [.711, .21, .767],
    0x2212: [.281, -0.221, .894],
    0x2213: [.537, .227, .894],
    0x2215: [.75, .25, .575],
    0x2216: [.75, .25, .575],
    0x2217: [.472, -0.028, .575],
    0x2218: [.474, -0.028, .575],
    0x2219: [.474, -0.028, .575],
    0x221A: [.82, .18, .958, { ic: .03 }],
    0x221D: [.451, .008, .894],
    0x221E: [.452, .008, 1.15],
    0x2220: [.714, 0, .722],
    0x2223: [.75, .249, .319],
    0x2224: [.75, .249, .319],
    0x2225: [.75, .248, .575],
    0x2226: [.75, .248, .575],
    0x2227: [.604, .017, .767],
    0x2228: [.604, .016, .767],
    0x2229: [.603, .016, .767],
    0x222A: [.604, .016, .767],
    0x222B: [.711, .211, .569, { ic: .063 }],
    0x223C: [.391, -0.109, .894],
    0x2240: [.583, .082, .319],
    0x2241: [.711, .21, .894],
    0x2243: [.502, 0, .894],
    0x2244: [.711, .21, .894],
    0x2245: [.638, .027, .894],
    0x2247: [.711, .21, .894],
    0x2248: [.524, -0.032, .894],
    0x2249: [.711, .21, .894],
    0x224D: [.533, .032, .894],
    0x2250: [.721, -0.109, .894],
    0x2260: [.711, .21, .894],
    0x2261: [.505, 0, .894],
    0x2262: [.711, .21, .894],
    0x2264: [.697, .199, .894],
    0x2265: [.697, .199, .894],
    0x226A: [.617, .116, 1.15],
    0x226B: [.618, .116, 1.15],
    0x226D: [.711, .21, .894],
    0x226E: [.711, .21, .894],
    0x226F: [.711, .21, .894],
    0x2270: [.711, .21, .894],
    0x2271: [.711, .21, .894],
    0x227A: [.585, .086, .894],
    0x227B: [.586, .086, .894],
    0x2280: [.711, .21, .894],
    0x2281: [.711, .21, .894],
    0x2282: [.587, .085, .894],
    0x2283: [.587, .086, .894],
    0x2284: [.711, .21, .894],
    0x2285: [.711, .21, .894],
    0x2286: [.697, .199, .894],
    0x2287: [.697, .199, .894],
    0x2288: [.711, .21, .894],
    0x2289: [.711, .21, .894],
    0x228E: [.604, .016, .767],
    0x2291: [.697, .199, .894],
    0x2292: [.697, .199, .894],
    0x2293: [.604, 0, .767],
    0x2294: [.604, 0, .767],
    0x2295: [.632, .132, .894],
    0x2296: [.632, .132, .894],
    0x2297: [.632, .132, .894],
    0x2298: [.632, .132, .894],
    0x2299: [.632, .132, .894],
    0x22A2: [.693, 0, .703],
    0x22A3: [.693, 0, .703],
    0x22A4: [.694, 0, .894],
    0x22A5: [.693, 0, .894],
    0x22A8: [.75, .249, .974],
    0x22AC: [.711, .21, .703],
    0x22AD: [.75, .249, .974],
    0x22C4: [.523, .021, .575],
    0x22C5: [.336, -0.166, .319],
    0x22C6: [.502, 0, .575],
    0x22C8: [.54, .039, 1],
    0x22E2: [.711, .21, .894],
    0x22E3: [.711, .21, .894],
    0x22EE: [.951, .029, .319],
    0x22EF: [.336, -0.166, 1.295],
    0x22F1: [.871, -0.101, 1.323],
    0x2308: [.75, .248, .511],
    0x2309: [.75, .248, .511],
    0x230A: [.749, .248, .511],
    0x230B: [.749, .248, .511],
    0x2322: [.405, -0.108, 1.15],
    0x2323: [.392, -0.126, 1.15],
    0x2329: [.75, .249, .447],
    0x232A: [.75, .249, .447],
    0x25B3: [.711, 0, 1.022],
    0x25B5: [.711, 0, 1.022],
    0x25B9: [.54, .039, .575],
    0x25BD: [.5, .21, 1.022],
    0x25BF: [.5, .21, 1.022],
    0x25C3: [.539, .038, .575],
    0x25EF: [.711, .211, 1.15],
    0x2660: [.719, .129, .894],
    0x2661: [.711, .024, .894],
    0x2662: [.719, .154, .894],
    0x2663: [.719, .129, .894],
    0x266D: [.75, .017, .447],
    0x266E: [.741, .223, .447],
    0x266F: [.724, .224, .447],
    0x2758: [.75, .249, .319],
    0x27E8: [.75, .249, .447],
    0x27E9: [.75, .249, .447],
    0x27F5: [.518, .017, 1.805],
    0x27F6: [.518, .017, 1.833],
    0x27F7: [.518, .017, 2.126],
    0x27F8: [.547, .046, 1.868],
    0x27F9: [.547, .046, 1.87],
    0x27FA: [.547, .046, 2.126],
    0x27FC: [.518, .017, 1.833],
    0x29F8: [.711, .21, .894],
    0x2A2F: [.53, .028, .894],
    0x2A3F: [.686, 0, .9],
    0x2AAF: [.696, .199, .894],
    0x2AB0: [.697, .199, .894],
    0x3008: [.75, .249, .447],
    0x3009: [.75, .249, .447],
};
//# sourceMappingURL=bold.js.map

/***/ }),

/***/ 46922:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.delimiters = exports.VSIZES = exports.HDW3 = exports.HDW2 = exports.HDW1 = void 0;
var FontData_js_1 = __webpack_require__(5549);
exports.HDW1 = [.75, .25, .875];
exports.HDW2 = [.85, .349, .667];
exports.HDW3 = [.583, .082, .5];
exports.VSIZES = [1, 1.2, 1.8, 2.4, 3];
var DELIM2F = { c: 0x2F, dir: FontData_js_1.V, sizes: exports.VSIZES };
var DELIMAF = { c: 0xAF, dir: FontData_js_1.H, sizes: [.5], stretch: [0, 0xAF], HDW: [.59, -0.544, .5] };
var DELIM2C6 = { c: 0x2C6, dir: FontData_js_1.H, sizes: [.5, .556, 1, 1.444, 1.889] };
var DELIM2DC = { c: 0x2DC, dir: FontData_js_1.H, sizes: [.5, .556, 1, 1.444, 1.889] };
var DELIM2013 = { c: 0x2013, dir: FontData_js_1.H, sizes: [.5], stretch: [0, 0x2013], HDW: [.285, -0.248, .5] };
var DELIM2190 = { c: 0x2190, dir: FontData_js_1.H, sizes: [1], stretch: [0x2190, 0x2212], HDW: exports.HDW3 };
var DELIM2192 = { c: 0x2192, dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x2212, 0x2192], HDW: exports.HDW3 };
var DELIM2194 = { c: 0x2194, dir: FontData_js_1.H, sizes: [1], stretch: [0x2190, 0x2212, 0x2192], HDW: exports.HDW3 };
var DELIM21A4 = { c: 0x21A4, dir: FontData_js_1.H, stretch: [0x2190, 0x2212, 0x2223], HDW: exports.HDW3, min: 1.278 };
var DELIM21A6 = { c: 0x21A6, dir: FontData_js_1.H, sizes: [1], stretch: [0x2223, 0x2212, 0x2192], HDW: exports.HDW3 };
var DELIM21D0 = { c: 0x21D0, dir: FontData_js_1.H, sizes: [1], stretch: [0x21D0, 0x3D], HDW: exports.HDW3 };
var DELIM21D2 = { c: 0x21D2, dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x3D, 0x21D2], HDW: exports.HDW3 };
var DELIM21D4 = { c: 0x21D4, dir: FontData_js_1.H, sizes: [1], stretch: [0x21D0, 0x3D, 0x21D2], HDW: exports.HDW3 };
var DELIM2212 = { c: 0x2212, dir: FontData_js_1.H, sizes: [.778], stretch: [0, 0x2212], HDW: exports.HDW3 };
var DELIM2223 = { c: 0x2223, dir: FontData_js_1.V, sizes: [1], stretch: [0, 0x2223], HDW: [.627, .015, .333] };
var DELIM23DC = { c: 0x23DC, dir: FontData_js_1.H, sizes: [.778, 1], schar: [0x2322, 0x2322], variants: [5, 0],
    stretch: [0xE150, 0xE154, 0xE151], HDW: [.32, .2, .5] };
var DELIM23DD = { c: 0x23DD, dir: FontData_js_1.H, sizes: [.778, 1], schar: [0x2323, 0x2323], variants: [5, 0],
    stretch: [0xE152, 0xE154, 0xE153], HDW: [.32, .2, .5] };
var DELIM23DE = { c: 0x23DE, dir: FontData_js_1.H, stretch: [0xE150, 0xE154, 0xE151, 0xE155], HDW: [.32, .2, .5], min: 1.8 };
var DELIM23DF = { c: 0x23DF, dir: FontData_js_1.H, stretch: [0xE152, 0xE154, 0xE153, 0xE156], HDW: [.32, .2, .5], min: 1.8 };
var DELIM27E8 = { c: 0x27E8, dir: FontData_js_1.V, sizes: exports.VSIZES };
var DELIM27E9 = { c: 0x27E9, dir: FontData_js_1.V, sizes: exports.VSIZES };
var DELIM2906 = { c: 0x2906, dir: FontData_js_1.H, stretch: [0x21D0, 0x3D, 0x2223], HDW: exports.HDW3, min: 1.278 };
var DELIM2907 = { c: 0x2907, dir: FontData_js_1.H, stretch: [0x22A8, 0x3D, 0x21D2], HDW: exports.HDW3, min: 1.278 };
exports.delimiters = {
    0x28: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x239B, 0x239C, 0x239D], HDW: [.85, .349, .875] },
    0x29: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x239E, 0x239F, 0x23A0], HDW: [.85, .349, .875] },
    0x2D: DELIM2212,
    0x2F: DELIM2F,
    0x3D: { dir: FontData_js_1.H, sizes: [.778], stretch: [0, 0x3D], HDW: exports.HDW3 },
    0x5B: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23A1, 0x23A2, 0x23A3], HDW: exports.HDW2 },
    0x5C: { dir: FontData_js_1.V, sizes: exports.VSIZES },
    0x5D: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23A4, 0x23A5, 0x23A6], HDW: exports.HDW2 },
    0x5E: DELIM2C6,
    0x5F: DELIM2013,
    0x7B: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23A7, 0x23AA, 0x23A9, 0x23A8], HDW: [.85, .349, .889] },
    0x7C: { dir: FontData_js_1.V, sizes: [1], stretch: [0, 0x2223], HDW: [.75, .25, .333] },
    0x7D: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23AB, 0x23AA, 0x23AD, 0x23AC], HDW: [.85, .349, .889] },
    0x7E: DELIM2DC,
    0xAF: DELIMAF,
    0x2C6: DELIM2C6,
    0x2C9: DELIMAF,
    0x2DC: DELIM2DC,
    0x302: DELIM2C6,
    0x303: DELIM2DC,
    0x332: DELIM2013,
    0x2013: DELIM2013,
    0x2014: DELIM2013,
    0x2015: DELIM2013,
    0x2016: { dir: FontData_js_1.V, sizes: [.602, 1], schar: [0, 0x2225], variants: [1, 0], stretch: [0, 0x2225], HDW: [.602, 0, .556] },
    0x2017: DELIM2013,
    0x203E: DELIMAF,
    0x20D7: DELIM2192,
    0x2190: DELIM2190,
    0x2191: { dir: FontData_js_1.V, sizes: [.888], stretch: [0x2191, 0x23D0], HDW: [.6, 0, .667] },
    0x2192: DELIM2192,
    0x2193: { dir: FontData_js_1.V, sizes: [.888], stretch: [0, 0x23D0, 0x2193], HDW: [.6, 0, .667] },
    0x2194: DELIM2194,
    0x2195: { dir: FontData_js_1.V, sizes: [1.044], stretch: [0x2191, 0x23D0, 0x2193], HDW: exports.HDW1 },
    0x219E: { dir: FontData_js_1.H, sizes: [1], stretch: [0x219E, 0x2212], HDW: exports.HDW3 },
    0x21A0: { dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x2212, 0x21A0], HDW: exports.HDW3 },
    0x21A4: DELIM21A4,
    0x21A5: { dir: FontData_js_1.V, stretch: [0x2191, 0x23D0, 0x22A5], HDW: exports.HDW1, min: 1.555 },
    0x21A6: DELIM21A6,
    0x21A7: { dir: FontData_js_1.V, stretch: [0x22A4, 0x23D0, 0x2193], HDW: exports.HDW1, min: 1.555 },
    0x21B0: { dir: FontData_js_1.V, sizes: [.722], stretch: [0x21B0, 0x23D0], HDW: exports.HDW1 },
    0x21B1: { dir: FontData_js_1.V, sizes: [.722], stretch: [0x21B1, 0x23D0], HDW: exports.HDW1 },
    0x21BC: { dir: FontData_js_1.H, sizes: [1], stretch: [0x21BC, 0x2212], HDW: exports.HDW3 },
    0x21BD: { dir: FontData_js_1.H, sizes: [1], stretch: [0x21BD, 0x2212], HDW: exports.HDW3 },
    0x21BE: { dir: FontData_js_1.V, sizes: [.888], stretch: [0x21BE, 0x23D0], HDW: exports.HDW1 },
    0x21BF: { dir: FontData_js_1.V, sizes: [.888], stretch: [0x21BF, 0x23D0], HDW: exports.HDW1 },
    0x21C0: { dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x2212, 0x21C0], HDW: exports.HDW3 },
    0x21C1: { dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x2212, 0x21C1], HDW: exports.HDW3 },
    0x21C2: { dir: FontData_js_1.V, sizes: [.888], stretch: [0, 0x23D0, 0x21C2], HDW: exports.HDW1 },
    0x21C3: { dir: FontData_js_1.V, sizes: [.888], stretch: [0, 0x23D0, 0x21C3], HDW: exports.HDW1 },
    0x21D0: DELIM21D0,
    0x21D1: { dir: FontData_js_1.V, sizes: [.888], stretch: [0x21D1, 0x2016], HDW: [.599, 0, .778] },
    0x21D2: DELIM21D2,
    0x21D3: { dir: FontData_js_1.V, sizes: [.888], stretch: [0, 0x2016, 0x21D3], HDW: [.6, 0, .778] },
    0x21D4: DELIM21D4,
    0x21D5: { dir: FontData_js_1.V, sizes: [1.044], stretch: [0x21D1, 0x2016, 0x21D3], HDW: [.75, .25, .778] },
    0x21DA: { dir: FontData_js_1.H, sizes: [1], stretch: [0x21DA, 0x2261], HDW: [.464, -0.036, .5] },
    0x21DB: { dir: FontData_js_1.H, sizes: [1], stretch: [0, 0x2261, 0x21DB], HDW: [.464, -0.036, .5] },
    0x2212: DELIM2212,
    0x2215: DELIM2F,
    0x221A: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0xE001, 0xE000, 0x23B7], fullExt: [.65, 2.3], HDW: [.85, .35, 1.056] },
    0x2223: DELIM2223,
    0x2225: { dir: FontData_js_1.V, sizes: [1], stretch: [0, 0x2225], HDW: [.627, .015, .556] },
    0x2308: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23A1, 0x23A2], HDW: exports.HDW2 },
    0x2309: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0x23A4, 0x23A5], HDW: exports.HDW2 },
    0x230A: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0, 0x23A2, 0x23A3], HDW: exports.HDW2 },
    0x230B: { dir: FontData_js_1.V, sizes: exports.VSIZES, stretch: [0, 0x23A5, 0x23A6], HDW: exports.HDW2 },
    0x2312: DELIM23DC,
    0x2322: DELIM23DC,
    0x2323: DELIM23DD,
    0x2329: DELIM27E8,
    0x232A: DELIM27E9,
    0x23AA: { dir: FontData_js_1.V, sizes: [.32], stretch: [0x23AA, 0x23AA, 0x23AA], HDW: [.29, .015, .889] },
    0x23AF: DELIM2013,
    0x23B0: { dir: FontData_js_1.V, sizes: [.989], stretch: [0x23A7, 0x23AA, 0x23AD], HDW: [.75, .25, .889] },
    0x23B1: { dir: FontData_js_1.V, sizes: [.989], stretch: [0x23AB, 0x23AA, 0x23A9], HDW: [.75, .25, .889] },
    0x23B4: { dir: FontData_js_1.H, stretch: [0x250C, 0x2212, 0x2510], HDW: exports.HDW3, min: 1 },
    0x23B5: { dir: FontData_js_1.H, stretch: [0x2514, 0x2212, 0x2518], HDW: exports.HDW3, min: 1 },
    0x23D0: { dir: FontData_js_1.V, sizes: [.602, 1], schar: [0, 0x2223], variants: [1, 0], stretch: [0, 0x2223], HDW: [.602, 0, .333] },
    0x23DC: DELIM23DC,
    0x23DD: DELIM23DD,
    0x23DE: DELIM23DE,
    0x23DF: DELIM23DF,
    0x23E0: { dir: FontData_js_1.H, stretch: [0x2CA, 0x2C9, 0x2CB], HDW: [.59, -0.544, .5], min: 1 },
    0x23E1: { dir: FontData_js_1.H, stretch: [0x2CB, 0x2C9, 0x2CA], HDW: [.59, -0.544, .5], min: 1 },
    0x2500: DELIM2013,
    0x2758: DELIM2223,
    0x27E8: DELIM27E8,
    0x27E9: DELIM27E9,
    0x27EE: { dir: FontData_js_1.V, sizes: [.989], stretch: [0x23A7, 0x23AA, 0x23A9], HDW: [.75, .25, .889] },
    0x27EF: { dir: FontData_js_1.V, sizes: [.989], stretch: [0x23AB, 0x23AA, 0x23AD], HDW: [.75, .25, .889] },
    0x27F5: DELIM2190,
    0x27F6: DELIM2192,
    0x27F7: DELIM2194,
    0x27F8: DELIM21D0,
    0x27F9: DELIM21D2,
    0x27FA: DELIM21D4,
    0x27FB: DELIM21A4,
    0x27FC: DELIM21A6,
    0x27FD: DELIM2906,
    0x27FE: DELIM2907,
    0x2906: DELIM2906,
    0x2907: DELIM2907,
    0x294E: { dir: FontData_js_1.H, stretch: [0x21BC, 0x2212, 0x21C0], HDW: exports.HDW3, min: 2 },
    0x294F: { dir: FontData_js_1.V, stretch: [0x21BE, 0x23D0, 0x21C2], HDW: exports.HDW1, min: 1.776 },
    0x2950: { dir: FontData_js_1.H, stretch: [0x21BD, 0x2212, 0x21C1], HDW: exports.HDW3, min: 2 },
    0x2951: { dir: FontData_js_1.V, stretch: [0x21BF, 0x23D0, 0x21C3], HDW: exports.HDW1, min: .5 },
    0x295A: { dir: FontData_js_1.H, stretch: [0x21BC, 0x2212, 0x2223], HDW: exports.HDW3, min: 1.278 },
    0x295B: { dir: FontData_js_1.H, stretch: [0x2223, 0x2212, 0x21C0], HDW: exports.HDW3, min: 1.278 },
    0x295C: { dir: FontData_js_1.V, stretch: [0x21BE, 0x23D0, 0x22A5], HDW: exports.HDW1, min: 1.556 },
    0x295D: { dir: FontData_js_1.V, stretch: [0x22A4, 0x23D0, 0x21C2], HDW: exports.HDW1, min: 1.556 },
    0x295E: { dir: FontData_js_1.H, stretch: [0x21BD, 0x2212, 0x2223], HDW: exports.HDW3, min: 1.278 },
    0x295F: { dir: FontData_js_1.H, stretch: [0x2223, 0x2212, 0x21C1], HDW: exports.HDW3, min: 1.278 },
    0x2960: { dir: FontData_js_1.V, stretch: [0x21BF, 0x23D0, 0x22A5], HDW: exports.HDW1, min: 1.776 },
    0x2961: { dir: FontData_js_1.V, stretch: [0x22A4, 0x23D0, 0x21C3], HDW: exports.HDW1, min: 1.776 },
    0x3008: DELIM27E8,
    0x3009: DELIM27E9,
    0xFE37: DELIM23DE,
    0xFE38: DELIM23DF,
};
//# sourceMappingURL=delimiters.js.map

/***/ }),

/***/ 94841:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.doubleStruck = void 0;
exports.doubleStruck = {};
//# sourceMappingURL=double-struck.js.map

/***/ }),

/***/ 19173:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.frakturBold = void 0;
exports.frakturBold = {
    0x21: [.689, .012, .349],
    0x22: [.695, -0.432, .254],
    0x26: [.696, .016, .871],
    0x27: [.695, -0.436, .25],
    0x28: [.737, .186, .459],
    0x29: [.735, .187, .459],
    0x2A: [.692, -0.449, .328],
    0x2B: [.598, .082, .893],
    0x2C: [.107, .191, .328],
    0x2D: [.275, -0.236, .893],
    0x2E: [.102, .015, .328],
    0x2F: [.721, .182, .593],
    0x30: [.501, .012, .593],
    0x31: [.489, 0, .593],
    0x32: [.491, 0, .593],
    0x33: [.487, .193, .593],
    0x34: [.495, .196, .593],
    0x35: [.481, .19, .593],
    0x36: [.704, .012, .593],
    0x37: [.479, .197, .593],
    0x38: [.714, .005, .593],
    0x39: [.487, .195, .593],
    0x3A: [.457, .012, .255],
    0x3B: [.458, .19, .255],
    0x3D: [.343, -0.168, .582],
    0x3F: [.697, .014, .428],
    0x5B: [.74, .13, .257],
    0x5D: [.738, .132, .257],
    0x5E: [.734, -0.452, .59],
    0x2018: [.708, -0.411, .254],
    0x2019: [.692, -0.394, .254],
    0x2044: [.721, .182, .593],
    0xE301: [.63, .027, .587],
    0xE302: [.693, .212, .394, { ic: .014 }],
    0xE303: [.681, .219, .387],
    0xE304: [.473, .212, .593],
    0xE305: [.684, .027, .393],
    0xE308: [.679, .22, .981],
    0xE309: [.717, .137, .727],
};
//# sourceMappingURL=fraktur-bold.js.map

/***/ }),

/***/ 57246:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.fraktur = void 0;
exports.fraktur = {
    0x21: [.689, .012, .296],
    0x22: [.695, -0.432, .215],
    0x26: [.698, .011, .738],
    0x27: [.695, -0.436, .212],
    0x28: [.737, .186, .389],
    0x29: [.735, .187, .389],
    0x2A: [.692, -0.449, .278],
    0x2B: [.598, .082, .756],
    0x2C: [.107, .191, .278],
    0x2D: [.275, -0.236, .756],
    0x2E: [.102, .015, .278],
    0x2F: [.721, .182, .502],
    0x30: [.492, .013, .502],
    0x31: [.468, 0, .502],
    0x32: [.474, 0, .502],
    0x33: [.473, .182, .502],
    0x34: [.476, .191, .502],
    0x35: [.458, .184, .502],
    0x36: [.7, .013, .502],
    0x37: [.468, .181, .502],
    0x38: [.705, .01, .502],
    0x39: [.469, .182, .502],
    0x3A: [.457, .012, .216],
    0x3B: [.458, .189, .216],
    0x3D: [.368, -0.132, .756],
    0x3F: [.693, .011, .362],
    0x5B: [.74, .13, .278],
    0x5D: [.738, .131, .278],
    0x5E: [.734, -0.452, .5],
    0x2018: [.708, -0.41, .215],
    0x2019: [.692, -0.395, .215],
    0x2044: [.721, .182, .502],
    0xE300: [.683, .032, .497],
    0xE301: [.616, .03, .498],
    0xE302: [.68, .215, .333],
    0xE303: [.679, .224, .329],
    0xE304: [.471, .214, .503],
    0xE305: [.686, .02, .333],
    0xE306: [.577, .021, .334, { ic: .013 }],
    0xE307: [.475, .022, .501, { ic: .013 }],
};
//# sourceMappingURL=fraktur.js.map

/***/ }),

/***/ 74629:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.italic = void 0;
exports.italic = {
    0x21: [.716, 0, .307, { ic: .073 }],
    0x22: [.694, -0.379, .514, { ic: .024 }],
    0x23: [.694, .194, .818, { ic: .01 }],
    0x25: [.75, .056, .818, { ic: .029 }],
    0x26: [.716, .022, .767, { ic: .035 }],
    0x27: [.694, -0.379, .307, { ic: .07 }],
    0x28: [.75, .25, .409, { ic: .108 }],
    0x29: [.75, .25, .409],
    0x2A: [.75, -0.32, .511, { ic: .073 }],
    0x2B: [.557, .057, .767],
    0x2C: [.121, .194, .307],
    0x2D: [.251, -0.18, .358],
    0x2E: [.121, 0, .307],
    0x2F: [.716, .215, .778],
    0x30: [.665, .021, .511, { ic: .051 }],
    0x31: [.666, 0, .511],
    0x32: [.666, .022, .511, { ic: .04 }],
    0x33: [.666, .022, .511, { ic: .051 }],
    0x34: [.666, .194, .511],
    0x35: [.666, .022, .511, { ic: .056 }],
    0x36: [.665, .022, .511, { ic: .054 }],
    0x37: [.666, .022, .511, { ic: .123 }],
    0x38: [.666, .021, .511, { ic: .042 }],
    0x39: [.666, .022, .511, { ic: .042 }],
    0x3A: [.431, 0, .307],
    0x3B: [.431, .194, .307],
    0x3D: [.367, -0.133, .767],
    0x3F: [.716, 0, .511, { ic: .04 }],
    0x40: [.705, .011, .767, { ic: .022 }],
    0x5B: [.75, .25, .307, { ic: .139 }],
    0x5D: [.75, .25, .307, { ic: .052 }],
    0x5E: [.694, -0.527, .511, { ic: .017 }],
    0x5F: [-0.025, .062, .511, { ic: .043 }],
    0x7E: [.318, -0.208, .511, { ic: .06 }],
    0x131: [.441, .01, .307, { ic: .033 }],
    0x237: [.442, .204, .332],
    0x300: [.697, -0.5, 0],
    0x301: [.697, -0.5, 0, { ic: .039 }],
    0x302: [.694, -0.527, 0, { ic: .017 }],
    0x303: [.668, -0.558, 0, { ic: .06 }],
    0x304: [.589, -0.544, 0, { ic: .054 }],
    0x306: [.694, -0.515, 0, { ic: .062 }],
    0x307: [.669, -0.548, 0],
    0x308: [.669, -0.554, 0, { ic: .045 }],
    0x30A: [.716, -0.542, 0],
    0x30B: [.697, -0.503, 0, { ic: .065 }],
    0x30C: [.638, -0.502, 0, { ic: .029 }],
    0x3DD: [.605, .085, .778],
    0x2013: [.285, -0.248, .511, { ic: .043 }],
    0x2014: [.285, -0.248, 1.022, { ic: .016 }],
    0x2015: [.285, -0.248, 1.022, { ic: .016 }],
    0x2017: [-0.025, .062, .511, { ic: .043 }],
    0x2018: [.694, -0.379, .307, { ic: .055 }],
    0x2019: [.694, -0.379, .307, { ic: .07 }],
    0x201C: [.694, -0.379, .514, { ic: .092 }],
    0x201D: [.694, -0.379, .514, { ic: .024 }],
    0x2044: [.716, .215, .778],
    0x210F: [.695, .013, .54, { ic: .022 }],
    0x2206: [.716, 0, .833, { sk: .167 }],
    0x29F8: [.716, .215, .778],
};
//# sourceMappingURL=italic.js.map

/***/ }),

/***/ 8709:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.largeop = void 0;
exports.largeop = {
    0x28: [1.15, .649, .597],
    0x29: [1.15, .649, .597],
    0x2F: [1.15, .649, .811],
    0x5B: [1.15, .649, .472],
    0x5C: [1.15, .649, .811],
    0x5D: [1.15, .649, .472],
    0x7B: [1.15, .649, .667],
    0x7D: [1.15, .649, .667],
    0x2C6: [.772, -0.565, 1],
    0x2DC: [.75, -0.611, 1],
    0x302: [.772, -0.565, 0],
    0x303: [.75, -0.611, 0],
    0x2016: [.602, 0, .778],
    0x2044: [1.15, .649, .811],
    0x2191: [.6, 0, .667],
    0x2193: [.6, 0, .667],
    0x21D1: [.599, 0, .778],
    0x21D3: [.6, 0, .778],
    0x220F: [.95, .45, 1.278],
    0x2210: [.95, .45, 1.278],
    0x2211: [.95, .45, 1.444],
    0x221A: [1.15, .65, 1, { ic: .02 }],
    0x2223: [.627, .015, .333],
    0x2225: [.627, .015, .556],
    0x222B: [1.36, .862, .556, { ic: .388 }],
    0x222C: [1.36, .862, 1.084, { ic: .388 }],
    0x222D: [1.36, .862, 1.592, { ic: .388 }],
    0x222E: [1.36, .862, .556, { ic: .388 }],
    0x22C0: [.95, .45, 1.111],
    0x22C1: [.95, .45, 1.111],
    0x22C2: [.949, .45, 1.111],
    0x22C3: [.95, .449, 1.111],
    0x2308: [1.15, .649, .528],
    0x2309: [1.15, .649, .528],
    0x230A: [1.15, .649, .528],
    0x230B: [1.15, .649, .528],
    0x2329: [1.15, .649, .611],
    0x232A: [1.15, .649, .611],
    0x23D0: [.602, 0, .667],
    0x2758: [.627, .015, .333],
    0x27E8: [1.15, .649, .611],
    0x27E9: [1.15, .649, .611],
    0x2A00: [.949, .449, 1.511],
    0x2A01: [.949, .449, 1.511],
    0x2A02: [.949, .449, 1.511],
    0x2A04: [.95, .449, 1.111],
    0x2A06: [.95, .45, 1.111],
    0x2A0C: [1.36, .862, 2.168, { ic: .388 }],
    0x3008: [1.15, .649, .611],
    0x3009: [1.15, .649, .611],
};
//# sourceMappingURL=largeop.js.map

/***/ }),

/***/ 60984:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.monospace = void 0;
exports.monospace = {
    0x20: [0, 0, .525],
    0x21: [.622, 0, .525],
    0x22: [.623, -0.333, .525],
    0x23: [.611, 0, .525],
    0x24: [.694, .082, .525],
    0x25: [.694, .083, .525],
    0x26: [.622, .011, .525],
    0x27: [.611, -0.287, .525],
    0x28: [.694, .082, .525],
    0x29: [.694, .082, .525],
    0x2A: [.52, -0.09, .525],
    0x2B: [.531, -0.081, .525],
    0x2C: [.14, .139, .525],
    0x2D: [.341, -0.271, .525],
    0x2E: [.14, 0, .525],
    0x2F: [.694, .083, .525],
    0x3A: [.431, 0, .525],
    0x3B: [.431, .139, .525],
    0x3C: [.557, -0.055, .525],
    0x3D: [.417, -0.195, .525],
    0x3E: [.557, -0.055, .525],
    0x3F: [.617, 0, .525],
    0x40: [.617, .006, .525],
    0x5B: [.694, .082, .525],
    0x5C: [.694, .083, .525],
    0x5D: [.694, .082, .525],
    0x5E: [.611, -0.46, .525],
    0x5F: [-0.025, .095, .525],
    0x60: [.681, -0.357, .525],
    0x7B: [.694, .083, .525],
    0x7C: [.694, .082, .525],
    0x7D: [.694, .083, .525],
    0x7E: [.611, -0.466, .525],
    0x7F: [.612, -0.519, .525],
    0xA0: [0, 0, .525],
    0x131: [.431, 0, .525],
    0x237: [.431, .228, .525],
    0x2B9: [.623, -0.334, .525],
    0x300: [.611, -0.485, 0],
    0x301: [.611, -0.485, 0],
    0x302: [.611, -0.46, 0],
    0x303: [.611, -0.466, 0],
    0x304: [.577, -0.5, 0],
    0x306: [.611, -0.504, 0],
    0x308: [.612, -0.519, 0],
    0x30A: [.619, -0.499, 0],
    0x30C: [.577, -0.449, 0],
    0x391: [.623, 0, .525],
    0x392: [.611, 0, .525],
    0x393: [.611, 0, .525],
    0x394: [.623, 0, .525],
    0x395: [.611, 0, .525],
    0x396: [.611, 0, .525],
    0x397: [.611, 0, .525],
    0x398: [.621, .01, .525],
    0x399: [.611, 0, .525],
    0x39A: [.611, 0, .525],
    0x39B: [.623, 0, .525],
    0x39C: [.611, 0, .525],
    0x39D: [.611, 0, .525],
    0x39E: [.611, 0, .525],
    0x39F: [.621, .01, .525],
    0x3A0: [.611, 0, .525],
    0x3A1: [.611, 0, .525],
    0x3A3: [.611, 0, .525],
    0x3A4: [.611, 0, .525],
    0x3A5: [.622, 0, .525],
    0x3A6: [.611, 0, .525],
    0x3A7: [.611, 0, .525],
    0x3A8: [.611, 0, .525],
    0x3A9: [.622, 0, .525],
    0x2017: [-0.025, .095, .525],
    0x2032: [.623, -0.334, .525],
    0x2033: [.623, 0, 1.05],
    0x2034: [.623, 0, 1.575],
    0x2044: [.694, .083, .525],
    0x2057: [.623, 0, 2.1],
    0x2206: [.623, 0, .525],
};
//# sourceMappingURL=monospace.js.map

/***/ }),

/***/ 65037:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.normal = void 0;
exports.normal = {
    0x20: [0, 0, .25],
    0x21: [.716, 0, .278],
    0x22: [.694, -0.379, .5],
    0x23: [.694, .194, .833],
    0x24: [.75, .056, .5],
    0x25: [.75, .056, .833],
    0x26: [.716, .022, .778],
    0x27: [.694, -0.379, .278],
    0x28: [.75, .25, .389],
    0x29: [.75, .25, .389],
    0x2A: [.75, -0.32, .5],
    0x2B: [.583, .082, .778],
    0x2C: [.121, .194, .278],
    0x2D: [.252, -0.179, .333],
    0x2E: [.12, 0, .278],
    0x2F: [.75, .25, .5],
    0x30: [.666, .022, .5],
    0x31: [.666, 0, .5],
    0x32: [.666, 0, .5],
    0x33: [.665, .022, .5],
    0x34: [.677, 0, .5],
    0x35: [.666, .022, .5],
    0x36: [.666, .022, .5],
    0x37: [.676, .022, .5],
    0x38: [.666, .022, .5],
    0x39: [.666, .022, .5],
    0x3A: [.43, 0, .278],
    0x3B: [.43, .194, .278],
    0x3C: [.54, .04, .778],
    0x3D: [.583, .082, .778],
    0x3E: [.54, .04, .778],
    0x3F: [.705, 0, .472],
    0x40: [.705, .011, .778],
    0x41: [.716, 0, .75],
    0x42: [.683, 0, .708],
    0x43: [.705, .021, .722],
    0x44: [.683, 0, .764],
    0x45: [.68, 0, .681],
    0x46: [.68, 0, .653],
    0x47: [.705, .022, .785],
    0x48: [.683, 0, .75],
    0x49: [.683, 0, .361],
    0x4A: [.683, .022, .514],
    0x4B: [.683, 0, .778],
    0x4C: [.683, 0, .625],
    0x4D: [.683, 0, .917],
    0x4E: [.683, 0, .75],
    0x4F: [.705, .022, .778],
    0x50: [.683, 0, .681],
    0x51: [.705, .193, .778],
    0x52: [.683, .022, .736],
    0x53: [.705, .022, .556],
    0x54: [.677, 0, .722],
    0x55: [.683, .022, .75],
    0x56: [.683, .022, .75],
    0x57: [.683, .022, 1.028],
    0x58: [.683, 0, .75],
    0x59: [.683, 0, .75],
    0x5A: [.683, 0, .611],
    0x5B: [.75, .25, .278],
    0x5C: [.75, .25, .5],
    0x5D: [.75, .25, .278],
    0x5E: [.694, -0.531, .5],
    0x5F: [-0.025, .062, .5],
    0x60: [.699, -0.505, .5],
    0x61: [.448, .011, .5],
    0x62: [.694, .011, .556],
    0x63: [.448, .011, .444],
    0x64: [.694, .011, .556],
    0x65: [.448, .011, .444],
    0x66: [.705, 0, .306, { ic: .066 }],
    0x67: [.453, .206, .5],
    0x68: [.694, 0, .556],
    0x69: [.669, 0, .278],
    0x6A: [.669, .205, .306],
    0x6B: [.694, 0, .528],
    0x6C: [.694, 0, .278],
    0x6D: [.442, 0, .833],
    0x6E: [.442, 0, .556],
    0x6F: [.448, .01, .5],
    0x70: [.442, .194, .556],
    0x71: [.442, .194, .528],
    0x72: [.442, 0, .392],
    0x73: [.448, .011, .394],
    0x74: [.615, .01, .389],
    0x75: [.442, .011, .556],
    0x76: [.431, .011, .528],
    0x77: [.431, .011, .722],
    0x78: [.431, 0, .528],
    0x79: [.431, .204, .528],
    0x7A: [.431, 0, .444],
    0x7B: [.75, .25, .5],
    0x7C: [.75, .249, .278],
    0x7D: [.75, .25, .5],
    0x7E: [.318, -0.215, .5],
    0xA0: [0, 0, .25],
    0xA3: [.714, .011, .769],
    0xA5: [.683, 0, .75],
    0xA8: [.669, -0.554, .5],
    0xAC: [.356, -0.089, .667],
    0xAE: [.709, .175, .947],
    0xAF: [.59, -0.544, .5],
    0xB0: [.715, -0.542, .5],
    0xB1: [.666, 0, .778],
    0xB4: [.699, -0.505, .5],
    0xB7: [.31, -0.19, .278],
    0xD7: [.491, -0.009, .778],
    0xF0: [.749, .021, .556],
    0xF7: [.537, .036, .778],
    0x131: [.442, 0, .278, { sk: .0278 }],
    0x237: [.442, .205, .306, { sk: .0833 }],
    0x2B9: [.56, -0.043, .275],
    0x2C6: [.694, -0.531, .5],
    0x2C7: [.644, -0.513, .5],
    0x2C9: [.59, -0.544, .5],
    0x2CA: [.699, -0.505, .5],
    0x2CB: [.699, -0.505, .5],
    0x2D8: [.694, -0.515, .5],
    0x2D9: [.669, -0.549, .5],
    0x2DA: [.715, -0.542, .5],
    0x2DC: [.668, -0.565, .5],
    0x300: [.699, -0.505, 0],
    0x301: [.699, -0.505, 0],
    0x302: [.694, -0.531, 0],
    0x303: [.668, -0.565, 0],
    0x304: [.59, -0.544, 0],
    0x306: [.694, -0.515, 0],
    0x307: [.669, -0.549, 0],
    0x308: [.669, -0.554, 0],
    0x30A: [.715, -0.542, 0],
    0x30B: [.701, -0.51, 0],
    0x30C: [.644, -0.513, 0],
    0x338: [.716, .215, 0],
    0x391: [.716, 0, .75],
    0x392: [.683, 0, .708],
    0x393: [.68, 0, .625],
    0x394: [.716, 0, .833],
    0x395: [.68, 0, .681],
    0x396: [.683, 0, .611],
    0x397: [.683, 0, .75],
    0x398: [.705, .022, .778],
    0x399: [.683, 0, .361],
    0x39A: [.683, 0, .778],
    0x39B: [.716, 0, .694],
    0x39C: [.683, 0, .917],
    0x39D: [.683, 0, .75],
    0x39E: [.677, 0, .667],
    0x39F: [.705, .022, .778],
    0x3A0: [.68, 0, .75],
    0x3A1: [.683, 0, .681],
    0x3A3: [.683, 0, .722],
    0x3A4: [.677, 0, .722],
    0x3A5: [.705, 0, .778],
    0x3A6: [.683, 0, .722],
    0x3A7: [.683, 0, .75],
    0x3A8: [.683, 0, .778],
    0x3A9: [.704, 0, .722],
    0x2000: [0, 0, .5],
    0x2001: [0, 0, 1],
    0x2002: [0, 0, .5],
    0x2003: [0, 0, 1],
    0x2004: [0, 0, .333],
    0x2005: [0, 0, .25],
    0x2006: [0, 0, .167],
    0x2009: [0, 0, .167],
    0x200A: [0, 0, .1],
    0x200B: [0, 0, 0],
    0x200C: [0, 0, 0],
    0x2013: [.285, -0.248, .5],
    0x2014: [.285, -0.248, 1],
    0x2015: [.285, -0.248, 1],
    0x2016: [.75, .25, .5],
    0x2017: [-0.025, .062, .5],
    0x2018: [.694, -0.379, .278],
    0x2019: [.694, -0.379, .278],
    0x201C: [.694, -0.379, .5],
    0x201D: [.694, -0.379, .5],
    0x2020: [.705, .216, .444],
    0x2021: [.705, .205, .444],
    0x2022: [.444, -0.055, .5],
    0x2026: [.12, 0, 1.172],
    0x2032: [.56, -0.043, .275],
    0x2033: [.56, 0, .55],
    0x2034: [.56, 0, .825],
    0x2035: [.56, -0.043, .275],
    0x2036: [.56, 0, .55],
    0x2037: [.56, 0, .825],
    0x203E: [.59, -0.544, .5],
    0x2044: [.75, .25, .5],
    0x2057: [.56, 0, 1.1],
    0x2060: [0, 0, 0],
    0x2061: [0, 0, 0],
    0x2062: [0, 0, 0],
    0x2063: [0, 0, 0],
    0x2064: [0, 0, 0],
    0x20D7: [.714, -0.516, .5],
    0x2102: [.702, .019, .722],
    0x210B: [.717, .036, .969, { ic: .272, sk: .333 }],
    0x210C: [.666, .133, .72],
    0x210D: [.683, 0, .778],
    0x210E: [.694, .011, .576, { sk: -0.0278 }],
    0x210F: [.695, .013, .54, { ic: .022 }],
    0x2110: [.717, .017, .809, { ic: .137, sk: .333 }],
    0x2111: [.686, .026, .554],
    0x2112: [.717, .017, .874, { ic: .161, sk: .306 }],
    0x2113: [.705, .02, .417, { sk: .111 }],
    0x2115: [.683, .02, .722],
    0x2118: [.453, .216, .636, { sk: .111 }],
    0x2119: [.683, 0, .611],
    0x211A: [.701, .181, .778],
    0x211B: [.717, .017, .85, { ic: .037, sk: .194 }],
    0x211C: [.686, .026, .828],
    0x211D: [.683, 0, .722],
    0x2124: [.683, 0, .667],
    0x2126: [.704, 0, .722],
    0x2127: [.684, .022, .722],
    0x2128: [.729, .139, .602],
    0x212C: [.708, .028, .908, { ic: .02, sk: .194 }],
    0x212D: [.685, .024, .613],
    0x2130: [.707, .008, .562, { ic: .156, sk: .139 }],
    0x2131: [.735, .036, .895, { ic: .095, sk: .222 }],
    0x2132: [.695, 0, .556],
    0x2133: [.721, .05, 1.08, { ic: .136, sk: .444 }],
    0x2135: [.694, 0, .611],
    0x2136: [.763, .021, .667, { ic: .02 }],
    0x2137: [.764, .043, .444],
    0x2138: [.764, .043, .667],
    0x2141: [.705, .023, .639],
    0x2190: [.511, .011, 1],
    0x2191: [.694, .193, .5],
    0x2192: [.511, .011, 1],
    0x2193: [.694, .194, .5],
    0x2194: [.511, .011, 1],
    0x2195: [.772, .272, .5],
    0x2196: [.72, .195, 1],
    0x2197: [.72, .195, 1],
    0x2198: [.695, .22, 1],
    0x2199: [.695, .22, 1],
    0x219A: [.437, -0.06, 1],
    0x219B: [.437, -0.06, 1],
    0x219E: [.417, -0.083, 1],
    0x21A0: [.417, -0.083, 1],
    0x21A2: [.417, -0.083, 1.111],
    0x21A3: [.417, -0.083, 1.111],
    0x21A6: [.511, .011, 1],
    0x21A9: [.511, .011, 1.126],
    0x21AA: [.511, .011, 1.126],
    0x21AB: [.575, .041, 1],
    0x21AC: [.575, .041, 1],
    0x21AD: [.417, -0.083, 1.389],
    0x21AE: [.437, -0.06, 1],
    0x21B0: [.722, 0, .5],
    0x21B1: [.722, 0, .5],
    0x21B6: [.461, 0, 1],
    0x21B7: [.46, 0, 1],
    0x21BA: [.65, .083, .778],
    0x21BB: [.65, .083, .778],
    0x21BC: [.511, -0.23, 1],
    0x21BD: [.27, .011, 1],
    0x21BE: [.694, .194, .417],
    0x21BF: [.694, .194, .417],
    0x21C0: [.511, -0.23, 1],
    0x21C1: [.27, .011, 1],
    0x21C2: [.694, .194, .417],
    0x21C3: [.694, .194, .417],
    0x21C4: [.667, 0, 1],
    0x21C6: [.667, 0, 1],
    0x21C7: [.583, .083, 1],
    0x21C8: [.694, .193, .833],
    0x21C9: [.583, .083, 1],
    0x21CA: [.694, .194, .833],
    0x21CB: [.514, .014, 1],
    0x21CC: [.671, .011, 1],
    0x21CD: [.534, .035, 1],
    0x21CE: [.534, .037, 1],
    0x21CF: [.534, .035, 1],
    0x21D0: [.525, .024, 1],
    0x21D1: [.694, .194, .611],
    0x21D2: [.525, .024, 1],
    0x21D3: [.694, .194, .611],
    0x21D4: [.526, .025, 1],
    0x21D5: [.772, .272, .611],
    0x21DA: [.611, .111, 1],
    0x21DB: [.611, .111, 1],
    0x21DD: [.417, -0.083, 1],
    0x21E0: [.437, -0.064, 1.334],
    0x21E2: [.437, -0.064, 1.334],
    0x2200: [.694, .022, .556],
    0x2201: [.846, .021, .5],
    0x2202: [.715, .022, .531, { ic: .035, sk: .0833 }],
    0x2203: [.694, 0, .556],
    0x2204: [.716, .215, .556],
    0x2205: [.772, .078, .5],
    0x2206: [.716, 0, .833],
    0x2207: [.683, .033, .833],
    0x2208: [.54, .04, .667],
    0x2209: [.716, .215, .667],
    0x220B: [.54, .04, .667],
    0x220C: [.716, .215, .667],
    0x220D: [.44, 0, .429, { ic: .027 }],
    0x220F: [.75, .25, .944],
    0x2210: [.75, .25, .944],
    0x2211: [.75, .25, 1.056],
    0x2212: [.583, .082, .778],
    0x2213: [.5, .166, .778],
    0x2214: [.766, .093, .778],
    0x2215: [.75, .25, .5],
    0x2216: [.75, .25, .5],
    0x2217: [.465, -0.035, .5],
    0x2218: [.444, -0.055, .5],
    0x2219: [.444, -0.055, .5],
    0x221A: [.8, .2, .833, { ic: .02 }],
    0x221D: [.442, .011, .778],
    0x221E: [.442, .011, 1],
    0x2220: [.694, 0, .722],
    0x2221: [.714, .02, .722],
    0x2222: [.551, .051, .722],
    0x2223: [.75, .249, .278],
    0x2224: [.75, .252, .278, { ic: .019 }],
    0x2225: [.75, .25, .5],
    0x2226: [.75, .25, .5, { ic: .018 }],
    0x2227: [.598, .022, .667],
    0x2228: [.598, .022, .667],
    0x2229: [.598, .022, .667],
    0x222A: [.598, .022, .667],
    0x222B: [.716, .216, .417, { ic: .055 }],
    0x222C: [.805, .306, .819, { ic: .138 }],
    0x222D: [.805, .306, 1.166, { ic: .138 }],
    0x222E: [.805, .306, .472, { ic: .138 }],
    0x2234: [.471, .082, .667],
    0x2235: [.471, .082, .667],
    0x223C: [.367, -0.133, .778],
    0x223D: [.367, -0.133, .778],
    0x2240: [.583, .083, .278],
    0x2241: [.467, -0.032, .778],
    0x2242: [.463, -0.034, .778],
    0x2243: [.464, -0.036, .778],
    0x2244: [.716, .215, .778],
    0x2245: [.589, -0.022, .778],
    0x2247: [.652, .155, .778],
    0x2248: [.483, -0.055, .778],
    0x2249: [.716, .215, .778],
    0x224A: [.579, .039, .778],
    0x224D: [.484, -0.016, .778],
    0x224E: [.492, -0.008, .778],
    0x224F: [.492, -0.133, .778],
    0x2250: [.67, -0.133, .778],
    0x2251: [.609, .108, .778],
    0x2252: [.601, .101, .778],
    0x2253: [.601, .102, .778],
    0x2256: [.367, -0.133, .778],
    0x2257: [.721, -0.133, .778],
    0x225C: [.859, -0.133, .778],
    0x2260: [.716, .215, .778],
    0x2261: [.464, -0.036, .778],
    0x2262: [.716, .215, .778],
    0x2264: [.636, .138, .778],
    0x2265: [.636, .138, .778],
    0x2266: [.753, .175, .778],
    0x2267: [.753, .175, .778],
    0x2268: [.752, .286, .778],
    0x2269: [.752, .286, .778],
    0x226A: [.568, .067, 1],
    0x226B: [.567, .067, 1],
    0x226C: [.75, .25, .5],
    0x226D: [.716, .215, .778],
    0x226E: [.708, .209, .778],
    0x226F: [.708, .209, .778],
    0x2270: [.801, .303, .778],
    0x2271: [.801, .303, .778],
    0x2272: [.732, .228, .778],
    0x2273: [.732, .228, .778],
    0x2274: [.732, .228, .778],
    0x2275: [.732, .228, .778],
    0x2276: [.681, .253, .778],
    0x2277: [.681, .253, .778],
    0x2278: [.716, .253, .778],
    0x2279: [.716, .253, .778],
    0x227A: [.539, .041, .778],
    0x227B: [.539, .041, .778],
    0x227C: [.58, .153, .778],
    0x227D: [.58, .154, .778],
    0x227E: [.732, .228, .778],
    0x227F: [.732, .228, .778],
    0x2280: [.705, .208, .778],
    0x2281: [.705, .208, .778],
    0x2282: [.54, .04, .778],
    0x2283: [.54, .04, .778],
    0x2284: [.716, .215, .778],
    0x2285: [.716, .215, .778],
    0x2286: [.636, .138, .778],
    0x2287: [.636, .138, .778],
    0x2288: [.801, .303, .778],
    0x2289: [.801, .303, .778],
    0x228A: [.635, .241, .778],
    0x228B: [.635, .241, .778],
    0x228E: [.598, .022, .667],
    0x228F: [.539, .041, .778],
    0x2290: [.539, .041, .778],
    0x2291: [.636, .138, .778],
    0x2292: [.636, .138, .778],
    0x2293: [.598, 0, .667],
    0x2294: [.598, 0, .667],
    0x2295: [.583, .083, .778],
    0x2296: [.583, .083, .778],
    0x2297: [.583, .083, .778],
    0x2298: [.583, .083, .778],
    0x2299: [.583, .083, .778],
    0x229A: [.582, .082, .778],
    0x229B: [.582, .082, .778],
    0x229D: [.582, .082, .778],
    0x229E: [.689, 0, .778],
    0x229F: [.689, 0, .778],
    0x22A0: [.689, 0, .778],
    0x22A1: [.689, 0, .778],
    0x22A2: [.694, 0, .611],
    0x22A3: [.694, 0, .611],
    0x22A4: [.668, 0, .778],
    0x22A5: [.668, 0, .778],
    0x22A8: [.75, .249, .867],
    0x22A9: [.694, 0, .722],
    0x22AA: [.694, 0, .889],
    0x22AC: [.695, 0, .611],
    0x22AD: [.695, 0, .611],
    0x22AE: [.695, 0, .722],
    0x22AF: [.695, 0, .722],
    0x22B2: [.539, .041, .778],
    0x22B3: [.539, .041, .778],
    0x22B4: [.636, .138, .778],
    0x22B5: [.636, .138, .778],
    0x22B8: [.408, -0.092, 1.111],
    0x22BA: [.431, .212, .556],
    0x22BB: [.716, 0, .611],
    0x22BC: [.716, 0, .611],
    0x22C0: [.75, .249, .833],
    0x22C1: [.75, .249, .833],
    0x22C2: [.75, .249, .833],
    0x22C3: [.75, .249, .833],
    0x22C4: [.488, -0.012, .5],
    0x22C5: [.31, -0.19, .278],
    0x22C6: [.486, -0.016, .5],
    0x22C7: [.545, .044, .778],
    0x22C8: [.505, .005, .9],
    0x22C9: [.492, -0.008, .778],
    0x22CA: [.492, -0.008, .778],
    0x22CB: [.694, .022, .778],
    0x22CC: [.694, .022, .778],
    0x22CD: [.464, -0.036, .778],
    0x22CE: [.578, .021, .76],
    0x22CF: [.578, .022, .76],
    0x22D0: [.54, .04, .778],
    0x22D1: [.54, .04, .778],
    0x22D2: [.598, .022, .667],
    0x22D3: [.598, .022, .667],
    0x22D4: [.736, .022, .667],
    0x22D6: [.541, .041, .778],
    0x22D7: [.541, .041, .778],
    0x22D8: [.568, .067, 1.333],
    0x22D9: [.568, .067, 1.333],
    0x22DA: [.886, .386, .778],
    0x22DB: [.886, .386, .778],
    0x22DE: [.734, 0, .778],
    0x22DF: [.734, 0, .778],
    0x22E0: [.801, .303, .778],
    0x22E1: [.801, .303, .778],
    0x22E2: [.716, .215, .778],
    0x22E3: [.716, .215, .778],
    0x22E6: [.73, .359, .778],
    0x22E7: [.73, .359, .778],
    0x22E8: [.73, .359, .778],
    0x22E9: [.73, .359, .778],
    0x22EA: [.706, .208, .778],
    0x22EB: [.706, .208, .778],
    0x22EC: [.802, .303, .778],
    0x22ED: [.801, .303, .778],
    0x22EE: [1.3, .03, .278],
    0x22EF: [.31, -0.19, 1.172],
    0x22F1: [1.52, -0.1, 1.282],
    0x2305: [.716, 0, .611],
    0x2306: [.813, .097, .611],
    0x2308: [.75, .25, .444],
    0x2309: [.75, .25, .444],
    0x230A: [.75, .25, .444],
    0x230B: [.75, .25, .444],
    0x231C: [.694, -0.306, .5],
    0x231D: [.694, -0.306, .5],
    0x231E: [.366, .022, .5],
    0x231F: [.366, .022, .5],
    0x2322: [.388, -0.122, 1],
    0x2323: [.378, -0.134, 1],
    0x2329: [.75, .25, .389],
    0x232A: [.75, .25, .389],
    0x23B0: [.744, .244, .412],
    0x23B1: [.744, .244, .412],
    0x23D0: [.602, 0, .667],
    0x24C8: [.709, .175, .902],
    0x250C: [.694, -0.306, .5],
    0x2510: [.694, -0.306, .5],
    0x2514: [.366, .022, .5],
    0x2518: [.366, .022, .5],
    0x2571: [.694, .195, .889],
    0x2572: [.694, .195, .889],
    0x25A0: [.689, 0, .778],
    0x25A1: [.689, 0, .778],
    0x25AA: [.689, 0, .778],
    0x25B2: [.575, .02, .722],
    0x25B3: [.716, 0, .889],
    0x25B4: [.575, .02, .722],
    0x25B5: [.716, 0, .889],
    0x25B6: [.539, .041, .778],
    0x25B8: [.539, .041, .778],
    0x25B9: [.505, .005, .5],
    0x25BC: [.576, .019, .722],
    0x25BD: [.5, .215, .889],
    0x25BE: [.576, .019, .722],
    0x25BF: [.5, .215, .889],
    0x25C0: [.539, .041, .778],
    0x25C2: [.539, .041, .778],
    0x25C3: [.505, .005, .5],
    0x25CA: [.716, .132, .667],
    0x25EF: [.715, .215, 1],
    0x25FB: [.689, 0, .778],
    0x25FC: [.689, 0, .778],
    0x2605: [.694, .111, .944],
    0x2660: [.727, .13, .778],
    0x2661: [.716, .033, .778],
    0x2662: [.727, .162, .778],
    0x2663: [.726, .13, .778],
    0x266D: [.75, .022, .389],
    0x266E: [.734, .223, .389],
    0x266F: [.723, .223, .389],
    0x2713: [.706, .034, .833],
    0x2720: [.716, .022, .833],
    0x2758: [.75, .249, .278],
    0x27E8: [.75, .25, .389],
    0x27E9: [.75, .25, .389],
    0x27EE: [.744, .244, .412],
    0x27EF: [.744, .244, .412],
    0x27F5: [.511, .011, 1.609],
    0x27F6: [.511, .011, 1.638],
    0x27F7: [.511, .011, 1.859],
    0x27F8: [.525, .024, 1.609],
    0x27F9: [.525, .024, 1.638],
    0x27FA: [.525, .024, 1.858],
    0x27FC: [.511, .011, 1.638],
    0x29EB: [.716, .132, .667],
    0x29F8: [.716, .215, .778],
    0x2A00: [.75, .25, 1.111],
    0x2A01: [.75, .25, 1.111],
    0x2A02: [.75, .25, 1.111],
    0x2A04: [.75, .249, .833],
    0x2A06: [.75, .249, .833],
    0x2A0C: [.805, .306, 1.638, { ic: .138 }],
    0x2A2F: [.491, -0.009, .778],
    0x2A3F: [.683, 0, .75],
    0x2A5E: [.813, .097, .611],
    0x2A7D: [.636, .138, .778],
    0x2A7E: [.636, .138, .778],
    0x2A85: [.762, .29, .778],
    0x2A86: [.762, .29, .778],
    0x2A87: [.635, .241, .778],
    0x2A88: [.635, .241, .778],
    0x2A89: [.761, .387, .778],
    0x2A8A: [.761, .387, .778],
    0x2A8B: [1.003, .463, .778],
    0x2A8C: [1.003, .463, .778],
    0x2A95: [.636, .138, .778],
    0x2A96: [.636, .138, .778],
    0x2AAF: [.636, .138, .778],
    0x2AB0: [.636, .138, .778],
    0x2AB5: [.752, .286, .778],
    0x2AB6: [.752, .286, .778],
    0x2AB7: [.761, .294, .778],
    0x2AB8: [.761, .294, .778],
    0x2AB9: [.761, .337, .778],
    0x2ABA: [.761, .337, .778],
    0x2AC5: [.753, .215, .778],
    0x2AC6: [.753, .215, .778],
    0x2ACB: [.783, .385, .778],
    0x2ACC: [.783, .385, .778],
    0x3008: [.75, .25, .389],
    0x3009: [.75, .25, .389],
    0xE006: [.43, .023, .222, { ic: .018 }],
    0xE007: [.431, .024, .389, { ic: .018 }],
    0xE008: [.605, .085, .778],
    0xE009: [.434, .006, .667, { ic: .067 }],
    0xE00C: [.752, .284, .778],
    0xE00D: [.752, .284, .778],
    0xE00E: [.919, .421, .778],
    0xE00F: [.801, .303, .778],
    0xE010: [.801, .303, .778],
    0xE011: [.919, .421, .778],
    0xE016: [.828, .33, .778],
    0xE017: [.752, .332, .778],
    0xE018: [.828, .33, .778],
    0xE019: [.752, .333, .778],
    0xE01A: [.634, .255, .778],
    0xE01B: [.634, .254, .778],
    0x1D400: [.698, 0, .869],
    0x1D401: [.686, 0, .818],
    0x1D402: [.697, .011, .831],
    0x1D403: [.686, 0, .882],
    0x1D404: [.68, 0, .756],
    0x1D405: [.68, 0, .724],
    0x1D406: [.697, .01, .904],
    0x1D407: [.686, 0, .9],
    0x1D408: [.686, 0, .436],
    0x1D409: [.686, .011, .594],
    0x1D40A: [.686, 0, .901],
    0x1D40B: [.686, 0, .692],
    0x1D40C: [.686, 0, 1.092],
    0x1D40D: [.686, 0, .9],
    0x1D40E: [.696, .01, .864],
    0x1D40F: [.686, 0, .786],
    0x1D410: [.696, .193, .864],
    0x1D411: [.686, .011, .862],
    0x1D412: [.697, .011, .639],
    0x1D413: [.675, 0, .8],
    0x1D414: [.686, .011, .885],
    0x1D415: [.686, .007, .869],
    0x1D416: [.686, .007, 1.189],
    0x1D417: [.686, 0, .869],
    0x1D418: [.686, 0, .869],
    0x1D419: [.686, 0, .703],
    0x1D41A: [.453, .006, .559],
    0x1D41B: [.694, .006, .639],
    0x1D41C: [.453, .006, .511],
    0x1D41D: [.694, .006, .639],
    0x1D41E: [.452, .006, .527],
    0x1D41F: [.7, 0, .351, { ic: .101 }],
    0x1D420: [.455, .201, .575],
    0x1D421: [.694, 0, .639],
    0x1D422: [.695, 0, .319],
    0x1D423: [.695, .2, .351],
    0x1D424: [.694, 0, .607],
    0x1D425: [.694, 0, .319],
    0x1D426: [.45, 0, .958],
    0x1D427: [.45, 0, .639],
    0x1D428: [.452, .005, .575],
    0x1D429: [.45, .194, .639],
    0x1D42A: [.45, .194, .607],
    0x1D42B: [.45, 0, .474],
    0x1D42C: [.453, .006, .454],
    0x1D42D: [.635, .005, .447],
    0x1D42E: [.45, .006, .639],
    0x1D42F: [.444, 0, .607],
    0x1D430: [.444, 0, .831],
    0x1D431: [.444, 0, .607],
    0x1D432: [.444, .2, .607],
    0x1D433: [.444, 0, .511],
    0x1D434: [.716, 0, .75, { sk: .139 }],
    0x1D435: [.683, 0, .759, { sk: .0833 }],
    0x1D436: [.705, .022, .715, { ic: .045, sk: .0833 }],
    0x1D437: [.683, 0, .828, { sk: .0556 }],
    0x1D438: [.68, 0, .738, { ic: .026, sk: .0833 }],
    0x1D439: [.68, 0, .643, { ic: .106, sk: .0833 }],
    0x1D43A: [.705, .022, .786, { sk: .0833 }],
    0x1D43B: [.683, 0, .831, { ic: .057, sk: .0556 }],
    0x1D43C: [.683, 0, .44, { ic: .064, sk: .111 }],
    0x1D43D: [.683, .022, .555, { ic: .078, sk: .167 }],
    0x1D43E: [.683, 0, .849, { ic: .04, sk: .0556 }],
    0x1D43F: [.683, 0, .681, { sk: .0278 }],
    0x1D440: [.683, 0, .97, { ic: .081, sk: .0833 }],
    0x1D441: [.683, 0, .803, { ic: .085, sk: .0833 }],
    0x1D442: [.704, .022, .763, { sk: .0833 }],
    0x1D443: [.683, 0, .642, { ic: .109, sk: .0833 }],
    0x1D444: [.704, .194, .791, { sk: .0833 }],
    0x1D445: [.683, .021, .759, { sk: .0833 }],
    0x1D446: [.705, .022, .613, { ic: .032, sk: .0833 }],
    0x1D447: [.677, 0, .584, { ic: .12, sk: .0833 }],
    0x1D448: [.683, .022, .683, { ic: .084, sk: .0278 }],
    0x1D449: [.683, .022, .583, { ic: .186 }],
    0x1D44A: [.683, .022, .944, { ic: .104 }],
    0x1D44B: [.683, 0, .828, { ic: .024, sk: .0833 }],
    0x1D44C: [.683, 0, .581, { ic: .182 }],
    0x1D44D: [.683, 0, .683, { ic: .04, sk: .0833 }],
    0x1D44E: [.441, .01, .529],
    0x1D44F: [.694, .011, .429],
    0x1D450: [.442, .011, .433, { sk: .0556 }],
    0x1D451: [.694, .01, .52, { sk: .167 }],
    0x1D452: [.442, .011, .466, { sk: .0556 }],
    0x1D453: [.705, .205, .49, { ic: .06, sk: .167 }],
    0x1D454: [.442, .205, .477, { sk: .0278 }],
    0x1D456: [.661, .011, .345],
    0x1D457: [.661, .204, .412],
    0x1D458: [.694, .011, .521],
    0x1D459: [.694, .011, .298, { sk: .0833 }],
    0x1D45A: [.442, .011, .878],
    0x1D45B: [.442, .011, .6],
    0x1D45C: [.441, .011, .485, { sk: .0556 }],
    0x1D45D: [.442, .194, .503, { sk: .0833 }],
    0x1D45E: [.442, .194, .446, { ic: .014, sk: .0833 }],
    0x1D45F: [.442, .011, .451, { sk: .0556 }],
    0x1D460: [.442, .01, .469, { sk: .0556 }],
    0x1D461: [.626, .011, .361, { sk: .0833 }],
    0x1D462: [.442, .011, .572, { sk: .0278 }],
    0x1D463: [.443, .011, .485, { sk: .0278 }],
    0x1D464: [.443, .011, .716, { sk: .0833 }],
    0x1D465: [.442, .011, .572, { sk: .0278 }],
    0x1D466: [.442, .205, .49, { sk: .0556 }],
    0x1D467: [.442, .011, .465, { sk: .0556 }],
    0x1D468: [.711, 0, .869, { sk: .16 }],
    0x1D469: [.686, 0, .866, { sk: .0958 }],
    0x1D46A: [.703, .017, .817, { ic: .038, sk: .0958 }],
    0x1D46B: [.686, 0, .938, { sk: .0639 }],
    0x1D46C: [.68, 0, .81, { ic: .015, sk: .0958 }],
    0x1D46D: [.68, 0, .689, { ic: .12, sk: .0958 }],
    0x1D46E: [.703, .016, .887, { sk: .0958 }],
    0x1D46F: [.686, 0, .982, { ic: .045, sk: .0639 }],
    0x1D470: [.686, 0, .511, { ic: .062, sk: .128 }],
    0x1D471: [.686, .017, .631, { ic: .063, sk: .192 }],
    0x1D472: [.686, 0, .971, { ic: .032, sk: .0639 }],
    0x1D473: [.686, 0, .756, { sk: .0319 }],
    0x1D474: [.686, 0, 1.142, { ic: .077, sk: .0958 }],
    0x1D475: [.686, 0, .95, { ic: .077, sk: .0958 }],
    0x1D476: [.703, .017, .837, { sk: .0958 }],
    0x1D477: [.686, 0, .723, { ic: .124, sk: .0958 }],
    0x1D478: [.703, .194, .869, { sk: .0958 }],
    0x1D479: [.686, .017, .872, { sk: .0958 }],
    0x1D47A: [.703, .017, .693, { ic: .021, sk: .0958 }],
    0x1D47B: [.675, 0, .637, { ic: .135, sk: .0958 }],
    0x1D47C: [.686, .016, .8, { ic: .077, sk: .0319 }],
    0x1D47D: [.686, .016, .678, { ic: .208 }],
    0x1D47E: [.686, .017, 1.093, { ic: .114 }],
    0x1D47F: [.686, 0, .947, { sk: .0958 }],
    0x1D480: [.686, 0, .675, { ic: .201 }],
    0x1D481: [.686, 0, .773, { ic: .032, sk: .0958 }],
    0x1D482: [.452, .008, .633],
    0x1D483: [.694, .008, .521],
    0x1D484: [.451, .008, .513, { sk: .0639 }],
    0x1D485: [.694, .008, .61, { sk: .192 }],
    0x1D486: [.452, .008, .554, { sk: .0639 }],
    0x1D487: [.701, .201, .568, { ic: .056, sk: .192 }],
    0x1D488: [.452, .202, .545, { sk: .0319 }],
    0x1D489: [.694, .008, .668, { sk: -0.0319 }],
    0x1D48A: [.694, .008, .405],
    0x1D48B: [.694, .202, .471],
    0x1D48C: [.694, .008, .604],
    0x1D48D: [.694, .008, .348, { sk: .0958 }],
    0x1D48E: [.452, .008, 1.032],
    0x1D48F: [.452, .008, .713],
    0x1D490: [.452, .008, .585, { sk: .0639 }],
    0x1D491: [.452, .194, .601, { sk: .0958 }],
    0x1D492: [.452, .194, .542, { sk: .0958 }],
    0x1D493: [.452, .008, .529, { sk: .0639 }],
    0x1D494: [.451, .008, .531, { sk: .0639 }],
    0x1D495: [.643, .007, .415, { sk: .0958 }],
    0x1D496: [.452, .008, .681, { sk: .0319 }],
    0x1D497: [.453, .008, .567, { sk: .0319 }],
    0x1D498: [.453, .008, .831, { sk: .0958 }],
    0x1D499: [.452, .008, .659, { sk: .0319 }],
    0x1D49A: [.452, .202, .59, { sk: .0639 }],
    0x1D49B: [.452, .008, .555, { sk: .0639 }],
    0x1D49C: [.717, .008, .803, { ic: .213, sk: .389 }],
    0x1D49E: [.728, .026, .666, { ic: .153, sk: .278 }],
    0x1D49F: [.708, .031, .774, { ic: .081, sk: .111 }],
    0x1D4A2: [.717, .037, .61, { ic: .128, sk: .25 }],
    0x1D4A5: [.717, .314, 1.052, { ic: .081, sk: .417 }],
    0x1D4A6: [.717, .037, .914, { ic: .29, sk: .361 }],
    0x1D4A9: [.726, .036, .902, { ic: .306, sk: .389 }],
    0x1D4AA: [.707, .008, .738, { ic: .067, sk: .167 }],
    0x1D4AB: [.716, .037, 1.013, { ic: .018, sk: .222 }],
    0x1D4AC: [.717, .017, .883, { sk: .278 }],
    0x1D4AE: [.708, .036, .868, { ic: .148, sk: .333 }],
    0x1D4AF: [.735, .037, .747, { ic: .249, sk: .222 }],
    0x1D4B0: [.717, .017, .8, { ic: .16, sk: .25 }],
    0x1D4B1: [.717, .017, .622, { ic: .228, sk: .222 }],
    0x1D4B2: [.717, .017, .805, { ic: .221, sk: .25 }],
    0x1D4B3: [.717, .017, .944, { ic: .187, sk: .278 }],
    0x1D4B4: [.716, .017, .71, { ic: .249, sk: .194 }],
    0x1D4B5: [.717, .016, .821, { ic: .211, sk: .306 }],
    0x1D504: [.696, .026, .718],
    0x1D505: [.691, .027, .884],
    0x1D507: [.685, .027, .832],
    0x1D508: [.685, .024, .663],
    0x1D509: [.686, .153, .611],
    0x1D50A: [.69, .026, .785],
    0x1D50D: [.686, .139, .552],
    0x1D50E: [.68, .027, .668, { ic: .014 }],
    0x1D50F: [.686, .026, .666],
    0x1D510: [.692, .027, 1.05],
    0x1D511: [.686, .025, .832],
    0x1D512: [.729, .027, .827],
    0x1D513: [.692, .218, .828],
    0x1D514: [.729, .069, .827],
    0x1D516: [.692, .027, .829],
    0x1D517: [.701, .027, .669],
    0x1D518: [.697, .027, .646, { ic: .019 }],
    0x1D519: [.686, .026, .831],
    0x1D51A: [.686, .027, 1.046],
    0x1D51B: [.688, .027, .719],
    0x1D51C: [.686, .218, .833],
    0x1D51E: [.47, .035, .5],
    0x1D51F: [.685, .031, .513],
    0x1D520: [.466, .029, .389],
    0x1D521: [.609, .033, .499],
    0x1D522: [.467, .03, .401],
    0x1D523: [.681, .221, .326],
    0x1D524: [.47, .209, .504],
    0x1D525: [.688, .205, .521],
    0x1D526: [.673, .02, .279],
    0x1D527: [.672, .208, .281],
    0x1D528: [.689, .025, .389],
    0x1D529: [.685, .02, .28],
    0x1D52A: [.475, .026, .767],
    0x1D52B: [.475, .022, .527],
    0x1D52C: [.48, .028, .489],
    0x1D52D: [.541, .212, .5],
    0x1D52E: [.479, .219, .489],
    0x1D52F: [.474, .021, .389],
    0x1D530: [.478, .029, .443],
    0x1D531: [.64, .02, .333, { ic: .015 }],
    0x1D532: [.474, .023, .517],
    0x1D533: [.53, .028, .512],
    0x1D534: [.532, .028, .774],
    0x1D535: [.472, .188, .389],
    0x1D536: [.528, .218, .499],
    0x1D537: [.471, .214, .391],
    0x1D538: [.701, 0, .722],
    0x1D539: [.683, 0, .667],
    0x1D53B: [.683, 0, .722],
    0x1D53C: [.683, 0, .667],
    0x1D53D: [.683, 0, .611],
    0x1D53E: [.702, .019, .778],
    0x1D540: [.683, 0, .389],
    0x1D541: [.683, .077, .5],
    0x1D542: [.683, 0, .778],
    0x1D543: [.683, 0, .667],
    0x1D544: [.683, 0, .944],
    0x1D546: [.701, .019, .778],
    0x1D54A: [.702, .012, .556],
    0x1D54B: [.683, 0, .667],
    0x1D54C: [.683, .019, .722],
    0x1D54D: [.683, .02, .722],
    0x1D54E: [.683, .019, 1],
    0x1D54F: [.683, 0, .722],
    0x1D550: [.683, 0, .722],
    0x1D56C: [.686, .031, .847],
    0x1D56D: [.684, .031, 1.044],
    0x1D56E: [.676, .032, .723],
    0x1D56F: [.683, .029, .982],
    0x1D570: [.686, .029, .783],
    0x1D571: [.684, .146, .722],
    0x1D572: [.687, .029, .927],
    0x1D573: [.683, .126, .851],
    0x1D574: [.681, .025, .655],
    0x1D575: [.68, .141, .652],
    0x1D576: [.681, .026, .789, { ic: .017 }],
    0x1D577: [.683, .028, .786],
    0x1D578: [.683, .032, 1.239],
    0x1D579: [.679, .03, .983],
    0x1D57A: [.726, .03, .976],
    0x1D57B: [.688, .223, .977],
    0x1D57C: [.726, .083, .976],
    0x1D57D: [.688, .028, .978],
    0x1D57E: [.685, .031, .978],
    0x1D57F: [.686, .03, .79, { ic: .012 }],
    0x1D580: [.688, .039, .851, { ic: .02 }],
    0x1D581: [.685, .029, .982],
    0x1D582: [.683, .03, 1.235],
    0x1D583: [.681, .035, .849],
    0x1D584: [.688, .214, .984],
    0x1D585: [.677, .148, .711],
    0x1D586: [.472, .032, .603],
    0x1D587: [.69, .032, .59],
    0x1D588: [.473, .026, .464],
    0x1D589: [.632, .028, .589],
    0x1D58A: [.471, .027, .472],
    0x1D58B: [.687, .222, .388],
    0x1D58C: [.472, .208, .595],
    0x1D58D: [.687, .207, .615],
    0x1D58E: [.686, .025, .331],
    0x1D58F: [.682, .203, .332],
    0x1D590: [.682, .025, .464],
    0x1D591: [.681, .024, .337],
    0x1D592: [.476, .031, .921],
    0x1D593: [.473, .028, .654],
    0x1D594: [.482, .034, .609],
    0x1D595: [.557, .207, .604],
    0x1D596: [.485, .211, .596],
    0x1D597: [.472, .026, .46],
    0x1D598: [.479, .034, .523],
    0x1D599: [.648, .027, .393, { ic: .014 }],
    0x1D59A: [.472, .032, .589, { ic: .014 }],
    0x1D59B: [.546, .027, .604],
    0x1D59C: [.549, .032, .918],
    0x1D59D: [.471, .188, .459],
    0x1D59E: [.557, .221, .589],
    0x1D59F: [.471, .214, .461],
    0x1D5A0: [.694, 0, .667],
    0x1D5A1: [.694, 0, .667],
    0x1D5A2: [.705, .011, .639],
    0x1D5A3: [.694, 0, .722],
    0x1D5A4: [.691, 0, .597],
    0x1D5A5: [.691, 0, .569],
    0x1D5A6: [.704, .011, .667],
    0x1D5A7: [.694, 0, .708],
    0x1D5A8: [.694, 0, .278],
    0x1D5A9: [.694, .022, .472],
    0x1D5AA: [.694, 0, .694],
    0x1D5AB: [.694, 0, .542],
    0x1D5AC: [.694, 0, .875],
    0x1D5AD: [.694, 0, .708],
    0x1D5AE: [.715, .022, .736],
    0x1D5AF: [.694, 0, .639],
    0x1D5B0: [.715, .125, .736],
    0x1D5B1: [.694, 0, .646],
    0x1D5B2: [.716, .022, .556],
    0x1D5B3: [.688, 0, .681],
    0x1D5B4: [.694, .022, .688],
    0x1D5B5: [.694, 0, .667],
    0x1D5B6: [.694, 0, .944],
    0x1D5B7: [.694, 0, .667],
    0x1D5B8: [.694, 0, .667],
    0x1D5B9: [.694, 0, .611],
    0x1D5BA: [.46, .01, .481],
    0x1D5BB: [.694, .011, .517],
    0x1D5BC: [.46, .01, .444],
    0x1D5BD: [.694, .01, .517],
    0x1D5BE: [.461, .01, .444],
    0x1D5BF: [.705, 0, .306, { ic: .041 }],
    0x1D5C0: [.455, .206, .5],
    0x1D5C1: [.694, 0, .517],
    0x1D5C2: [.68, 0, .239],
    0x1D5C3: [.68, .205, .267],
    0x1D5C4: [.694, 0, .489],
    0x1D5C5: [.694, 0, .239],
    0x1D5C6: [.455, 0, .794],
    0x1D5C7: [.455, 0, .517],
    0x1D5C8: [.46, .01, .5],
    0x1D5C9: [.455, .194, .517],
    0x1D5CA: [.455, .194, .517],
    0x1D5CB: [.455, 0, .342],
    0x1D5CC: [.46, .01, .383],
    0x1D5CD: [.571, .01, .361],
    0x1D5CE: [.444, .01, .517],
    0x1D5CF: [.444, 0, .461],
    0x1D5D0: [.444, 0, .683],
    0x1D5D1: [.444, 0, .461],
    0x1D5D2: [.444, .204, .461],
    0x1D5D3: [.444, 0, .435],
    0x1D5D4: [.694, 0, .733],
    0x1D5D5: [.694, 0, .733],
    0x1D5D6: [.704, .011, .703],
    0x1D5D7: [.694, 0, .794],
    0x1D5D8: [.691, 0, .642],
    0x1D5D9: [.691, 0, .611],
    0x1D5DA: [.705, .011, .733],
    0x1D5DB: [.694, 0, .794],
    0x1D5DC: [.694, 0, .331],
    0x1D5DD: [.694, .022, .519],
    0x1D5DE: [.694, 0, .764],
    0x1D5DF: [.694, 0, .581],
    0x1D5E0: [.694, 0, .978],
    0x1D5E1: [.694, 0, .794],
    0x1D5E2: [.716, .022, .794],
    0x1D5E3: [.694, 0, .703],
    0x1D5E4: [.716, .106, .794],
    0x1D5E5: [.694, 0, .703],
    0x1D5E6: [.716, .022, .611],
    0x1D5E7: [.688, 0, .733],
    0x1D5E8: [.694, .022, .764],
    0x1D5E9: [.694, 0, .733],
    0x1D5EA: [.694, 0, 1.039],
    0x1D5EB: [.694, 0, .733],
    0x1D5EC: [.694, 0, .733],
    0x1D5ED: [.694, 0, .672],
    0x1D5EE: [.475, .011, .525],
    0x1D5EF: [.694, .01, .561],
    0x1D5F0: [.475, .011, .489],
    0x1D5F1: [.694, .011, .561],
    0x1D5F2: [.474, .01, .511],
    0x1D5F3: [.705, 0, .336, { ic: .045 }],
    0x1D5F4: [.469, .206, .55],
    0x1D5F5: [.694, 0, .561],
    0x1D5F6: [.695, 0, .256],
    0x1D5F7: [.695, .205, .286],
    0x1D5F8: [.694, 0, .531],
    0x1D5F9: [.694, 0, .256],
    0x1D5FA: [.469, 0, .867],
    0x1D5FB: [.468, 0, .561],
    0x1D5FC: [.474, .011, .55],
    0x1D5FD: [.469, .194, .561],
    0x1D5FE: [.469, .194, .561],
    0x1D5FF: [.469, 0, .372],
    0x1D600: [.474, .01, .422],
    0x1D601: [.589, .01, .404],
    0x1D602: [.458, .011, .561],
    0x1D603: [.458, 0, .5],
    0x1D604: [.458, 0, .744],
    0x1D605: [.458, 0, .5],
    0x1D606: [.458, .205, .5],
    0x1D607: [.458, 0, .476],
    0x1D608: [.694, 0, .667],
    0x1D609: [.694, 0, .667, { ic: .029 }],
    0x1D60A: [.705, .01, .639, { ic: .08 }],
    0x1D60B: [.694, 0, .722, { ic: .025 }],
    0x1D60C: [.691, 0, .597, { ic: .091 }],
    0x1D60D: [.691, 0, .569, { ic: .104 }],
    0x1D60E: [.705, .011, .667, { ic: .063 }],
    0x1D60F: [.694, 0, .708, { ic: .06 }],
    0x1D610: [.694, 0, .278, { ic: .06 }],
    0x1D611: [.694, .022, .472, { ic: .063 }],
    0x1D612: [.694, 0, .694, { ic: .091 }],
    0x1D613: [.694, 0, .542],
    0x1D614: [.694, 0, .875, { ic: .054 }],
    0x1D615: [.694, 0, .708, { ic: .058 }],
    0x1D616: [.716, .022, .736, { ic: .027 }],
    0x1D617: [.694, 0, .639, { ic: .051 }],
    0x1D618: [.716, .125, .736, { ic: .027 }],
    0x1D619: [.694, 0, .646, { ic: .052 }],
    0x1D61A: [.716, .022, .556, { ic: .053 }],
    0x1D61B: [.688, 0, .681, { ic: .109 }],
    0x1D61C: [.694, .022, .688, { ic: .059 }],
    0x1D61D: [.694, 0, .667, { ic: .132 }],
    0x1D61E: [.694, 0, .944, { ic: .132 }],
    0x1D61F: [.694, 0, .667, { ic: .091 }],
    0x1D620: [.694, 0, .667, { ic: .143 }],
    0x1D621: [.694, 0, .611, { ic: .091 }],
    0x1D622: [.461, .01, .481],
    0x1D623: [.694, .011, .517, { ic: .022 }],
    0x1D624: [.46, .011, .444, { ic: .055 }],
    0x1D625: [.694, .01, .517, { ic: .071 }],
    0x1D626: [.46, .011, .444, { ic: .028 }],
    0x1D627: [.705, 0, .306, { ic: .188 }],
    0x1D628: [.455, .206, .5, { ic: .068 }],
    0x1D629: [.694, 0, .517],
    0x1D62A: [.68, 0, .239, { ic: .076 }],
    0x1D62B: [.68, .204, .267, { ic: .069 }],
    0x1D62C: [.694, 0, .489, { ic: .054 }],
    0x1D62D: [.694, 0, .239, { ic: .072 }],
    0x1D62E: [.455, 0, .794],
    0x1D62F: [.454, 0, .517],
    0x1D630: [.461, .011, .5, { ic: .023 }],
    0x1D631: [.455, .194, .517, { ic: .021 }],
    0x1D632: [.455, .194, .517, { ic: .021 }],
    0x1D633: [.455, 0, .342, { ic: .082 }],
    0x1D634: [.461, .011, .383, { ic: .053 }],
    0x1D635: [.571, .011, .361, { ic: .049 }],
    0x1D636: [.444, .01, .517, { ic: .02 }],
    0x1D637: [.444, 0, .461, { ic: .079 }],
    0x1D638: [.444, 0, .683, { ic: .079 }],
    0x1D639: [.444, 0, .461, { ic: .076 }],
    0x1D63A: [.444, .205, .461, { ic: .079 }],
    0x1D63B: [.444, 0, .435, { ic: .059 }],
    0x1D670: [.623, 0, .525],
    0x1D671: [.611, 0, .525],
    0x1D672: [.622, .011, .525],
    0x1D673: [.611, 0, .525],
    0x1D674: [.611, 0, .525],
    0x1D675: [.611, 0, .525],
    0x1D676: [.622, .011, .525],
    0x1D677: [.611, 0, .525],
    0x1D678: [.611, 0, .525],
    0x1D679: [.611, .011, .525],
    0x1D67A: [.611, 0, .525],
    0x1D67B: [.611, 0, .525],
    0x1D67C: [.611, 0, .525],
    0x1D67D: [.611, 0, .525],
    0x1D67E: [.621, .01, .525],
    0x1D67F: [.611, 0, .525],
    0x1D680: [.621, .138, .525],
    0x1D681: [.611, .011, .525],
    0x1D682: [.622, .011, .525],
    0x1D683: [.611, 0, .525],
    0x1D684: [.611, .011, .525],
    0x1D685: [.611, .007, .525],
    0x1D686: [.611, .007, .525],
    0x1D687: [.611, 0, .525],
    0x1D688: [.611, 0, .525],
    0x1D689: [.611, 0, .525],
    0x1D68A: [.439, .006, .525],
    0x1D68B: [.611, .006, .525],
    0x1D68C: [.44, .006, .525],
    0x1D68D: [.611, .006, .525],
    0x1D68E: [.44, .006, .525],
    0x1D68F: [.617, 0, .525],
    0x1D690: [.442, .229, .525],
    0x1D691: [.611, 0, .525],
    0x1D692: [.612, 0, .525],
    0x1D693: [.612, .228, .525],
    0x1D694: [.611, 0, .525],
    0x1D695: [.611, 0, .525],
    0x1D696: [.436, 0, .525, { ic: .011 }],
    0x1D697: [.436, 0, .525],
    0x1D698: [.44, .006, .525],
    0x1D699: [.437, .221, .525],
    0x1D69A: [.437, .221, .525, { ic: .02 }],
    0x1D69B: [.437, 0, .525],
    0x1D69C: [.44, .006, .525],
    0x1D69D: [.554, .006, .525],
    0x1D69E: [.431, .005, .525],
    0x1D69F: [.431, 0, .525],
    0x1D6A0: [.431, 0, .525],
    0x1D6A1: [.431, 0, .525],
    0x1D6A2: [.431, .228, .525],
    0x1D6A3: [.431, 0, .525],
    0x1D6A8: [.698, 0, .869],
    0x1D6A9: [.686, 0, .818],
    0x1D6AA: [.68, 0, .692],
    0x1D6AB: [.698, 0, .958],
    0x1D6AC: [.68, 0, .756],
    0x1D6AD: [.686, 0, .703],
    0x1D6AE: [.686, 0, .9],
    0x1D6AF: [.696, .01, .894],
    0x1D6B0: [.686, 0, .436],
    0x1D6B1: [.686, 0, .901],
    0x1D6B2: [.698, 0, .806],
    0x1D6B3: [.686, 0, 1.092],
    0x1D6B4: [.686, 0, .9],
    0x1D6B5: [.675, 0, .767],
    0x1D6B6: [.696, .01, .864],
    0x1D6B7: [.68, 0, .9],
    0x1D6B8: [.686, 0, .786],
    0x1D6BA: [.686, 0, .831],
    0x1D6BB: [.675, 0, .8],
    0x1D6BC: [.697, 0, .894],
    0x1D6BD: [.686, 0, .831],
    0x1D6BE: [.686, 0, .869],
    0x1D6BF: [.686, 0, .894],
    0x1D6C0: [.696, 0, .831],
    0x1D6C1: [.686, .024, .958],
    0x1D6E2: [.716, 0, .75, { sk: .139 }],
    0x1D6E3: [.683, 0, .759, { sk: .0833 }],
    0x1D6E4: [.68, 0, .615, { ic: .106, sk: .0833 }],
    0x1D6E5: [.716, 0, .833, { sk: .167 }],
    0x1D6E6: [.68, 0, .738, { ic: .026, sk: .0833 }],
    0x1D6E7: [.683, 0, .683, { ic: .04, sk: .0833 }],
    0x1D6E8: [.683, 0, .831, { ic: .057, sk: .0556 }],
    0x1D6E9: [.704, .022, .763, { sk: .0833 }],
    0x1D6EA: [.683, 0, .44, { ic: .064, sk: .111 }],
    0x1D6EB: [.683, 0, .849, { ic: .04, sk: .0556 }],
    0x1D6EC: [.716, 0, .694, { sk: .167 }],
    0x1D6ED: [.683, 0, .97, { ic: .081, sk: .0833 }],
    0x1D6EE: [.683, 0, .803, { ic: .085, sk: .0833 }],
    0x1D6EF: [.677, 0, .742, { ic: .035, sk: .0833 }],
    0x1D6F0: [.704, .022, .763, { sk: .0833 }],
    0x1D6F1: [.68, 0, .831, { ic: .056, sk: .0556 }],
    0x1D6F2: [.683, 0, .642, { ic: .109, sk: .0833 }],
    0x1D6F4: [.683, 0, .78, { ic: .026, sk: .0833 }],
    0x1D6F5: [.677, 0, .584, { ic: .12, sk: .0833 }],
    0x1D6F6: [.705, 0, .583, { ic: .117, sk: .0556 }],
    0x1D6F7: [.683, 0, .667, { sk: .0833 }],
    0x1D6F8: [.683, 0, .828, { ic: .024, sk: .0833 }],
    0x1D6F9: [.683, 0, .612, { ic: .08, sk: .0556 }],
    0x1D6FA: [.704, 0, .772, { ic: .014, sk: .0833 }],
    0x1D6FC: [.442, .011, .64, { sk: .0278 }],
    0x1D6FD: [.705, .194, .566, { sk: .0833 }],
    0x1D6FE: [.441, .216, .518, { ic: .025 }],
    0x1D6FF: [.717, .01, .444, { sk: .0556 }],
    0x1D700: [.452, .022, .466, { sk: .0833 }],
    0x1D701: [.704, .204, .438, { ic: .033, sk: .0833 }],
    0x1D702: [.442, .216, .497, { sk: .0556 }],
    0x1D703: [.705, .01, .469, { sk: .0833 }],
    0x1D704: [.442, .01, .354, { sk: .0556 }],
    0x1D705: [.442, .011, .576],
    0x1D706: [.694, .012, .583],
    0x1D707: [.442, .216, .603, { sk: .0278 }],
    0x1D708: [.442, 0, .494, { ic: .036, sk: .0278 }],
    0x1D709: [.704, .205, .438, { sk: .111 }],
    0x1D70A: [.441, .011, .485, { sk: .0556 }],
    0x1D70B: [.431, .011, .57],
    0x1D70C: [.442, .216, .517, { sk: .0833 }],
    0x1D70D: [.442, .107, .363, { ic: .042, sk: .0833 }],
    0x1D70E: [.431, .011, .571],
    0x1D70F: [.431, .013, .437, { ic: .08, sk: .0278 }],
    0x1D710: [.443, .01, .54, { sk: .0278 }],
    0x1D711: [.442, .218, .654, { sk: .0833 }],
    0x1D712: [.442, .204, .626, { sk: .0556 }],
    0x1D713: [.694, .205, .651, { sk: .111 }],
    0x1D714: [.443, .011, .622],
    0x1D715: [.715, .022, .531, { ic: .035, sk: .0833 }],
    0x1D716: [.431, .011, .406, { sk: .0556 }],
    0x1D717: [.705, .011, .591, { sk: .0833 }],
    0x1D718: [.434, .006, .667, { ic: .067 }],
    0x1D719: [.694, .205, .596, { sk: .0833 }],
    0x1D71A: [.442, .194, .517, { sk: .0833 }],
    0x1D71B: [.431, .01, .828],
    0x1D71C: [.711, 0, .869, { sk: .16 }],
    0x1D71D: [.686, 0, .866, { sk: .0958 }],
    0x1D71E: [.68, 0, .657, { ic: .12, sk: .0958 }],
    0x1D71F: [.711, 0, .958, { sk: .192 }],
    0x1D720: [.68, 0, .81, { ic: .015, sk: .0958 }],
    0x1D721: [.686, 0, .773, { ic: .032, sk: .0958 }],
    0x1D722: [.686, 0, .982, { ic: .045, sk: .0639 }],
    0x1D723: [.702, .017, .867, { sk: .0958 }],
    0x1D724: [.686, 0, .511, { ic: .062, sk: .128 }],
    0x1D725: [.686, 0, .971, { ic: .032, sk: .0639 }],
    0x1D726: [.711, 0, .806, { sk: .192 }],
    0x1D727: [.686, 0, 1.142, { ic: .077, sk: .0958 }],
    0x1D728: [.686, 0, .95, { ic: .077, sk: .0958 }],
    0x1D729: [.675, 0, .841, { ic: .026, sk: .0958 }],
    0x1D72A: [.703, .017, .837, { sk: .0958 }],
    0x1D72B: [.68, 0, .982, { ic: .044, sk: .0639 }],
    0x1D72C: [.686, 0, .723, { ic: .124, sk: .0958 }],
    0x1D72E: [.686, 0, .885, { ic: .017, sk: .0958 }],
    0x1D72F: [.675, 0, .637, { ic: .135, sk: .0958 }],
    0x1D730: [.703, 0, .671, { ic: .131, sk: .0639 }],
    0x1D731: [.686, 0, .767, { sk: .0958 }],
    0x1D732: [.686, 0, .947, { sk: .0958 }],
    0x1D733: [.686, 0, .714, { ic: .076, sk: .0639 }],
    0x1D734: [.703, 0, .879, { sk: .0958 }],
    0x1D736: [.452, .008, .761, { sk: .0319 }],
    0x1D737: [.701, .194, .66, { sk: .0958 }],
    0x1D738: [.451, .211, .59, { ic: .027 }],
    0x1D739: [.725, .008, .522, { sk: .0639 }],
    0x1D73A: [.461, .017, .529, { sk: .0958 }],
    0x1D73B: [.711, .202, .508, { ic: .013, sk: .0958 }],
    0x1D73C: [.452, .211, .6, { sk: .0639 }],
    0x1D73D: [.702, .008, .562, { sk: .0958 }],
    0x1D73E: [.452, .008, .412, { sk: .0639 }],
    0x1D73F: [.452, .008, .668],
    0x1D740: [.694, .013, .671],
    0x1D741: [.452, .211, .708, { sk: .0319 }],
    0x1D742: [.452, 0, .577, { ic: .031, sk: .0319 }],
    0x1D743: [.711, .201, .508, { sk: .128 }],
    0x1D744: [.452, .008, .585, { sk: .0639 }],
    0x1D745: [.444, .008, .682],
    0x1D746: [.451, .211, .612, { sk: .0958 }],
    0x1D747: [.451, .105, .424, { ic: .033, sk: .0958 }],
    0x1D748: [.444, .008, .686],
    0x1D749: [.444, .013, .521, { ic: .089, sk: .0319 }],
    0x1D74A: [.453, .008, .631, { sk: .0319 }],
    0x1D74B: [.452, .216, .747, { sk: .0958 }],
    0x1D74C: [.452, .201, .718, { sk: .0639 }],
    0x1D74D: [.694, .202, .758, { sk: .128 }],
    0x1D74E: [.453, .008, .718],
    0x1D74F: [.71, .017, .628, { ic: .029, sk: .0958 }],
    0x1D750: [.444, .007, .483, { sk: .0639 }],
    0x1D751: [.701, .008, .692, { sk: .0958 }],
    0x1D752: [.434, .006, .667, { ic: .067 }],
    0x1D753: [.694, .202, .712, { sk: .0958 }],
    0x1D754: [.451, .194, .612, { sk: .0958 }],
    0x1D755: [.444, .008, .975],
    0x1D756: [.694, 0, .733],
    0x1D757: [.694, 0, .733],
    0x1D758: [.691, 0, .581],
    0x1D759: [.694, 0, .917],
    0x1D75A: [.691, 0, .642],
    0x1D75B: [.694, 0, .672],
    0x1D75C: [.694, 0, .794],
    0x1D75D: [.716, .022, .856],
    0x1D75E: [.694, 0, .331],
    0x1D75F: [.694, 0, .764],
    0x1D760: [.694, 0, .672],
    0x1D761: [.694, 0, .978],
    0x1D762: [.694, 0, .794],
    0x1D763: [.688, 0, .733],
    0x1D764: [.716, .022, .794],
    0x1D765: [.691, 0, .794],
    0x1D766: [.694, 0, .703],
    0x1D768: [.694, 0, .794],
    0x1D769: [.688, 0, .733],
    0x1D76A: [.715, 0, .856],
    0x1D76B: [.694, 0, .794],
    0x1D76C: [.694, 0, .733],
    0x1D76D: [.694, 0, .856],
    0x1D76E: [.716, 0, .794],
    0x1D7CE: [.654, .01, .575],
    0x1D7CF: [.655, 0, .575],
    0x1D7D0: [.654, 0, .575],
    0x1D7D1: [.655, .011, .575],
    0x1D7D2: [.656, 0, .575],
    0x1D7D3: [.655, .011, .575],
    0x1D7D4: [.655, .011, .575],
    0x1D7D5: [.676, .011, .575],
    0x1D7D6: [.654, .011, .575],
    0x1D7D7: [.654, .011, .575],
    0x1D7E2: [.678, .022, .5],
    0x1D7E3: [.678, 0, .5],
    0x1D7E4: [.677, 0, .5],
    0x1D7E5: [.678, .022, .5],
    0x1D7E6: [.656, 0, .5],
    0x1D7E7: [.656, .021, .5],
    0x1D7E8: [.677, .022, .5],
    0x1D7E9: [.656, .011, .5],
    0x1D7EA: [.678, .022, .5],
    0x1D7EB: [.677, .022, .5],
    0x1D7EC: [.715, .022, .55],
    0x1D7ED: [.716, 0, .55],
    0x1D7EE: [.716, 0, .55],
    0x1D7EF: [.716, .022, .55],
    0x1D7F0: [.694, 0, .55],
    0x1D7F1: [.694, .022, .55],
    0x1D7F2: [.716, .022, .55],
    0x1D7F3: [.695, .011, .55],
    0x1D7F4: [.715, .022, .55],
    0x1D7F5: [.716, .022, .55],
    0x1D7F6: [.621, .01, .525],
    0x1D7F7: [.622, 0, .525],
    0x1D7F8: [.622, 0, .525],
    0x1D7F9: [.622, .011, .525],
    0x1D7FA: [.624, 0, .525],
    0x1D7FB: [.611, .01, .525],
    0x1D7FC: [.622, .011, .525],
    0x1D7FD: [.627, .01, .525],
    0x1D7FE: [.621, .01, .525],
    0x1D7FF: [.622, .011, .525],
};
//# sourceMappingURL=normal.js.map

/***/ }),

/***/ 75656:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifBoldItalic = void 0;
exports.sansSerifBoldItalic = {
    0x131: [.458, 0, .256],
    0x237: [.458, .205, .286],
};
//# sourceMappingURL=sans-serif-bold-italic.js.map

/***/ }),

/***/ 86650:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifBold = void 0;
exports.sansSerifBold = {
    0x21: [.694, 0, .367],
    0x22: [.694, -0.442, .558],
    0x23: [.694, .193, .917],
    0x24: [.75, .056, .55],
    0x25: [.75, .056, 1.029],
    0x26: [.716, .022, .831],
    0x27: [.694, -0.442, .306],
    0x28: [.75, .249, .428],
    0x29: [.75, .25, .428],
    0x2A: [.75, -0.293, .55],
    0x2B: [.617, .116, .856],
    0x2C: [.146, .106, .306],
    0x2D: [.273, -0.186, .367],
    0x2E: [.146, 0, .306],
    0x2F: [.75, .249, .55],
    0x3A: [.458, 0, .306],
    0x3B: [.458, .106, .306],
    0x3D: [.407, -0.094, .856],
    0x3F: [.705, 0, .519],
    0x40: [.704, .011, .733],
    0x5B: [.75, .25, .343],
    0x5D: [.75, .25, .343],
    0x5E: [.694, -0.537, .55],
    0x5F: [-0.023, .11, .55],
    0x7E: [.344, -0.198, .55],
    0x131: [.458, 0, .256],
    0x237: [.458, .205, .286],
    0x300: [.694, -0.537, 0],
    0x301: [.694, -0.537, 0],
    0x302: [.694, -0.537, 0],
    0x303: [.694, -0.548, 0],
    0x304: [.66, -0.56, 0],
    0x306: [.694, -0.552, 0],
    0x307: [.695, -0.596, 0],
    0x308: [.695, -0.595, 0],
    0x30A: [.694, -0.538, 0],
    0x30B: [.694, -0.537, 0],
    0x30C: [.657, -0.5, 0],
    0x2013: [.327, -0.24, .55],
    0x2014: [.327, -0.24, 1.1],
    0x2015: [.327, -0.24, 1.1],
    0x2017: [-0.023, .11, .55],
    0x2018: [.694, -0.443, .306],
    0x2019: [.694, -0.442, .306],
    0x201C: [.694, -0.443, .558],
    0x201D: [.694, -0.442, .558],
    0x2044: [.75, .249, .55],
    0x2206: [.694, 0, .917],
};
//# sourceMappingURL=sans-serif-bold.js.map

/***/ }),

/***/ 51518:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerifItalic = void 0;
exports.sansSerifItalic = {
    0x21: [.694, 0, .319, { ic: .036 }],
    0x22: [.694, -0.471, .5],
    0x23: [.694, .194, .833, { ic: .018 }],
    0x24: [.75, .056, .5, { ic: .065 }],
    0x25: [.75, .056, .833],
    0x26: [.716, .022, .758],
    0x27: [.694, -0.471, .278, { ic: .057 }],
    0x28: [.75, .25, .389, { ic: .102 }],
    0x29: [.75, .25, .389],
    0x2A: [.75, -0.306, .5, { ic: .068 }],
    0x2B: [.583, .083, .778],
    0x2C: [.098, .125, .278],
    0x2D: [.259, -0.186, .333],
    0x2E: [.098, 0, .278],
    0x2F: [.75, .25, .5, { ic: .1 }],
    0x30: [.678, .022, .5, { ic: .049 }],
    0x31: [.678, 0, .5],
    0x32: [.678, 0, .5, { ic: .051 }],
    0x33: [.678, .022, .5, { ic: .044 }],
    0x34: [.656, 0, .5, { ic: .021 }],
    0x35: [.656, .022, .5, { ic: .055 }],
    0x36: [.678, .022, .5, { ic: .048 }],
    0x37: [.656, .011, .5, { ic: .096 }],
    0x38: [.678, .022, .5, { ic: .054 }],
    0x39: [.677, .022, .5, { ic: .045 }],
    0x3A: [.444, 0, .278],
    0x3B: [.444, .125, .278],
    0x3D: [.37, -0.13, .778, { ic: .018 }],
    0x3F: [.704, 0, .472, { ic: .064 }],
    0x40: [.705, .01, .667, { ic: .04 }],
    0x5B: [.75, .25, .289, { ic: .136 }],
    0x5D: [.75, .25, .289, { ic: .064 }],
    0x5E: [.694, -0.527, .5, { ic: .033 }],
    0x5F: [-0.038, .114, .5, { ic: .065 }],
    0x7E: [.327, -0.193, .5, { ic: .06 }],
    0x131: [.444, 0, .239, { ic: .019 }],
    0x237: [.444, .204, .267, { ic: .019 }],
    0x300: [.694, -0.527, 0],
    0x301: [.694, -0.527, 0, { ic: .063 }],
    0x302: [.694, -0.527, 0, { ic: .033 }],
    0x303: [.677, -0.543, 0, { ic: .06 }],
    0x304: [.631, -0.552, 0, { ic: .064 }],
    0x306: [.694, -0.508, 0, { ic: .073 }],
    0x307: [.68, -0.576, 0],
    0x308: [.68, -0.582, 0, { ic: .04 }],
    0x30A: [.693, -0.527, 0],
    0x30B: [.694, -0.527, 0, { ic: .063 }],
    0x30C: [.654, -0.487, 0, { ic: .06 }],
    0x391: [.694, 0, .667],
    0x392: [.694, 0, .667, { ic: .029 }],
    0x393: [.691, 0, .542, { ic: .104 }],
    0x394: [.694, 0, .833],
    0x395: [.691, 0, .597, { ic: .091 }],
    0x396: [.694, 0, .611, { ic: .091 }],
    0x397: [.694, 0, .708, { ic: .06 }],
    0x398: [.715, .022, .778, { ic: .026 }],
    0x399: [.694, 0, .278, { ic: .06 }],
    0x39A: [.694, 0, .694, { ic: .091 }],
    0x39B: [.694, 0, .611],
    0x39C: [.694, 0, .875, { ic: .054 }],
    0x39D: [.694, 0, .708, { ic: .058 }],
    0x39E: [.688, 0, .667, { ic: .098 }],
    0x39F: [.716, .022, .736, { ic: .027 }],
    0x3A0: [.691, 0, .708, { ic: .06 }],
    0x3A1: [.694, 0, .639, { ic: .051 }],
    0x3A3: [.694, 0, .722, { ic: .091 }],
    0x3A4: [.688, 0, .681, { ic: .109 }],
    0x3A5: [.716, 0, .778, { ic: .065 }],
    0x3A6: [.694, 0, .722, { ic: .021 }],
    0x3A7: [.694, 0, .667, { ic: .091 }],
    0x3A8: [.694, 0, .778, { ic: .076 }],
    0x3A9: [.716, 0, .722, { ic: .047 }],
    0x2013: [.312, -0.236, .5, { ic: .065 }],
    0x2014: [.312, -0.236, 1, { ic: .065 }],
    0x2015: [.312, -0.236, 1, { ic: .065 }],
    0x2017: [-0.038, .114, .5, { ic: .065 }],
    0x2018: [.694, -0.471, .278, { ic: .058 }],
    0x2019: [.694, -0.471, .278, { ic: .057 }],
    0x201C: [.694, -0.471, .5, { ic: .114 }],
    0x201D: [.694, -0.471, .5],
    0x2044: [.75, .25, .5, { ic: .1 }],
    0x2206: [.694, 0, .833],
};
//# sourceMappingURL=sans-serif-italic.js.map

/***/ }),

/***/ 71290:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.sansSerif = void 0;
exports.sansSerif = {
    0x21: [.694, 0, .319],
    0x22: [.694, -0.471, .5],
    0x23: [.694, .194, .833],
    0x24: [.75, .056, .5],
    0x25: [.75, .056, .833],
    0x26: [.716, .022, .758],
    0x27: [.694, -0.471, .278],
    0x28: [.75, .25, .389],
    0x29: [.75, .25, .389],
    0x2A: [.75, -0.306, .5],
    0x2B: [.583, .082, .778],
    0x2C: [.098, .125, .278],
    0x2D: [.259, -0.186, .333],
    0x2E: [.098, 0, .278],
    0x2F: [.75, .25, .5],
    0x3A: [.444, 0, .278],
    0x3B: [.444, .125, .278],
    0x3D: [.37, -0.13, .778],
    0x3F: [.704, 0, .472],
    0x40: [.704, .011, .667],
    0x5B: [.75, .25, .289],
    0x5D: [.75, .25, .289],
    0x5E: [.694, -0.527, .5],
    0x5F: [-0.038, .114, .5],
    0x7E: [.327, -0.193, .5],
    0x131: [.444, 0, .239],
    0x237: [.444, .205, .267],
    0x300: [.694, -0.527, 0],
    0x301: [.694, -0.527, 0],
    0x302: [.694, -0.527, 0],
    0x303: [.677, -0.543, 0],
    0x304: [.631, -0.552, 0],
    0x306: [.694, -0.508, 0],
    0x307: [.68, -0.576, 0],
    0x308: [.68, -0.582, 0],
    0x30A: [.694, -0.527, 0],
    0x30B: [.694, -0.527, 0],
    0x30C: [.654, -0.487, 0],
    0x391: [.694, 0, .667],
    0x392: [.694, 0, .667],
    0x393: [.691, 0, .542],
    0x394: [.694, 0, .833],
    0x395: [.691, 0, .597],
    0x396: [.694, 0, .611],
    0x397: [.694, 0, .708],
    0x398: [.716, .021, .778],
    0x399: [.694, 0, .278],
    0x39A: [.694, 0, .694],
    0x39B: [.694, 0, .611],
    0x39C: [.694, 0, .875],
    0x39D: [.694, 0, .708],
    0x39E: [.688, 0, .667],
    0x39F: [.715, .022, .736],
    0x3A0: [.691, 0, .708],
    0x3A1: [.694, 0, .639],
    0x3A3: [.694, 0, .722],
    0x3A4: [.688, 0, .681],
    0x3A5: [.716, 0, .778],
    0x3A6: [.694, 0, .722],
    0x3A7: [.694, 0, .667],
    0x3A8: [.694, 0, .778],
    0x3A9: [.716, 0, .722],
    0x2013: [.312, -0.236, .5],
    0x2014: [.312, -0.236, 1],
    0x2015: [.312, -0.236, 1],
    0x2017: [-0.038, .114, .5],
    0x2018: [.694, -0.471, .278],
    0x2019: [.694, -0.471, .278],
    0x201C: [.694, -0.471, .5],
    0x201D: [.694, -0.471, .5],
    0x2044: [.75, .25, .5],
    0x2206: [.694, 0, .833],
};
//# sourceMappingURL=sans-serif.js.map

/***/ }),

/***/ 89520:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.scriptBold = void 0;
exports.scriptBold = {};
//# sourceMappingURL=script-bold.js.map

/***/ }),

/***/ 10970:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.script = void 0;
exports.script = {};
//# sourceMappingURL=script.js.map

/***/ }),

/***/ 34664:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.smallop = void 0;
exports.smallop = {
    0x28: [.85, .349, .458],
    0x29: [.85, .349, .458],
    0x2F: [.85, .349, .578],
    0x5B: [.85, .349, .417],
    0x5C: [.85, .349, .578],
    0x5D: [.85, .349, .417],
    0x7B: [.85, .349, .583],
    0x7D: [.85, .349, .583],
    0x2C6: [.744, -0.551, .556],
    0x2DC: [.722, -0.597, .556],
    0x302: [.744, -0.551, 0],
    0x303: [.722, -0.597, 0],
    0x2016: [.602, 0, .778],
    0x2044: [.85, .349, .578],
    0x2191: [.6, 0, .667],
    0x2193: [.6, 0, .667],
    0x21D1: [.599, 0, .778],
    0x21D3: [.6, 0, .778],
    0x220F: [.75, .25, .944],
    0x2210: [.75, .25, .944],
    0x2211: [.75, .25, 1.056],
    0x221A: [.85, .35, 1, { ic: .02 }],
    0x2223: [.627, .015, .333],
    0x2225: [.627, .015, .556],
    0x222B: [.805, .306, .472, { ic: .138 }],
    0x222C: [.805, .306, .819, { ic: .138 }],
    0x222D: [.805, .306, 1.166, { ic: .138 }],
    0x222E: [.805, .306, .472, { ic: .138 }],
    0x22C0: [.75, .249, .833],
    0x22C1: [.75, .249, .833],
    0x22C2: [.75, .249, .833],
    0x22C3: [.75, .249, .833],
    0x2308: [.85, .349, .472],
    0x2309: [.85, .349, .472],
    0x230A: [.85, .349, .472],
    0x230B: [.85, .349, .472],
    0x2329: [.85, .35, .472],
    0x232A: [.85, .35, .472],
    0x23D0: [.602, 0, .667],
    0x2758: [.627, .015, .333],
    0x27E8: [.85, .35, .472],
    0x27E9: [.85, .35, .472],
    0x2A00: [.75, .25, 1.111],
    0x2A01: [.75, .25, 1.111],
    0x2A02: [.75, .25, 1.111],
    0x2A04: [.75, .249, .833],
    0x2A06: [.75, .249, .833],
    0x2A0C: [.805, .306, 1.638, { ic: .138 }],
    0x3008: [.85, .35, .472],
    0x3009: [.85, .35, .472],
};
//# sourceMappingURL=smallop.js.map

/***/ }),

/***/ 43545:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texCalligraphicBold = void 0;
exports.texCalligraphicBold = {
    0x41: [.751, .049, .921, { ic: .068, sk: .224 }],
    0x42: [.705, .017, .748, { sk: .16 }],
    0x43: [.703, .02, .613, { sk: .16 }],
    0x44: [.686, 0, .892, { sk: .0958 }],
    0x45: [.703, .016, .607, { ic: .02, sk: .128 }],
    0x46: [.686, .03, .814, { ic: .116, sk: .128 }],
    0x47: [.703, .113, .682, { sk: .128 }],
    0x48: [.686, .048, .987, { sk: .128 }],
    0x49: [.686, 0, .642, { ic: .104, sk: .0319 }],
    0x4A: [.686, .114, .779, { ic: .158, sk: .192 }],
    0x4B: [.703, .017, .871, { sk: .0639 }],
    0x4C: [.703, .017, .788, { sk: .16 }],
    0x4D: [.703, .049, 1.378, { sk: .16 }],
    0x4E: [.84, .049, .937, { ic: .168, sk: .0958 }],
    0x4F: [.703, .017, .906, { sk: .128 }],
    0x50: [.686, .067, .81, { ic: .036, sk: .0958 }],
    0x51: [.703, .146, .939, { sk: .128 }],
    0x52: [.686, .017, .99, { sk: .0958 }],
    0x53: [.703, .016, .696, { ic: .025, sk: .16 }],
    0x54: [.72, .069, .644, { ic: .303, sk: .0319 }],
    0x55: [.686, .024, .715, { ic: .056, sk: .0958 }],
    0x56: [.686, .077, .737, { ic: .037, sk: .0319 }],
    0x57: [.686, .077, 1.169, { ic: .037, sk: .0958 }],
    0x58: [.686, 0, .817, { ic: .089, sk: .16 }],
    0x59: [.686, .164, .759, { ic: .038, sk: .0958 }],
    0x5A: [.686, 0, .818, { ic: .035, sk: .16 }],
    0x131: [.452, .008, .394, { sk: .0319 }],
    0x237: [.451, .201, .439, { sk: .0958 }],
};
//# sourceMappingURL=tex-calligraphic-bold.js.map

/***/ }),

/***/ 55389:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texCalligraphic = void 0;
exports.texCalligraphic = {
    0x41: [.728, .05, .798, { ic: .021, sk: .194 }],
    0x42: [.705, .022, .657, { sk: .139 }],
    0x43: [.705, .025, .527, { sk: .139 }],
    0x44: [.683, 0, .771, { sk: .0833 }],
    0x45: [.705, .022, .528, { ic: .036, sk: .111 }],
    0x46: [.683, .032, .719, { ic: .11, sk: .111 }],
    0x47: [.704, .119, .595, { sk: .111 }],
    0x48: [.683, .048, .845, { sk: .111 }],
    0x49: [.683, 0, .545, { ic: .097, sk: .0278 }],
    0x4A: [.683, .119, .678, { ic: .161, sk: .167 }],
    0x4B: [.705, .022, .762, { sk: .0556 }],
    0x4C: [.705, .022, .69, { sk: .139 }],
    0x4D: [.705, .05, 1.201, { sk: .139 }],
    0x4E: [.789, .05, .82, { ic: .159, sk: .0833 }],
    0x4F: [.705, .022, .796, { sk: .111 }],
    0x50: [.683, .057, .696, { ic: .037, sk: .0833 }],
    0x51: [.705, .131, .817, { sk: .111 }],
    0x52: [.682, .022, .848, { sk: .0833 }],
    0x53: [.705, .022, .606, { ic: .036, sk: .139 }],
    0x54: [.717, .068, .545, { ic: .288, sk: .0278 }],
    0x55: [.683, .028, .626, { ic: .061, sk: .0833 }],
    0x56: [.683, .052, .613, { ic: .045, sk: .0278 }],
    0x57: [.683, .053, .988, { ic: .046, sk: .0833 }],
    0x58: [.683, 0, .713, { ic: .094, sk: .139 }],
    0x59: [.683, .143, .668, { ic: .046, sk: .0833 }],
    0x5A: [.683, 0, .725, { ic: .042, sk: .139 }],
};
//# sourceMappingURL=tex-calligraphic.js.map

/***/ }),

/***/ 1641:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texMathit = void 0;
exports.texMathit = {
    0x41: [.716, 0, .743],
    0x42: [.683, 0, .704],
    0x43: [.705, .021, .716],
    0x44: [.683, 0, .755],
    0x45: [.68, 0, .678],
    0x46: [.68, 0, .653],
    0x47: [.705, .022, .774],
    0x48: [.683, 0, .743],
    0x49: [.683, 0, .386],
    0x4A: [.683, .021, .525],
    0x4B: [.683, 0, .769],
    0x4C: [.683, 0, .627],
    0x4D: [.683, 0, .897],
    0x4E: [.683, 0, .743],
    0x4F: [.704, .022, .767],
    0x50: [.683, 0, .678],
    0x51: [.704, .194, .767],
    0x52: [.683, .022, .729],
    0x53: [.705, .022, .562],
    0x54: [.677, 0, .716],
    0x55: [.683, .022, .743],
    0x56: [.683, .022, .743],
    0x57: [.683, .022, .999],
    0x58: [.683, 0, .743],
    0x59: [.683, 0, .743],
    0x5A: [.683, 0, .613],
    0x61: [.442, .011, .511],
    0x62: [.694, .011, .46],
    0x63: [.441, .01, .46],
    0x64: [.694, .011, .511],
    0x65: [.442, .01, .46],
    0x66: [.705, .204, .307],
    0x67: [.442, .205, .46],
    0x68: [.694, .011, .511],
    0x69: [.656, .01, .307],
    0x6A: [.656, .204, .307],
    0x6B: [.694, .011, .46],
    0x6C: [.694, .011, .256],
    0x6D: [.442, .011, .818],
    0x6E: [.442, .011, .562],
    0x6F: [.442, .011, .511],
    0x70: [.442, .194, .511],
    0x71: [.442, .194, .46],
    0x72: [.442, .011, .422],
    0x73: [.442, .011, .409],
    0x74: [.626, .011, .332],
    0x75: [.441, .011, .537],
    0x76: [.443, .01, .46],
    0x77: [.443, .011, .664],
    0x78: [.442, .011, .464],
    0x79: [.441, .205, .486],
    0x7A: [.442, .011, .409],
};
//# sourceMappingURL=tex-mathit.js.map

/***/ }),

/***/ 33944:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texOldstyleBold = void 0;
exports.texOldstyleBold = {
    0x30: [.46, .017, .575],
    0x31: [.461, 0, .575],
    0x32: [.46, 0, .575],
    0x33: [.461, .211, .575],
    0x34: [.469, .194, .575],
    0x35: [.461, .211, .575],
    0x36: [.66, .017, .575],
    0x37: [.476, .211, .575],
    0x38: [.661, .017, .575],
    0x39: [.461, .21, .575],
    0x41: [.751, .049, .921, { ic: .068, sk: .224 }],
    0x42: [.705, .017, .748, { sk: .16 }],
    0x43: [.703, .02, .613, { sk: .16 }],
    0x44: [.686, 0, .892, { sk: .0958 }],
    0x45: [.703, .016, .607, { ic: .02, sk: .128 }],
    0x46: [.686, .03, .814, { ic: .116, sk: .128 }],
    0x47: [.703, .113, .682, { sk: .128 }],
    0x48: [.686, .048, .987, { sk: .128 }],
    0x49: [.686, 0, .642, { ic: .104, sk: .0319 }],
    0x4A: [.686, .114, .779, { ic: .158, sk: .192 }],
    0x4B: [.703, .017, .871, { sk: .0639 }],
    0x4C: [.703, .017, .788, { sk: .16 }],
    0x4D: [.703, .049, 1.378, { sk: .16 }],
    0x4E: [.84, .049, .937, { ic: .168, sk: .0958 }],
    0x4F: [.703, .017, .906, { sk: .128 }],
    0x50: [.686, .067, .81, { ic: .036, sk: .0958 }],
    0x51: [.703, .146, .939, { sk: .128 }],
    0x52: [.686, .017, .99, { sk: .0958 }],
    0x53: [.703, .016, .696, { ic: .025, sk: .16 }],
    0x54: [.72, .069, .644, { ic: .303, sk: .0319 }],
    0x55: [.686, .024, .715, { ic: .056, sk: .0958 }],
    0x56: [.686, .077, .737, { ic: .037, sk: .0319 }],
    0x57: [.686, .077, 1.169, { ic: .037, sk: .0958 }],
    0x58: [.686, 0, .817, { ic: .089, sk: .16 }],
    0x59: [.686, .164, .759, { ic: .038, sk: .0958 }],
    0x5A: [.686, 0, .818, { ic: .035, sk: .16 }],
};
//# sourceMappingURL=tex-oldstyle-bold.js.map

/***/ }),

/***/ 98690:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texOldstyle = void 0;
exports.texOldstyle = {
    0x30: [.452, .022, .5],
    0x31: [.453, 0, .5],
    0x32: [.453, 0, .5],
    0x33: [.452, .216, .5],
    0x34: [.464, .194, .5],
    0x35: [.453, .216, .5],
    0x36: [.665, .022, .5],
    0x37: [.463, .216, .5],
    0x38: [.666, .021, .5],
    0x39: [.453, .216, .5],
    0x41: [.728, .05, .798, { ic: .021, sk: .194 }],
    0x42: [.705, .022, .657, { sk: .139 }],
    0x43: [.705, .025, .527, { sk: .139 }],
    0x44: [.683, 0, .771, { sk: .0833 }],
    0x45: [.705, .022, .528, { ic: .036, sk: .111 }],
    0x46: [.683, .032, .719, { ic: .11, sk: .111 }],
    0x47: [.704, .119, .595, { sk: .111 }],
    0x48: [.683, .048, .845, { sk: .111 }],
    0x49: [.683, 0, .545, { ic: .097, sk: .0278 }],
    0x4A: [.683, .119, .678, { ic: .161, sk: .167 }],
    0x4B: [.705, .022, .762, { sk: .0556 }],
    0x4C: [.705, .022, .69, { sk: .139 }],
    0x4D: [.705, .05, 1.201, { sk: .139 }],
    0x4E: [.789, .05, .82, { ic: .159, sk: .0833 }],
    0x4F: [.705, .022, .796, { sk: .111 }],
    0x50: [.683, .057, .696, { ic: .037, sk: .0833 }],
    0x51: [.705, .131, .817, { sk: .111 }],
    0x52: [.682, .022, .848, { sk: .0833 }],
    0x53: [.705, .022, .606, { ic: .036, sk: .139 }],
    0x54: [.717, .068, .545, { ic: .288, sk: .0278 }],
    0x55: [.683, .028, .626, { ic: .061, sk: .0833 }],
    0x56: [.683, .052, .613, { ic: .045, sk: .0278 }],
    0x57: [.683, .053, .988, { ic: .046, sk: .0833 }],
    0x58: [.683, 0, .713, { ic: .094, sk: .139 }],
    0x59: [.683, .143, .668, { ic: .046, sk: .0833 }],
    0x5A: [.683, 0, .725, { ic: .042, sk: .139 }],
};
//# sourceMappingURL=tex-oldstyle.js.map

/***/ }),

/***/ 96464:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texSize3 = void 0;
exports.texSize3 = {
    0x28: [1.45, .949, .736],
    0x29: [1.45, .949, .736],
    0x2F: [1.45, .949, 1.044],
    0x5B: [1.45, .949, .528],
    0x5C: [1.45, .949, 1.044],
    0x5D: [1.45, .949, .528],
    0x7B: [1.45, .949, .75],
    0x7D: [1.45, .949, .75],
    0x2C6: [.772, -0.564, 1.444],
    0x2DC: [.749, -0.61, 1.444],
    0x302: [.772, -0.564, 0],
    0x303: [.749, -0.61, 0],
    0x2044: [1.45, .949, 1.044],
    0x221A: [1.45, .95, 1, { ic: .02 }],
    0x2308: [1.45, .949, .583],
    0x2309: [1.45, .949, .583],
    0x230A: [1.45, .949, .583],
    0x230B: [1.45, .949, .583],
    0x2329: [1.45, .95, .75],
    0x232A: [1.45, .949, .75],
    0x27E8: [1.45, .95, .75],
    0x27E9: [1.45, .949, .75],
    0x3008: [1.45, .95, .75],
    0x3009: [1.45, .949, .75],
};
//# sourceMappingURL=tex-size3.js.map

/***/ }),

/***/ 58905:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texSize4 = void 0;
exports.texSize4 = {
    0x28: [1.75, 1.249, .792],
    0x29: [1.75, 1.249, .792],
    0x2F: [1.75, 1.249, 1.278],
    0x5B: [1.75, 1.249, .583],
    0x5C: [1.75, 1.249, 1.278],
    0x5D: [1.75, 1.249, .583],
    0x7B: [1.75, 1.249, .806],
    0x7D: [1.75, 1.249, .806],
    0x2C6: [.845, -0.561, 1.889, { ic: .013 }],
    0x2DC: [.823, -0.583, 1.889],
    0x302: [.845, -0.561, 0, { ic: .013 }],
    0x303: [.823, -0.583, 0],
    0x2044: [1.75, 1.249, 1.278],
    0x221A: [1.75, 1.25, 1, { ic: .02 }],
    0x2308: [1.75, 1.249, .639],
    0x2309: [1.75, 1.249, .639],
    0x230A: [1.75, 1.249, .639],
    0x230B: [1.75, 1.249, .639],
    0x2329: [1.75, 1.248, .806],
    0x232A: [1.75, 1.248, .806],
    0x239B: [1.154, .655, .875],
    0x239C: [.61, .01, .875],
    0x239D: [1.165, .644, .875],
    0x239E: [1.154, .655, .875],
    0x239F: [.61, .01, .875],
    0x23A0: [1.165, .644, .875],
    0x23A1: [1.154, .645, .667],
    0x23A2: [.602, 0, .667],
    0x23A3: [1.155, .644, .667],
    0x23A4: [1.154, .645, .667],
    0x23A5: [.602, 0, .667],
    0x23A6: [1.155, .644, .667],
    0x23A7: [.899, .01, .889],
    0x23A8: [1.16, .66, .889],
    0x23A9: [.01, .899, .889],
    0x23AA: [.29, .015, .889],
    0x23AB: [.899, .01, .889],
    0x23AC: [1.16, .66, .889],
    0x23AD: [.01, .899, .889],
    0x23B7: [.935, .885, 1.056],
    0x27E8: [1.75, 1.248, .806],
    0x27E9: [1.75, 1.248, .806],
    0x3008: [1.75, 1.248, .806],
    0x3009: [1.75, 1.248, .806],
    0xE000: [.625, .014, 1.056],
    0xE001: [.605, .014, 1.056, { ic: .02 }],
    0xE150: [.12, .213, .45, { ic: .01 }],
    0xE151: [.12, .213, .45, { ic: .024 }],
    0xE152: [.333, 0, .45, { ic: .01 }],
    0xE153: [.333, 0, .45, { ic: .024 }],
    0xE154: [.32, .2, .4, { ic: .01 }],
    0xE155: [.333, 0, .9, { ic: .01 }],
    0xE156: [.12, .213, .9, { ic: .01 }],
};
//# sourceMappingURL=tex-size4.js.map

/***/ }),

/***/ 46527:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.texVariant = void 0;
exports.texVariant = {
    0x2C6: [.845, -0.561, 2.333, { ic: .013 }],
    0x2DC: [.899, -0.628, 2.333],
    0x302: [.845, -0.561, 0, { ic: .013 }],
    0x303: [.899, -0.628, 0],
    0x3F0: [.434, .006, .667, { ic: .067 }],
    0x210F: [.695, .013, .54, { ic: .022 }],
    0x2190: [.437, -0.064, .5],
    0x2192: [.437, -0.064, .5],
    0x21CC: [.514, .014, 1],
    0x2204: [.86, .166, .556],
    0x2205: [.587, 0, .778],
    0x2212: [.27, -0.23, .5],
    0x2216: [.43, .023, .778],
    0x221D: [.472, -0.028, .778],
    0x2223: [.43, .023, .222],
    0x2224: [.43, .023, .222, { ic: .018 }],
    0x2225: [.431, .023, .389],
    0x2226: [.431, .024, .389, { ic: .018 }],
    0x223C: [.365, -0.132, .778],
    0x2248: [.481, -0.05, .778],
    0x2268: [.752, .284, .778],
    0x2269: [.752, .284, .778],
    0x2270: [.919, .421, .778],
    0x2271: [.919, .421, .778],
    0x2288: [.828, .33, .778],
    0x2289: [.828, .33, .778],
    0x228A: [.634, .255, .778],
    0x228B: [.634, .254, .778],
    0x22A8: [.694, 0, .611],
    0x22C5: [.189, 0, .278],
    0x2322: [.378, -0.122, .778],
    0x2323: [.378, -0.143, .778],
    0x25B3: [.575, .02, .722],
    0x25BD: [.576, .019, .722],
    0x2A87: [.801, .303, .778],
    0x2A88: [.801, .303, .778],
    0x2ACB: [.752, .332, .778],
    0x2ACC: [.752, .333, .778],
};
//# sourceMappingURL=tex-variant.js.map

/***/ }),

/***/ 56780:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.px = exports.emRounded = exports.em = exports.percent = exports.length2em = exports.MATHSPACE = exports.RELUNITS = exports.UNITS = exports.BIGDIMEN = void 0;
exports.BIGDIMEN = 1000000;
exports.UNITS = {
    px: 1,
    'in': 96,
    cm: 96 / 2.54,
    mm: 96 / 25.4
};
exports.RELUNITS = {
    em: 1,
    ex: .431,
    pt: 1 / 10,
    pc: 12 / 10,
    mu: 1 / 18
};
exports.MATHSPACE = {
    veryverythinmathspace: 1 / 18,
    verythinmathspace: 2 / 18,
    thinmathspace: 3 / 18,
    mediummathspace: 4 / 18,
    thickmathspace: 5 / 18,
    verythickmathspace: 6 / 18,
    veryverythickmathspace: 7 / 18,
    negativeveryverythinmathspace: -1 / 18,
    negativeverythinmathspace: -2 / 18,
    negativethinmathspace: -3 / 18,
    negativemediummathspace: -4 / 18,
    negativethickmathspace: -5 / 18,
    negativeverythickmathspace: -6 / 18,
    negativeveryverythickmathspace: -7 / 18,
    thin: .04,
    medium: .06,
    thick: .1,
    normal: 1,
    big: 2,
    small: 1 / Math.sqrt(2),
    infinity: exports.BIGDIMEN
};
function length2em(length, size, scale, em) {
    if (size === void 0) { size = 0; }
    if (scale === void 0) { scale = 1; }
    if (em === void 0) { em = 16; }
    if (typeof length !== 'string') {
        length = String(length);
    }
    if (length === '' || length == null) {
        return size;
    }
    if (exports.MATHSPACE[length]) {
        return exports.MATHSPACE[length];
    }
    var match = length.match(/^\s*([-+]?(?:\.\d+|\d+(?:\.\d*)?))?(pt|em|ex|mu|px|pc|in|mm|cm|%)?/);
    if (!match) {
        return size;
    }
    var m = parseFloat(match[1] || '1'), unit = match[2];
    if (exports.UNITS.hasOwnProperty(unit)) {
        return m * exports.UNITS[unit] / em / scale;
    }
    if (exports.RELUNITS.hasOwnProperty(unit)) {
        return m * exports.RELUNITS[unit];
    }
    if (unit === '%') {
        return m / 100 * size;
    }
    return m * size;
}
exports.length2em = length2em;
function percent(m) {
    return (100 * m).toFixed(1).replace(/\.?0+$/, '') + '%';
}
exports.percent = percent;
function em(m) {
    if (Math.abs(m) < .001)
        return '0';
    return (m.toFixed(3).replace(/\.?0+$/, '')) + 'em';
}
exports.em = em;
function emRounded(m, em) {
    if (em === void 0) { em = 16; }
    m = (Math.round(m * em) + .05) / em;
    if (Math.abs(m) < .001)
        return '0em';
    return m.toFixed(3).replace(/\.?0+$/, '') + 'em';
}
exports.emRounded = emRounded;
function px(m, M, em) {
    if (M === void 0) { M = -exports.BIGDIMEN; }
    if (em === void 0) { em = 16; }
    m *= em;
    if (M && m < M)
        m = M;
    if (Math.abs(m) < .1)
        return '0';
    return m.toFixed(1).replace(/\.0$/, '') + 'px';
}
exports.px = px;
//# sourceMappingURL=lengths.js.map

/***/ })

}]);
//# sourceMappingURL=2065.e9b5d8d0a8bec3304454.js.map?v=e9b5d8d0a8bec3304454