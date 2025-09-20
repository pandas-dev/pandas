"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8285,4498],{

/***/ 78285:
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
exports.SafeHandler = exports.SafeMathDocumentMixin = void 0;
var safe_js_1 = __webpack_require__(19290);
function SafeMathDocumentMixin(BaseDocument) {
    var _a;
    return _a = (function (_super) {
            __extends(class_1, _super);
            function class_1() {
                var e_1, _a;
                var args = [];
                for (var _i = 0; _i < arguments.length; _i++) {
                    args[_i] = arguments[_i];
                }
                var _this = _super.apply(this, __spreadArray([], __read(args), false)) || this;
                _this.safe = new _this.options.SafeClass(_this, _this.options.safeOptions);
                var ProcessBits = _this.constructor.ProcessBits;
                if (!ProcessBits.has('safe')) {
                    ProcessBits.allocate('safe');
                }
                try {
                    for (var _b = __values(_this.inputJax), _c = _b.next(); !_c.done; _c = _b.next()) {
                        var jax = _c.value;
                        if (jax.name.match(/MathML/)) {
                            jax.mathml.filterAttribute = _this.safe.mmlAttribute.bind(_this.safe);
                            jax.mathml.filterClassList = _this.safe.mmlClassList.bind(_this.safe);
                        }
                        else if (jax.name.match(/TeX/)) {
                            jax.postFilters.add(_this.sanitize.bind(jax), -5.5);
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
                return _this;
            }
            class_1.prototype.sanitize = function (data) {
                data.math.root = this.parseOptions.root;
                data.document.safe.sanitize(data.math, data.document);
            };
            return class_1;
        }(BaseDocument)),
        _a.OPTIONS = __assign(__assign({}, BaseDocument.OPTIONS), { safeOptions: __assign({}, safe_js_1.Safe.OPTIONS), SafeClass: safe_js_1.Safe }),
        _a;
}
exports.SafeMathDocumentMixin = SafeMathDocumentMixin;
function SafeHandler(handler) {
    handler.documentClass = SafeMathDocumentMixin(handler.documentClass);
    return handler;
}
exports.SafeHandler = SafeHandler;
//# sourceMappingURL=SafeHandler.js.map

/***/ }),

/***/ 35230:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
exports.SafeMethods = void 0;
var lengths_js_1 = __webpack_require__(56780);
exports.SafeMethods = {
    filterURL: function (safe, url) {
        var protocol = (url.match(/^\s*([a-z]+):/i) || [null, ''])[1].toLowerCase();
        var allow = safe.allow.URLs;
        return (allow === 'all' || (allow === 'safe' &&
            (safe.options.safeProtocols[protocol] || !protocol))) ? url : null;
    },
    filterClassList: function (safe, list) {
        var _this = this;
        var classes = list.trim().replace(/\s\s+/g, ' ').split(/ /);
        return classes.map(function (name) { return _this.filterClass(safe, name) || ''; }).join(' ').trim().replace(/\s\s+/g, '');
    },
    filterClass: function (safe, CLASS) {
        var allow = safe.allow.classes;
        return (allow === 'all' || (allow === 'safe' && CLASS.match(safe.options.classPattern))) ? CLASS : null;
    },
    filterID: function (safe, id) {
        var allow = safe.allow.cssIDs;
        return (allow === 'all' || (allow === 'safe' && id.match(safe.options.idPattern))) ? id : null;
    },
    filterStyles: function (safe, styles) {
        var e_1, _a, e_2, _b;
        if (safe.allow.styles === 'all')
            return styles;
        if (safe.allow.styles !== 'safe')
            return null;
        var adaptor = safe.adaptor;
        var options = safe.options;
        try {
            var div1 = adaptor.node('div', { style: styles });
            var div2 = adaptor.node('div');
            try {
                for (var _c = __values(Object.keys(options.safeStyles)), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var style = _d.value;
                    if (options.styleParts[style]) {
                        try {
                            for (var _e = (e_2 = void 0, __values(['Top', 'Right', 'Bottom', 'Left'])), _f = _e.next(); !_f.done; _f = _e.next()) {
                                var sufix = _f.value;
                                var name_1 = style + sufix;
                                var value = this.filterStyle(safe, name_1, div1);
                                if (value) {
                                    adaptor.setStyle(div2, name_1, value);
                                }
                            }
                        }
                        catch (e_2_1) { e_2 = { error: e_2_1 }; }
                        finally {
                            try {
                                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
                            }
                            finally { if (e_2) throw e_2.error; }
                        }
                    }
                    else {
                        var value = this.filterStyle(safe, style, div1);
                        if (value) {
                            adaptor.setStyle(div2, style, value);
                        }
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
            styles = adaptor.allStyles(div2);
        }
        catch (err) {
            styles = '';
        }
        return styles;
    },
    filterStyle: function (safe, style, div) {
        var value = safe.adaptor.getStyle(div, style);
        if (typeof value !== 'string' || value === '' || value.match(/^\s*calc/) ||
            (value.match(/javascript:/) && !safe.options.safeProtocols.javascript) ||
            (value.match(/data:/) && !safe.options.safeProtocols.data)) {
            return null;
        }
        var name = style.replace(/Top|Right|Left|Bottom/, '');
        if (!safe.options.safeStyles[style] && !safe.options.safeStyles[name]) {
            return null;
        }
        return this.filterStyleValue(safe, style, value, div);
    },
    filterStyleValue: function (safe, style, value, div) {
        var name = safe.options.styleLengths[style];
        if (!name) {
            return value;
        }
        if (typeof name !== 'string') {
            return this.filterStyleLength(safe, style, value);
        }
        var length = this.filterStyleLength(safe, name, safe.adaptor.getStyle(div, name));
        if (!length) {
            return null;
        }
        safe.adaptor.setStyle(div, name, length);
        return safe.adaptor.getStyle(div, style);
    },
    filterStyleLength: function (safe, style, value) {
        if (!value.match(/^(.+)(em|ex|ch|rem|px|mm|cm|in|pt|pc|%)$/))
            return null;
        var em = (0, lengths_js_1.length2em)(value, 1);
        var lengths = safe.options.styleLengths[style];
        var _a = __read((Array.isArray(lengths) ? lengths : [-safe.options.lengthMax, safe.options.lengthMax]), 2), m = _a[0], M = _a[1];
        return (m <= em && em <= M ? value : (em < m ? m : M).toFixed(3).replace(/\.?0+$/, '') + 'em');
    },
    filterFontSize: function (safe, size) {
        return this.filterStyleLength(safe, 'fontSize', size);
    },
    filterSizeMultiplier: function (safe, size) {
        var _a = __read(safe.options.scriptsizemultiplierRange || [-Infinity, Infinity], 2), m = _a[0], M = _a[1];
        return Math.min(M, Math.max(m, parseFloat(size))).toString();
    },
    filterScriptLevel: function (safe, level) {
        var _a = __read(safe.options.scriptlevelRange || [-Infinity, Infinity], 2), m = _a[0], M = _a[1];
        return Math.min(M, Math.max(m, parseInt(level))).toString();
    },
    filterData: function (safe, value, id) {
        return (id.match(safe.options.dataPattern) ? value : null);
    }
};
//# sourceMappingURL=SafeMethods.js.map

/***/ }),

/***/ 19290:
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
exports.Safe = void 0;
var Options_js_1 = __webpack_require__(4498);
var SafeMethods_js_1 = __webpack_require__(35230);
var Safe = (function () {
    function Safe(document, options) {
        this.filterAttributes = new Map([
            ['href', 'filterURL'],
            ['src', 'filterURL'],
            ['altimg', 'filterURL'],
            ['class', 'filterClassList'],
            ['style', 'filterStyles'],
            ['id', 'filterID'],
            ['fontsize', 'filterFontSize'],
            ['mathsize', 'filterFontSize'],
            ['scriptminsize', 'filterFontSize'],
            ['scriptsizemultiplier', 'filterSizeMultiplier'],
            ['scriptlevel', 'filterScriptLevel'],
            ['data-', 'filterData']
        ]);
        this.filterMethods = __assign({}, SafeMethods_js_1.SafeMethods);
        this.adaptor = document.adaptor;
        this.options = options;
        this.allow = this.options.allow;
    }
    Safe.prototype.sanitize = function (math, document) {
        try {
            math.root.walkTree(this.sanitizeNode.bind(this));
        }
        catch (err) {
            document.options.compileError(document, math, err);
        }
    };
    Safe.prototype.sanitizeNode = function (node) {
        var e_1, _a;
        var attributes = node.attributes.getAllAttributes();
        try {
            for (var _b = __values(Object.keys(attributes)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var id = _c.value;
                var method = this.filterAttributes.get(id);
                if (method) {
                    var value = this.filterMethods[method](this, attributes[id]);
                    if (value) {
                        if (value !== (typeof value === 'number' ? parseFloat(attributes[id]) : attributes[id])) {
                            attributes[id] = value;
                        }
                    }
                    else {
                        delete attributes[id];
                    }
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
    Safe.prototype.mmlAttribute = function (id, value) {
        if (id === 'class')
            return null;
        var method = this.filterAttributes.get(id);
        var filter = (method || (id.substr(0, 5) === 'data-' ? this.filterAttributes.get('data-') : null));
        if (!filter) {
            return value;
        }
        var result = this.filterMethods[filter](this, value, id);
        return (typeof result === 'number' || typeof result === 'boolean' ? String(result) : result);
    };
    Safe.prototype.mmlClassList = function (list) {
        var _this = this;
        return list.map(function (name) { return _this.filterMethods.filterClass(_this, name); })
            .filter(function (value) { return value !== null; });
    };
    Safe.OPTIONS = {
        allow: {
            URLs: 'safe',
            classes: 'safe',
            cssIDs: 'safe',
            styles: 'safe'
        },
        lengthMax: 3,
        scriptsizemultiplierRange: [.6, 1],
        scriptlevelRange: [-2, 2],
        classPattern: /^mjx-[-a-zA-Z0-9_.]+$/,
        idPattern: /^mjx-[-a-zA-Z0-9_.]+$/,
        dataPattern: /^data-mjx-/,
        safeProtocols: (0, Options_js_1.expandable)({
            http: true,
            https: true,
            file: true,
            javascript: false,
            data: false
        }),
        safeStyles: (0, Options_js_1.expandable)({
            color: true,
            backgroundColor: true,
            border: true,
            cursor: true,
            margin: true,
            padding: true,
            textShadow: true,
            fontFamily: true,
            fontSize: true,
            fontStyle: true,
            fontWeight: true,
            opacity: true,
            outline: true
        }),
        styleParts: (0, Options_js_1.expandable)({
            border: true,
            padding: true,
            margin: true,
            outline: true
        }),
        styleLengths: (0, Options_js_1.expandable)({
            borderTop: 'borderTopWidth',
            borderRight: 'borderRightWidth',
            borderBottom: 'borderBottomWidth',
            borderLeft: 'borderLeftWidth',
            paddingTop: true,
            paddingRight: true,
            paddingBottom: true,
            paddingLeft: true,
            marginTop: true,
            marginRight: true,
            marginBottom: true,
            marginLeft: true,
            outlineTop: true,
            outlineRight: true,
            outlineBottom: true,
            outlineLeft: true,
            fontSize: [.707, 1.44]
        })
    };
    return Safe;
}());
exports.Safe = Safe;
//# sourceMappingURL=safe.js.map

/***/ }),

/***/ 4498:
/***/ (function(__unused_webpack_module, exports) {


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
exports.lookup = exports.separateOptions = exports.selectOptionsFromKeys = exports.selectOptions = exports.userOptions = exports.defaultOptions = exports.insert = exports.copy = exports.keys = exports.makeArray = exports.expandable = exports.Expandable = exports.OPTIONS = exports.REMOVE = exports.APPEND = exports.isObject = void 0;
var OBJECT = {}.constructor;
function isObject(obj) {
    return typeof obj === 'object' && obj !== null &&
        (obj.constructor === OBJECT || obj.constructor === Expandable);
}
exports.isObject = isObject;
exports.APPEND = '[+]';
exports.REMOVE = '[-]';
exports.OPTIONS = {
    invalidOption: 'warn',
    optionError: function (message, _key) {
        if (exports.OPTIONS.invalidOption === 'fatal') {
            throw new Error(message);
        }
        console.warn('MathJax: ' + message);
    }
};
var Expandable = (function () {
    function Expandable() {
    }
    return Expandable;
}());
exports.Expandable = Expandable;
function expandable(def) {
    return Object.assign(Object.create(Expandable.prototype), def);
}
exports.expandable = expandable;
function makeArray(x) {
    return Array.isArray(x) ? x : [x];
}
exports.makeArray = makeArray;
function keys(def) {
    if (!def) {
        return [];
    }
    return Object.keys(def).concat(Object.getOwnPropertySymbols(def));
}
exports.keys = keys;
function copy(def) {
    var e_1, _a;
    var props = {};
    try {
        for (var _b = __values(keys(def)), _c = _b.next(); !_c.done; _c = _b.next()) {
            var key = _c.value;
            var prop = Object.getOwnPropertyDescriptor(def, key);
            var value = prop.value;
            if (Array.isArray(value)) {
                prop.value = insert([], value, false);
            }
            else if (isObject(value)) {
                prop.value = copy(value);
            }
            if (prop.enumerable) {
                props[key] = prop;
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
    return Object.defineProperties(def.constructor === Expandable ? expandable({}) : {}, props);
}
exports.copy = copy;
function insert(dst, src, warn) {
    var e_2, _a;
    if (warn === void 0) { warn = true; }
    var _loop_1 = function (key) {
        if (warn && dst[key] === undefined && dst.constructor !== Expandable) {
            if (typeof key === 'symbol') {
                key = key.toString();
            }
            exports.OPTIONS.optionError("Invalid option \"".concat(key, "\" (no default value)."), key);
            return "continue";
        }
        var sval = src[key], dval = dst[key];
        if (isObject(sval) && dval !== null &&
            (typeof dval === 'object' || typeof dval === 'function')) {
            var ids = keys(sval);
            if (Array.isArray(dval) &&
                ((ids.length === 1 && (ids[0] === exports.APPEND || ids[0] === exports.REMOVE) && Array.isArray(sval[ids[0]])) ||
                    (ids.length === 2 && ids.sort().join(',') === exports.APPEND + ',' + exports.REMOVE &&
                        Array.isArray(sval[exports.APPEND]) && Array.isArray(sval[exports.REMOVE])))) {
                if (sval[exports.REMOVE]) {
                    dval = dst[key] = dval.filter(function (x) { return sval[exports.REMOVE].indexOf(x) < 0; });
                }
                if (sval[exports.APPEND]) {
                    dst[key] = __spreadArray(__spreadArray([], __read(dval), false), __read(sval[exports.APPEND]), false);
                }
            }
            else {
                insert(dval, sval, warn);
            }
        }
        else if (Array.isArray(sval)) {
            dst[key] = [];
            insert(dst[key], sval, false);
        }
        else if (isObject(sval)) {
            dst[key] = copy(sval);
        }
        else {
            dst[key] = sval;
        }
    };
    try {
        for (var _b = __values(keys(src)), _c = _b.next(); !_c.done; _c = _b.next()) {
            var key = _c.value;
            _loop_1(key);
        }
    }
    catch (e_2_1) { e_2 = { error: e_2_1 }; }
    finally {
        try {
            if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
        }
        finally { if (e_2) throw e_2.error; }
    }
    return dst;
}
exports.insert = insert;
function defaultOptions(options) {
    var defs = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        defs[_i - 1] = arguments[_i];
    }
    defs.forEach(function (def) { return insert(options, def, false); });
    return options;
}
exports.defaultOptions = defaultOptions;
function userOptions(options) {
    var defs = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        defs[_i - 1] = arguments[_i];
    }
    defs.forEach(function (def) { return insert(options, def, true); });
    return options;
}
exports.userOptions = userOptions;
function selectOptions(options) {
    var e_3, _a;
    var keys = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        keys[_i - 1] = arguments[_i];
    }
    var subset = {};
    try {
        for (var keys_1 = __values(keys), keys_1_1 = keys_1.next(); !keys_1_1.done; keys_1_1 = keys_1.next()) {
            var key = keys_1_1.value;
            if (options.hasOwnProperty(key)) {
                subset[key] = options[key];
            }
        }
    }
    catch (e_3_1) { e_3 = { error: e_3_1 }; }
    finally {
        try {
            if (keys_1_1 && !keys_1_1.done && (_a = keys_1.return)) _a.call(keys_1);
        }
        finally { if (e_3) throw e_3.error; }
    }
    return subset;
}
exports.selectOptions = selectOptions;
function selectOptionsFromKeys(options, object) {
    return selectOptions.apply(void 0, __spreadArray([options], __read(Object.keys(object)), false));
}
exports.selectOptionsFromKeys = selectOptionsFromKeys;
function separateOptions(options) {
    var e_4, _a, e_5, _b;
    var objects = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        objects[_i - 1] = arguments[_i];
    }
    var results = [];
    try {
        for (var objects_1 = __values(objects), objects_1_1 = objects_1.next(); !objects_1_1.done; objects_1_1 = objects_1.next()) {
            var object = objects_1_1.value;
            var exists = {}, missing = {};
            try {
                for (var _c = (e_5 = void 0, __values(Object.keys(options || {}))), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var key = _d.value;
                    (object[key] === undefined ? missing : exists)[key] = options[key];
                }
            }
            catch (e_5_1) { e_5 = { error: e_5_1 }; }
            finally {
                try {
                    if (_d && !_d.done && (_b = _c.return)) _b.call(_c);
                }
                finally { if (e_5) throw e_5.error; }
            }
            results.push(exists);
            options = missing;
        }
    }
    catch (e_4_1) { e_4 = { error: e_4_1 }; }
    finally {
        try {
            if (objects_1_1 && !objects_1_1.done && (_a = objects_1.return)) _a.call(objects_1);
        }
        finally { if (e_4) throw e_4.error; }
    }
    results.unshift(options);
    return results;
}
exports.separateOptions = separateOptions;
function lookup(name, lookup, def) {
    if (def === void 0) { def = null; }
    return (lookup.hasOwnProperty(name) ? lookup[name] : def);
}
exports.lookup = lookup;
//# sourceMappingURL=Options.js.map

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
//# sourceMappingURL=8285.8bade38c361d9af60b43.js.map?v=8bade38c361d9af60b43