"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7471],{

/***/ 35663:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractHandler = void 0;
var MathDocument_js_1 = __webpack_require__(13139);
var DefaultMathDocument = (function (_super) {
    __extends(DefaultMathDocument, _super);
    function DefaultMathDocument() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return DefaultMathDocument;
}(MathDocument_js_1.AbstractMathDocument));
var AbstractHandler = (function () {
    function AbstractHandler(adaptor, priority) {
        if (priority === void 0) { priority = 5; }
        this.documentClass = DefaultMathDocument;
        this.adaptor = adaptor;
        this.priority = priority;
    }
    Object.defineProperty(AbstractHandler.prototype, "name", {
        get: function () {
            return this.constructor.NAME;
        },
        enumerable: false,
        configurable: true
    });
    AbstractHandler.prototype.handlesDocument = function (_document) {
        return false;
    };
    AbstractHandler.prototype.create = function (document, options) {
        return new this.documentClass(document, this.adaptor, options);
    };
    AbstractHandler.NAME = 'generic';
    return AbstractHandler;
}());
exports.AbstractHandler = AbstractHandler;
//# sourceMappingURL=Handler.js.map

/***/ }),

/***/ 76771:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractInputJax = void 0;
var Options_js_1 = __webpack_require__(4498);
var FunctionList_js_1 = __webpack_require__(18341);
var AbstractInputJax = (function () {
    function AbstractInputJax(options) {
        if (options === void 0) { options = {}; }
        this.adaptor = null;
        this.mmlFactory = null;
        var CLASS = this.constructor;
        this.options = (0, Options_js_1.userOptions)((0, Options_js_1.defaultOptions)({}, CLASS.OPTIONS), options);
        this.preFilters = new FunctionList_js_1.FunctionList();
        this.postFilters = new FunctionList_js_1.FunctionList();
    }
    Object.defineProperty(AbstractInputJax.prototype, "name", {
        get: function () {
            return this.constructor.NAME;
        },
        enumerable: false,
        configurable: true
    });
    AbstractInputJax.prototype.setAdaptor = function (adaptor) {
        this.adaptor = adaptor;
    };
    AbstractInputJax.prototype.setMmlFactory = function (mmlFactory) {
        this.mmlFactory = mmlFactory;
    };
    AbstractInputJax.prototype.initialize = function () {
    };
    AbstractInputJax.prototype.reset = function () {
        var _args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            _args[_i] = arguments[_i];
        }
    };
    Object.defineProperty(AbstractInputJax.prototype, "processStrings", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    AbstractInputJax.prototype.findMath = function (_node, _options) {
        return [];
    };
    AbstractInputJax.prototype.executeFilters = function (filters, math, document, data) {
        var args = { math: math, document: document, data: data };
        filters.execute(args);
        return args.data;
    };
    AbstractInputJax.NAME = 'generic';
    AbstractInputJax.OPTIONS = {};
    return AbstractInputJax;
}());
exports.AbstractInputJax = AbstractInputJax;
//# sourceMappingURL=InputJax.js.map

/***/ }),

/***/ 13139:
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
exports.AbstractMathDocument = exports.resetAllOptions = exports.resetOptions = exports.RenderList = void 0;
var Options_js_1 = __webpack_require__(4498);
var InputJax_js_1 = __webpack_require__(76771);
var OutputJax_js_1 = __webpack_require__(44005);
var MathList_js_1 = __webpack_require__(80426);
var MathItem_js_1 = __webpack_require__(21605);
var MmlFactory_js_1 = __webpack_require__(72666);
var BitField_js_1 = __webpack_require__(62844);
var PrioritizedList_js_1 = __webpack_require__(98721);
var RenderList = (function (_super) {
    __extends(RenderList, _super);
    function RenderList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    RenderList.create = function (actions) {
        var e_1, _a;
        var list = new this();
        try {
            for (var _b = __values(Object.keys(actions)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var id = _c.value;
                var _d = __read(this.action(id, actions[id]), 2), action = _d[0], priority = _d[1];
                if (priority) {
                    list.add(action, priority);
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
        return list;
    };
    RenderList.action = function (id, action) {
        var _a, _b, _c, _d;
        var renderDoc, renderMath;
        var convert = true;
        var priority = action[0];
        if (action.length === 1 || typeof action[1] === 'boolean') {
            action.length === 2 && (convert = action[1]);
            _a = __read(this.methodActions(id), 2), renderDoc = _a[0], renderMath = _a[1];
        }
        else if (typeof action[1] === 'string') {
            if (typeof action[2] === 'string') {
                action.length === 4 && (convert = action[3]);
                var _e = __read(action.slice(1), 2), method1 = _e[0], method2 = _e[1];
                _b = __read(this.methodActions(method1, method2), 2), renderDoc = _b[0], renderMath = _b[1];
            }
            else {
                action.length === 3 && (convert = action[2]);
                _c = __read(this.methodActions(action[1]), 2), renderDoc = _c[0], renderMath = _c[1];
            }
        }
        else {
            action.length === 4 && (convert = action[3]);
            _d = __read(action.slice(1), 2), renderDoc = _d[0], renderMath = _d[1];
        }
        return [{ id: id, renderDoc: renderDoc, renderMath: renderMath, convert: convert }, priority];
    };
    RenderList.methodActions = function (method1, method2) {
        if (method2 === void 0) { method2 = method1; }
        return [
            function (document) { method1 && document[method1](); return false; },
            function (math, document) { method2 && math[method2](document); return false; }
        ];
    };
    RenderList.prototype.renderDoc = function (document, start) {
        var e_2, _a;
        if (start === void 0) { start = MathItem_js_1.STATE.UNPROCESSED; }
        try {
            for (var _b = __values(this.items), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                if (item.priority >= start) {
                    if (item.item.renderDoc(document))
                        return;
                }
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_2) throw e_2.error; }
        }
    };
    RenderList.prototype.renderMath = function (math, document, start) {
        var e_3, _a;
        if (start === void 0) { start = MathItem_js_1.STATE.UNPROCESSED; }
        try {
            for (var _b = __values(this.items), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                if (item.priority >= start) {
                    if (item.item.renderMath(math, document))
                        return;
                }
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_3) throw e_3.error; }
        }
    };
    RenderList.prototype.renderConvert = function (math, document, end) {
        var e_4, _a;
        if (end === void 0) { end = MathItem_js_1.STATE.LAST; }
        try {
            for (var _b = __values(this.items), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                if (item.priority > end)
                    return;
                if (item.item.convert) {
                    if (item.item.renderMath(math, document))
                        return;
                }
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_4) throw e_4.error; }
        }
    };
    RenderList.prototype.findID = function (id) {
        var e_5, _a;
        try {
            for (var _b = __values(this.items), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                if (item.item.id === id) {
                    return item.item;
                }
            }
        }
        catch (e_5_1) { e_5 = { error: e_5_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_5) throw e_5.error; }
        }
        return null;
    };
    return RenderList;
}(PrioritizedList_js_1.PrioritizedList));
exports.RenderList = RenderList;
exports.resetOptions = {
    all: false,
    processed: false,
    inputJax: null,
    outputJax: null
};
exports.resetAllOptions = {
    all: true,
    processed: true,
    inputJax: [],
    outputJax: []
};
var DefaultInputJax = (function (_super) {
    __extends(DefaultInputJax, _super);
    function DefaultInputJax() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    DefaultInputJax.prototype.compile = function (_math) {
        return null;
    };
    return DefaultInputJax;
}(InputJax_js_1.AbstractInputJax));
var DefaultOutputJax = (function (_super) {
    __extends(DefaultOutputJax, _super);
    function DefaultOutputJax() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    DefaultOutputJax.prototype.typeset = function (_math, _document) {
        if (_document === void 0) { _document = null; }
        return null;
    };
    DefaultOutputJax.prototype.escaped = function (_math, _document) {
        return null;
    };
    return DefaultOutputJax;
}(OutputJax_js_1.AbstractOutputJax));
var DefaultMathList = (function (_super) {
    __extends(DefaultMathList, _super);
    function DefaultMathList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return DefaultMathList;
}(MathList_js_1.AbstractMathList));
var DefaultMathItem = (function (_super) {
    __extends(DefaultMathItem, _super);
    function DefaultMathItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return DefaultMathItem;
}(MathItem_js_1.AbstractMathItem));
var AbstractMathDocument = (function () {
    function AbstractMathDocument(document, adaptor, options) {
        var _this = this;
        var CLASS = this.constructor;
        this.document = document;
        this.options = (0, Options_js_1.userOptions)((0, Options_js_1.defaultOptions)({}, CLASS.OPTIONS), options);
        this.math = new (this.options['MathList'] || DefaultMathList)();
        this.renderActions = RenderList.create(this.options['renderActions']);
        this.processed = new AbstractMathDocument.ProcessBits();
        this.outputJax = this.options['OutputJax'] || new DefaultOutputJax();
        var inputJax = this.options['InputJax'] || [new DefaultInputJax()];
        if (!Array.isArray(inputJax)) {
            inputJax = [inputJax];
        }
        this.inputJax = inputJax;
        this.adaptor = adaptor;
        this.outputJax.setAdaptor(adaptor);
        this.inputJax.map(function (jax) { return jax.setAdaptor(adaptor); });
        this.mmlFactory = this.options['MmlFactory'] || new MmlFactory_js_1.MmlFactory();
        this.inputJax.map(function (jax) { return jax.setMmlFactory(_this.mmlFactory); });
        this.outputJax.initialize();
        this.inputJax.map(function (jax) { return jax.initialize(); });
    }
    Object.defineProperty(AbstractMathDocument.prototype, "kind", {
        get: function () {
            return this.constructor.KIND;
        },
        enumerable: false,
        configurable: true
    });
    AbstractMathDocument.prototype.addRenderAction = function (id) {
        var action = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            action[_i - 1] = arguments[_i];
        }
        var _a = __read(RenderList.action(id, action), 2), fn = _a[0], p = _a[1];
        this.renderActions.add(fn, p);
    };
    AbstractMathDocument.prototype.removeRenderAction = function (id) {
        var action = this.renderActions.findID(id);
        if (action) {
            this.renderActions.remove(action);
        }
    };
    AbstractMathDocument.prototype.render = function () {
        this.renderActions.renderDoc(this);
        return this;
    };
    AbstractMathDocument.prototype.rerender = function (start) {
        if (start === void 0) { start = MathItem_js_1.STATE.RERENDER; }
        this.state(start - 1);
        this.render();
        return this;
    };
    AbstractMathDocument.prototype.convert = function (math, options) {
        if (options === void 0) { options = {}; }
        var _a = (0, Options_js_1.userOptions)({
            format: this.inputJax[0].name, display: true, end: MathItem_js_1.STATE.LAST,
            em: 16, ex: 8, containerWidth: null, lineWidth: 1000000, scale: 1, family: ''
        }, options), format = _a.format, display = _a.display, end = _a.end, ex = _a.ex, em = _a.em, containerWidth = _a.containerWidth, lineWidth = _a.lineWidth, scale = _a.scale, family = _a.family;
        if (containerWidth === null) {
            containerWidth = 80 * ex;
        }
        var jax = this.inputJax.reduce(function (jax, ijax) { return (ijax.name === format ? ijax : jax); }, null);
        var mitem = new this.options.MathItem(math, jax, display);
        mitem.start.node = this.adaptor.body(this.document);
        mitem.setMetrics(em, ex, containerWidth, lineWidth, scale);
        if (this.outputJax.options.mtextInheritFont) {
            mitem.outputData.mtextFamily = family;
        }
        if (this.outputJax.options.merrorInheritFont) {
            mitem.outputData.merrorFamily = family;
        }
        mitem.convert(this, end);
        return (mitem.typesetRoot || mitem.root);
    };
    AbstractMathDocument.prototype.findMath = function (_options) {
        if (_options === void 0) { _options = null; }
        this.processed.set('findMath');
        return this;
    };
    AbstractMathDocument.prototype.compile = function () {
        var e_6, _a, e_7, _b;
        if (!this.processed.isSet('compile')) {
            var recompile = [];
            try {
                for (var _c = __values(this.math), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var math = _d.value;
                    this.compileMath(math);
                    if (math.inputData.recompile !== undefined) {
                        recompile.push(math);
                    }
                }
            }
            catch (e_6_1) { e_6 = { error: e_6_1 }; }
            finally {
                try {
                    if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
                }
                finally { if (e_6) throw e_6.error; }
            }
            try {
                for (var recompile_1 = __values(recompile), recompile_1_1 = recompile_1.next(); !recompile_1_1.done; recompile_1_1 = recompile_1.next()) {
                    var math = recompile_1_1.value;
                    var data = math.inputData.recompile;
                    math.state(data.state);
                    math.inputData.recompile = data;
                    this.compileMath(math);
                }
            }
            catch (e_7_1) { e_7 = { error: e_7_1 }; }
            finally {
                try {
                    if (recompile_1_1 && !recompile_1_1.done && (_b = recompile_1.return)) _b.call(recompile_1);
                }
                finally { if (e_7) throw e_7.error; }
            }
            this.processed.set('compile');
        }
        return this;
    };
    AbstractMathDocument.prototype.compileMath = function (math) {
        try {
            math.compile(this);
        }
        catch (err) {
            if (err.retry || err.restart) {
                throw err;
            }
            this.options['compileError'](this, math, err);
            math.inputData['error'] = err;
        }
    };
    AbstractMathDocument.prototype.compileError = function (math, err) {
        math.root = this.mmlFactory.create('math', null, [
            this.mmlFactory.create('merror', { 'data-mjx-error': err.message, title: err.message }, [
                this.mmlFactory.create('mtext', null, [
                    this.mmlFactory.create('text').setText('Math input error')
                ])
            ])
        ]);
        if (math.display) {
            math.root.attributes.set('display', 'block');
        }
        math.inputData.error = err.message;
    };
    AbstractMathDocument.prototype.typeset = function () {
        var e_8, _a;
        if (!this.processed.isSet('typeset')) {
            try {
                for (var _b = __values(this.math), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var math = _c.value;
                    try {
                        math.typeset(this);
                    }
                    catch (err) {
                        if (err.retry || err.restart) {
                            throw err;
                        }
                        this.options['typesetError'](this, math, err);
                        math.outputData['error'] = err;
                    }
                }
            }
            catch (e_8_1) { e_8 = { error: e_8_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_8) throw e_8.error; }
            }
            this.processed.set('typeset');
        }
        return this;
    };
    AbstractMathDocument.prototype.typesetError = function (math, err) {
        math.typesetRoot = this.adaptor.node('mjx-container', {
            class: 'MathJax mjx-output-error',
            jax: this.outputJax.name,
        }, [
            this.adaptor.node('span', {
                'data-mjx-error': err.message,
                title: err.message,
                style: {
                    color: 'red',
                    'background-color': 'yellow',
                    'line-height': 'normal'
                }
            }, [
                this.adaptor.text('Math output error')
            ])
        ]);
        if (math.display) {
            this.adaptor.setAttributes(math.typesetRoot, {
                style: {
                    display: 'block',
                    margin: '1em 0',
                    'text-align': 'center'
                }
            });
        }
        math.outputData.error = err.message;
    };
    AbstractMathDocument.prototype.getMetrics = function () {
        if (!this.processed.isSet('getMetrics')) {
            this.outputJax.getMetrics(this);
            this.processed.set('getMetrics');
        }
        return this;
    };
    AbstractMathDocument.prototype.updateDocument = function () {
        var e_9, _a;
        if (!this.processed.isSet('updateDocument')) {
            try {
                for (var _b = __values(this.math.reversed()), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var math = _c.value;
                    math.updateDocument(this);
                }
            }
            catch (e_9_1) { e_9 = { error: e_9_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_9) throw e_9.error; }
            }
            this.processed.set('updateDocument');
        }
        return this;
    };
    AbstractMathDocument.prototype.removeFromDocument = function (_restore) {
        if (_restore === void 0) { _restore = false; }
        return this;
    };
    AbstractMathDocument.prototype.state = function (state, restore) {
        var e_10, _a;
        if (restore === void 0) { restore = false; }
        try {
            for (var _b = __values(this.math), _c = _b.next(); !_c.done; _c = _b.next()) {
                var math = _c.value;
                math.state(state, restore);
            }
        }
        catch (e_10_1) { e_10 = { error: e_10_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_10) throw e_10.error; }
        }
        if (state < MathItem_js_1.STATE.INSERTED) {
            this.processed.clear('updateDocument');
        }
        if (state < MathItem_js_1.STATE.TYPESET) {
            this.processed.clear('typeset');
            this.processed.clear('getMetrics');
        }
        if (state < MathItem_js_1.STATE.COMPILED) {
            this.processed.clear('compile');
        }
        return this;
    };
    AbstractMathDocument.prototype.reset = function (options) {
        var _a;
        if (options === void 0) { options = { processed: true }; }
        options = (0, Options_js_1.userOptions)(Object.assign({}, exports.resetOptions), options);
        options.all && Object.assign(options, exports.resetAllOptions);
        options.processed && this.processed.reset();
        options.inputJax && this.inputJax.forEach(function (jax) { return jax.reset.apply(jax, __spreadArray([], __read(options.inputJax), false)); });
        options.outputJax && (_a = this.outputJax).reset.apply(_a, __spreadArray([], __read(options.outputJax), false));
        return this;
    };
    AbstractMathDocument.prototype.clear = function () {
        this.reset();
        this.math.clear();
        return this;
    };
    AbstractMathDocument.prototype.concat = function (list) {
        this.math.merge(list);
        return this;
    };
    AbstractMathDocument.prototype.clearMathItemsWithin = function (containers) {
        var _a;
        var items = this.getMathItemsWithin(containers);
        (_a = this.math).remove.apply(_a, __spreadArray([], __read(items), false));
        return items;
    };
    AbstractMathDocument.prototype.getMathItemsWithin = function (elements) {
        var e_11, _a, e_12, _b;
        if (!Array.isArray(elements)) {
            elements = [elements];
        }
        var adaptor = this.adaptor;
        var items = [];
        var containers = adaptor.getElements(elements, this.document);
        try {
            ITEMS: for (var _c = __values(this.math), _d = _c.next(); !_d.done; _d = _c.next()) {
                var item = _d.value;
                try {
                    for (var containers_1 = (e_12 = void 0, __values(containers)), containers_1_1 = containers_1.next(); !containers_1_1.done; containers_1_1 = containers_1.next()) {
                        var container = containers_1_1.value;
                        if (item.start.node && adaptor.contains(container, item.start.node)) {
                            items.push(item);
                            continue ITEMS;
                        }
                    }
                }
                catch (e_12_1) { e_12 = { error: e_12_1 }; }
                finally {
                    try {
                        if (containers_1_1 && !containers_1_1.done && (_b = containers_1.return)) _b.call(containers_1);
                    }
                    finally { if (e_12) throw e_12.error; }
                }
            }
        }
        catch (e_11_1) { e_11 = { error: e_11_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_11) throw e_11.error; }
        }
        return items;
    };
    AbstractMathDocument.KIND = 'MathDocument';
    AbstractMathDocument.OPTIONS = {
        OutputJax: null,
        InputJax: null,
        MmlFactory: null,
        MathList: DefaultMathList,
        MathItem: DefaultMathItem,
        compileError: function (doc, math, err) {
            doc.compileError(math, err);
        },
        typesetError: function (doc, math, err) {
            doc.typesetError(math, err);
        },
        renderActions: (0, Options_js_1.expandable)({
            find: [MathItem_js_1.STATE.FINDMATH, 'findMath', '', false],
            compile: [MathItem_js_1.STATE.COMPILED],
            metrics: [MathItem_js_1.STATE.METRICS, 'getMetrics', '', false],
            typeset: [MathItem_js_1.STATE.TYPESET],
            update: [MathItem_js_1.STATE.INSERTED, 'updateDocument', false]
        })
    };
    AbstractMathDocument.ProcessBits = (0, BitField_js_1.BitFieldClass)('findMath', 'compile', 'getMetrics', 'typeset', 'updateDocument');
    return AbstractMathDocument;
}());
exports.AbstractMathDocument = AbstractMathDocument;
//# sourceMappingURL=MathDocument.js.map

/***/ }),

/***/ 80426:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractMathList = void 0;
var LinkedList_js_1 = __webpack_require__(6811);
var AbstractMathList = (function (_super) {
    __extends(AbstractMathList, _super);
    function AbstractMathList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    AbstractMathList.prototype.isBefore = function (a, b) {
        return (a.start.i < b.start.i || (a.start.i === b.start.i && a.start.n < b.start.n));
    };
    return AbstractMathList;
}(LinkedList_js_1.LinkedList));
exports.AbstractMathList = AbstractMathList;
//# sourceMappingURL=MathList.js.map

/***/ }),

/***/ 44005:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractOutputJax = void 0;
var Options_js_1 = __webpack_require__(4498);
var FunctionList_js_1 = __webpack_require__(18341);
var AbstractOutputJax = (function () {
    function AbstractOutputJax(options) {
        if (options === void 0) { options = {}; }
        this.adaptor = null;
        var CLASS = this.constructor;
        this.options = (0, Options_js_1.userOptions)((0, Options_js_1.defaultOptions)({}, CLASS.OPTIONS), options);
        this.postFilters = new FunctionList_js_1.FunctionList();
    }
    Object.defineProperty(AbstractOutputJax.prototype, "name", {
        get: function () {
            return this.constructor.NAME;
        },
        enumerable: false,
        configurable: true
    });
    AbstractOutputJax.prototype.setAdaptor = function (adaptor) {
        this.adaptor = adaptor;
    };
    AbstractOutputJax.prototype.initialize = function () {
    };
    AbstractOutputJax.prototype.reset = function () {
        var _args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            _args[_i] = arguments[_i];
        }
    };
    AbstractOutputJax.prototype.getMetrics = function (_document) {
    };
    AbstractOutputJax.prototype.styleSheet = function (_document) {
        return null;
    };
    AbstractOutputJax.prototype.pageElements = function (_document) {
        return null;
    };
    AbstractOutputJax.prototype.executeFilters = function (filters, math, document, data) {
        var args = { math: math, document: document, data: data };
        filters.execute(args);
        return args.data;
    };
    AbstractOutputJax.NAME = 'generic';
    AbstractOutputJax.OPTIONS = {};
    return AbstractOutputJax;
}());
exports.AbstractOutputJax = AbstractOutputJax;
//# sourceMappingURL=OutputJax.js.map

/***/ }),

/***/ 99331:
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
exports.HTMLDocument = void 0;
var MathDocument_js_1 = __webpack_require__(13139);
var Options_js_1 = __webpack_require__(4498);
var HTMLMathItem_js_1 = __webpack_require__(92870);
var HTMLMathList_js_1 = __webpack_require__(65337);
var HTMLDomStrings_js_1 = __webpack_require__(43751);
var MathItem_js_1 = __webpack_require__(21605);
var HTMLDocument = (function (_super) {
    __extends(HTMLDocument, _super);
    function HTMLDocument(document, adaptor, options) {
        var _this = this;
        var _a = __read((0, Options_js_1.separateOptions)(options, HTMLDomStrings_js_1.HTMLDomStrings.OPTIONS), 2), html = _a[0], dom = _a[1];
        _this = _super.call(this, document, adaptor, html) || this;
        _this.domStrings = _this.options['DomStrings'] || new HTMLDomStrings_js_1.HTMLDomStrings(dom);
        _this.domStrings.adaptor = adaptor;
        _this.styles = [];
        return _this;
    }
    HTMLDocument.prototype.findPosition = function (N, index, delim, nodes) {
        var e_1, _a;
        var adaptor = this.adaptor;
        try {
            for (var _b = __values(nodes[N]), _c = _b.next(); !_c.done; _c = _b.next()) {
                var list = _c.value;
                var _d = __read(list, 2), node = _d[0], n = _d[1];
                if (index <= n && adaptor.kind(node) === '#text') {
                    return { node: node, n: Math.max(index, 0), delim: delim };
                }
                index -= n;
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return { node: null, n: 0, delim: delim };
    };
    HTMLDocument.prototype.mathItem = function (item, jax, nodes) {
        var math = item.math;
        var start = this.findPosition(item.n, item.start.n, item.open, nodes);
        var end = this.findPosition(item.n, item.end.n, item.close, nodes);
        return new this.options.MathItem(math, jax, item.display, start, end);
    };
    HTMLDocument.prototype.findMath = function (options) {
        var e_2, _a, e_3, _b, _c, e_4, _d, e_5, _e;
        if (!this.processed.isSet('findMath')) {
            this.adaptor.document = this.document;
            options = (0, Options_js_1.userOptions)({ elements: this.options.elements || [this.adaptor.body(this.document)] }, options);
            try {
                for (var _f = __values(this.adaptor.getElements(options['elements'], this.document)), _g = _f.next(); !_g.done; _g = _f.next()) {
                    var container = _g.value;
                    var _h = __read([null, null], 2), strings = _h[0], nodes = _h[1];
                    try {
                        for (var _j = (e_3 = void 0, __values(this.inputJax)), _k = _j.next(); !_k.done; _k = _j.next()) {
                            var jax = _k.value;
                            var list = new (this.options['MathList'])();
                            if (jax.processStrings) {
                                if (strings === null) {
                                    _c = __read(this.domStrings.find(container), 2), strings = _c[0], nodes = _c[1];
                                }
                                try {
                                    for (var _l = (e_4 = void 0, __values(jax.findMath(strings))), _m = _l.next(); !_m.done; _m = _l.next()) {
                                        var math = _m.value;
                                        list.push(this.mathItem(math, jax, nodes));
                                    }
                                }
                                catch (e_4_1) { e_4 = { error: e_4_1 }; }
                                finally {
                                    try {
                                        if (_m && !_m.done && (_d = _l.return)) _d.call(_l);
                                    }
                                    finally { if (e_4) throw e_4.error; }
                                }
                            }
                            else {
                                try {
                                    for (var _o = (e_5 = void 0, __values(jax.findMath(container))), _p = _o.next(); !_p.done; _p = _o.next()) {
                                        var math = _p.value;
                                        var item = new this.options.MathItem(math.math, jax, math.display, math.start, math.end);
                                        list.push(item);
                                    }
                                }
                                catch (e_5_1) { e_5 = { error: e_5_1 }; }
                                finally {
                                    try {
                                        if (_p && !_p.done && (_e = _o.return)) _e.call(_o);
                                    }
                                    finally { if (e_5) throw e_5.error; }
                                }
                            }
                            this.math.merge(list);
                        }
                    }
                    catch (e_3_1) { e_3 = { error: e_3_1 }; }
                    finally {
                        try {
                            if (_k && !_k.done && (_b = _j.return)) _b.call(_j);
                        }
                        finally { if (e_3) throw e_3.error; }
                    }
                }
            }
            catch (e_2_1) { e_2 = { error: e_2_1 }; }
            finally {
                try {
                    if (_g && !_g.done && (_a = _f.return)) _a.call(_f);
                }
                finally { if (e_2) throw e_2.error; }
            }
            this.processed.set('findMath');
        }
        return this;
    };
    HTMLDocument.prototype.updateDocument = function () {
        if (!this.processed.isSet('updateDocument')) {
            this.addPageElements();
            this.addStyleSheet();
            _super.prototype.updateDocument.call(this);
            this.processed.set('updateDocument');
        }
        return this;
    };
    HTMLDocument.prototype.addPageElements = function () {
        var body = this.adaptor.body(this.document);
        var node = this.documentPageElements();
        if (node) {
            this.adaptor.append(body, node);
        }
    };
    HTMLDocument.prototype.addStyleSheet = function () {
        var sheet = this.documentStyleSheet();
        var adaptor = this.adaptor;
        if (sheet && !adaptor.parent(sheet)) {
            var head = adaptor.head(this.document);
            var styles = this.findSheet(head, adaptor.getAttribute(sheet, 'id'));
            if (styles) {
                adaptor.replace(sheet, styles);
            }
            else {
                adaptor.append(head, sheet);
            }
        }
    };
    HTMLDocument.prototype.findSheet = function (head, id) {
        var e_6, _a;
        if (id) {
            try {
                for (var _b = __values(this.adaptor.tags(head, 'style')), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var sheet = _c.value;
                    if (this.adaptor.getAttribute(sheet, 'id') === id) {
                        return sheet;
                    }
                }
            }
            catch (e_6_1) { e_6 = { error: e_6_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_6) throw e_6.error; }
            }
        }
        return null;
    };
    HTMLDocument.prototype.removeFromDocument = function (restore) {
        var e_7, _a;
        if (restore === void 0) { restore = false; }
        if (this.processed.isSet('updateDocument')) {
            try {
                for (var _b = __values(this.math), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var math = _c.value;
                    if (math.state() >= MathItem_js_1.STATE.INSERTED) {
                        math.state(MathItem_js_1.STATE.TYPESET, restore);
                    }
                }
            }
            catch (e_7_1) { e_7 = { error: e_7_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_7) throw e_7.error; }
            }
        }
        this.processed.clear('updateDocument');
        return this;
    };
    HTMLDocument.prototype.documentStyleSheet = function () {
        return this.outputJax.styleSheet(this);
    };
    HTMLDocument.prototype.documentPageElements = function () {
        return this.outputJax.pageElements(this);
    };
    HTMLDocument.prototype.addStyles = function (styles) {
        this.styles.push(styles);
    };
    HTMLDocument.prototype.getStyles = function () {
        return this.styles;
    };
    HTMLDocument.KIND = 'HTML';
    HTMLDocument.OPTIONS = __assign(__assign({}, MathDocument_js_1.AbstractMathDocument.OPTIONS), { renderActions: (0, Options_js_1.expandable)(__assign(__assign({}, MathDocument_js_1.AbstractMathDocument.OPTIONS.renderActions), { styles: [MathItem_js_1.STATE.INSERTED + 1, '', 'updateStyleSheet', false] })), MathList: HTMLMathList_js_1.HTMLMathList, MathItem: HTMLMathItem_js_1.HTMLMathItem, DomStrings: null });
    return HTMLDocument;
}(MathDocument_js_1.AbstractMathDocument));
exports.HTMLDocument = HTMLDocument;
//# sourceMappingURL=HTMLDocument.js.map

/***/ }),

/***/ 43751:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
exports.HTMLDomStrings = void 0;
var Options_js_1 = __webpack_require__(4498);
var HTMLDomStrings = (function () {
    function HTMLDomStrings(options) {
        if (options === void 0) { options = null; }
        var CLASS = this.constructor;
        this.options = (0, Options_js_1.userOptions)((0, Options_js_1.defaultOptions)({}, CLASS.OPTIONS), options);
        this.init();
        this.getPatterns();
    }
    HTMLDomStrings.prototype.init = function () {
        this.strings = [];
        this.string = '';
        this.snodes = [];
        this.nodes = [];
        this.stack = [];
    };
    HTMLDomStrings.prototype.getPatterns = function () {
        var skip = (0, Options_js_1.makeArray)(this.options['skipHtmlTags']);
        var ignore = (0, Options_js_1.makeArray)(this.options['ignoreHtmlClass']);
        var process = (0, Options_js_1.makeArray)(this.options['processHtmlClass']);
        this.skipHtmlTags = new RegExp('^(?:' + skip.join('|') + ')$', 'i');
        this.ignoreHtmlClass = new RegExp('(?:^| )(?:' + ignore.join('|') + ')(?: |$)');
        this.processHtmlClass = new RegExp('(?:^| )(?:' + process + ')(?: |$)');
    };
    HTMLDomStrings.prototype.pushString = function () {
        if (this.string.match(/\S/)) {
            this.strings.push(this.string);
            this.nodes.push(this.snodes);
        }
        this.string = '';
        this.snodes = [];
    };
    HTMLDomStrings.prototype.extendString = function (node, text) {
        this.snodes.push([node, text.length]);
        this.string += text;
    };
    HTMLDomStrings.prototype.handleText = function (node, ignore) {
        if (!ignore) {
            this.extendString(node, this.adaptor.value(node));
        }
        return this.adaptor.next(node);
    };
    HTMLDomStrings.prototype.handleTag = function (node, ignore) {
        if (!ignore) {
            var text = this.options['includeHtmlTags'][this.adaptor.kind(node)];
            this.extendString(node, text);
        }
        return this.adaptor.next(node);
    };
    HTMLDomStrings.prototype.handleContainer = function (node, ignore) {
        this.pushString();
        var cname = this.adaptor.getAttribute(node, 'class') || '';
        var tname = this.adaptor.kind(node) || '';
        var process = this.processHtmlClass.exec(cname);
        var next = node;
        if (this.adaptor.firstChild(node) && !this.adaptor.getAttribute(node, 'data-MJX') &&
            (process || !this.skipHtmlTags.exec(tname))) {
            if (this.adaptor.next(node)) {
                this.stack.push([this.adaptor.next(node), ignore]);
            }
            next = this.adaptor.firstChild(node);
            ignore = (ignore || this.ignoreHtmlClass.exec(cname)) && !process;
        }
        else {
            next = this.adaptor.next(node);
        }
        return [next, ignore];
    };
    HTMLDomStrings.prototype.handleOther = function (node, _ignore) {
        this.pushString();
        return this.adaptor.next(node);
    };
    HTMLDomStrings.prototype.find = function (node) {
        var _a, _b;
        this.init();
        var stop = this.adaptor.next(node);
        var ignore = false;
        var include = this.options['includeHtmlTags'];
        while (node && node !== stop) {
            var kind = this.adaptor.kind(node);
            if (kind === '#text') {
                node = this.handleText(node, ignore);
            }
            else if (include.hasOwnProperty(kind)) {
                node = this.handleTag(node, ignore);
            }
            else if (kind) {
                _a = __read(this.handleContainer(node, ignore), 2), node = _a[0], ignore = _a[1];
            }
            else {
                node = this.handleOther(node, ignore);
            }
            if (!node && this.stack.length) {
                this.pushString();
                _b = __read(this.stack.pop(), 2), node = _b[0], ignore = _b[1];
            }
        }
        this.pushString();
        var result = [this.strings, this.nodes];
        this.init();
        return result;
    };
    HTMLDomStrings.OPTIONS = {
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code', 'annotation', 'annotation-xml'],
        includeHtmlTags: { br: '\n', wbr: '', '#comment': '' },
        ignoreHtmlClass: 'mathjax_ignore',
        processHtmlClass: 'mathjax_process'
    };
    return HTMLDomStrings;
}());
exports.HTMLDomStrings = HTMLDomStrings;
//# sourceMappingURL=HTMLDomStrings.js.map

/***/ }),

/***/ 97471:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.HTMLHandler = void 0;
var Handler_js_1 = __webpack_require__(35663);
var HTMLDocument_js_1 = __webpack_require__(99331);
var HTMLHandler = (function (_super) {
    __extends(HTMLHandler, _super);
    function HTMLHandler() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.documentClass = HTMLDocument_js_1.HTMLDocument;
        return _this;
    }
    HTMLHandler.prototype.handlesDocument = function (document) {
        var adaptor = this.adaptor;
        if (typeof (document) === 'string') {
            try {
                document = adaptor.parse(document, 'text/html');
            }
            catch (err) { }
        }
        if (document instanceof adaptor.window.Document ||
            document instanceof adaptor.window.HTMLElement ||
            document instanceof adaptor.window.DocumentFragment) {
            return true;
        }
        return false;
    };
    HTMLHandler.prototype.create = function (document, options) {
        var adaptor = this.adaptor;
        if (typeof (document) === 'string') {
            document = adaptor.parse(document, 'text/html');
        }
        else if (document instanceof adaptor.window.HTMLElement ||
            document instanceof adaptor.window.DocumentFragment) {
            var child = document;
            document = adaptor.parse('', 'text/html');
            adaptor.append(adaptor.body(document), child);
        }
        return _super.prototype.create.call(this, document, options);
    };
    return HTMLHandler;
}(Handler_js_1.AbstractHandler));
exports.HTMLHandler = HTMLHandler;
//# sourceMappingURL=HTMLHandler.js.map

/***/ }),

/***/ 92870:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.HTMLMathItem = void 0;
var MathItem_js_1 = __webpack_require__(21605);
var HTMLMathItem = (function (_super) {
    __extends(HTMLMathItem, _super);
    function HTMLMathItem(math, jax, display, start, end) {
        if (display === void 0) { display = true; }
        if (start === void 0) { start = { node: null, n: 0, delim: '' }; }
        if (end === void 0) { end = { node: null, n: 0, delim: '' }; }
        return _super.call(this, math, jax, display, start, end) || this;
    }
    Object.defineProperty(HTMLMathItem.prototype, "adaptor", {
        get: function () {
            return this.inputJax.adaptor;
        },
        enumerable: false,
        configurable: true
    });
    HTMLMathItem.prototype.updateDocument = function (_html) {
        if (this.state() < MathItem_js_1.STATE.INSERTED) {
            if (this.inputJax.processStrings) {
                var node = this.start.node;
                if (node === this.end.node) {
                    if (this.end.n && this.end.n < this.adaptor.value(this.end.node).length) {
                        this.adaptor.split(this.end.node, this.end.n);
                    }
                    if (this.start.n) {
                        node = this.adaptor.split(this.start.node, this.start.n);
                    }
                    this.adaptor.replace(this.typesetRoot, node);
                }
                else {
                    if (this.start.n) {
                        node = this.adaptor.split(node, this.start.n);
                    }
                    while (node !== this.end.node) {
                        var next = this.adaptor.next(node);
                        this.adaptor.remove(node);
                        node = next;
                    }
                    this.adaptor.insert(this.typesetRoot, node);
                    if (this.end.n < this.adaptor.value(node).length) {
                        this.adaptor.split(node, this.end.n);
                    }
                    this.adaptor.remove(node);
                }
            }
            else {
                this.adaptor.replace(this.typesetRoot, this.start.node);
            }
            this.start.node = this.end.node = this.typesetRoot;
            this.start.n = this.end.n = 0;
            this.state(MathItem_js_1.STATE.INSERTED);
        }
    };
    HTMLMathItem.prototype.updateStyleSheet = function (document) {
        document.addStyleSheet();
    };
    HTMLMathItem.prototype.removeFromDocument = function (restore) {
        if (restore === void 0) { restore = false; }
        if (this.state() >= MathItem_js_1.STATE.TYPESET) {
            var adaptor = this.adaptor;
            var node = this.start.node;
            var math = adaptor.text('');
            if (restore) {
                var text = this.start.delim + this.math + this.end.delim;
                if (this.inputJax.processStrings) {
                    math = adaptor.text(text);
                }
                else {
                    var doc = adaptor.parse(text, 'text/html');
                    math = adaptor.firstChild(adaptor.body(doc));
                }
            }
            if (adaptor.parent(node)) {
                adaptor.replace(math, node);
            }
            this.start.node = this.end.node = math;
            this.start.n = this.end.n = 0;
        }
    };
    return HTMLMathItem;
}(MathItem_js_1.AbstractMathItem));
exports.HTMLMathItem = HTMLMathItem;
//# sourceMappingURL=HTMLMathItem.js.map

/***/ }),

/***/ 65337:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.HTMLMathList = void 0;
var MathList_js_1 = __webpack_require__(80426);
var HTMLMathList = (function (_super) {
    __extends(HTMLMathList, _super);
    function HTMLMathList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return HTMLMathList;
}(MathList_js_1.AbstractMathList));
exports.HTMLMathList = HTMLMathList;
//# sourceMappingURL=HTMLMathList.js.map

/***/ }),

/***/ 62844:
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
exports.BitFieldClass = exports.BitField = void 0;
var BitField = (function () {
    function BitField() {
        this.bits = 0;
    }
    BitField.allocate = function () {
        var e_1, _a;
        var names = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            names[_i] = arguments[_i];
        }
        try {
            for (var names_1 = __values(names), names_1_1 = names_1.next(); !names_1_1.done; names_1_1 = names_1.next()) {
                var name_1 = names_1_1.value;
                if (this.has(name_1)) {
                    throw new Error('Bit already allocated for ' + name_1);
                }
                if (this.next === BitField.MAXBIT) {
                    throw new Error('Maximum number of bits already allocated');
                }
                this.names.set(name_1, this.next);
                this.next <<= 1;
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (names_1_1 && !names_1_1.done && (_a = names_1.return)) _a.call(names_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
    };
    BitField.has = function (name) {
        return this.names.has(name);
    };
    BitField.prototype.set = function (name) {
        this.bits |= this.getBit(name);
    };
    BitField.prototype.clear = function (name) {
        this.bits &= ~this.getBit(name);
    };
    BitField.prototype.isSet = function (name) {
        return !!(this.bits & this.getBit(name));
    };
    BitField.prototype.reset = function () {
        this.bits = 0;
    };
    BitField.prototype.getBit = function (name) {
        var bit = this.constructor.names.get(name);
        if (!bit) {
            throw new Error('Unknown bit-field name: ' + name);
        }
        return bit;
    };
    BitField.MAXBIT = 1 << 31;
    BitField.next = 1;
    BitField.names = new Map();
    return BitField;
}());
exports.BitField = BitField;
function BitFieldClass() {
    var names = [];
    for (var _i = 0; _i < arguments.length; _i++) {
        names[_i] = arguments[_i];
    }
    var Bits = (function (_super) {
        __extends(Bits, _super);
        function Bits() {
            return _super !== null && _super.apply(this, arguments) || this;
        }
        return Bits;
    }(BitField));
    Bits.allocate.apply(Bits, __spreadArray([], __read(names), false));
    return Bits;
}
exports.BitFieldClass = BitFieldClass;
//# sourceMappingURL=BitField.js.map

/***/ }),

/***/ 18341:
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
exports.FunctionList = void 0;
var PrioritizedList_js_1 = __webpack_require__(98721);
var FunctionList = (function (_super) {
    __extends(FunctionList, _super);
    function FunctionList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    FunctionList.prototype.execute = function () {
        var e_1, _a;
        var data = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            data[_i] = arguments[_i];
        }
        try {
            for (var _b = __values(this), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                var result = item.item.apply(item, __spreadArray([], __read(data), false));
                if (result === false) {
                    return false;
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
        return true;
    };
    FunctionList.prototype.asyncExecute = function () {
        var data = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            data[_i] = arguments[_i];
        }
        var i = -1;
        var items = this.items;
        return new Promise(function (ok, fail) {
            (function execute() {
                var _a;
                while (++i < items.length) {
                    var result = (_a = items[i]).item.apply(_a, __spreadArray([], __read(data), false));
                    if (result instanceof Promise) {
                        result.then(execute).catch(function (err) { return fail(err); });
                        return;
                    }
                    if (result === false) {
                        ok(false);
                        return;
                    }
                }
                ok(true);
            })();
        });
    };
    return FunctionList;
}(PrioritizedList_js_1.PrioritizedList));
exports.FunctionList = FunctionList;
//# sourceMappingURL=FunctionList.js.map

/***/ }),

/***/ 6811:
/***/ (function(__unused_webpack_module, exports) {


var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = op[0] & 2 ? y["return"] : op[0] ? y["throw"] || ((t = y["return"]) && t.call(y), 0) : y.next) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [op[0] & 2, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
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
exports.LinkedList = exports.ListItem = exports.END = void 0;
exports.END = Symbol();
var ListItem = (function () {
    function ListItem(data) {
        if (data === void 0) { data = null; }
        this.next = null;
        this.prev = null;
        this.data = data;
    }
    return ListItem;
}());
exports.ListItem = ListItem;
var LinkedList = (function () {
    function LinkedList() {
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        this.list = new ListItem(exports.END);
        this.list.next = this.list.prev = this.list;
        this.push.apply(this, __spreadArray([], __read(args), false));
    }
    LinkedList.prototype.isBefore = function (a, b) {
        return a < b;
    };
    LinkedList.prototype.push = function () {
        var e_1, _a;
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        try {
            for (var args_1 = __values(args), args_1_1 = args_1.next(); !args_1_1.done; args_1_1 = args_1.next()) {
                var data = args_1_1.value;
                var item = new ListItem(data);
                item.next = this.list;
                item.prev = this.list.prev;
                this.list.prev = item;
                item.prev.next = item;
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (args_1_1 && !args_1_1.done && (_a = args_1.return)) _a.call(args_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return this;
    };
    LinkedList.prototype.pop = function () {
        var item = this.list.prev;
        if (item.data === exports.END) {
            return null;
        }
        this.list.prev = item.prev;
        item.prev.next = this.list;
        item.next = item.prev = null;
        return item.data;
    };
    LinkedList.prototype.unshift = function () {
        var e_2, _a;
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        try {
            for (var _b = __values(args.slice(0).reverse()), _c = _b.next(); !_c.done; _c = _b.next()) {
                var data = _c.value;
                var item = new ListItem(data);
                item.next = this.list.next;
                item.prev = this.list;
                this.list.next = item;
                item.next.prev = item;
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_2) throw e_2.error; }
        }
        return this;
    };
    LinkedList.prototype.shift = function () {
        var item = this.list.next;
        if (item.data === exports.END) {
            return null;
        }
        this.list.next = item.next;
        item.next.prev = this.list;
        item.next = item.prev = null;
        return item.data;
    };
    LinkedList.prototype.remove = function () {
        var e_3, _a;
        var items = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            items[_i] = arguments[_i];
        }
        var map = new Map();
        try {
            for (var items_1 = __values(items), items_1_1 = items_1.next(); !items_1_1.done; items_1_1 = items_1.next()) {
                var item_1 = items_1_1.value;
                map.set(item_1, true);
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (items_1_1 && !items_1_1.done && (_a = items_1.return)) _a.call(items_1);
            }
            finally { if (e_3) throw e_3.error; }
        }
        var item = this.list.next;
        while (item.data !== exports.END) {
            var next = item.next;
            if (map.has(item.data)) {
                item.prev.next = item.next;
                item.next.prev = item.prev;
                item.next = item.prev = null;
            }
            item = next;
        }
    };
    LinkedList.prototype.clear = function () {
        this.list.next.prev = this.list.prev.next = null;
        this.list.next = this.list.prev = this.list;
        return this;
    };
    LinkedList.prototype[Symbol.iterator] = function () {
        var current;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    current = this.list.next;
                    _a.label = 1;
                case 1:
                    if (!(current.data !== exports.END)) return [3, 3];
                    return [4, current.data];
                case 2:
                    _a.sent();
                    current = current.next;
                    return [3, 1];
                case 3: return [2];
            }
        });
    };
    LinkedList.prototype.reversed = function () {
        var current;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    current = this.list.prev;
                    _a.label = 1;
                case 1:
                    if (!(current.data !== exports.END)) return [3, 3];
                    return [4, current.data];
                case 2:
                    _a.sent();
                    current = current.prev;
                    return [3, 1];
                case 3: return [2];
            }
        });
    };
    LinkedList.prototype.insert = function (data, isBefore) {
        if (isBefore === void 0) { isBefore = null; }
        if (isBefore === null) {
            isBefore = this.isBefore.bind(this);
        }
        var item = new ListItem(data);
        var cur = this.list.next;
        while (cur.data !== exports.END && isBefore(cur.data, item.data)) {
            cur = cur.next;
        }
        item.prev = cur.prev;
        item.next = cur;
        cur.prev.next = cur.prev = item;
        return this;
    };
    LinkedList.prototype.sort = function (isBefore) {
        var e_4, _a;
        if (isBefore === void 0) { isBefore = null; }
        if (isBefore === null) {
            isBefore = this.isBefore.bind(this);
        }
        var lists = [];
        try {
            for (var _b = __values(this), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                lists.push(new LinkedList(item));
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_4) throw e_4.error; }
        }
        this.list.next = this.list.prev = this.list;
        while (lists.length > 1) {
            var l1 = lists.shift();
            var l2 = lists.shift();
            l1.merge(l2, isBefore);
            lists.push(l1);
        }
        if (lists.length) {
            this.list = lists[0].list;
        }
        return this;
    };
    LinkedList.prototype.merge = function (list, isBefore) {
        var _a, _b, _c, _d, _e;
        if (isBefore === void 0) { isBefore = null; }
        if (isBefore === null) {
            isBefore = this.isBefore.bind(this);
        }
        var lcur = this.list.next;
        var mcur = list.list.next;
        while (lcur.data !== exports.END && mcur.data !== exports.END) {
            if (isBefore(mcur.data, lcur.data)) {
                _a = __read([lcur, mcur], 2), mcur.prev.next = _a[0], lcur.prev.next = _a[1];
                _b = __read([lcur.prev, mcur.prev], 2), mcur.prev = _b[0], lcur.prev = _b[1];
                _c = __read([list.list, this.list], 2), this.list.prev.next = _c[0], list.list.prev.next = _c[1];
                _d = __read([list.list.prev, this.list.prev], 2), this.list.prev = _d[0], list.list.prev = _d[1];
                _e = __read([mcur.next, lcur], 2), lcur = _e[0], mcur = _e[1];
            }
            else {
                lcur = lcur.next;
            }
        }
        if (mcur.data !== exports.END) {
            this.list.prev.next = list.list.next;
            list.list.next.prev = this.list.prev;
            list.list.prev.next = this.list;
            this.list.prev = list.list.prev;
            list.list.next = list.list.prev = list.list;
        }
        return this;
    };
    return LinkedList;
}());
exports.LinkedList = LinkedList;
//# sourceMappingURL=LinkedList.js.map

/***/ }),

/***/ 98721:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.PrioritizedList = void 0;
var PrioritizedList = (function () {
    function PrioritizedList() {
        this.items = [];
        this.items = [];
    }
    PrioritizedList.prototype[Symbol.iterator] = function () {
        var i = 0;
        var items = this.items;
        return {
            next: function () {
                return { value: items[i++], done: (i > items.length) };
            }
        };
    };
    PrioritizedList.prototype.add = function (item, priority) {
        if (priority === void 0) { priority = PrioritizedList.DEFAULTPRIORITY; }
        var i = this.items.length;
        do {
            i--;
        } while (i >= 0 && priority < this.items[i].priority);
        this.items.splice(i + 1, 0, { item: item, priority: priority });
        return item;
    };
    PrioritizedList.prototype.remove = function (item) {
        var i = this.items.length;
        do {
            i--;
        } while (i >= 0 && this.items[i].item !== item);
        if (i >= 0) {
            this.items.splice(i, 1);
        }
    };
    PrioritizedList.DEFAULTPRIORITY = 5;
    return PrioritizedList;
}());
exports.PrioritizedList = PrioritizedList;
//# sourceMappingURL=PrioritizedList.js.map

/***/ })

}]);
//# sourceMappingURL=7471.27c6037e2917dcd9958a.js.map?v=27c6037e2917dcd9958a