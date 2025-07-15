"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4971],{

/***/ 27573:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.VERSION = void 0;
exports.VERSION = '3.2.2';
//# sourceMappingURL=version.js.map

/***/ }),

/***/ 12514:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.HandlerList = void 0;
var PrioritizedList_js_1 = __webpack_require__(98721);
var HandlerList = (function (_super) {
    __extends(HandlerList, _super);
    function HandlerList() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    HandlerList.prototype.register = function (handler) {
        return this.add(handler, handler.priority);
    };
    HandlerList.prototype.unregister = function (handler) {
        this.remove(handler);
    };
    HandlerList.prototype.handlesDocument = function (document) {
        var e_1, _a;
        try {
            for (var _b = __values(this), _c = _b.next(); !_c.done; _c = _b.next()) {
                var item = _c.value;
                var handler = item.item;
                if (handler.handlesDocument(document)) {
                    return handler;
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
        throw new Error("Can't find handler for document");
    };
    HandlerList.prototype.document = function (document, options) {
        if (options === void 0) { options = null; }
        return this.handlesDocument(document).create(document, options);
    };
    return HandlerList;
}(PrioritizedList_js_1.PrioritizedList));
exports.HandlerList = HandlerList;
//# sourceMappingURL=HandlerList.js.map

/***/ }),

/***/ 44971:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.mathjax = void 0;
var version_js_1 = __webpack_require__(27573);
var HandlerList_js_1 = __webpack_require__(12514);
var Retries_js_1 = __webpack_require__(10956);
exports.mathjax = {
    version: version_js_1.VERSION,
    handlers: new HandlerList_js_1.HandlerList(),
    document: function (document, options) {
        return exports.mathjax.handlers.document(document, options);
    },
    handleRetriesFor: Retries_js_1.handleRetriesFor,
    retryAfter: Retries_js_1.retryAfter,
    asyncLoad: null,
};
//# sourceMappingURL=mathjax.js.map

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

/***/ }),

/***/ 10956:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.retryAfter = exports.handleRetriesFor = void 0;
function handleRetriesFor(code) {
    return new Promise(function run(ok, fail) {
        try {
            ok(code());
        }
        catch (err) {
            if (err.retry && err.retry instanceof Promise) {
                err.retry.then(function () { return run(ok, fail); })
                    .catch(function (perr) { return fail(perr); });
            }
            else if (err.restart && err.restart.isCallback) {
                MathJax.Callback.After(function () { return run(ok, fail); }, err.restart);
            }
            else {
                fail(err);
            }
        }
    });
}
exports.handleRetriesFor = handleRetriesFor;
function retryAfter(promise) {
    var err = new Error('MathJax retry');
    err.retry = promise;
    throw err;
}
exports.retryAfter = retryAfter;
//# sourceMappingURL=Retries.js.map

/***/ })

}]);
//# sourceMappingURL=4971.e850b0a1dcb6d3fce7a4.js.map?v=e850b0a1dcb6d3fce7a4