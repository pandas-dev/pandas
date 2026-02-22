"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[792],{

/***/ 53821:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MathJax = exports.combineWithMathJax = exports.combineDefaults = exports.combineConfig = exports.isObject = void 0;
var version_js_1 = __webpack_require__(27573);
function isObject(x) {
    return typeof x === 'object' && x !== null;
}
exports.isObject = isObject;
function combineConfig(dst, src) {
    var e_1, _a;
    try {
        for (var _b = __values(Object.keys(src)), _c = _b.next(); !_c.done; _c = _b.next()) {
            var id = _c.value;
            if (id === '__esModule')
                continue;
            if (isObject(dst[id]) && isObject(src[id]) &&
                !(src[id] instanceof Promise)) {
                combineConfig(dst[id], src[id]);
            }
            else if (src[id] !== null && src[id] !== undefined) {
                dst[id] = src[id];
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
    return dst;
}
exports.combineConfig = combineConfig;
function combineDefaults(dst, name, src) {
    var e_2, _a;
    if (!dst[name]) {
        dst[name] = {};
    }
    dst = dst[name];
    try {
        for (var _b = __values(Object.keys(src)), _c = _b.next(); !_c.done; _c = _b.next()) {
            var id = _c.value;
            if (isObject(dst[id]) && isObject(src[id])) {
                combineDefaults(dst, id, src[id]);
            }
            else if (dst[id] == null && src[id] != null) {
                dst[id] = src[id];
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
    return dst;
}
exports.combineDefaults = combineDefaults;
function combineWithMathJax(config) {
    return combineConfig(exports.MathJax, config);
}
exports.combineWithMathJax = combineWithMathJax;
if (typeof __webpack_require__.g.MathJax === 'undefined') {
    __webpack_require__.g.MathJax = {};
}
if (!__webpack_require__.g.MathJax.version) {
    __webpack_require__.g.MathJax = {
        version: version_js_1.VERSION,
        _: {},
        config: __webpack_require__.g.MathJax
    };
}
exports.MathJax = __webpack_require__.g.MathJax;
//# sourceMappingURL=global.js.map

/***/ }),

/***/ 14187:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {

var __dirname = "/";

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
var e_1, _a;
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.CONFIG = exports.MathJax = exports.Loader = exports.PathFilters = exports.PackageError = exports.Package = void 0;
var global_js_1 = __webpack_require__(53821);
var package_js_1 = __webpack_require__(91593);
var package_js_2 = __webpack_require__(91593);
Object.defineProperty(exports, "Package", ({ enumerable: true, get: function () { return package_js_2.Package; } }));
Object.defineProperty(exports, "PackageError", ({ enumerable: true, get: function () { return package_js_2.PackageError; } }));
var FunctionList_js_1 = __webpack_require__(18341);
exports.PathFilters = {
    source: function (data) {
        if (exports.CONFIG.source.hasOwnProperty(data.name)) {
            data.name = exports.CONFIG.source[data.name];
        }
        return true;
    },
    normalize: function (data) {
        var name = data.name;
        if (!name.match(/^(?:[a-z]+:\/)?\/|[a-z]:\\|\[/i)) {
            data.name = '[mathjax]/' + name.replace(/^\.\//, '');
        }
        if (data.addExtension && !name.match(/\.[^\/]+$/)) {
            data.name += '.js';
        }
        return true;
    },
    prefix: function (data) {
        var match;
        while ((match = data.name.match(/^\[([^\]]*)\]/))) {
            if (!exports.CONFIG.paths.hasOwnProperty(match[1]))
                break;
            data.name = exports.CONFIG.paths[match[1]] + data.name.substr(match[0].length);
        }
        return true;
    }
};
var Loader;
(function (Loader) {
    var VERSION = global_js_1.MathJax.version;
    Loader.versions = new Map();
    function ready() {
        var e_2, _a;
        var names = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            names[_i] = arguments[_i];
        }
        if (names.length === 0) {
            names = Array.from(package_js_1.Package.packages.keys());
        }
        var promises = [];
        try {
            for (var names_1 = __values(names), names_1_1 = names_1.next(); !names_1_1.done; names_1_1 = names_1.next()) {
                var name_1 = names_1_1.value;
                var extension = package_js_1.Package.packages.get(name_1) || new package_js_1.Package(name_1, true);
                promises.push(extension.promise);
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (names_1_1 && !names_1_1.done && (_a = names_1.return)) _a.call(names_1);
            }
            finally { if (e_2) throw e_2.error; }
        }
        return Promise.all(promises);
    }
    Loader.ready = ready;
    function load() {
        var e_3, _a;
        var names = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            names[_i] = arguments[_i];
        }
        if (names.length === 0) {
            return Promise.resolve();
        }
        var promises = [];
        var _loop_1 = function (name_2) {
            var extension = package_js_1.Package.packages.get(name_2);
            if (!extension) {
                extension = new package_js_1.Package(name_2);
                extension.provides(exports.CONFIG.provides[name_2]);
            }
            extension.checkNoLoad();
            promises.push(extension.promise.then(function () {
                if (!exports.CONFIG.versionWarnings)
                    return;
                if (extension.isLoaded && !Loader.versions.has(package_js_1.Package.resolvePath(name_2))) {
                    console.warn("No version information available for component ".concat(name_2));
                }
            }));
        };
        try {
            for (var names_2 = __values(names), names_2_1 = names_2.next(); !names_2_1.done; names_2_1 = names_2.next()) {
                var name_2 = names_2_1.value;
                _loop_1(name_2);
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (names_2_1 && !names_2_1.done && (_a = names_2.return)) _a.call(names_2);
            }
            finally { if (e_3) throw e_3.error; }
        }
        package_js_1.Package.loadAll();
        return Promise.all(promises);
    }
    Loader.load = load;
    function preLoad() {
        var e_4, _a;
        var names = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            names[_i] = arguments[_i];
        }
        try {
            for (var names_3 = __values(names), names_3_1 = names_3.next(); !names_3_1.done; names_3_1 = names_3.next()) {
                var name_3 = names_3_1.value;
                var extension = package_js_1.Package.packages.get(name_3);
                if (!extension) {
                    extension = new package_js_1.Package(name_3, true);
                    extension.provides(exports.CONFIG.provides[name_3]);
                }
                extension.loaded();
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (names_3_1 && !names_3_1.done && (_a = names_3.return)) _a.call(names_3);
            }
            finally { if (e_4) throw e_4.error; }
        }
    }
    Loader.preLoad = preLoad;
    function defaultReady() {
        if (typeof exports.MathJax.startup !== 'undefined') {
            exports.MathJax.config.startup.ready();
        }
    }
    Loader.defaultReady = defaultReady;
    function getRoot() {
        var root = __dirname + '/../../es5';
        if (typeof document !== 'undefined') {
            var script = document.currentScript || document.getElementById('MathJax-script');
            if (script) {
                root = script.src.replace(/\/[^\/]*$/, '');
            }
        }
        return root;
    }
    Loader.getRoot = getRoot;
    function checkVersion(name, version, _type) {
        Loader.versions.set(package_js_1.Package.resolvePath(name), VERSION);
        if (exports.CONFIG.versionWarnings && version !== VERSION) {
            console.warn("Component ".concat(name, " uses ").concat(version, " of MathJax; version in use is ").concat(VERSION));
            return true;
        }
        return false;
    }
    Loader.checkVersion = checkVersion;
    Loader.pathFilters = new FunctionList_js_1.FunctionList();
    Loader.pathFilters.add(exports.PathFilters.source, 0);
    Loader.pathFilters.add(exports.PathFilters.normalize, 10);
    Loader.pathFilters.add(exports.PathFilters.prefix, 20);
})(Loader = exports.Loader || (exports.Loader = {}));
exports.MathJax = global_js_1.MathJax;
if (typeof exports.MathJax.loader === 'undefined') {
    (0, global_js_1.combineDefaults)(exports.MathJax.config, 'loader', {
        paths: {
            mathjax: Loader.getRoot()
        },
        source: {},
        dependencies: {},
        provides: {},
        load: [],
        ready: Loader.defaultReady.bind(Loader),
        failed: function (error) { return console.log("MathJax(".concat(error.package || '?', "): ").concat(error.message)); },
        require: null,
        pathFilters: [],
        versionWarnings: true
    });
    (0, global_js_1.combineWithMathJax)({
        loader: Loader
    });
    try {
        for (var _b = __values(exports.MathJax.config.loader.pathFilters), _c = _b.next(); !_c.done; _c = _b.next()) {
            var filter = _c.value;
            if (Array.isArray(filter)) {
                Loader.pathFilters.add(filter[0], filter[1]);
            }
            else {
                Loader.pathFilters.add(filter);
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
}
exports.CONFIG = exports.MathJax.config.loader;
//# sourceMappingURL=loader.js.map

/***/ }),

/***/ 91593:
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
exports.Package = exports.PackageError = void 0;
var loader_js_1 = __webpack_require__(14187);
var PackageError = (function (_super) {
    __extends(PackageError, _super);
    function PackageError(message, name) {
        var _this = _super.call(this, message) || this;
        _this.package = name;
        return _this;
    }
    return PackageError;
}(Error));
exports.PackageError = PackageError;
var Package = (function () {
    function Package(name, noLoad) {
        if (noLoad === void 0) { noLoad = false; }
        this.isLoaded = false;
        this.isLoading = false;
        this.hasFailed = false;
        this.dependents = [];
        this.dependencies = [];
        this.dependencyCount = 0;
        this.provided = [];
        this.name = name;
        this.noLoad = noLoad;
        Package.packages.set(name, this);
        this.promise = this.makePromise(this.makeDependencies());
    }
    Object.defineProperty(Package.prototype, "canLoad", {
        get: function () {
            return this.dependencyCount === 0 && !this.noLoad && !this.isLoading && !this.hasFailed;
        },
        enumerable: false,
        configurable: true
    });
    Package.resolvePath = function (name, addExtension) {
        if (addExtension === void 0) { addExtension = true; }
        var data = { name: name, original: name, addExtension: addExtension };
        loader_js_1.Loader.pathFilters.execute(data);
        return data.name;
    };
    Package.loadAll = function () {
        var e_1, _a;
        try {
            for (var _b = __values(this.packages.values()), _c = _b.next(); !_c.done; _c = _b.next()) {
                var extension = _c.value;
                if (extension.canLoad) {
                    extension.load();
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
    Package.prototype.makeDependencies = function () {
        var e_2, _a;
        var promises = [];
        var map = Package.packages;
        var noLoad = this.noLoad;
        var name = this.name;
        var dependencies = [];
        if (loader_js_1.CONFIG.dependencies.hasOwnProperty(name)) {
            dependencies.push.apply(dependencies, __spreadArray([], __read(loader_js_1.CONFIG.dependencies[name]), false));
        }
        else if (name !== 'core') {
            dependencies.push('core');
        }
        try {
            for (var dependencies_1 = __values(dependencies), dependencies_1_1 = dependencies_1.next(); !dependencies_1_1.done; dependencies_1_1 = dependencies_1.next()) {
                var dependent = dependencies_1_1.value;
                var extension = map.get(dependent) || new Package(dependent, noLoad);
                if (this.dependencies.indexOf(extension) < 0) {
                    extension.addDependent(this, noLoad);
                    this.dependencies.push(extension);
                    if (!extension.isLoaded) {
                        this.dependencyCount++;
                        promises.push(extension.promise);
                    }
                }
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (dependencies_1_1 && !dependencies_1_1.done && (_a = dependencies_1.return)) _a.call(dependencies_1);
            }
            finally { if (e_2) throw e_2.error; }
        }
        return promises;
    };
    Package.prototype.makePromise = function (promises) {
        var _this = this;
        var promise = new Promise((function (resolve, reject) {
            _this.resolve = resolve;
            _this.reject = reject;
        }));
        var config = (loader_js_1.CONFIG[this.name] || {});
        if (config.ready) {
            promise = promise.then(function (_name) { return config.ready(_this.name); });
        }
        if (promises.length) {
            promises.push(promise);
            promise = Promise.all(promises).then(function (names) { return names.join(', '); });
        }
        if (config.failed) {
            promise.catch(function (message) { return config.failed(new PackageError(message, _this.name)); });
        }
        return promise;
    };
    Package.prototype.load = function () {
        if (!this.isLoaded && !this.isLoading && !this.noLoad) {
            this.isLoading = true;
            var url = Package.resolvePath(this.name);
            if (loader_js_1.CONFIG.require) {
                this.loadCustom(url);
            }
            else {
                this.loadScript(url);
            }
        }
    };
    Package.prototype.loadCustom = function (url) {
        var _this = this;
        try {
            var result = loader_js_1.CONFIG.require(url);
            if (result instanceof Promise) {
                result.then(function () { return _this.checkLoad(); })
                    .catch(function (err) { return _this.failed('Can\'t load "' + url + '"\n' + err.message.trim()); });
            }
            else {
                this.checkLoad();
            }
        }
        catch (err) {
            this.failed(err.message);
        }
    };
    Package.prototype.loadScript = function (url) {
        var _this = this;
        var script = document.createElement('script');
        script.src = url;
        script.charset = 'UTF-8';
        script.onload = function (_event) { return _this.checkLoad(); };
        script.onerror = function (_event) { return _this.failed('Can\'t load "' + url + '"'); };
        document.head.appendChild(script);
    };
    Package.prototype.loaded = function () {
        var e_3, _a, e_4, _b;
        this.isLoaded = true;
        this.isLoading = false;
        try {
            for (var _c = __values(this.dependents), _d = _c.next(); !_d.done; _d = _c.next()) {
                var dependent = _d.value;
                dependent.requirementSatisfied();
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_3) throw e_3.error; }
        }
        try {
            for (var _e = __values(this.provided), _f = _e.next(); !_f.done; _f = _e.next()) {
                var provided = _f.value;
                provided.loaded();
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
            }
            finally { if (e_4) throw e_4.error; }
        }
        this.resolve(this.name);
    };
    Package.prototype.failed = function (message) {
        this.hasFailed = true;
        this.isLoading = false;
        this.reject(new PackageError(message, this.name));
    };
    Package.prototype.checkLoad = function () {
        var _this = this;
        var config = (loader_js_1.CONFIG[this.name] || {});
        var checkReady = config.checkReady || (function () { return Promise.resolve(); });
        checkReady().then(function () { return _this.loaded(); })
            .catch(function (message) { return _this.failed(message); });
    };
    Package.prototype.requirementSatisfied = function () {
        if (this.dependencyCount) {
            this.dependencyCount--;
            if (this.canLoad) {
                this.load();
            }
        }
    };
    Package.prototype.provides = function (names) {
        var e_5, _a;
        if (names === void 0) { names = []; }
        try {
            for (var names_1 = __values(names), names_1_1 = names_1.next(); !names_1_1.done; names_1_1 = names_1.next()) {
                var name_1 = names_1_1.value;
                var provided = Package.packages.get(name_1);
                if (!provided) {
                    if (!loader_js_1.CONFIG.dependencies[name_1]) {
                        loader_js_1.CONFIG.dependencies[name_1] = [];
                    }
                    loader_js_1.CONFIG.dependencies[name_1].push(name_1);
                    provided = new Package(name_1, true);
                    provided.isLoading = true;
                }
                this.provided.push(provided);
            }
        }
        catch (e_5_1) { e_5 = { error: e_5_1 }; }
        finally {
            try {
                if (names_1_1 && !names_1_1.done && (_a = names_1.return)) _a.call(names_1);
            }
            finally { if (e_5) throw e_5.error; }
        }
    };
    Package.prototype.addDependent = function (extension, noLoad) {
        this.dependents.push(extension);
        if (!noLoad) {
            this.checkNoLoad();
        }
    };
    Package.prototype.checkNoLoad = function () {
        var e_6, _a;
        if (this.noLoad) {
            this.noLoad = false;
            try {
                for (var _b = __values(this.dependencies), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var dependency = _c.value;
                    dependency.checkNoLoad();
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
    };
    Package.packages = new Map();
    return Package;
}());
exports.Package = Package;
//# sourceMappingURL=package.js.map

/***/ }),

/***/ 20792:
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
var __spreadArray = (this && this.__spreadArray) || function (to, from, pack) {
    if (pack || arguments.length === 2) for (var i = 0, l = from.length, ar; i < l; i++) {
        if (ar || !(i in from)) {
            if (!ar) ar = Array.prototype.slice.call(from, 0, i);
            ar[i] = from[i];
        }
    }
    return to.concat(ar || Array.prototype.slice.call(from));
};
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.RequireConfiguration = exports.options = exports.RequireMethods = exports.RequireLoad = void 0;
var Configuration_js_1 = __webpack_require__(63401);
var SymbolMap_js_1 = __webpack_require__(65695);
var TexError_js_1 = __importDefault(__webpack_require__(54420));
var global_js_1 = __webpack_require__(53821);
var package_js_1 = __webpack_require__(91593);
var loader_js_1 = __webpack_require__(14187);
var mathjax_js_1 = __webpack_require__(44971);
var Options_js_1 = __webpack_require__(4498);
var MJCONFIG = global_js_1.MathJax.config;
function RegisterExtension(jax, name) {
    var _a;
    var require = jax.parseOptions.options.require;
    var required = jax.parseOptions.packageData.get('require').required;
    var extension = name.substr(require.prefix.length);
    if (required.indexOf(extension) < 0) {
        required.push(extension);
        RegisterDependencies(jax, loader_js_1.CONFIG.dependencies[name]);
        var handler = Configuration_js_1.ConfigurationHandler.get(extension);
        if (handler) {
            var options_1 = MJCONFIG[name] || {};
            if (handler.options && Object.keys(handler.options).length === 1 && handler.options[extension]) {
                options_1 = (_a = {}, _a[extension] = options_1, _a);
            }
            jax.configuration.add(extension, jax, options_1);
            var configured = jax.parseOptions.packageData.get('require').configured;
            if (handler.preprocessors.length && !configured.has(extension)) {
                configured.set(extension, true);
                mathjax_js_1.mathjax.retryAfter(Promise.resolve());
            }
        }
    }
}
function RegisterDependencies(jax, names) {
    var e_1, _a;
    if (names === void 0) { names = []; }
    var prefix = jax.parseOptions.options.require.prefix;
    try {
        for (var names_1 = __values(names), names_1_1 = names_1.next(); !names_1_1.done; names_1_1 = names_1.next()) {
            var name_1 = names_1_1.value;
            if (name_1.substr(0, prefix.length) === prefix) {
                RegisterExtension(jax, name_1);
            }
        }
    }
    catch (e_1_1) { e_1 = { error: e_1_1 }; }
    finally {
        try {
            if (names_1_1 && !names_1_1.done && (_a = names_1.return)) _a.call(names_1);
        }
        finally { if (e_1) throw e_1.error; }
    }
}
function RequireLoad(parser, name) {
    var options = parser.options.require;
    var allow = options.allow;
    var extension = (name.substr(0, 1) === '[' ? '' : options.prefix) + name;
    var allowed = (allow.hasOwnProperty(extension) ? allow[extension] :
        allow.hasOwnProperty(name) ? allow[name] : options.defaultAllow);
    if (!allowed) {
        throw new TexError_js_1.default('BadRequire', 'Extension "%1" is not allowed to be loaded', extension);
    }
    if (package_js_1.Package.packages.has(extension)) {
        RegisterExtension(parser.configuration.packageData.get('require').jax, extension);
    }
    else {
        mathjax_js_1.mathjax.retryAfter(loader_js_1.Loader.load(extension));
    }
}
exports.RequireLoad = RequireLoad;
function config(_config, jax) {
    jax.parseOptions.packageData.set('require', {
        jax: jax,
        required: __spreadArray([], __read(jax.options.packages), false),
        configured: new Map()
    });
    var options = jax.parseOptions.options.require;
    var prefix = options.prefix;
    if (prefix.match(/[^_a-zA-Z0-9]/)) {
        throw Error('Illegal characters used in \\require prefix');
    }
    if (!loader_js_1.CONFIG.paths[prefix]) {
        loader_js_1.CONFIG.paths[prefix] = '[mathjax]/input/tex/extensions';
    }
    options.prefix = '[' + prefix + ']/';
}
exports.RequireMethods = {
    Require: function (parser, name) {
        var required = parser.GetArgument(name);
        if (required.match(/[^_a-zA-Z0-9]/) || required === '') {
            throw new TexError_js_1.default('BadPackageName', 'Argument for %1 is not a valid package name', name);
        }
        RequireLoad(parser, required);
    }
};
exports.options = {
    require: {
        allow: (0, Options_js_1.expandable)({
            base: false,
            'all-packages': false,
            autoload: false,
            configmacros: false,
            tagformat: false,
            setoptions: false
        }),
        defaultAllow: true,
        prefix: 'tex'
    }
};
new SymbolMap_js_1.CommandMap('require', { require: 'Require' }, exports.RequireMethods);
exports.RequireConfiguration = Configuration_js_1.Configuration.create('require', { handler: { macro: ['require'] }, config: config, options: exports.options });
//# sourceMappingURL=RequireConfiguration.js.map

/***/ })

}]);
//# sourceMappingURL=792.050c0efb8da8e633f900.js.map?v=050c0efb8da8e633f900