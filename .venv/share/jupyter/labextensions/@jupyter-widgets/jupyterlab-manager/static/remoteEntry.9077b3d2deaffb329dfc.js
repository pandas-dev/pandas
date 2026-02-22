var _JUPYTERLAB;
/******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "webpack/container/entry/@jupyter-widgets/jupyterlab-manager":
/*!***********************!*\
  !*** container entry ***!
  \***********************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

var moduleMap = {
	"./index": () => {
		return Promise.all([__webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling"), __webpack_require__.e("lib_index_js")]).then(() => (() => ((__webpack_require__(/*! ./lib/index.js */ "./lib/index.js")))));
	},
	"./extension": () => {
		return Promise.all([__webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling"), __webpack_require__.e("lib_index_js")]).then(() => (() => ((__webpack_require__(/*! ./lib/index.js */ "./lib/index.js")))));
	}
};
var get = (module, getScope) => {
	__webpack_require__.R = getScope;
	getScope = (
		__webpack_require__.o(moduleMap, module)
			? moduleMap[module]()
			: Promise.resolve().then(() => {
				throw new Error('Module "' + module + '" does not exist in container.');
			})
	);
	__webpack_require__.R = undefined;
	return getScope;
};
var init = (shareScope, initScope) => {
	if (!__webpack_require__.S) return;
	var name = "default"
	var oldScope = __webpack_require__.S[name];
	if(oldScope && oldScope !== shareScope) throw new Error("Container initialization failed as it has already been initialized with a different share scope");
	__webpack_require__.S[name] = shareScope;
	return __webpack_require__.I(name, initScope);
};

// This exports getters to disallow modifications
__webpack_require__.d(exports, {
	get: () => (get),
	init: () => (init)
});

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			id: moduleId,
/******/ 			loaded: false,
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = __webpack_modules__;
/******/ 	
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = __webpack_module_cache__;
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/compat get default export */
/******/ 	(() => {
/******/ 		// getDefaultExport function for compatibility with non-harmony modules
/******/ 		__webpack_require__.n = (module) => {
/******/ 			var getter = module && module.__esModule ?
/******/ 				() => (module['default']) :
/******/ 				() => (module);
/******/ 			__webpack_require__.d(getter, { a: getter });
/******/ 			return getter;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/define property getters */
/******/ 	(() => {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = (exports, definition) => {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/ensure chunk */
/******/ 	(() => {
/******/ 		__webpack_require__.f = {};
/******/ 		// This file contains only the entry chunk.
/******/ 		// The chunk loading function for additional chunks
/******/ 		__webpack_require__.e = (chunkId) => {
/******/ 			return Promise.all(Object.keys(__webpack_require__.f).reduce((promises, key) => {
/******/ 				__webpack_require__.f[key](chunkId, promises);
/******/ 				return promises;
/******/ 			}, []));
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/get javascript chunk filename */
/******/ 	(() => {
/******/ 		// This function allow to reference async chunks
/******/ 		__webpack_require__.u = (chunkId) => {
/******/ 			// return url for filenames based on template
/******/ 			return "" + chunkId + "." + {"webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery":"65fdeedf445afbd3c750","webpack_sharing_consume_default_lumino_coreutils":"3cdb0a7dfa5b1685ad03","webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base":"14f8904ad66e5de4794a","webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling":"121fc1d57584a4cad646","lib_index_js":"56f37fb50492d0d63a45","vendors-node_modules_base64-js_index_js-node_modules_sanitize-html_index_js":"c79fc2bcc3f69676beda","packages_base-manager_lib_index_js":"13941b475e5b5a1ab133","vendors-node_modules_backbone_backbone_js-node_modules_lodash_isEqual_js":"ab532abfb6832a50f59e","webpack_sharing_consume_default_lumino_messaging":"5d5527fdcb258b5816a0","packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery":"5dd13f8e980fa3c50bfe","vendors-node_modules_d3-color_src_color_js-node_modules_d3-format_src_defaultLocale_js-node_m-09b215":"2643c43f22ad111f4f82","packages_controls_lib_index_js":"6d5d05e0ec5e0c9f8185","packages_controls_lib_version_js":"105c0edf5e497d01e132","packages_output_lib_index_js":"49c9e4037a3b9e9e3b18","vendors-node_modules_jquery_dist_jquery_js":"3afdf979a9caffdf15a7","vendors-node_modules_semver_index_js":"7a90c08b9ce3dac425c6","@jupyter-widgets/controls":"9535f14e043c94385218"}[chunkId] + ".js";
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/global */
/******/ 	(() => {
/******/ 		__webpack_require__.g = (function() {
/******/ 			if (typeof globalThis === 'object') return globalThis;
/******/ 			try {
/******/ 				return this || new Function('return this')();
/******/ 			} catch (e) {
/******/ 				if (typeof window === 'object') return window;
/******/ 			}
/******/ 		})();
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	(() => {
/******/ 		__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/load script */
/******/ 	(() => {
/******/ 		var inProgress = {};
/******/ 		var dataWebpackPrefix = "@jupyter-widgets/jupyterlab-manager:";
/******/ 		// loadScript function to load a script via script tag
/******/ 		__webpack_require__.l = (url, done, key, chunkId) => {
/******/ 			if(inProgress[url]) { inProgress[url].push(done); return; }
/******/ 			var script, needAttach;
/******/ 			if(key !== undefined) {
/******/ 				var scripts = document.getElementsByTagName("script");
/******/ 				for(var i = 0; i < scripts.length; i++) {
/******/ 					var s = scripts[i];
/******/ 					if(s.getAttribute("src") == url || s.getAttribute("data-webpack") == dataWebpackPrefix + key) { script = s; break; }
/******/ 				}
/******/ 			}
/******/ 			if(!script) {
/******/ 				needAttach = true;
/******/ 				script = document.createElement('script');
/******/ 		
/******/ 				script.charset = 'utf-8';
/******/ 				script.timeout = 120;
/******/ 				if (__webpack_require__.nc) {
/******/ 					script.setAttribute("nonce", __webpack_require__.nc);
/******/ 				}
/******/ 				script.setAttribute("data-webpack", dataWebpackPrefix + key);
/******/ 		
/******/ 				script.src = url;
/******/ 			}
/******/ 			inProgress[url] = [done];
/******/ 			var onScriptComplete = (prev, event) => {
/******/ 				// avoid mem leaks in IE.
/******/ 				script.onerror = script.onload = null;
/******/ 				clearTimeout(timeout);
/******/ 				var doneFns = inProgress[url];
/******/ 				delete inProgress[url];
/******/ 				script.parentNode && script.parentNode.removeChild(script);
/******/ 				doneFns && doneFns.forEach((fn) => (fn(event)));
/******/ 				if(prev) return prev(event);
/******/ 			}
/******/ 			var timeout = setTimeout(onScriptComplete.bind(null, undefined, { type: 'timeout', target: script }), 120000);
/******/ 			script.onerror = onScriptComplete.bind(null, script.onerror);
/******/ 			script.onload = onScriptComplete.bind(null, script.onload);
/******/ 			needAttach && document.head.appendChild(script);
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/node module decorator */
/******/ 	(() => {
/******/ 		__webpack_require__.nmd = (module) => {
/******/ 			module.paths = [];
/******/ 			if (!module.children) module.children = [];
/******/ 			return module;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/sharing */
/******/ 	(() => {
/******/ 		__webpack_require__.S = {};
/******/ 		var initPromises = {};
/******/ 		var initTokens = {};
/******/ 		__webpack_require__.I = (name, initScope) => {
/******/ 			if(!initScope) initScope = [];
/******/ 			// handling circular init calls
/******/ 			var initToken = initTokens[name];
/******/ 			if(!initToken) initToken = initTokens[name] = {};
/******/ 			if(initScope.indexOf(initToken) >= 0) return;
/******/ 			initScope.push(initToken);
/******/ 			// only runs once
/******/ 			if(initPromises[name]) return initPromises[name];
/******/ 			// creates a new share scope if needed
/******/ 			if(!__webpack_require__.o(__webpack_require__.S, name)) __webpack_require__.S[name] = {};
/******/ 			// runs all init snippets from all modules reachable
/******/ 			var scope = __webpack_require__.S[name];
/******/ 			var warn = (msg) => {
/******/ 				if (typeof console !== "undefined" && console.warn) console.warn(msg);
/******/ 			};
/******/ 			var uniqueName = "@jupyter-widgets/jupyterlab-manager";
/******/ 			var register = (name, version, factory, eager) => {
/******/ 				var versions = scope[name] = scope[name] || {};
/******/ 				var activeVersion = versions[version];
/******/ 				if(!activeVersion || (!activeVersion.loaded && (!eager != !activeVersion.eager ? eager : uniqueName > activeVersion.from))) versions[version] = { get: factory, from: uniqueName, eager: !!eager };
/******/ 			};
/******/ 			var initExternal = (id) => {
/******/ 				var handleError = (err) => (warn("Initialization of sharing external failed: " + err));
/******/ 				try {
/******/ 					var module = __webpack_require__(id);
/******/ 					if(!module) return;
/******/ 					var initFn = (module) => (module && module.init && module.init(__webpack_require__.S[name], initScope))
/******/ 					if(module.then) return promises.push(module.then(initFn, handleError));
/******/ 					var initResult = initFn(module);
/******/ 					if(initResult && initResult.then) return promises.push(initResult['catch'](handleError));
/******/ 				} catch(err) { handleError(err); }
/******/ 			}
/******/ 			var promises = [];
/******/ 			switch(name) {
/******/ 				case "default": {
/******/ 					register("@jupyter-widgets/base-manager", "1.0.12", () => (Promise.all([__webpack_require__.e("vendors-node_modules_base64-js_index_js-node_modules_sanitize-html_index_js"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("packages_base-manager_lib_index_js")]).then(() => (() => (__webpack_require__(/*! ../../packages/base-manager/lib/index.js */ "../../packages/base-manager/lib/index.js"))))));
/******/ 					register("@jupyter-widgets/base", "6.0.11", () => (Promise.all([__webpack_require__.e("vendors-node_modules_backbone_backbone_js-node_modules_lodash_isEqual_js"), __webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_lumino_messaging"), __webpack_require__.e("packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery")]).then(() => (() => (__webpack_require__(/*! ../../packages/base/lib/index.js */ "../../packages/base/lib/index.js"))))));
/******/ 					register("@jupyter-widgets/controls", "5.0.12", () => (Promise.all([__webpack_require__.e("vendors-node_modules_d3-color_src_color_js-node_modules_d3-format_src_defaultLocale_js-node_m-09b215"), __webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling"), __webpack_require__.e("webpack_sharing_consume_default_lumino_messaging"), __webpack_require__.e("packages_controls_lib_index_js"), __webpack_require__.e("packages_controls_lib_version_js")]).then(() => (() => (__webpack_require__(/*! ../../packages/controls/lib/index.js */ "../../packages/controls/lib/index.js"))))));
/******/ 					register("@jupyter-widgets/jupyterlab-manager", "5.0.15", () => (Promise.all([__webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling"), __webpack_require__.e("lib_index_js")]).then(() => (() => (__webpack_require__(/*! ./lib/index.js */ "./lib/index.js"))))));
/******/ 					register("@jupyter-widgets/output", "6.0.11", () => (Promise.all([__webpack_require__.e("webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base"), __webpack_require__.e("packages_output_lib_index_js")]).then(() => (() => (__webpack_require__(/*! ../../packages/output/lib/index.js */ "../../packages/output/lib/index.js"))))));
/******/ 					register("jquery", "3.7.1", () => (__webpack_require__.e("vendors-node_modules_jquery_dist_jquery_js").then(() => (() => (__webpack_require__(/*! ../../node_modules/jquery/dist/jquery.js */ "../../node_modules/jquery/dist/jquery.js"))))));
/******/ 					register("semver", "7.6.3", () => (__webpack_require__.e("vendors-node_modules_semver_index_js").then(() => (() => (__webpack_require__(/*! ../../node_modules/semver/index.js */ "../../node_modules/semver/index.js"))))));
/******/ 				}
/******/ 				break;
/******/ 			}
/******/ 			if(!promises.length) return initPromises[name] = 1;
/******/ 			return initPromises[name] = Promise.all(promises).then(() => (initPromises[name] = 1));
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/publicPath */
/******/ 	(() => {
/******/ 		var scriptUrl;
/******/ 		if (__webpack_require__.g.importScripts) scriptUrl = __webpack_require__.g.location + "";
/******/ 		var document = __webpack_require__.g.document;
/******/ 		if (!scriptUrl && document) {
/******/ 			if (document.currentScript && document.currentScript.tagName.toUpperCase() === 'SCRIPT')
/******/ 				scriptUrl = document.currentScript.src;
/******/ 			if (!scriptUrl) {
/******/ 				var scripts = document.getElementsByTagName("script");
/******/ 				if(scripts.length) {
/******/ 					var i = scripts.length - 1;
/******/ 					while (i > -1 && (!scriptUrl || !/^http(s?):/.test(scriptUrl))) scriptUrl = scripts[i--].src;
/******/ 				}
/******/ 			}
/******/ 		}
/******/ 		// When supporting browsers where an automatic publicPath is not supported you must specify an output.publicPath manually via configuration
/******/ 		// or pass an empty string ("") and set the __webpack_public_path__ variable from your code to use your own logic.
/******/ 		if (!scriptUrl) throw new Error("Automatic publicPath is not supported in this browser");
/******/ 		scriptUrl = scriptUrl.replace(/#.*$/, "").replace(/\?.*$/, "").replace(/\/[^\/]+$/, "/");
/******/ 		__webpack_require__.p = scriptUrl;
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/consumes */
/******/ 	(() => {
/******/ 		var parseVersion = (str) => {
/******/ 			// see webpack/lib/util/semver.js for original code
/******/ 			var p=p=>{return p.split(".").map((p=>{return+p==p?+p:p}))},n=/^([^-+]+)?(?:-([^+]+))?(?:\+(.+))?$/.exec(str),r=n[1]?p(n[1]):[];return n[2]&&(r.length++,r.push.apply(r,p(n[2]))),n[3]&&(r.push([]),r.push.apply(r,p(n[3]))),r;
/******/ 		}
/******/ 		var versionLt = (a, b) => {
/******/ 			// see webpack/lib/util/semver.js for original code
/******/ 			a=parseVersion(a),b=parseVersion(b);for(var r=0;;){if(r>=a.length)return r<b.length&&"u"!=(typeof b[r])[0];var e=a[r],n=(typeof e)[0];if(r>=b.length)return"u"==n;var t=b[r],f=(typeof t)[0];if(n!=f)return"o"==n&&"n"==f||("s"==f||"u"==n);if("o"!=n&&"u"!=n&&e!=t)return e<t;r++}
/******/ 		}
/******/ 		var rangeToString = (range) => {
/******/ 			// see webpack/lib/util/semver.js for original code
/******/ 			var r=range[0],n="";if(1===range.length)return"*";if(r+.5){n+=0==r?">=":-1==r?"<":1==r?"^":2==r?"~":r>0?"=":"!=";for(var e=1,a=1;a<range.length;a++){e--,n+="u"==(typeof(t=range[a]))[0]?"-":(e>0?".":"")+(e=2,t)}return n}var g=[];for(a=1;a<range.length;a++){var t=range[a];g.push(0===t?"not("+o()+")":1===t?"("+o()+" || "+o()+")":2===t?g.pop()+" "+g.pop():rangeToString(t))}return o();function o(){return g.pop().replace(/^\((.+)\)$/,"$1")}
/******/ 		}
/******/ 		var satisfy = (range, version) => {
/******/ 			// see webpack/lib/util/semver.js for original code
/******/ 			if(0 in range){version=parseVersion(version);var e=range[0],r=e<0;r&&(e=-e-1);for(var n=0,i=1,a=!0;;i++,n++){var f,s,g=i<range.length?(typeof range[i])[0]:"";if(n>=version.length||"o"==(s=(typeof(f=version[n]))[0]))return!a||("u"==g?i>e&&!r:""==g!=r);if("u"==s){if(!a||"u"!=g)return!1}else if(a)if(g==s)if(i<=e){if(f!=range[i])return!1}else{if(r?f>range[i]:f<range[i])return!1;f!=range[i]&&(a=!1)}else if("s"!=g&&"n"!=g){if(r||i<=e)return!1;a=!1,i--}else{if(i<=e||s<g!=r)return!1;a=!1}else"s"!=g&&"n"!=g&&(a=!1,i--)}}var t=[],o=t.pop.bind(t);for(n=1;n<range.length;n++){var u=range[n];t.push(1==u?o()|o():2==u?o()&o():u?satisfy(u,version):!o())}return!!o();
/******/ 		}
/******/ 		var exists = (scope, key) => {
/******/ 			return scope && __webpack_require__.o(scope, key);
/******/ 		}
/******/ 		var get = (entry) => {
/******/ 			entry.loaded = 1;
/******/ 			return entry.get()
/******/ 		};
/******/ 		var eagerOnly = (versions) => {
/******/ 			return Object.keys(versions).reduce((filtered, version) => {
/******/ 					if (versions[version].eager) {
/******/ 						filtered[version] = versions[version];
/******/ 					}
/******/ 					return filtered;
/******/ 			}, {});
/******/ 		};
/******/ 		var findLatestVersion = (scope, key, eager) => {
/******/ 			var versions = eager ? eagerOnly(scope[key]) : scope[key];
/******/ 			var key = Object.keys(versions).reduce((a, b) => {
/******/ 				return !a || versionLt(a, b) ? b : a;
/******/ 			}, 0);
/******/ 			return key && versions[key];
/******/ 		};
/******/ 		var findSatisfyingVersion = (scope, key, requiredVersion, eager) => {
/******/ 			var versions = eager ? eagerOnly(scope[key]) : scope[key];
/******/ 			var key = Object.keys(versions).reduce((a, b) => {
/******/ 				if (!satisfy(requiredVersion, b)) return a;
/******/ 				return !a || versionLt(a, b) ? b : a;
/******/ 			}, 0);
/******/ 			return key && versions[key]
/******/ 		};
/******/ 		var findSingletonVersionKey = (scope, key, eager) => {
/******/ 			var versions = eager ? eagerOnly(scope[key]) : scope[key];
/******/ 			return Object.keys(versions).reduce((a, b) => {
/******/ 				return !a || (!versions[a].loaded && versionLt(a, b)) ? b : a;
/******/ 			}, 0);
/******/ 		};
/******/ 		var getInvalidSingletonVersionMessage = (scope, key, version, requiredVersion) => {
/******/ 			return "Unsatisfied version " + version + " from " + (version && scope[key][version].from) + " of shared singleton module " + key + " (required " + rangeToString(requiredVersion) + ")"
/******/ 		};
/******/ 		var getInvalidVersionMessage = (scope, scopeName, key, requiredVersion, eager) => {
/******/ 			var versions = scope[key];
/******/ 			return "No satisfying version (" + rangeToString(requiredVersion) + ")" + (eager ? " for eager consumption" : "") + " of shared module " + key + " found in shared scope " + scopeName + ".\n" +
/******/ 				"Available versions: " + Object.keys(versions).map((key) => {
/******/ 				return key + " from " + versions[key].from;
/******/ 			}).join(", ");
/******/ 		};
/******/ 		var fail = (msg) => {
/******/ 			throw new Error(msg);
/******/ 		}
/******/ 		var failAsNotExist = (scopeName, key) => {
/******/ 			return fail("Shared module " + key + " doesn't exist in shared scope " + scopeName);
/******/ 		}
/******/ 		var warn = /*#__PURE__*/ (msg) => {
/******/ 			if (typeof console !== "undefined" && console.warn) console.warn(msg);
/******/ 		};
/******/ 		var init = (fn) => (function(scopeName, key, eager, c, d) {
/******/ 			var promise = __webpack_require__.I(scopeName);
/******/ 			if (promise && promise.then && !eager) {
/******/ 				return promise.then(fn.bind(fn, scopeName, __webpack_require__.S[scopeName], key, false, c, d));
/******/ 			}
/******/ 			return fn(scopeName, __webpack_require__.S[scopeName], key, eager, c, d);
/******/ 		});
/******/ 		
/******/ 		var useFallback = (scopeName, key, fallback) => {
/******/ 			return fallback ? fallback() : failAsNotExist(scopeName, key);
/******/ 		}
/******/ 		var load = /*#__PURE__*/ init((scopeName, scope, key, eager, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			return get(findLatestVersion(scope, key, eager));
/******/ 		});
/******/ 		var loadVersion = /*#__PURE__*/ init((scopeName, scope, key, eager, requiredVersion, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			var satisfyingVersion = findSatisfyingVersion(scope, key, requiredVersion, eager);
/******/ 			if (satisfyingVersion) return get(satisfyingVersion);
/******/ 			warn(getInvalidVersionMessage(scope, scopeName, key, requiredVersion, eager))
/******/ 			return get(findLatestVersion(scope, key, eager));
/******/ 		});
/******/ 		var loadStrictVersion = /*#__PURE__*/ init((scopeName, scope, key, eager, requiredVersion, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			var satisfyingVersion = findSatisfyingVersion(scope, key, requiredVersion, eager);
/******/ 			if (satisfyingVersion) return get(satisfyingVersion);
/******/ 			if (fallback) return fallback();
/******/ 			fail(getInvalidVersionMessage(scope, scopeName, key, requiredVersion, eager));
/******/ 		});
/******/ 		var loadSingleton = /*#__PURE__*/ init((scopeName, scope, key, eager, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			var version = findSingletonVersionKey(scope, key, eager);
/******/ 			return get(scope[key][version]);
/******/ 		});
/******/ 		var loadSingletonVersion = /*#__PURE__*/ init((scopeName, scope, key, eager, requiredVersion, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			var version = findSingletonVersionKey(scope, key, eager);
/******/ 			if (!satisfy(requiredVersion, version)) {
/******/ 				warn(getInvalidSingletonVersionMessage(scope, key, version, requiredVersion));
/******/ 			}
/******/ 			return get(scope[key][version]);
/******/ 		});
/******/ 		var loadStrictSingletonVersion = /*#__PURE__*/ init((scopeName, scope, key, eager, requiredVersion, fallback) => {
/******/ 			if (!exists(scope, key)) return useFallback(scopeName, key, fallback);
/******/ 			var version = findSingletonVersionKey(scope, key, eager);
/******/ 			if (!satisfy(requiredVersion, version)) {
/******/ 				fail(getInvalidSingletonVersionMessage(scope, key, version, requiredVersion));
/******/ 			}
/******/ 			return get(scope[key][version]);
/******/ 		});
/******/ 		var installedModules = {};
/******/ 		var moduleToHandlerMapping = {
/******/ 			"webpack/sharing/consume/default/@lumino/widgets": () => (loadSingletonVersion("default", "@lumino/widgets", false, [1,2,3,1,,"alpha",0])),
/******/ 			"webpack/sharing/consume/default/jquery/jquery?123b": () => (loadStrictVersion("default", "jquery", false, [1,3,1,1], () => (__webpack_require__.e("vendors-node_modules_jquery_dist_jquery_js").then(() => (() => (__webpack_require__(/*! jquery */ "../../node_modules/jquery/dist/jquery.js"))))))),
/******/ 			"webpack/sharing/consume/default/@lumino/coreutils": () => (loadSingletonVersion("default", "@lumino/coreutils", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base": () => (loadStrictVersion("default", "@jupyter-widgets/base", false, [1,6,0,11], () => (Promise.all([__webpack_require__.e("vendors-node_modules_backbone_backbone_js-node_modules_lodash_isEqual_js"), __webpack_require__.e("webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery"), __webpack_require__.e("webpack_sharing_consume_default_lumino_coreutils"), __webpack_require__.e("webpack_sharing_consume_default_lumino_messaging"), __webpack_require__.e("packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery")]).then(() => (() => (__webpack_require__(/*! @jupyter-widgets/base */ "../../packages/base/lib/index.js"))))))),
/******/ 			"webpack/sharing/consume/default/@lumino/algorithm": () => (loadSingletonVersion("default", "@lumino/algorithm", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/@lumino/signaling": () => (loadSingletonVersion("default", "@lumino/signaling", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/settingregistry": () => (loadSingletonVersion("default", "@jupyterlab/settingregistry", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/console": () => (loadSingletonVersion("default", "@jupyterlab/console", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/notebook": () => (loadSingletonVersion("default", "@jupyterlab/notebook", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/mainmenu": () => (loadSingletonVersion("default", "@jupyterlab/mainmenu", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/rendermime": () => (loadSingletonVersion("default", "@jupyterlab/rendermime", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/logconsole": () => (loadSingletonVersion("default", "@jupyterlab/logconsole", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@lumino/disposable": () => (loadSingletonVersion("default", "@lumino/disposable", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/@jupyter-widgets/base-manager/@jupyter-widgets/base-manager": () => (loadStrictVersion("default", "@jupyter-widgets/base-manager", false, [1,1,0,12], () => (Promise.all([__webpack_require__.e("vendors-node_modules_base64-js_index_js-node_modules_sanitize-html_index_js"), __webpack_require__.e("packages_base-manager_lib_index_js")]).then(() => (() => (__webpack_require__(/*! @jupyter-widgets/base-manager */ "../../packages/base-manager/lib/index.js"))))))),
/******/ 			"webpack/sharing/consume/default/semver/semver": () => (loadStrictVersion("default", "semver", false, [1,7,3,5], () => (__webpack_require__.e("vendors-node_modules_semver_index_js").then(() => (() => (__webpack_require__(/*! semver */ "../../node_modules/semver/index.js"))))))),
/******/ 			"webpack/sharing/consume/default/@jupyter-widgets/output/@jupyter-widgets/output": () => (loadStrictVersion("default", "@jupyter-widgets/output", false, [1,6,0,11], () => (__webpack_require__.e("packages_output_lib_index_js").then(() => (() => (__webpack_require__(/*! @jupyter-widgets/output */ "../../packages/output/lib/index.js"))))))),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/outputarea": () => (loadVersion("default", "@jupyterlab/outputarea", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/services": () => (loadSingletonVersion("default", "@jupyterlab/services", false, [1,7,4,6])),
/******/ 			"webpack/sharing/consume/default/@jupyterlab/translation": () => (loadSingletonVersion("default", "@jupyterlab/translation", false, [1,4,4,6])),
/******/ 			"webpack/sharing/consume/default/@lumino/messaging": () => (loadSingletonVersion("default", "@lumino/messaging", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/jquery/jquery?5edd": () => (load("default", "jquery", false, () => (__webpack_require__.e("vendors-node_modules_jquery_dist_jquery_js").then(() => (() => (__webpack_require__(/*! jquery */ "../../node_modules/jquery/dist/jquery.js"))))))),
/******/ 			"webpack/sharing/consume/default/@lumino/domutils": () => (loadSingletonVersion("default", "@lumino/domutils", false, [1,2,0,0])),
/******/ 			"webpack/sharing/consume/default/@jupyter-widgets/controls/@jupyter-widgets/controls": () => (loadStrictVersion("default", "@jupyter-widgets/controls", false, [1,5,0,12], () => (Promise.all([__webpack_require__.e("vendors-node_modules_d3-color_src_color_js-node_modules_d3-format_src_defaultLocale_js-node_m-09b215"), __webpack_require__.e("webpack_sharing_consume_default_lumino_messaging"), __webpack_require__.e("packages_controls_lib_index_js")]).then(() => (() => (__webpack_require__(/*! @jupyter-widgets/controls */ "../../packages/controls/lib/index.js")))))))
/******/ 		};
/******/ 		// no consumes in initial chunks
/******/ 		var chunkMapping = {
/******/ 			"webpack_sharing_consume_default_lumino_widgets-webpack_sharing_consume_default_jquery_jquery": [
/******/ 				"webpack/sharing/consume/default/@lumino/widgets",
/******/ 				"webpack/sharing/consume/default/jquery/jquery?123b"
/******/ 			],
/******/ 			"webpack_sharing_consume_default_lumino_coreutils": [
/******/ 				"webpack/sharing/consume/default/@lumino/coreutils"
/******/ 			],
/******/ 			"webpack_sharing_consume_default_jupyter-widgets_base_jupyter-widgets_base": [
/******/ 				"webpack/sharing/consume/default/@jupyter-widgets/base/@jupyter-widgets/base"
/******/ 			],
/******/ 			"webpack_sharing_consume_default_lumino_algorithm-webpack_sharing_consume_default_lumino_signaling": [
/******/ 				"webpack/sharing/consume/default/@lumino/algorithm",
/******/ 				"webpack/sharing/consume/default/@lumino/signaling"
/******/ 			],
/******/ 			"lib_index_js": [
/******/ 				"webpack/sharing/consume/default/@jupyterlab/settingregistry",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/console",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/notebook",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/mainmenu",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/rendermime",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/logconsole",
/******/ 				"webpack/sharing/consume/default/@lumino/disposable",
/******/ 				"webpack/sharing/consume/default/@jupyter-widgets/base-manager/@jupyter-widgets/base-manager",
/******/ 				"webpack/sharing/consume/default/semver/semver",
/******/ 				"webpack/sharing/consume/default/@jupyter-widgets/output/@jupyter-widgets/output",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/outputarea",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/services",
/******/ 				"webpack/sharing/consume/default/@jupyterlab/translation"
/******/ 			],
/******/ 			"webpack_sharing_consume_default_lumino_messaging": [
/******/ 				"webpack/sharing/consume/default/@lumino/messaging"
/******/ 			],
/******/ 			"packages_base_lib_index_js-webpack_sharing_consume_default_jquery_jquery": [
/******/ 				"webpack/sharing/consume/default/jquery/jquery?5edd"
/******/ 			],
/******/ 			"packages_controls_lib_index_js": [
/******/ 				"webpack/sharing/consume/default/@lumino/domutils"
/******/ 			],
/******/ 			"@jupyter-widgets/controls": [
/******/ 				"webpack/sharing/consume/default/@jupyter-widgets/controls/@jupyter-widgets/controls"
/******/ 			]
/******/ 		};
/******/ 		var startedInstallModules = {};
/******/ 		__webpack_require__.f.consumes = (chunkId, promises) => {
/******/ 			if(__webpack_require__.o(chunkMapping, chunkId)) {
/******/ 				chunkMapping[chunkId].forEach((id) => {
/******/ 					if(__webpack_require__.o(installedModules, id)) return promises.push(installedModules[id]);
/******/ 					if(!startedInstallModules[id]) {
/******/ 					var onFactory = (factory) => {
/******/ 						installedModules[id] = 0;
/******/ 						__webpack_require__.m[id] = (module) => {
/******/ 							delete __webpack_require__.c[id];
/******/ 							module.exports = factory();
/******/ 						}
/******/ 					};
/******/ 					startedInstallModules[id] = true;
/******/ 					var onError = (error) => {
/******/ 						delete installedModules[id];
/******/ 						__webpack_require__.m[id] = (module) => {
/******/ 							delete __webpack_require__.c[id];
/******/ 							throw error;
/******/ 						}
/******/ 					};
/******/ 					try {
/******/ 						var promise = moduleToHandlerMapping[id]();
/******/ 						if(promise.then) {
/******/ 							promises.push(installedModules[id] = promise.then(onFactory)['catch'](onError));
/******/ 						} else onFactory(promise);
/******/ 					} catch(e) { onError(e); }
/******/ 					}
/******/ 				});
/******/ 			}
/******/ 		}
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/jsonp chunk loading */
/******/ 	(() => {
/******/ 		__webpack_require__.b = document.baseURI || self.location.href;
/******/ 		
/******/ 		// object to store loaded and loading chunks
/******/ 		// undefined = chunk not loaded, null = chunk preloaded/prefetched
/******/ 		// [resolve, reject, Promise] = chunk loading, 0 = chunk loaded
/******/ 		var installedChunks = {
/******/ 			"@jupyter-widgets/jupyterlab-manager": 0
/******/ 		};
/******/ 		
/******/ 		__webpack_require__.f.j = (chunkId, promises) => {
/******/ 				// JSONP chunk loading for javascript
/******/ 				var installedChunkData = __webpack_require__.o(installedChunks, chunkId) ? installedChunks[chunkId] : undefined;
/******/ 				if(installedChunkData !== 0) { // 0 means "already installed".
/******/ 		
/******/ 					// a Promise means "currently loading".
/******/ 					if(installedChunkData) {
/******/ 						promises.push(installedChunkData[2]);
/******/ 					} else {
/******/ 						if(!/^(webpack_sharing_consume_default_(lumino_((algorithm\-webpack_sharing_consume_default_lumino_signal|messag)ing|coreutils|widgets\-webpack_sharing_consume_default_jquery_jquery)|jupyter\-widgets_base_jupyter\-widgets_base)|@jupyter\-widgets\/controls)$/.test(chunkId)) {
/******/ 							// setup Promise in chunk cache
/******/ 							var promise = new Promise((resolve, reject) => (installedChunkData = installedChunks[chunkId] = [resolve, reject]));
/******/ 							promises.push(installedChunkData[2] = promise);
/******/ 		
/******/ 							// start chunk loading
/******/ 							var url = __webpack_require__.p + __webpack_require__.u(chunkId);
/******/ 							// create error before stack unwound to get useful stacktrace later
/******/ 							var error = new Error();
/******/ 							var loadingEnded = (event) => {
/******/ 								if(__webpack_require__.o(installedChunks, chunkId)) {
/******/ 									installedChunkData = installedChunks[chunkId];
/******/ 									if(installedChunkData !== 0) installedChunks[chunkId] = undefined;
/******/ 									if(installedChunkData) {
/******/ 										var errorType = event && (event.type === 'load' ? 'missing' : event.type);
/******/ 										var realSrc = event && event.target && event.target.src;
/******/ 										error.message = 'Loading chunk ' + chunkId + ' failed.\n(' + errorType + ': ' + realSrc + ')';
/******/ 										error.name = 'ChunkLoadError';
/******/ 										error.type = errorType;
/******/ 										error.request = realSrc;
/******/ 										installedChunkData[1](error);
/******/ 									}
/******/ 								}
/******/ 							};
/******/ 							__webpack_require__.l(url, loadingEnded, "chunk-" + chunkId, chunkId);
/******/ 						} else installedChunks[chunkId] = 0;
/******/ 					}
/******/ 				}
/******/ 		};
/******/ 		
/******/ 		// no prefetching
/******/ 		
/******/ 		// no preloaded
/******/ 		
/******/ 		// no HMR
/******/ 		
/******/ 		// no HMR manifest
/******/ 		
/******/ 		// no on chunks loaded
/******/ 		
/******/ 		// install a JSONP callback for chunk loading
/******/ 		var webpackJsonpCallback = (parentChunkLoadingFunction, data) => {
/******/ 			var [chunkIds, moreModules, runtime] = data;
/******/ 			// add "moreModules" to the modules object,
/******/ 			// then flag all "chunkIds" as loaded and fire callback
/******/ 			var moduleId, chunkId, i = 0;
/******/ 			if(chunkIds.some((id) => (installedChunks[id] !== 0))) {
/******/ 				for(moduleId in moreModules) {
/******/ 					if(__webpack_require__.o(moreModules, moduleId)) {
/******/ 						__webpack_require__.m[moduleId] = moreModules[moduleId];
/******/ 					}
/******/ 				}
/******/ 				if(runtime) var result = runtime(__webpack_require__);
/******/ 			}
/******/ 			if(parentChunkLoadingFunction) parentChunkLoadingFunction(data);
/******/ 			for(;i < chunkIds.length; i++) {
/******/ 				chunkId = chunkIds[i];
/******/ 				if(__webpack_require__.o(installedChunks, chunkId) && installedChunks[chunkId]) {
/******/ 					installedChunks[chunkId][0]();
/******/ 				}
/******/ 				installedChunks[chunkId] = 0;
/******/ 			}
/******/ 		
/******/ 		}
/******/ 		
/******/ 		var chunkLoadingGlobal = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] = self["webpackChunk_jupyter_widgets_jupyterlab_manager"] || [];
/******/ 		chunkLoadingGlobal.forEach(webpackJsonpCallback.bind(null, 0));
/******/ 		chunkLoadingGlobal.push = webpackJsonpCallback.bind(null, chunkLoadingGlobal.push.bind(chunkLoadingGlobal));
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/nonce */
/******/ 	(() => {
/******/ 		__webpack_require__.nc = undefined;
/******/ 	})();
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// module cache are used so entry inlining is disabled
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	var __webpack_exports__ = __webpack_require__("webpack/container/entry/@jupyter-widgets/jupyterlab-manager");
/******/ 	(_JUPYTERLAB = typeof _JUPYTERLAB === "undefined" ? {} : _JUPYTERLAB)["@jupyter-widgets/jupyterlab-manager"] = __webpack_exports__;
/******/ 	
/******/ })()
;
//# sourceMappingURL=remoteEntry.9077b3d2deaffb329dfc.js.map