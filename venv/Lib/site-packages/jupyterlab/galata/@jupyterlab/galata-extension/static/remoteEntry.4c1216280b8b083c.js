var _JUPYTERLAB;
(() => {
var __webpack_modules__ = ({
"../../node_modules/process/browser.js"(module) {
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


},
"@module-federation/runtime/rspack.js!=!data:text/javascript,import%20__module_federation_bundler_runtime__%20from%20%22%2Fhome%2Frunner%2Fwork%2Fjupyterlab%2Fjupyterlab%2F.jupyter_releaser_checkout%2Fnode_modules%2F%40module-federation%2Fwebpack-bundler-runtime%2Fdist%2Findex.cjs%22%3Bconst%20__module_federation_runtime_plugins__%20%3D%20%5B%5D.filter((%7B%20plugin%20%7D)%20%3D%3E%20plugin).map((%7B%20plugin%2C%20params%20%7D)%20%3D%3E%20plugin(params))%3Bconst%20__module_federation_remote_infos__%20%3D%20%7B%7D%3Bconst%20__module_federation_container_name__%20%3D%20%22%40jupyterlab%2Fgalata-extension%22%3Bconst%20__module_federation_share_strategy__%20%3D%20%22version-first%22%3Bconst%20__module_federation_share_fallbacks__%20%3D%20undefined%3Bconst%20__module_federation_library_type__%20%3D%20%22var%22%3Bif((__webpack_require__.initializeSharingData%7C%7C__webpack_require__.initializeExposesData)%26%26__webpack_require__.federation)%7Bvar%20_ref%2C_ref1%2C_ref2%2C_ref3%2C_ref4%3Bvar%20__webpack_require___remotesLoadingData%2C__webpack_require___remotesLoadingData1%2C__webpack_require___initializeSharingData%2C__webpack_require___consumesLoadingData%2C__webpack_require___consumesLoadingData1%2C__webpack_require___initializeExposesData%2C__webpack_require___consumesLoadingData2%3Bconst%20override%3D(obj%2Ckey%2Cvalue)%3D%3E%7Bif(!obj)return%3Bif(obj%5Bkey%5D)obj%5Bkey%5D%3Dvalue%7D%3Bconst%20merge%3D(obj%2Ckey%2Cfn)%3D%3E%7Bconst%20value%3Dfn()%3Bif(Array.isArray(value))%7Bvar%20_obj%2C_key%2C_%3B(_%3D(_obj%3Dobj)%5B_key%3Dkey%5D)!%3D%3Dnull%26%26_!%3D%3Dvoid%200%3F_%3A_obj%5B_key%5D%3D%5B%5D%3Bobj%5Bkey%5D.push(...value)%7Delse%20if(typeof%20value%3D%3D%3D%22object%22%26%26value!%3D%3Dnull)%7Bvar%20_obj1%2C_key1%2C_1%3B(_1%3D(_obj1%3Dobj)%5B_key1%3Dkey%5D)!%3D%3Dnull%26%26_1!%3D%3Dvoid%200%3F_1%3A_obj1%5B_key1%5D%3D%7B%7D%3BObject.assign(obj%5Bkey%5D%2Cvalue)%7D%7D%3Bconst%20early%3D(obj%2Ckey%2Cinitial)%3D%3E%7Bvar%20_obj%2C_key%2C_%3B(_%3D(_obj%3Dobj)%5B_key%3Dkey%5D)!%3D%3Dnull%26%26_!%3D%3Dvoid%200%3F_%3A_obj%5B_key%5D%3Dinitial()%7D%3Bconst%20remotesLoadingChunkMapping%3D(_ref%3D(__webpack_require___remotesLoadingData%3D__webpack_require__.remotesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___remotesLoadingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___remotesLoadingData.chunkMapping)!%3D%3Dnull%26%26_ref!%3D%3Dvoid%200%3F_ref%3A%7B%7D%3Bconst%20remotesLoadingModuleIdToRemoteDataMapping%3D(_ref1%3D(__webpack_require___remotesLoadingData1%3D__webpack_require__.remotesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___remotesLoadingData1%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___remotesLoadingData1.moduleIdToRemoteDataMapping)!%3D%3Dnull%26%26_ref1!%3D%3Dvoid%200%3F_ref1%3A%7B%7D%3Bconst%20initializeSharingScopeToInitDataMapping%3D(_ref2%3D(__webpack_require___initializeSharingData%3D__webpack_require__.initializeSharingData)%3D%3D%3Dnull%7C%7C__webpack_require___initializeSharingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___initializeSharingData.scopeToSharingDataMapping)!%3D%3Dnull%26%26_ref2!%3D%3Dvoid%200%3F_ref2%3A%7B%7D%3Bconst%20consumesLoadingChunkMapping%3D(_ref3%3D(__webpack_require___consumesLoadingData%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData.chunkMapping)!%3D%3Dnull%26%26_ref3!%3D%3Dvoid%200%3F_ref3%3A%7B%7D%3Bconst%20consumesLoadingModuleToConsumeDataMapping%3D(_ref4%3D(__webpack_require___consumesLoadingData1%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData1%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData1.moduleIdToConsumeDataMapping)!%3D%3Dnull%26%26_ref4!%3D%3Dvoid%200%3F_ref4%3A%7B%7D%3Bconst%20consumesLoadinginstalledModules%3D%7B%7D%3Bconst%20initializeSharingInitPromises%3D%5B%5D%3Bconst%20initializeSharingInitTokens%3D%7B%7D%3Bconst%20containerShareScope%3D(__webpack_require___initializeExposesData%3D__webpack_require__.initializeExposesData)%3D%3D%3Dnull%7C%7C__webpack_require___initializeExposesData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___initializeExposesData.shareScope%3Bfor(const%20key%20in%20__module_federation_bundler_runtime__)%7B__webpack_require__.federation%5Bkey%5D%3D__module_federation_bundler_runtime__%5Bkey%5D%7Dearly(__webpack_require__.federation%2C%22libraryType%22%2C()%3D%3E__module_federation_library_type__)%3Bearly(__webpack_require__.federation%2C%22sharedFallback%22%2C()%3D%3E__module_federation_share_fallbacks__)%3Bconst%20sharedFallback%3D__webpack_require__.federation.sharedFallback%3Bearly(__webpack_require__.federation%2C%22consumesLoadingModuleToHandlerMapping%22%2C()%3D%3E%7Bconst%20consumesLoadingModuleToHandlerMapping%3D%7B%7D%3Bfor(let%5BmoduleId%2Cdata%5Dof%20Object.entries(consumesLoadingModuleToConsumeDataMapping))%7Bvar%20__webpack_require___federation_bundlerRuntime%3BconsumesLoadingModuleToHandlerMapping%5BmoduleId%5D%3D%7Bgetter%3AsharedFallback%3F(__webpack_require___federation_bundlerRuntime%3D__webpack_require__.federation.bundlerRuntime)%3D%3D%3Dnull%7C%7C__webpack_require___federation_bundlerRuntime%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___federation_bundlerRuntime.getSharedFallbackGetter(%7BshareKey%3Adata.shareKey%2Cfactory%3Adata.fallback%2CwebpackRequire%3A__webpack_require__%2ClibraryType%3A__webpack_require__.federation.libraryType%7D)%3Adata.fallback%2CtreeShakingGetter%3AsharedFallback%3Fdata.fallback%3Aundefined%2CshareInfo%3A%7BshareConfig%3A%7BfixedDependencies%3Afalse%2CrequiredVersion%3Adata.requiredVersion%2CstrictVersion%3Adata.strictVersion%2Csingleton%3Adata.singleton%2Ceager%3Adata.eager%7D%2Cscope%3A%5Bdata.shareScope%5D%7D%2CshareKey%3Adata.shareKey%2CtreeShaking%3A__webpack_require__.federation.sharedFallback%3F%7Bget%3Adata.fallback%2Cmode%3Adata.treeShakingMode%7D%3Aundefined%7D%7Dreturn%20consumesLoadingModuleToHandlerMapping%7D)%3Bearly(__webpack_require__.federation%2C%22initOptions%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.initOptions%2C%22name%22%2C()%3D%3E__module_federation_container_name__)%3Bearly(__webpack_require__.federation.initOptions%2C%22shareStrategy%22%2C()%3D%3E__module_federation_share_strategy__)%3Bearly(__webpack_require__.federation.initOptions%2C%22shared%22%2C()%3D%3E%7Bconst%20shared%3D%7B%7D%3Bfor(let%5Bscope%2Cstages%5Dof%20Object.entries(initializeSharingScopeToInitDataMapping))%7Bfor(let%20stage%20of%20stages)%7Bif(typeof%20stage%3D%3D%3D%22object%22%26%26stage!%3D%3Dnull)%7Bconst%7Bname%2Cversion%2Cfactory%2Ceager%2Csingleton%2CrequiredVersion%2CstrictVersion%2CtreeShakingMode%7D%3Dstage%3Bconst%20shareConfig%3D%7B%7D%3Bconst%20isValidValue%3Dfunction(val)%7Breturn%20typeof%20val!%3D%3D%22undefined%22%7D%3Bif(isValidValue(singleton))%7BshareConfig.singleton%3Dsingleton%7Dif(isValidValue(requiredVersion))%7BshareConfig.requiredVersion%3DrequiredVersion%7Dif(isValidValue(eager))%7BshareConfig.eager%3Deager%7Dif(isValidValue(strictVersion))%7BshareConfig.strictVersion%3DstrictVersion%7Dconst%20options%3D%7Bversion%2Cscope%3A%5Bscope%5D%2CshareConfig%2Cget%3Afactory%2CtreeShaking%3AtreeShakingMode%3F%7Bmode%3AtreeShakingMode%7D%3Aundefined%7D%3Bif(shared%5Bname%5D)%7Bshared%5Bname%5D.push(options)%7Delse%7Bshared%5Bname%5D%3D%5Boptions%5D%7D%7D%7D%7Dreturn%20shared%7D)%3Bmerge(__webpack_require__.federation.initOptions%2C%22remotes%22%2C()%3D%3EObject.values(__module_federation_remote_infos__).flat().filter(remote%3D%3Eremote.externalType%3D%3D%3D%22script%22))%3Bmerge(__webpack_require__.federation.initOptions%2C%22plugins%22%2C()%3D%3E__module_federation_runtime_plugins__)%3Bearly(__webpack_require__.federation%2C%22bundlerRuntimeOptions%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions%2C%22remotes%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22chunkMapping%22%2C()%3D%3EremotesLoadingChunkMapping)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22remoteInfos%22%2C()%3D%3E__module_federation_remote_infos__)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22idToExternalAndNameMapping%22%2C()%3D%3E%7Bconst%20remotesLoadingIdToExternalAndNameMappingMapping%3D%7B%7D%3Bfor(let%5BmoduleId%2Cdata%5Dof%20Object.entries(remotesLoadingModuleIdToRemoteDataMapping))%7BremotesLoadingIdToExternalAndNameMappingMapping%5BmoduleId%5D%3D%5Bdata.shareScope%2Cdata.name%2Cdata.externalModuleId%2Cdata.remoteName%5D%7Dreturn%20remotesLoadingIdToExternalAndNameMappingMapping%7D)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22webpackRequire%22%2C()%3D%3E__webpack_require__)%3Bmerge(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22idToRemoteMap%22%2C()%3D%3E%7Bconst%20idToRemoteMap%3D%7B%7D%3Bfor(let%5Bid%2CremoteData%5Dof%20Object.entries(remotesLoadingModuleIdToRemoteDataMapping))%7Bconst%20info%3D__module_federation_remote_infos__%5BremoteData.remoteName%5D%3Bif(info)idToRemoteMap%5Bid%5D%3Dinfo%7Dreturn%20idToRemoteMap%7D)%3Boverride(__webpack_require__%2C%22S%22%2C__webpack_require__.federation.bundlerRuntime.S)%3Bif(__webpack_require__.federation.attachShareScopeMap)%7B__webpack_require__.federation.attachShareScopeMap(__webpack_require__)%7Doverride(__webpack_require__.f%2C%22remotes%22%2C(chunkId%2Cpromises)%3D%3E__webpack_require__.federation.bundlerRuntime.remotes(%7BchunkId%2Cpromises%2CchunkMapping%3AremotesLoadingChunkMapping%2CidToExternalAndNameMapping%3A__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToExternalAndNameMapping%2CidToRemoteMap%3A__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToRemoteMap%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__.f%2C%22consumes%22%2C(chunkId%2Cpromises)%3D%3E__webpack_require__.federation.bundlerRuntime.consumes(%7BchunkId%2Cpromises%2CchunkMapping%3AconsumesLoadingChunkMapping%2CmoduleToHandlerMapping%3A__webpack_require__.federation.consumesLoadingModuleToHandlerMapping%2CinstalledModules%3AconsumesLoadinginstalledModules%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22I%22%2C(name%2CinitScope)%3D%3E__webpack_require__.federation.bundlerRuntime.I(%7BshareScopeName%3Aname%2CinitScope%2CinitPromises%3AinitializeSharingInitPromises%2CinitTokens%3AinitializeSharingInitTokens%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22initContainer%22%2C(shareScope%2CinitScope%2CremoteEntryInitOptions)%3D%3E__webpack_require__.federation.bundlerRuntime.initContainerEntry(%7BshareScope%2CinitScope%2CremoteEntryInitOptions%2CshareScopeKey%3AcontainerShareScope%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22getContainer%22%2C(module%2CgetScope)%3D%3E%7Bvar%20moduleMap%3D__webpack_require__.initializeExposesData.moduleMap%3B__webpack_require__.R%3DgetScope%3BgetScope%3DObject.prototype.hasOwnProperty.call(moduleMap%2Cmodule)%3FmoduleMap%5Bmodule%5D()%3APromise.resolve().then(()%3D%3E%7Bthrow%20new%20Error('Module%20%22'%2Bmodule%2B'%22%20does%20not%20exist%20in%20container.')%7D)%3B__webpack_require__.R%3Dundefined%3Breturn%20getScope%7D)%3B__webpack_require__.federation.instance%3D__webpack_require__.federation.bundlerRuntime.init(%7BwebpackRequire%3A__webpack_require__%7D)%3Bif((__webpack_require___consumesLoadingData2%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData2%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData2.initialConsumes)%7B__webpack_require__.federation.bundlerRuntime.installInitialConsumes(%7BwebpackRequire%3A__webpack_require__%2CinstalledModules%3AconsumesLoadinginstalledModules%2CinitialConsumes%3A__webpack_require__.consumesLoadingData.initialConsumes%2CmoduleToHandlerMapping%3A__webpack_require__.federation.consumesLoadingModuleToHandlerMapping%7D)%7D%7D"(__unused_rspack_module, __webpack_exports__, __webpack_require__) {
"use strict";
__webpack_require__.r(__webpack_exports__);
/* import */ var _home_runner_work_jupyterlab_jupyterlab_jupyter_releaser_checkout_node_modules_module_federation_webpack_bundler_runtime_dist_index_cjs__rspack_import_0 = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/index.cjs");
/* import */ var _home_runner_work_jupyterlab_jupyterlab_jupyter_releaser_checkout_node_modules_module_federation_webpack_bundler_runtime_dist_index_cjs__rspack_import_0_default = /*#__PURE__*/__webpack_require__.n(_home_runner_work_jupyterlab_jupyterlab_jupyter_releaser_checkout_node_modules_module_federation_webpack_bundler_runtime_dist_index_cjs__rspack_import_0);
const __module_federation_runtime_plugins__ = [].filter(({ plugin }) => plugin).map(({ plugin, params }) => plugin(params));const __module_federation_remote_infos__ = {};const __module_federation_container_name__ = "@jupyterlab/galata-extension";const __module_federation_share_strategy__ = "version-first";const __module_federation_share_fallbacks__ = undefined;const __module_federation_library_type__ = "var";if((__webpack_require__.initializeSharingData||__webpack_require__.initializeExposesData)&&__webpack_require__.federation){var _ref,_ref1,_ref2,_ref3,_ref4;var __webpack_require___remotesLoadingData,__webpack_require___remotesLoadingData1,__webpack_require___initializeSharingData,__webpack_require___consumesLoadingData,__webpack_require___consumesLoadingData1,__webpack_require___initializeExposesData,__webpack_require___consumesLoadingData2;const override=(obj,key,value)=>{if(!obj)return;if(obj[key])obj[key]=value};const merge=(obj,key,fn)=>{const value=fn();if(Array.isArray(value)){var _obj,_key,_;(_=(_obj=obj)[_key=key])!==null&&_!==void 0?_:_obj[_key]=[];obj[key].push(...value)}else if(typeof value==="object"&&value!==null){var _obj1,_key1,_1;(_1=(_obj1=obj)[_key1=key])!==null&&_1!==void 0?_1:_obj1[_key1]={};Object.assign(obj[key],value)}};const early=(obj,key,initial)=>{var _obj,_key,_;(_=(_obj=obj)[_key=key])!==null&&_!==void 0?_:_obj[_key]=initial()};const remotesLoadingChunkMapping=(_ref=(__webpack_require___remotesLoadingData=__webpack_require__.remotesLoadingData)===null||__webpack_require___remotesLoadingData===void 0?void 0:__webpack_require___remotesLoadingData.chunkMapping)!==null&&_ref!==void 0?_ref:{};const remotesLoadingModuleIdToRemoteDataMapping=(_ref1=(__webpack_require___remotesLoadingData1=__webpack_require__.remotesLoadingData)===null||__webpack_require___remotesLoadingData1===void 0?void 0:__webpack_require___remotesLoadingData1.moduleIdToRemoteDataMapping)!==null&&_ref1!==void 0?_ref1:{};const initializeSharingScopeToInitDataMapping=(_ref2=(__webpack_require___initializeSharingData=__webpack_require__.initializeSharingData)===null||__webpack_require___initializeSharingData===void 0?void 0:__webpack_require___initializeSharingData.scopeToSharingDataMapping)!==null&&_ref2!==void 0?_ref2:{};const consumesLoadingChunkMapping=(_ref3=(__webpack_require___consumesLoadingData=__webpack_require__.consumesLoadingData)===null||__webpack_require___consumesLoadingData===void 0?void 0:__webpack_require___consumesLoadingData.chunkMapping)!==null&&_ref3!==void 0?_ref3:{};const consumesLoadingModuleToConsumeDataMapping=(_ref4=(__webpack_require___consumesLoadingData1=__webpack_require__.consumesLoadingData)===null||__webpack_require___consumesLoadingData1===void 0?void 0:__webpack_require___consumesLoadingData1.moduleIdToConsumeDataMapping)!==null&&_ref4!==void 0?_ref4:{};const consumesLoadinginstalledModules={};const initializeSharingInitPromises=[];const initializeSharingInitTokens={};const containerShareScope=(__webpack_require___initializeExposesData=__webpack_require__.initializeExposesData)===null||__webpack_require___initializeExposesData===void 0?void 0:__webpack_require___initializeExposesData.shareScope;for(const key in (_home_runner_work_jupyterlab_jupyterlab_jupyter_releaser_checkout_node_modules_module_federation_webpack_bundler_runtime_dist_index_cjs__rspack_import_0_default())){__webpack_require__.federation[key]=(_home_runner_work_jupyterlab_jupyterlab_jupyter_releaser_checkout_node_modules_module_federation_webpack_bundler_runtime_dist_index_cjs__rspack_import_0_default())[key]}early(__webpack_require__.federation,"libraryType",()=>__module_federation_library_type__);early(__webpack_require__.federation,"sharedFallback",()=>__module_federation_share_fallbacks__);const sharedFallback=__webpack_require__.federation.sharedFallback;early(__webpack_require__.federation,"consumesLoadingModuleToHandlerMapping",()=>{const consumesLoadingModuleToHandlerMapping={};for(let[moduleId,data]of Object.entries(consumesLoadingModuleToConsumeDataMapping)){var __webpack_require___federation_bundlerRuntime;consumesLoadingModuleToHandlerMapping[moduleId]={getter:sharedFallback?(__webpack_require___federation_bundlerRuntime=__webpack_require__.federation.bundlerRuntime)===null||__webpack_require___federation_bundlerRuntime===void 0?void 0:__webpack_require___federation_bundlerRuntime.getSharedFallbackGetter({shareKey:data.shareKey,factory:data.fallback,webpackRequire:__webpack_require__,libraryType:__webpack_require__.federation.libraryType}):data.fallback,treeShakingGetter:sharedFallback?data.fallback:undefined,shareInfo:{shareConfig:{fixedDependencies:false,requiredVersion:data.requiredVersion,strictVersion:data.strictVersion,singleton:data.singleton,eager:data.eager},scope:[data.shareScope]},shareKey:data.shareKey,treeShaking:__webpack_require__.federation.sharedFallback?{get:data.fallback,mode:data.treeShakingMode}:undefined}}return consumesLoadingModuleToHandlerMapping});early(__webpack_require__.federation,"initOptions",()=>({}));early(__webpack_require__.federation.initOptions,"name",()=>__module_federation_container_name__);early(__webpack_require__.federation.initOptions,"shareStrategy",()=>__module_federation_share_strategy__);early(__webpack_require__.federation.initOptions,"shared",()=>{const shared={};for(let[scope,stages]of Object.entries(initializeSharingScopeToInitDataMapping)){for(let stage of stages){if(typeof stage==="object"&&stage!==null){const{name,version,factory,eager,singleton,requiredVersion,strictVersion,treeShakingMode}=stage;const shareConfig={};const isValidValue=function(val){return typeof val!=="undefined"};if(isValidValue(singleton)){shareConfig.singleton=singleton}if(isValidValue(requiredVersion)){shareConfig.requiredVersion=requiredVersion}if(isValidValue(eager)){shareConfig.eager=eager}if(isValidValue(strictVersion)){shareConfig.strictVersion=strictVersion}const options={version,scope:[scope],shareConfig,get:factory,treeShaking:treeShakingMode?{mode:treeShakingMode}:undefined};if(shared[name]){shared[name].push(options)}else{shared[name]=[options]}}}}return shared});merge(__webpack_require__.federation.initOptions,"remotes",()=>Object.values(__module_federation_remote_infos__).flat().filter(remote=>remote.externalType==="script"));merge(__webpack_require__.federation.initOptions,"plugins",()=>__module_federation_runtime_plugins__);early(__webpack_require__.federation,"bundlerRuntimeOptions",()=>({}));early(__webpack_require__.federation.bundlerRuntimeOptions,"remotes",()=>({}));early(__webpack_require__.federation.bundlerRuntimeOptions.remotes,"chunkMapping",()=>remotesLoadingChunkMapping);early(__webpack_require__.federation.bundlerRuntimeOptions.remotes,"remoteInfos",()=>__module_federation_remote_infos__);early(__webpack_require__.federation.bundlerRuntimeOptions.remotes,"idToExternalAndNameMapping",()=>{const remotesLoadingIdToExternalAndNameMappingMapping={};for(let[moduleId,data]of Object.entries(remotesLoadingModuleIdToRemoteDataMapping)){remotesLoadingIdToExternalAndNameMappingMapping[moduleId]=[data.shareScope,data.name,data.externalModuleId,data.remoteName]}return remotesLoadingIdToExternalAndNameMappingMapping});early(__webpack_require__.federation.bundlerRuntimeOptions.remotes,"webpackRequire",()=>__webpack_require__);merge(__webpack_require__.federation.bundlerRuntimeOptions.remotes,"idToRemoteMap",()=>{const idToRemoteMap={};for(let[id,remoteData]of Object.entries(remotesLoadingModuleIdToRemoteDataMapping)){const info=__module_federation_remote_infos__[remoteData.remoteName];if(info)idToRemoteMap[id]=info}return idToRemoteMap});override(__webpack_require__,"S",__webpack_require__.federation.bundlerRuntime.S);if(__webpack_require__.federation.attachShareScopeMap){__webpack_require__.federation.attachShareScopeMap(__webpack_require__)}override(__webpack_require__.f,"remotes",(chunkId,promises)=>__webpack_require__.federation.bundlerRuntime.remotes({chunkId,promises,chunkMapping:remotesLoadingChunkMapping,idToExternalAndNameMapping:__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToExternalAndNameMapping,idToRemoteMap:__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToRemoteMap,webpackRequire:__webpack_require__}));override(__webpack_require__.f,"consumes",(chunkId,promises)=>__webpack_require__.federation.bundlerRuntime.consumes({chunkId,promises,chunkMapping:consumesLoadingChunkMapping,moduleToHandlerMapping:__webpack_require__.federation.consumesLoadingModuleToHandlerMapping,installedModules:consumesLoadinginstalledModules,webpackRequire:__webpack_require__}));override(__webpack_require__,"I",(name,initScope)=>__webpack_require__.federation.bundlerRuntime.I({shareScopeName:name,initScope,initPromises:initializeSharingInitPromises,initTokens:initializeSharingInitTokens,webpackRequire:__webpack_require__}));override(__webpack_require__,"initContainer",(shareScope,initScope,remoteEntryInitOptions)=>__webpack_require__.federation.bundlerRuntime.initContainerEntry({shareScope,initScope,remoteEntryInitOptions,shareScopeKey:containerShareScope,webpackRequire:__webpack_require__}));override(__webpack_require__,"getContainer",(module,getScope)=>{var moduleMap=__webpack_require__.initializeExposesData.moduleMap;__webpack_require__.R=getScope;getScope=Object.prototype.hasOwnProperty.call(moduleMap,module)?moduleMap[module]():Promise.resolve().then(()=>{throw new Error('Module "'+module+'" does not exist in container.')});__webpack_require__.R=undefined;return getScope});__webpack_require__.federation.instance=__webpack_require__.federation.bundlerRuntime.init({webpackRequire:__webpack_require__});if((__webpack_require___consumesLoadingData2=__webpack_require__.consumesLoadingData)===null||__webpack_require___consumesLoadingData2===void 0?void 0:__webpack_require___consumesLoadingData2.initialConsumes){__webpack_require__.federation.bundlerRuntime.installInitialConsumes({webpackRequire:__webpack_require__,installedModules:consumesLoadinginstalledModules,initialConsumes:__webpack_require__.consumesLoadingData.initialConsumes,moduleToHandlerMapping:__webpack_require__.federation.consumesLoadingModuleToHandlerMapping})}}

},
"webpack/container/entry/@jupyterlab/galata-extension"(__unused_rspack_module, exports, __webpack_require__) {
"use strict";

__webpack_require__.d(exports, {
	get: () => (__webpack_require__.getContainer),
	init: () => (__webpack_require__.initContainer)
});

},
"../../node_modules/@module-federation/error-codes/dist/browser.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
const require_getShortErrorMsg = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/getShortErrorMsg.cjs");

//#region src/browser.ts
function logAndReport(code, descMap, args, logger, originalErrorMsg, context) {
	return logger(require_getShortErrorMsg.getShortErrorMsg(code, descMap, args, originalErrorMsg));
}

//#endregion
exports.logAndReport = logAndReport;
//# sourceMappingURL=browser.cjs.map

},
"../../node_modules/@module-federation/error-codes/dist/desc.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/error-codes.cjs");

//#region src/desc.ts
const runtimeDescMap = {
	[require_error_codes.RUNTIME_001]: "Failed to get remoteEntry exports.",
	[require_error_codes.RUNTIME_002]: "The remote entry interface does not contain \"init\"",
	[require_error_codes.RUNTIME_003]: "Failed to get manifest.",
	[require_error_codes.RUNTIME_004]: "Failed to locate remote.",
	[require_error_codes.RUNTIME_005]: "Invalid loadShareSync function call from bundler runtime",
	[require_error_codes.RUNTIME_006]: "Invalid loadShareSync function call from runtime",
	[require_error_codes.RUNTIME_007]: "Failed to get remote snapshot.",
	[require_error_codes.RUNTIME_008]: "Failed to load script resources.",
	[require_error_codes.RUNTIME_009]: "Please call createInstance first.",
	[require_error_codes.RUNTIME_010]: "The name option cannot be changed after initialization. If you want to create a new instance with a different name, please use \"createInstance\" api.",
	[require_error_codes.RUNTIME_011]: "The remoteEntry URL is missing from the remote snapshot.",
	[require_error_codes.RUNTIME_012]: "The getter for the shared module is not a function. This may be caused by setting \"shared.import: false\" without the host providing the corresponding lib."
};
const typeDescMap = { [require_error_codes.TYPE_001]: "Failed to generate type declaration. Execute the below cmd to reproduce and fix the error." };
const buildDescMap = {
	[require_error_codes.BUILD_001]: "Failed to find expose module.",
	[require_error_codes.BUILD_002]: "PublicPath is required in prod mode."
};
const errorDescMap = {
	...runtimeDescMap,
	...typeDescMap,
	...buildDescMap
};

//#endregion
exports.buildDescMap = buildDescMap;
exports.errorDescMap = errorDescMap;
exports.runtimeDescMap = runtimeDescMap;
exports.typeDescMap = typeDescMap;
//# sourceMappingURL=desc.cjs.map

},
"../../node_modules/@module-federation/error-codes/dist/error-codes.cjs"(__unused_rspack_module, exports) {

//#region src/error-codes.ts
const RUNTIME_001 = "RUNTIME-001";
const RUNTIME_002 = "RUNTIME-002";
const RUNTIME_003 = "RUNTIME-003";
const RUNTIME_004 = "RUNTIME-004";
const RUNTIME_005 = "RUNTIME-005";
const RUNTIME_006 = "RUNTIME-006";
const RUNTIME_007 = "RUNTIME-007";
const RUNTIME_008 = "RUNTIME-008";
const RUNTIME_009 = "RUNTIME-009";
const RUNTIME_010 = "RUNTIME-010";
const RUNTIME_011 = "RUNTIME-011";
const RUNTIME_012 = "RUNTIME-012";
const TYPE_001 = "TYPE-001";
const BUILD_001 = "BUILD-001";
const BUILD_002 = "BUILD-002";

//#endregion
exports.BUILD_001 = BUILD_001;
exports.BUILD_002 = BUILD_002;
exports.RUNTIME_001 = RUNTIME_001;
exports.RUNTIME_002 = RUNTIME_002;
exports.RUNTIME_003 = RUNTIME_003;
exports.RUNTIME_004 = RUNTIME_004;
exports.RUNTIME_005 = RUNTIME_005;
exports.RUNTIME_006 = RUNTIME_006;
exports.RUNTIME_007 = RUNTIME_007;
exports.RUNTIME_008 = RUNTIME_008;
exports.RUNTIME_009 = RUNTIME_009;
exports.RUNTIME_010 = RUNTIME_010;
exports.RUNTIME_011 = RUNTIME_011;
exports.RUNTIME_012 = RUNTIME_012;
exports.TYPE_001 = TYPE_001;
//# sourceMappingURL=error-codes.cjs.map

},
"../../node_modules/@module-federation/error-codes/dist/getShortErrorMsg.cjs"(__unused_rspack_module, exports) {

//#region src/getShortErrorMsg.ts
const getDocsUrl = (errorCode) => {
	return `View the docs to see how to solve: https://module-federation.io/guide/troubleshooting/${errorCode.split("-")[0].toLowerCase()}#${errorCode.toLowerCase()}`;
};
const getShortErrorMsg = (errorCode, errorDescMap, args, originalErrorMsg) => {
	const msg = [`${[errorDescMap[errorCode]]} #${errorCode}`];
	args && msg.push(`args: ${JSON.stringify(args)}`);
	msg.push(getDocsUrl(errorCode));
	originalErrorMsg && msg.push(`Original Error Message:\n ${originalErrorMsg}`);
	return msg.join("\n");
};

//#endregion
exports.getShortErrorMsg = getShortErrorMsg;
//# sourceMappingURL=getShortErrorMsg.cjs.map

},
"../../node_modules/@module-federation/error-codes/dist/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
const require_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/error-codes.cjs");
const require_getShortErrorMsg = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/getShortErrorMsg.cjs");
const require_desc = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/desc.cjs");

exports.BUILD_001 = require_error_codes.BUILD_001;
exports.BUILD_002 = require_error_codes.BUILD_002;
exports.RUNTIME_001 = require_error_codes.RUNTIME_001;
exports.RUNTIME_002 = require_error_codes.RUNTIME_002;
exports.RUNTIME_003 = require_error_codes.RUNTIME_003;
exports.RUNTIME_004 = require_error_codes.RUNTIME_004;
exports.RUNTIME_005 = require_error_codes.RUNTIME_005;
exports.RUNTIME_006 = require_error_codes.RUNTIME_006;
exports.RUNTIME_007 = require_error_codes.RUNTIME_007;
exports.RUNTIME_008 = require_error_codes.RUNTIME_008;
exports.RUNTIME_009 = require_error_codes.RUNTIME_009;
exports.RUNTIME_010 = require_error_codes.RUNTIME_010;
exports.RUNTIME_011 = require_error_codes.RUNTIME_011;
exports.RUNTIME_012 = require_error_codes.RUNTIME_012;
exports.TYPE_001 = require_error_codes.TYPE_001;
exports.buildDescMap = require_desc.buildDescMap;
exports.errorDescMap = require_desc.errorDescMap;
exports.getShortErrorMsg = require_getShortErrorMsg.getShortErrorMsg;
exports.runtimeDescMap = require_desc.runtimeDescMap;
exports.typeDescMap = require_desc.typeDescMap;

},
"../../node_modules/@module-federation/runtime-core/dist/_virtual/_rolldown/runtime.cjs"(__unused_rspack_module, exports) {
//#region \0rolldown/runtime.js
var __defProp = Object.defineProperty;
var __exportAll = (all, no_symbols) => {
	let target = {};
	for (var name in all) {
		__defProp(target, name, {
			get: all[name],
			enumerable: true
		});
	}
	if (!no_symbols) {
		__defProp(target, Symbol.toStringTag, { value: "Module" });
	}
	return target;
};

//#endregion

exports.__exportAll = __exportAll;

},
"../../node_modules/@module-federation/runtime-core/dist/constant.cjs"(__unused_rspack_module, exports) {

//#region src/constant.ts
const DEFAULT_SCOPE = "default";
const DEFAULT_REMOTE_TYPE = "global";

//#endregion
exports.DEFAULT_REMOTE_TYPE = DEFAULT_REMOTE_TYPE;
exports.DEFAULT_SCOPE = DEFAULT_SCOPE;
//# sourceMappingURL=constant.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/core.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/constant.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
const require_env = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/env.cjs");
const require_plugin = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/plugin.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_index$1 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/module/index.cjs");
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");
const require_asyncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs");
const require_syncWaterfallHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs");
const require_asyncWaterfallHooks = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs");
const require_pluginSystem = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/index.cjs");
const require_index$3 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/index.cjs");
const require_generate_preload_assets = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/plugins/generate-preload-assets.cjs");
const require_SnapshotHandler = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/SnapshotHandler.cjs");
const require_index$4 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/shared/index.cjs");
const require_index$5 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/remote/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/core.ts
const USE_SNAPSHOT = typeof FEDERATION_OPTIMIZE_NO_SNAPSHOT_PLUGIN === "boolean" ? !FEDERATION_OPTIMIZE_NO_SNAPSHOT_PLUGIN : true;
var ModuleFederation = class {
	constructor(userOptions) {
		this.hooks = new require_pluginSystem.PluginSystem({
			beforeInit: new require_syncWaterfallHook.SyncWaterfallHook("beforeInit"),
			init: new require_syncHook.SyncHook(),
			beforeInitContainer: new require_asyncWaterfallHooks.AsyncWaterfallHook("beforeInitContainer"),
			initContainer: new require_asyncWaterfallHooks.AsyncWaterfallHook("initContainer")
		});
		this.version = "2.4.0";
		this.moduleCache = /* @__PURE__ */ new Map();
		this.loaderHook = new require_pluginSystem.PluginSystem({
			getModuleInfo: new require_syncHook.SyncHook(),
			createScript: new require_syncHook.SyncHook(),
			createLink: new require_syncHook.SyncHook(),
			fetch: new require_asyncHook.AsyncHook(),
			loadEntryError: new require_asyncHook.AsyncHook(),
			getModuleFactory: new require_asyncHook.AsyncHook()
		});
		this.bridgeHook = new require_pluginSystem.PluginSystem({
			beforeBridgeRender: new require_syncHook.SyncHook(),
			afterBridgeRender: new require_syncHook.SyncHook(),
			beforeBridgeDestroy: new require_syncHook.SyncHook(),
			afterBridgeDestroy: new require_syncHook.SyncHook()
		});
		const plugins = USE_SNAPSHOT ? [require_index$3.snapshotPlugin(), require_generate_preload_assets.generatePreloadAssetsPlugin()] : [];
		const defaultOptions = {
			id: require_env.getBuilderId(),
			name: userOptions.name,
			plugins,
			remotes: [],
			shared: {},
			inBrowser: _module_federation_sdk.isBrowserEnvValue
		};
		this.name = userOptions.name;
		this.options = defaultOptions;
		this.snapshotHandler = new require_SnapshotHandler.SnapshotHandler(this);
		this.sharedHandler = new require_index$4.SharedHandler(this);
		this.remoteHandler = new require_index$5.RemoteHandler(this);
		this.shareScopeMap = this.sharedHandler.shareScopeMap;
		this.registerPlugins([...defaultOptions.plugins, ...userOptions.plugins || []]);
		this.options = this.formatOptions(defaultOptions, userOptions);
	}
	initOptions(userOptions) {
		if (userOptions.name && userOptions.name !== this.options.name) require_logger.error((0, _module_federation_error_codes.getShortErrorMsg)(_module_federation_error_codes.RUNTIME_010, _module_federation_error_codes.runtimeDescMap));
		this.registerPlugins(userOptions.plugins);
		const options = this.formatOptions(this.options, userOptions);
		this.options = options;
		return options;
	}
	async loadShare(pkgName, extraOptions) {
		return this.sharedHandler.loadShare(pkgName, extraOptions);
	}
	loadShareSync(pkgName, extraOptions) {
		return this.sharedHandler.loadShareSync(pkgName, extraOptions);
	}
	initializeSharing(shareScopeName = require_constant.DEFAULT_SCOPE, extraOptions) {
		return this.sharedHandler.initializeSharing(shareScopeName, extraOptions);
	}
	initRawContainer(name, url, container) {
		const remoteInfo = require_load.getRemoteInfo({
			name,
			entry: url
		});
		const module = new require_index$1.Module({
			host: this,
			remoteInfo
		});
		module.remoteEntryExports = container;
		this.moduleCache.set(name, module);
		return module;
	}
	async loadRemote(id, options) {
		return this.remoteHandler.loadRemote(id, options);
	}
	async preloadRemote(preloadOptions) {
		return this.remoteHandler.preloadRemote(preloadOptions);
	}
	initShareScopeMap(scopeName, shareScope, extraOptions = {}) {
		this.sharedHandler.initShareScopeMap(scopeName, shareScope, extraOptions);
	}
	formatOptions(globalOptions, userOptions) {
		const { allShareInfos: shared } = require_share.formatShareConfigs(globalOptions, userOptions);
		const { userOptions: userOptionsRes, options: globalOptionsRes } = this.hooks.lifecycle.beforeInit.emit({
			origin: this,
			userOptions,
			options: globalOptions,
			shareInfo: shared
		});
		const remotes = this.remoteHandler.formatAndRegisterRemote(globalOptionsRes, userOptionsRes);
		const { allShareInfos } = this.sharedHandler.registerShared(globalOptionsRes, userOptionsRes);
		const plugins = [...globalOptionsRes.plugins];
		if (userOptionsRes.plugins) userOptionsRes.plugins.forEach((plugin) => {
			if (!plugins.includes(plugin)) plugins.push(plugin);
		});
		const optionsRes = {
			...globalOptions,
			...userOptions,
			plugins,
			remotes,
			shared: allShareInfos,
			id: userOptionsRes.id || globalOptions.id
		};
		this.hooks.lifecycle.init.emit({
			origin: this,
			options: optionsRes
		});
		return optionsRes;
	}
	registerPlugins(plugins) {
		const pluginRes = require_plugin.registerPlugins(plugins, this);
		this.options.plugins = this.options.plugins.reduce((res, plugin) => {
			if (!plugin) return res;
			if (res && !res.find((item) => item.name === plugin.name)) res.push(plugin);
			return res;
		}, pluginRes || []);
	}
	registerRemotes(remotes, options) {
		return this.remoteHandler.registerRemotes(remotes, options);
	}
	registerShared(shared) {
		this.sharedHandler.registerShared(this.options, {
			...this.options,
			shared
		});
	}
};

//#endregion
exports.ModuleFederation = ModuleFederation;
//# sourceMappingURL=core.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/global.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/global.ts
const CurrentGlobal = typeof globalThis === "object" ? globalThis : window;
const nativeGlobal = (() => {
	try {
		return document.defaultView;
	} catch {
		return CurrentGlobal;
	}
})();
const Global = nativeGlobal;
function definePropertyGlobalVal(target, key, val) {
	Object.defineProperty(target, key, {
		value: val,
		configurable: false,
		writable: true
	});
}
function includeOwnProperty(target, key) {
	return Object.hasOwnProperty.call(target, key);
}
if (!includeOwnProperty(CurrentGlobal, "__GLOBAL_LOADING_REMOTE_ENTRY__")) definePropertyGlobalVal(CurrentGlobal, "__GLOBAL_LOADING_REMOTE_ENTRY__", {});
const globalLoading = CurrentGlobal.__GLOBAL_LOADING_REMOTE_ENTRY__;
function setGlobalDefaultVal(target) {
	if (includeOwnProperty(target, "__VMOK__") && !includeOwnProperty(target, "__FEDERATION__")) definePropertyGlobalVal(target, "__FEDERATION__", target.__VMOK__);
	if (!includeOwnProperty(target, "__FEDERATION__")) {
		definePropertyGlobalVal(target, "__FEDERATION__", {
			__GLOBAL_PLUGIN__: [],
			__INSTANCES__: [],
			moduleInfo: {},
			__SHARE__: {},
			__MANIFEST_LOADING__: {},
			__PRELOADED_MAP__: /* @__PURE__ */ new Map()
		});
		definePropertyGlobalVal(target, "__VMOK__", target.__FEDERATION__);
	}
	target.__FEDERATION__.__GLOBAL_PLUGIN__ ??= [];
	target.__FEDERATION__.__INSTANCES__ ??= [];
	target.__FEDERATION__.moduleInfo ??= {};
	target.__FEDERATION__.__SHARE__ ??= {};
	target.__FEDERATION__.__MANIFEST_LOADING__ ??= {};
	target.__FEDERATION__.__PRELOADED_MAP__ ??= /* @__PURE__ */ new Map();
}
setGlobalDefaultVal(CurrentGlobal);
setGlobalDefaultVal(nativeGlobal);
function resetFederationGlobalInfo() {
	CurrentGlobal.__FEDERATION__.__GLOBAL_PLUGIN__ = [];
	CurrentGlobal.__FEDERATION__.__INSTANCES__ = [];
	CurrentGlobal.__FEDERATION__.moduleInfo = {};
	CurrentGlobal.__FEDERATION__.__SHARE__ = {};
	CurrentGlobal.__FEDERATION__.__MANIFEST_LOADING__ = {};
	Object.keys(globalLoading).forEach((key) => {
		delete globalLoading[key];
	});
}
function setGlobalFederationInstance(FederationInstance) {
	CurrentGlobal.__FEDERATION__.__INSTANCES__.push(FederationInstance);
}
function getGlobalFederationConstructor() {
	return CurrentGlobal.__FEDERATION__.__DEBUG_CONSTRUCTOR__;
}
function setGlobalFederationConstructor(FederationConstructor, isDebug = (0, _module_federation_sdk.isDebugMode)()) {
	if (isDebug) {
		CurrentGlobal.__FEDERATION__.__DEBUG_CONSTRUCTOR__ = FederationConstructor;
		CurrentGlobal.__FEDERATION__.__DEBUG_CONSTRUCTOR_VERSION__ = "2.4.0";
	}
}
function getInfoWithoutType(target, key) {
	if (typeof key === "string") if (target[key]) return {
		value: target[key],
		key
	};
	else {
		const targetKeys = Object.keys(target);
		for (const targetKey of targetKeys) {
			const [targetTypeOrName, _] = targetKey.split(":");
			const nKey = `${targetTypeOrName}:${key}`;
			const typeWithKeyRes = target[nKey];
			if (typeWithKeyRes) return {
				value: typeWithKeyRes,
				key: nKey
			};
		}
		return {
			value: void 0,
			key
		};
	}
	else require_logger.error(`getInfoWithoutType: "key" must be a string, got ${typeof key} (${JSON.stringify(key)}).`);
}
const getGlobalSnapshot = () => nativeGlobal.__FEDERATION__.moduleInfo;
const getTargetSnapshotInfoByModuleInfo = (moduleInfo, snapshot) => {
	const getModuleInfo = getInfoWithoutType(snapshot, require_tool.getFMId(moduleInfo)).value;
	if (getModuleInfo && !getModuleInfo.version && "version" in moduleInfo && moduleInfo["version"]) getModuleInfo.version = moduleInfo["version"];
	if (getModuleInfo) return getModuleInfo;
	if ("version" in moduleInfo && moduleInfo["version"]) {
		const { version, ...resModuleInfo } = moduleInfo;
		const moduleKeyWithoutVersion = require_tool.getFMId(resModuleInfo);
		const getModuleInfoWithoutVersion = getInfoWithoutType(nativeGlobal.__FEDERATION__.moduleInfo, moduleKeyWithoutVersion).value;
		if (getModuleInfoWithoutVersion?.version === version) return getModuleInfoWithoutVersion;
	}
};
const getGlobalSnapshotInfoByModuleInfo = (moduleInfo) => getTargetSnapshotInfoByModuleInfo(moduleInfo, nativeGlobal.__FEDERATION__.moduleInfo);
const setGlobalSnapshotInfoByModuleInfo = (remoteInfo, moduleDetailInfo) => {
	const moduleKey = require_tool.getFMId(remoteInfo);
	nativeGlobal.__FEDERATION__.moduleInfo[moduleKey] = moduleDetailInfo;
	return nativeGlobal.__FEDERATION__.moduleInfo;
};
const addGlobalSnapshot = (moduleInfos) => {
	nativeGlobal.__FEDERATION__.moduleInfo = {
		...nativeGlobal.__FEDERATION__.moduleInfo,
		...moduleInfos
	};
	return () => {
		const keys = Object.keys(moduleInfos);
		for (const key of keys) delete nativeGlobal.__FEDERATION__.moduleInfo[key];
	};
};
const getRemoteEntryExports = (name, globalName) => {
	const remoteEntryKey = globalName || `__FEDERATION_${name}:custom__`;
	return {
		remoteEntryKey,
		entryExports: CurrentGlobal[remoteEntryKey]
	};
};
const registerGlobalPlugins = (plugins) => {
	const { __GLOBAL_PLUGIN__ } = nativeGlobal.__FEDERATION__;
	plugins.forEach((plugin) => {
		if (__GLOBAL_PLUGIN__.findIndex((p) => p.name === plugin.name) === -1) __GLOBAL_PLUGIN__.push(plugin);
		else require_logger.warn(`The plugin ${plugin.name} has been registered.`);
	});
};
const getGlobalHostPlugins = () => nativeGlobal.__FEDERATION__.__GLOBAL_PLUGIN__;
const getPreloaded = (id) => CurrentGlobal.__FEDERATION__.__PRELOADED_MAP__.get(id);
const setPreloaded = (id) => CurrentGlobal.__FEDERATION__.__PRELOADED_MAP__.set(id, true);

//#endregion
exports.CurrentGlobal = CurrentGlobal;
exports.Global = Global;
exports.addGlobalSnapshot = addGlobalSnapshot;
exports.getGlobalFederationConstructor = getGlobalFederationConstructor;
exports.getGlobalHostPlugins = getGlobalHostPlugins;
exports.getGlobalSnapshot = getGlobalSnapshot;
exports.getGlobalSnapshotInfoByModuleInfo = getGlobalSnapshotInfoByModuleInfo;
exports.getInfoWithoutType = getInfoWithoutType;
exports.getPreloaded = getPreloaded;
exports.getRemoteEntryExports = getRemoteEntryExports;
exports.getTargetSnapshotInfoByModuleInfo = getTargetSnapshotInfoByModuleInfo;
exports.globalLoading = globalLoading;
exports.nativeGlobal = nativeGlobal;
exports.registerGlobalPlugins = registerGlobalPlugins;
exports.resetFederationGlobalInfo = resetFederationGlobalInfo;
exports.setGlobalFederationConstructor = setGlobalFederationConstructor;
exports.setGlobalFederationInstance = setGlobalFederationInstance;
exports.setGlobalSnapshotInfoByModuleInfo = setGlobalSnapshotInfoByModuleInfo;
exports.setPreloaded = setPreloaded;
//# sourceMappingURL=global.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/helpers.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
const require_manifest = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_preload = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/preload.cjs");

//#region src/helpers.ts
const ShareUtils = {
	getRegisteredShare: require_share.getRegisteredShare,
	getGlobalShareScope: require_share.getGlobalShareScope
};
const GlobalUtils = {
	Global: require_global.Global,
	nativeGlobal: require_global.nativeGlobal,
	resetFederationGlobalInfo: require_global.resetFederationGlobalInfo,
	setGlobalFederationInstance: require_global.setGlobalFederationInstance,
	getGlobalFederationConstructor: require_global.getGlobalFederationConstructor,
	setGlobalFederationConstructor: require_global.setGlobalFederationConstructor,
	getInfoWithoutType: require_global.getInfoWithoutType,
	getGlobalSnapshot: require_global.getGlobalSnapshot,
	getTargetSnapshotInfoByModuleInfo: require_global.getTargetSnapshotInfoByModuleInfo,
	getGlobalSnapshotInfoByModuleInfo: require_global.getGlobalSnapshotInfoByModuleInfo,
	setGlobalSnapshotInfoByModuleInfo: require_global.setGlobalSnapshotInfoByModuleInfo,
	addGlobalSnapshot: require_global.addGlobalSnapshot,
	getRemoteEntryExports: require_global.getRemoteEntryExports,
	registerGlobalPlugins: require_global.registerGlobalPlugins,
	getGlobalHostPlugins: require_global.getGlobalHostPlugins,
	getPreloaded: require_global.getPreloaded,
	setPreloaded: require_global.setPreloaded
};
var helpers_default = {
	global: GlobalUtils,
	share: ShareUtils,
	utils: {
		matchRemoteWithNameAndExpose: require_manifest.matchRemoteWithNameAndExpose,
		preloadAssets: require_preload.preloadAssets,
		getRemoteInfo: require_load.getRemoteInfo
	}
};

//#endregion
exports["default"] = helpers_default;
//# sourceMappingURL=helpers.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_index = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/index.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
const require_manifest = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_helpers = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/helpers.cjs");
const require_index$2 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/module/index.cjs");
const require_core = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/core.cjs");
const require_index$3 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/type/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/index.ts
const helpers = require_helpers.default;

//#endregion
exports.CurrentGlobal = require_global.CurrentGlobal;
exports.Global = require_global.Global;
exports.Module = require_index$2.Module;
exports.ModuleFederation = require_core.ModuleFederation;
exports.addGlobalSnapshot = require_global.addGlobalSnapshot;
exports.assert = require_logger.assert;
exports.error = require_logger.error;
exports.getGlobalFederationConstructor = require_global.getGlobalFederationConstructor;
exports.getGlobalSnapshot = require_global.getGlobalSnapshot;
exports.getInfoWithoutType = require_global.getInfoWithoutType;
exports.getRegisteredShare = require_share.getRegisteredShare;
exports.getRemoteEntry = require_load.getRemoteEntry;
exports.getRemoteInfo = require_load.getRemoteInfo;
exports.helpers = helpers;
exports.isStaticResourcesEqual = require_tool.isStaticResourcesEqual;
Object.defineProperty(exports, "loadScript", ({
  enumerable: true,
  get: function () {
    return _module_federation_sdk.loadScript;
  }
}));
Object.defineProperty(exports, "loadScriptNode", ({
  enumerable: true,
  get: function () {
    return _module_federation_sdk.loadScriptNode;
  }
}));
exports.matchRemoteWithNameAndExpose = require_manifest.matchRemoteWithNameAndExpose;
exports.registerGlobalPlugins = require_global.registerGlobalPlugins;
exports.resetFederationGlobalInfo = require_global.resetFederationGlobalInfo;
exports.safeWrapper = require_tool.safeWrapper;
exports.satisfy = require_index.satisfy;
exports.setGlobalFederationConstructor = require_global.setGlobalFederationConstructor;
exports.setGlobalFederationInstance = require_global.setGlobalFederationInstance;
Object.defineProperty(exports, "types", ({
  enumerable: true,
  get: function () {
    return require_index$3.type_exports;
  }
}));
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/module/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
const require_context = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/module/index.ts
function createRemoteEntryInitOptions(remoteInfo, hostShareScopeMap, rawInitScope) {
	const localShareScopeMap = hostShareScopeMap;
	const shareScopeKeys = Array.isArray(remoteInfo.shareScope) ? remoteInfo.shareScope : [remoteInfo.shareScope];
	if (!shareScopeKeys.length) shareScopeKeys.push("default");
	shareScopeKeys.forEach((shareScopeKey) => {
		if (!localShareScopeMap[shareScopeKey]) localShareScopeMap[shareScopeKey] = {};
	});
	const remoteEntryInitOptions = {
		version: remoteInfo.version || "",
		shareScopeKeys: Array.isArray(remoteInfo.shareScope) ? shareScopeKeys : remoteInfo.shareScope || "default"
	};
	Object.defineProperty(remoteEntryInitOptions, "shareScopeMap", {
		value: localShareScopeMap,
		enumerable: false
	});
	return {
		remoteEntryInitOptions,
		shareScope: localShareScopeMap[shareScopeKeys[0]],
		initScope: rawInitScope ?? []
	};
}
var Module = class {
	constructor({ remoteInfo, host }) {
		this.inited = false;
		this.initing = false;
		this.lib = void 0;
		this.remoteInfo = remoteInfo;
		this.host = host;
	}
	async getEntry() {
		if (this.remoteEntryExports) return this.remoteEntryExports;
		const remoteEntryExports = await require_load.getRemoteEntry({
			origin: this.host,
			remoteInfo: this.remoteInfo,
			remoteEntryExports: this.remoteEntryExports
		});
		require_logger.assert(remoteEntryExports, `remoteEntryExports is undefined \n ${(0, _module_federation_sdk.safeToString)(this.remoteInfo)}`);
		this.remoteEntryExports = remoteEntryExports;
		return this.remoteEntryExports;
	}
	async init(id, remoteSnapshot, rawInitScope) {
		const remoteEntryExports = await this.getEntry();
		if (this.inited) return remoteEntryExports;
		if (this.initPromise) {
			await this.initPromise;
			return remoteEntryExports;
		}
		this.initing = true;
		this.initPromise = (async () => {
			const { remoteEntryInitOptions, shareScope, initScope } = createRemoteEntryInitOptions(this.remoteInfo, this.host.shareScopeMap, rawInitScope);
			const initContainerOptions = await this.host.hooks.lifecycle.beforeInitContainer.emit({
				shareScope,
				remoteEntryInitOptions,
				initScope,
				remoteInfo: this.remoteInfo,
				origin: this.host
			});
			if (typeof remoteEntryExports?.init === "undefined") require_logger.error(_module_federation_error_codes.RUNTIME_002, _module_federation_error_codes.runtimeDescMap, {
				hostName: this.host.name,
				remoteName: this.remoteInfo.name,
				remoteEntryUrl: this.remoteInfo.entry,
				remoteEntryKey: this.remoteInfo.entryGlobalName
			}, void 0, require_context.optionsToMFContext(this.host.options));
			await remoteEntryExports.init(initContainerOptions.shareScope, initContainerOptions.initScope, initContainerOptions.remoteEntryInitOptions);
			await this.host.hooks.lifecycle.initContainer.emit({
				...initContainerOptions,
				id,
				remoteSnapshot,
				remoteEntryExports
			});
			this.inited = true;
		})();
		try {
			await this.initPromise;
		} finally {
			this.initing = false;
			this.initPromise = void 0;
		}
		return remoteEntryExports;
	}
	async get(id, expose, options, remoteSnapshot) {
		const { loadFactory = true } = options || { loadFactory: true };
		const remoteEntryExports = await this.init(id, remoteSnapshot);
		this.lib = remoteEntryExports;
		let moduleFactory;
		moduleFactory = await this.host.loaderHook.lifecycle.getModuleFactory.emit({
			remoteEntryExports,
			expose,
			moduleInfo: this.remoteInfo
		});
		if (!moduleFactory) moduleFactory = await remoteEntryExports.get(expose);
		require_logger.assert(moduleFactory, `${require_tool.getFMId(this.remoteInfo)} remote don't export ${expose}.`);
		const symbolName = require_tool.processModuleAlias(this.remoteInfo.name, expose);
		const wrapModuleFactory = this.wraperFactory(moduleFactory, symbolName);
		if (!loadFactory) return wrapModuleFactory;
		return await wrapModuleFactory();
	}
	wraperFactory(moduleFactory, id) {
		function defineModuleId(res, id) {
			if (res && typeof res === "object" && Object.isExtensible(res) && !Object.getOwnPropertyDescriptor(res, Symbol.for("mf_module_id"))) Object.defineProperty(res, Symbol.for("mf_module_id"), {
				value: id,
				enumerable: false
			});
		}
		if (moduleFactory instanceof Promise) return async () => {
			const res = await moduleFactory();
			defineModuleId(res, id);
			return res;
		};
		else return () => {
			const res = moduleFactory();
			defineModuleId(res, id);
			return res;
		};
	}
};

//#endregion
exports.Module = Module;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/plugins/generate-preload-assets.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_preload = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/preload.cjs");
const require_index$1 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/plugins/generate-preload-assets.ts
function splitId(id) {
	const splitInfo = id.split(":");
	if (splitInfo.length === 1) return {
		name: splitInfo[0],
		version: void 0
	};
	else if (splitInfo.length === 2) return {
		name: splitInfo[0],
		version: splitInfo[1]
	};
	else return {
		name: splitInfo[1],
		version: splitInfo[2]
	};
}
function traverseModuleInfo(globalSnapshot, remoteInfo, traverse, isRoot, memo = {}, remoteSnapshot) {
	const { value: snapshotValue } = require_global.getInfoWithoutType(globalSnapshot, require_tool.getFMId(remoteInfo));
	const effectiveRemoteSnapshot = remoteSnapshot || snapshotValue;
	if (effectiveRemoteSnapshot && !(0, _module_federation_sdk.isManifestProvider)(effectiveRemoteSnapshot)) {
		traverse(effectiveRemoteSnapshot, remoteInfo, isRoot);
		if (effectiveRemoteSnapshot.remotesInfo) {
			const remoteKeys = Object.keys(effectiveRemoteSnapshot.remotesInfo);
			for (const key of remoteKeys) {
				if (memo[key]) continue;
				memo[key] = true;
				const subRemoteInfo = splitId(key);
				const remoteValue = effectiveRemoteSnapshot.remotesInfo[key];
				traverseModuleInfo(globalSnapshot, {
					name: subRemoteInfo.name,
					version: remoteValue.matchedVersion
				}, traverse, false, memo, void 0);
			}
		}
	}
}
const isExisted = (type, url) => {
	return document.querySelector(`${type}[${type === "link" ? "href" : "src"}="${url}"]`);
};
function generatePreloadAssets(origin, preloadOptions, remote, globalSnapshot, remoteSnapshot) {
	const cssAssets = [];
	const jsAssets = [];
	const entryAssets = [];
	const loadedSharedJsAssets = /* @__PURE__ */ new Set();
	const loadedSharedCssAssets = /* @__PURE__ */ new Set();
	const { options } = origin;
	const { preloadConfig: rootPreloadConfig } = preloadOptions;
	const { depsRemote } = rootPreloadConfig;
	traverseModuleInfo(globalSnapshot, remote, (moduleInfoSnapshot, remoteInfo, isRoot) => {
		let preloadConfig;
		if (isRoot) preloadConfig = rootPreloadConfig;
		else if (Array.isArray(depsRemote)) {
			const findPreloadConfig = depsRemote.find((remoteConfig) => {
				if (remoteConfig.nameOrAlias === remoteInfo.name || remoteConfig.nameOrAlias === remoteInfo.alias) return true;
				return false;
			});
			if (!findPreloadConfig) return;
			preloadConfig = require_preload.defaultPreloadArgs(findPreloadConfig);
		} else if (depsRemote === true) preloadConfig = rootPreloadConfig;
		else return;
		const remoteEntryUrl = (0, _module_federation_sdk.getResourceUrl)(moduleInfoSnapshot, require_tool.getRemoteEntryInfoFromSnapshot(moduleInfoSnapshot).url);
		if (remoteEntryUrl) entryAssets.push({
			name: remoteInfo.name,
			moduleInfo: {
				name: remoteInfo.name,
				entry: remoteEntryUrl,
				type: "remoteEntryType" in moduleInfoSnapshot ? moduleInfoSnapshot.remoteEntryType : "global",
				entryGlobalName: "globalName" in moduleInfoSnapshot ? moduleInfoSnapshot.globalName : remoteInfo.name,
				shareScope: "",
				version: "version" in moduleInfoSnapshot ? moduleInfoSnapshot.version : void 0
			},
			url: remoteEntryUrl
		});
		let moduleAssetsInfo = "modules" in moduleInfoSnapshot ? moduleInfoSnapshot.modules : [];
		const normalizedPreloadExposes = require_preload.normalizePreloadExposes(preloadConfig.exposes);
		if (normalizedPreloadExposes.length && "modules" in moduleInfoSnapshot) moduleAssetsInfo = moduleInfoSnapshot?.modules?.reduce((assets, moduleAssetInfo) => {
			if (normalizedPreloadExposes?.indexOf(moduleAssetInfo.moduleName) !== -1) assets.push(moduleAssetInfo);
			return assets;
		}, []);
		function handleAssets(assets) {
			const assetsRes = assets.map((asset) => (0, _module_federation_sdk.getResourceUrl)(moduleInfoSnapshot, asset));
			if (preloadConfig.filter) return assetsRes.filter(preloadConfig.filter);
			return assetsRes;
		}
		if (moduleAssetsInfo) {
			const assetsLength = moduleAssetsInfo.length;
			for (let index = 0; index < assetsLength; index++) {
				const assetsInfo = moduleAssetsInfo[index];
				const exposeFullPath = `${remoteInfo.name}/${assetsInfo.moduleName}`;
				origin.remoteHandler.hooks.lifecycle.handlePreloadModule.emit({
					id: assetsInfo.moduleName === "." ? remoteInfo.name : exposeFullPath,
					name: remoteInfo.name,
					remoteSnapshot: moduleInfoSnapshot,
					preloadConfig,
					remote: remoteInfo,
					origin
				});
				if (require_global.getPreloaded(exposeFullPath)) continue;
				if (preloadConfig.resourceCategory === "all") {
					cssAssets.push(...handleAssets(assetsInfo.assets.css.async));
					cssAssets.push(...handleAssets(assetsInfo.assets.css.sync));
					jsAssets.push(...handleAssets(assetsInfo.assets.js.async));
					jsAssets.push(...handleAssets(assetsInfo.assets.js.sync));
				} else if (preloadConfig.resourceCategory === "sync") {
					cssAssets.push(...handleAssets(assetsInfo.assets.css.sync));
					jsAssets.push(...handleAssets(assetsInfo.assets.js.sync));
				}
				require_global.setPreloaded(exposeFullPath);
			}
		}
	}, true, {}, remoteSnapshot);
	if (remoteSnapshot.shared && remoteSnapshot.shared.length > 0) {
		const collectSharedAssets = (shareInfo, snapshotShared) => {
			const { shared: registeredShared } = require_share.getRegisteredShare(origin.shareScopeMap, snapshotShared.sharedName, shareInfo, origin.sharedHandler.hooks.lifecycle.resolveShare) || {};
			if (registeredShared && typeof registeredShared.lib === "function") {
				snapshotShared.assets.js.sync.forEach((asset) => {
					loadedSharedJsAssets.add(asset);
				});
				snapshotShared.assets.css.sync.forEach((asset) => {
					loadedSharedCssAssets.add(asset);
				});
			}
		};
		remoteSnapshot.shared.forEach((shared) => {
			const shareInfos = options.shared?.[shared.sharedName];
			if (!shareInfos) return;
			const sharedOptions = shared.version ? shareInfos.find((s) => s.version === shared.version) : shareInfos;
			if (!sharedOptions) return;
			require_tool.arrayOptions(sharedOptions).forEach((s) => {
				collectSharedAssets(s, shared);
			});
		});
	}
	const needPreloadJsAssets = jsAssets.filter((asset) => !loadedSharedJsAssets.has(asset) && !isExisted("script", asset));
	return {
		cssAssets: cssAssets.filter((asset) => !loadedSharedCssAssets.has(asset) && !isExisted("link", asset)),
		jsAssetsWithoutEntry: needPreloadJsAssets,
		entryAssets: entryAssets.filter((entry) => !isExisted("script", entry.url))
	};
}
const generatePreloadAssetsPlugin = function() {
	return {
		name: "generate-preload-assets-plugin",
		async generatePreloadAssets(args) {
			const { origin, preloadOptions, remoteInfo, remote, globalSnapshot, remoteSnapshot } = args;
			if (!_module_federation_sdk.isBrowserEnvValue) return {
				cssAssets: [],
				jsAssetsWithoutEntry: [],
				entryAssets: []
			};
			if (require_tool.isRemoteInfoWithEntry(remote) && require_tool.isPureRemoteEntry(remote)) return {
				cssAssets: [],
				jsAssetsWithoutEntry: [],
				entryAssets: [{
					name: remote.name,
					url: remote.entry,
					moduleInfo: {
						name: remoteInfo.name,
						entry: remote.entry,
						type: remoteInfo.type || "global",
						entryGlobalName: "",
						shareScope: ""
					}
				}]
			};
			require_index$1.assignRemoteInfo(remoteInfo, remoteSnapshot);
			return generatePreloadAssets(origin, preloadOptions, remoteInfo, globalSnapshot, remoteSnapshot);
		}
	};
};

//#endregion
exports.generatePreloadAssetsPlugin = generatePreloadAssetsPlugin;
//# sourceMappingURL=generate-preload-assets.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/SnapshotHandler.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
const require_context = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_asyncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs");
const require_asyncWaterfallHooks = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs");
const require_pluginSystem = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/plugins/snapshot/SnapshotHandler.ts
function getGlobalRemoteInfo(moduleInfo, origin) {
	const hostGlobalSnapshot = require_global.getGlobalSnapshotInfoByModuleInfo({
		name: origin.name,
		version: origin.options.version
	});
	const globalRemoteInfo = hostGlobalSnapshot && "remotesInfo" in hostGlobalSnapshot && hostGlobalSnapshot.remotesInfo && require_global.getInfoWithoutType(hostGlobalSnapshot.remotesInfo, moduleInfo.name).value;
	if (globalRemoteInfo && globalRemoteInfo.matchedVersion) return {
		hostGlobalSnapshot,
		globalSnapshot: require_global.getGlobalSnapshot(),
		remoteSnapshot: require_global.getGlobalSnapshotInfoByModuleInfo({
			name: moduleInfo.name,
			version: globalRemoteInfo.matchedVersion
		})
	};
	return {
		hostGlobalSnapshot: void 0,
		globalSnapshot: require_global.getGlobalSnapshot(),
		remoteSnapshot: require_global.getGlobalSnapshotInfoByModuleInfo({
			name: moduleInfo.name,
			version: "version" in moduleInfo ? moduleInfo.version : void 0
		})
	};
}
var SnapshotHandler = class {
	constructor(HostInstance) {
		this.loadingHostSnapshot = null;
		this.manifestCache = /* @__PURE__ */ new Map();
		this.hooks = new require_pluginSystem.PluginSystem({
			beforeLoadRemoteSnapshot: new require_asyncHook.AsyncHook("beforeLoadRemoteSnapshot"),
			loadSnapshot: new require_asyncWaterfallHooks.AsyncWaterfallHook("loadGlobalSnapshot"),
			loadRemoteSnapshot: new require_asyncWaterfallHooks.AsyncWaterfallHook("loadRemoteSnapshot"),
			afterLoadSnapshot: new require_asyncWaterfallHooks.AsyncWaterfallHook("afterLoadSnapshot")
		});
		this.manifestLoading = require_global.Global.__FEDERATION__.__MANIFEST_LOADING__;
		this.HostInstance = HostInstance;
		this.loaderHook = HostInstance.loaderHook;
	}
	async loadRemoteSnapshotInfo({ moduleInfo, id, expose }) {
		const { options } = this.HostInstance;
		await this.hooks.lifecycle.beforeLoadRemoteSnapshot.emit({
			options,
			moduleInfo
		});
		let hostSnapshot = require_global.getGlobalSnapshotInfoByModuleInfo({
			name: this.HostInstance.options.name,
			version: this.HostInstance.options.version
		});
		if (!hostSnapshot) {
			hostSnapshot = {
				version: this.HostInstance.options.version || "",
				remoteEntry: "",
				remotesInfo: {}
			};
			require_global.addGlobalSnapshot({ [this.HostInstance.options.name]: hostSnapshot });
		}
		if (hostSnapshot && "remotesInfo" in hostSnapshot && !require_global.getInfoWithoutType(hostSnapshot.remotesInfo, moduleInfo.name).value) {
			if ("version" in moduleInfo || "entry" in moduleInfo) hostSnapshot.remotesInfo = {
				...hostSnapshot?.remotesInfo,
				[moduleInfo.name]: { matchedVersion: "version" in moduleInfo ? moduleInfo.version : moduleInfo.entry }
			};
		}
		const { hostGlobalSnapshot, remoteSnapshot, globalSnapshot } = this.getGlobalRemoteInfo(moduleInfo);
		const { remoteSnapshot: globalRemoteSnapshot, globalSnapshot: globalSnapshotRes } = await this.hooks.lifecycle.loadSnapshot.emit({
			options,
			moduleInfo,
			hostGlobalSnapshot,
			remoteSnapshot,
			globalSnapshot
		});
		let mSnapshot;
		let gSnapshot;
		if (globalRemoteSnapshot) if ((0, _module_federation_sdk.isManifestProvider)(globalRemoteSnapshot)) {
			const remoteEntry = _module_federation_sdk.isBrowserEnvValue ? globalRemoteSnapshot.remoteEntry : globalRemoteSnapshot.ssrRemoteEntry || globalRemoteSnapshot.remoteEntry || "";
			const moduleSnapshot = await this.getManifestJson(remoteEntry, moduleInfo, {});
			const globalSnapshotRes = require_global.setGlobalSnapshotInfoByModuleInfo({
				...moduleInfo,
				entry: remoteEntry
			}, moduleSnapshot);
			mSnapshot = moduleSnapshot;
			gSnapshot = globalSnapshotRes;
		} else {
			const { remoteSnapshot: remoteSnapshotRes } = await this.hooks.lifecycle.loadRemoteSnapshot.emit({
				options: this.HostInstance.options,
				moduleInfo,
				remoteSnapshot: globalRemoteSnapshot,
				from: "global"
			});
			mSnapshot = remoteSnapshotRes;
			gSnapshot = globalSnapshotRes;
		}
		else if (require_tool.isRemoteInfoWithEntry(moduleInfo)) {
			const moduleSnapshot = await this.getManifestJson(moduleInfo.entry, moduleInfo, {});
			const globalSnapshotRes = require_global.setGlobalSnapshotInfoByModuleInfo(moduleInfo, moduleSnapshot);
			const { remoteSnapshot: remoteSnapshotRes } = await this.hooks.lifecycle.loadRemoteSnapshot.emit({
				options: this.HostInstance.options,
				moduleInfo,
				remoteSnapshot: moduleSnapshot,
				from: "global"
			});
			mSnapshot = remoteSnapshotRes;
			gSnapshot = globalSnapshotRes;
		} else require_logger.error(_module_federation_error_codes.RUNTIME_007, _module_federation_error_codes.runtimeDescMap, {
			remoteName: moduleInfo.name,
			remoteVersion: moduleInfo.version,
			hostName: this.HostInstance.options.name,
			globalSnapshot: JSON.stringify(globalSnapshotRes)
		}, void 0, require_context.optionsToMFContext(this.HostInstance.options));
		await this.hooks.lifecycle.afterLoadSnapshot.emit({
			id,
			host: this.HostInstance,
			options,
			moduleInfo,
			remoteSnapshot: mSnapshot
		});
		return {
			remoteSnapshot: mSnapshot,
			globalSnapshot: gSnapshot
		};
	}
	getGlobalRemoteInfo(moduleInfo) {
		return getGlobalRemoteInfo(moduleInfo, this.HostInstance);
	}
	async getManifestJson(manifestUrl, moduleInfo, extraOptions) {
		const getManifest = async () => {
			let manifestJson = this.manifestCache.get(manifestUrl);
			if (manifestJson) return manifestJson;
			try {
				let res = await this.loaderHook.lifecycle.fetch.emit(manifestUrl, {}, require_load.getRemoteInfo(moduleInfo));
				if (!res || !(res instanceof Response)) res = await fetch(manifestUrl, {});
				manifestJson = await res.json();
			} catch (err) {
				manifestJson = await this.HostInstance.remoteHandler.hooks.lifecycle.errorLoadRemote.emit({
					id: manifestUrl,
					error: err,
					from: "runtime",
					lifecycle: "afterResolve",
					origin: this.HostInstance
				});
				if (!manifestJson) {
					delete this.manifestLoading[manifestUrl];
					require_logger.error(_module_federation_error_codes.RUNTIME_003, _module_federation_error_codes.runtimeDescMap, {
						manifestUrl,
						moduleName: moduleInfo.name,
						hostName: this.HostInstance.options.name
					}, `${err}`, require_context.optionsToMFContext(this.HostInstance.options));
				}
			}
			require_logger.assert(manifestJson.metaData && manifestJson.exposes && manifestJson.shared, `"${manifestUrl}" is not a valid federation manifest for remote "${moduleInfo.name}". Missing required fields: ${[
				!manifestJson.metaData && "metaData",
				!manifestJson.exposes && "exposes",
				!manifestJson.shared && "shared"
			].filter(Boolean).join(", ")}.`);
			this.manifestCache.set(manifestUrl, manifestJson);
			return manifestJson;
		};
		const asyncLoadProcess = async () => {
			const manifestJson = await getManifest();
			const remoteSnapshot = (0, _module_federation_sdk.generateSnapshotFromManifest)(manifestJson, { version: manifestUrl });
			const { remoteSnapshot: remoteSnapshotRes } = await this.hooks.lifecycle.loadRemoteSnapshot.emit({
				options: this.HostInstance.options,
				moduleInfo,
				manifestJson,
				remoteSnapshot,
				manifestUrl,
				from: "manifest"
			});
			return remoteSnapshotRes;
		};
		if (!this.manifestLoading[manifestUrl]) this.manifestLoading[manifestUrl] = asyncLoadProcess().then((res) => res);
		return this.manifestLoading[manifestUrl];
	}
};

//#endregion
exports.SnapshotHandler = SnapshotHandler;
exports.getGlobalRemoteInfo = getGlobalRemoteInfo;
//# sourceMappingURL=SnapshotHandler.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_preload = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/preload.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/plugins/snapshot/index.ts
function assignRemoteInfo(remoteInfo, remoteSnapshot) {
	const remoteEntryInfo = require_tool.getRemoteEntryInfoFromSnapshot(remoteSnapshot);
	if (!remoteEntryInfo.url) require_logger.error(_module_federation_error_codes.RUNTIME_011, _module_federation_error_codes.runtimeDescMap, { remoteName: remoteInfo.name });
	let entryUrl = (0, _module_federation_sdk.getResourceUrl)(remoteSnapshot, remoteEntryInfo.url);
	if (!_module_federation_sdk.isBrowserEnvValue && !entryUrl.startsWith("http")) entryUrl = `https:${entryUrl}`;
	remoteInfo.type = remoteEntryInfo.type;
	remoteInfo.entryGlobalName = remoteEntryInfo.globalName;
	remoteInfo.entry = entryUrl;
	remoteInfo.version = remoteSnapshot.version;
	remoteInfo.buildVersion = remoteSnapshot.buildVersion;
}
function snapshotPlugin() {
	return {
		name: "snapshot-plugin",
		async afterResolve(args) {
			const { remote, pkgNameOrAlias, expose, origin, remoteInfo, id } = args;
			if (!require_tool.isRemoteInfoWithEntry(remote) || !require_tool.isPureRemoteEntry(remote)) {
				const { remoteSnapshot, globalSnapshot } = await origin.snapshotHandler.loadRemoteSnapshotInfo({
					moduleInfo: remote,
					id
				});
				assignRemoteInfo(remoteInfo, remoteSnapshot);
				const preloadOptions = {
					remote,
					preloadConfig: {
						nameOrAlias: pkgNameOrAlias,
						exposes: [expose],
						resourceCategory: "sync",
						share: false,
						depsRemote: false
					}
				};
				const assets = await origin.remoteHandler.hooks.lifecycle.generatePreloadAssets.emit({
					origin,
					preloadOptions,
					remoteInfo,
					remote,
					remoteSnapshot,
					globalSnapshot
				});
				if (assets) require_preload.preloadAssets(remoteInfo, origin, assets, false);
				return {
					...args,
					remoteSnapshot
				};
			}
			return args;
		}
	};
}

//#endregion
exports.assignRemoteInfo = assignRemoteInfo;
exports.snapshotPlugin = snapshotPlugin;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/remote/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/constant.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
const require_manifest = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
const require_context = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_preload = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/preload.cjs");
const require_index$1 = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/module/index.cjs");
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");
const require_asyncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs");
const require_syncWaterfallHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs");
const require_asyncWaterfallHooks = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs");
const require_pluginSystem = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/index.cjs");
const require_SnapshotHandler = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/plugins/snapshot/SnapshotHandler.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/remote/index.ts
var RemoteHandler = class {
	constructor(host) {
		this.hooks = new require_pluginSystem.PluginSystem({
			beforeRegisterRemote: new require_syncWaterfallHook.SyncWaterfallHook("beforeRegisterRemote"),
			registerRemote: new require_syncWaterfallHook.SyncWaterfallHook("registerRemote"),
			beforeRequest: new require_asyncWaterfallHooks.AsyncWaterfallHook("beforeRequest"),
			onLoad: new require_asyncHook.AsyncHook("onLoad"),
			handlePreloadModule: new require_syncHook.SyncHook("handlePreloadModule"),
			errorLoadRemote: new require_asyncHook.AsyncHook("errorLoadRemote"),
			beforePreloadRemote: new require_asyncHook.AsyncHook("beforePreloadRemote"),
			generatePreloadAssets: new require_asyncHook.AsyncHook("generatePreloadAssets"),
			afterPreloadRemote: new require_asyncHook.AsyncHook(),
			loadEntry: new require_asyncHook.AsyncHook()
		});
		this.host = host;
		this.idToRemoteMap = {};
	}
	formatAndRegisterRemote(globalOptions, userOptions) {
		return (userOptions.remotes || []).reduce((res, remote) => {
			this.registerRemote(remote, res, { force: false });
			return res;
		}, globalOptions.remotes);
	}
	setIdToRemoteMap(id, remoteMatchInfo) {
		const { remote, expose } = remoteMatchInfo;
		const { name, alias } = remote;
		this.idToRemoteMap[id] = {
			name: remote.name,
			expose
		};
		if (alias && id.startsWith(name)) {
			const idWithAlias = id.replace(name, alias);
			this.idToRemoteMap[idWithAlias] = {
				name: remote.name,
				expose
			};
			return;
		}
		if (alias && id.startsWith(alias)) {
			const idWithName = id.replace(alias, name);
			this.idToRemoteMap[idWithName] = {
				name: remote.name,
				expose
			};
		}
	}
	async loadRemote(id, options) {
		const { host } = this;
		try {
			const { loadFactory = true } = options || { loadFactory: true };
			const { module, moduleOptions, remoteMatchInfo } = await this.getRemoteModuleAndOptions({ id });
			const { pkgNameOrAlias, remote, expose, id: idRes, remoteSnapshot } = remoteMatchInfo;
			const moduleOrFactory = await module.get(idRes, expose, options, remoteSnapshot);
			const moduleWrapper = await this.hooks.lifecycle.onLoad.emit({
				id: idRes,
				pkgNameOrAlias,
				expose,
				exposeModule: loadFactory ? moduleOrFactory : void 0,
				exposeModuleFactory: loadFactory ? void 0 : moduleOrFactory,
				remote,
				options: moduleOptions,
				moduleInstance: module,
				origin: host
			});
			this.setIdToRemoteMap(id, remoteMatchInfo);
			if (typeof moduleWrapper === "function") return moduleWrapper;
			return moduleOrFactory;
		} catch (error) {
			const { from = "runtime" } = options || { from: "runtime" };
			const failOver = await this.hooks.lifecycle.errorLoadRemote.emit({
				id,
				error,
				from,
				lifecycle: "onLoad",
				origin: host
			});
			if (!failOver) throw error;
			return failOver;
		}
	}
	async preloadRemote(preloadOptions) {
		const { host } = this;
		await this.hooks.lifecycle.beforePreloadRemote.emit({
			preloadOps: preloadOptions,
			options: host.options,
			origin: host
		});
		const preloadOps = require_preload.formatPreloadArgs(host.options.remotes, preloadOptions);
		await Promise.all(preloadOps.map(async (ops) => {
			const { remote } = ops;
			const remoteInfo = require_load.getRemoteInfo(remote);
			const { globalSnapshot, remoteSnapshot } = await host.snapshotHandler.loadRemoteSnapshotInfo({ moduleInfo: remote });
			const assets = await this.hooks.lifecycle.generatePreloadAssets.emit({
				origin: host,
				preloadOptions: ops,
				remote,
				remoteInfo,
				globalSnapshot,
				remoteSnapshot
			});
			if (!assets) return;
			require_preload.preloadAssets(remoteInfo, host, assets);
		}));
	}
	registerRemotes(remotes, options) {
		const { host } = this;
		remotes.forEach((remote) => {
			this.registerRemote(remote, host.options.remotes, { force: options?.force });
		});
	}
	async getRemoteModuleAndOptions(options) {
		const { host } = this;
		const { id } = options;
		let loadRemoteArgs;
		try {
			loadRemoteArgs = await this.hooks.lifecycle.beforeRequest.emit({
				id,
				options: host.options,
				origin: host
			});
		} catch (error) {
			loadRemoteArgs = await this.hooks.lifecycle.errorLoadRemote.emit({
				id,
				options: host.options,
				origin: host,
				from: "runtime",
				error,
				lifecycle: "beforeRequest"
			});
			if (!loadRemoteArgs) throw error;
		}
		const { id: idRes } = loadRemoteArgs;
		const remoteSplitInfo = require_manifest.matchRemoteWithNameAndExpose(host.options.remotes, idRes);
		if (!remoteSplitInfo) require_logger.error(_module_federation_error_codes.RUNTIME_004, _module_federation_error_codes.runtimeDescMap, {
			hostName: host.options.name,
			requestId: idRes
		}, void 0, require_context.optionsToMFContext(host.options));
		const { remote: rawRemote } = remoteSplitInfo;
		const remoteInfo = require_load.getRemoteInfo(rawRemote);
		const matchInfo = await host.sharedHandler.hooks.lifecycle.afterResolve.emit({
			id: idRes,
			...remoteSplitInfo,
			options: host.options,
			origin: host,
			remoteInfo
		});
		const { remote, expose } = matchInfo;
		require_logger.assert(remote && expose, `The 'beforeRequest' hook was executed, but it failed to return the correct 'remote' and 'expose' values while loading ${idRes}.`);
		let module = host.moduleCache.get(remote.name);
		const moduleOptions = {
			host,
			remoteInfo
		};
		if (!module) {
			module = new require_index$1.Module(moduleOptions);
			host.moduleCache.set(remote.name, module);
		}
		return {
			module,
			moduleOptions,
			remoteMatchInfo: matchInfo
		};
	}
	registerRemote(remote, targetRemotes, options) {
		const { host } = this;
		const normalizeRemote = () => {
			if (remote.alias) {
				const findEqual = targetRemotes.find((item) => remote.alias && (item.name.startsWith(remote.alias) || item.alias?.startsWith(remote.alias)));
				require_logger.assert(!findEqual, `The alias ${remote.alias} of remote ${remote.name} is not allowed to be the prefix of ${findEqual && findEqual.name} name or alias`);
			}
			if ("entry" in remote) {
				if (_module_federation_sdk.isBrowserEnvValue && typeof window !== "undefined" && !remote.entry.startsWith("http")) remote.entry = new URL(remote.entry, window.location.origin).href;
			}
			if (!remote.shareScope) remote.shareScope = require_constant.DEFAULT_SCOPE;
			if (!remote.type) remote.type = require_constant.DEFAULT_REMOTE_TYPE;
		};
		this.hooks.lifecycle.beforeRegisterRemote.emit({
			remote,
			origin: host
		});
		const registeredRemote = targetRemotes.find((item) => item.name === remote.name);
		if (!registeredRemote) {
			normalizeRemote();
			targetRemotes.push(remote);
			this.hooks.lifecycle.registerRemote.emit({
				remote,
				origin: host
			});
		} else {
			const messages = [`The remote "${remote.name}" is already registered.`, "Please note that overriding it may cause unexpected errors."];
			if (options?.force) {
				this.removeRemote(registeredRemote);
				normalizeRemote();
				targetRemotes.push(remote);
				this.hooks.lifecycle.registerRemote.emit({
					remote,
					origin: host
				});
				(0, _module_federation_sdk.warn)(messages.join(" "));
			}
		}
	}
	removeRemote(remote) {
		try {
			const { host } = this;
			const { name } = remote;
			const remoteIndex = host.options.remotes.findIndex((item) => item.name === name);
			if (remoteIndex !== -1) host.options.remotes.splice(remoteIndex, 1);
			const loadedModule = host.moduleCache.get(remote.name);
			if (loadedModule) {
				const remoteInfo = loadedModule.remoteInfo;
				const key = remoteInfo.entryGlobalName;
				if (require_global.CurrentGlobal[key]) if (Object.getOwnPropertyDescriptor(require_global.CurrentGlobal, key)?.configurable) delete require_global.CurrentGlobal[key];
				else require_global.CurrentGlobal[key] = void 0;
				const remoteEntryUniqueKey = require_load.getRemoteEntryUniqueKey(loadedModule.remoteInfo);
				if (require_global.globalLoading[remoteEntryUniqueKey]) delete require_global.globalLoading[remoteEntryUniqueKey];
				host.snapshotHandler.manifestCache.delete(remoteInfo.entry);
				let remoteInsId = remoteInfo.buildVersion ? (0, _module_federation_sdk.composeKeyWithSeparator)(remoteInfo.name, remoteInfo.buildVersion) : remoteInfo.name;
				const remoteInsIndex = require_global.CurrentGlobal.__FEDERATION__.__INSTANCES__.findIndex((ins) => {
					if (remoteInfo.buildVersion) return ins.options.id === remoteInsId;
					else return ins.name === remoteInsId;
				});
				if (remoteInsIndex !== -1) {
					const remoteIns = require_global.CurrentGlobal.__FEDERATION__.__INSTANCES__[remoteInsIndex];
					remoteInsId = remoteIns.options.id || remoteInsId;
					const globalShareScopeMap = require_share.getGlobalShareScope();
					let isAllSharedNotUsed = true;
					const needDeleteKeys = [];
					Object.keys(globalShareScopeMap).forEach((instId) => {
						const shareScopeMap = globalShareScopeMap[instId];
						shareScopeMap && Object.keys(shareScopeMap).forEach((shareScope) => {
							const shareScopeVal = shareScopeMap[shareScope];
							shareScopeVal && Object.keys(shareScopeVal).forEach((shareName) => {
								const sharedPkgs = shareScopeVal[shareName];
								sharedPkgs && Object.keys(sharedPkgs).forEach((shareVersion) => {
									const shared = sharedPkgs[shareVersion];
									if (shared && typeof shared === "object" && shared.from === remoteInfo.name) if (shared.loaded || shared.loading) {
										shared.useIn = shared.useIn.filter((usedHostName) => usedHostName !== remoteInfo.name);
										if (shared.useIn.length) isAllSharedNotUsed = false;
										else needDeleteKeys.push([
											instId,
											shareScope,
											shareName,
											shareVersion
										]);
									} else needDeleteKeys.push([
										instId,
										shareScope,
										shareName,
										shareVersion
									]);
								});
							});
						});
					});
					if (isAllSharedNotUsed) {
						remoteIns.shareScopeMap = {};
						delete globalShareScopeMap[remoteInsId];
					}
					needDeleteKeys.forEach(([insId, shareScope, shareName, shareVersion]) => {
						delete globalShareScopeMap[insId]?.[shareScope]?.[shareName]?.[shareVersion];
					});
					require_global.CurrentGlobal.__FEDERATION__.__INSTANCES__.splice(remoteInsIndex, 1);
				}
				const { hostGlobalSnapshot } = require_SnapshotHandler.getGlobalRemoteInfo(remote, host);
				if (hostGlobalSnapshot) {
					const remoteKey = hostGlobalSnapshot && "remotesInfo" in hostGlobalSnapshot && hostGlobalSnapshot.remotesInfo && require_global.getInfoWithoutType(hostGlobalSnapshot.remotesInfo, remote.name).key;
					if (remoteKey) {
						delete hostGlobalSnapshot.remotesInfo[remoteKey];
						if (Boolean(require_global.Global.__FEDERATION__.__MANIFEST_LOADING__[remoteKey])) delete require_global.Global.__FEDERATION__.__MANIFEST_LOADING__[remoteKey];
					}
				}
				host.moduleCache.delete(remote.name);
			}
		} catch (err) {
			require_logger.logger.error(`removeRemote failed: ${err instanceof Error ? err.message : String(err)}`);
		}
	}
};

//#endregion
exports.RemoteHandler = RemoteHandler;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/shared/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/constant.cjs");
const require_share = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs");
const require_context = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");
const require_asyncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs");
const require_syncWaterfallHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs");
const require_asyncWaterfallHooks = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs");
const require_pluginSystem = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/shared/index.ts
var SharedHandler = class {
	constructor(host) {
		this.hooks = new require_pluginSystem.PluginSystem({
			beforeRegisterShare: new require_syncWaterfallHook.SyncWaterfallHook("beforeRegisterShare"),
			afterResolve: new require_asyncWaterfallHooks.AsyncWaterfallHook("afterResolve"),
			beforeLoadShare: new require_asyncWaterfallHooks.AsyncWaterfallHook("beforeLoadShare"),
			loadShare: new require_asyncHook.AsyncHook(),
			resolveShare: new require_syncWaterfallHook.SyncWaterfallHook("resolveShare"),
			initContainerShareScopeMap: new require_syncWaterfallHook.SyncWaterfallHook("initContainerShareScopeMap")
		});
		this.host = host;
		this.shareScopeMap = {};
		this.initTokens = {};
		this._setGlobalShareScopeMap(host.options);
	}
	registerShared(globalOptions, userOptions) {
		const { newShareInfos, allShareInfos } = require_share.formatShareConfigs(globalOptions, userOptions);
		Object.keys(newShareInfos).forEach((sharedKey) => {
			newShareInfos[sharedKey].forEach((sharedVal) => {
				sharedVal.scope.forEach((sc) => {
					this.hooks.lifecycle.beforeRegisterShare.emit({
						origin: this.host,
						pkgName: sharedKey,
						shared: sharedVal
					});
					if (!this.shareScopeMap[sc]?.[sharedKey]) this.setShared({
						pkgName: sharedKey,
						lib: sharedVal.lib,
						get: sharedVal.get,
						loaded: sharedVal.loaded || Boolean(sharedVal.lib),
						shared: sharedVal,
						from: userOptions.name
					});
				});
			});
		});
		return {
			newShareInfos,
			allShareInfos
		};
	}
	async loadShare(pkgName, extraOptions) {
		const { host } = this;
		const shareOptions = require_share.getTargetSharedOptions({
			pkgName,
			extraOptions,
			shareInfos: host.options.shared
		});
		if (shareOptions?.scope) await Promise.all(shareOptions.scope.map(async (shareScope) => {
			await Promise.all(this.initializeSharing(shareScope, { strategy: shareOptions.strategy }));
		}));
		const { shareInfo: shareOptionsRes } = await this.hooks.lifecycle.beforeLoadShare.emit({
			pkgName,
			shareInfo: shareOptions,
			shared: host.options.shared,
			origin: host
		});
		require_logger.assert(shareOptionsRes, `Cannot find shared "${pkgName}" in host "${host.options.name}". Ensure the shared config for "${pkgName}" is declared in the federation plugin options and the host has been initialized before loading shares.`);
		const { shared: registeredShared, useTreesShaking } = require_share.getRegisteredShare(this.shareScopeMap, pkgName, shareOptionsRes, this.hooks.lifecycle.resolveShare) || {};
		if (registeredShared) {
			const targetShared = require_share.directShare(registeredShared, useTreesShaking);
			if (targetShared.lib) {
				require_share.addUseIn(targetShared, host.options.name);
				return targetShared.lib;
			} else if (targetShared.loading && !targetShared.loaded) {
				const factory = await targetShared.loading;
				targetShared.loaded = true;
				if (!targetShared.lib) targetShared.lib = factory;
				require_share.addUseIn(targetShared, host.options.name);
				return factory;
			} else {
				const asyncLoadProcess = async () => {
					const factory = await targetShared.get();
					require_share.addUseIn(targetShared, host.options.name);
					targetShared.loaded = true;
					targetShared.lib = factory;
					return factory;
				};
				const loading = asyncLoadProcess();
				this.setShared({
					pkgName,
					loaded: false,
					shared: registeredShared,
					from: host.options.name,
					lib: null,
					loading,
					treeShaking: useTreesShaking ? targetShared : void 0
				});
				return loading;
			}
		} else {
			if (extraOptions?.customShareInfo) return false;
			const _useTreeShaking = require_share.shouldUseTreeShaking(shareOptionsRes.treeShaking);
			const targetShared = require_share.directShare(shareOptionsRes, _useTreeShaking);
			const asyncLoadProcess = async () => {
				const factory = await targetShared.get();
				targetShared.lib = factory;
				targetShared.loaded = true;
				require_share.addUseIn(targetShared, host.options.name);
				const { shared: gShared, useTreesShaking: gUseTreeShaking } = require_share.getRegisteredShare(this.shareScopeMap, pkgName, shareOptionsRes, this.hooks.lifecycle.resolveShare) || {};
				if (gShared) {
					const targetGShared = require_share.directShare(gShared, gUseTreeShaking);
					targetGShared.lib = factory;
					targetGShared.loaded = true;
					gShared.from = shareOptionsRes.from;
				}
				return factory;
			};
			const loading = asyncLoadProcess();
			this.setShared({
				pkgName,
				loaded: false,
				shared: shareOptionsRes,
				from: host.options.name,
				lib: null,
				loading,
				treeShaking: _useTreeShaking ? targetShared : void 0
			});
			return loading;
		}
	}
	/**
	* This function initializes the sharing sequence (executed only once per share scope).
	* It accepts one argument, the name of the share scope.
	* If the share scope does not exist, it creates one.
	*/
	initializeSharing(shareScopeName = require_constant.DEFAULT_SCOPE, extraOptions) {
		const { host } = this;
		const from = extraOptions?.from;
		const strategy = extraOptions?.strategy;
		let initScope = extraOptions?.initScope;
		const promises = [];
		if (from !== "build") {
			const { initTokens } = this;
			if (!initScope) initScope = [];
			let initToken = initTokens[shareScopeName];
			if (!initToken) initToken = initTokens[shareScopeName] = { from: this.host.name };
			if (initScope.indexOf(initToken) >= 0) return promises;
			initScope.push(initToken);
		}
		const shareScope = this.shareScopeMap;
		const hostName = host.options.name;
		if (!shareScope[shareScopeName]) shareScope[shareScopeName] = {};
		const scope = shareScope[shareScopeName];
		const register = (name, shared) => {
			const { version, eager } = shared;
			scope[name] = scope[name] || {};
			const versions = scope[name];
			const activeVersion = versions[version] && require_share.directShare(versions[version]);
			const activeVersionEager = Boolean(activeVersion && ("eager" in activeVersion && activeVersion.eager || "shareConfig" in activeVersion && activeVersion.shareConfig?.eager));
			if (!activeVersion || activeVersion.strategy !== "loaded-first" && !activeVersion.loaded && (Boolean(!eager) !== !activeVersionEager ? eager : hostName > versions[version].from)) versions[version] = shared;
		};
		const initRemoteModule = async (key) => {
			const { module } = await host.remoteHandler.getRemoteModuleAndOptions({ id: key });
			let remoteEntryExports = void 0;
			try {
				remoteEntryExports = await module.getEntry();
			} catch (error) {
				remoteEntryExports = await host.remoteHandler.hooks.lifecycle.errorLoadRemote.emit({
					id: key,
					error,
					from: "runtime",
					lifecycle: "beforeLoadShare",
					origin: host
				});
				if (!remoteEntryExports) return;
			} finally {
				if (remoteEntryExports?.init && !module.initing) {
					module.remoteEntryExports = remoteEntryExports;
					await module.init(void 0, void 0, initScope);
				}
			}
		};
		Object.keys(host.options.shared).forEach((shareName) => {
			host.options.shared[shareName].forEach((shared) => {
				if (shared.scope.includes(shareScopeName)) register(shareName, shared);
			});
		});
		if (host.options.shareStrategy === "version-first" || strategy === "version-first") host.options.remotes.forEach((remote) => {
			if (remote.shareScope === shareScopeName) promises.push(initRemoteModule(remote.name));
		});
		return promises;
	}
	loadShareSync(pkgName, extraOptions) {
		const { host } = this;
		const shareOptions = require_share.getTargetSharedOptions({
			pkgName,
			extraOptions,
			shareInfos: host.options.shared
		});
		if (shareOptions?.scope) shareOptions.scope.forEach((shareScope) => {
			this.initializeSharing(shareScope, { strategy: shareOptions.strategy });
		});
		const { shared: registeredShared, useTreesShaking } = require_share.getRegisteredShare(this.shareScopeMap, pkgName, shareOptions, this.hooks.lifecycle.resolveShare) || {};
		if (registeredShared) {
			if (typeof registeredShared.lib === "function") {
				require_share.addUseIn(registeredShared, host.options.name);
				if (!registeredShared.loaded) {
					registeredShared.loaded = true;
					if (registeredShared.from === host.options.name) shareOptions.loaded = true;
				}
				return registeredShared.lib;
			}
			if (typeof registeredShared.get === "function") {
				const module = registeredShared.get();
				if (!(module instanceof Promise)) {
					require_share.addUseIn(registeredShared, host.options.name);
					this.setShared({
						pkgName,
						loaded: true,
						from: host.options.name,
						lib: module,
						shared: registeredShared
					});
					return module;
				}
			}
		}
		if (shareOptions.lib) {
			if (!shareOptions.loaded) shareOptions.loaded = true;
			return shareOptions.lib;
		}
		if (shareOptions.get) {
			const module = shareOptions.get();
			if (module instanceof Promise) require_logger.error(extraOptions?.from === "build" ? _module_federation_error_codes.RUNTIME_005 : _module_federation_error_codes.RUNTIME_006, _module_federation_error_codes.runtimeDescMap, {
				hostName: host.options.name,
				sharedPkgName: pkgName
			}, void 0, require_context.optionsToMFContext(host.options));
			shareOptions.lib = module;
			this.setShared({
				pkgName,
				loaded: true,
				from: host.options.name,
				lib: shareOptions.lib,
				shared: shareOptions
			});
			return shareOptions.lib;
		}
		require_logger.error(_module_federation_error_codes.RUNTIME_006, _module_federation_error_codes.runtimeDescMap, {
			hostName: host.options.name,
			sharedPkgName: pkgName
		}, void 0, require_context.optionsToMFContext(host.options));
	}
	initShareScopeMap(scopeName, shareScope, extraOptions = {}) {
		const { host } = this;
		this.shareScopeMap[scopeName] = shareScope;
		this.hooks.lifecycle.initContainerShareScopeMap.emit({
			shareScope,
			options: host.options,
			origin: host,
			scopeName,
			hostShareScopeMap: extraOptions.hostShareScopeMap
		});
	}
	setShared({ pkgName, shared, from, lib, loading, loaded, get, treeShaking }) {
		const { version, scope = "default", ...shareInfo } = shared;
		const scopes = Array.isArray(scope) ? scope : [scope];
		const mergeAttrs = (shared) => {
			const merge = (s, key, val) => {
				if (val && !s[key]) s[key] = val;
			};
			const targetShared = treeShaking ? shared.treeShaking : shared;
			merge(targetShared, "loaded", loaded);
			merge(targetShared, "loading", loading);
			merge(targetShared, "get", get);
		};
		scopes.forEach((sc) => {
			if (!this.shareScopeMap[sc]) this.shareScopeMap[sc] = {};
			if (!this.shareScopeMap[sc][pkgName]) this.shareScopeMap[sc][pkgName] = {};
			if (!this.shareScopeMap[sc][pkgName][version]) this.shareScopeMap[sc][pkgName][version] = {
				version,
				scope: [sc],
				...shareInfo,
				lib
			};
			const registeredShared = this.shareScopeMap[sc][pkgName][version];
			mergeAttrs(registeredShared);
			if (from && registeredShared.from !== from) registeredShared.from = from;
		});
	}
	_setGlobalShareScopeMap(hostOptions) {
		const globalShareScopeMap = require_share.getGlobalShareScope();
		const identifier = hostOptions.id || hostOptions.name;
		if (identifier && !globalShareScopeMap[identifier]) globalShareScopeMap[identifier] = this.shareScopeMap;
	}
};

//#endregion
exports.SharedHandler = SharedHandler;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/type/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/_virtual/_rolldown/runtime.cjs");

//#region src/type/index.ts
var type_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "type_exports", ({
  enumerable: true,
  get: function () {
    return type_exports;
  }
}));
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs"(__unused_rspack_module, exports) {

//#region src/utils/context.ts
function remoteToEntry(r) {
	return {
		name: r.name,
		alias: r.alias,
		entry: "entry" in r ? r.entry : void 0,
		version: "version" in r ? r.version : void 0,
		type: r.type,
		entryGlobalName: r.entryGlobalName,
		shareScope: r.shareScope
	};
}
/**
* Build a partial MFContext from runtime Options.
* Used to enrich diagnostic entries with host context at error sites.
*/
function optionsToMFContext(options) {
	const shared = {};
	for (const [pkgName, versions] of Object.entries(options.shared)) {
		const first = versions[0];
		if (first) shared[pkgName] = {
			version: first.version,
			singleton: first.shareConfig?.singleton,
			requiredVersion: first.shareConfig?.requiredVersion === false ? false : first.shareConfig?.requiredVersion,
			eager: first.eager,
			strictVersion: first.shareConfig?.strictVersion
		};
	}
	return {
		project: {
			name: options.name,
			mfRole: options.remotes?.length > 0 ? "host" : "unknown"
		},
		mfConfig: {
			name: options.name,
			remotes: options.remotes?.map(remoteToEntry) ?? [],
			shared
		}
	};
}

//#endregion
exports.optionsToMFContext = optionsToMFContext;
//# sourceMappingURL=context.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/env.cjs"(__unused_rspack_module, exports, __webpack_require__) {
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/utils/env.ts
function getBuilderId() {
	return typeof FEDERATION_BUILD_IDENTIFIER !== "undefined" ? FEDERATION_BUILD_IDENTIFIER : "";
}

//#endregion
exports.getBuilderId = getBuilderId;
//# sourceMappingURL=env.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");

//#region src/utils/hooks/asyncHook.ts
var AsyncHook = class extends require_syncHook.SyncHook {
	emit(...data) {
		let result;
		const ls = Array.from(this.listeners);
		if (ls.length > 0) {
			let i = 0;
			const call = (prev) => {
				if (prev === false) return false;
				else if (i < ls.length) return Promise.resolve(ls[i++].apply(null, data)).then(call);
				else return prev;
			};
			result = call();
		}
		return Promise.resolve(result);
	}
};

//#endregion
exports.AsyncHook = AsyncHook;
//# sourceMappingURL=asyncHook.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");
const require_syncWaterfallHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs");

//#region src/utils/hooks/asyncWaterfallHooks.ts
var AsyncWaterfallHook = class extends require_syncHook.SyncHook {
	constructor(type) {
		super();
		this.onerror = require_logger.error;
		this.type = type;
	}
	emit(data) {
		if (!require_tool.isObject(data)) require_logger.error(`The response data for the "${this.type}" hook must be an object.`);
		const ls = Array.from(this.listeners);
		if (ls.length > 0) {
			let i = 0;
			const processError = (e) => {
				require_logger.warn(e);
				this.onerror(e);
				return data;
			};
			const call = (prevData) => {
				if (require_syncWaterfallHook.checkReturnData(data, prevData)) {
					data = prevData;
					if (i < ls.length) try {
						return Promise.resolve(ls[i++](data)).then(call, processError);
					} catch (e) {
						return processError(e);
					}
				} else this.onerror(`A plugin returned an incorrect value for the "${this.type}" type.`);
				return data;
			};
			return Promise.resolve(call(data));
		}
		return Promise.resolve(data);
	}
};

//#endregion
exports.AsyncWaterfallHook = AsyncWaterfallHook;
//# sourceMappingURL=asyncWaterfallHooks.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/index.cjs"(__unused_rspack_module, __unused_rspack_exports, __webpack_require__) {
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");
const require_asyncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncHook.cjs");
const require_syncWaterfallHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs");
const require_asyncWaterfallHooks = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/asyncWaterfallHooks.cjs");
const require_pluginSystem = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs");


},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/pluginSystem.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
__webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs");

//#region src/utils/hooks/pluginSystem.ts
var PluginSystem = class {
	constructor(lifecycle) {
		this.registerPlugins = {};
		this.lifecycle = lifecycle;
		this.lifecycleKeys = Object.keys(lifecycle);
	}
	applyPlugin(plugin, instance) {
		require_logger.assert(require_tool.isPlainObject(plugin), "Plugin configuration is invalid.");
		const pluginName = plugin.name;
		require_logger.assert(pluginName, "A name must be provided by the plugin.");
		if (!this.registerPlugins[pluginName]) {
			this.registerPlugins[pluginName] = plugin;
			plugin.apply?.(instance);
			Object.keys(this.lifecycle).forEach((key) => {
				const pluginLife = plugin[key];
				if (pluginLife) this.lifecycle[key].on(pluginLife);
			});
		}
	}
	removePlugin(pluginName) {
		require_logger.assert(pluginName, "A name is required.");
		const plugin = this.registerPlugins[pluginName];
		require_logger.assert(plugin, `The plugin "${pluginName}" is not registered.`);
		Object.keys(plugin).forEach((key) => {
			if (key !== "name") this.lifecycle[key].remove(plugin[key]);
		});
	}
};

//#endregion
exports.PluginSystem = PluginSystem;
//# sourceMappingURL=pluginSystem.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs"(__unused_rspack_module, exports) {

//#region src/utils/hooks/syncHook.ts
var SyncHook = class {
	constructor(type) {
		this.type = "";
		this.listeners = /* @__PURE__ */ new Set();
		if (type) this.type = type;
	}
	on(fn) {
		if (typeof fn === "function") this.listeners.add(fn);
	}
	once(fn) {
		const self = this;
		this.on(function wrapper(...args) {
			self.remove(wrapper);
			return fn.apply(null, args);
		});
	}
	emit(...data) {
		let result;
		if (this.listeners.size > 0) this.listeners.forEach((fn) => {
			result = fn(...data);
		});
		return result;
	}
	remove(fn) {
		this.listeners.delete(fn);
	}
	removeAll() {
		this.listeners.clear();
	}
};

//#endregion
exports.SyncHook = SyncHook;
//# sourceMappingURL=syncHook.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncWaterfallHook.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_syncHook = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/hooks/syncHook.cjs");

//#region src/utils/hooks/syncWaterfallHook.ts
function checkReturnData(originalData, returnedData) {
	if (!require_tool.isObject(returnedData)) return false;
	if (originalData !== returnedData) {
		for (const key in originalData) if (!(key in returnedData)) return false;
	}
	return true;
}
var SyncWaterfallHook = class extends require_syncHook.SyncHook {
	constructor(type) {
		super();
		this.onerror = require_logger.error;
		this.type = type;
	}
	emit(data) {
		if (!require_tool.isObject(data)) require_logger.error(`The data for the "${this.type}" hook should be an object.`);
		for (const fn of this.listeners) try {
			const tempData = fn(data);
			if (checkReturnData(data, tempData)) data = tempData;
			else {
				this.onerror(`A plugin returned an unacceptable value for the "${this.type}" type.`);
				break;
			}
		} catch (e) {
			require_logger.warn(e);
			this.onerror(e);
		}
		return data;
	}
};

//#endregion
exports.SyncWaterfallHook = SyncWaterfallHook;
exports.checkReturnData = checkReturnData;
//# sourceMappingURL=syncWaterfallHook.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/index.cjs"(__unused_rspack_module, __unused_rspack_exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_env = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/env.cjs");
const require_manifest = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs");
const require_plugin = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/plugin.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
const require_context = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/context.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");


},
"../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/constant.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/utils/load.ts
const importCallback = ".then(callbacks[0]).catch(callbacks[1])";
async function loadEsmEntry({ entry, remoteEntryExports }) {
	return new Promise((resolve, reject) => {
		try {
			if (!remoteEntryExports) if (typeof FEDERATION_ALLOW_NEW_FUNCTION !== "undefined") new Function("callbacks", `import("${entry}")${importCallback}`)([resolve, reject]);
			else import(
				/* webpackIgnore: true */
				/* @vite-ignore */
				entry
).then(resolve).catch(reject);
			else resolve(remoteEntryExports);
		} catch (e) {
			require_logger.error(`Failed to load ESM entry from "${entry}". ${e instanceof Error ? e.message : String(e)}`);
		}
	});
}
async function loadSystemJsEntry({ entry, remoteEntryExports }) {
	return new Promise((resolve, reject) => {
		try {
			if (!remoteEntryExports) if (false) {}
			else new Function("callbacks", `System.import("${entry}")${importCallback}`)([resolve, reject]);
			else resolve(remoteEntryExports);
		} catch (e) {
			require_logger.error(`Failed to load SystemJS entry from "${entry}". ${e instanceof Error ? e.message : String(e)}`);
		}
	});
}
function handleRemoteEntryLoaded(name, globalName, entry) {
	const { remoteEntryKey, entryExports } = require_global.getRemoteEntryExports(name, globalName);
	if (!entryExports) require_logger.error(_module_federation_error_codes.RUNTIME_001, _module_federation_error_codes.runtimeDescMap, {
		remoteName: name,
		remoteEntryUrl: entry,
		remoteEntryKey
	});
	return entryExports;
}
async function loadEntryScript({ name, globalName, entry, remoteInfo, loaderHook, getEntryUrl }) {
	const { entryExports: remoteEntryExports } = require_global.getRemoteEntryExports(name, globalName);
	if (remoteEntryExports) return remoteEntryExports;
	const url = getEntryUrl ? getEntryUrl(entry) : entry;
	return (0, _module_federation_sdk.loadScript)(url, {
		attrs: {},
		createScriptHook: (url, attrs) => {
			const res = loaderHook.lifecycle.createScript.emit({
				url,
				attrs,
				remoteInfo
			});
			if (!res) return;
			if (res instanceof HTMLScriptElement) return res;
			if ("script" in res || "timeout" in res) return res;
		}
	}).then(() => {
		return handleRemoteEntryLoaded(name, globalName, entry);
	}, (loadError) => {
		const originalMsg = loadError instanceof Error ? loadError.message : String(loadError);
		require_logger.error(_module_federation_error_codes.RUNTIME_008, _module_federation_error_codes.runtimeDescMap, {
			remoteName: name,
			resourceUrl: url
		}, originalMsg);
	});
}
async function loadEntryDom({ remoteInfo, remoteEntryExports, loaderHook, getEntryUrl }) {
	const { entry, entryGlobalName: globalName, name, type } = remoteInfo;
	switch (type) {
		case "esm":
		case "module": return loadEsmEntry({
			entry,
			remoteEntryExports
		});
		case "system": return loadSystemJsEntry({
			entry,
			remoteEntryExports
		});
		default: return loadEntryScript({
			entry,
			globalName,
			name,
			remoteInfo,
			loaderHook,
			getEntryUrl
		});
	}
}
async function loadEntryNode({ remoteInfo, loaderHook }) {
	const { entry, entryGlobalName: globalName, name, type } = remoteInfo;
	const { entryExports: remoteEntryExports } = require_global.getRemoteEntryExports(name, globalName);
	if (remoteEntryExports) return remoteEntryExports;
	return (0, _module_federation_sdk.loadScriptNode)(entry, {
		attrs: {
			name,
			globalName,
			type
		},
		loaderHook: { createScriptHook: (url, attrs = {}) => {
			const res = loaderHook.lifecycle.createScript.emit({
				url,
				attrs,
				remoteInfo
			});
			if (!res) return;
			if ("url" in res) return res;
		} }
	}).then(() => {
		return handleRemoteEntryLoaded(name, globalName, entry);
	}).catch((e) => {
		require_logger.error(`Failed to load Node.js entry for remote "${name}" from "${entry}". ${e instanceof Error ? e.message : String(e)}`);
	});
}
function getRemoteEntryUniqueKey(remoteInfo) {
	const { entry, name } = remoteInfo;
	return (0, _module_federation_sdk.composeKeyWithSeparator)(name, entry);
}
async function getRemoteEntry(params) {
	const { origin, remoteEntryExports, remoteInfo, getEntryUrl, _inErrorHandling = false } = params;
	const uniqueKey = getRemoteEntryUniqueKey(remoteInfo);
	if (remoteEntryExports) return remoteEntryExports;
	if (!require_global.globalLoading[uniqueKey]) {
		const loadEntryHook = origin.remoteHandler.hooks.lifecycle.loadEntry;
		const loaderHook = origin.loaderHook;
		require_global.globalLoading[uniqueKey] = loadEntryHook.emit({
			loaderHook,
			remoteInfo,
			remoteEntryExports
		}).then((res) => {
			if (res) return res;
			return (typeof ENV_TARGET !== "undefined" ? ENV_TARGET === "web" : _module_federation_sdk.isBrowserEnvValue) ? loadEntryDom({
				remoteInfo,
				remoteEntryExports,
				loaderHook,
				getEntryUrl
			}) : loadEntryNode({
				remoteInfo,
				loaderHook
			});
		}).catch(async (err) => {
			const uniqueKey = getRemoteEntryUniqueKey(remoteInfo);
			const isScriptExecutionError = err instanceof Error && err.message.includes("ScriptExecutionError");
			if (err instanceof Error && err.message.includes(_module_federation_error_codes.RUNTIME_008) && !isScriptExecutionError && !_inErrorHandling) {
				const wrappedGetRemoteEntry = (params) => {
					return getRemoteEntry({
						...params,
						_inErrorHandling: true
					});
				};
				const RemoteEntryExports = await origin.loaderHook.lifecycle.loadEntryError.emit({
					getRemoteEntry: wrappedGetRemoteEntry,
					origin,
					remoteInfo,
					remoteEntryExports,
					globalLoading: require_global.globalLoading,
					uniqueKey
				});
				if (RemoteEntryExports) return RemoteEntryExports;
			}
			throw err;
		});
	}
	return require_global.globalLoading[uniqueKey];
}
function getRemoteInfo(remote) {
	return {
		...remote,
		entry: "entry" in remote ? remote.entry : "",
		type: remote.type || require_constant.DEFAULT_REMOTE_TYPE,
		entryGlobalName: remote.entryGlobalName || remote.name,
		shareScope: remote.shareScope || require_constant.DEFAULT_SCOPE
	};
}

//#endregion
exports.getRemoteEntry = getRemoteEntry;
exports.getRemoteEntryUniqueKey = getRemoteEntryUniqueKey;
exports.getRemoteInfo = getRemoteInfo;
//# sourceMappingURL=load.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs"(__unused_rspack_module, exports, __webpack_require__) {
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");
let _module_federation_error_codes_browser = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/browser.cjs");

//#region src/utils/logger.ts
const LOG_CATEGORY = "[ Federation Runtime ]";
const logger = (0, _module_federation_sdk.createLogger)(LOG_CATEGORY);
function assert(condition, msgOrCode, descMap, args, context) {
	if (!condition) if (descMap !== void 0) error(msgOrCode, descMap, args, void 0, context);
	else error(msgOrCode);
}
function error(msgOrCode, descMap, args, originalErrorMsg, context) {
	if (descMap !== void 0) return (0, _module_federation_error_codes_browser.logAndReport)(msgOrCode, descMap, args ?? {}, (msg) => {
		throw new Error(`${LOG_CATEGORY}: ${msg}`);
	}, originalErrorMsg, context);
	const msg = msgOrCode;
	if (msg instanceof Error) {
		if (!msg.message.startsWith(LOG_CATEGORY)) msg.message = `${LOG_CATEGORY}: ${msg.message}`;
		throw msg;
	}
	throw new Error(`${LOG_CATEGORY}: ${msg}`);
}
function warn(msg) {
	if (msg instanceof Error) {
		if (!msg.message.startsWith(LOG_CATEGORY)) msg.message = `${LOG_CATEGORY}: ${msg.message}`;
		logger.warn(msg);
	} else logger.warn(msg);
}

//#endregion
exports.assert = assert;
exports.error = error;
exports.logger = logger;
exports.warn = warn;
//# sourceMappingURL=logger.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs"(__unused_rspack_module, exports) {

//#region src/utils/manifest.ts
function matchRemoteWithNameAndExpose(remotes, id) {
	for (const remote of remotes) {
		const isNameMatched = id.startsWith(remote.name);
		let expose = id.replace(remote.name, "");
		if (isNameMatched) {
			if (expose.startsWith("/")) {
				const pkgNameOrAlias = remote.name;
				expose = `.${expose}`;
				return {
					pkgNameOrAlias,
					expose,
					remote
				};
			} else if (expose === "") return {
				pkgNameOrAlias: remote.name,
				expose: ".",
				remote
			};
		}
		const isAliasMatched = remote.alias && id.startsWith(remote.alias);
		let exposeWithAlias = remote.alias && id.replace(remote.alias, "");
		if (remote.alias && isAliasMatched) {
			if (exposeWithAlias && exposeWithAlias.startsWith("/")) {
				const pkgNameOrAlias = remote.alias;
				exposeWithAlias = `.${exposeWithAlias}`;
				return {
					pkgNameOrAlias,
					expose: exposeWithAlias,
					remote
				};
			} else if (exposeWithAlias === "") return {
				pkgNameOrAlias: remote.alias,
				expose: ".",
				remote
			};
		}
	}
}
function matchRemote(remotes, nameOrAlias) {
	for (const remote of remotes) {
		if (nameOrAlias === remote.name) return remote;
		if (remote.alias && nameOrAlias === remote.alias) return remote;
	}
}

//#endregion
exports.matchRemote = matchRemote;
exports.matchRemoteWithNameAndExpose = matchRemoteWithNameAndExpose;
//# sourceMappingURL=manifest.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/plugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");

//#region src/utils/plugin.ts
function registerPlugins(plugins, instance) {
	const globalPlugins = require_global.getGlobalHostPlugins();
	const hookInstances = [
		instance.hooks,
		instance.remoteHandler.hooks,
		instance.sharedHandler.hooks,
		instance.snapshotHandler.hooks,
		instance.loaderHook,
		instance.bridgeHook
	];
	if (globalPlugins.length > 0) globalPlugins.forEach((plugin) => {
		if (plugins?.find((item) => item.name !== plugin.name)) plugins.push(plugin);
	});
	if (plugins && plugins.length > 0) plugins.forEach((plugin) => {
		hookInstances.forEach((hookInstance) => {
			hookInstance.applyPlugin(plugin, instance);
		});
	});
	return plugins;
}

//#endregion
exports.registerPlugins = registerPlugins;
//# sourceMappingURL=plugin.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/preload.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_manifest = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/manifest.cjs");
const require_load = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/load.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/utils/preload.ts
function defaultPreloadArgs(preloadConfig) {
	return {
		resourceCategory: "sync",
		share: true,
		depsRemote: true,
		...preloadConfig
	};
}
function formatPreloadArgs(remotes, preloadArgs) {
	return preloadArgs.map((args) => {
		const remoteInfo = require_manifest.matchRemote(remotes, args.nameOrAlias);
		require_logger.assert(remoteInfo, `Unable to preload ${args.nameOrAlias} as it is not included in ${!remoteInfo && (0, _module_federation_sdk.safeToString)({
			remoteInfo,
			remotes
		})}`);
		return {
			remote: remoteInfo,
			preloadConfig: defaultPreloadArgs(args)
		};
	});
}
function normalizePreloadExposes(exposes) {
	if (!exposes) return [];
	return exposes.map((expose) => {
		if (expose === ".") return expose;
		if (expose.startsWith("./")) return expose.replace("./", "");
		return expose;
	});
}
function preloadAssets(remoteInfo, host, assets, useLinkPreload = true) {
	const { cssAssets, jsAssetsWithoutEntry, entryAssets } = assets;
	if (host.options.inBrowser) {
		entryAssets.forEach((asset) => {
			const { moduleInfo } = asset;
			const module = host.moduleCache.get(remoteInfo.name);
			if (module) require_load.getRemoteEntry({
				origin: host,
				remoteInfo: moduleInfo,
				remoteEntryExports: module.remoteEntryExports
			});
			else require_load.getRemoteEntry({
				origin: host,
				remoteInfo: moduleInfo,
				remoteEntryExports: void 0
			});
		});
		if (useLinkPreload) {
			const defaultAttrs = {
				rel: "preload",
				as: "style"
			};
			cssAssets.forEach((cssUrl) => {
				const { link: cssEl, needAttach } = (0, _module_federation_sdk.createLink)({
					url: cssUrl,
					cb: () => {},
					attrs: defaultAttrs,
					createLinkHook: (url, attrs) => {
						const res = host.loaderHook.lifecycle.createLink.emit({
							url,
							attrs,
							remoteInfo
						});
						if (res instanceof HTMLLinkElement) return res;
					}
				});
				needAttach && document.head.appendChild(cssEl);
			});
		} else {
			const defaultAttrs = {
				rel: "stylesheet",
				type: "text/css"
			};
			cssAssets.forEach((cssUrl) => {
				const { link: cssEl, needAttach } = (0, _module_federation_sdk.createLink)({
					url: cssUrl,
					cb: () => {},
					attrs: defaultAttrs,
					createLinkHook: (url, attrs) => {
						const res = host.loaderHook.lifecycle.createLink.emit({
							url,
							attrs,
							remoteInfo
						});
						if (res instanceof HTMLLinkElement) return res;
					},
					needDeleteLink: false
				});
				needAttach && document.head.appendChild(cssEl);
			});
		}
		if (useLinkPreload) {
			const defaultAttrs = {
				rel: "preload",
				as: "script"
			};
			jsAssetsWithoutEntry.forEach((jsUrl) => {
				const { link: linkEl, needAttach } = (0, _module_federation_sdk.createLink)({
					url: jsUrl,
					cb: () => {},
					attrs: defaultAttrs,
					createLinkHook: (url, attrs) => {
						const res = host.loaderHook.lifecycle.createLink.emit({
							url,
							attrs,
							remoteInfo
						});
						if (res instanceof HTMLLinkElement) return res;
					}
				});
				needAttach && document.head.appendChild(linkEl);
			});
		} else {
			const defaultAttrs = {
				fetchpriority: "high",
				type: remoteInfo?.type === "module" ? "module" : "text/javascript"
			};
			jsAssetsWithoutEntry.forEach((jsUrl) => {
				const { script: scriptEl, needAttach } = (0, _module_federation_sdk.createScript)({
					url: jsUrl,
					cb: () => {},
					attrs: defaultAttrs,
					createScriptHook: (url, attrs) => {
						const res = host.loaderHook.lifecycle.createScript.emit({
							url,
							attrs,
							remoteInfo
						});
						if (res instanceof HTMLScriptElement) return res;
					},
					needDeleteScript: true
				});
				needAttach && document.head.appendChild(scriptEl);
			});
		}
	}
}

//#endregion
exports.defaultPreloadArgs = defaultPreloadArgs;
exports.formatPreloadArgs = formatPreloadArgs;
exports.normalizePreloadExposes = normalizePreloadExposes;
exports.preloadAssets = preloadAssets;
//# sourceMappingURL=preload.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/semver/compare.cjs"(__unused_rspack_module, exports) {

//#region src/utils/semver/compare.ts
function compareAtom(rangeAtom, versionAtom) {
	rangeAtom = Number(rangeAtom) || rangeAtom;
	versionAtom = Number(versionAtom) || versionAtom;
	if (rangeAtom > versionAtom) return 1;
	if (rangeAtom === versionAtom) return 0;
	return -1;
}
function comparePreRelease(rangeAtom, versionAtom) {
	const { preRelease: rangePreRelease } = rangeAtom;
	const { preRelease: versionPreRelease } = versionAtom;
	if (rangePreRelease === void 0 && Boolean(versionPreRelease)) return 1;
	if (Boolean(rangePreRelease) && versionPreRelease === void 0) return -1;
	if (rangePreRelease === void 0 && versionPreRelease === void 0) return 0;
	for (let i = 0, n = rangePreRelease.length; i <= n; i++) {
		const rangeElement = rangePreRelease[i];
		const versionElement = versionPreRelease[i];
		if (rangeElement === versionElement) continue;
		if (rangeElement === void 0 && versionElement === void 0) return 0;
		if (!rangeElement) return 1;
		if (!versionElement) return -1;
		return compareAtom(rangeElement, versionElement);
	}
	return 0;
}
function compareVersion(rangeAtom, versionAtom) {
	return compareAtom(rangeAtom.major, versionAtom.major) || compareAtom(rangeAtom.minor, versionAtom.minor) || compareAtom(rangeAtom.patch, versionAtom.patch) || comparePreRelease(rangeAtom, versionAtom);
}
function eq(rangeAtom, versionAtom) {
	return rangeAtom.version === versionAtom.version;
}
function compare(rangeAtom, versionAtom) {
	switch (rangeAtom.operator) {
		case "":
		case "=": return eq(rangeAtom, versionAtom);
		case ">": return compareVersion(rangeAtom, versionAtom) < 0;
		case ">=": return eq(rangeAtom, versionAtom) || compareVersion(rangeAtom, versionAtom) < 0;
		case "<": return compareVersion(rangeAtom, versionAtom) > 0;
		case "<=": return eq(rangeAtom, versionAtom) || compareVersion(rangeAtom, versionAtom) > 0;
		case void 0: return true;
		default: return false;
	}
}

//#endregion
exports.compare = compare;
//# sourceMappingURL=compare.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/semver/constants.cjs"(__unused_rspack_module, exports) {

//#region src/utils/semver/constants.ts
const buildIdentifier = "[0-9A-Za-z-]+";
const build = `(?:\\+(${buildIdentifier}(?:\\.${buildIdentifier})*))`;
const numericIdentifier = "0|[1-9]\\d*";
const numericIdentifierLoose = "[0-9]+";
const nonNumericIdentifier = "\\d*[a-zA-Z-][a-zA-Z0-9-]*";
const preReleaseIdentifierLoose = `(?:${numericIdentifierLoose}|${nonNumericIdentifier})`;
const preReleaseLoose = `(?:-?(${preReleaseIdentifierLoose}(?:\\.${preReleaseIdentifierLoose})*))`;
const preReleaseIdentifier = `(?:${numericIdentifier}|${nonNumericIdentifier})`;
const preRelease = `(?:-(${preReleaseIdentifier}(?:\\.${preReleaseIdentifier})*))`;
const xRangeIdentifier = `${numericIdentifier}|x|X|\\*`;
const xRangePlain = `[v=\\s]*(${xRangeIdentifier})(?:\\.(${xRangeIdentifier})(?:\\.(${xRangeIdentifier})(?:${preRelease})?${build}?)?)?`;
const hyphenRange = `^\\s*(${xRangePlain})\\s+-\\s+(${xRangePlain})\\s*$`;
const loosePlain = `[v=\\s]*${`(${numericIdentifierLoose})\\.(${numericIdentifierLoose})\\.(${numericIdentifierLoose})`}${preReleaseLoose}?${build}?`;
const gtlt = "((?:<|>)?=?)";
const comparatorTrim = `(\\s*)${gtlt}\\s*(${loosePlain}|${xRangePlain})`;
const loneTilde = "(?:~>?)";
const tildeTrim = `(\\s*)${loneTilde}\\s+`;
const loneCaret = "(?:\\^)";
const caretTrim = `(\\s*)${loneCaret}\\s+`;
const star = "(<|>)?=?\\s*\\*";
const caret = `^${loneCaret}${xRangePlain}$`;
const fullPlain = `v?${`(${numericIdentifier})\\.(${numericIdentifier})\\.(${numericIdentifier})`}${preRelease}?${build}?`;
const tilde = `^${loneTilde}${xRangePlain}$`;
const xRange = `^${gtlt}\\s*${xRangePlain}$`;
const comparator = `^${gtlt}\\s*(${fullPlain})$|^$`;
const gte0 = "^\\s*>=\\s*0.0.0\\s*$";

//#endregion
exports.caret = caret;
exports.caretTrim = caretTrim;
exports.comparator = comparator;
exports.comparatorTrim = comparatorTrim;
exports.gte0 = gte0;
exports.hyphenRange = hyphenRange;
exports.star = star;
exports.tilde = tilde;
exports.tildeTrim = tildeTrim;
exports.xRange = xRange;
//# sourceMappingURL=constants.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/semver/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_utils = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/utils.cjs");
const require_parser = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/parser.cjs");
const require_compare = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/compare.cjs");

//#region src/utils/semver/index.ts
function parseComparatorString(range) {
	return require_utils.pipe(require_parser.parseCarets, require_parser.parseTildes, require_parser.parseXRanges, require_parser.parseStar)(range);
}
function parseRange(range) {
	return require_utils.pipe(require_parser.parseHyphen, require_parser.parseComparatorTrim, require_parser.parseTildeTrim, require_parser.parseCaretTrim)(range.trim()).split(/\s+/).join(" ");
}
function satisfy(version, range) {
	if (!version) return false;
	const extractedVersion = require_utils.extractComparator(version);
	if (!extractedVersion) return false;
	const [, versionOperator, , versionMajor, versionMinor, versionPatch, versionPreRelease] = extractedVersion;
	const versionAtom = {
		operator: versionOperator,
		version: require_utils.combineVersion(versionMajor, versionMinor, versionPatch, versionPreRelease),
		major: versionMajor,
		minor: versionMinor,
		patch: versionPatch,
		preRelease: versionPreRelease?.split(".")
	};
	const orRanges = range.split("||");
	for (const orRange of orRanges) {
		const trimmedOrRange = orRange.trim();
		if (!trimmedOrRange) return true;
		if (trimmedOrRange === "*" || trimmedOrRange === "x") return true;
		try {
			const parsedSubRange = parseRange(trimmedOrRange);
			if (!parsedSubRange.trim()) return true;
			const parsedComparatorString = parsedSubRange.split(" ").map((rangeVersion) => parseComparatorString(rangeVersion)).join(" ");
			if (!parsedComparatorString.trim()) return true;
			const comparators = parsedComparatorString.split(/\s+/).map((comparator) => require_parser.parseGTE0(comparator)).filter(Boolean);
			if (comparators.length === 0) continue;
			let subRangeSatisfied = true;
			for (const comparator of comparators) {
				const extractedComparator = require_utils.extractComparator(comparator);
				if (!extractedComparator) {
					subRangeSatisfied = false;
					break;
				}
				const [, rangeOperator, , rangeMajor, rangeMinor, rangePatch, rangePreRelease] = extractedComparator;
				if (!require_compare.compare({
					operator: rangeOperator,
					version: require_utils.combineVersion(rangeMajor, rangeMinor, rangePatch, rangePreRelease),
					major: rangeMajor,
					minor: rangeMinor,
					patch: rangePatch,
					preRelease: rangePreRelease?.split(".")
				}, versionAtom)) {
					subRangeSatisfied = false;
					break;
				}
			}
			if (subRangeSatisfied) return true;
		} catch (e) {
			console.error(`[semver] Error processing range part "${trimmedOrRange}":`, e);
			continue;
		}
	}
	return false;
}

//#endregion
exports.satisfy = satisfy;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/semver/parser.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_constants = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/constants.cjs");
const require_utils = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/utils.cjs");

//#region src/utils/semver/parser.ts
function parseHyphen(range) {
	return range.replace(require_utils.parseRegex(require_constants.hyphenRange), (_range, from, fromMajor, fromMinor, fromPatch, _fromPreRelease, _fromBuild, to, toMajor, toMinor, toPatch, toPreRelease) => {
		if (require_utils.isXVersion(fromMajor)) from = "";
		else if (require_utils.isXVersion(fromMinor)) from = `>=${fromMajor}.0.0`;
		else if (require_utils.isXVersion(fromPatch)) from = `>=${fromMajor}.${fromMinor}.0`;
		else from = `>=${from}`;
		if (require_utils.isXVersion(toMajor)) to = "";
		else if (require_utils.isXVersion(toMinor)) to = `<${Number(toMajor) + 1}.0.0-0`;
		else if (require_utils.isXVersion(toPatch)) to = `<${toMajor}.${Number(toMinor) + 1}.0-0`;
		else if (toPreRelease) to = `<=${toMajor}.${toMinor}.${toPatch}-${toPreRelease}`;
		else to = `<=${to}`;
		return `${from} ${to}`.trim();
	});
}
function parseComparatorTrim(range) {
	return range.replace(require_utils.parseRegex(require_constants.comparatorTrim), "$1$2$3");
}
function parseTildeTrim(range) {
	return range.replace(require_utils.parseRegex(require_constants.tildeTrim), "$1~");
}
function parseCaretTrim(range) {
	return range.replace(require_utils.parseRegex(require_constants.caretTrim), "$1^");
}
function parseCarets(range) {
	return range.trim().split(/\s+/).map((rangeVersion) => rangeVersion.replace(require_utils.parseRegex(require_constants.caret), (_, major, minor, patch, preRelease) => {
		if (require_utils.isXVersion(major)) return "";
		else if (require_utils.isXVersion(minor)) return `>=${major}.0.0 <${Number(major) + 1}.0.0-0`;
		else if (require_utils.isXVersion(patch)) if (major === "0") return `>=${major}.${minor}.0 <${major}.${Number(minor) + 1}.0-0`;
		else return `>=${major}.${minor}.0 <${Number(major) + 1}.0.0-0`;
		else if (preRelease) if (major === "0") if (minor === "0") return `>=${major}.${minor}.${patch}-${preRelease} <${major}.${minor}.${Number(patch) + 1}-0`;
		else return `>=${major}.${minor}.${patch}-${preRelease} <${major}.${Number(minor) + 1}.0-0`;
		else return `>=${major}.${minor}.${patch}-${preRelease} <${Number(major) + 1}.0.0-0`;
		else {
			if (major === "0") if (minor === "0") return `>=${major}.${minor}.${patch} <${major}.${minor}.${Number(patch) + 1}-0`;
			else return `>=${major}.${minor}.${patch} <${major}.${Number(minor) + 1}.0-0`;
			return `>=${major}.${minor}.${patch} <${Number(major) + 1}.0.0-0`;
		}
	})).join(" ");
}
function parseTildes(range) {
	return range.trim().split(/\s+/).map((rangeVersion) => rangeVersion.replace(require_utils.parseRegex(require_constants.tilde), (_, major, minor, patch, preRelease) => {
		if (require_utils.isXVersion(major)) return "";
		else if (require_utils.isXVersion(minor)) return `>=${major}.0.0 <${Number(major) + 1}.0.0-0`;
		else if (require_utils.isXVersion(patch)) return `>=${major}.${minor}.0 <${major}.${Number(minor) + 1}.0-0`;
		else if (preRelease) return `>=${major}.${minor}.${patch}-${preRelease} <${major}.${Number(minor) + 1}.0-0`;
		return `>=${major}.${minor}.${patch} <${major}.${Number(minor) + 1}.0-0`;
	})).join(" ");
}
function parseXRanges(range) {
	return range.split(/\s+/).map((rangeVersion) => rangeVersion.trim().replace(require_utils.parseRegex(require_constants.xRange), (ret, gtlt, major, minor, patch, preRelease) => {
		const isXMajor = require_utils.isXVersion(major);
		const isXMinor = isXMajor || require_utils.isXVersion(minor);
		const isXPatch = isXMinor || require_utils.isXVersion(patch);
		if (gtlt === "=" && isXPatch) gtlt = "";
		preRelease = "";
		if (isXMajor) if (gtlt === ">" || gtlt === "<") return "<0.0.0-0";
		else return "*";
		else if (gtlt && isXPatch) {
			if (isXMinor) minor = 0;
			patch = 0;
			if (gtlt === ">") {
				gtlt = ">=";
				if (isXMinor) {
					major = Number(major) + 1;
					minor = 0;
					patch = 0;
				} else {
					minor = Number(minor) + 1;
					patch = 0;
				}
			} else if (gtlt === "<=") {
				gtlt = "<";
				if (isXMinor) major = Number(major) + 1;
				else minor = Number(minor) + 1;
			}
			if (gtlt === "<") preRelease = "-0";
			return `${gtlt + major}.${minor}.${patch}${preRelease}`;
		} else if (isXMinor) return `>=${major}.0.0${preRelease} <${Number(major) + 1}.0.0-0`;
		else if (isXPatch) return `>=${major}.${minor}.0${preRelease} <${major}.${Number(minor) + 1}.0-0`;
		return ret;
	})).join(" ");
}
function parseStar(range) {
	return range.trim().replace(require_utils.parseRegex(require_constants.star), "");
}
function parseGTE0(comparatorString) {
	return comparatorString.trim().replace(require_utils.parseRegex(require_constants.gte0), "");
}

//#endregion
exports.parseCaretTrim = parseCaretTrim;
exports.parseCarets = parseCarets;
exports.parseComparatorTrim = parseComparatorTrim;
exports.parseGTE0 = parseGTE0;
exports.parseHyphen = parseHyphen;
exports.parseStar = parseStar;
exports.parseTildeTrim = parseTildeTrim;
exports.parseTildes = parseTildes;
exports.parseXRanges = parseXRanges;
//# sourceMappingURL=parser.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/semver/utils.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_constants = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/constants.cjs");

//#region src/utils/semver/utils.ts
function parseRegex(source) {
	return new RegExp(source);
}
function isXVersion(version) {
	return !version || version.toLowerCase() === "x" || version === "*";
}
function pipe(...fns) {
	return (x) => fns.reduce((v, f) => f(v), x);
}
function extractComparator(comparatorString) {
	return comparatorString.match(parseRegex(require_constants.comparator));
}
function combineVersion(major, minor, patch, preRelease) {
	const mainVersion = `${major}.${minor}.${patch}`;
	if (preRelease) return `${mainVersion}-${preRelease}`;
	return mainVersion;
}

//#endregion
exports.combineVersion = combineVersion;
exports.extractComparator = extractComparator;
exports.isXVersion = isXVersion;
exports.parseRegex = parseRegex;
exports.pipe = pipe;
//# sourceMappingURL=utils.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/share.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
const require_tool = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs");
const require_global = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/global.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/constant.cjs");
const require_index = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/semver/index.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/utils/share.ts
function formatShare(shareArgs, from, name, shareStrategy) {
	let get;
	if ("get" in shareArgs) get = shareArgs.get;
	else if ("lib" in shareArgs) get = () => Promise.resolve(shareArgs.lib);
	else get = () => Promise.resolve(() => {
		require_logger.error(`Cannot get shared "${name}" from "${from}": neither "get" nor "lib" is provided in the share config.`);
	});
	if (shareArgs.shareConfig?.eager && shareArgs.treeShaking?.mode) require_logger.error(`Invalid shared config for "${name}" from "${from}": cannot use both "eager: true" and "treeShaking.mode" simultaneously. Choose one strategy.`);
	return {
		deps: [],
		useIn: [],
		from,
		loading: null,
		...shareArgs,
		shareConfig: {
			requiredVersion: `^${shareArgs.version}`,
			singleton: false,
			eager: false,
			strictVersion: false,
			...shareArgs.shareConfig
		},
		get,
		loaded: shareArgs?.loaded || "lib" in shareArgs ? true : void 0,
		version: shareArgs.version ?? "0",
		scope: Array.isArray(shareArgs.scope) ? shareArgs.scope : [shareArgs.scope ?? "default"],
		strategy: (shareArgs.strategy ?? shareStrategy) || "version-first",
		treeShaking: shareArgs.treeShaking ? {
			...shareArgs.treeShaking,
			mode: shareArgs.treeShaking.mode ?? "server-calc",
			status: shareArgs.treeShaking.status ?? _module_federation_sdk.TreeShakingStatus.UNKNOWN,
			useIn: []
		} : void 0
	};
}
function formatShareConfigs(prevOptions, newOptions) {
	const shareArgs = newOptions.shared || {};
	const from = newOptions.name;
	const newShareInfos = Object.keys(shareArgs).reduce((res, pkgName) => {
		const arrayShareArgs = require_tool.arrayOptions(shareArgs[pkgName]);
		res[pkgName] = res[pkgName] || [];
		arrayShareArgs.forEach((shareConfig) => {
			res[pkgName].push(formatShare(shareConfig, from, pkgName, newOptions.shareStrategy));
		});
		return res;
	}, {});
	const allShareInfos = { ...prevOptions.shared };
	Object.keys(newShareInfos).forEach((shareKey) => {
		if (!allShareInfos[shareKey]) allShareInfos[shareKey] = newShareInfos[shareKey];
		else newShareInfos[shareKey].forEach((newUserSharedOptions) => {
			if (!allShareInfos[shareKey].find((sharedVal) => sharedVal.version === newUserSharedOptions.version)) allShareInfos[shareKey].push(newUserSharedOptions);
		});
	});
	return {
		allShareInfos,
		newShareInfos
	};
}
function shouldUseTreeShaking(treeShaking, usedExports) {
	if (!treeShaking) return false;
	const { status, mode } = treeShaking;
	if (status === _module_federation_sdk.TreeShakingStatus.NO_USE) return false;
	if (status === _module_federation_sdk.TreeShakingStatus.CALCULATED) return true;
	if (mode === "runtime-infer") {
		if (!usedExports) return true;
		return isMatchUsedExports(treeShaking, usedExports);
	}
	return false;
}
/**
* compare version a and b, return true if a is less than b
*/
function versionLt(a, b) {
	const transformInvalidVersion = (version) => {
		if (!Number.isNaN(Number(version))) {
			const splitArr = version.split(".");
			let validVersion = version;
			for (let i = 0; i < 3 - splitArr.length; i++) validVersion += ".0";
			return validVersion;
		}
		return version;
	};
	if (require_index.satisfy(transformInvalidVersion(a), `<=${transformInvalidVersion(b)}`)) return true;
	else return false;
}
const findVersion = (shareVersionMap, cb) => {
	const callback = cb || function(prev, cur) {
		return versionLt(prev, cur);
	};
	return Object.keys(shareVersionMap).reduce((prev, cur) => {
		if (!prev) return cur;
		if (callback(prev, cur)) return cur;
		if (prev === "0") return cur;
		return prev;
	}, 0);
};
const isLoaded = (shared) => {
	return Boolean(shared.loaded) || typeof shared.lib === "function";
};
const isLoading = (shared) => {
	return Boolean(shared.loading);
};
const isMatchUsedExports = (treeShaking, usedExports) => {
	if (!treeShaking || !usedExports) return false;
	const { usedExports: treeShakingUsedExports } = treeShaking;
	if (!treeShakingUsedExports) return false;
	if (usedExports.every((e) => treeShakingUsedExports.includes(e))) return true;
	return false;
};
function findSingletonVersionOrderByVersion(shareScopeMap, scope, pkgName, treeShaking) {
	const versions = shareScopeMap[scope][pkgName];
	let version = "";
	let useTreesShaking = shouldUseTreeShaking(treeShaking);
	const callback = function(prev, cur) {
		if (useTreesShaking) {
			if (!versions[prev].treeShaking) return true;
			if (!versions[cur].treeShaking) return false;
			return !isLoaded(versions[prev].treeShaking) && versionLt(prev, cur);
		}
		return !isLoaded(versions[prev]) && versionLt(prev, cur);
	};
	if (useTreesShaking) {
		version = findVersion(shareScopeMap[scope][pkgName], callback);
		if (version) return {
			version,
			useTreesShaking
		};
		useTreesShaking = false;
	}
	return {
		version: findVersion(shareScopeMap[scope][pkgName], callback),
		useTreesShaking
	};
}
const isLoadingOrLoaded = (shared) => {
	return isLoaded(shared) || isLoading(shared);
};
function findSingletonVersionOrderByLoaded(shareScopeMap, scope, pkgName, treeShaking) {
	const versions = shareScopeMap[scope][pkgName];
	let version = "";
	let useTreesShaking = shouldUseTreeShaking(treeShaking);
	const callback = function(prev, cur) {
		if (useTreesShaking) {
			if (!versions[prev].treeShaking) return true;
			if (!versions[cur].treeShaking) return false;
			if (isLoadingOrLoaded(versions[cur].treeShaking)) if (isLoadingOrLoaded(versions[prev].treeShaking)) return Boolean(versionLt(prev, cur));
			else return true;
			if (isLoadingOrLoaded(versions[prev].treeShaking)) return false;
		}
		if (isLoadingOrLoaded(versions[cur])) if (isLoadingOrLoaded(versions[prev])) return Boolean(versionLt(prev, cur));
		else return true;
		if (isLoadingOrLoaded(versions[prev])) return false;
		return versionLt(prev, cur);
	};
	if (useTreesShaking) {
		version = findVersion(shareScopeMap[scope][pkgName], callback);
		if (version) return {
			version,
			useTreesShaking
		};
		useTreesShaking = false;
	}
	return {
		version: findVersion(shareScopeMap[scope][pkgName], callback),
		useTreesShaking
	};
}
function getFindShareFunction(strategy) {
	if (strategy === "loaded-first") return findSingletonVersionOrderByLoaded;
	return findSingletonVersionOrderByVersion;
}
function getRegisteredShare(localShareScopeMap, pkgName, shareInfo, resolveShare) {
	if (!localShareScopeMap) return;
	const { shareConfig, scope = require_constant.DEFAULT_SCOPE, strategy, treeShaking } = shareInfo;
	const scopes = Array.isArray(scope) ? scope : [scope];
	for (const sc of scopes) if (shareConfig && localShareScopeMap[sc] && localShareScopeMap[sc][pkgName]) {
		const { requiredVersion } = shareConfig;
		const { version: maxOrSingletonVersion, useTreesShaking } = getFindShareFunction(strategy)(localShareScopeMap, sc, pkgName, treeShaking);
		const defaultResolver = () => {
			const shared = localShareScopeMap[sc][pkgName][maxOrSingletonVersion];
			if (shareConfig.singleton) {
				if (typeof requiredVersion === "string" && !require_index.satisfy(maxOrSingletonVersion, requiredVersion)) {
					const msg = `Version ${maxOrSingletonVersion} from ${maxOrSingletonVersion && shared.from} of shared singleton module ${pkgName} does not satisfy the requirement of ${shareInfo.from} which needs ${requiredVersion})`;
					if (shareConfig.strictVersion) require_logger.error(msg);
					else require_logger.warn(msg);
				}
				return {
					shared,
					useTreesShaking
				};
			} else {
				if (requiredVersion === false || requiredVersion === "*") return {
					shared,
					useTreesShaking
				};
				if (require_index.satisfy(maxOrSingletonVersion, requiredVersion)) return {
					shared,
					useTreesShaking
				};
				const _usedTreeShaking = shouldUseTreeShaking(treeShaking);
				if (_usedTreeShaking) for (const [versionKey, versionValue] of Object.entries(localShareScopeMap[sc][pkgName])) {
					if (!shouldUseTreeShaking(versionValue.treeShaking, treeShaking?.usedExports)) continue;
					if (require_index.satisfy(versionKey, requiredVersion)) return {
						shared: versionValue,
						useTreesShaking: _usedTreeShaking
					};
				}
				for (const [versionKey, versionValue] of Object.entries(localShareScopeMap[sc][pkgName])) if (require_index.satisfy(versionKey, requiredVersion)) return {
					shared: versionValue,
					useTreesShaking: false
				};
			}
		};
		const params = {
			shareScopeMap: localShareScopeMap,
			scope: sc,
			pkgName,
			version: maxOrSingletonVersion,
			GlobalFederation: require_global.Global.__FEDERATION__,
			shareInfo,
			resolver: defaultResolver
		};
		return (resolveShare.emit(params) || params).resolver();
	}
}
function getGlobalShareScope() {
	return require_global.Global.__FEDERATION__.__SHARE__;
}
function getTargetSharedOptions(options) {
	const { pkgName, extraOptions, shareInfos } = options;
	const defaultResolver = (sharedOptions) => {
		if (!sharedOptions) return;
		const shareVersionMap = {};
		sharedOptions.forEach((shared) => {
			shareVersionMap[shared.version] = shared;
		});
		const callback = function(prev, cur) {
			return !isLoaded(shareVersionMap[prev]) && versionLt(prev, cur);
		};
		return shareVersionMap[findVersion(shareVersionMap, callback)];
	};
	const resolver = extraOptions?.resolver ?? defaultResolver;
	const isPlainObject = (val) => {
		return val !== null && typeof val === "object" && !Array.isArray(val);
	};
	const merge = (...sources) => {
		const out = {};
		for (const src of sources) {
			if (!src) continue;
			for (const [key, value] of Object.entries(src)) {
				const prev = out[key];
				if (isPlainObject(prev) && isPlainObject(value)) out[key] = merge(prev, value);
				else if (value !== void 0) out[key] = value;
			}
		}
		return out;
	};
	return merge(resolver(shareInfos[pkgName]), extraOptions?.customShareInfo);
}
const addUseIn = (shared, from) => {
	if (!shared.useIn) shared.useIn = [];
	require_tool.addUniqueItem(shared.useIn, from);
};
function directShare(shared, useTreesShaking) {
	if (useTreesShaking && shared.treeShaking) return shared.treeShaking;
	return shared;
}

//#endregion
exports.addUseIn = addUseIn;
exports.directShare = directShare;
exports.formatShareConfigs = formatShareConfigs;
exports.getGlobalShareScope = getGlobalShareScope;
exports.getRegisteredShare = getRegisteredShare;
exports.getTargetSharedOptions = getTargetSharedOptions;
exports.shouldUseTreeShaking = shouldUseTreeShaking;
//# sourceMappingURL=share.cjs.map

},
"../../node_modules/@module-federation/runtime-core/dist/utils/tool.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_logger = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/utils/logger.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/utils/tool.ts
function addUniqueItem(arr, item) {
	if (arr.findIndex((name) => name === item) === -1) arr.push(item);
	return arr;
}
function getFMId(remoteInfo) {
	if ("version" in remoteInfo && remoteInfo.version) return `${remoteInfo.name}:${remoteInfo.version}`;
	else if ("entry" in remoteInfo && remoteInfo.entry) return `${remoteInfo.name}:${remoteInfo.entry}`;
	else return `${remoteInfo.name}`;
}
function isRemoteInfoWithEntry(remote) {
	return typeof remote.entry !== "undefined";
}
function isPureRemoteEntry(remote) {
	return !remote.entry.includes(".json");
}
async function safeWrapper(callback, disableWarn) {
	try {
		return await callback();
	} catch (e) {
		!disableWarn && require_logger.warn(e);
		return;
	}
}
function isObject(val) {
	return val && typeof val === "object";
}
const objectToString = Object.prototype.toString;
function isPlainObject(val) {
	return objectToString.call(val) === "[object Object]";
}
function isStaticResourcesEqual(url1, url2) {
	const REG_EXP = /^(https?:)?\/\//i;
	return url1.replace(REG_EXP, "").replace(/\/$/, "") === url2.replace(REG_EXP, "").replace(/\/$/, "");
}
function arrayOptions(options) {
	return Array.isArray(options) ? options : [options];
}
function getRemoteEntryInfoFromSnapshot(snapshot) {
	const defaultRemoteEntryInfo = {
		url: "",
		type: "global",
		globalName: ""
	};
	if (_module_federation_sdk.isBrowserEnvValue || (0, _module_federation_sdk.isReactNativeEnv)() || !("ssrRemoteEntry" in snapshot)) return "remoteEntry" in snapshot ? {
		url: snapshot.remoteEntry,
		type: snapshot.remoteEntryType,
		globalName: snapshot.globalName
	} : defaultRemoteEntryInfo;
	if ("ssrRemoteEntry" in snapshot) return {
		url: snapshot.ssrRemoteEntry || defaultRemoteEntryInfo.url,
		type: snapshot.ssrRemoteEntryType || defaultRemoteEntryInfo.type,
		globalName: snapshot.globalName
	};
	return defaultRemoteEntryInfo;
}
const processModuleAlias = (name, subPath) => {
	let moduleName;
	if (name.endsWith("/")) moduleName = name.slice(0, -1);
	else moduleName = name;
	if (subPath.startsWith(".")) subPath = subPath.slice(1);
	moduleName = moduleName + subPath;
	return moduleName;
};

//#endregion
exports.addUniqueItem = addUniqueItem;
exports.arrayOptions = arrayOptions;
exports.getFMId = getFMId;
exports.getRemoteEntryInfoFromSnapshot = getRemoteEntryInfoFromSnapshot;
exports.isObject = isObject;
exports.isPlainObject = isPlainObject;
exports.isPureRemoteEntry = isPureRemoteEntry;
exports.isRemoteInfoWithEntry = isRemoteInfoWithEntry;
exports.isStaticResourcesEqual = isStaticResourcesEqual;
exports.objectToString = objectToString;
exports.processModuleAlias = processModuleAlias;
exports.safeWrapper = safeWrapper;
//# sourceMappingURL=tool.cjs.map

},
"../../node_modules/@module-federation/runtime/dist/_virtual/_rolldown/runtime.cjs"(__unused_rspack_module, exports) {
//#region \0rolldown/runtime.js
var __create = Object.create;
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __getProtoOf = Object.getPrototypeOf;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __copyProps = (to, from, except, desc) => {
	if (from && typeof from === "object" || typeof from === "function") {
		for (var keys = __getOwnPropNames(from), i = 0, n = keys.length, key; i < n; i++) {
			key = keys[i];
			if (!__hasOwnProp.call(to, key) && key !== except) {
				__defProp(to, key, {
					get: ((k) => from[k]).bind(null, key),
					enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable
				});
			}
		}
	}
	return to;
};
var __toESM = (mod, isNodeMode, target) => (target = mod != null ? __create(__getProtoOf(mod)) : {}, __copyProps(isNodeMode || !mod || !mod.__esModule ? __defProp(target, "default", {
	value: mod,
	enumerable: true
}) : target, mod));

//#endregion

exports.__toESM = __toESM;

},
"../../node_modules/@module-federation/runtime/dist/helpers.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperties(exports, { __esModule: { value: true }, [Symbol.toStringTag]: { value: 'Module' } });
const require_runtime = __webpack_require__("../../node_modules/@module-federation/runtime/dist/_virtual/_rolldown/runtime.cjs");
const require_utils = __webpack_require__("../../node_modules/@module-federation/runtime/dist/utils.cjs");
let _module_federation_runtime_core = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/index.cjs");

//#region src/helpers.ts
const global = {
	..._module_federation_runtime_core.helpers.global,
	getGlobalFederationInstance: require_utils.getGlobalFederationInstance
};
const share = _module_federation_runtime_core.helpers.share;
const utils = _module_federation_runtime_core.helpers.utils;
const runtimeHelpers = {
	global,
	share,
	utils
};

//#endregion
exports["default"] = runtimeHelpers;
exports.global = global;
exports.share = share;
exports.utils = utils;
//# sourceMappingURL=helpers.cjs.map

},
"../../node_modules/@module-federation/runtime/dist/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
const require_runtime = __webpack_require__("../../node_modules/@module-federation/runtime/dist/_virtual/_rolldown/runtime.cjs");
const require_utils = __webpack_require__("../../node_modules/@module-federation/runtime/dist/utils.cjs");
let _module_federation_runtime_core = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/index.cjs");
let _module_federation_error_codes = __webpack_require__("../../node_modules/@module-federation/error-codes/dist/index.cjs");

//#region src/index.ts
function createInstance(options) {
	const instance = new (((0, _module_federation_runtime_core.getGlobalFederationConstructor)()) || _module_federation_runtime_core.ModuleFederation)({
		id: `${options.name}@${options.version || Date.now()}`,
		...options
	});
	(0, _module_federation_runtime_core.setGlobalFederationInstance)(instance);
	return instance;
}
let FederationInstance = null;
/**
* @deprecated Use createInstance or getInstance instead
*/
function init(options) {
	const instance = require_utils.getGlobalFederationInstance(options.name, options.version);
	const normalizedOptions = {
		...options,
		id: options.id || ""
	};
	if (!instance) {
		FederationInstance = createInstance(normalizedOptions);
		return FederationInstance;
	} else {
		instance.initOptions(normalizedOptions);
		if (!FederationInstance) FederationInstance = instance;
		return instance;
	}
}
function loadRemote(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.loadRemote.apply(FederationInstance, args);
}
function loadShare(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.loadShare.apply(FederationInstance, args);
}
function loadShareSync(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.loadShareSync.apply(FederationInstance, args);
}
function preloadRemote(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.preloadRemote.apply(FederationInstance, args);
}
function registerRemotes(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.registerRemotes.apply(FederationInstance, args);
}
function registerPlugins(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.registerPlugins.apply(FederationInstance, args);
}
function getInstance() {
	return FederationInstance;
}
function registerShared(...args) {
	(0, _module_federation_runtime_core.assert)(FederationInstance, _module_federation_error_codes.RUNTIME_009, _module_federation_error_codes.runtimeDescMap);
	return FederationInstance.registerShared.apply(FederationInstance, args);
}
(0, _module_federation_runtime_core.setGlobalFederationConstructor)(_module_federation_runtime_core.ModuleFederation);

//#endregion
Object.defineProperty(exports, "Module", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.Module;
  }
}));
exports.ModuleFederation = _module_federation_runtime_core.ModuleFederation;
exports.createInstance = createInstance;
exports.getInstance = getInstance;
Object.defineProperty(exports, "getRemoteEntry", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.getRemoteEntry;
  }
}));
Object.defineProperty(exports, "getRemoteInfo", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.getRemoteInfo;
  }
}));
exports.init = init;
exports.loadRemote = loadRemote;
Object.defineProperty(exports, "loadScript", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.loadScript;
  }
}));
Object.defineProperty(exports, "loadScriptNode", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.loadScriptNode;
  }
}));
exports.loadShare = loadShare;
exports.loadShareSync = loadShareSync;
exports.preloadRemote = preloadRemote;
Object.defineProperty(exports, "registerGlobalPlugins", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime_core.registerGlobalPlugins;
  }
}));
exports.registerPlugins = registerPlugins;
exports.registerRemotes = registerRemotes;
exports.registerShared = registerShared;
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/runtime/dist/utils.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/runtime/dist/_virtual/_rolldown/runtime.cjs");
let _module_federation_runtime_core = __webpack_require__("../../node_modules/@module-federation/runtime-core/dist/index.cjs");

//#region src/utils.ts
function getBuilderId() {
	return typeof FEDERATION_BUILD_IDENTIFIER !== "undefined" ? FEDERATION_BUILD_IDENTIFIER : "";
}
function getGlobalFederationInstance(name, version) {
	const buildId = getBuilderId();
	return _module_federation_runtime_core.CurrentGlobal.__FEDERATION__.__INSTANCES__.find((GMInstance) => {
		if (buildId && GMInstance.options.id === buildId) return true;
		if (GMInstance.options.name === name && !GMInstance.options.version && !version) return true;
		if (GMInstance.options.name === name && version && GMInstance.options.version === version) return true;
		return false;
	});
}

//#endregion
exports.getGlobalFederationInstance = getGlobalFederationInstance;
//# sourceMappingURL=utils.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs"(__unused_rspack_module, exports) {
//#region \0rolldown/runtime.js
var __defProp = Object.defineProperty;
var __exportAll = (all, no_symbols) => {
	let target = {};
	for (var name in all) {
		__defProp(target, name, {
			get: all[name],
			enumerable: true
		});
	}
	if (!no_symbols) {
		__defProp(target, Symbol.toStringTag, { value: "Module" });
	}
	return target;
};

//#endregion

exports.__exportAll = __exportAll;

},
"../../node_modules/@module-federation/sdk/dist/constant.cjs"(__unused_rspack_module, exports) {

//#region src/constant.ts
const FederationModuleManifest = "federation-manifest.json";
const MANIFEST_EXT = ".json";
const BROWSER_LOG_KEY = "FEDERATION_DEBUG";
const NameTransformSymbol = {
	AT: "@",
	HYPHEN: "-",
	SLASH: "/"
};
const NameTransformMap = {
	[NameTransformSymbol.AT]: "scope_",
	[NameTransformSymbol.HYPHEN]: "_",
	[NameTransformSymbol.SLASH]: "__"
};
const EncodedNameTransformMap = {
	[NameTransformMap[NameTransformSymbol.AT]]: NameTransformSymbol.AT,
	[NameTransformMap[NameTransformSymbol.HYPHEN]]: NameTransformSymbol.HYPHEN,
	[NameTransformMap[NameTransformSymbol.SLASH]]: NameTransformSymbol.SLASH
};
const SEPARATOR = ":";
const ManifestFileName = "mf-manifest.json";
const StatsFileName = "mf-stats.json";
const MFModuleType = {
	NPM: "npm",
	APP: "app"
};
const MODULE_DEVTOOL_IDENTIFIER = "__MF_DEVTOOLS_MODULE_INFO__";
const ENCODE_NAME_PREFIX = "ENCODE_NAME_PREFIX";
const TEMP_DIR = ".federation";
let TreeShakingStatus = /* @__PURE__ */ function(TreeShakingStatus) {
	/**
	* Not handled by deploy server, needs to infer by the real runtime period.
	*/
	TreeShakingStatus[TreeShakingStatus["UNKNOWN"] = 1] = "UNKNOWN";
	/**
	* It means the shared has been calculated , runtime should take this shared as first choice.
	*/
	TreeShakingStatus[TreeShakingStatus["CALCULATED"] = 2] = "CALCULATED";
	/**
	* It means the shared has been calculated, and marked as no used
	*/
	TreeShakingStatus[TreeShakingStatus["NO_USE"] = 0] = "NO_USE";
	return TreeShakingStatus;
}({});

//#endregion
exports.BROWSER_LOG_KEY = BROWSER_LOG_KEY;
exports.ENCODE_NAME_PREFIX = ENCODE_NAME_PREFIX;
exports.EncodedNameTransformMap = EncodedNameTransformMap;
exports.FederationModuleManifest = FederationModuleManifest;
exports.MANIFEST_EXT = MANIFEST_EXT;
exports.MFModuleType = MFModuleType;
exports.MODULE_DEVTOOL_IDENTIFIER = MODULE_DEVTOOL_IDENTIFIER;
exports.ManifestFileName = ManifestFileName;
exports.NameTransformMap = NameTransformMap;
exports.NameTransformSymbol = NameTransformSymbol;
exports.SEPARATOR = SEPARATOR;
exports.StatsFileName = StatsFileName;
exports.TEMP_DIR = TEMP_DIR;
exports.TreeShakingStatus = TreeShakingStatus;
//# sourceMappingURL=constant.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/createModuleFederationConfig.cjs"(__unused_rspack_module, exports) {

//#region src/createModuleFederationConfig.ts
const createModuleFederationConfig = (options) => {
	return options;
};

//#endregion
exports.createModuleFederationConfig = createModuleFederationConfig;
//# sourceMappingURL=createModuleFederationConfig.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/dom.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_utils = __webpack_require__("../../node_modules/@module-federation/sdk/dist/utils.cjs");

//#region src/dom.ts
async function safeWrapper(callback, disableWarn) {
	try {
		return await callback();
	} catch (e) {
		!disableWarn && require_utils.warn(e);
		return;
	}
}
function isStaticResourcesEqual(url1, url2) {
	const REG_EXP = /^(https?:)?\/\//i;
	return url1.replace(REG_EXP, "").replace(/\/$/, "") === url2.replace(REG_EXP, "").replace(/\/$/, "");
}
function createScript(info) {
	let script = null;
	let needAttach = true;
	let timeout = 2e4;
	let timeoutId;
	const scripts = document.getElementsByTagName("script");
	for (let i = 0; i < scripts.length; i++) {
		const s = scripts[i];
		const scriptSrc = s.getAttribute("src");
		if (scriptSrc && isStaticResourcesEqual(scriptSrc, info.url)) {
			script = s;
			needAttach = false;
			break;
		}
	}
	if (!script) {
		const attrs = info.attrs;
		script = document.createElement("script");
		script.type = attrs?.["type"] === "module" ? "module" : "text/javascript";
		let createScriptRes = void 0;
		if (info.createScriptHook) {
			createScriptRes = info.createScriptHook(info.url, info.attrs);
			if (createScriptRes instanceof HTMLScriptElement) script = createScriptRes;
			else if (typeof createScriptRes === "object") {
				if ("script" in createScriptRes && createScriptRes.script) script = createScriptRes.script;
				if ("timeout" in createScriptRes && createScriptRes.timeout) timeout = createScriptRes.timeout;
			}
		}
		if (!script.src) script.src = info.url;
		if (attrs && !createScriptRes) Object.keys(attrs).forEach((name) => {
			if (script) {
				if (name === "async" || name === "defer") script[name] = attrs[name];
				else if (!script.getAttribute(name)) script.setAttribute(name, attrs[name]);
			}
		});
	}
	let executionError = null;
	const executionErrorHandler = typeof window !== "undefined" ? (evt) => {
		if (evt.filename && isStaticResourcesEqual(evt.filename, info.url)) {
			const err = /* @__PURE__ */ new Error(`ScriptExecutionError: Script "${info.url}" loaded but threw a runtime error during execution: ${evt.message} (${evt.filename}:${evt.lineno}:${evt.colno})`);
			err.name = "ScriptExecutionError";
			executionError = err;
		}
	} : null;
	if (executionErrorHandler) window.addEventListener("error", executionErrorHandler);
	const onScriptComplete = async (prev, event) => {
		clearTimeout(timeoutId);
		if (executionErrorHandler) window.removeEventListener("error", executionErrorHandler);
		const onScriptCompleteCallback = () => {
			if (event?.type === "error") {
				const networkError = /* @__PURE__ */ new Error(event?.isTimeout ? `ScriptNetworkError: Script "${info.url}" timed out.` : `ScriptNetworkError: Failed to load script "${info.url}" - the script URL is unreachable or the server returned an error (network failure, 404, CORS, etc.)`);
				networkError.name = "ScriptNetworkError";
				info?.onErrorCallback && info?.onErrorCallback(networkError);
			} else if (executionError) info?.onErrorCallback && info?.onErrorCallback(executionError);
			else info?.cb && info?.cb();
		};
		if (script) {
			script.onerror = null;
			script.onload = null;
			safeWrapper(() => {
				const { needDeleteScript = true } = info;
				if (needDeleteScript) script?.parentNode && script.parentNode.removeChild(script);
			});
			if (prev && typeof prev === "function") {
				const result = prev(event);
				if (result instanceof Promise) {
					const res = await result;
					onScriptCompleteCallback();
					return res;
				}
				onScriptCompleteCallback();
				return result;
			}
		}
		onScriptCompleteCallback();
	};
	script.onerror = onScriptComplete.bind(null, script.onerror);
	script.onload = onScriptComplete.bind(null, script.onload);
	timeoutId = setTimeout(() => {
		onScriptComplete(null, {
			type: "error",
			isTimeout: true
		});
	}, timeout);
	return {
		script,
		needAttach
	};
}
function createLink(info) {
	let link = null;
	let needAttach = true;
	const links = document.getElementsByTagName("link");
	for (let i = 0; i < links.length; i++) {
		const l = links[i];
		const linkHref = l.getAttribute("href");
		const linkRel = l.getAttribute("rel");
		if (linkHref && isStaticResourcesEqual(linkHref, info.url) && linkRel === info.attrs["rel"]) {
			link = l;
			needAttach = false;
			break;
		}
	}
	if (!link) {
		link = document.createElement("link");
		link.setAttribute("href", info.url);
		let createLinkRes = void 0;
		const attrs = info.attrs;
		if (info.createLinkHook) {
			createLinkRes = info.createLinkHook(info.url, attrs);
			if (createLinkRes instanceof HTMLLinkElement) link = createLinkRes;
		}
		if (attrs && !createLinkRes) Object.keys(attrs).forEach((name) => {
			if (link && !link.getAttribute(name)) link.setAttribute(name, attrs[name]);
		});
	}
	const onLinkComplete = (prev, event) => {
		const onLinkCompleteCallback = () => {
			if (event?.type === "error") info?.onErrorCallback && info?.onErrorCallback(event);
			else info?.cb && info?.cb();
		};
		if (link) {
			link.onerror = null;
			link.onload = null;
			safeWrapper(() => {
				const { needDeleteLink = true } = info;
				if (needDeleteLink) link?.parentNode && link.parentNode.removeChild(link);
			});
			if (prev) {
				const res = prev(event);
				onLinkCompleteCallback();
				return res;
			}
		}
		onLinkCompleteCallback();
	};
	link.onerror = onLinkComplete.bind(null, link.onerror);
	link.onload = onLinkComplete.bind(null, link.onload);
	return {
		link,
		needAttach
	};
}
function loadScript(url, info) {
	const { attrs = {}, createScriptHook } = info;
	return new Promise((resolve, reject) => {
		const { script, needAttach } = createScript({
			url,
			cb: resolve,
			onErrorCallback: reject,
			attrs: {
				fetchpriority: "high",
				...attrs
			},
			createScriptHook,
			needDeleteScript: true
		});
		needAttach && document.head.appendChild(script);
	});
}

//#endregion
exports.createLink = createLink;
exports.createScript = createScript;
exports.isStaticResourcesEqual = isStaticResourcesEqual;
exports.loadScript = loadScript;
exports.safeWrapper = safeWrapper;
//# sourceMappingURL=dom.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/env.cjs"(__unused_rspack_module, exports, __webpack_require__) {
/* provided dependency */ var process = __webpack_require__("../../node_modules/process/browser.js");
const require_constant = __webpack_require__("../../node_modules/@module-federation/sdk/dist/constant.cjs");

//#region src/env.ts
const isBrowserEnvValue = typeof ENV_TARGET !== "undefined" ? ENV_TARGET === "web" : typeof window !== "undefined" && typeof window.document !== "undefined";
function isBrowserEnv() {
	return isBrowserEnvValue;
}
function isReactNativeEnv() {
	return typeof navigator !== "undefined" && navigator?.product === "ReactNative";
}
function isBrowserDebug() {
	try {
		if (isBrowserEnv() && window.localStorage) return Boolean(localStorage.getItem(require_constant.BROWSER_LOG_KEY));
	} catch (error) {
		return false;
	}
	return false;
}
function isDebugMode() {
	if (typeof process !== "undefined" && process.env && process.env["FEDERATION_DEBUG"]) return Boolean(process.env["FEDERATION_DEBUG"]);
	if (typeof FEDERATION_DEBUG !== "undefined" && Boolean(FEDERATION_DEBUG)) return true;
	return isBrowserDebug();
}
const getProcessEnv = function() {
	return typeof process !== "undefined" && process.env ? process.env : {};
};

//#endregion
exports.getProcessEnv = getProcessEnv;
exports.isBrowserEnv = isBrowserEnv;
exports.isBrowserEnvValue = isBrowserEnvValue;
exports.isDebugMode = isDebugMode;
exports.isReactNativeEnv = isReactNativeEnv;
//# sourceMappingURL=env.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/generateSnapshotFromManifest.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_constant = __webpack_require__("../../node_modules/@module-federation/sdk/dist/constant.cjs");

//#region src/generateSnapshotFromManifest.ts
const simpleJoinRemoteEntry = (rPath, rName) => {
	if (!rPath) return rName;
	const transformPath = (str) => {
		if (str === ".") return "";
		if (str.startsWith("./")) return str.replace("./", "");
		if (str.startsWith("/")) {
			const strWithoutSlash = str.slice(1);
			if (strWithoutSlash.endsWith("/")) return strWithoutSlash.slice(0, -1);
			return strWithoutSlash;
		}
		return str;
	};
	const transformedPath = transformPath(rPath);
	if (!transformedPath) return rName;
	if (transformedPath.endsWith("/")) return `${transformedPath}${rName}`;
	return `${transformedPath}/${rName}`;
};
function inferAutoPublicPath(url) {
	return url.replace(/#.*$/, "").replace(/\?.*$/, "").replace(/\/[^\/]+$/, "/");
}
function generateSnapshotFromManifest(manifest, options = {}) {
	const { remotes = {}, overrides = {}, version } = options;
	let remoteSnapshot;
	const getPublicPath = () => {
		if ("publicPath" in manifest.metaData) {
			if ((manifest.metaData.publicPath === "auto" || manifest.metaData.publicPath === "") && version) return inferAutoPublicPath(version);
			return manifest.metaData.publicPath;
		} else return manifest.metaData.getPublicPath;
	};
	const overridesKeys = Object.keys(overrides);
	let remotesInfo = {};
	if (!Object.keys(remotes).length) remotesInfo = manifest.remotes?.reduce((res, next) => {
		let matchedVersion;
		const name = next.federationContainerName;
		if (overridesKeys.includes(name)) matchedVersion = overrides[name];
		else if ("version" in next) matchedVersion = next.version;
		else matchedVersion = next.entry;
		res[name] = { matchedVersion };
		return res;
	}, {}) || {};
	Object.keys(remotes).forEach((key) => remotesInfo[key] = { matchedVersion: overridesKeys.includes(key) ? overrides[key] : remotes[key] });
	const { remoteEntry: { path: remoteEntryPath, name: remoteEntryName, type: remoteEntryType }, types: remoteTypes = {
		path: "",
		name: "",
		zip: "",
		api: ""
	}, buildInfo: { buildVersion }, globalName, ssrRemoteEntry } = manifest.metaData;
	const { exposes } = manifest;
	let basicRemoteSnapshot = {
		version: version ? version : "",
		buildVersion,
		globalName,
		remoteEntry: simpleJoinRemoteEntry(remoteEntryPath, remoteEntryName),
		remoteEntryType,
		remoteTypes: simpleJoinRemoteEntry(remoteTypes.path, remoteTypes.name),
		remoteTypesZip: remoteTypes.zip || "",
		remoteTypesAPI: remoteTypes.api || "",
		remotesInfo,
		shared: manifest?.shared.map((item) => ({
			assets: item.assets,
			sharedName: item.name,
			version: item.version,
			usedExports: item.referenceExports || []
		})),
		modules: exposes?.map((expose) => ({
			moduleName: expose.name,
			modulePath: expose.path,
			assets: expose.assets
		}))
	};
	if ("publicPath" in manifest.metaData) {
		remoteSnapshot = {
			...basicRemoteSnapshot,
			publicPath: getPublicPath()
		};
		if (typeof manifest.metaData.ssrPublicPath === "string") remoteSnapshot.ssrPublicPath = manifest.metaData.ssrPublicPath;
	} else remoteSnapshot = {
		...basicRemoteSnapshot,
		getPublicPath: getPublicPath()
	};
	if (ssrRemoteEntry) {
		const fullSSRRemoteEntry = simpleJoinRemoteEntry(ssrRemoteEntry.path, ssrRemoteEntry.name);
		remoteSnapshot.ssrRemoteEntry = fullSSRRemoteEntry;
		remoteSnapshot.ssrRemoteEntryType = ssrRemoteEntry.type || "commonjs-module";
	}
	return remoteSnapshot;
}
function isManifestProvider(moduleInfo) {
	if ("remoteEntry" in moduleInfo && moduleInfo.remoteEntry.includes(require_constant.MANIFEST_EXT)) return true;
	else return false;
}
function getManifestFileName(manifestOptions) {
	if (!manifestOptions) return {
		statsFileName: require_constant.StatsFileName,
		manifestFileName: require_constant.ManifestFileName
	};
	let filePath = typeof manifestOptions === "boolean" ? "" : manifestOptions.filePath || "";
	let fileName = typeof manifestOptions === "boolean" ? "" : manifestOptions.fileName || "";
	const JSON_EXT = ".json";
	const addExt = (name) => {
		if (name.endsWith(JSON_EXT)) return name;
		return `${name}${JSON_EXT}`;
	};
	const insertSuffix = (name, suffix) => {
		return name.replace(JSON_EXT, `${suffix}${JSON_EXT}`);
	};
	const manifestFileName = fileName ? addExt(fileName) : require_constant.ManifestFileName;
	return {
		statsFileName: simpleJoinRemoteEntry(filePath, fileName ? insertSuffix(manifestFileName, "-stats") : require_constant.StatsFileName),
		manifestFileName: simpleJoinRemoteEntry(filePath, manifestFileName)
	};
}

//#endregion
exports.generateSnapshotFromManifest = generateSnapshotFromManifest;
exports.getManifestFileName = getManifestFileName;
exports.inferAutoPublicPath = inferAutoPublicPath;
exports.isManifestProvider = isManifestProvider;
exports.simpleJoinRemoteEntry = simpleJoinRemoteEntry;
//# sourceMappingURL=generateSnapshotFromManifest.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
const require_constant = __webpack_require__("../../node_modules/@module-federation/sdk/dist/constant.cjs");
const require_ContainerPlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/ContainerPlugin.cjs");
const require_ContainerReferencePlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/ContainerReferencePlugin.cjs");
const require_ModuleFederationPlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/ModuleFederationPlugin.cjs");
const require_SharePlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/SharePlugin.cjs");
const require_ConsumeSharedPlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/ConsumeSharedPlugin.cjs");
const require_ProvideSharedPlugin = __webpack_require__("../../node_modules/@module-federation/sdk/dist/types/plugins/ProvideSharedPlugin.cjs");
const require_env = __webpack_require__("../../node_modules/@module-federation/sdk/dist/env.cjs");
const require_utils = __webpack_require__("../../node_modules/@module-federation/sdk/dist/utils.cjs");
const require_generateSnapshotFromManifest = __webpack_require__("../../node_modules/@module-federation/sdk/dist/generateSnapshotFromManifest.cjs");
const require_logger = __webpack_require__("../../node_modules/@module-federation/sdk/dist/logger.cjs");
const require_dom = __webpack_require__("../../node_modules/@module-federation/sdk/dist/dom.cjs");
const require_node = __webpack_require__("../../node_modules/@module-federation/sdk/dist/node.cjs");
const require_normalizeOptions = __webpack_require__("../../node_modules/@module-federation/sdk/dist/normalizeOptions.cjs");
const require_createModuleFederationConfig = __webpack_require__("../../node_modules/@module-federation/sdk/dist/createModuleFederationConfig.cjs");

exports.BROWSER_LOG_KEY = require_constant.BROWSER_LOG_KEY;
exports.ENCODE_NAME_PREFIX = require_constant.ENCODE_NAME_PREFIX;
exports.EncodedNameTransformMap = require_constant.EncodedNameTransformMap;
exports.FederationModuleManifest = require_constant.FederationModuleManifest;
exports.MANIFEST_EXT = require_constant.MANIFEST_EXT;
exports.MFModuleType = require_constant.MFModuleType;
exports.MODULE_DEVTOOL_IDENTIFIER = require_constant.MODULE_DEVTOOL_IDENTIFIER;
exports.ManifestFileName = require_constant.ManifestFileName;
exports.NameTransformMap = require_constant.NameTransformMap;
exports.NameTransformSymbol = require_constant.NameTransformSymbol;
exports.SEPARATOR = require_constant.SEPARATOR;
exports.StatsFileName = require_constant.StatsFileName;
exports.TEMP_DIR = require_constant.TEMP_DIR;
exports.TreeShakingStatus = require_constant.TreeShakingStatus;
exports.assert = require_utils.assert;
exports.bindLoggerToCompiler = require_logger.bindLoggerToCompiler;
exports.composeKeyWithSeparator = require_utils.composeKeyWithSeparator;
Object.defineProperty(exports, "consumeSharedPlugin", ({
  enumerable: true,
  get: function () {
    return require_ConsumeSharedPlugin.ConsumeSharedPlugin_exports;
  }
}));
Object.defineProperty(exports, "containerPlugin", ({
  enumerable: true,
  get: function () {
    return require_ContainerPlugin.ContainerPlugin_exports;
  }
}));
Object.defineProperty(exports, "containerReferencePlugin", ({
  enumerable: true,
  get: function () {
    return require_ContainerReferencePlugin.ContainerReferencePlugin_exports;
  }
}));
exports.createInfrastructureLogger = require_logger.createInfrastructureLogger;
exports.createLink = require_dom.createLink;
exports.createLogger = require_logger.createLogger;
exports.createModuleFederationConfig = require_createModuleFederationConfig.createModuleFederationConfig;
exports.createScript = require_dom.createScript;
exports.createScriptNode = require_node.createScriptNode;
exports.decodeName = require_utils.decodeName;
exports.encodeName = require_utils.encodeName;
exports.error = require_utils.error;
exports.generateExposeFilename = require_utils.generateExposeFilename;
exports.generateShareFilename = require_utils.generateShareFilename;
exports.generateSnapshotFromManifest = require_generateSnapshotFromManifest.generateSnapshotFromManifest;
exports.getManifestFileName = require_generateSnapshotFromManifest.getManifestFileName;
exports.getProcessEnv = require_env.getProcessEnv;
exports.getResourceUrl = require_utils.getResourceUrl;
exports.inferAutoPublicPath = require_generateSnapshotFromManifest.inferAutoPublicPath;
exports.infrastructureLogger = require_logger.infrastructureLogger;
exports.isBrowserEnv = require_env.isBrowserEnv;
exports.isBrowserEnvValue = require_env.isBrowserEnvValue;
exports.isDebugMode = require_env.isDebugMode;
exports.isManifestProvider = require_generateSnapshotFromManifest.isManifestProvider;
exports.isReactNativeEnv = require_env.isReactNativeEnv;
exports.isRequiredVersion = require_utils.isRequiredVersion;
exports.isStaticResourcesEqual = require_dom.isStaticResourcesEqual;
exports.loadScript = require_dom.loadScript;
exports.loadScriptNode = require_node.loadScriptNode;
exports.logger = require_logger.logger;
Object.defineProperty(exports, "moduleFederationPlugin", ({
  enumerable: true,
  get: function () {
    return require_ModuleFederationPlugin.ModuleFederationPlugin_exports;
  }
}));
exports.normalizeOptions = require_normalizeOptions.normalizeOptions;
exports.parseEntry = require_utils.parseEntry;
Object.defineProperty(exports, "provideSharedPlugin", ({
  enumerable: true,
  get: function () {
    return require_ProvideSharedPlugin.ProvideSharedPlugin_exports;
  }
}));
exports.safeToString = require_utils.safeToString;
exports.safeWrapper = require_dom.safeWrapper;
Object.defineProperty(exports, "sharePlugin", ({
  enumerable: true,
  get: function () {
    return require_SharePlugin.SharePlugin_exports;
  }
}));
exports.simpleJoinRemoteEntry = require_generateSnapshotFromManifest.simpleJoinRemoteEntry;
exports.warn = require_utils.warn;

},
"../../node_modules/@module-federation/sdk/dist/logger.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_env = __webpack_require__("../../node_modules/@module-federation/sdk/dist/env.cjs");

//#region src/logger.ts
const PREFIX = "[ Module Federation ]";
const DEFAULT_DELEGATE = console;
const LOGGER_STACK_SKIP_TOKENS = [
	"logger.ts",
	"logger.js",
	"captureStackTrace",
	"Logger.emit",
	"Logger.log",
	"Logger.info",
	"Logger.warn",
	"Logger.error",
	"Logger.debug"
];
function captureStackTrace() {
	try {
		const stack = (/* @__PURE__ */ new Error()).stack;
		if (!stack) return;
		const [, ...rawLines] = stack.split("\n");
		const filtered = rawLines.filter((line) => !LOGGER_STACK_SKIP_TOKENS.some((token) => line.includes(token)));
		if (!filtered.length) return;
		return `Stack trace:\n${filtered.slice(0, 5).join("\n")}`;
	} catch {
		return;
	}
}
var Logger = class {
	constructor(prefix, delegate = DEFAULT_DELEGATE) {
		this.prefix = prefix;
		this.delegate = delegate ?? DEFAULT_DELEGATE;
	}
	setPrefix(prefix) {
		this.prefix = prefix;
	}
	setDelegate(delegate) {
		this.delegate = delegate ?? DEFAULT_DELEGATE;
	}
	emit(method, args) {
		const delegate = this.delegate;
		const stackTrace = require_env.isDebugMode() ? captureStackTrace() : void 0;
		const enrichedArgs = stackTrace ? [...args, stackTrace] : args;
		const order = (() => {
			switch (method) {
				case "log": return ["log", "info"];
				case "info": return ["info", "log"];
				case "warn": return [
					"warn",
					"info",
					"log"
				];
				case "error": return [
					"error",
					"warn",
					"log"
				];
				default: return ["debug", "log"];
			}
		})();
		for (const candidate of order) {
			const handler = delegate[candidate];
			if (typeof handler === "function") {
				handler.call(delegate, this.prefix, ...enrichedArgs);
				return;
			}
		}
		for (const candidate of order) {
			const handler = DEFAULT_DELEGATE[candidate];
			if (typeof handler === "function") {
				handler.call(DEFAULT_DELEGATE, this.prefix, ...enrichedArgs);
				return;
			}
		}
	}
	log(...args) {
		this.emit("log", args);
	}
	warn(...args) {
		this.emit("warn", args);
	}
	error(...args) {
		this.emit("error", args);
	}
	success(...args) {
		this.emit("info", args);
	}
	info(...args) {
		this.emit("info", args);
	}
	ready(...args) {
		this.emit("info", args);
	}
	debug(...args) {
		if (require_env.isDebugMode()) this.emit("debug", args);
	}
};
function createLogger(prefix) {
	return new Logger(prefix);
}
function createInfrastructureLogger(prefix) {
	const infrastructureLogger = new Logger(prefix);
	Object.defineProperty(infrastructureLogger, "__mf_infrastructure_logger__", {
		value: true,
		enumerable: false,
		configurable: false
	});
	return infrastructureLogger;
}
function bindLoggerToCompiler(loggerInstance, compiler, name) {
	if (!loggerInstance.__mf_infrastructure_logger__) return;
	if (!compiler?.getInfrastructureLogger) return;
	try {
		const infrastructureLogger = compiler.getInfrastructureLogger(name);
		if (infrastructureLogger && typeof infrastructureLogger === "object" && (typeof infrastructureLogger.log === "function" || typeof infrastructureLogger.info === "function" || typeof infrastructureLogger.warn === "function" || typeof infrastructureLogger.error === "function")) loggerInstance.setDelegate(infrastructureLogger);
	} catch {
		loggerInstance.setDelegate(void 0);
	}
}
const logger = createLogger(PREFIX);
const infrastructureLogger = createInfrastructureLogger(PREFIX);

//#endregion
exports.bindLoggerToCompiler = bindLoggerToCompiler;
exports.createInfrastructureLogger = createInfrastructureLogger;
exports.createLogger = createLogger;
exports.infrastructureLogger = infrastructureLogger;
exports.logger = logger;
//# sourceMappingURL=logger.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/node.cjs"(__unused_rspack_module, exports) {

//#region src/node.ts
const sdkImportCache = /* @__PURE__ */ new Map();
function importNodeModule(name) {
	if (!name) throw new Error("import specifier is required");
	if (sdkImportCache.has(name)) return sdkImportCache.get(name);
	const promise = new Function("name", `return import(name)`)(name).then((res) => res).catch((error) => {
		console.error(`Error importing module ${name}:`, error);
		sdkImportCache.delete(name);
		throw error;
	});
	sdkImportCache.set(name, promise);
	return promise;
}
const loadNodeFetch = async () => {
	const fetchModule = await importNodeModule("node-fetch");
	return fetchModule.default || fetchModule;
};
const lazyLoaderHookFetch = async (input, init, loaderHook) => {
	const hook = (url, init) => {
		return loaderHook.lifecycle.fetch.emit(url, init);
	};
	const res = await hook(input, init || {});
	if (!res || !(res instanceof Response)) return (typeof fetch === "undefined" ? await loadNodeFetch() : fetch)(input, init || {});
	return res;
};
const createScriptNode = typeof ENV_TARGET === "undefined" || ENV_TARGET !== "web" ? (url, cb, attrs, loaderHook) => {
	if (loaderHook?.createScriptHook) {
		const hookResult = loaderHook.createScriptHook(url);
		if (hookResult && typeof hookResult === "object" && "url" in hookResult) url = hookResult.url;
	}
	let urlObj;
	try {
		urlObj = new URL(url);
	} catch (e) {
		console.error("Error constructing URL:", e);
		cb(/* @__PURE__ */ new Error(`Invalid URL: ${e}`));
		return;
	}
	const getFetch = async () => {
		if (loaderHook?.fetch) return (input, init) => lazyLoaderHookFetch(input, init, loaderHook);
		return typeof fetch === "undefined" ? loadNodeFetch() : fetch;
	};
	const handleScriptFetch = async (f, urlObj) => {
		try {
			const res = await f(urlObj.href);
			const data = await res.text();
			const [path, vm] = await Promise.all([importNodeModule("path"), importNodeModule("vm")]);
			const scriptContext = {
				exports: {},
				module: { exports: {} }
			};
			const urlDirname = urlObj.pathname.split("/").slice(0, -1).join("/");
			const filename = path.basename(urlObj.pathname);
			const script = new vm.Script(`(function(exports, module, require, __dirname, __filename) {${data}\n})`, {
				filename,
				importModuleDynamically: vm.constants?.USE_MAIN_CONTEXT_DEFAULT_LOADER ?? importNodeModule
			});
			let requireFn;
			requireFn = eval("require");
			script.runInThisContext()(scriptContext.exports, scriptContext.module, requireFn, urlDirname, filename);
			const exportedInterface = scriptContext.module.exports || scriptContext.exports;
			if (attrs && exportedInterface && attrs["globalName"]) {
				cb(void 0, exportedInterface[attrs["globalName"]] || exportedInterface);
				return;
			}
			cb(void 0, exportedInterface);
		} catch (e) {
			cb(e instanceof Error ? e : /* @__PURE__ */ new Error(`Script execution error: ${e}`));
		}
	};
	getFetch().then(async (f) => {
		if (attrs?.["type"] === "esm" || attrs?.["type"] === "module") return loadModule(urlObj.href, {
			fetch: f,
			vm: await importNodeModule("vm")
		}).then(async (module) => {
			await module.evaluate();
			cb(void 0, module.namespace);
		}).catch((e) => {
			cb(e instanceof Error ? e : /* @__PURE__ */ new Error(`Script execution error: ${e}`));
		});
		handleScriptFetch(f, urlObj);
	}).catch((err) => {
		cb(err);
	});
} : (url, cb, attrs, loaderHook) => {
	cb(/* @__PURE__ */ new Error("createScriptNode is disabled in non-Node.js environment"));
};
const loadScriptNode = typeof ENV_TARGET === "undefined" || ENV_TARGET !== "web" ? (url, info) => {
	return new Promise((resolve, reject) => {
		createScriptNode(url, (error, scriptContext) => {
			if (error) reject(error);
			else {
				const remoteEntryKey = info?.attrs?.["globalName"] || `__FEDERATION_${info?.attrs?.["name"]}:custom__`;
				resolve(globalThis[remoteEntryKey] = scriptContext);
			}
		}, info.attrs, info.loaderHook);
	});
} : (url, info) => {
	throw new Error("loadScriptNode is disabled in non-Node.js environment");
};
const esmModuleCache = /* @__PURE__ */ new Map();
async function loadModule(url, options) {
	if (esmModuleCache.has(url)) return esmModuleCache.get(url);
	const { fetch, vm } = options;
	const code = await (await fetch(url)).text();
	const module = new vm.SourceTextModule(code, { importModuleDynamically: async (specifier, script) => {
		const resolvedUrl = new URL(specifier, url).href;
		return loadModule(resolvedUrl, options);
	} });
	esmModuleCache.set(url, module);
	await module.link(async (specifier) => {
		const resolvedUrl = new URL(specifier, url).href;
		return await loadModule(resolvedUrl, options);
	});
	return module;
}

//#endregion
exports.createScriptNode = createScriptNode;
exports.loadScriptNode = loadScriptNode;
//# sourceMappingURL=node.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/normalizeOptions.cjs"(__unused_rspack_module, exports) {

//#region src/normalizeOptions.ts
function normalizeOptions(enableDefault, defaultOptions, key) {
	return function(options) {
		if (options === false) return false;
		if (typeof options === "undefined") if (enableDefault) return defaultOptions;
		else return false;
		if (options === true) return defaultOptions;
		if (options && typeof options === "object") return {
			...defaultOptions,
			...options
		};
		throw new Error(`Unexpected type for \`${key}\`, expect boolean/undefined/object, got: ${typeof options}`);
	};
}

//#endregion
exports.normalizeOptions = normalizeOptions;
//# sourceMappingURL=normalizeOptions.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/ConsumeSharedPlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/ConsumeSharedPlugin.ts
var ConsumeSharedPlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "ConsumeSharedPlugin_exports", ({
  enumerable: true,
  get: function () {
    return ConsumeSharedPlugin_exports;
  }
}));
//# sourceMappingURL=ConsumeSharedPlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/ContainerPlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/ContainerPlugin.ts
var ContainerPlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "ContainerPlugin_exports", ({
  enumerable: true,
  get: function () {
    return ContainerPlugin_exports;
  }
}));
//# sourceMappingURL=ContainerPlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/ContainerReferencePlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/ContainerReferencePlugin.ts
var ContainerReferencePlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "ContainerReferencePlugin_exports", ({
  enumerable: true,
  get: function () {
    return ContainerReferencePlugin_exports;
  }
}));
//# sourceMappingURL=ContainerReferencePlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/ModuleFederationPlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/ModuleFederationPlugin.ts
var ModuleFederationPlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "ModuleFederationPlugin_exports", ({
  enumerable: true,
  get: function () {
    return ModuleFederationPlugin_exports;
  }
}));
//# sourceMappingURL=ModuleFederationPlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/ProvideSharedPlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/ProvideSharedPlugin.ts
var ProvideSharedPlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "ProvideSharedPlugin_exports", ({
  enumerable: true,
  get: function () {
    return ProvideSharedPlugin_exports;
  }
}));
//# sourceMappingURL=ProvideSharedPlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/types/plugins/SharePlugin.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/sdk/dist/_virtual/_rolldown/runtime.cjs");

//#region src/types/plugins/SharePlugin.ts
var SharePlugin_exports = /* @__PURE__ */ require_runtime.__exportAll({});

//#endregion
Object.defineProperty(exports, "SharePlugin_exports", ({
  enumerable: true,
  get: function () {
    return SharePlugin_exports;
  }
}));
//# sourceMappingURL=SharePlugin.cjs.map

},
"../../node_modules/@module-federation/sdk/dist/utils.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_constant = __webpack_require__("../../node_modules/@module-federation/sdk/dist/constant.cjs");
const require_env = __webpack_require__("../../node_modules/@module-federation/sdk/dist/env.cjs");

//#region src/utils.ts
const LOG_CATEGORY = "[ Federation Runtime ]";
const parseEntry = (str, devVerOrUrl, separator = require_constant.SEPARATOR) => {
	const strSplit = str.split(separator);
	const devVersionOrUrl = require_env.getProcessEnv()["NODE_ENV"] === "development" && devVerOrUrl;
	const defaultVersion = "*";
	const isEntry = (s) => s.startsWith("http") || s.includes(require_constant.MANIFEST_EXT);
	if (strSplit.length >= 2) {
		let [name, ...versionOrEntryArr] = strSplit;
		if (str.startsWith(separator)) {
			name = strSplit.slice(0, 2).join(separator);
			versionOrEntryArr = [devVersionOrUrl || strSplit.slice(2).join(separator)];
		}
		let versionOrEntry = devVersionOrUrl || versionOrEntryArr.join(separator);
		if (isEntry(versionOrEntry)) return {
			name,
			entry: versionOrEntry
		};
		else return {
			name,
			version: versionOrEntry || defaultVersion
		};
	} else if (strSplit.length === 1) {
		const [name] = strSplit;
		if (devVersionOrUrl && isEntry(devVersionOrUrl)) return {
			name,
			entry: devVersionOrUrl
		};
		return {
			name,
			version: devVersionOrUrl || defaultVersion
		};
	} else throw `Invalid entry value: ${str}`;
};
const composeKeyWithSeparator = function(...args) {
	if (!args.length) return "";
	return args.reduce((sum, cur) => {
		if (!cur) return sum;
		if (!sum) return cur;
		return `${sum}${require_constant.SEPARATOR}${cur}`;
	}, "");
};
const encodeName = function(name, prefix = "", withExt = false) {
	try {
		const ext = withExt ? ".js" : "";
		return `${prefix}${name.replace(new RegExp(`${require_constant.NameTransformSymbol.AT}`, "g"), require_constant.NameTransformMap[require_constant.NameTransformSymbol.AT]).replace(new RegExp(`${require_constant.NameTransformSymbol.HYPHEN}`, "g"), require_constant.NameTransformMap[require_constant.NameTransformSymbol.HYPHEN]).replace(new RegExp(`${require_constant.NameTransformSymbol.SLASH}`, "g"), require_constant.NameTransformMap[require_constant.NameTransformSymbol.SLASH])}${ext}`;
	} catch (err) {
		throw err;
	}
};
const decodeName = function(name, prefix, withExt) {
	try {
		let decodedName = name;
		if (prefix) {
			if (!decodedName.startsWith(prefix)) return decodedName;
			decodedName = decodedName.replace(new RegExp(prefix, "g"), "");
		}
		decodedName = decodedName.replace(new RegExp(`${require_constant.NameTransformMap[require_constant.NameTransformSymbol.AT]}`, "g"), require_constant.EncodedNameTransformMap[require_constant.NameTransformMap[require_constant.NameTransformSymbol.AT]]).replace(new RegExp(`${require_constant.NameTransformMap[require_constant.NameTransformSymbol.SLASH]}`, "g"), require_constant.EncodedNameTransformMap[require_constant.NameTransformMap[require_constant.NameTransformSymbol.SLASH]]).replace(new RegExp(`${require_constant.NameTransformMap[require_constant.NameTransformSymbol.HYPHEN]}`, "g"), require_constant.EncodedNameTransformMap[require_constant.NameTransformMap[require_constant.NameTransformSymbol.HYPHEN]]);
		if (withExt) decodedName = decodedName.replace(".js", "");
		return decodedName;
	} catch (err) {
		throw err;
	}
};
const generateExposeFilename = (exposeName, withExt) => {
	if (!exposeName) return "";
	let expose = exposeName;
	if (expose === ".") expose = "default_export";
	if (expose.startsWith("./")) expose = expose.replace("./", "");
	return encodeName(expose, "__federation_expose_", withExt);
};
const generateShareFilename = (pkgName, withExt) => {
	if (!pkgName) return "";
	return encodeName(pkgName, "__federation_shared_", withExt);
};
const getResourceUrl = (module, sourceUrl) => {
	if ("getPublicPath" in module) {
		let publicPath;
		if (!module.getPublicPath.startsWith("function")) publicPath = new Function(module.getPublicPath)();
		else publicPath = new Function("return " + module.getPublicPath)()();
		return `${publicPath}${sourceUrl}`;
	} else if ("publicPath" in module) {
		if (!require_env.isBrowserEnv() && !require_env.isReactNativeEnv() && "ssrPublicPath" in module && typeof module.ssrPublicPath === "string") return `${module.ssrPublicPath}${sourceUrl}`;
		return `${module.publicPath}${sourceUrl}`;
	} else {
		console.warn("Cannot get resource URL. If in debug mode, please ignore.", module, sourceUrl);
		return "";
	}
};
const assert = (condition, msg) => {
	if (!condition) error(msg);
};
const error = (msg) => {
	throw new Error(`${LOG_CATEGORY}: ${msg}`);
};
const warn = (msg) => {
	console.warn(`${LOG_CATEGORY}: ${msg}`);
};
function safeToString(info) {
	try {
		return JSON.stringify(info, null, 2);
	} catch (e) {
		return "";
	}
}
const VERSION_PATTERN_REGEXP = /^([\d^=v<>~]|[*xX]$)/;
function isRequiredVersion(str) {
	return VERSION_PATTERN_REGEXP.test(str);
}

//#endregion
exports.assert = assert;
exports.composeKeyWithSeparator = composeKeyWithSeparator;
exports.decodeName = decodeName;
exports.encodeName = encodeName;
exports.error = error;
exports.generateExposeFilename = generateExposeFilename;
exports.generateShareFilename = generateShareFilename;
exports.getResourceUrl = getResourceUrl;
exports.isRequiredVersion = isRequiredVersion;
exports.parseEntry = parseEntry;
exports.safeToString = safeToString;
exports.warn = warn;
//# sourceMappingURL=utils.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/_virtual/_rolldown/runtime.cjs"(__unused_rspack_module, exports) {
//#region \0rolldown/runtime.js
var __create = Object.create;
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __getProtoOf = Object.getPrototypeOf;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __copyProps = (to, from, except, desc) => {
	if (from && typeof from === "object" || typeof from === "function") {
		for (var keys = __getOwnPropNames(from), i = 0, n = keys.length, key; i < n; i++) {
			key = keys[i];
			if (!__hasOwnProp.call(to, key) && key !== except) {
				__defProp(to, key, {
					get: ((k) => from[k]).bind(null, key),
					enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable
				});
			}
		}
	}
	return to;
};
var __toESM = (mod, isNodeMode, target) => (target = mod != null ? __create(__getProtoOf(mod)) : {}, __copyProps(isNodeMode || !mod || !mod.__esModule ? __defProp(target, "default", {
	value: mod,
	enumerable: true
}) : target, mod));

//#endregion

exports.__toESM = __toESM;

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/attachShareScopeMap.cjs"(__unused_rspack_module, exports) {

//#region src/attachShareScopeMap.ts
function attachShareScopeMap(webpackRequire) {
	if (!webpackRequire.S || webpackRequire.federation.hasAttachShareScopeMap || !webpackRequire.federation.instance || !webpackRequire.federation.instance.shareScopeMap) return;
	webpackRequire.S = webpackRequire.federation.instance.shareScopeMap;
	webpackRequire.federation.hasAttachShareScopeMap = true;
}

//#endregion
exports.attachShareScopeMap = attachShareScopeMap;
//# sourceMappingURL=attachShareScopeMap.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/constant.cjs"(__unused_rspack_module, exports) {
Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });

//#region src/constant.ts
const FEDERATION_SUPPORTED_TYPES = ["script"];

//#endregion
exports.FEDERATION_SUPPORTED_TYPES = FEDERATION_SUPPORTED_TYPES;
//# sourceMappingURL=constant.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/consumes.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_attachShareScopeMap = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/attachShareScopeMap.cjs");
const require_updateOptions = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/updateOptions.cjs");
const require_getUsedExports = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/getUsedExports.cjs");
const require_error_codes = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/error-codes/dist/error-codes.cjs");
const require_getShortErrorMsg = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/error-codes/dist/getShortErrorMsg.cjs");

//#region src/consumes.ts
function consumes(options) {
	require_updateOptions.updateConsumeOptions(options);
	const { chunkId, promises, installedModules, webpackRequire, chunkMapping, moduleToHandlerMapping } = options;
	require_attachShareScopeMap.attachShareScopeMap(webpackRequire);
	if (webpackRequire.o(chunkMapping, chunkId)) chunkMapping[chunkId].forEach((id) => {
		if (webpackRequire.o(installedModules, id)) return promises.push(installedModules[id]);
		const onFactory = (factory) => {
			installedModules[id] = 0;
			webpackRequire.m[id] = (module) => {
				delete webpackRequire.c[id];
				const result = factory();
				const { shareInfo } = moduleToHandlerMapping[id];
				if (shareInfo?.shareConfig?.layer && result && typeof result === "object") try {
					if (!result.hasOwnProperty("layer") || result.layer === void 0) result.layer = shareInfo.shareConfig.layer;
				} catch (e) {}
				module.exports = result;
			};
		};
		const onError = (error) => {
			delete installedModules[id];
			webpackRequire.m[id] = (module) => {
				delete webpackRequire.c[id];
				throw error;
			};
		};
		try {
			const federationInstance = webpackRequire.federation.instance;
			if (!federationInstance) throw new Error("Federation instance not found!");
			const { shareKey, getter, shareInfo, treeShakingGetter } = moduleToHandlerMapping[id];
			const usedExports = require_getUsedExports.getUsedExports(webpackRequire, shareKey);
			const customShareInfo = { ...shareInfo };
			if (Array.isArray(customShareInfo.scope) && Array.isArray(customShareInfo.scope[0])) customShareInfo.scope = customShareInfo.scope[0];
			if (usedExports) customShareInfo.treeShaking = {
				usedExports,
				useIn: [federationInstance.options.name]
			};
			const promise = federationInstance.loadShare(shareKey, { customShareInfo }).then((factory) => {
				if (factory === false) {
					if (typeof getter !== "function") throw new Error(require_getShortErrorMsg.getShortErrorMsg(require_error_codes.RUNTIME_012, { [require_error_codes.RUNTIME_012]: "The getter for the shared module is not a function. This may be caused by setting \"shared.import: false\" without the host providing the corresponding lib." }, { shareKey }));
					return treeShakingGetter?.() || getter();
				}
				return factory;
			});
			if (promise.then) promises.push(installedModules[id] = promise.then(onFactory).catch(onError));
			else onFactory(promise);
		} catch (e) {
			onError(e);
		}
	});
}

//#endregion
exports.consumes = consumes;
//# sourceMappingURL=consumes.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/error-codes/dist/error-codes.cjs"(__unused_rspack_module, exports) {

//#region ../error-codes/dist/error-codes.mjs
const RUNTIME_012 = "RUNTIME-012";

//#endregion
exports.RUNTIME_012 = RUNTIME_012;
//# sourceMappingURL=error-codes.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/error-codes/dist/getShortErrorMsg.cjs"(__unused_rspack_module, exports) {

//#region ../error-codes/dist/getShortErrorMsg.mjs
const getDocsUrl = (errorCode) => {
	return `View the docs to see how to solve: https://module-federation.io/guide/troubleshooting/${errorCode.split("-")[0].toLowerCase()}#${errorCode.toLowerCase()}`;
};
const getShortErrorMsg = (errorCode, errorDescMap, args, originalErrorMsg) => {
	const msg = [`${[errorDescMap[errorCode]]} #${errorCode}`];
	args && msg.push(`args: ${JSON.stringify(args)}`);
	msg.push(getDocsUrl(errorCode));
	originalErrorMsg && msg.push(`Original Error Message:\n ${originalErrorMsg}`);
	return msg.join("\n");
};

//#endregion
exports.getShortErrorMsg = getShortErrorMsg;
//# sourceMappingURL=getShortErrorMsg.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/getSharedFallbackGetter.cjs"(__unused_rspack_module, exports) {

//#region src/getSharedFallbackGetter.ts
const getSharedFallbackGetter = ({ shareKey, factory, version, webpackRequire, libraryType = "global" }) => {
	const { runtime, instance, bundlerRuntime, sharedFallback } = webpackRequire.federation;
	if (!sharedFallback) return factory;
	const fallbackItems = sharedFallback[shareKey];
	if (!fallbackItems) return factory;
	const fallbackItem = version ? fallbackItems.find((item) => item[1] === version) : fallbackItems[0];
	if (!fallbackItem) throw new Error(`No fallback item found for shareKey: ${shareKey} and version: ${version}`);
	return () => runtime.getRemoteEntry({
		origin: webpackRequire.federation.instance,
		remoteInfo: {
			name: fallbackItem[2],
			entry: `${webpackRequire.p}${fallbackItem[0]}`,
			type: libraryType,
			entryGlobalName: fallbackItem[2],
			shareScope: "default"
		}
	}).then((shareEntry) => {
		if (!shareEntry) throw new Error(`Failed to load fallback entry for shareKey: ${shareKey} and version: ${version}`);
		return shareEntry.init(webpackRequire.federation.instance, bundlerRuntime).then(() => shareEntry.get());
	});
};

//#endregion
exports.getSharedFallbackGetter = getSharedFallbackGetter;
//# sourceMappingURL=getSharedFallbackGetter.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/getUsedExports.cjs"(__unused_rspack_module, exports) {

//#region src/getUsedExports.ts
function getUsedExports(webpackRequire, sharedName) {
	const usedExports = webpackRequire.federation.usedExports;
	if (!usedExports) return;
	return usedExports[sharedName];
}

//#endregion
exports.getUsedExports = getUsedExports;
//# sourceMappingURL=getUsedExports.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/index.cjs"(__unused_rspack_module, exports, __webpack_require__) {
Object.defineProperties(exports, { __esModule: { value: true }, [Symbol.toStringTag]: { value: 'Module' } });
const require_runtime = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/_virtual/_rolldown/runtime.cjs");
const require_attachShareScopeMap = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/attachShareScopeMap.cjs");
const require_remotes = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/remotes.cjs");
const require_consumes = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/consumes.cjs");
const require_initializeSharing = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/initializeSharing.cjs");
const require_installInitialConsumes = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/installInitialConsumes.cjs");
const require_initContainerEntry = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/initContainerEntry.cjs");
const require_init = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/init.cjs");
const require_getSharedFallbackGetter = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/getSharedFallbackGetter.cjs");
let _module_federation_runtime = __webpack_require__("../../node_modules/@module-federation/runtime/dist/index.cjs");
_module_federation_runtime = require_runtime.__toESM(_module_federation_runtime);

//#region src/index.ts
const federation = {
	runtime: _module_federation_runtime,
	instance: void 0,
	initOptions: void 0,
	bundlerRuntime: {
		remotes: require_remotes.remotes,
		consumes: require_consumes.consumes,
		I: require_initializeSharing.initializeSharing,
		S: {},
		installInitialConsumes: require_installInitialConsumes.installInitialConsumes,
		initContainerEntry: require_initContainerEntry.initContainerEntry,
		init: require_init.init,
		getSharedFallbackGetter: require_getSharedFallbackGetter.getSharedFallbackGetter
	},
	attachShareScopeMap: require_attachShareScopeMap.attachShareScopeMap,
	bundlerRuntimeOptions: {}
};
const instance = federation.instance;
const initOptions = federation.initOptions;
const bundlerRuntime = federation.bundlerRuntime;
const bundlerRuntimeOptions = federation.bundlerRuntimeOptions;

//#endregion
exports.attachShareScopeMap = require_attachShareScopeMap.attachShareScopeMap;
exports.bundlerRuntime = bundlerRuntime;
exports.bundlerRuntimeOptions = bundlerRuntimeOptions;
exports["default"] = federation;
exports.initOptions = initOptions;
exports.instance = instance;
Object.defineProperty(exports, "runtime", ({
  enumerable: true,
  get: function () {
    return _module_federation_runtime;
  }
}));
//# sourceMappingURL=index.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/init.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/_virtual/_rolldown/runtime.cjs");
let _module_federation_runtime = __webpack_require__("../../node_modules/@module-federation/runtime/dist/index.cjs");
let _module_federation_runtime_helpers = __webpack_require__("../../node_modules/@module-federation/runtime/dist/helpers.cjs");
_module_federation_runtime_helpers = require_runtime.__toESM(_module_federation_runtime_helpers);

//#region src/init.ts
function init({ webpackRequire }) {
	const { initOptions, runtime, sharedFallback, bundlerRuntime, libraryType } = webpackRequire.federation;
	if (!initOptions) throw new Error("initOptions is required!");
	const treeShakingSharePlugin = function() {
		return {
			name: "tree-shake-plugin",
			beforeInit(args) {
				const { userOptions, origin, options: registeredOptions } = args;
				const version = userOptions.version || registeredOptions.version;
				if (!sharedFallback) return args;
				const currentShared = userOptions.shared || {};
				const shared = [];
				Object.keys(currentShared).forEach((sharedName) => {
					(Array.isArray(currentShared[sharedName]) ? currentShared[sharedName] : [currentShared[sharedName]]).forEach((sharedArg) => {
						shared.push([sharedName, sharedArg]);
						if ("get" in sharedArg) {
							sharedArg.treeShaking ||= {};
							sharedArg.treeShaking.get = sharedArg.get;
							sharedArg.get = bundlerRuntime.getSharedFallbackGetter({
								shareKey: sharedName,
								factory: sharedArg.get,
								webpackRequire,
								libraryType,
								version: sharedArg.version
							});
						}
					});
				});
				const hostGlobalSnapshot = _module_federation_runtime_helpers.default.global.getGlobalSnapshotInfoByModuleInfo({
					name: origin.name,
					version
				});
				if (!hostGlobalSnapshot || !("shared" in hostGlobalSnapshot)) return args;
				Object.keys(registeredOptions.shared || {}).forEach((pkgName) => {
					registeredOptions.shared[pkgName].forEach((sharedArg) => {
						shared.push([pkgName, sharedArg]);
					});
				});
				const patchShared = (pkgName, shared) => {
					const shareSnapshot = hostGlobalSnapshot.shared.find((item) => item.sharedName === pkgName);
					if (!shareSnapshot) return;
					const { treeShaking } = shared;
					if (!treeShaking) return;
					const { secondarySharedTreeShakingName, secondarySharedTreeShakingEntry, treeShakingStatus } = shareSnapshot;
					if (treeShaking.status === treeShakingStatus) return;
					treeShaking.status = treeShakingStatus;
					if (secondarySharedTreeShakingEntry && libraryType && secondarySharedTreeShakingName) treeShaking.get = async () => {
						const shareEntry = await (0, _module_federation_runtime.getRemoteEntry)({
							origin,
							remoteInfo: {
								name: secondarySharedTreeShakingName,
								entry: secondarySharedTreeShakingEntry,
								type: libraryType,
								entryGlobalName: secondarySharedTreeShakingName,
								shareScope: "default"
							}
						});
						await shareEntry.init(origin, __webpack_require__.federation.bundlerRuntime);
						return shareEntry.get();
					};
				};
				shared.forEach(([pkgName, sharedArg]) => {
					patchShared(pkgName, sharedArg);
				});
				return args;
			}
		};
	};
	initOptions.plugins ||= [];
	initOptions.plugins.push(treeShakingSharePlugin());
	return runtime.init(initOptions);
}

//#endregion
exports.init = init;
//# sourceMappingURL=init.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/initContainerEntry.cjs"(__unused_rspack_module, exports) {

//#region src/initContainerEntry.ts
function initContainerEntry(options) {
	const { webpackRequire, shareScope, initScope, shareScopeKey, remoteEntryInitOptions } = options;
	if (!webpackRequire.S) return;
	if (!webpackRequire.federation || !webpackRequire.federation.instance || !webpackRequire.federation.initOptions) return;
	const federationInstance = webpackRequire.federation.instance;
	federationInstance.initOptions({
		name: webpackRequire.federation.initOptions.name,
		remotes: [],
		...remoteEntryInitOptions
	});
	const hostShareScopeKeys = remoteEntryInitOptions?.shareScopeKeys;
	const hostShareScopeMap = remoteEntryInitOptions?.shareScopeMap;
	if (!shareScopeKey || typeof shareScopeKey === "string") {
		const key = shareScopeKey || "default";
		if (Array.isArray(hostShareScopeKeys)) hostShareScopeKeys.forEach((hostKey) => {
			if (!hostShareScopeMap[hostKey]) hostShareScopeMap[hostKey] = {};
			const sc = hostShareScopeMap[hostKey];
			federationInstance.initShareScopeMap(hostKey, sc, { hostShareScopeMap: remoteEntryInitOptions?.shareScopeMap || {} });
		});
		else federationInstance.initShareScopeMap(key, shareScope, { hostShareScopeMap: remoteEntryInitOptions?.shareScopeMap || {} });
	} else shareScopeKey.forEach((key) => {
		if (!hostShareScopeKeys || !hostShareScopeMap) {
			federationInstance.initShareScopeMap(key, shareScope, { hostShareScopeMap: remoteEntryInitOptions?.shareScopeMap || {} });
			return;
		}
		if (!hostShareScopeMap[key]) hostShareScopeMap[key] = {};
		const sc = hostShareScopeMap[key];
		federationInstance.initShareScopeMap(key, sc, { hostShareScopeMap: remoteEntryInitOptions?.shareScopeMap || {} });
	});
	if (webpackRequire.federation.attachShareScopeMap) webpackRequire.federation.attachShareScopeMap(webpackRequire);
	if (!Array.isArray(shareScopeKey)) return webpackRequire.I(shareScopeKey || "default", initScope);
	if (Boolean(webpackRequire.federation.initOptions.shared)) return webpackRequire.I(shareScopeKey, initScope);
	return Promise.all(shareScopeKey.map((key) => {
		return webpackRequire.I(key, initScope);
	})).then(() => true);
}

//#endregion
exports.initContainerEntry = initContainerEntry;
//# sourceMappingURL=initContainerEntry.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/initializeSharing.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_attachShareScopeMap = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/attachShareScopeMap.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/constant.cjs");

//#region src/initializeSharing.ts
function initializeSharing({ shareScopeName, webpackRequire, initPromises, initTokens, initScope }) {
	const shareScopeKeys = Array.isArray(shareScopeName) ? shareScopeName : [shareScopeName];
	var initializeSharingPromises = [];
	var _initializeSharing = function(shareScopeKey) {
		if (!initScope) initScope = [];
		const mfInstance = webpackRequire.federation.instance;
		var initToken = initTokens[shareScopeKey];
		if (!initToken) initToken = initTokens[shareScopeKey] = { from: mfInstance.name };
		if (initScope.indexOf(initToken) >= 0) return;
		initScope.push(initToken);
		const promise = initPromises[shareScopeKey];
		if (promise) return promise;
		var warn = (msg) => typeof console !== "undefined" && console.warn && console.warn(msg);
		var initExternal = (id) => {
			var handleError = (err) => warn("Initialization of sharing external failed: " + err);
			try {
				var module = webpackRequire(id);
				if (!module) return;
				var initFn = (module) => module && module.init && module.init(webpackRequire.S[shareScopeKey], initScope, {
					shareScopeMap: webpackRequire.S || {},
					shareScopeKeys: shareScopeName
				});
				if (module.then) return promises.push(module.then(initFn, handleError));
				var initResult = initFn(module);
				if (initResult && typeof initResult !== "boolean" && initResult.then) return promises.push(initResult["catch"](handleError));
			} catch (err) {
				handleError(err);
			}
		};
		const promises = mfInstance.initializeSharing(shareScopeKey, {
			strategy: mfInstance.options.shareStrategy,
			initScope,
			from: "build"
		});
		require_attachShareScopeMap.attachShareScopeMap(webpackRequire);
		const bundlerRuntimeRemotesOptions = webpackRequire.federation.bundlerRuntimeOptions.remotes;
		if (bundlerRuntimeRemotesOptions) Object.keys(bundlerRuntimeRemotesOptions.idToRemoteMap).forEach((moduleId) => {
			const info = bundlerRuntimeRemotesOptions.idToRemoteMap[moduleId];
			const externalModuleId = bundlerRuntimeRemotesOptions.idToExternalAndNameMapping[moduleId][2];
			if (info.length > 1) initExternal(externalModuleId);
			else if (info.length === 1) {
				const remoteInfo = info[0];
				if (!require_constant.FEDERATION_SUPPORTED_TYPES.includes(remoteInfo.externalType)) initExternal(externalModuleId);
			}
		});
		if (!promises.length) return initPromises[shareScopeKey] = true;
		return initPromises[shareScopeKey] = Promise.all(promises).then(() => initPromises[shareScopeKey] = true);
	};
	shareScopeKeys.forEach((key) => {
		initializeSharingPromises.push(_initializeSharing(key));
	});
	return Promise.all(initializeSharingPromises).then(() => true);
}

//#endregion
exports.initializeSharing = initializeSharing;
//# sourceMappingURL=initializeSharing.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/installInitialConsumes.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_updateOptions = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/updateOptions.cjs");
const require_getUsedExports = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/getUsedExports.cjs");

//#region src/installInitialConsumes.ts
function handleInitialConsumes(options) {
	const { moduleId, moduleToHandlerMapping, webpackRequire, asyncLoad } = options;
	const federationInstance = webpackRequire.federation.instance;
	if (!federationInstance) throw new Error("Federation instance not found!");
	const { shareKey, shareInfo } = moduleToHandlerMapping[moduleId];
	try {
		const usedExports = require_getUsedExports.getUsedExports(webpackRequire, shareKey);
		const customShareInfo = { ...shareInfo };
		if (usedExports) customShareInfo.treeShaking = {
			usedExports,
			useIn: [federationInstance.options.name]
		};
		if (asyncLoad) return federationInstance.loadShare(shareKey, { customShareInfo });
		return federationInstance.loadShareSync(shareKey, { customShareInfo });
	} catch (err) {
		console.error("loadShareSync failed! The function should not be called unless you set \"eager:true\". If you do not set it, and encounter this issue, you can check whether an async boundary is implemented.");
		console.error("The original error message is as follows: ");
		throw err;
	}
}
function installInitialConsumes(options) {
	require_updateOptions.updateConsumeOptions(options);
	const { moduleToHandlerMapping, webpackRequire, installedModules, initialConsumes, asyncLoad } = options;
	const factoryIdSets = [];
	initialConsumes.forEach((id) => {
		const factoryGetter = () => handleInitialConsumes({
			moduleId: id,
			moduleToHandlerMapping,
			webpackRequire,
			asyncLoad
		});
		factoryIdSets.push([id, factoryGetter]);
	});
	const setModule = (id, factoryGetter) => {
		webpackRequire.m[id] = (module) => {
			installedModules[id] = 0;
			delete webpackRequire.c[id];
			const factory = factoryGetter();
			if (typeof factory !== "function") throw new Error(`Shared module is not available for eager consumption: ${id}`);
			const result = factory();
			const { shareInfo } = moduleToHandlerMapping[id];
			if (shareInfo?.shareConfig?.layer && result && typeof result === "object") try {
				if (!result.hasOwnProperty("layer") || result.layer === void 0) result.layer = shareInfo.shareConfig.layer;
			} catch (e) {}
			module.exports = result;
		};
	};
	if (asyncLoad) return Promise.all(factoryIdSets.map(async ([id, factoryGetter]) => {
		const result = await factoryGetter();
		setModule(id, () => result);
	}));
	factoryIdSets.forEach(([id, factoryGetter]) => {
		setModule(id, factoryGetter);
	});
}

//#endregion
exports.installInitialConsumes = installInitialConsumes;
//# sourceMappingURL=installInitialConsumes.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/remotes.cjs"(__unused_rspack_module, exports, __webpack_require__) {
const require_runtime = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/_virtual/_rolldown/runtime.cjs");
const require_attachShareScopeMap = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/attachShareScopeMap.cjs");
const require_constant = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/constant.cjs");
const require_updateOptions = __webpack_require__("../../node_modules/@module-federation/webpack-bundler-runtime/dist/updateOptions.cjs");
let _module_federation_sdk = __webpack_require__("../../node_modules/@module-federation/sdk/dist/index.cjs");

//#region src/remotes.ts
function remotes(options) {
	require_updateOptions.updateRemoteOptions(options);
	const { chunkId, promises, webpackRequire, chunkMapping, idToExternalAndNameMapping, idToRemoteMap } = options;
	require_attachShareScopeMap.attachShareScopeMap(webpackRequire);
	if (webpackRequire.o(chunkMapping, chunkId)) chunkMapping[chunkId].forEach((id) => {
		let getScope = webpackRequire.R;
		if (!getScope) getScope = [];
		const data = idToExternalAndNameMapping[id];
		const remoteInfos = idToRemoteMap[id] || [];
		if (getScope.indexOf(data) >= 0) return;
		getScope.push(data);
		if (data.p) return promises.push(data.p);
		const onError = (error) => {
			if (!error) error = /* @__PURE__ */ new Error("Container missing");
			if (typeof error.message === "string") error.message += `\nwhile loading "${data[1]}" from ${data[2]}`;
			webpackRequire.m[id] = () => {
				throw error;
			};
			data.p = 0;
		};
		const handleFunction = (fn, arg1, arg2, d, next, first) => {
			try {
				const promise = fn(arg1, arg2);
				if (promise && promise.then) {
					const p = promise.then((result) => next(result, d), onError);
					if (first) promises.push(data.p = p);
					else return p;
				} else return next(promise, d, first);
			} catch (error) {
				onError(error);
			}
		};
		const onExternal = (external, _, first) => external ? handleFunction(webpackRequire.I, data[0], 0, external, onInitialized, first) : onError();
		var onInitialized = (_, external, first) => handleFunction(external.get, data[1], getScope, 0, onFactory, first);
		var onFactory = (factory) => {
			data.p = 1;
			webpackRequire.m[id] = (module) => {
				module.exports = factory();
			};
		};
		const onRemoteLoaded = () => {
			try {
				const remoteModuleName = (0, _module_federation_sdk.decodeName)(remoteInfos[0].name, _module_federation_sdk.ENCODE_NAME_PREFIX) + data[1].slice(1);
				const instance = webpackRequire.federation.instance;
				const loadRemote = () => webpackRequire.federation.instance.loadRemote(remoteModuleName, {
					loadFactory: false,
					from: "build"
				});
				if (instance.options.shareStrategy === "version-first") {
					const shareScopes = Array.isArray(data[0]) ? data[0] : [data[0]];
					return Promise.all(shareScopes.map((shareScope) => instance.sharedHandler.initializeSharing(shareScope))).then(() => {
						return loadRemote();
					});
				}
				return loadRemote();
			} catch (error) {
				onError(error);
			}
		};
		if (remoteInfos.length === 1 && require_constant.FEDERATION_SUPPORTED_TYPES.includes(remoteInfos[0].externalType) && remoteInfos[0].name) handleFunction(onRemoteLoaded, data[2], 0, 0, onFactory, 1);
		else handleFunction(webpackRequire, data[2], 0, 0, onExternal, 1);
	});
}

//#endregion
exports.remotes = remotes;
//# sourceMappingURL=remotes.cjs.map

},
"../../node_modules/@module-federation/webpack-bundler-runtime/dist/updateOptions.cjs"(__unused_rspack_module, exports) {

//#region src/updateOptions.ts
function updateConsumeOptions(options) {
	const { webpackRequire, moduleToHandlerMapping } = options;
	const { consumesLoadingData, initializeSharingData } = webpackRequire;
	const { sharedFallback, bundlerRuntime, libraryType } = webpackRequire.federation;
	if (consumesLoadingData && !consumesLoadingData._updated) {
		const { moduleIdToConsumeDataMapping: updatedModuleIdToConsumeDataMapping = {}, initialConsumes: updatedInitialConsumes = [], chunkMapping: updatedChunkMapping = {} } = consumesLoadingData;
		Object.entries(updatedModuleIdToConsumeDataMapping).forEach(([id, data]) => {
			if (!moduleToHandlerMapping[id]) moduleToHandlerMapping[id] = {
				getter: sharedFallback ? bundlerRuntime?.getSharedFallbackGetter({
					shareKey: data.shareKey,
					factory: data.fallback,
					webpackRequire,
					libraryType
				}) : data.fallback,
				treeShakingGetter: sharedFallback ? data.fallback : void 0,
				shareInfo: {
					shareConfig: {
						requiredVersion: data.requiredVersion,
						strictVersion: data.strictVersion,
						singleton: data.singleton,
						eager: data.eager,
						layer: data.layer
					},
					scope: Array.isArray(data.shareScope) ? data.shareScope : [data.shareScope || "default"],
					treeShaking: sharedFallback ? {
						get: data.fallback,
						mode: data.treeShakingMode
					} : void 0
				},
				shareKey: data.shareKey
			};
		});
		if ("initialConsumes" in options) {
			const { initialConsumes = [] } = options;
			updatedInitialConsumes.forEach((id) => {
				if (!initialConsumes.includes(id)) initialConsumes.push(id);
			});
		}
		if ("chunkMapping" in options) {
			const { chunkMapping = {} } = options;
			Object.entries(updatedChunkMapping).forEach(([id, chunkModules]) => {
				if (!chunkMapping[id]) chunkMapping[id] = [];
				chunkModules.forEach((moduleId) => {
					if (!chunkMapping[id].includes(moduleId)) chunkMapping[id].push(moduleId);
				});
			});
		}
		consumesLoadingData._updated = 1;
	}
	if (initializeSharingData && !initializeSharingData._updated) {
		const { federation } = webpackRequire;
		if (!federation.instance || !initializeSharingData.scopeToSharingDataMapping) return;
		const shared = {};
		for (let [scope, stages] of Object.entries(initializeSharingData.scopeToSharingDataMapping)) for (let stage of stages) if (typeof stage === "object" && stage !== null) {
			const { name, version, factory, eager, singleton, requiredVersion, strictVersion } = stage;
			const shareConfig = { requiredVersion: `^${version}` };
			const isValidValue = function(val) {
				return typeof val !== "undefined";
			};
			if (isValidValue(singleton)) shareConfig.singleton = singleton;
			if (isValidValue(requiredVersion)) shareConfig.requiredVersion = requiredVersion;
			if (isValidValue(eager)) shareConfig.eager = eager;
			if (isValidValue(strictVersion)) shareConfig.strictVersion = strictVersion;
			const options = {
				version,
				scope: [scope],
				shareConfig,
				get: factory
			};
			if (shared[name]) shared[name].push(options);
			else shared[name] = [options];
		}
		federation.instance.registerShared(shared);
		initializeSharingData._updated = 1;
	}
}
function updateRemoteOptions(options) {
	const { webpackRequire, idToExternalAndNameMapping = {}, idToRemoteMap = {}, chunkMapping = {} } = options;
	const { remotesLoadingData } = webpackRequire;
	const remoteInfos = webpackRequire.federation?.bundlerRuntimeOptions?.remotes?.remoteInfos;
	if (!remotesLoadingData || remotesLoadingData._updated || !remoteInfos) return;
	const { chunkMapping: updatedChunkMapping, moduleIdToRemoteDataMapping } = remotesLoadingData;
	if (!updatedChunkMapping || !moduleIdToRemoteDataMapping) return;
	for (let [moduleId, data] of Object.entries(moduleIdToRemoteDataMapping)) {
		if (!idToExternalAndNameMapping[moduleId]) idToExternalAndNameMapping[moduleId] = [
			data.shareScope,
			data.name,
			data.externalModuleId
		];
		if (!idToRemoteMap[moduleId] && remoteInfos[data.remoteName]) {
			const items = remoteInfos[data.remoteName];
			idToRemoteMap[moduleId] ||= [];
			items.forEach((item) => {
				if (!idToRemoteMap[moduleId].includes(item)) idToRemoteMap[moduleId].push(item);
			});
		}
	}
	if (chunkMapping) Object.entries(updatedChunkMapping).forEach(([id, chunkModules]) => {
		if (!chunkMapping[id]) chunkMapping[id] = [];
		chunkModules.forEach((moduleId) => {
			if (!chunkMapping[id].includes(moduleId)) chunkMapping[id].push(moduleId);
		});
	});
	remotesLoadingData._updated = 1;
}

//#endregion
exports.updateConsumeOptions = updateConsumeOptions;
exports.updateRemoteOptions = updateRemoteOptions;
//# sourceMappingURL=updateOptions.cjs.map

},

});
// The module cache
var __webpack_module_cache__ = {};

// The require function
function __webpack_require__(moduleId) {

// Check if module is in cache
var cachedModule = __webpack_module_cache__[moduleId];
if (cachedModule !== undefined) {
return cachedModule.exports;
}
// Create a new module (and put it into the cache)
var module = (__webpack_module_cache__[moduleId] = {
id: moduleId,
exports: {}
});
// Execute the module function
__webpack_modules__[moduleId](module, module.exports, __webpack_require__);

// Return the exports of the module
return module.exports;

}

// expose the modules object (__webpack_modules__)
__webpack_require__.m = __webpack_modules__;

// expose the module cache
__webpack_require__.c = __webpack_module_cache__;

// the startup function
__webpack_require__.x = () => {
// Load entry module and return exports
var __webpack_exports__ = __webpack_require__("webpack/container/entry/@jupyterlab/galata-extension");
return __webpack_exports__
};

// webpack/runtime/compat_get_default_export
(() => {
// getDefaultExport function for compatibility with non-ESM modules
__webpack_require__.n = (module) => {
	var getter = module && module.__esModule ?
		() => (module['default']) :
		() => (module);
	__webpack_require__.d(getter, { a: getter });
	return getter;
};

})();
// webpack/runtime/define_property_getters
(() => {
__webpack_require__.d = (exports, definition) => {
	for(var key in definition) {
        if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
            Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
        }
    }
};
})();
// webpack/runtime/ensure_chunk
(() => {
__webpack_require__.f = {};
// This file contains only the entry chunk.
// The chunk loading function for additional chunks
__webpack_require__.e = (chunkId) => {
	return Promise.all(
		Object.keys(__webpack_require__.f).reduce((promises, key) => {
			__webpack_require__.f[key](chunkId, promises);
			return promises;
		}, [])
	);
};
})();
// webpack/runtime/get javascript chunk filename
(() => {
// This function allow to reference chunks
__webpack_require__.u = (chunkId) => {
  // return url for filenames not based on template
  
  // return url for filenames based on template
  return "" + chunkId + "." + {"lib_extension_index_js": "573abb5f68d369fe","vendors-node_modules_fontsource-variable_noto-sans-sc_index_css": "386da7677f8b1f52","vendors-node_modules_fontsource_dejavu-mono_index_css": "d98af1128b7590b9","vendors-node_modules_fontsource_dejavu-sans_index_css": "bf5a7bbac03d188c",}[chunkId] + ".js"
}
})();
// webpack/runtime/global
(() => {
__webpack_require__.g = (() => {
	if (typeof globalThis === 'object') return globalThis;
	try {
		return this || new Function('return this')();
	} catch (e) {
		if (typeof window === 'object') return window;
	}
})();
})();
// webpack/runtime/has_own_property
(() => {
__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
})();
// webpack/runtime/load_script
(() => {
var inProgress = {};

var uniqueName = "@jupyterlab/galata-extension:";
// loadScript function to load a script via script tag
__webpack_require__.l = function (url, done, key, chunkId) {
	if (inProgress[url]) {
		inProgress[url].push(done);
		return;
	}
	var script, needAttach;
	if (key !== undefined) {
		var scripts = document.getElementsByTagName("script");
		for (var i = 0; i < scripts.length; i++) {
			var s = scripts[i];
			if (s.getAttribute("src") == url || s.getAttribute("data-rspack") == uniqueName + key) {
				script = s;
				break;
			}
		}
	}
	if (!script) {
		needAttach = true;
		script = document.createElement('script');

script.timeout = 120;
if (__webpack_require__.nc) {
  script.setAttribute("nonce", __webpack_require__.nc);
}

script.setAttribute("data-rspack", uniqueName + key);



script.src = url;


	}
	inProgress[url] = [done];
	var onScriptComplete = function (prev, event) {
		script.onerror = script.onload = null;
		clearTimeout(timeout);
		var doneFns = inProgress[url];
		delete inProgress[url];
		script.parentNode && script.parentNode.removeChild(script);
		doneFns &&
			doneFns.forEach(function (fn) {
				return fn(event);
			});
		if (prev) return prev(event);
	};
	var timeout = setTimeout(
		onScriptComplete.bind(null, undefined, {
			type: 'timeout',
			target: script
		}),
		120000
	);
	script.onerror = onScriptComplete.bind(null, script.onerror);
	script.onload = onScriptComplete.bind(null, script.onload);
	needAttach && document.head.appendChild(script);
};

})();
// webpack/runtime/make_namespace_object
(() => {
// define __esModule on exports
__webpack_require__.r = (exports) => {
	if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
		Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
	}
	Object.defineProperty(exports, '__esModule', { value: true });
};
})();
// webpack/runtime/module_federation/runtime
(() => {

if(!__webpack_require__.federation){
    __webpack_require__.federation = {
        
chunkMatcher: function(chunkId) {
    return true;
},
rootOutputDir: "",

    };
}

})();
// webpack/runtime/nonce
(() => {
__webpack_require__.nc = undefined;
})();
// webpack/runtime/sharing
(() => {

__webpack_require__.S = {};
__webpack_require__.initializeSharingData = { scopeToSharingDataMapping: { "default": [{ name: "@fontsource-variable/noto-sans-sc", version: "5.2.10", factory: () => (__webpack_require__.e("vendors-node_modules_fontsource-variable_noto-sans-sc_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource-variable/noto-sans-sc/index.css"))))), eager: 0, treeShakingMode: null }, { name: "@fontsource/dejavu-mono", version: "5.2.5", factory: () => (__webpack_require__.e("vendors-node_modules_fontsource_dejavu-mono_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource/dejavu-mono/index.css"))))), eager: 0, treeShakingMode: null }, { name: "@fontsource/dejavu-sans", version: "5.2.5", factory: () => (__webpack_require__.e("vendors-node_modules_fontsource_dejavu-sans_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource/dejavu-sans/index.css"))))), eager: 0, treeShakingMode: null }, { name: "@jupyterlab/galata-extension", version: "5.6.1", factory: () => (__webpack_require__.e("lib_extension_index_js").then(() => (() => (__webpack_require__("../lib/extension/index.js"))))), eager: 0, treeShakingMode: null, singleton: 1 }] }, uniqueName: "@jupyterlab/galata-extension" };
__webpack_require__.I = __webpack_require__.I || function() { throw new Error("should have __webpack_require__.I") }

})();
// webpack/runtime/auto_public_path
(() => {
var scriptUrl;

if (__webpack_require__.g.importScripts) scriptUrl = __webpack_require__.g.location + "";
var document = __webpack_require__.g.document;
if (!scriptUrl && document) {
  // Technically we could use `document.currentScript instanceof window.HTMLScriptElement`,
  // but an attacker could try to inject `<script>HTMLScriptElement = HTMLImageElement</script>`
  // and use `<img name="currentScript" src="https://attacker.controlled.server/"></img>`
  if (document.currentScript && document.currentScript.tagName.toUpperCase() === 'SCRIPT') scriptUrl = document.currentScript.src;
  if (!scriptUrl) {
    var scripts = document.getElementsByTagName("script");
    if (scripts.length) {
      var i = scripts.length - 1;
      while (i > -1 && (!scriptUrl || !/^http(s?):/.test(scriptUrl))) scriptUrl = scripts[i--].src;
    }
  }
}

// When supporting browsers where an automatic publicPath is not supported you must specify an output.publicPath manually via configuration",
// or pass an empty string ("") and set the __webpack_public_path__ variable from your code to use your own logic.',
if (!scriptUrl) throw new Error("Automatic publicPath is not supported in this browser");
scriptUrl = scriptUrl.replace(/^blob:/, "").replace(/#.*$/, "").replace(/\?.*$/, "").replace(/\/[^\/]+$/, "/");

__webpack_require__.p = scriptUrl;

})();
// webpack/runtime/consumes_loading
(() => {

__webpack_require__.consumesLoadingData = { chunkMapping: {"lib_extension_index_js":["webpack/sharing/consume/default/@jupyterlab/terminal","webpack/sharing/consume/default/@fontsource/dejavu-sans/@fontsource/dejavu-sans","webpack/sharing/consume/default/@jupyterlab/apputils","webpack/sharing/consume/default/@fontsource/dejavu-mono/@fontsource/dejavu-mono","webpack/sharing/consume/default/@fontsource-variable/noto-sans-sc/@fontsource-variable/noto-sans-sc","webpack/sharing/consume/default/@lumino/signaling","webpack/sharing/consume/default/@lumino/algorithm","webpack/sharing/consume/default/@lumino/coreutils"]}, moduleIdToConsumeDataMapping: {"webpack/sharing/consume/default/@fontsource/dejavu-sans/@fontsource/dejavu-sans": { shareScope: "default", shareKey: "@fontsource/dejavu-sans", import: "@fontsource/dejavu-sans", requiredVersion: "*", strictVersion: true, singleton: false, eager: false, fallback: () => (__webpack_require__.e("vendors-node_modules_fontsource_dejavu-sans_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource/dejavu-sans/index.css"))))), treeShakingMode: null }, "webpack/sharing/consume/default/@jupyterlab/terminal": { shareScope: "default", shareKey: "@jupyterlab/terminal", import: null, requiredVersion: "^4.6.1", strictVersion: false, singleton: true, eager: false, fallback: undefined, treeShakingMode: null }, "webpack/sharing/consume/default/@jupyterlab/apputils": { shareScope: "default", shareKey: "@jupyterlab/apputils", import: null, requiredVersion: "^4.7.1", strictVersion: false, singleton: true, eager: false, fallback: undefined, treeShakingMode: null }, "webpack/sharing/consume/default/@fontsource-variable/noto-sans-sc/@fontsource-variable/noto-sans-sc": { shareScope: "default", shareKey: "@fontsource-variable/noto-sans-sc", import: "@fontsource-variable/noto-sans-sc", requiredVersion: "*", strictVersion: true, singleton: false, eager: false, fallback: () => (__webpack_require__.e("vendors-node_modules_fontsource-variable_noto-sans-sc_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource-variable/noto-sans-sc/index.css"))))), treeShakingMode: null }, "webpack/sharing/consume/default/@fontsource/dejavu-mono/@fontsource/dejavu-mono": { shareScope: "default", shareKey: "@fontsource/dejavu-mono", import: "@fontsource/dejavu-mono", requiredVersion: "*", strictVersion: true, singleton: false, eager: false, fallback: () => (__webpack_require__.e("vendors-node_modules_fontsource_dejavu-mono_index_css").then(() => (() => (__webpack_require__("../../node_modules/@fontsource/dejavu-mono/index.css"))))), treeShakingMode: null }, "webpack/sharing/consume/default/@lumino/signaling": { shareScope: "default", shareKey: "@lumino/signaling", import: null, requiredVersion: "^2.0.0", strictVersion: false, singleton: true, eager: false, fallback: undefined, treeShakingMode: null }, "webpack/sharing/consume/default/@lumino/coreutils": { shareScope: "default", shareKey: "@lumino/coreutils", import: null, requiredVersion: "^2.0.0", strictVersion: false, singleton: true, eager: false, fallback: undefined, treeShakingMode: null }, "webpack/sharing/consume/default/@lumino/algorithm": { shareScope: "default", shareKey: "@lumino/algorithm", import: null, requiredVersion: "^2.0.0", strictVersion: false, singleton: true, eager: false, fallback: undefined, treeShakingMode: null }}, initialConsumes: [] };
__webpack_require__.f.consumes = __webpack_require__.f.consumes || function() { throw new Error("should have __webpack_require__.f.consumes") }
})();
// webpack/runtime/initialize_exposes
(() => {

    __webpack_require__.initializeExposesData = {
  moduleMap: {
  "./index": () => {
return __webpack_require__.e("lib_extension_index_js").then(() => (() => (__webpack_require__("../lib/extension/index.js"))));
},
"./extension": () => {
return __webpack_require__.e("lib_extension_index_js").then(() => (() => (__webpack_require__("../lib/extension/index.js"))));
},
},
  shareScope: "default",
};
__webpack_require__.getContainer = __webpack_require__.getContainer || function() { throw new Error("should have __webpack_require__.getContainer") };__webpack_require__.initContainer = __webpack_require__.initContainer || function() { throw new Error("should have __webpack_require__.initContainer") };
})();
// webpack/runtime/jsonp_chunk_loading
(() => {
__webpack_require__.b = document.baseURI || self.location.href;

      // object to store loaded and loading chunks
      // undefined = chunk not loaded, null = chunk preloaded/prefetched
      // [resolve, reject, Promise] = chunk loading, 0 = chunk loaded
      var installedChunks = {"@jupyterlab/galata-extension": 0,};
      
        __webpack_require__.f.j = function (chunkId, promises) {
          // JSONP chunk loading for javascript
var installedChunkData = __webpack_require__.o(installedChunks, chunkId)
	? installedChunks[chunkId]
	: undefined;
if (installedChunkData !== 0) {
	// 0 means "already installed".

	// a Promise means "currently loading".
	if (installedChunkData) {
		promises.push(installedChunkData[2]);
	} else {
		if (true) {
			// setup Promise in chunk cache
			var promise = new Promise((resolve, reject) => (installedChunkData = installedChunks[chunkId] = [resolve, reject]));
			promises.push((installedChunkData[2] = promise));

			// start chunk loading
			var url = __webpack_require__.p + __webpack_require__.u(chunkId);
			// create error before stack unwound to get useful stacktrace later
			var error = new Error();
			var loadingEnded = function (event) {
				if (__webpack_require__.o(installedChunks, chunkId)) {
					installedChunkData = installedChunks[chunkId];
					if (installedChunkData !== 0) installedChunks[chunkId] = undefined;
					if (installedChunkData) {
						var errorType =
							event && (event.type === 'load' ? 'missing' : event.type);
						var realSrc = event && event.target && event.target.src;
						error.message =
							'Loading chunk ' +
							chunkId +
							' failed.\n(' +
							errorType +
							': ' +
							realSrc +
							')';
						error.name = 'ChunkLoadError';
						error.type = errorType;
						error.request = realSrc;
						installedChunkData[1](error);
					}
				}
			};
			__webpack_require__.l(url, loadingEnded, "chunk-" + chunkId, chunkId);
		} 
	}
}

        }
        // install a JSONP callback for chunk loading
var __rspack_jsonp = (parentChunkLoadingFunction, data) => {
	var [chunkIds, moreModules, runtime] = data;
	// add "moreModules" to the modules object,
	// then flag all "chunkIds" as loaded and fire callback
	var moduleId, chunkId, i = 0;
	if (chunkIds.some((id) => (installedChunks[id] !== 0))) {
		for (moduleId in moreModules) {
			if (__webpack_require__.o(moreModules, moduleId)) {
				__webpack_require__.m[moduleId] = moreModules[moduleId];
			}
		}
		if (runtime) var result = runtime(__webpack_require__);
	}
	if (parentChunkLoadingFunction) parentChunkLoadingFunction(data);
	for (; i < chunkIds.length; i++) {
		chunkId = chunkIds[i];
		if (
			__webpack_require__.o(installedChunks, chunkId) &&
			installedChunks[chunkId]
		) {
			installedChunks[chunkId][0]();
		}
		installedChunks[chunkId] = 0;
	}
	
};

var chunkLoadingGlobal = self["rspackChunk_jupyterlab_galata_extension"] = self["rspackChunk_jupyterlab_galata_extension"] || [];
chunkLoadingGlobal.forEach(__rspack_jsonp.bind(null, 0));
chunkLoadingGlobal.push = __rspack_jsonp.bind(null, chunkLoadingGlobal.push.bind(chunkLoadingGlobal));

})();
// webpack/runtime/embed_federation_runtime
(() => {
var prevStartup = __webpack_require__.x;
var hasRun = false;
__webpack_require__.x = function () {
	if (!hasRun) {
		hasRun = true;
		__webpack_require__("@module-federation/runtime/rspack.js!=!data:text/javascript,import%20__module_federation_bundler_runtime__%20from%20%22%2Fhome%2Frunner%2Fwork%2Fjupyterlab%2Fjupyterlab%2F.jupyter_releaser_checkout%2Fnode_modules%2F%40module-federation%2Fwebpack-bundler-runtime%2Fdist%2Findex.cjs%22%3Bconst%20__module_federation_runtime_plugins__%20%3D%20%5B%5D.filter((%7B%20plugin%20%7D)%20%3D%3E%20plugin).map((%7B%20plugin%2C%20params%20%7D)%20%3D%3E%20plugin(params))%3Bconst%20__module_federation_remote_infos__%20%3D%20%7B%7D%3Bconst%20__module_federation_container_name__%20%3D%20%22%40jupyterlab%2Fgalata-extension%22%3Bconst%20__module_federation_share_strategy__%20%3D%20%22version-first%22%3Bconst%20__module_federation_share_fallbacks__%20%3D%20undefined%3Bconst%20__module_federation_library_type__%20%3D%20%22var%22%3Bif((__webpack_require__.initializeSharingData%7C%7C__webpack_require__.initializeExposesData)%26%26__webpack_require__.federation)%7Bvar%20_ref%2C_ref1%2C_ref2%2C_ref3%2C_ref4%3Bvar%20__webpack_require___remotesLoadingData%2C__webpack_require___remotesLoadingData1%2C__webpack_require___initializeSharingData%2C__webpack_require___consumesLoadingData%2C__webpack_require___consumesLoadingData1%2C__webpack_require___initializeExposesData%2C__webpack_require___consumesLoadingData2%3Bconst%20override%3D(obj%2Ckey%2Cvalue)%3D%3E%7Bif(!obj)return%3Bif(obj%5Bkey%5D)obj%5Bkey%5D%3Dvalue%7D%3Bconst%20merge%3D(obj%2Ckey%2Cfn)%3D%3E%7Bconst%20value%3Dfn()%3Bif(Array.isArray(value))%7Bvar%20_obj%2C_key%2C_%3B(_%3D(_obj%3Dobj)%5B_key%3Dkey%5D)!%3D%3Dnull%26%26_!%3D%3Dvoid%200%3F_%3A_obj%5B_key%5D%3D%5B%5D%3Bobj%5Bkey%5D.push(...value)%7Delse%20if(typeof%20value%3D%3D%3D%22object%22%26%26value!%3D%3Dnull)%7Bvar%20_obj1%2C_key1%2C_1%3B(_1%3D(_obj1%3Dobj)%5B_key1%3Dkey%5D)!%3D%3Dnull%26%26_1!%3D%3Dvoid%200%3F_1%3A_obj1%5B_key1%5D%3D%7B%7D%3BObject.assign(obj%5Bkey%5D%2Cvalue)%7D%7D%3Bconst%20early%3D(obj%2Ckey%2Cinitial)%3D%3E%7Bvar%20_obj%2C_key%2C_%3B(_%3D(_obj%3Dobj)%5B_key%3Dkey%5D)!%3D%3Dnull%26%26_!%3D%3Dvoid%200%3F_%3A_obj%5B_key%5D%3Dinitial()%7D%3Bconst%20remotesLoadingChunkMapping%3D(_ref%3D(__webpack_require___remotesLoadingData%3D__webpack_require__.remotesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___remotesLoadingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___remotesLoadingData.chunkMapping)!%3D%3Dnull%26%26_ref!%3D%3Dvoid%200%3F_ref%3A%7B%7D%3Bconst%20remotesLoadingModuleIdToRemoteDataMapping%3D(_ref1%3D(__webpack_require___remotesLoadingData1%3D__webpack_require__.remotesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___remotesLoadingData1%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___remotesLoadingData1.moduleIdToRemoteDataMapping)!%3D%3Dnull%26%26_ref1!%3D%3Dvoid%200%3F_ref1%3A%7B%7D%3Bconst%20initializeSharingScopeToInitDataMapping%3D(_ref2%3D(__webpack_require___initializeSharingData%3D__webpack_require__.initializeSharingData)%3D%3D%3Dnull%7C%7C__webpack_require___initializeSharingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___initializeSharingData.scopeToSharingDataMapping)!%3D%3Dnull%26%26_ref2!%3D%3Dvoid%200%3F_ref2%3A%7B%7D%3Bconst%20consumesLoadingChunkMapping%3D(_ref3%3D(__webpack_require___consumesLoadingData%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData.chunkMapping)!%3D%3Dnull%26%26_ref3!%3D%3Dvoid%200%3F_ref3%3A%7B%7D%3Bconst%20consumesLoadingModuleToConsumeDataMapping%3D(_ref4%3D(__webpack_require___consumesLoadingData1%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData1%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData1.moduleIdToConsumeDataMapping)!%3D%3Dnull%26%26_ref4!%3D%3Dvoid%200%3F_ref4%3A%7B%7D%3Bconst%20consumesLoadinginstalledModules%3D%7B%7D%3Bconst%20initializeSharingInitPromises%3D%5B%5D%3Bconst%20initializeSharingInitTokens%3D%7B%7D%3Bconst%20containerShareScope%3D(__webpack_require___initializeExposesData%3D__webpack_require__.initializeExposesData)%3D%3D%3Dnull%7C%7C__webpack_require___initializeExposesData%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___initializeExposesData.shareScope%3Bfor(const%20key%20in%20__module_federation_bundler_runtime__)%7B__webpack_require__.federation%5Bkey%5D%3D__module_federation_bundler_runtime__%5Bkey%5D%7Dearly(__webpack_require__.federation%2C%22libraryType%22%2C()%3D%3E__module_federation_library_type__)%3Bearly(__webpack_require__.federation%2C%22sharedFallback%22%2C()%3D%3E__module_federation_share_fallbacks__)%3Bconst%20sharedFallback%3D__webpack_require__.federation.sharedFallback%3Bearly(__webpack_require__.federation%2C%22consumesLoadingModuleToHandlerMapping%22%2C()%3D%3E%7Bconst%20consumesLoadingModuleToHandlerMapping%3D%7B%7D%3Bfor(let%5BmoduleId%2Cdata%5Dof%20Object.entries(consumesLoadingModuleToConsumeDataMapping))%7Bvar%20__webpack_require___federation_bundlerRuntime%3BconsumesLoadingModuleToHandlerMapping%5BmoduleId%5D%3D%7Bgetter%3AsharedFallback%3F(__webpack_require___federation_bundlerRuntime%3D__webpack_require__.federation.bundlerRuntime)%3D%3D%3Dnull%7C%7C__webpack_require___federation_bundlerRuntime%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___federation_bundlerRuntime.getSharedFallbackGetter(%7BshareKey%3Adata.shareKey%2Cfactory%3Adata.fallback%2CwebpackRequire%3A__webpack_require__%2ClibraryType%3A__webpack_require__.federation.libraryType%7D)%3Adata.fallback%2CtreeShakingGetter%3AsharedFallback%3Fdata.fallback%3Aundefined%2CshareInfo%3A%7BshareConfig%3A%7BfixedDependencies%3Afalse%2CrequiredVersion%3Adata.requiredVersion%2CstrictVersion%3Adata.strictVersion%2Csingleton%3Adata.singleton%2Ceager%3Adata.eager%7D%2Cscope%3A%5Bdata.shareScope%5D%7D%2CshareKey%3Adata.shareKey%2CtreeShaking%3A__webpack_require__.federation.sharedFallback%3F%7Bget%3Adata.fallback%2Cmode%3Adata.treeShakingMode%7D%3Aundefined%7D%7Dreturn%20consumesLoadingModuleToHandlerMapping%7D)%3Bearly(__webpack_require__.federation%2C%22initOptions%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.initOptions%2C%22name%22%2C()%3D%3E__module_federation_container_name__)%3Bearly(__webpack_require__.federation.initOptions%2C%22shareStrategy%22%2C()%3D%3E__module_federation_share_strategy__)%3Bearly(__webpack_require__.federation.initOptions%2C%22shared%22%2C()%3D%3E%7Bconst%20shared%3D%7B%7D%3Bfor(let%5Bscope%2Cstages%5Dof%20Object.entries(initializeSharingScopeToInitDataMapping))%7Bfor(let%20stage%20of%20stages)%7Bif(typeof%20stage%3D%3D%3D%22object%22%26%26stage!%3D%3Dnull)%7Bconst%7Bname%2Cversion%2Cfactory%2Ceager%2Csingleton%2CrequiredVersion%2CstrictVersion%2CtreeShakingMode%7D%3Dstage%3Bconst%20shareConfig%3D%7B%7D%3Bconst%20isValidValue%3Dfunction(val)%7Breturn%20typeof%20val!%3D%3D%22undefined%22%7D%3Bif(isValidValue(singleton))%7BshareConfig.singleton%3Dsingleton%7Dif(isValidValue(requiredVersion))%7BshareConfig.requiredVersion%3DrequiredVersion%7Dif(isValidValue(eager))%7BshareConfig.eager%3Deager%7Dif(isValidValue(strictVersion))%7BshareConfig.strictVersion%3DstrictVersion%7Dconst%20options%3D%7Bversion%2Cscope%3A%5Bscope%5D%2CshareConfig%2Cget%3Afactory%2CtreeShaking%3AtreeShakingMode%3F%7Bmode%3AtreeShakingMode%7D%3Aundefined%7D%3Bif(shared%5Bname%5D)%7Bshared%5Bname%5D.push(options)%7Delse%7Bshared%5Bname%5D%3D%5Boptions%5D%7D%7D%7D%7Dreturn%20shared%7D)%3Bmerge(__webpack_require__.federation.initOptions%2C%22remotes%22%2C()%3D%3EObject.values(__module_federation_remote_infos__).flat().filter(remote%3D%3Eremote.externalType%3D%3D%3D%22script%22))%3Bmerge(__webpack_require__.federation.initOptions%2C%22plugins%22%2C()%3D%3E__module_federation_runtime_plugins__)%3Bearly(__webpack_require__.federation%2C%22bundlerRuntimeOptions%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions%2C%22remotes%22%2C()%3D%3E(%7B%7D))%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22chunkMapping%22%2C()%3D%3EremotesLoadingChunkMapping)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22remoteInfos%22%2C()%3D%3E__module_federation_remote_infos__)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22idToExternalAndNameMapping%22%2C()%3D%3E%7Bconst%20remotesLoadingIdToExternalAndNameMappingMapping%3D%7B%7D%3Bfor(let%5BmoduleId%2Cdata%5Dof%20Object.entries(remotesLoadingModuleIdToRemoteDataMapping))%7BremotesLoadingIdToExternalAndNameMappingMapping%5BmoduleId%5D%3D%5Bdata.shareScope%2Cdata.name%2Cdata.externalModuleId%2Cdata.remoteName%5D%7Dreturn%20remotesLoadingIdToExternalAndNameMappingMapping%7D)%3Bearly(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22webpackRequire%22%2C()%3D%3E__webpack_require__)%3Bmerge(__webpack_require__.federation.bundlerRuntimeOptions.remotes%2C%22idToRemoteMap%22%2C()%3D%3E%7Bconst%20idToRemoteMap%3D%7B%7D%3Bfor(let%5Bid%2CremoteData%5Dof%20Object.entries(remotesLoadingModuleIdToRemoteDataMapping))%7Bconst%20info%3D__module_federation_remote_infos__%5BremoteData.remoteName%5D%3Bif(info)idToRemoteMap%5Bid%5D%3Dinfo%7Dreturn%20idToRemoteMap%7D)%3Boverride(__webpack_require__%2C%22S%22%2C__webpack_require__.federation.bundlerRuntime.S)%3Bif(__webpack_require__.federation.attachShareScopeMap)%7B__webpack_require__.federation.attachShareScopeMap(__webpack_require__)%7Doverride(__webpack_require__.f%2C%22remotes%22%2C(chunkId%2Cpromises)%3D%3E__webpack_require__.federation.bundlerRuntime.remotes(%7BchunkId%2Cpromises%2CchunkMapping%3AremotesLoadingChunkMapping%2CidToExternalAndNameMapping%3A__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToExternalAndNameMapping%2CidToRemoteMap%3A__webpack_require__.federation.bundlerRuntimeOptions.remotes.idToRemoteMap%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__.f%2C%22consumes%22%2C(chunkId%2Cpromises)%3D%3E__webpack_require__.federation.bundlerRuntime.consumes(%7BchunkId%2Cpromises%2CchunkMapping%3AconsumesLoadingChunkMapping%2CmoduleToHandlerMapping%3A__webpack_require__.federation.consumesLoadingModuleToHandlerMapping%2CinstalledModules%3AconsumesLoadinginstalledModules%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22I%22%2C(name%2CinitScope)%3D%3E__webpack_require__.federation.bundlerRuntime.I(%7BshareScopeName%3Aname%2CinitScope%2CinitPromises%3AinitializeSharingInitPromises%2CinitTokens%3AinitializeSharingInitTokens%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22initContainer%22%2C(shareScope%2CinitScope%2CremoteEntryInitOptions)%3D%3E__webpack_require__.federation.bundlerRuntime.initContainerEntry(%7BshareScope%2CinitScope%2CremoteEntryInitOptions%2CshareScopeKey%3AcontainerShareScope%2CwebpackRequire%3A__webpack_require__%7D))%3Boverride(__webpack_require__%2C%22getContainer%22%2C(module%2CgetScope)%3D%3E%7Bvar%20moduleMap%3D__webpack_require__.initializeExposesData.moduleMap%3B__webpack_require__.R%3DgetScope%3BgetScope%3DObject.prototype.hasOwnProperty.call(moduleMap%2Cmodule)%3FmoduleMap%5Bmodule%5D()%3APromise.resolve().then(()%3D%3E%7Bthrow%20new%20Error('Module%20%22'%2Bmodule%2B'%22%20does%20not%20exist%20in%20container.')%7D)%3B__webpack_require__.R%3Dundefined%3Breturn%20getScope%7D)%3B__webpack_require__.federation.instance%3D__webpack_require__.federation.bundlerRuntime.init(%7BwebpackRequire%3A__webpack_require__%7D)%3Bif((__webpack_require___consumesLoadingData2%3D__webpack_require__.consumesLoadingData)%3D%3D%3Dnull%7C%7C__webpack_require___consumesLoadingData2%3D%3D%3Dvoid%200%3Fvoid%200%3A__webpack_require___consumesLoadingData2.initialConsumes)%7B__webpack_require__.federation.bundlerRuntime.installInitialConsumes(%7BwebpackRequire%3A__webpack_require__%2CinstalledModules%3AconsumesLoadinginstalledModules%2CinitialConsumes%3A__webpack_require__.consumesLoadingData.initialConsumes%2CmoduleToHandlerMapping%3A__webpack_require__.federation.consumesLoadingModuleToHandlerMapping%7D)%7D%7D")
	}
	if (typeof prevStartup === "function") {
		return prevStartup();
	}
	console.warn("[MF] Invalid prevStartup");
};

})();
// module factories are used so entry inlining is disabled
// run startup
var __webpack_exports__ = __webpack_require__.x();
(_JUPYTERLAB = typeof _JUPYTERLAB === 'undefined' ? {} : _JUPYTERLAB)["@jupyterlab/galata-extension"] = __webpack_exports__;
})()
;
//# sourceMappingURL=remoteEntry.4c1216280b8b083c.js.map