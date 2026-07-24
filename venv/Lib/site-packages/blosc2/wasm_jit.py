#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import os
from pathlib import Path

_HELPERS_REGISTERED = False

_REGISTER_HELPERS_JS = r"""
(() => {
  const g = globalThis;
  if (g.__blosc2_me_jit_helper_ptrs) {
    return g.__blosc2_me_jit_helper_ptrs;
  }

  const candidates = [];
  const addCandidate = (name, obj) => {
    if (!obj || (typeof obj !== "object" && typeof obj !== "function")) {
      return;
    }
    candidates.push({ name, obj });
  };
  const addDerivedCandidates = (baseName, obj) => {
    if (!obj || (typeof obj !== "object" && typeof obj !== "function")) {
      return;
    }
    addCandidate(`${baseName}._module`, obj._module);
    addCandidate(`${baseName}.module`, obj.module);
    addCandidate(`${baseName}.Module`, obj.Module);
    addCandidate(`${baseName}.asm`, obj.asm);
    addCandidate(`${baseName}.wasmExports`, obj.wasmExports);
    addCandidate(`${baseName}.wasm`, obj.wasm);
    addCandidate(`${baseName}.__wasm`, obj.__wasm);
    addCandidate(`${baseName}.pyodide`, obj.pyodide);
    addCandidate(`${baseName}._api`, obj._api);
  };

  addCandidate("globalThis", g);
  addCandidate("globalThis.Module", g.Module);
  addCandidate("globalThis.__blosc2_pyodide_module", g.__blosc2_pyodide_module);
  addCandidate("globalThis.__blosc2_pyodide_api", g.__blosc2_pyodide_api);
  addCandidate("globalThis.pyodide", g.pyodide);
  addCandidate("globalThis.pyodide._module", g.pyodide && g.pyodide._module);
  addCandidate("globalThis.pyodide.module", g.pyodide && g.pyodide.module);
  addCandidate("globalThis.pyodide.Module", g.pyodide && g.pyodide.Module);
  addCandidate("globalThis.pyodide._api", g.pyodide && g.pyodide._api);
  addCandidate("globalThis.pyodide._api._module", g.pyodide && g.pyodide._api && g.pyodide._api._module);
  addCandidate("globalThis.pyodide._api.Module", g.pyodide && g.pyodide._api && g.pyodide._api.Module);
  addDerivedCandidates("globalThis", g);
  addDerivedCandidates("globalThis.pyodide", g.pyodide);
  addDerivedCandidates("globalThis.__blosc2_pyodide_module", g.__blosc2_pyodide_module);
  addDerivedCandidates("globalThis.__blosc2_pyodide_api", g.__blosc2_pyodide_api);

  const resolve = (name) => {
    for (const cand of candidates) {
      let value;
      try {
        value = cand.obj[name];
      } catch (_e) {
        value = undefined;
      }
      if (value !== undefined && value !== null) {
        if (typeof value === "function") {
          return value.bind(cand.obj);
        }
        return value;
      }
    }
    if (g[name] !== undefined && g[name] !== null) {
      return g[name];
    }
    return null;
  };

  const wasmExports = resolve("wasmExports") || resolve("exports");
  const asmObj = resolve("asm");

  const isWasmMemory = (value) =>
    typeof WebAssembly !== "undefined" &&
    typeof WebAssembly.Memory !== "undefined" &&
    value instanceof WebAssembly.Memory;
  const isWasmTable = (value) =>
    typeof WebAssembly !== "undefined" &&
    typeof WebAssembly.Table !== "undefined" &&
    value instanceof WebAssembly.Table;
  const heapU8ForProbe = resolve("HEAPU8");
  const heapBufferForProbe = heapU8ForProbe && heapU8ForProbe.buffer ? heapU8ForProbe.buffer : null;
  const heapBufferLenForProbe =
    heapBufferForProbe && typeof heapBufferForProbe.byteLength === "number"
      ? heapBufferForProbe.byteLength
      : -1;

  const isMemoryLike = (value) => {
    if (!value) {
      return false;
    }
    if (isWasmMemory(value)) {
      return true;
    }
    let buf = null;
    try {
      buf = value.buffer;
    } catch (_e) {
      buf = null;
    }
    if (!buf || typeof buf.byteLength !== "number") {
      return false;
    }
    if (typeof value.grow !== "function") {
      return false;
    }
    if (heapBufferForProbe && buf !== heapBufferForProbe) {
      const bufLen = typeof buf.byteLength === "number" ? buf.byteLength : -1;
      if (heapBufferLenForProbe > 0 && bufLen > 0 && bufLen < heapBufferLenForProbe) {
        return false;
      }
    }
    return true;
  };

  const isTableLike = (value) => {
    if (!value) {
      return false;
    }
    if (isWasmTable(value)) {
      return true;
    }
    return (
      typeof value.get === "function" &&
      typeof value.grow === "function" &&
      typeof value.length === "number"
    );
  };

  const findMemoryOrTableByType = (wantMemory) => {
    const isObj = (v) => v && (typeof v === "object" || typeof v === "function");
    const seen = new Set();
    const queue = [];
    const maxDepth = 6;
    const maxVisited = 5000;

    for (const cand of candidates) {
      if (isObj(cand.obj)) {
        queue.push({ value: cand.obj, depth: 0 });
      }
    }

    while (queue.length > 0 && seen.size < maxVisited) {
      const node = queue.shift();
      const obj = node.value;
      const depth = node.depth;
      if (!isObj(obj) || seen.has(obj)) {
        continue;
      }
      seen.add(obj);

      if (wantMemory && isMemoryLike(obj)) {
        return obj;
      }
      if (!wantMemory && isTableLike(obj)) {
        return obj;
      }
      if (depth >= maxDepth) {
        continue;
      }

      let keys = [];
      try {
        keys = Object.getOwnPropertyNames(obj);
      } catch (_e) {
        keys = [];
      }
      let symKeys = [];
      try {
        symKeys = Object.getOwnPropertySymbols(obj);
      } catch (_e) {
        symKeys = [];
      }
      const allKeys = keys.concat(symKeys);

      for (const key of allKeys) {
        let value;
        try {
          value = obj[key];
        } catch (_e) {
          continue;
        }

        if (wantMemory && isMemoryLike(value)) {
          return value;
        }
        if (!wantMemory && isTableLike(value)) {
          return value;
        }
        if (isObj(value)) {
          if (wantMemory && isMemoryLike(value.memory)) {
            return value.memory;
          }
          if (!wantMemory && isTableLike(value.__indirect_function_table)) {
            return value.__indirect_function_table;
          }
          queue.push({ value, depth: depth + 1 });
        }
      }

      let proto = null;
      try {
        proto = Object.getPrototypeOf(obj);
      } catch (_e) {
        proto = null;
      }
      if (isObj(proto)) {
        queue.push({ value: proto, depth: depth + 1 });
      }
    }

    return null;
  };

  const captureMemoryViaGrowHook = () => {
    if (
      typeof WebAssembly === "undefined" ||
      typeof WebAssembly.Memory === "undefined" ||
      !WebAssembly.Memory.prototype ||
      typeof WebAssembly.Memory.prototype.grow !== "function"
    ) {
      return null;
    }

    const growMemory = resolve("growMemory");
    const resizeHeap = resolve("_emscripten_resize_heap");
    if (typeof growMemory !== "function" && typeof resizeHeap !== "function") {
      return null;
    }

    const heapU8 = resolve("HEAPU8");
    const currentBytes =
      heapU8 && heapU8.buffer && typeof heapU8.buffer.byteLength === "number"
        ? heapU8.buffer.byteLength
        : 0;
    if (currentBytes <= 0) {
      return null;
    }
    const onePage = 64 * 1024;
    let targetBytes = currentBytes + onePage;
    const getHeapMax = resolve("getHeapMax");
    if (typeof getHeapMax === "function") {
      try {
        const maxBytes = getHeapMax();
        if (typeof maxBytes === "number" && maxBytes > 0) {
          targetBytes = Math.min(targetBytes, maxBytes);
        }
      } catch (_e) {
        /* ignore */
      }
    }
    if (targetBytes <= currentBytes) {
      return null;
    }

    let captured = null;
    const originalGrow = WebAssembly.Memory.prototype.grow;
    WebAssembly.Memory.prototype.grow = function patchedGrow(pages) {
      captured = this;
      return originalGrow.call(this, pages);
    };

    try {
      if (typeof growMemory === "function") {
        growMemory(targetBytes);
      } else if (typeof resizeHeap === "function") {
        resizeHeap(targetBytes);
      }
    } catch (_e) {
      /* best effort only */
    } finally {
      WebAssembly.Memory.prototype.grow = originalGrow;
    }

    if (captured && isMemoryLike(captured)) {
      return captured;
    }
    return null;
  };

  const deriveRuntimeFromAdjustedImports = () => {
    for (const cand of candidates) {
      const obj = cand.obj;
      if (!obj || typeof obj.adjustWasmImports !== "function") {
        continue;
      }
      try {
        const importsObj = { env: {} };
        const adjustedMaybe = obj.adjustWasmImports(importsObj);
        const adjusted =
          adjustedMaybe && (typeof adjustedMaybe === "object" || typeof adjustedMaybe === "function")
            ? adjustedMaybe
            : importsObj;
        const env =
          (adjusted && adjusted.env) ||
          (importsObj && importsObj.env) ||
          null;
        if (!env) {
          continue;
        }
        const mem =
          env.memory ||
          env.wasmMemory ||
          (adjusted && (adjusted.memory || adjusted.wasmMemory)) ||
          null;
        const tbl =
          env.__indirect_function_table ||
          env.wasmTable ||
          (adjusted && (adjusted.__indirect_function_table || adjusted.wasmTable)) ||
          null;
        if (mem || tbl) {
          return { memory: mem, table: tbl };
        }
      } catch (_e) {
        continue;
      }
    }
    return null;
  };

  const adjustedRuntime = deriveRuntimeFromAdjustedImports();

  const wasmMemory =
    resolve("wasmMemory") ||
    resolve("memory") ||
    resolve("wasmMemoryObject") ||
    resolve("__wasmMemory") ||
    (asmObj && asmObj.memory) ||
    (asmObj && asmObj.wasmMemory) ||
    (wasmExports && wasmExports.memory) ||
    (adjustedRuntime && adjustedRuntime.memory) ||
    captureMemoryViaGrowHook() ||
    findMemoryOrTableByType(true) ||
    null;
  const wasmTable =
    resolve("wasmTable") ||
    resolve("__indirect_function_table") ||
    (asmObj && asmObj.__indirect_function_table) ||
    (asmObj && asmObj.wasmTable) ||
    (wasmExports && wasmExports.__indirect_function_table) ||
    (adjustedRuntime && adjustedRuntime.table) ||
    findMemoryOrTableByType(false) ||
    null;
  const runtime = {
    HEAPF32: resolve("HEAPF32"),
    HEAPF64: resolve("HEAPF64"),
    HEAPU8: heapU8ForProbe,
    wasmMemory,
    wasmTable,
    addFunction: resolve("addFunction"),
    removeFunction: resolve("removeFunction"),
    stackSave: resolve("stackSave"),
    stackAlloc: resolve("stackAlloc"),
    stackRestore: resolve("stackRestore"),
    lengthBytesUTF8: resolve("lengthBytesUTF8"),
    stringToUTF8: resolve("stringToUTF8"),
    err: resolve("err"),
  };

  const required = [
    "HEAPF32",
    "HEAPF64",
    "HEAPU8",
    "wasmMemory",
    "wasmTable",
    "addFunction",
    "removeFunction",
    "stackSave",
    "stackAlloc",
    "stackRestore",
    "lengthBytesUTF8",
    "stringToUTF8",
  ];
  const missing = required.filter((name) => !runtime[name]);
  if (missing.length > 0) {
    const aliasKeys = [
      "wasmMemory",
      "memory",
      "wasmExports",
      "asm",
      "__indirect_function_table",
      "wasmTable",
      "adjustWasmImports",
    ];
    const keyRegex = /(mem|wasm|asm|module|heap)/i;
    const diag = candidates.map((cand) => {
      const have = required.filter((name) => {
        try {
          return !!cand.obj[name];
        } catch (_e) {
          return false;
        }
      });
      const aliases = aliasKeys.filter((name) => {
        try {
          return cand.obj[name] !== undefined && cand.obj[name] !== null;
        } catch (_e) {
          return false;
        }
      });
      let ownKeys = [];
      try {
        ownKeys = Object.getOwnPropertyNames(cand.obj);
      } catch (_e) {
        ownKeys = [];
      }
      const interesting = ownKeys.filter((k) => keyRegex.test(k)).slice(0, 20);
      return `${cand.name}=[${have.join(",")}],aliases=[${aliases.join(",")}],keys=[${interesting.join(",")}]`;
    }).join(" | ");
    return {
      instantiatePtr: 0,
      freePtr: 0,
      error: `missing runtime members: ${missing.join(", ")}; candidates: ${diag}`,
    };
  }

  if (typeof g._meJitInstantiate !== "function" || typeof g._meJitFreeFn !== "function") {
    return { instantiatePtr: 0, freePtr: 0, error: "me_jit_glue exports unavailable" };
  }

  const refreshRuntimeViews = () => {
    const updater = resolve("updateMemoryViews");
    if (typeof updater === "function") {
      try {
        updater();
      } catch (_e) {
        /* best effort only */
      }
      runtime.HEAPU8 = resolve("HEAPU8") || runtime.HEAPU8;
      runtime.HEAPF32 = resolve("HEAPF32") || runtime.HEAPF32;
      runtime.HEAPF64 = resolve("HEAPF64") || runtime.HEAPF64;
    }

    const mem = runtime.wasmMemory;
    const buffer = mem && mem.buffer ? mem.buffer : null;
    if (!buffer || typeof buffer.byteLength !== "number" || buffer.byteLength === 0) {
      return null;
    }

    const heapU8 = runtime.HEAPU8;
    if (!heapU8 || heapU8.buffer !== buffer || heapU8.byteLength === 0) {
      runtime.HEAPU8 = new Uint8Array(buffer);
    }
    const heapF32 = runtime.HEAPF32;
    if (!heapF32 || heapF32.buffer !== buffer || heapF32.byteLength === 0) {
      runtime.HEAPF32 = new Float32Array(buffer);
    }
    const heapF64 = runtime.HEAPF64;
    if (!heapF64 || heapF64.buffer !== buffer || heapF64.byteLength === 0) {
      runtime.HEAPF64 = new Float64Array(buffer);
    }

    return runtime.HEAPU8;
  };

  const instantiateWrapper = (wasmPtr, wasmLen, bridgeLookupFnIdx) => {
    const start = wasmPtr >>> 0;
    const len = wasmLen >>> 0;
    if (start === 0 || len === 0) {
      return 0;
    }
    const heapU8 = refreshRuntimeViews();
    if (!heapU8) {
      return 0;
    }
    const end = (start + len) >>> 0;
    if (end > heapU8.byteLength || end < start) {
      return 0;
    }
    const wasmBytes = new Uint8Array(len);
    wasmBytes.set(heapU8.subarray(start, end));
    return g._meJitInstantiate(runtime, wasmBytes, bridgeLookupFnIdx | 0) | 0;
  };
  const freeWrapper = (fnIdx) => {
    g._meJitFreeFn(runtime, fnIdx | 0);
  };

  const instantiatePtr = runtime.addFunction(instantiateWrapper, "iiii");
  const freePtr = runtime.addFunction(freeWrapper, "vi");
  g.__blosc2_me_jit_helper_ptrs = {
    instantiatePtr,
    freePtr,
    instantiateWrapper,
    freeWrapper,
    runtime,
  };
  return g.__blosc2_me_jit_helper_ptrs;
})()
"""


def _trace_enabled() -> bool:
    value = os.environ.get("ME_DSL_TRACE", "")
    return value.lower() in {"1", "true", "on", "yes"}


def _trace(message: str) -> None:
    if _trace_enabled():
        print(f"[blosc2.wasm-jit] {message}")


def _js_eval(js_mod, source: str):
    evaluator = getattr(js_mod, "eval", None)
    if evaluator is not None:
        return evaluator(source)
    return js_mod.globalThis.eval(source)


def _load_glue_once(js_mod) -> bool:
    has_exports = _js_eval(
        js_mod,
        "typeof globalThis._meJitInstantiate === 'function' && "
        "typeof globalThis._meJitFreeFn === 'function'",
    )
    if bool(has_exports):
        return True

    glue_path = Path(__file__).with_name("me_jit_glue.js")
    try:
        glue_source = glue_path.read_text(encoding="utf-8")
    except OSError as exc:
        _trace(f"could not read {glue_path.name}: {exc}")
        return False

    try:
        _js_eval(js_mod, glue_source)
    except Exception as exc:  # pragma: no cover - pyodide-specific error path
        _trace(f"failed to evaluate {glue_path.name}: {exc}")
        return False

    has_exports = _js_eval(
        js_mod,
        "typeof globalThis._meJitInstantiate === 'function' && "
        "typeof globalThis._meJitFreeFn === 'function'",
    )
    return bool(has_exports)


def _inject_pyodide_runtime_handles(js_mod) -> None:
    try:
        import pyodide_js
    except ImportError:
        return

    module_obj = None
    for name in ("_module", "module", "Module"):
        module_obj = getattr(pyodide_js, name, None)
        if module_obj is not None:
            break
    if module_obj is not None:
        js_mod.globalThis.__blosc2_pyodide_module = module_obj
        _trace("captured pyodide_js module handle")

    api_obj = getattr(pyodide_js, "_api", None)
    if api_obj is not None:
        js_mod.globalThis.__blosc2_pyodide_api = api_obj
        _trace("captured pyodide_js API handle")


def _create_helper_ptrs(js_mod) -> tuple[int, int] | None:
    try:
        result = _js_eval(js_mod, _REGISTER_HELPERS_JS)
    except Exception as exc:  # pragma: no cover - pyodide-specific error path
        _trace(f"helper setup JS failed: {exc}")
        return None

    try:
        instantiate_ptr = int(result.instantiatePtr)
        free_ptr = int(result.freePtr)
    except Exception as exc:  # pragma: no cover - pyodide-specific error path
        _trace(f"unexpected helper setup result: {exc}")
        return None

    if instantiate_ptr == 0 or free_ptr == 0:
        with_error = getattr(result, "error", None)
        if with_error:
            _trace(str(with_error))
        return None
    return instantiate_ptr, free_ptr


def init_wasm_jit_helpers() -> bool:
    global _HELPERS_REGISTERED
    if _HELPERS_REGISTERED:
        return True

    try:
        import js
    except ImportError:
        return False

    from . import blosc2_ext

    if not hasattr(blosc2_ext, "register_wasm_jit_helpers"):
        _trace("extension does not expose register_wasm_jit_helpers")
        return False

    _inject_pyodide_runtime_handles(js)
    if not _load_glue_once(js):
        _trace("me_jit_glue.js was not loaded")
        return False

    helper_ptrs = _create_helper_ptrs(js)
    if helper_ptrs is None:
        _trace("could not allocate addFunction helper pointers")
        return False

    instantiate_ptr, free_ptr = helper_ptrs
    try:
        blosc2_ext.register_wasm_jit_helpers(instantiate_ptr, free_ptr)
    except Exception as exc:  # pragma: no cover - pyodide-specific error path
        _trace(f"C helper registration failed: {exc}")
        return False
    _HELPERS_REGISTERED = True
    _trace(f"registered wasm JIT helper pointers instantiate={instantiate_ptr} free={free_ptr}")
    return True
