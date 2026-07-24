/* Runtime-agnostic wasm32 JIT JS glue for miniexpr.
 * Callers provide runtime facilities via the `runtime` object.
 */
(function(root) {
    'use strict';

    function _meJitInstantiate(runtime, wasmBytes, bridgeLookupFnIdx, emitWarnings) {
        if (!runtime || !wasmBytes) {
            return 0;
        }

        var HEAPF64 = runtime.HEAPF64;
        var HEAPF32 = runtime.HEAPF32;
        var wasmMemory = runtime.wasmMemory;
        var wasmTable = runtime.wasmTable;
        var stackSave = runtime.stackSave;
        var stackAlloc = runtime.stackAlloc;
        var stackRestore = runtime.stackRestore;
        var lengthBytesUTF8 = runtime.lengthBytesUTF8;
        var stringToUTF8 = runtime.stringToUTF8;
        var addFunction = runtime.addFunction;
        var err = runtime.err || function(message) {
            if (typeof console !== 'undefined' && typeof console.error === 'function') {
                console.error(message);
            }
        };
        var meJitEmitWarnings = true;

        if (typeof emitWarnings === 'boolean') {
            meJitEmitWarnings = emitWarnings;
        } else if (typeof runtime.meJitEmitWarnings === 'boolean') {
            meJitEmitWarnings = runtime.meJitEmitWarnings;
        }

        if (!HEAPF64 || !HEAPF32 || !wasmMemory || !wasmTable ||
            typeof stackSave !== 'function' || typeof stackAlloc !== 'function' ||
            typeof stackRestore !== 'function' || typeof lengthBytesUTF8 !== 'function' ||
            typeof stringToUTF8 !== 'function' || typeof addFunction !== 'function') {
            err('[me-wasm-jit] invalid runtime object');
            return 0;
        }

        var src = wasmBytes;
    var enc = new TextEncoder();
    var dec = new TextDecoder();
    /* --- LEB128 helpers ------------------------------------------------- */
    function readULEB(buf, pos) {
        var r = 0, s = 0, b;
        do { b = buf[pos++]; r |= (b & 0x7f) << s; s += 7; } while (b & 0x80);
        return [r, pos];
    }
    function encULEB(v) {
        var a = [];
        do { var b = v & 0x7f; v >>>= 7; if (v) b |= 0x80; a.push(b); } while (v);
        return a;
    }
    function encStr(s) {
        var b = enc.encode(s);
        return encULEB(b.length).concat(Array.from(b));
    }
    function readName(buf, pos) {
        var t = readULEB(buf, pos);
        var n = t[0];
        pos = t[1];
        var s = dec.decode(buf.subarray(pos, pos + n));
        return [s, pos + n];
    }
    function skipLimits(buf, pos) {
        var t = readULEB(buf, pos);
        var flags = t[0];
        pos = t[1];
        t = readULEB(buf, pos);
        pos = t[1];
        if (flags & 0x01) {
            t = readULEB(buf, pos);
            pos = t[1];
        }
        return pos;
    }
    function encMemoryImport() {
        var imp = [];
        imp = imp.concat(encStr("env"), encStr("memory"));
        imp.push(0x02, 0x00); /* memory, limits-flag: no-max */
        imp = imp.concat(encULEB(256));
        return imp;
    }
    function buildImportSecWithMemory() {
        var body = encULEB(1);
        body = body.concat(encMemoryImport());
        var sec = [0x02];
        sec = sec.concat(encULEB(body.length));
        return sec.concat(body);
    }
    function patchImportSec(secData) {
        var pos = 0;
        var t = readULEB(secData, pos);
        var nimports = t[0];
        pos = t[1];
        var entries = [];
        var hasEnvMemory = false;
        for (var i = 0; i < nimports; i++) {
            var start = pos;
            var moduleName = "";
            var fieldName = "";
            t = readName(secData, pos);
            moduleName = t[0];
            pos = t[1];
            t = readName(secData, pos);
            fieldName = t[0];
            pos = t[1];
            var kind = secData[pos++];
            if (kind === 0x00) {
                t = readULEB(secData, pos);
                pos = t[1];
            }
            else if (kind === 0x01) {
                pos++; /* elem type */
                pos = skipLimits(secData, pos);
            }
            else if (kind === 0x02) {
                pos = skipLimits(secData, pos);
                if (moduleName === "env" && fieldName === "memory") {
                    hasEnvMemory = true;
                }
            }
            else if (kind === 0x03) {
                pos += 2; /* valtype + mutability */
            }
            else {
                throw new Error("unsupported wasm import kind " + kind);
            }
            entries.push(Array.from(secData.subarray(start, pos)));
        }
        if (!hasEnvMemory) {
            entries.push(encMemoryImport());
        }
        var body = encULEB(entries.length);
        for (var ei = 0; ei < entries.length; ei++) {
            body = body.concat(entries[ei]);
        }
        var sec = [0x02];
        sec = sec.concat(encULEB(body.length));
        return sec.concat(body);
    }
    function buildEnvImports() {
        var bridgeLookup = null;
        var bridgeCache = Object.create(null);
        if (bridgeLookupFnIdx) {
            bridgeLookup = wasmTable.get(bridgeLookupFnIdx);
        }
        function lookupBridge(name) {
            if (!bridgeLookup) {
                return null;
            }
            if (Object.prototype.hasOwnProperty.call(bridgeCache, name)) {
                return bridgeCache[name];
            }
            var sp = stackSave();
            try {
                var nbytes = lengthBytesUTF8(name) + 1;
                var namePtr = stackAlloc(nbytes);
                stringToUTF8(name, namePtr, nbytes);
                var fnIdx = bridgeLookup(namePtr) | 0;
                bridgeCache[name] = fnIdx ? wasmTable.get(fnIdx) : null;
            } finally {
                stackRestore(sp);
            }
            return bridgeCache[name];
        }
        function bindBridge(name, fallback) {
            var fn = lookupBridge(name);
            return fn ? fn : fallback;
        }
        function fdim(x, y) { return x > y ? (x - y) : 0.0; }
        function copysign(x, y) {
            if (y === 0) {
                return (1 / y === -Infinity) ? -Math.abs(x) : Math.abs(x);
            }
            return y < 0 ? -Math.abs(x) : Math.abs(x);
        }
        function ldexp(x, e) { return x * Math.pow(2.0, e); }
        function rint(x) {
            if (!isFinite(x)) {
                return x;
            }
            var n = Math.round(x);
            if (Math.abs(x - n) === 0.5) {
                n = 2 * Math.round(x / 2);
            }
            return n;
        }
        function remainder(x, y) {
            if (!isFinite(x) || !isFinite(y) || y === 0.0) {
                return NaN;
            }
            return x - y * Math.round(x / y);
        }
        function erfApprox(x) {
            var sign = x < 0 ? -1.0 : 1.0;
            x = Math.abs(x);
            var a1 = 0.254829592;
            var a2 = -0.284496736;
            var a3 = 1.421413741;
            var a4 = -1.453152027;
            var a5 = 1.061405429;
            var p = 0.3275911;
            var t = 1.0 / (1.0 + p * x);
            var y = 1.0 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-x * x);
            return sign * y;
        }
        function erfcApprox(x) { return 1.0 - erfApprox(x); }
        function tgammaApprox(z) {
            var p = [
                676.5203681218851, -1259.1392167224028, 771.32342877765313,
                -176.61502916214059, 12.507343278686905, -0.13857109526572012,
                9.9843695780195716e-6, 1.5056327351493116e-7
            ];
            if (z < 0.5) {
                return Math.PI / (Math.sin(Math.PI * z) * tgammaApprox(1.0 - z));
            }
            z -= 1.0;
            var x = 0.99999999999980993;
            for (var i = 0; i < p.length; i++) {
                x += p[i] / (z + i + 1.0);
            }
            var t = z + p.length - 0.5;
            return Math.sqrt(2.0 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
        }
        function lgammaApprox(x) {
            var g = tgammaApprox(x);
            return Math.log(Math.abs(g));
        }
        function nextafterApprox(x, y) {
            if (isNaN(x) || isNaN(y)) {
                return NaN;
            }
            if (x === y) {
                return y;
            }
            if (x === 0.0) {
                return y > 0.0 ? Number.MIN_VALUE : -Number.MIN_VALUE;
            }
            var buf = new ArrayBuffer(8);
            var dv = new DataView(buf);
            dv.setFloat64(0, x, true);
            var bits = dv.getBigUint64(0, true);
            if ((y > x) === (x > 0.0)) {
                bits += 1n;
            }
            else {
                bits -= 1n;
            }
            dv.setBigUint64(0, bits, true);
            return dv.getFloat64(0, true);
        }
        function meJitExp10(x) { return Math.pow(10.0, x); }
        function meJitSinpi(x) { return Math.sin(Math.PI * x); }
        function meJitCospi(x) { return Math.cos(Math.PI * x); }
        var mathExp2 = Math.exp2 ? Math.exp2 : function(x) { return Math.pow(2.0, x); };
        function meJitLogaddexp(a, b) {
            var hi = a > b ? a : b;
            var lo = a > b ? b : a;
            return hi + Math.log1p(Math.exp(lo - hi));
        }
        function meJitWhere(c, x, y) { return c !== 0.0 ? x : y; }
        function vecUnaryF64(inPtr, outPtr, n, fn) {
            var ii = inPtr >> 3;
            var oo = outPtr >> 3;
            for (var i = 0; i < n; i++) {
                HEAPF64[oo + i] = fn(HEAPF64[ii + i]);
            }
        }
        function vecBinaryF64(aPtr, bPtr, outPtr, n, fn) {
            var aa = aPtr >> 3;
            var bb = bPtr >> 3;
            var oo = outPtr >> 3;
            for (var i = 0; i < n; i++) {
                HEAPF64[oo + i] = fn(HEAPF64[aa + i], HEAPF64[bb + i]);
            }
        }
        function vecUnaryF32(inPtr, outPtr, n, fn) {
            var ii = inPtr >> 2;
            var oo = outPtr >> 2;
            for (var i = 0; i < n; i++) {
                HEAPF32[oo + i] = fn(HEAPF32[ii + i]);
            }
        }
        function vecBinaryF32(aPtr, bPtr, outPtr, n, fn) {
            var aa = aPtr >> 2;
            var bb = bPtr >> 2;
            var oo = outPtr >> 2;
            for (var i = 0; i < n; i++) {
                HEAPF32[oo + i] = fn(HEAPF32[aa + i], HEAPF32[bb + i]);
            }
        }
        var env = {
            memory: wasmMemory,
            acos: Math.acos, acosh: Math.acosh, asin: Math.asin, asinh: Math.asinh,
            atan: Math.atan, atan2: Math.atan2, atanh: Math.atanh, cbrt: Math.cbrt,
            ceil: Math.ceil, copysign: copysign, cos: Math.cos, cosh: Math.cosh,
            erf: erfApprox, erfc: erfcApprox, exp: Math.exp, exp2: mathExp2,
            expm1: Math.expm1, fabs: Math.abs, fdim: fdim, floor: Math.floor,
            fma: function(a, b, c) { return a * b + c; }, fmax: Math.max, fmin: Math.min,
            fmod: function(a, b) { return a % b; }, hypot: Math.hypot, ldexp: ldexp,
            lgamma: lgammaApprox, log: Math.log, log10: Math.log10, log1p: Math.log1p,
            log2: Math.log2, nextafter: nextafterApprox, pow: Math.pow, remainder: remainder,
            rint: rint, round: Math.round, sin: Math.sin, sinh: Math.sinh, sqrt: Math.sqrt,
            tan: Math.tan, tanh: Math.tanh, tgamma: tgammaApprox, trunc: Math.trunc,
            me_jit_exp10: meJitExp10, me_jit_sinpi: meJitSinpi, me_jit_cospi: meJitCospi,
            me_jit_logaddexp: meJitLogaddexp, me_jit_where: meJitWhere
        };
        env.me_wasm32_cast_int = function(x) {
            return x < 0 ? Math.ceil(x) : Math.floor(x);
        };
        env.me_wasm32_cast_float = function(x) {
            return x;
        };
        env.me_wasm32_cast_bool = function(x) {
            return x !== 0 ? 1 : 0;
        };
        env.memset = bindBridge("memset", function(ptr, value, n) {
            if (n > 0) {
                HEAPU8.fill(value & 255, ptr, ptr + n);
            }
            return ptr | 0;
        });
        /* Prefer host wasm bridge symbols; keep JS fallbacks for robustness. */
        env.me_jit_exp10 = bindBridge("me_jit_exp10", env.me_jit_exp10);
        env.me_jit_sinpi = bindBridge("me_jit_sinpi", env.me_jit_sinpi);
        env.me_jit_cospi = bindBridge("me_jit_cospi", env.me_jit_cospi);
        env.me_jit_logaddexp = bindBridge("me_jit_logaddexp", env.me_jit_logaddexp);
        env.me_jit_where = bindBridge("me_jit_where", env.me_jit_where);
        env.me_jit_vec_sin_f64 = bindBridge("me_jit_vec_sin_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.sin); });
        env.me_jit_vec_cos_f64 = bindBridge("me_jit_vec_cos_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.cos); });
        env.me_jit_vec_exp_f64 = bindBridge("me_jit_vec_exp_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.exp); });
        env.me_jit_vec_log_f64 = bindBridge("me_jit_vec_log_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.log); });
        env.me_jit_vec_exp10_f64 = bindBridge("me_jit_vec_exp10_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, meJitExp10); });
        env.me_jit_vec_sinpi_f64 = bindBridge("me_jit_vec_sinpi_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, meJitSinpi); });
        env.me_jit_vec_cospi_f64 = bindBridge("me_jit_vec_cospi_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, meJitCospi); });
        env.me_jit_vec_atan2_f64 = bindBridge("me_jit_vec_atan2_f64", function(aPtr, bPtr, outPtr, n) { vecBinaryF64(aPtr, bPtr, outPtr, n, Math.atan2); });
        env.me_jit_vec_hypot_f64 = bindBridge("me_jit_vec_hypot_f64", function(aPtr, bPtr, outPtr, n) { vecBinaryF64(aPtr, bPtr, outPtr, n, Math.hypot); });
        env.me_jit_vec_pow_f64 = bindBridge("me_jit_vec_pow_f64", function(aPtr, bPtr, outPtr, n) { vecBinaryF64(aPtr, bPtr, outPtr, n, Math.pow); });
        env.me_jit_vec_fmax_f64 = bindBridge("me_jit_vec_fmax_f64", function(aPtr, bPtr, outPtr, n) { vecBinaryF64(aPtr, bPtr, outPtr, n, Math.max); });
        env.me_jit_vec_fmin_f64 = bindBridge("me_jit_vec_fmin_f64", function(aPtr, bPtr, outPtr, n) { vecBinaryF64(aPtr, bPtr, outPtr, n, Math.min); });
        env.me_jit_vec_expm1_f64 = bindBridge("me_jit_vec_expm1_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.expm1); });
        env.me_jit_vec_log10_f64 = bindBridge("me_jit_vec_log10_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.log10); });
        env.me_jit_vec_sinh_f64 = bindBridge("me_jit_vec_sinh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.sinh); });
        env.me_jit_vec_cosh_f64 = bindBridge("me_jit_vec_cosh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.cosh); });
        env.me_jit_vec_tanh_f64 = bindBridge("me_jit_vec_tanh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.tanh); });
        env.me_jit_vec_asinh_f64 = bindBridge("me_jit_vec_asinh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.asinh); });
        env.me_jit_vec_acosh_f64 = bindBridge("me_jit_vec_acosh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.acosh); });
        env.me_jit_vec_atanh_f64 = bindBridge("me_jit_vec_atanh_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.atanh); });
        env.me_jit_vec_abs_f64 = bindBridge("me_jit_vec_abs_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.abs); });
        env.me_jit_vec_sqrt_f64 = bindBridge("me_jit_vec_sqrt_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.sqrt); });
        env.me_jit_vec_log1p_f64 = bindBridge("me_jit_vec_log1p_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.log1p); });
        env.me_jit_vec_exp2_f64 = bindBridge("me_jit_vec_exp2_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, mathExp2); });
        env.me_jit_vec_log2_f64 = bindBridge("me_jit_vec_log2_f64", function(inPtr, outPtr, n) { vecUnaryF64(inPtr, outPtr, n, Math.log2); });
        env.me_jit_vec_sin_f32 = bindBridge("me_jit_vec_sin_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.sin); });
        env.me_jit_vec_cos_f32 = bindBridge("me_jit_vec_cos_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.cos); });
        env.me_jit_vec_exp_f32 = bindBridge("me_jit_vec_exp_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.exp); });
        env.me_jit_vec_log_f32 = bindBridge("me_jit_vec_log_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.log); });
        env.me_jit_vec_exp10_f32 = bindBridge("me_jit_vec_exp10_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, meJitExp10); });
        env.me_jit_vec_sinpi_f32 = bindBridge("me_jit_vec_sinpi_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, meJitSinpi); });
        env.me_jit_vec_cospi_f32 = bindBridge("me_jit_vec_cospi_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, meJitCospi); });
        env.me_jit_vec_atan2_f32 = bindBridge("me_jit_vec_atan2_f32", function(aPtr, bPtr, outPtr, n) { vecBinaryF32(aPtr, bPtr, outPtr, n, Math.atan2); });
        env.me_jit_vec_hypot_f32 = bindBridge("me_jit_vec_hypot_f32", function(aPtr, bPtr, outPtr, n) { vecBinaryF32(aPtr, bPtr, outPtr, n, Math.hypot); });
        env.me_jit_vec_pow_f32 = bindBridge("me_jit_vec_pow_f32", function(aPtr, bPtr, outPtr, n) { vecBinaryF32(aPtr, bPtr, outPtr, n, Math.pow); });
        env.me_jit_vec_fmax_f32 = bindBridge("me_jit_vec_fmax_f32", function(aPtr, bPtr, outPtr, n) { vecBinaryF32(aPtr, bPtr, outPtr, n, Math.max); });
        env.me_jit_vec_fmin_f32 = bindBridge("me_jit_vec_fmin_f32", function(aPtr, bPtr, outPtr, n) { vecBinaryF32(aPtr, bPtr, outPtr, n, Math.min); });
        env.me_jit_vec_expm1_f32 = bindBridge("me_jit_vec_expm1_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.expm1); });
        env.me_jit_vec_log10_f32 = bindBridge("me_jit_vec_log10_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.log10); });
        env.me_jit_vec_sinh_f32 = bindBridge("me_jit_vec_sinh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.sinh); });
        env.me_jit_vec_cosh_f32 = bindBridge("me_jit_vec_cosh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.cosh); });
        env.me_jit_vec_tanh_f32 = bindBridge("me_jit_vec_tanh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.tanh); });
        env.me_jit_vec_asinh_f32 = bindBridge("me_jit_vec_asinh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.asinh); });
        env.me_jit_vec_acosh_f32 = bindBridge("me_jit_vec_acosh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.acosh); });
        env.me_jit_vec_atanh_f32 = bindBridge("me_jit_vec_atanh_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.atanh); });
        env.me_jit_vec_abs_f32 = bindBridge("me_jit_vec_abs_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.abs); });
        env.me_jit_vec_sqrt_f32 = bindBridge("me_jit_vec_sqrt_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.sqrt); });
        env.me_jit_vec_log1p_f32 = bindBridge("me_jit_vec_log1p_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.log1p); });
        env.me_jit_vec_exp2_f32 = bindBridge("me_jit_vec_exp2_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, mathExp2); });
        env.me_jit_vec_log2_f32 = bindBridge("me_jit_vec_log2_f32", function(inPtr, outPtr, n) { vecUnaryF32(inPtr, outPtr, n, Math.log2); });
        return env;
    }
    /* --- parse sections ------------------------------------------------- */
    var pos = 8, sections = [];
    while (pos < src.length) {
        var id = src[pos++];
        var tmp = readULEB(src, pos), len = tmp[0]; pos = tmp[1];
        sections.push({ id: id, data: src.subarray(pos, pos + len) });
        pos += len;
    }
    /* --- reassemble with patched memory -------------------------------- */
    var out = [0x00,0x61,0x73,0x6d, 0x01,0x00,0x00,0x00];
    var impDone = false;
    for (var i = 0; i < sections.length; i++) {
        var s = sections[i];
        if (s.id === 5) continue; /* drop memory section */
        if (s.id === 2) {
            out = out.concat(patchImportSec(s.data));
            impDone = true;
            continue;
        }
        if (!impDone && s.id > 2) {
            out = out.concat(buildImportSecWithMemory());
            impDone = true;
        }
        if (s.id === 7) { /* strip memory export from export section */
            var ep = 0, et = readULEB(s.data, ep), ne = et[0]; ep = et[1];
            var exps = [];
            for (var e = 0; e < ne; e++) {
                var nt = readULEB(s.data, ep), nl = nt[0]; ep = nt[1];
                var nm = dec.decode(s.data.subarray(ep, ep + nl)); ep += nl;
                var kd = s.data[ep++];
                var xt = readULEB(s.data, ep), xi = xt[0]; ep = xt[1];
                if (nm === "memory" && kd === 0x02) continue;
                exps.push({ n: nm, k: kd, i: xi });
            }
            var eb = encULEB(exps.length);
            for (var e = 0; e < exps.length; e++) {
                eb = eb.concat(encStr(exps[e].n));
                eb.push(exps[e].k);
                eb = eb.concat(encULEB(exps[e].i));
            }
            out.push(0x07);
            out = out.concat(encULEB(eb.length));
            out = out.concat(eb);
            continue;
        }
        out.push(s.id);
        out = out.concat(encULEB(s.data.length));
        out = out.concat(Array.from(s.data));
    }
    if (!impDone) {
        out = out.concat(buildImportSecWithMemory());
    }
    /* --- instantiate with shared memory -------------------------------- */
    var patched = new Uint8Array(out);
    try {
        var mod = new WebAssembly.Module(patched);
        var inst = new WebAssembly.Instance(mod, { env: buildEnvImports() });
    } catch (e) {
        if (meJitEmitWarnings) {
            err("[me-wasm-jit] " + e.message);
        }
        return 0;
    }
    var fn = inst.exports["me_dsl_jit_kernel"];
    if (!fn) {
        if (meJitEmitWarnings) {
            err("[me-wasm-jit] missing export");
        }
        return 0;
    }
    return addFunction(fn, "iiii");
    }

    function _meJitFreeFn(runtime, idx) {
        if (!runtime || typeof runtime.removeFunction !== 'function') {
            return;
        }
        if (idx) {
            runtime.removeFunction(idx);
        }
    }

    root._meJitInstantiate = _meJitInstantiate;
    root._meJitFreeFn = _meJitFreeFn;

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = {
            _meJitInstantiate: _meJitInstantiate,
            _meJitFreeFn: _meJitFreeFn
        };
    }
})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this));
