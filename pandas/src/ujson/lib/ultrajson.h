/*
Copyright (c) 2011, Jonas Tarnstrom and ESN Social Software AB
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
   must display the following acknowledgement:
   This product includes software developed by ESN Social Software AB (www.esn.me).
4. Neither the name of the ESN Social Software AB nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY ESN SOCIAL SOFTWARE AB ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ESN SOCIAL SOFTWARE AB BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Portions of code from:
MODP_ASCII - Ascii transformations (upper/lower, etc)
http://code.google.com/p/stringencoders/
Copyright (c) 2007  Nick Galbreath -- nickg [at] modp [dot] com. All rights reserved.

*/

/*
Ultra fast JSON encoder and decoder
Developed by Jonas Tarnstrom (jonas@esn.me).

Encoder notes:
------------------

:: Cyclic references ::
Cyclic referenced objects are not detected.
Set JSONObjectEncoder.recursionMax to suitable value or make sure input object
tree doesn't have cyclic references.

*/

#ifndef __ULTRAJSON_H__
#define __ULTRAJSON_H__

#include <stdio.h>
#include <wchar.h>

//#define JSON_DECODE_NUMERIC_AS_DOUBLE

// Don't output any extra whitespaces when encoding
#define JSON_NO_EXTRA_WHITESPACE

// Max decimals to encode double floating point numbers with
#ifndef JSON_DOUBLE_MAX_DECIMALS
#define JSON_DOUBLE_MAX_DECIMALS 15
#endif

// Max recursion depth, default for encoder
#ifndef JSON_MAX_RECURSION_DEPTH
#define JSON_MAX_RECURSION_DEPTH 1024
#endif

/*
Dictates and limits how much stack space for buffers UltraJSON will use before resorting to provided heap functions */
#ifndef JSON_MAX_STACK_BUFFER_SIZE
#define JSON_MAX_STACK_BUFFER_SIZE 131072
#endif

#ifdef _WIN32

typedef __int64 JSINT64;
typedef unsigned __int64 JSUINT64;

typedef unsigned __int32 uint32_t;
typedef __int32 JSINT32;
typedef uint32_t JSUINT32;
typedef unsigned __int8 JSUINT8;
typedef unsigned __int16 JSUTF16;
typedef unsigned __int32 JSUTF32;
typedef __int64 JSLONG;

#define EXPORTFUNCTION __declspec(dllexport)

#define FASTCALL_MSVC __fastcall
#define FASTCALL_ATTR
#define INLINE_PREFIX __inline

#else

#include <sys/types.h>
typedef int64_t JSINT64;
typedef u_int64_t JSUINT64;

typedef int32_t JSINT32;
typedef u_int32_t JSUINT32;

#define FASTCALL_MSVC
#define FASTCALL_ATTR __attribute__((fastcall))
#define INLINE_PREFIX inline

typedef u_int32_t uint32_t;

typedef u_int8_t JSUINT8;
typedef u_int16_t JSUTF16;
typedef u_int32_t JSUTF32;

typedef int64_t JSLONG;

#define EXPORTFUNCTION
#endif

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define __LITTLE_ENDIAN__
#else

#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define __BIG_ENDIAN__
#endif

#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
#error "Endianess not supported"
#endif

enum JSTYPES
{
    JT_NULL,        // NULL
    JT_TRUE,        //boolean true
    JT_FALSE,       //boolean false
    JT_INT,         //(JSINT32 (signed 32-bit))
    JT_LONG,        //(JSINT64 (signed 64-bit))
    JT_DOUBLE,  //(double)
    JT_UTF8,        //(char 8-bit)
    JT_ARRAY,       // Array structure
    JT_OBJECT,  // Key/Value structure
    JT_INVALID, // Internal, do not return nor expect
};

typedef void * JSOBJ;
typedef void * JSITER;

typedef struct __JSONTypeContext
{
    int type;
    void *encoder;
    void *prv[32];
} JSONTypeContext;

/*
Function pointer declarations, suitable for implementing UltraJSON */
typedef void (*JSPFN_ITERBEGIN)(JSOBJ obj, JSONTypeContext *tc);
typedef int (*JSPFN_ITERNEXT)(JSOBJ obj, JSONTypeContext *tc);
typedef void (*JSPFN_ITEREND)(JSOBJ obj, JSONTypeContext *tc);
typedef JSOBJ (*JSPFN_ITERGETVALUE)(JSOBJ obj, JSONTypeContext *tc);
typedef char *(*JSPFN_ITERGETNAME)(JSOBJ obj, JSONTypeContext *tc, size_t *outLen);
typedef void *(*JSPFN_MALLOC)(size_t size);
typedef void (*JSPFN_FREE)(void *pptr);
typedef void *(*JSPFN_REALLOC)(void *base, size_t size);

typedef struct __JSONObjectEncoder
{
    void (*beginTypeContext)(JSOBJ obj, JSONTypeContext *tc);
    void (*endTypeContext)(JSOBJ obj, JSONTypeContext *tc);
    const char *(*getStringValue)(JSOBJ obj, JSONTypeContext *tc, size_t *_outLen);
    JSINT64 (*getLongValue)(JSOBJ obj, JSONTypeContext *tc);
    JSINT32 (*getIntValue)(JSOBJ obj, JSONTypeContext *tc);
    double (*getDoubleValue)(JSOBJ obj, JSONTypeContext *tc);

    /*
    Begin iteration of an iteratable object (JS_ARRAY or JS_OBJECT)
    Implementor should setup iteration state in ti->prv
    */
    JSPFN_ITERBEGIN iterBegin;

    /*
    Retrieve next object in an iteration. Should return 0 to indicate iteration has reached end or 1 if there are more items.
    Implementor is responsible for keeping state of the iteration. Use ti->prv fields for this
    */
    JSPFN_ITERNEXT iterNext;

    /*
    Ends the iteration of an iteratable object.
    Any iteration state stored in ti->prv can be freed here
    */
    JSPFN_ITEREND iterEnd;

    /*
    Returns a reference to the value object of an iterator
    The is responsible for the life-cycle of the returned string. Use iterNext/iterEnd and ti->prv to keep track of current object
    */
    JSPFN_ITERGETVALUE iterGetValue;

    /*
    Return name of iterator.
    The is responsible for the life-cycle of the returned string. Use iterNext/iterEnd and ti->prv to keep track of current object
    */
    JSPFN_ITERGETNAME iterGetName;

    /*
    Release a value as indicated by setting ti->release = 1 in the previous getValue call.
    The ti->prv array should contain the necessary context to release the value
    */
    void (*releaseObject)(JSOBJ obj);

    /* Library functions
    Set to NULL to use STDLIB malloc,realloc,free */
    JSPFN_MALLOC malloc;
    JSPFN_REALLOC realloc;
    JSPFN_FREE free;

    /*
    Configuration for max recursion, set to 0 to use default (see JSON_MAX_RECURSION_DEPTH)*/
    int recursionMax;

    /*
    Configuration for max decimals of double floating poiunt numbers to encode (0-9) */
    int doublePrecision;

    /*
    If true output will be ASCII with all characters above 127 encoded as \uXXXX. If false output will be UTF-8 or what ever charset strings are brought as */
    int forceASCII;


    /*
    Set to an error message if error occured */
    const char *errorMsg;
    JSOBJ errorObj;

    /* Buffer stuff */
    char *start;
    char *offset;
    char *end;
    int heap;
    int level;

} JSONObjectEncoder;


/*
Encode an object structure into JSON.

Arguments:
obj - An anonymous type representing the object
enc - Function definitions for querying JSOBJ type
buffer - Preallocated buffer to store result in. If NULL function allocates own buffer
cbBuffer - Length of buffer (ignored if buffer is NULL)

Returns:
Encoded JSON object as a null terminated char string.

NOTE:
If the supplied buffer wasn't enough to hold the result the function will allocate a new buffer.
Life cycle of the provided buffer must still be handled by caller.

If the return value doesn't equal the specified buffer caller must release the memory using
JSONObjectEncoder.free or free() as specified when calling this function.
*/
EXPORTFUNCTION char *JSON_EncodeObject(JSOBJ obj, JSONObjectEncoder *enc, char *buffer, size_t cbBuffer);



typedef struct __JSONObjectDecoder
{
    JSOBJ (*newString)(wchar_t *start, wchar_t *end);
    int (*objectAddKey)(JSOBJ obj, JSOBJ name, JSOBJ value);
    int (*arrayAddItem)(JSOBJ obj, JSOBJ value);
    JSOBJ (*newTrue)(void);
    JSOBJ (*newFalse)(void);
    JSOBJ (*newNull)(void);
    JSOBJ (*newObject)(void *decoder);
    JSOBJ (*endObject)(JSOBJ obj);
    JSOBJ (*newArray)(void *decoder);
    JSOBJ (*endArray)(JSOBJ obj);
    JSOBJ (*newInt)(JSINT32 value);
    JSOBJ (*newLong)(JSINT64 value);
    JSOBJ (*newDouble)(double value);
    void (*releaseObject)(JSOBJ obj, void *decoder);
    JSPFN_MALLOC malloc;
    JSPFN_FREE free;
    JSPFN_REALLOC realloc;

    char *errorStr;
    char *errorOffset;



} JSONObjectDecoder;

EXPORTFUNCTION JSOBJ JSON_DecodeObject(JSONObjectDecoder *dec, const char *buffer, size_t cbBuffer);

#endif
