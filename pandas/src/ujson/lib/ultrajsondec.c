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

#include "ultrajson.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <wchar.h>

struct DecoderState
{
	char *start;
	char *end;
	wchar_t *escStart;
	wchar_t *escEnd;
	int escHeap;
	int lastType;
	JSONObjectDecoder *dec;
};

JSOBJ FASTCALL_MSVC decode_any( struct DecoderState *ds) FASTCALL_ATTR;
typedef JSOBJ (*PFN_DECODER)( struct DecoderState *ds);
#define RETURN_JSOBJ_NULLCHECK(_expr) return(_expr);

double createDouble(double intNeg, double intValue, double frcValue, int frcDecimalCount)
{
	static const double g_pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 1000000000000000};

	return (intValue + (frcValue / g_pow10[frcDecimalCount])) * intNeg;
}

static JSOBJ SetError( struct DecoderState *ds, int offset, const char *message)
{
	ds->dec->errorOffset = ds->start + offset;
	ds->dec->errorStr = (char *) message;
	return NULL;
}


FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_numeric ( struct DecoderState *ds)
{
#ifdef JSON_DECODE_NUMERIC_AS_DOUBLE
	double intNeg = 1;
	double intValue;
#else
	int intNeg = 1;
	JSLONG intValue;
#endif

	double expNeg;
	int chr;
	int decimalCount = 0;
	double frcValue = 0.0;
	double expValue;
	char *offset = ds->start;

	if (*(offset) == '-')
	{
		offset ++;
		intNeg = -1;
	}

	// Scan integer part
	intValue = 0;

	while (1)
	{
		chr = (int) (unsigned char) *(offset);

		switch (chr)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			//FIXME: Check for arithemtic overflow here
			//PERF: Don't do 64-bit arithmetic here unless we know we have to
#ifdef JSON_DECODE_NUMERIC_AS_DOUBLE
			intValue = intValue * 10.0 + (double) (chr - 48);
#else
			intValue = intValue * 10LL + (JSLONG) (chr - 48);
#endif
			offset ++;
			break;

		case '.':
			offset ++;
			goto DECODE_FRACTION;
			break;

		case 'e':
		case 'E':
			offset ++;
			goto DECODE_EXPONENT;
			break;

		default:
			goto BREAK_INT_LOOP;
			break;
		}
	}

BREAK_INT_LOOP:

	ds->lastType = JT_INT;
	ds->start = offset;

	//If input string is LONGLONG_MIN here the value is already negative so we should not flip it

#ifdef JSON_DECODE_NUMERIC_AS_DOUBLE
#else
	if (intValue < 0)
	{
		intNeg = 1;
	}
#endif

	//dbg1 = (intValue * intNeg);
	//dbg2 = (JSLONG) dbg1;

#ifdef JSON_DECODE_NUMERIC_AS_DOUBLE
	if (intValue > (double) INT_MAX || intValue < (double) INT_MIN)
#else
	if ( (intValue >> 31))
#endif
	{	
		RETURN_JSOBJ_NULLCHECK(ds->dec->newLong( (JSINT64) (intValue * (JSINT64) intNeg)));
	}
	else
	{
		RETURN_JSOBJ_NULLCHECK(ds->dec->newInt( (JSINT32) (intValue * intNeg)));
	}



DECODE_FRACTION:

	// Scan fraction part
	frcValue = 0.0;
	while (1)
	{
		chr = (int) (unsigned char) *(offset);

		switch (chr)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			if (decimalCount < JSON_DOUBLE_MAX_DECIMALS)
			{
				frcValue = frcValue * 10.0 + (double) (chr - 48);
				decimalCount ++;
			}
			offset ++;
			break;

		case 'e':
		case 'E':
			offset ++;
			goto DECODE_EXPONENT;
			break;

		default:
			goto BREAK_FRC_LOOP;
		}
	}

BREAK_FRC_LOOP:

	if (intValue < 0)
	{
		intNeg = 1;
	}

	//FIXME: Check for arithemtic overflow here
	ds->lastType = JT_DOUBLE;
	ds->start = offset;
	RETURN_JSOBJ_NULLCHECK(ds->dec->newDouble (createDouble( (double) intNeg, (double) intValue, frcValue, decimalCount)));

DECODE_EXPONENT:
	expNeg = 1.0;

	if (*(offset) == '-')
	{
		expNeg = -1.0;
		offset ++;
	}
	else
	if (*(offset) == '+')
	{
		expNeg = +1.0;
		offset ++;
	}

	expValue = 0.0;

	while (1)
	{
		chr = (int) (unsigned char) *(offset);

		switch (chr)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			expValue = expValue * 10.0 + (double) (chr - 48);
			offset ++;
			break;

		default:
			goto BREAK_EXP_LOOP;

		}
	}

BREAK_EXP_LOOP:

#ifdef JSON_DECODE_NUMERIC_AS_DOUBLE
#else
	if (intValue < 0)
	{
		intNeg = 1;
	}
#endif
	
	//FIXME: Check for arithemtic overflow here
	ds->lastType = JT_DOUBLE;
	ds->start = offset;
	RETURN_JSOBJ_NULLCHECK(ds->dec->newDouble (createDouble( (double) intNeg, (double) intValue , frcValue, decimalCount) * pow(10.0, expValue * expNeg)));
}

FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_true ( struct DecoderState *ds) 
{
	char *offset = ds->start;
	offset ++;

	if (*(offset++) != 'r')
		goto SETERROR;
	if (*(offset++) != 'u')
		goto SETERROR;
	if (*(offset++) != 'e')
		goto SETERROR;

	ds->lastType = JT_TRUE;
	ds->start = offset;
	RETURN_JSOBJ_NULLCHECK(ds->dec->newTrue());

SETERROR:
	return SetError(ds, -1, "Unexpected character found when decoding 'true'");
}

FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_false ( struct DecoderState *ds) 
{
	char *offset = ds->start;
	offset ++;

	if (*(offset++) != 'a')
		goto SETERROR;
	if (*(offset++) != 'l')
		goto SETERROR;
	if (*(offset++) != 's')
		goto SETERROR;
	if (*(offset++) != 'e')
		goto SETERROR;

	ds->lastType = JT_FALSE;
	ds->start = offset;
	RETURN_JSOBJ_NULLCHECK(ds->dec->newFalse());

SETERROR:
	return SetError(ds, -1, "Unexpected character found when decoding 'false'");

}


FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_null ( struct DecoderState *ds) 
{
	char *offset = ds->start;
	offset ++;

	if (*(offset++) != 'u')
		goto SETERROR;
	if (*(offset++) != 'l')
		goto SETERROR;
	if (*(offset++) != 'l')
		goto SETERROR;

	ds->lastType = JT_NULL;
	ds->start = offset;
	RETURN_JSOBJ_NULLCHECK(ds->dec->newNull());

SETERROR:
	return SetError(ds, -1, "Unexpected character found when decoding 'null'");
}

FASTCALL_ATTR void FASTCALL_MSVC SkipWhitespace(struct DecoderState *ds) 
{
	char *offset = ds->start;

	while (1)
	{
		switch (*offset)
		{
		case ' ':
		case '\t':
		case '\r':
		case '\n':
			offset ++;
			break;

		default:
			ds->start = offset;
			return;
		}
	}
}


enum DECODESTRINGSTATE
{
	DS_ISNULL = 0x32,
	DS_ISQUOTE,
	DS_ISESCAPE,
	DS_UTFLENERROR,

};

static const JSUINT8 g_decoderLookup[256] = 
{
/* 0x00 */ DS_ISNULL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0x10 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
/* 0x20 */ 1, 1, DS_ISQUOTE, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0x30 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
/* 0x40 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0x50 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, DS_ISESCAPE, 1, 1, 1,
/* 0x60 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0x70 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
/* 0x80 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0x90 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
/* 0xa0 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
/* 0xb0 */ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
/* 0xc0 */ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
/* 0xd0 */ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
/* 0xe0 */ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
/* 0xf0 */ 4, 4, 4, 4, 4, 4, 4, 4, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, DS_UTFLENERROR, 
};


FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_string ( struct DecoderState *ds)
{
	JSUTF16 sur[2] = { 0 };
	int iSur = 0;
	int index;
	wchar_t *escOffset;
	size_t escLen = (ds->escEnd - ds->escStart);
	JSUINT8 *inputOffset;
	JSUINT8 oct;
	JSUTF32 ucs;
	ds->lastType = JT_INVALID;
	ds->start ++;

	if ( (ds->end - ds->start) > escLen)
	{
		size_t newSize = (ds->end - ds->start);

		if (ds->escHeap)
		{
			ds->escStart = (wchar_t *) ds->dec->realloc (ds->escStart, newSize * sizeof(wchar_t));
		}
		else
		{
			wchar_t *oldStart = ds->escStart;
			ds->escHeap = 1;
			ds->escStart = (wchar_t *) ds->dec->malloc (newSize * sizeof(wchar_t));
			memcpy (ds->escStart, oldStart, escLen * sizeof(wchar_t));
		}

		ds->escEnd = ds->escStart + newSize;
	}

	escOffset = ds->escStart;
	inputOffset = ds->start;

	while(1)
	{
		switch (g_decoderLookup[(JSUINT8)(*inputOffset)])
		{
		case DS_ISNULL:
			return SetError(ds, -1, "Unmatched ''\"' when when decoding 'string'");

		case DS_ISQUOTE:
			ds->lastType = JT_UTF8;
			inputOffset ++;
			ds->start += ( (char *) inputOffset - (ds->start));
			RETURN_JSOBJ_NULLCHECK(ds->dec->newString(ds->escStart, escOffset));

		case DS_UTFLENERROR:
			return SetError (ds, -1, "Invalid UTF-8 sequence length when decoding 'string'");

		case DS_ISESCAPE:
			inputOffset ++;
			switch (*inputOffset)
			{
			case '\\': *(escOffset++) = L'\\'; inputOffset++; continue;
			case '\"': *(escOffset++) = L'\"'; inputOffset++; continue;
			case '/':  *(escOffset++) = L'/';  inputOffset++; continue;
			case 'b':  *(escOffset++) = L'\b'; inputOffset++; continue;
			case 'f':  *(escOffset++) = L'\f'; inputOffset++; continue;
			case 'n':  *(escOffset++) = L'\n'; inputOffset++; continue;
			case 'r':  *(escOffset++) = L'\r'; inputOffset++; continue;
			case 't':  *(escOffset++) = L'\t'; inputOffset++; continue;

			case 'u':
				{
					int index;
					inputOffset ++;

					for (index = 0; index < 4; index ++)
					{
						switch (*inputOffset)
						{
						case '\0':	return SetError (ds, -1, "Unterminated unicode escape sequence when decoding 'string'");
						default:		return SetError (ds, -1, "Unexpected character in unicode escape sequence when decoding 'string'");

						case '0':
						case '1':
						case '2':
						case '3':
						case '4':
						case '5':
						case '6':
						case '7':
						case '8':
						case '9':
							sur[iSur] = (sur[iSur] << 4) + (JSUTF16) (*inputOffset - '0');
							break;

						case 'a':
						case 'b':
						case 'c':
						case 'd':
						case 'e':
						case 'f':
							sur[iSur] = (sur[iSur] << 4) + 10 + (JSUTF16) (*inputOffset - 'a');
							break;

						case 'A':
						case 'B':
						case 'C':
						case 'D':
						case 'E':
						case 'F':
							sur[iSur] = (sur[iSur] << 4) + 10 + (JSUTF16) (*inputOffset - 'A');
							break;
						}

						inputOffset ++;
					}


					if (iSur == 0)
					{
						if((sur[iSur] & 0xfc00) == 0xd800)
						{
							// First of a surrogate pair, continue parsing
							iSur ++;
							break;
						} 
						(*escOffset++) = (wchar_t) sur[iSur];
						iSur = 0;
					}
					else
					{
						// Decode pair
						if ((sur[1] & 0xfc00) != 0xdc00)
						{
							return SetError (ds, -1, "Unpaired high surrogate when decoding 'string'");
						}

#if WCHAR_MAX == 0xffff
						(*escOffset++) = (wchar_t) sur[0];
						(*escOffset++) = (wchar_t) sur[1];
#else
						(*escOffset++) = (wchar_t) 0x10000 + (((sur[0] - 0xd800) << 10) | (sur[1] - 0xdc00));
#endif
						iSur = 0;
					}
					break;
				}

			case '\0': return SetError(ds, -1, "Unterminated escape sequence when decoding 'string'");
			default: return SetError(ds, -1, "Unrecognized escape sequence when decoding 'string'");
			}
			break;

		case 1:
			*(escOffset++) = (wchar_t) (*inputOffset++); 
			break;

		case 2:
		{
			ucs = (*inputOffset++) & 0x1f;
			ucs <<= 6;
			if (((*inputOffset) & 0x80) != 0x80)
			{
				return SetError(ds, -1, "Invalid octet in UTF-8 sequence when decoding 'string'");
			}
			ucs |= (*inputOffset++) & 0x3f;
			if (ucs < 0x80)	return SetError (ds, -1, "Overlong 2 byte UTF-8 sequence detected when decoding 'string'");
			*(escOffset++) = (wchar_t) ucs;
			break;
		}

		case 3:
		{
			JSUTF32 ucs = 0;
			ucs |= (*inputOffset++) & 0x0f;

			for (index = 0; index < 2; index ++)
			{
				ucs <<= 6;
				oct = (*inputOffset++);

				if ((oct & 0x80) != 0x80)
				{
					return SetError(ds, -1, "Invalid octet in UTF-8 sequence when decoding 'string'");
				}

				ucs |= oct & 0x3f;
			}

			if (ucs < 0x800) return SetError (ds, -1, "Overlong 3 byte UTF-8 sequence detected when encoding string");
			*(escOffset++) = (wchar_t) ucs;
			break;
		}

		case 4:
		{
			JSUTF32 ucs = 0;
			ucs |= (*inputOffset++) & 0x07;

			for (index = 0; index < 3; index ++)
			{
				ucs <<= 6;
				oct = (*inputOffset++);

				if ((oct & 0x80) != 0x80)
				{
					return SetError(ds, -1, "Invalid octet in UTF-8 sequence when decoding 'string'");
				}

				ucs |= oct & 0x3f;
			}

			if (ucs < 0x10000) return SetError (ds, -1, "Overlong 4 byte UTF-8 sequence detected when decoding 'string'");

			#if WCHAR_MAX == 0xffff
			if (ucs >= 0x10000)
			{
				ucs -= 0x10000;
				*(escOffset++) = (ucs >> 10) + 0xd800;
				*(escOffset++) = (ucs & 0x3ff) + 0xdc00;
			}
			else
			{
				*(escOffset++) = (wchar_t) ucs;
			}
			#else
			*(escOffset++) = (wchar_t) ucs;
			#endif
			break;
		}
		}
	}
}

FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_array( struct DecoderState *ds)
{
	JSOBJ itemValue;
	JSOBJ newObj = ds->dec->newArray(ds->dec);

	ds->lastType = JT_INVALID;
	ds->start ++;

	while (1)//(*ds->start) != '\0')
	{
		SkipWhitespace(ds);

		if ((*ds->start) == ']')
		{
			ds->start++;
			return ds->dec->endArray(newObj);
		}

		itemValue = decode_any(ds);

		if (itemValue == NULL)
		{
			ds->dec->releaseObject(newObj, ds->dec);
			return NULL;
		}

		if (!ds->dec->arrayAddItem (newObj, itemValue)) 
		{
			ds->dec->releaseObject(newObj, ds->dec);
			return NULL;
		}

		SkipWhitespace(ds);

		switch (*(ds->start++))
		{
			case ']':
				return ds->dec->endArray(newObj);

			case ',':
				break;

			default:
				ds->dec->releaseObject(newObj, ds->dec);
				return SetError(ds, -1, "Unexpected character in found when decoding array value");
		}
	}

	ds->dec->releaseObject(newObj, ds->dec);
	return SetError(ds, -1, "Unmatched ']' when decoding 'array'");
}



FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_object( struct DecoderState *ds)
{
	JSOBJ itemName;
	JSOBJ itemValue;
	JSOBJ newObj = ds->dec->newObject(ds->dec);

	ds->start ++;

	while (1)
	{
		SkipWhitespace(ds);

		if ((*ds->start) == '}')
		{
			ds->start ++;
			return ds->dec->endObject(newObj);
		}

		ds->lastType = JT_INVALID;
		itemName = decode_any(ds);

		if (itemName == NULL)
		{
			ds->dec->releaseObject(newObj, ds->dec);
			return NULL;
		}

		if (ds->lastType != JT_UTF8)
		{
			ds->dec->releaseObject(newObj, ds->dec);
			ds->dec->releaseObject(itemName, ds->dec);
			return SetError(ds, -1, "Key name of object must be 'string' when decoding 'object'");
		}

		SkipWhitespace(ds);

		if (*(ds->start++) != ':')
		{
			ds->dec->releaseObject(newObj, ds->dec);
			ds->dec->releaseObject(itemName, ds->dec);
			return SetError(ds, -1, "No ':' found when decoding object value");
		}

		SkipWhitespace(ds);

		itemValue = decode_any(ds);

		if (itemValue == NULL)
		{
			ds->dec->releaseObject(newObj, ds->dec);
			ds->dec->releaseObject(itemName, ds->dec);
			return NULL;
		}

		if (!ds->dec->objectAddKey (newObj, itemName, itemValue)) 
		{
			ds->dec->releaseObject(newObj, ds->dec);
			ds->dec->releaseObject(itemName, ds->dec);
			ds->dec->releaseObject(itemValue, ds->dec);
			return NULL;
		}

		SkipWhitespace(ds);

		switch (*(ds->start++))
		{
			case '}':
				return ds->dec->endObject(newObj);

			case ',':
				break;

			default:
				ds->dec->releaseObject(newObj, ds->dec);
				return SetError(ds, -1, "Unexpected character in found when decoding object value");
		}
	}

	ds->dec->releaseObject(newObj, ds->dec);
	return SetError(ds, -1, "Unmatched '}' when decoding object");
}

FASTCALL_ATTR JSOBJ FASTCALL_MSVC decode_any(struct DecoderState *ds)
{
	while (1)
	{
		switch (*ds->start)
		{
			case '\"': 
				return decode_string (ds);
			case '0': 
			case '1':
			case '2': 
			case '3': 
			case '4': 
			case '5':
			case '6': 
			case '7': 
			case '8': 
			case '9': 
			case '-': 
				return decode_numeric (ds);

			case '[':	return decode_array (ds);
			case '{': return decode_object (ds);
			case 't': return decode_true (ds);
			case 'f': return decode_false (ds);
			case 'n': return decode_null (ds);

			case ' ':
			case '\t':
			case '\r':
			case '\n':
				// White space
				ds->start ++;
				break;

			default:
				return SetError(ds, -1, "Expected object or value");
		}
	}
}


JSOBJ JSON_DecodeObject(JSONObjectDecoder *dec, const char *buffer, size_t cbBuffer)
{

	/*
	FIXME: Base the size of escBuffer of that of cbBuffer so that the unicode escaping doesn't run into the wall each time */
	struct DecoderState ds;
	wchar_t escBuffer[(JSON_MAX_STACK_BUFFER_SIZE / sizeof(wchar_t))];
	JSOBJ ret;
	
	ds.start = (char *) buffer;
	ds.end = ds.start + cbBuffer;

	ds.escStart = escBuffer;
	ds.escEnd = ds.escStart + (JSON_MAX_STACK_BUFFER_SIZE / sizeof(wchar_t));
	ds.escHeap = 0;
	ds.dec = dec;
	ds.dec->errorStr = NULL;
	ds.dec->errorOffset = NULL;

	ds.dec = dec;

	ret = decode_any (&ds);
	
	if (ds.escHeap)
	{
		dec->free(ds.escStart);
	}
	return ret;
}
