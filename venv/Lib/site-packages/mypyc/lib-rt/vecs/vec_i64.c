#define VEC VecI64
#define VEC_TYPE VecI64Type
#define VEC_OBJECT VecI64Object
#define BUF_OBJECT VecI64BufObject
#define BUF_TYPE VecI64BufType
#define NAME(suffix) VecI64##suffix
#define FUNC(suffix) VecI64_##suffix
#define ITEM_TYPE_STR "i64"
#define ITEM_TYPE_MAGIC VEC_ITEM_TYPE_I64
#define ITEM_C_TYPE int64_t
#define FEATURES Vec_I64API

#define BOX_ITEM VecI64_BoxItem
#define UNBOX_ITEM VecI64_UnboxItem
#define IS_UNBOX_ERROR VecI64_IsUnboxError
#define BUFFER_FORMAT_CHAR_OK(c) ((c) == 'q' || ((c) == 'l' && sizeof(long) == 8))
#define BUFFER_FORMAT "q"

#include "vec_template.c"
