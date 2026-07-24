#define VEC VecI32
#define VEC_TYPE VecI32Type
#define VEC_OBJECT VecI32Object
#define BUF_OBJECT VecI32BufObject
#define BUF_TYPE VecI32BufType
#define NAME(suffix) VecI32##suffix
#define FUNC(suffix) VecI32_##suffix
#define ITEM_TYPE_STR "i32"
#define ITEM_TYPE_MAGIC VEC_ITEM_TYPE_I32
#define ITEM_C_TYPE int32_t
#define FEATURES Vec_I32API

#define BOX_ITEM VecI32_BoxItem
#define UNBOX_ITEM VecI32_UnboxItem
#define IS_UNBOX_ERROR VecI32_IsUnboxError
#define BUFFER_FORMAT_CHAR_OK(c) ((c) == 'i' || ((c) == 'l' && sizeof(long) == 4))
#define BUFFER_FORMAT "i"

#include "vec_template.c"
