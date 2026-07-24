#define VEC VecI16
#define VEC_TYPE VecI16Type
#define VEC_OBJECT VecI16Object
#define BUF_OBJECT VecI16BufObject
#define BUF_TYPE VecI16BufType
#define NAME(suffix) VecI16##suffix
#define FUNC(suffix) VecI16_##suffix
#define ITEM_TYPE_STR "i16"
#define ITEM_TYPE_MAGIC VEC_ITEM_TYPE_I16
#define ITEM_C_TYPE int16_t
#define FEATURES Vec_I16API

#define BOX_ITEM VecI16_BoxItem
#define UNBOX_ITEM VecI16_UnboxItem
#define IS_UNBOX_ERROR VecI16_IsUnboxError
#define BUFFER_FORMAT_CHAR_OK(c) ((c) == 'h')
#define BUFFER_FORMAT "h"

#include "vec_template.c"
