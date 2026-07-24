#define VEC VecU8
#define VEC_TYPE VecU8Type
#define VEC_OBJECT VecU8Object
#define BUF_OBJECT VecU8BufObject
#define BUF_TYPE VecU8BufType
#define NAME(suffix) VecU8##suffix
#define FUNC(suffix) VecU8_##suffix
#define ITEM_TYPE_STR "u8"
#define ITEM_TYPE_MAGIC VEC_ITEM_TYPE_U8
#define ITEM_C_TYPE uint8_t
#define FEATURES Vec_U8API

#define BOX_ITEM VecU8_BoxItem
#define UNBOX_ITEM VecU8_UnboxItem
#define IS_UNBOX_ERROR VecU8_IsUnboxError
#define BUFFER_FORMAT_CHAR_OK(c) ((c) == 'B' || (c) == 'c')
#define BUFFER_FORMAT "B"

#include "vec_template.c"
