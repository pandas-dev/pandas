#define VEC VecBool
#define VEC_TYPE VecBoolType
#define VEC_OBJECT VecBoolObject
#define BUF_OBJECT VecBoolBufObject
#define BUF_TYPE VecBoolBufType
#define NAME(suffix) VecBool##suffix
#define FUNC(suffix) VecBool_##suffix
#define ITEM_TYPE_STR "bool"
#define ITEM_TYPE_MAGIC VEC_ITEM_TYPE_BOOL
#define ITEM_C_TYPE char
#define FEATURES Vec_BoolAPI

#define BOX_ITEM VecBool_BoxItem
#define UNBOX_ITEM VecBool_UnboxItem
#define IS_UNBOX_ERROR VecBool_IsUnboxError
#define BUFFER_FORMAT "b"

#include "vec_template.c"
