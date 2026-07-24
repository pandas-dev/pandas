#ifndef LIBRT_TIME_API_H
#define LIBRT_TIME_API_H

#include "librt_time.h"

int
import_librt_time(void);

extern void *LibRTTime_API[LIBRT_TIME_API_LEN];

#define LibRTTime_ABIVersion (*(int (*)(void)) LibRTTime_API[0])
#define LibRTTime_APIVersion (*(int (*)(void)) LibRTTime_API[1])
#define LibRTTime_time (*(double (*)(void)) LibRTTime_API[2])

#endif  // LIBRT_TIME_API_H
