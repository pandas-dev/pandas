cdef extern from "datetime.h":
    
    ctypedef extern class datetime.datetime [object PyDateTime_DateTime]:
        cdef int *data
        cdef long hashcode
        cdef char hastzinfo
        
    int PyDateTime_GET_YEAR(datetime o)
    int PyDateTime_GET_MONTH(datetime o)
    int PyDateTime_GET_DAY(datetime o)
    int PyDateTime_TIME_GET_HOUR(datetime o)
    int PyDateTime_TIME_GET_MINUTE(datetime o)
    int PyDateTime_TIME_GET_SECOND(datetime o)
    int PyDateTime_TIME_GET_MICROSECOND(datetime o)