import os
import resource

import psutil
from logzero import logger

RESOURCE_LOG = """---RESOURCE---
User time: {0}
System time: {1}
Max resident size: {2}
Block input operations: {3}
Block output operations: {4}
---MEMORY_INFO---
RSS: {5}
VMS: {6}
Data: {7}
"""


def log_resource_usage(step: str):
    mem = psutil.Process(os.getpid()).memory_info()
    res = resource.getrusage(resource.RUSAGE_SELF)
    # MacOs only
    try:
        data = mem.data
    except AttributeError:
        data = 0
    res_log = RESOURCE_LOG.format(
        res.ru_utime,
        res.ru_stime,
        res.ru_maxrss,
        res.ru_inblock,
        res.ru_oublock,
        mem.rss,
        mem.vms,
        data,
    )
    logger.info(f"[resource][{step}] {res_log}")
