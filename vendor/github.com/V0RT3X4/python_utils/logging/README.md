# Vortexa Utils Logging Helpers

Small helper scripts to abstract logging-related boilerplate code.  


## log_unhandled_exceptions

Function decorator, designed to be wrapped around any `main()` (or equivalent) function, to capture errors, prefix them with `ERROR`, and raise them in-line, when executed in AWS Batch.

### Problem:

AWS Batch jobs all output logs onto CloudWatch Log Group (`/aws/batch/job`).  Therefore, to raise specific alarms, python jobs should use logging, with the logger pattern containing a unique identifier for the job (such as the job/repo name), so the CloudWatch can filter logs and look for specific exceptions.

When Errors are raised by a python program logging to CloudWatch, the loger pattern and the Error/stacktrace are output on 2 consecutive lines.  CloudWatch Alarm triggers can only look for patterns combinations which are in-line, therefore, for a CloudWatch Alarm to be raised when a job fails, the logger pattern and some form of identifiable error key most be printed in-line.


### Solution:

`log_unhandled_exceptions` decorator, can be wrapped around main executing functions, and if any errors are raised during run-time, will capture these errors, and raise them in-line with the logging pattern, using the common pattern `ERROR: <error>`.  CloudWatch alerts can now be set to look for (1) the unique logging pattern of the project (i.e. name) and (2) the key `ERROR`, to raise targeted alerts.  The full stacktrace will still be output to Cloudwatch logs.

### Usage:

```python
from vortexa_utils.logging import log_unhandled_exceptions

# The following is the logger set-up boilerplate code.
# This can be done as below, or imported from a project-logger dir.  
# The following is only intended as a sample and should not be copied without understanding what is happening.
import logging

logger = logging.getLogger(__name__)
log_format = logging.Formatter(
    f"PROJECT_NAME:%(name)s:%(message)s"
) # Only a sample format, can be designed at will, as long as unique identifier (e.g. PROJECT_NAME) is included
logger.setFormatter(log_format)
logger.setLevel(logging.INFO)

@log_unhandled_exceptions(logger)
def main():
    return int(1) + str('two')
    
if __name__ == "__main__":
    main()
```

Code snippet above would return:

```
PROJECT_NAME:__main__:ERROR: unsupported operan types(s) for +: 'int' and 'str'
    Traceback (most recent call last):
        ... <remaining stacktrace> ...  
    TypeError: unsupported operand type(s) for +: 'int' and 'str'
```

As a result, a cloudwatch alarm can now be set on the pattern `PROJECT_NAME ERROR`
