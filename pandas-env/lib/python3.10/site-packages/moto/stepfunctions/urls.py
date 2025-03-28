from .responses import StepFunctionResponse

url_bases = [r"https?://states\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": StepFunctionResponse.dispatch}
