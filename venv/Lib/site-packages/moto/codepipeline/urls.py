from .responses import CodePipelineResponse

url_bases = [r"https?://codepipeline\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": CodePipelineResponse.dispatch}
