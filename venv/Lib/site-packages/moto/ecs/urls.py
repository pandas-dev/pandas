from .responses import EC2ContainerServiceResponse

url_bases = [r"https?://ecs\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": EC2ContainerServiceResponse.dispatch}
