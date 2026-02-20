from .responses import Ec2InstanceConnectResponse

url_bases = [r"https?://ec2-instance-connect\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": Ec2InstanceConnectResponse.dispatch}
