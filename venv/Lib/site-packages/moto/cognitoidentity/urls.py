from .responses import CognitoIdentityResponse

url_bases = [r"https?://cognito-identity\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": CognitoIdentityResponse.dispatch}
