from .responses import SageMakerResponse

url_bases = [
    r"https?://api\.sagemaker\.(.+)\.amazonaws.com",
]

url_paths = {
    "{0}/$": SageMakerResponse.dispatch,
}
