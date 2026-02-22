# AWS has a 20 second idle timeout:
#   https://web.archive.org/web/20150926192339/https://forums.aws.amazon.com/message.jspa?messageID=215367
# and aiohttp default timeout is 30s so we set it to something
# reasonable here
DEFAULT_KEEPALIVE_TIMEOUT = 12
