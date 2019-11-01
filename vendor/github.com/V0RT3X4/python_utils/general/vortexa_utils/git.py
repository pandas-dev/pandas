import os

from git import Git, Repo
from logzero import logger
from tenacity import retry, wait_fixed, stop_after_attempt


@retry(wait=wait_fixed(10), stop=stop_after_attempt(3))
def clone_repo(repo_url: str, path: str, ssh_key: str):
    os.environ['GIT_SSH_COMMAND'] = f'ssh -i {ssh_key}'
    with Git().custom_environment():
        logger.info('Cloning git repo %s to %s', repo_url, path)
        Repo.clone_from(repo_url, path, branch='master')
        logger.info('Repo cloned successfully')
