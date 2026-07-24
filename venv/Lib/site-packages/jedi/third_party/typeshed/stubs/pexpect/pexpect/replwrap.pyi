import sys
from _typeshed import Incomplete

PY3: Incomplete
basestring = str
PEXPECT_PROMPT: str
PEXPECT_CONTINUATION_PROMPT: str

class REPLWrapper:
    child: Incomplete
    prompt: Incomplete
    continuation_prompt: Incomplete
    def __init__(
        self,
        cmd_or_spawn,
        orig_prompt,
        prompt_change,
        new_prompt="[PEXPECT_PROMPT>",
        continuation_prompt="[PEXPECT_PROMPT+",
        extra_init_cmd=None,
    ) -> None: ...
    def set_prompt(self, orig_prompt, prompt_change) -> None: ...
    def run_command(self, command, timeout: float | None = -1, async_: bool = False): ...

def python(command: str = sys.executable): ...
def bash(command: str = "bash"): ...
def zsh(command: str = "zsh", args=("--no-rcs", "-V", "+Z")): ...
