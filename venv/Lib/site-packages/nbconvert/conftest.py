"""pytest configuration."""

import asyncio
import os

if os.name == "nt":
    asyncio.set_event_loop_policy(
        asyncio.WindowsSelectorEventLoopPolicy()  # type:ignore[attr-defined]
    )
