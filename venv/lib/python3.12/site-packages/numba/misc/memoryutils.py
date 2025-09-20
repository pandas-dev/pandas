"""
Memory monitoring utilities for measuring memory usage.

Example usage:
    tracker = MemoryTracker("my_function")
    with tracker.monitor():
        my_function()
    # Access data: tracker.rss_delta, tracker.duration, etc.
    # Get formatted string: tracker.get_summary()
"""
from __future__ import annotations
import os
import contextlib
import time
from typing import Dict, Optional


try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


IS_SUPPORTED = _HAS_PSUTIL


def get_available_memory() -> Optional[int]:
    """
    Get current available system memory in bytes.

    Used for memory threshold checking in parallel test execution.

    Returns:
        int or None: Available memory in bytes, or None if unavailable
    """
    if _HAS_PSUTIL:
        try:
            sys_mem = psutil.virtual_memory()
            return sys_mem.available
        except Exception:
            pass

    return None


def get_memory_usage() -> Dict[str, Optional[int]]:
    """
    Get memory usage information needed for monitoring.

    Returns only RSS and available memory which are the fields
    actually used by the MemoryTracker.

    Returns:
        dict: Memory usage information including:
            - rss: Current process RSS (physical memory currently used)
            - available: Available system memory
    """
    memory_info = {}

    if _HAS_PSUTIL:
        try:
            # Get current process RSS
            process = psutil.Process(os.getpid())
            mem_info = process.memory_info()
            memory_info["rss"] = mem_info.rss

            # Get system available memory
            sys_mem = psutil.virtual_memory()
            memory_info["available"] = sys_mem.available

        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

    # Set defaults if unavailable
    if "rss" not in memory_info:
        memory_info["rss"] = None
    if "available" not in memory_info:
        memory_info["available"] = None

    return memory_info


class MemoryTracker:
    """
    A simple memory monitor that tracks RSS delta and timing.

    Stores monitoring data in instance attributes for later access.
    Each instance is typically used for monitoring a single operation.
    """
    pid: int
    name: str
    start_time: float | None
    end_time: float | None
    start_memory: Dict[str, int | None] | None
    end_memory: Dict[str, int | None] | None
    duration: float | None
    rss_delta: int | None

    def __init__(self, name: str):
        """Initialize a MemoryTracker with empty monitoring data."""
        self.pid = os.getpid()
        self.name = name
        self.start_time = None
        self.end_time = None
        self.start_memory = None
        self.end_memory = None
        self.duration = None
        self.rss_delta = None

    @contextlib.contextmanager
    def monitor(self):
        """
        Context manager to monitor memory usage during function execution.

        Records start/end memory usage and timing, calculates RSS delta,
        and stores all data in instance attributes.

        Args:
            name (str): Name/identifier for the function or operation being
                        monitored

        Yields:
            self: The MemoryTracker instance for accessing stored data
        """
        # Store data in self and record start time and memory usage
        self.start_time = time.time()
        self.start_memory = get_memory_usage()

        try:
            yield self
        finally:
            # Record end time and memory usage
            self.end_time = time.time()
            self.end_memory = get_memory_usage()
            self.duration = self.end_time - self.start_time

            # Calculate RSS delta
            start_rss = self.start_memory.get("rss", 0)
            end_rss = self.end_memory.get("rss", 0)
            self.rss_delta = ((end_rss - start_rss)
                              if start_rss and end_rss else 0)

    def get_summary(self) -> str:
        """
        Return a formatted summary of the memory monitoring data.

        Formats the stored monitoring data into a human-readable string
        containing name, PID, RSS delta, available memory, duration,
        and start time.

        Returns:
            str: Formatted summary string with monitoring results

        Note:
            Should be called after monitor() context has completed
            to ensure all data is available.
        """
        if self.start_memory is None or self.end_memory is None:
            raise ValueError("Memory monitoring data not available")

        current_available = self.end_memory.get("available")

        def format_bytes(bytes_val, show_sign=False):
            """Convert bytes to human readable format"""
            if bytes_val is None:
                return "N/A"
            if bytes_val == 0:
                return "0 B"

            sign = ""
            if show_sign:
                sign = "-" if bytes_val < 0 else "+"
            bytes_val = abs(bytes_val)

            for unit in ["B", "KB", "MB", "GB"]:
                if bytes_val < 1024.0:
                    return f"{sign}{bytes_val:.2f} {unit}"
                bytes_val /= 1024.0
            return f"{sign}{bytes_val:.2f} TB"

        start_ts = time.strftime("%H:%M:%S", time.localtime(self.start_time))
        start_rss = self.start_memory.get("rss", 0)
        end_rss = self.end_memory.get("rss", 0)

        buf = [
            f"Name: {self.name}",
            f"PID: {self.pid}",
            f"Start: {start_ts}",
            f"Duration: {self.duration:.3f}s",
            f"Start RSS: {format_bytes(start_rss)}",
            f"End RSS: {format_bytes(end_rss)}",
            f"RSS delta: {format_bytes(self.rss_delta, show_sign=True)}",
            f"Avail memory: {format_bytes(current_available)}",
        ]
        return ' | '.join(buf)
