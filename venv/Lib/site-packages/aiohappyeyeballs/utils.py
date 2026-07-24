"""Utility functions for aiohappyeyeballs."""

import ipaddress
import socket

from .types import AddrInfoType


def addr_to_addr_infos(
    addr: tuple[str, int, int, int] | tuple[str, int, int] | tuple[str, int] | None,
) -> list[AddrInfoType] | None:
    """Convert an address tuple to a list of addr_info tuples."""
    if addr is None:
        return None
    host = addr[0]
    port = addr[1]
    is_ipv6 = ":" in host
    if is_ipv6:
        flowinfo = 0
        scopeid = 0
        addr_len = len(addr)
        if addr_len >= 4:
            scopeid = addr[3]  # type: ignore[misc]
        if addr_len >= 3:
            flowinfo = addr[2]  # type: ignore[misc]
        addr = (host, port, flowinfo, scopeid)
        family = socket.AF_INET6
    else:
        addr = (host, port)
        family = socket.AF_INET
    return [(family, socket.SOCK_STREAM, socket.IPPROTO_TCP, "", addr)]


def pop_addr_infos_interleave(
    addr_infos: list[AddrInfoType], interleave: int | None = None
) -> None:
    """
    Pop addr_info from the list of addr_infos by family up to interleave times.

    The interleave parameter is used to know how many addr_infos for
    each family should be popped of the top of the list.
    """
    if interleave is None:
        interleave = 1
    seen: dict[int, int] = {}
    kept: list[AddrInfoType] = []
    for addr_info in addr_infos:
        family = addr_info[0]
        count = seen.get(family, 0)
        if count >= interleave:
            kept.append(addr_info)
        seen[family] = count + 1
    addr_infos[:] = kept


def _addr_tuple_to_ip_address(
    addr: tuple[str, int] | tuple[str, int, int, int],
) -> tuple[ipaddress.IPv4Address, int] | tuple[ipaddress.IPv6Address, int, int, int]:
    """Convert an address tuple to an IPv4Address."""
    return (ipaddress.ip_address(addr[0]), *addr[1:])


def remove_addr_infos(
    addr_infos: list[AddrInfoType],
    addr: tuple[str, int] | tuple[str, int, int, int],
) -> None:
    """
    Remove an address from the list of addr_infos.

    The addr value is typically the return value of
    sock.getpeername().
    """
    kept = [ai for ai in addr_infos if ai[-1] != addr]
    if len(kept) == len(addr_infos):
        # Slow path in case addr is formatted differently
        match_addr = _addr_tuple_to_ip_address(addr)
        kept = [
            ai for ai in addr_infos if _addr_tuple_to_ip_address(ai[-1]) != match_addr
        ]
    if len(kept) == len(addr_infos):
        raise ValueError(f"Address {addr} not found in addr_infos")
    addr_infos[:] = kept
