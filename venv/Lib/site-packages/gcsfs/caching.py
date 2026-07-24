from collections import deque

from fsspec.caching import BaseCache, register_cache


class ReadAheadChunked(BaseCache):
    """
    An optimized ReadAhead cache that fetches multiple chunks in a single
    HTTP request but manages them as separate bytes objects to avoid
    expensive memory slicing.

    While this approach primarily optimizes for CPU and memory allocation overhead,
    it strictly maintains the same semantics as the existing readahead cache.
    For example, if a user requests 5MB and the cache fetches 10MB, it serves the
    requested 5MB but retains that data in memory to handle potential backward seeks.
    This mirrors the standard readahead behavior, which does not eagerly discard served
    chunks until a new fetch is required.
    """

    name = "readahead_chunked"

    def __init__(self, blocksize: int, fetcher, size: int) -> None:
        super().__init__(blocksize, fetcher, size)
        self.chunks = deque()  # Entries: (start, end, data_bytes)

    @property
    def cache(self):
        """
        Compatibility property for tests/legacy code that expects 'cache'
        to be a single bytestring.

        WARNING: Accessing this property forces a memory copy of the
        entire current buffer, negating the Zero-Copy optimization
        of ReadAheadChunked. Use for debugging/testing only.
        """
        if not self.chunks:
            return b""
        return b"".join(chunk[2] for chunk in self.chunks)

    def _fetch(self, start: int | None, end: int | None) -> bytes:
        if start is None:
            start = 0
        if end is None or end > self.size:
            end = self.size
        if start >= self.size:
            return b""

        # Handle backward seeks that go beyond the start of our cache window
        if self.chunks and self.chunks[0][0] > start:
            self.chunks.clear()

        parts = []
        current_pos = start

        # Satisfy as much as possible from the existing cache (Zero-Copy)
        for c_start, c_end, c_data in self.chunks:
            if c_end <= start:
                continue  # Skip chunks completely before our window

            if c_start >= end:
                break  # If we've reached chunks completely past our window, stop

            if c_end > current_pos:
                slice_start = max(0, current_pos - c_start)
                slice_end = min(len(c_data), end - c_start)

                if slice_start == 0 and slice_end == len(c_data):
                    # Zero-copy: Direct reference to the full object
                    parts.append(c_data)
                else:
                    # Slicing creates a copy, but it's unavoidable for partials
                    parts.append(c_data[slice_start:slice_end])

                current_pos += slice_end - slice_start

        # Fetch missing data if necessary
        should_fetch_backend = current_pos < end
        if should_fetch_backend:
            # On a cache miss, we replace the entire window (standard readahead behavior)
            self.chunks.clear()

            missing_len = min(self.size - current_pos, end - current_pos)
            readahead_block = min(
                self.size - (current_pos + missing_len), self.blocksize
            )

            self.miss_count += 1
            chunk_lengths = [missing_len]
            if readahead_block > 0:
                chunk_lengths.append(readahead_block)

            # Vector read call
            new_chunks = self.fetcher(start=current_pos, chunk_lengths=chunk_lengths)

            # Process the requested data
            req_data = new_chunks[0]
            self.chunks.append((current_pos, current_pos + len(req_data), req_data))
            self.total_requested_bytes += len(req_data)
            parts.append(req_data)

            # Process the readahead data (if any)
            if len(new_chunks) > 1:
                ra_data = new_chunks[1]
                ra_start = current_pos + len(req_data)
                self.chunks.append((ra_start, ra_start + len(ra_data), ra_data))
                self.total_requested_bytes += len(ra_data)

        if not parts:
            return b""

        if not should_fetch_backend:
            self.hit_count += 1

        # Optimization: return the single object directly if possible
        if len(parts) == 1:
            return parts[0]

        return b"".join(parts)


register_cache(ReadAheadChunked, clobber=True)
