import io


socket_bytes = io.BytesIO(b"So\x01me\r\nbytes\rto\nparsB")

byte_joiner = b''.join

list(socket_bytes)

def split_on(buffer, *spliters):
    if not spliters:
        spliters = {b'\n', b'\r'}
    else:
        spliters = set(spliters)
    line = []
    while True:
        b = buffer.read(1)
        split = b in {b'\n', b'\r'}

        if split or not b:
            if line:
                yield byte_joiner(line)
        if split:
            line = []
        elif not b:
            return
        else:
            line.append(b)

gen = split_on(socket_bytes)
list(gen)
