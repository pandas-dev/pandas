import socket

HOST = '127.0.0.1'  # Standard loopback interface address (localhost)
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)

with

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

sock.bind((HOST, PORT))
sock.listen()

while True:
    conn, addr = sock.accept()
    with conn:
        while True:
            conn.sendall(b'some\rdata\nbyt\1\xffest\r\nadslfkja\n\raslkdj')
