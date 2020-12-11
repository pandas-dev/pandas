class DNS_Resolver:
    def __init__(self):
        self.ip_address = "11.111.111.11"
        self.port = "8502"
        self.cache = {
        “.”: “199.7.83.42” #root dns server
        }

    def do_lookup(self, recv_socket, host):
        domain = host.split('.')
        domain.reverse()
        temp = resolve_dns(self.cache["."], domain[0])
        name = domain[0]
        domain.pop(0)
        for dom in domain:
            name = dom + "." + name
            temp = resolve_dns(temp, name)
        self.cache[host] = temp
        self.send_response(recv_socket, host, temp)


    def send_response(self. recv_socket, host, ip):
        send_msg = {"host": host, "resolved_ip": ip}
        msg = json.dumps(send_msg)
        recv_socket.sendall(msg.encode("utf-8"))



    def main(self):
        threads = []
        while True:
            recv_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            recv_socket.bind((self.ip_address, self.port))
            message = []
            while True:
                try:
                    data = clientsocket.recv(4096)
                except socket.timeout:
                    continue
                if not data:
                    break
                message.append(data)
            # parsing msg received into dict
            message_bytes = b''.join(message)
            message_str = message_bytes.decode("utf-8")
            message_dict = json.loads(message_str)
            host = message_dict["host"]
            if host in self.cache:
                self.send_response(recv_socket, host, self.cache[host])
            else:
                threads.append(threading.Thread(target=self.do_lookup, args=(recv_socket, host)))
                threads[-1].start()

        for thread in threads:
            thread.join()

        recv_socket.close()