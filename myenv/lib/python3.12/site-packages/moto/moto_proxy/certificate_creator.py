import os
import threading
import time
from subprocess import PIPE, Popen
from uuid import uuid4

from . import debug, info


def join_with_script_dir(path: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


class CertificateCreator:
    cakey = join_with_script_dir("ca.key")
    cacert = join_with_script_dir("ca.crt")
    certkey = join_with_script_dir("cert.key")
    certdir = join_with_script_dir("certs/")

    lock = threading.Lock()

    def validate(self) -> None:
        # Verify the CertificateAuthority files exist
        if not os.path.isfile(CertificateCreator.cakey):
            raise Exception(f"Cannot find {CertificateCreator.cakey}")
        if not os.path.isfile(CertificateCreator.cacert):
            raise Exception(f"Cannot find {CertificateCreator.cacert}")
        if not os.path.isfile(CertificateCreator.certkey):
            raise Exception(f"Cannot find {CertificateCreator.certkey}")
        if not os.path.isdir(CertificateCreator.certdir):
            raise Exception(f"Cannot find {CertificateCreator.certdir}")
        # Verify the `certs` dir is reachable
        try:
            test_file_location = f"{CertificateCreator.certdir}/{uuid4()}.txt"
            debug(
                f"Writing test file to {test_file_location} to verify the directory is writable..."
            )
            with open(test_file_location, "w") as file:
                file.write("test")
            os.remove(test_file_location)
        except Exception:
            info("Failed to write test file")
            info(
                f"The directory {CertificateCreator.certdir} does not seem to be writable"
            )
            raise

    def create(self, path: str) -> str:
        """
        Create an SSL certificate for the supplied hostname.
        This method will return a path to the certificate.
        """
        full_name = path.split(":")[0]

        with CertificateCreator.lock:
            # We don't want to create certificates for every possible endpoint
            # Especially with randomly named S3-buckets

            # We can create certificates that match wildcards to reduce the total number
            # For example:
            # Hostname:      somebucket.s3.amazonaws.com
            # Certificate:   *.s3.amazonaws.com
            #
            # All requests that match this wildcard certificate will reuse it

            wildcard_name = f"*.{'.'.join(full_name.split('.')[1:])}"
            server_csr = f"{self.certdir.rstrip('/')}/{wildcard_name}.csr"

            # Verify if the certificate already exists
            certpath = f"{self.certdir.rstrip('/')}/{wildcard_name}.crt"
            if not os.path.isfile(certpath):
                # Create a Config-file that contains the wildcard-name
                with open(f"{self.certdir.rstrip('/')}/req.conf.tmpl", "r") as f:
                    config_template = f.read()
                config_template = config_template.replace("{{full_name}}", full_name)
                config_template = config_template.replace(
                    "{{wildcard_name}}", wildcard_name
                )
                config_template_name = (
                    f"{self.certdir.rstrip('/')}/{wildcard_name}.conf"
                )
                with open(config_template_name, "w") as f:
                    f.write(config_template)

                # Create an Certificate Signing Request
                #
                subject = f"/CN={full_name}"[0:64]
                commands = [
                    "openssl",
                    "req",
                    "-new",
                    "-key",
                    self.certkey,
                    "-out",
                    server_csr,
                ]
                commands.extend(["-subj", subject, "-config", config_template_name])

                p1 = Popen(commands)
                p1.communicate()
                debug(f"Created CSR in {server_csr}")

                # Create the actual certificate used by the requests
                p2 = Popen(
                    [
                        "openssl",
                        "x509",
                        "-req",
                        "-in",
                        server_csr,
                        "-days",
                        "3650",
                        "-CA",
                        self.cacert,
                        "-CAkey",
                        self.cakey,
                        "-set_serial",
                        f"{int(time.time() * 1000)}",
                        "-out",
                        certpath,
                        "-extensions",
                        "req_ext",
                        "-extfile",
                        config_template_name,
                    ],
                    stderr=PIPE,
                )
                p2.communicate()
                debug(f"Created certificate for {path} called {certpath}")
                os.remove(server_csr)
                os.remove(config_template_name)
                debug(f"Removed intermediate certificates for {certpath}")
        return certpath
