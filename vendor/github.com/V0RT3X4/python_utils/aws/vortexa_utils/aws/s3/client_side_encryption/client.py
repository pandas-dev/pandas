# @Author: richard
# @Date:   2018-11-28T15:15:54+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T18:07:33+00:00
import boto3
from .kms_cipher_provider import KMSCipherProvider
from .decrypt_handeler import DecryptHandeler


class Client(object):
    """
    Client Side Encryption S3 Client.

    Attributes
    ----------
    s3 : botocore.client.S3
    cipher_provider : .cipher_provider.CipherProvider

    Methods
    -------
    get_object(Bucket, Key)
        get and decrypt an object from s3

    """

    def __init__(
            self,
            client=None,
            cipher_provider=None,
            key_id=None,
            strict=None,
            **kwargs):
        """Initilises the client side encryption s3 client.

        Parameters
        ----------
        client : botocore.client.S3
            Optional S3 client to use for s3 interaction
            Will create client if not set.

        cipher_provider : CipherProvider
            Optional `CipherProvider` to provide encryption cipher
            Will default to `KMSCipherProvider()` if not set.

        key_id : str
            The kms `key id`, `alias` or `aws::arn`
            for the `KMSCipherProvider`.

        region_name : str
             The region for the kms and s3 client resources.

        """
        region_name = kwargs.get('region')
        self.s3 = client or boto3.client('s3', **kwargs)
        self.cipher_provider = (
            cipher_provider or KMSCipherProvider(
                key_id=key_id,
                region_name=region_name
            )
        )
        self.strict = strict

    def get_object(self, Bucket, Key):
        """Retrieve object from Amazon S3.

        See also:
        `AWS API Documentation <https://docs.aws.amazon.com/goto/WebAPI/s3-2006-03-01/GetObject>`_

        `AWS Client Side Encryption <https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingClientSideEncryption.html>`_

        Parameters
        ----------
        Bucket : str
            **[REQUIRED]** The Bucket
        Key : str
            **[REQUIRED]** The Path Key in the Bucket

        """
        # location_info = self.s3.get_bucket_location(Bucket=Bucket)
        # bucket_region = location_info['LocationConstraint']

        obj = self.s3.get_object(Bucket=Bucket, Key=Key)
        handeler = DecryptHandeler(obj, self, self.strict)
        return handeler.decrypt()

    def object_encrypted(self, Bucket, Key) -> bool:
        """Check if object has encryption envelope.

        Parameters
        ----------
        Bucket : str
            **[REQUIRED]** The Bucket
        Key : str
            **[REQUIRED]** The Path Key in the Bucket

        Returns
        -------
        bool

        """
        obj = self.s3.head_object(Bucket=Bucket, Key=Key)
        handeler = DecryptHandeler(obj, self)
        return handeler.extract_envelop() is not None
