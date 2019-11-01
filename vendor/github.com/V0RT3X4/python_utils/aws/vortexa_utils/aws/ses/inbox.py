# @Author: richard
# @Date:   2018-12-06T18:06:25+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T19:36:09+00:00
from typing import Iterable
import logging
from datetime import datetime
from functools import wraps
import boto3
# cd aws/vortexa_utils
# import aws.s3.client_side_encryption.client as client
import vortexa_utils.aws.s3.client_side_encryption.client as client
import email
import email.policy
import email.parser
from email.message import EmailMessage

logger = logging.getLogger(__name__)


class Inbox(object):
    """Short summary.

    Parameters
    ----------
    default_bucket : str
        Default s3  Bucket to assosiate the inbox with.

    """

    def __init__(self, default_bucket: str = None, **kwargs):
        """Short summary.

        Parameters
        ----------
        default_bucket : str
            Default s3  Bucket to assosiate the inbox with.
        strict : bool
            When True will not fetch unencrypted emails. Defaults to False.
        **kwargs : dict
            **`kwargs` to pass to `s3.client`.

        """
        self.bucket = default_bucket
        self.s3crypto = client.Client(**kwargs)
        self.s3 = self.s3crypto.s3
        # Specify the default policy for email parsing else Parser defaults to
        # email.policy.compat32 for python 3 and 2 compatibility
        self.parser = email.parser.BytesParser(policy=email.policy.default)

    def get_email(self, Key: str, Bucket: str = None) -> EmailMessage:
        """Get `EmailMessage` Object from `Bucket`.

        Parameters
        ----------
        Key : str
            `Key` name of email in s3.
        Bucket : str
            s3 `Bucket` to look for email, will search `self.bucket` if `None`.

        Returns
        -------
        email.message.EmailMessage
            Email object.

        """
        Bucket = Bucket or self.bucket
        if Bucket is None:
            raise ValueError("Bucket not set")
        with self.s3crypto.get_object(Bucket=Bucket, Key=Key) as io:
            return self.parser.parse(io)

    def list_objects(
            self,
            Bucket: str = None,
            Path: str = None,
            Begin: datetime = None,
            Until: datetime = None):
        # type:  (...) -> Iterable['boto3.resources.factory.s3.ObjectSummary']
        """List all objects in `Bucket` prefixed by `Path`.

        Parameters
        ----------
        Bucket : str
            S3 `Bucket` to look for emails will search `self.bucket` if `None`.
        Path : str
            The `Path` prefix to filter the emails by, no filter if `None`.
        Begin : datetime
            Filter object from this datetime.
        Until : datetime = None
            Filter objects untill this datetime.

        Returns
        -------
        iterable boto3.resources.factory.s3.ObjectSummary
            List of matching email objects.

        """
        bucket = boto3.resource('s3').Bucket(Bucket or self.bucket)
        objs = bucket.objects.filter(Prefix=Path)
        if Begin:
            objs = (obj for obj in objs if obj.last_modified >= Begin)
        if Until:
            objs = (obj for obj in objs if obj.last_modified <= Until)

        if Begin is None and Until is None:
            # if no timestamps dont bother sorting
            return objs
        return sorted(objs, key=lambda o: o.last_modified)

    @wraps(list_objects, assigned=('__annotations__',))
    def list_emails(self, **kwargs) -> Iterable[EmailMessage]:
        """List all emails in `Bucket` prefixed by `Path`.

        Parameters
        ----------
        Bucket : str
            S3 `Bucket` to look for emails will search `self.bucket` if `None`.
        Path : str
            The `Path` prefix to filter the emails by, no filter if `None`.
        Begin : datetime
            Filter object from this datetime.
        Until : datetime = None
            Filter objects untill this datetime.

        Returns
        -------
        iterable emails
            List of matching email objects.

        Examples
        -------
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>> inbox = Inbox()
        >>> inboc.list_emails('/some/sub/folder')

        """
        objects = self.list_objects(**kwargs)
        for obj in objects:
            yield self.get_email(obj.key, obj.bucket_name)
