from tempfile import NamedTemporaryFile

import boto3
from logzero import logger
from pandas import DataFrame, read_hdf, read_csv


class S3Client:
    def __init__(self, s3_bucket: str):
        self.s3 = boto3.client("s3")
        self.s3_bucket = s3_bucket

    def upload(self, filename: str, s3_key: str, owner_acl: bool = True):
        logger.info("[s3] Started uploading: %s", s3_key)
        self.s3.upload_file(filename, self.s3_bucket, s3_key)
        logger.info("[s3] Finished uploading: %s", s3_key)
        if owner_acl:
            self.s3.put_object_acl(
                ACL="bucket-owner-full-control", Bucket=self.s3_bucket, Key=s3_key
            )
            logger.info("[s3] bucket-owner-full-control ACL set")

    def hdf_pd(self, filename) -> DataFrame:
        return self.__s3_pd__(filename, "hdf")

    def csv_pd(self, filename) -> DataFrame:
        return self.__s3_pd__(filename, "csv")

    def copy(self, src, dest, owner_acl: bool = True):
        copy_source = {"Bucket": self.s3_bucket, "Key": src}
        self.s3.copy_object(CopySource=copy_source, Bucket=self.s3_bucket, Key=dest)
        if owner_acl:
            self.s3.put_object_acl(
                ACL="bucket-owner-full-control", Bucket=self.s3_bucket, Key=dest
            )
            logger.info("[s3] bucket-owner-full-control ACL set")

    def __s3_pd__(self, filename, filetype) -> DataFrame:
        with NamedTemporaryFile("wb") as f:
            logger.info(f"[s3] Started downloading: s3://{self.s3_bucket}/{filename}")
            self.s3.download_fileobj(self.s3_bucket, filename, f)
            f.flush()
            logger.info(f"[s3] Finished downloading: s3://{self.s3_bucket}/{filename}")
            logger.info("[pandas] Started loading: %s", filename)
            if filetype == "hdf":
                df: DataFrame = read_hdf(f.name)
            elif filetype == "csv":
                df: DataFrame = read_csv(f.name)
            logger.info("[pandas] Finished loading: %s", filename)
        return df
