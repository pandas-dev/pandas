import json

from moto.core.responses import BaseResponse

from .models import DataSyncBackend, Location, datasync_backends


class DataSyncResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="datasync")

    @property
    def datasync_backend(self) -> DataSyncBackend:
        return datasync_backends[self.current_account][self.region]

    def list_locations(self) -> str:
        locations = list()
        for arn, location in self.datasync_backend.locations.items():
            locations.append({"LocationArn": arn, "LocationUri": location.uri})
        return json.dumps({"Locations": locations})

    def _get_location(self, location_arn: str, typ: str) -> Location:
        return self.datasync_backend._get_location(location_arn, typ)

    def create_location_s3(self) -> str:
        # s3://bucket_name/folder/
        s3_bucket_arn = self._get_param("S3BucketArn")
        subdirectory = self._get_param("Subdirectory")
        metadata = {"S3Config": self._get_param("S3Config")}
        location_uri_elts = ["s3:/", s3_bucket_arn.split(":")[-1]]
        if subdirectory:
            location_uri_elts.append(subdirectory)
        location_uri = "/".join(location_uri_elts)
        arn = self.datasync_backend.create_location(
            location_uri, metadata=metadata, typ="S3"
        )
        return json.dumps({"LocationArn": arn})

    def describe_location_s3(self) -> str:
        location_arn = self._get_param("LocationArn")
        location = self._get_location(location_arn, typ="S3")
        return json.dumps(
            {
                "LocationArn": location.arn,
                "LocationUri": location.uri,
                "S3Config": location.metadata["S3Config"],
            }
        )

    def create_location_smb(self) -> str:
        # smb://smb.share.fqdn/AWS_Test/
        subdirectory = self._get_param("Subdirectory")
        server_hostname = self._get_param("ServerHostname")
        metadata = {
            "AgentArns": self._get_param("AgentArns"),
            "User": self._get_param("User"),
            "Domain": self._get_param("Domain"),
            "MountOptions": self._get_param("MountOptions"),
        }

        location_uri = "/".join(["smb:/", server_hostname, subdirectory])
        arn = self.datasync_backend.create_location(
            location_uri, metadata=metadata, typ="SMB"
        )
        return json.dumps({"LocationArn": arn})

    def describe_location_smb(self) -> str:
        location_arn = self._get_param("LocationArn")
        location = self._get_location(location_arn, typ="SMB")
        return json.dumps(
            {
                "LocationArn": location.arn,
                "LocationUri": location.uri,
                "AgentArns": location.metadata["AgentArns"],
                "User": location.metadata["User"],
                "Domain": location.metadata["Domain"],
                "MountOptions": location.metadata["MountOptions"],
            }
        )

    def delete_location(self) -> str:
        location_arn = self._get_param("LocationArn")
        self.datasync_backend.delete_location(location_arn)
        return json.dumps({})

    def create_task(self) -> str:
        destination_location_arn = self._get_param("DestinationLocationArn")
        source_location_arn = self._get_param("SourceLocationArn")
        name = self._get_param("Name")
        metadata = {
            "CloudWatchLogGroupArn": self._get_param("CloudWatchLogGroupArn"),
            "Options": self._get_param("Options"),
            "Excludes": self._get_param("Excludes"),
            "Tags": self._get_param("Tags"),
        }
        arn = self.datasync_backend.create_task(
            source_location_arn, destination_location_arn, name, metadata=metadata
        )
        return json.dumps({"TaskArn": arn})

    def update_task(self) -> str:
        task_arn = self._get_param("TaskArn")
        self.datasync_backend.update_task(
            task_arn,
            name=self._get_param("Name"),
            metadata={
                "CloudWatchLogGroupArn": self._get_param("CloudWatchLogGroupArn"),
                "Options": self._get_param("Options"),
                "Excludes": self._get_param("Excludes"),
                "Tags": self._get_param("Tags"),
            },
        )
        return json.dumps({})

    def list_tasks(self) -> str:
        tasks = list()
        for arn, task in self.datasync_backend.tasks.items():
            tasks.append({"Name": task.name, "Status": task.status, "TaskArn": arn})
        return json.dumps({"Tasks": tasks})

    def delete_task(self) -> str:
        task_arn = self._get_param("TaskArn")
        self.datasync_backend.delete_task(task_arn)
        return json.dumps({})

    def describe_task(self) -> str:
        task_arn = self._get_param("TaskArn")
        task = self.datasync_backend._get_task(task_arn)
        return json.dumps(
            {
                "TaskArn": task.arn,
                "Status": task.status,
                "Name": task.name,
                "CurrentTaskExecutionArn": task.current_task_execution_arn,
                "SourceLocationArn": task.source_location_arn,
                "DestinationLocationArn": task.destination_location_arn,
                "CloudWatchLogGroupArn": task.metadata["CloudWatchLogGroupArn"],
                "Options": task.metadata["Options"],
                "Excludes": task.metadata["Excludes"],
            }
        )

    def start_task_execution(self) -> str:
        task_arn = self._get_param("TaskArn")
        arn = self.datasync_backend.start_task_execution(task_arn)
        return json.dumps({"TaskExecutionArn": arn})

    def cancel_task_execution(self) -> str:
        task_execution_arn = self._get_param("TaskExecutionArn")
        self.datasync_backend.cancel_task_execution(task_execution_arn)
        return json.dumps({})

    def describe_task_execution(self) -> str:
        task_execution_arn = self._get_param("TaskExecutionArn")
        task_execution = self.datasync_backend._get_task_execution(task_execution_arn)
        result = json.dumps(
            {"TaskExecutionArn": task_execution.arn, "Status": task_execution.status}
        )
        if task_execution.status == "SUCCESS":
            self.datasync_backend.tasks[task_execution.task_arn].status = "AVAILABLE"
        # Simulate task being executed
        task_execution.iterate_status()
        return result
