from moto.core.responses import BaseResponse
from moto.core.utils import tags_from_query_string

from .exceptions import InvalidParameterValueError
from .models import EBBackend, eb_backends


class EBResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="elasticbeanstalk")

    @property
    def elasticbeanstalk_backend(self) -> EBBackend:
        """
        :rtype: EBBackend
        """
        return eb_backends[self.current_account][self.region]

    def create_application(self) -> str:
        app = self.elasticbeanstalk_backend.create_application(
            application_name=self._get_param("ApplicationName")
        )

        template = self.response_template(EB_CREATE_APPLICATION)
        return template.render(
            region_name=self.elasticbeanstalk_backend.region_name, application=app
        )

    def describe_applications(self) -> str:
        template = self.response_template(EB_DESCRIBE_APPLICATIONS)
        return template.render(
            applications=self.elasticbeanstalk_backend.applications.values()
        )

    def create_environment(self) -> str:
        application_name = self._get_param("ApplicationName")
        try:
            app = self.elasticbeanstalk_backend.applications[application_name]
        except KeyError:
            raise InvalidParameterValueError(
                f"No Application named '{application_name}' found."
            )

        tags = tags_from_query_string(self.querystring, prefix="Tags.member")
        env = self.elasticbeanstalk_backend.create_environment(
            app,
            environment_name=self._get_param("EnvironmentName"),
            stack_name=self._get_param("SolutionStackName"),
            tags=tags,
        )

        template = self.response_template(EB_CREATE_ENVIRONMENT)
        return template.render(
            environment=env, region=self.elasticbeanstalk_backend.region_name
        )

    def describe_environments(self) -> str:
        envs = self.elasticbeanstalk_backend.describe_environments()

        template = self.response_template(EB_DESCRIBE_ENVIRONMENTS)
        return template.render(environments=envs)

    def list_available_solution_stacks(self) -> str:
        return EB_LIST_AVAILABLE_SOLUTION_STACKS

    def update_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        tags_to_add = tags_from_query_string(
            self.querystring, prefix="TagsToAdd.member"
        )
        tags_to_remove = self._get_multi_param("TagsToRemove.member")
        self.elasticbeanstalk_backend.update_tags_for_resource(
            resource_arn, tags_to_add, tags_to_remove
        )

        return EB_UPDATE_TAGS_FOR_RESOURCE

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        tags = self.elasticbeanstalk_backend.list_tags_for_resource(resource_arn)

        template = self.response_template(EB_LIST_TAGS_FOR_RESOURCE)
        return template.render(tags=tags, arn=resource_arn)

    def delete_application(self) -> str:
        application_name = self._get_param("ApplicationName")
        self.elasticbeanstalk_backend.delete_application(
            application_name=application_name,
        )
        return DELETE_APPLICATION_TEMPLATE


EB_CREATE_APPLICATION = """
<CreateApplicationResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <CreateApplicationResult>
    <Application>
      <ConfigurationTemplates/>
      <DateCreated>2019-09-03T13:08:29.049Z</DateCreated>
      <ResourceLifecycleConfig>
        <VersionLifecycleConfig>
          <MaxAgeRule>
            <DeleteSourceFromS3>false</DeleteSourceFromS3>
            <MaxAgeInDays>180</MaxAgeInDays>
            <Enabled>false</Enabled>
          </MaxAgeRule>
          <MaxCountRule>
            <DeleteSourceFromS3>false</DeleteSourceFromS3>
            <MaxCount>200</MaxCount>
            <Enabled>false</Enabled>
          </MaxCountRule>
        </VersionLifecycleConfig>
      </ResourceLifecycleConfig>
      <ApplicationArn>{{ application.arn }}</ApplicationArn>
      <ApplicationName>{{ application.application_name }}</ApplicationName>
      <DateUpdated>2019-09-03T13:08:29.049Z</DateUpdated>
    </Application>
  </CreateApplicationResult>
  <ResponseMetadata>
    <RequestId>1b6173c8-13aa-4b0a-99e9-eb36a1fb2778</RequestId>
  </ResponseMetadata>
</CreateApplicationResponse>
"""

EB_DESCRIBE_APPLICATIONS = """
<DescribeApplicationsResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <DescribeApplicationsResult>
    <Applications>
      {% for application in applications %}
      <member>
        <ConfigurationTemplates/>
        <DateCreated>2019-09-03T13:08:29.049Z</DateCreated>
        <ResourceLifecycleConfig>
          <VersionLifecycleConfig>
            <MaxAgeRule>
              <MaxAgeInDays>180</MaxAgeInDays>
              <DeleteSourceFromS3>false</DeleteSourceFromS3>
              <Enabled>false</Enabled>
            </MaxAgeRule>
            <MaxCountRule>
              <DeleteSourceFromS3>false</DeleteSourceFromS3>
              <MaxCount>200</MaxCount>
              <Enabled>false</Enabled>
            </MaxCountRule>
          </VersionLifecycleConfig>
        </ResourceLifecycleConfig>
        <ApplicationArn>{{ application.arn }}</ApplicationArn>
        <ApplicationName>{{ application.application_name }}</ApplicationName>
        <DateUpdated>2019-09-03T13:08:29.049Z</DateUpdated>
      </member>
      {% endfor %}
    </Applications>
  </DescribeApplicationsResult>
  <ResponseMetadata>
    <RequestId>015a05eb-282e-4b76-bd18-663fdfaf42e4</RequestId>
  </ResponseMetadata>
</DescribeApplicationsResponse>
"""

DELETE_APPLICATION_TEMPLATE = """
<DeleteApplicationResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <ResponseMetadata>
    <RequestId>015a05eb-282e-4b76-bd18-663fdfaf42e4</RequestId>
  </ResponseMetadata>
</DeleteApplicationResponse>
"""

EB_CREATE_ENVIRONMENT = """
<CreateEnvironmentResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <CreateEnvironmentResult>
    <SolutionStackName>{{ environment.solution_stack_name }}</SolutionStackName>
    <Health>Grey</Health>
    <EnvironmentArn>{{ environment.environment_arn }}</EnvironmentArn>
    <DateUpdated>2019-09-04T09:41:24.222Z</DateUpdated>
    <DateCreated>2019-09-04T09:41:24.222Z</DateCreated>
    <EnvironmentId>{{ environment_id }}</EnvironmentId>
    <PlatformArn>{{ environment.platform_arn }}</PlatformArn>
    <Tier>
      <Name>WebServer</Name>
      <Type>Standard</Type>
      <Version>1.0</Version>
    </Tier>
    <EnvironmentName>{{ environment.environment_name }}</EnvironmentName>
    <ApplicationName>{{ environment.application_name }}</ApplicationName>
    <Status>Launching</Status>
  </CreateEnvironmentResult>
  <ResponseMetadata>
    <RequestId>18dc8158-f5d7-4d5a-82ef-07fcaadf81c6</RequestId>
  </ResponseMetadata>
</CreateEnvironmentResponse>
"""

EB_DESCRIBE_ENVIRONMENTS = """
<DescribeEnvironmentsResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <DescribeEnvironmentsResult>
    <Environments>
      {% for env in environments %}
      <member>
        <SolutionStackName>{{ env.solution_stack_name }}</SolutionStackName>
        <Health>Grey</Health>
        <EnvironmentArn>{{ env.environment_arn }}</EnvironmentArn>
        <MinCapacityEnabled>false</MinCapacityEnabled>
        <DateUpdated>2019-08-30T09:35:10.913Z</DateUpdated>
        <AbortableOperationInProgress>false</AbortableOperationInProgress>
        <Alerts/>
        <DateCreated>2019-08-22T07:02:47.332Z</DateCreated>
        <EnvironmentId>{{ env.environment_id }}</EnvironmentId>
        <VersionLabel>1</VersionLabel>
        <PlatformArn>{{ env.platform_arn }}</PlatformArn>
        <Tier>
          <Name>WebServer</Name>
          <Type>Standard</Type>
          <Version>1.0</Version>
        </Tier>
        <HealthStatus>No Data</HealthStatus>
        <EnvironmentName>{{ env.environment_name }}</EnvironmentName>
        <EndpointURL></EndpointURL>
        <CNAME></CNAME>
        <EnvironmentLinks/>
        <ApplicationName>{{ env.application_name }}</ApplicationName>
        <Status>Ready</Status>
      </member>
      {% endfor %}
    </Environments>
  </DescribeEnvironmentsResult>
  <ResponseMetadata>
    <RequestId>dd56b215-01a0-40b2-bd1e-57589c39424f</RequestId>
  </ResponseMetadata>
</DescribeEnvironmentsResponse>
"""

# Current list as of 2019-09-04
EB_LIST_AVAILABLE_SOLUTION_STACKS = """
<ListAvailableSolutionStacksResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <ListAvailableSolutionStacksResult>
    <SolutionStacks>
      <member>64bit Amazon Linux 2018.03 v4.10.1 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.9.2 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.8.0 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.6.0 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.5.3 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.5.1 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v4.5.0 running Node.js</member>
      <member>64bit Amazon Linux 2017.09 v4.4.6 running Node.js</member>
      <member>64bit Amazon Linux 2017.09 v4.4.5 running Node.js</member>
      <member>64bit Amazon Linux 2017.09 v4.4.4 running Node.js</member>
      <member>64bit Amazon Linux 2017.09 v4.4.2 running Node.js</member>
      <member>64bit Amazon Linux 2017.09 v4.4.0 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.3.0 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.2.2 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.2.1 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.2.0 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.1.1 running Node.js</member>
      <member>64bit Amazon Linux 2017.03 v4.1.0 running Node.js</member>
      <member>64bit Amazon Linux 2016.09 v4.0.1 running Node.js</member>
      <member>64bit Amazon Linux 2016.09 v4.0.0 running Node.js</member>
      <member>64bit Amazon Linux 2016.09 v3.3.1 running Node.js</member>
      <member>64bit Amazon Linux 2016.09 v3.1.0 running Node.js</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.4</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.5</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.6</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.0</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.1</member>
      <member>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.12 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.7 running PHP 7.1</member>
      <member>64bit Amazon Linux 2018.03 v2.8.6 running PHP 7.1</member>
      <member>64bit Amazon Linux 2018.03 v2.8.6 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.5 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.4 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.3 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.2 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.1 running PHP 7.2</member>
      <member>64bit Amazon Linux 2018.03 v2.8.0 running PHP 7.1</member>
      <member>64bit Amazon Linux 2018.03 v2.7.1 running PHP 5.6</member>
      <member>64bit Amazon Linux 2018.03 v2.7.1 running PHP 7.0</member>
      <member>64bit Amazon Linux 2018.03 v2.7.1 running PHP 7.1</member>
      <member>64bit Amazon Linux 2018.03 v2.7.0 running PHP 7.0</member>
      <member>64bit Amazon Linux 2018.03 v2.7.0 running PHP 7.1</member>
      <member>64bit Amazon Linux 2017.09 v2.6.6 running PHP 5.4</member>
      <member>64bit Amazon Linux 2017.09 v2.6.6 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.6 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.5 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.4</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.5</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running PHP 7.1</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.4</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.5</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running PHP 7.1</member>
      <member>64bit Amazon Linux 2017.09 v2.6.2 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.2 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.1 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.4</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.5</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running PHP 7.1</member>
      <member>64bit Amazon Linux 2017.03 v2.5.0 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.03 v2.5.0 running PHP 7.1</member>
      <member>64bit Amazon Linux 2017.03 v2.4.4 running PHP 5.5</member>
      <member>64bit Amazon Linux 2017.03 v2.4.4 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.03 v2.4.4 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.03 v2.4.3 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.4</member>
      <member>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.5</member>
      <member>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.6</member>
      <member>64bit Amazon Linux 2017.03 v2.4.2 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.03 v2.4.1 running PHP 7.0</member>
      <member>64bit Amazon Linux 2017.03 v2.4.0 running PHP 7.0</member>
      <member>64bit Amazon Linux 2016.09 v2.3.2 running PHP 7.0</member>
      <member>64bit Amazon Linux 2016.09 v2.3.1 running PHP 7.0</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Python 3.6</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Python 3.4</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Python</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Python 2.7</member>
      <member>64bit Amazon Linux 2018.03 v2.7.5 running Python 3.6</member>
      <member>64bit Amazon Linux 2018.03 v2.7.1 running Python 3.6</member>
      <member>64bit Amazon Linux 2018.03 v2.7.0 running Python 3.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running Python 3.6</member>
      <member>64bit Amazon Linux 2017.09 v2.6.1 running Python 3.6</member>
      <member>64bit Amazon Linux 2017.03 v2.4.0 running Python 3.4</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.6 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.5 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.4 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.3 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.2 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.1 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.0 (Puma)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.6 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.5 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.4 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.3 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.2 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.1 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.0 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 1.9.3</member>
      <member>64bit Amazon Linux 2018.03 v2.8.0 running Ruby 2.5 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2017.03 v2.4.4 running Ruby 2.3 (Puma)</member>
      <member>64bit Amazon Linux 2017.03 v2.4.4 running Ruby 2.3 (Passenger Standalone)</member>
      <member>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 8.5 Java 8</member>
      <member>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 7 Java 7</member>
      <member>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 7 Java 6</member>
      <member>64bit Amazon Linux 2018.03 v3.1.1 running Tomcat 8.5 Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.6.5 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.6.2 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.6.1 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.6.0 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2016.09 v2.5.4 running Tomcat 8 Java 8</member>
      <member>64bit Amazon Linux 2016.03 v2.1.0 running Tomcat 8 Java 8</member>
      <member>64bit Windows Server Core 2016 v2.2.1 running IIS 10.0</member>
      <member>64bit Windows Server 2016 v2.2.1 running IIS 10.0</member>
      <member>64bit Windows Server Core 2012 R2 v2.2.1 running IIS 8.5</member>
      <member>64bit Windows Server 2012 R2 v2.2.1 running IIS 8.5</member>
      <member>64bit Windows Server Core 2016 v1.2.0 running IIS 10.0</member>
      <member>64bit Windows Server 2016 v1.2.0 running IIS 10.0</member>
      <member>64bit Windows Server Core 2012 R2 v1.2.0 running IIS 8.5</member>
      <member>64bit Windows Server 2012 R2 v1.2.0 running IIS 8.5</member>
      <member>64bit Windows Server 2012 v1.2.0 running IIS 8</member>
      <member>64bit Windows Server 2008 R2 v1.2.0 running IIS 7.5</member>
      <member>64bit Windows Server Core 2012 R2 running IIS 8.5</member>
      <member>64bit Windows Server 2012 R2 running IIS 8.5</member>
      <member>64bit Windows Server 2012 running IIS 8</member>
      <member>64bit Windows Server 2008 R2 running IIS 7.5</member>
      <member>64bit Amazon Linux 2018.03 v2.12.16 running Docker 18.06.1-ce</member>
      <member>64bit Amazon Linux 2016.09 v2.5.2 running Docker 1.12.6</member>
      <member>64bit Amazon Linux 2018.03 v2.15.2 running Multi-container Docker 18.06.1-ce (Generic)</member>
      <member>64bit Debian jessie v2.12.16 running Go 1.4 (Preconfigured - Docker)</member>
      <member>64bit Debian jessie v2.12.16 running Go 1.3 (Preconfigured - Docker)</member>
      <member>64bit Debian jessie v2.12.16 running Python 3.4 (Preconfigured - Docker)</member>
      <member>64bit Debian jessie v2.10.0 running Python 3.4 (Preconfigured - Docker)</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.9.1 running Java 7</member>
      <member>64bit Amazon Linux 2018.03 v2.8.0 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.7.6 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.7.5 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.7.4 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.7.2 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.7.1 running Java 8</member>
      <member>64bit Amazon Linux 2017.09 v2.6.8 running Java 8</member>
      <member>64bit Amazon Linux 2017.09 v2.6.5 running Java 8</member>
      <member>64bit Amazon Linux 2017.09 v2.6.4 running Java 8</member>
      <member>64bit Amazon Linux 2017.09 v2.6.3 running Java 8</member>
      <member>64bit Amazon Linux 2017.09 v2.6.0 running Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.5.4 running Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.5.3 running Java 8</member>
      <member>64bit Amazon Linux 2017.03 v2.5.2 running Java 8</member>
      <member>64bit Amazon Linux 2016.09 v2.4.4 running Java 8</member>
      <member>64bit Amazon Linux 2018.03 v2.12.1 running Go 1.12.7</member>
      <member>64bit Amazon Linux 2018.03 v2.6.14 running Packer 1.0.3</member>
      <member>64bit Amazon Linux 2018.03 v2.12.16 running GlassFish 5.0 Java 8 (Preconfigured - Docker)</member>
    </SolutionStacks>
    <SolutionStackDetails>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.10.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.9.2 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.8.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.6.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.5.3 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.5.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v4.5.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v4.4.6 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v4.4.5 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v4.4.4 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v4.4.2 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v4.4.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.3.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.2.2 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.2.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.2.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.1.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v4.1.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v4.0.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v4.0.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v3.3.1 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v3.1.0 running Node.js</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.14 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.12 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.7 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.6 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.6 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.5 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.4 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.3 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.2 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.1 running PHP 7.2</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.0 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.1 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.1 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.1 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.0 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.0 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.6 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.6 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.6 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.5 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.2 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.2 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.1 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.5.0 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.5.0 running PHP 7.1</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.4 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.4 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.4 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.3 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.2 running PHP 5.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.2 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.1 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.0 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v2.3.2 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v2.3.1 running PHP 7.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Python 3.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Python</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Python 2.7</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.5 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.1 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.0 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.1 running Python 3.6</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.0 running Python 3.4</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.6 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.5 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.4 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.3 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.2 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.1 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.0 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.6 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.5 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.4 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.3 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.2 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.1 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 2.0 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.10.1 running Ruby 1.9.3</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.0 running Ruby 2.5 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.4 running Ruby 2.3 (Puma)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.4.4 running Ruby 2.3 (Passenger Standalone)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 8.5 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 7 Java 7</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v3.2.1 running Tomcat 7 Java 6</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v3.1.1 running Tomcat 8.5 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.6.5 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.6.2 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.6.1 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.6.0 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v2.5.4 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.03 v2.1.0 running Tomcat 8 Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>war</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server Core 2016 v2.2.1 running IIS 10.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2016 v2.2.1 running IIS 10.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server Core 2012 R2 v2.2.1 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2012 R2 v2.2.1 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server Core 2016 v1.2.0 running IIS 10.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2016 v1.2.0 running IIS 10.0</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server Core 2012 R2 v1.2.0 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2012 R2 v1.2.0 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2012 v1.2.0 running IIS 8</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2008 R2 v1.2.0 running IIS 7.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server Core 2012 R2 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2012 R2 running IIS 8.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2012 running IIS 8</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Windows Server 2008 R2 running IIS 7.5</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.12.16 running Docker 18.06.1-ce</SolutionStackName>
        <PermittedFileTypes/>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v2.5.2 running Docker 1.12.6</SolutionStackName>
        <PermittedFileTypes/>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.15.2 running Multi-container Docker 18.06.1-ce (Generic)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
          <member>json</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Debian jessie v2.12.16 running Go 1.4 (Preconfigured - Docker)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Debian jessie v2.12.16 running Go 1.3 (Preconfigured - Docker)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Debian jessie v2.12.16 running Python 3.4 (Preconfigured - Docker)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Debian jessie v2.10.0 running Python 3.4 (Preconfigured - Docker)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.9.1 running Java 7</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.8.0 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.6 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.5 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.4 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.2 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.7.1 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.8 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.5 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.4 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.3 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.09 v2.6.0 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.5.4 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.5.3 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2017.03 v2.5.2 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2016.09 v2.4.4 running Java 8</SolutionStackName>
        <PermittedFileTypes>
          <member>jar</member>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.12.1 running Go 1.12.7</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.6.14 running Packer 1.0.3</SolutionStackName>
        <PermittedFileTypes/>
      </member>
      <member>
        <SolutionStackName>64bit Amazon Linux 2018.03 v2.12.16 running GlassFish 5.0 Java 8 (Preconfigured - Docker)</SolutionStackName>
        <PermittedFileTypes>
          <member>zip</member>
        </PermittedFileTypes>
      </member>
    </SolutionStackDetails>
  </ListAvailableSolutionStacksResult>
  <ResponseMetadata>
    <RequestId>bd6bd2b2-9983-4845-b53b-fe53e8a5e1e7</RequestId>
  </ResponseMetadata>
</ListAvailableSolutionStacksResponse>
"""

EB_UPDATE_TAGS_FOR_RESOURCE = """
<UpdateTagsForResourceResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <ResponseMetadata>
    <RequestId>f355d788-e67e-440f-b915-99e35254ffee</RequestId>
  </ResponseMetadata>
</UpdateTagsForResourceResponse>
"""

EB_LIST_TAGS_FOR_RESOURCE = """
<ListTagsForResourceResponse xmlns="http://elasticbeanstalk.amazonaws.com/docs/2010-12-01/">
  <ListTagsForResourceResult>
    <ResourceTags>
      {% for key, value in tags.items() %}
      <member>
        <Key>{{ key }}</Key>
        <Value>{{ value }}</Value>
      </member>
      {% endfor %}
    </ResourceTags>
    <ResourceArn>{{ arn }}</ResourceArn>
  </ListTagsForResourceResult>
  <ResponseMetadata>
    <RequestId>178e410f-3b57-456f-a64c-a3b6a16da9ab</RequestId>
  </ResponseMetadata>
</ListTagsForResourceResponse>
"""
