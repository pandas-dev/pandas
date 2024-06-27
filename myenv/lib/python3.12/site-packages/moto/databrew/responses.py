import json
from typing import Any, Dict
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import DataBrewBackend, databrew_backends


class DataBrewResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="databrew")

    @property
    def databrew_backend(self) -> DataBrewBackend:
        """Return backend instance specific for this region."""
        return databrew_backends[self.current_account][self.region]

    # region Recipes
    @property
    def parameters(self) -> Dict[str, Any]:  # type: ignore[misc]
        return json.loads(self.body)

    def create_recipe(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_CreateRecipe.html
        recipe_description = self.parameters.get("Description")
        recipe_steps = self.parameters.get("Steps")
        recipe_name = self.parameters.get("Name")
        tags = self.parameters.get("Tags")
        return json.dumps(
            self.databrew_backend.create_recipe(
                recipe_name,  # type: ignore[arg-type]
                recipe_description,  # type: ignore[arg-type]
                recipe_steps,  # type: ignore[arg-type]
                tags,  # type: ignore[arg-type]
            ).as_dict()
        )

    def delete_recipe_version(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_DeleteRecipeVersion.html
        split_path = self._get_path().strip("/").split("/")
        recipe_name = split_path[1]
        recipe_version = split_path[3]
        self.databrew_backend.delete_recipe_version(recipe_name, recipe_version)
        return json.dumps({"Name": recipe_name, "RecipeVersion": recipe_version})

    def _get_path(self) -> str:
        return unquote(self.parsed_url.path)

    def list_recipes(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_ListRecipes.html
        next_token = self._get_param("NextToken", self._get_param("nextToken"))
        max_results = self._get_int_param(
            "MaxResults", self._get_int_param("maxResults")
        )
        recipe_version = self._get_param(
            "RecipeVersion", self._get_param("recipeVersion")
        )

        # pylint: disable=unexpected-keyword-arg, unbalanced-tuple-unpacking
        recipe_list, next_token = self.databrew_backend.list_recipes(
            next_token=next_token,
            max_results=max_results,
            recipe_version=recipe_version,
        )
        return json.dumps(
            {
                "Recipes": [recipe.as_dict() for recipe in recipe_list],
                "NextToken": next_token,
            }
        )

    def list_recipe_versions(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_ListRecipeVersions.html
        recipe_name = self._get_param("Name", self._get_param("name"))
        next_token = self._get_param("NextToken", self._get_param("nextToken"))
        max_results = self._get_int_param(
            "MaxResults", self._get_int_param("maxResults")
        )

        # pylint: disable=unexpected-keyword-arg, unbalanced-tuple-unpacking
        recipe_list, next_token = self.databrew_backend.list_recipe_versions(
            recipe_name=recipe_name, next_token=next_token, max_results=max_results
        )
        return json.dumps(
            {
                "Recipes": [recipe.as_dict() for recipe in recipe_list],
                "NextToken": next_token,
            }
        )

    def publish_recipe(self) -> str:
        recipe_name = self._get_path().strip("/").split("/", 2)[1]
        recipe_description = self.parameters.get("Description")
        self.databrew_backend.publish_recipe(recipe_name, recipe_description)
        return json.dumps({"Name": recipe_name})

    def update_recipe(self) -> str:
        recipe_name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        recipe_description = self.parameters.get("Description")
        recipe_steps = self.parameters.get("Steps")

        self.databrew_backend.update_recipe(
            recipe_name,
            recipe_description,  # type: ignore[arg-type]
            recipe_steps,  # type: ignore[arg-type]
        )
        return json.dumps({"Name": recipe_name})

    def describe_recipe(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_DescribeRecipe.html
        recipe_name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        recipe_version = self._get_param(
            "RecipeVersion", self._get_param("recipeVersion")
        )
        recipe = self.databrew_backend.describe_recipe(
            recipe_name, recipe_version=recipe_version
        )
        return json.dumps(recipe.as_dict())

    # endregion

    # region Rulesets

    def create_ruleset(self) -> str:
        ruleset_description = self.parameters.get("Description")
        ruleset_rules = self.parameters.get("Rules")
        ruleset_name = self.parameters.get("Name")
        ruleset_target_arn = self.parameters.get("TargetArn")
        tags = self.parameters.get("Tags")

        return json.dumps(
            self.databrew_backend.create_ruleset(
                ruleset_name,  # type: ignore[arg-type]
                ruleset_description,  # type: ignore[arg-type]
                ruleset_rules,  # type: ignore[arg-type]
                ruleset_target_arn,  # type: ignore[arg-type]
                tags,  # type: ignore[arg-type]
            ).as_dict()
        )

    def update_ruleset(self) -> str:
        ruleset_name = self._get_path().split("/")[-1]
        ruleset_description = self.parameters.get("Description")
        ruleset_rules = self.parameters.get("Rules")
        tags = self.parameters.get("Tags")

        ruleset = self.databrew_backend.update_ruleset(
            ruleset_name,
            ruleset_description,  # type: ignore[arg-type]
            ruleset_rules,  # type: ignore[arg-type]
            tags,  # type: ignore[arg-type]
        )
        return json.dumps(ruleset.as_dict())

    def describe_ruleset(self) -> str:
        ruleset_name = self._get_path().split("/")[-1]
        ruleset = self.databrew_backend.describe_ruleset(ruleset_name)
        return json.dumps(ruleset.as_dict())

    def delete_ruleset(self) -> str:
        ruleset_name = self._get_path().split("/")[-1]
        self.databrew_backend.delete_ruleset(ruleset_name)
        return json.dumps({"Name": ruleset_name})

    def list_rulesets(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_ListRulesets.html
        next_token = self._get_param("NextToken", self._get_param("nextToken"))
        max_results = self._get_int_param(
            "MaxResults", self._get_int_param("maxResults")
        )

        # pylint: disable=unexpected-keyword-arg, unbalanced-tuple-unpacking
        ruleset_list, next_token = self.databrew_backend.list_rulesets(
            next_token=next_token, max_results=max_results
        )
        return json.dumps(
            {
                "Rulesets": [ruleset.as_dict() for ruleset in ruleset_list],
                "NextToken": next_token,
            }
        )

    # endregion

    # region Datasets

    def create_dataset(self) -> str:
        dataset_name = self.parameters.get("Name")
        dataset_format = self.parameters.get("Format")
        dataset_format_options = self.parameters.get("FormatOptions")
        dataset_input = self.parameters.get("Input")
        dataset_path_otions = self.parameters.get("PathOptions")
        dataset_tags = self.parameters.get("Tags")

        return json.dumps(
            self.databrew_backend.create_dataset(
                dataset_name,  # type: ignore[arg-type]
                dataset_format,  # type: ignore[arg-type]
                dataset_format_options,  # type: ignore[arg-type]
                dataset_input,  # type: ignore[arg-type]
                dataset_path_otions,  # type: ignore[arg-type]
                dataset_tags,  # type: ignore[arg-type]
            ).as_dict()
        )

    def list_datasets(self) -> str:
        next_token = self._get_param("NextToken", self._get_param("nextToken"))
        max_results = self._get_int_param(
            "MaxResults", self._get_int_param("maxResults")
        )

        # pylint: disable=unexpected-keyword-arg, unbalanced-tuple-unpacking
        dataset_list, next_token = self.databrew_backend.list_datasets(
            next_token=next_token, max_results=max_results
        )

        return json.dumps(
            {
                "Datasets": [dataset.as_dict() for dataset in dataset_list],
                "NextToken": next_token,
            }
        )

    def update_dataset(self) -> str:
        dataset_name = self._get_path().split("/")[-1]
        dataset_format = self.parameters.get("Format")
        dataset_format_options = self.parameters.get("FormatOptions")
        dataset_input = self.parameters.get("Input")
        dataset_path_otions = self.parameters.get("PathOptions")
        dataset_tags = self.parameters.get("Tags")

        dataset = self.databrew_backend.update_dataset(
            dataset_name,
            dataset_format,  # type: ignore[arg-type]
            dataset_format_options,  # type: ignore[arg-type]
            dataset_input,  # type: ignore[arg-type]
            dataset_path_otions,  # type: ignore[arg-type]
            dataset_tags,  # type: ignore[arg-type]
        )
        return json.dumps(dataset.as_dict())

    def delete_dataset(self) -> str:
        dataset_name = self._get_path().split("/")[-1]
        self.databrew_backend.delete_dataset(dataset_name)
        return json.dumps({"Name": dataset_name})

    def describe_dataset(self) -> str:
        dataset_name = self._get_path().split("/")[-1]
        dataset = self.databrew_backend.describe_dataset(dataset_name)
        return json.dumps(dataset.as_dict())

    # endregion

    # region Jobs
    def list_jobs(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_ListJobs.html
        dataset_name = self._get_param("datasetName")
        project_name = self._get_param("projectName")
        next_token = self._get_param("NextToken", self._get_param("nextToken"))
        max_results = self._get_int_param(
            "MaxResults", self._get_int_param("maxResults")
        )

        job_list, next_token = self.databrew_backend.list_jobs(
            dataset_name=dataset_name,
            project_name=project_name,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(
            {
                "Jobs": [job.as_dict() for job in job_list],
                "NextToken": next_token,
            }
        )

    def describe_job(self) -> str:
        job_name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        job = self.databrew_backend.describe_job(job_name)
        return json.dumps(job.as_dict())

    def delete_job(self) -> str:
        job_name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        self.databrew_backend.delete_job(job_name)
        return json.dumps({"Name": job_name})

    def create_profile_job(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_CreateProfileJob.html
        kwargs = {
            "dataset_name": self._get_param("DatasetName"),
            "name": self._get_param("Name"),
            "output_location": self._get_param("OutputLocation"),
            "role_arn": self._get_param("RoleArn"),
            "configuration": self._get_param("Configuration"),
            "encryption_key_arn": self._get_param("EncryptionKeyArn"),
            "encryption_mode": self._get_param("EncryptionMode"),
            "job_sample": self._get_param("JobSample"),
            "log_subscription": self._get_param("LogSubscription"),
            "max_capacity": self._get_int_param("MaxCapacity"),
            "max_retries": self._get_int_param("MaxRetries"),
            "tags": self._get_param("Tags"),
            "timeout": self._get_int_param("Timeout"),
            "validation_configurations": self._get_param("ValidationConfigurations"),
        }
        return json.dumps(self.databrew_backend.create_profile_job(**kwargs).as_dict())

    def update_profile_job(self) -> str:
        name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        # https://docs.aws.amazon.com/databrew/latest/dg/API_UpdateProfileJob.html
        kwargs = {
            "name": name,
            "output_location": self._get_param("OutputLocation"),
            "role_arn": self._get_param("RoleArn"),
            "configuration": self._get_param("Configuration"),
            "encryption_key_arn": self._get_param("EncryptionKeyArn"),
            "encryption_mode": self._get_param("EncryptionMode"),
            "job_sample": self._get_param("JobSample"),
            "log_subscription": self._get_param("LogSubscription"),
            "max_capacity": self._get_int_param("MaxCapacity"),
            "max_retries": self._get_int_param("MaxRetries"),
            "timeout": self._get_int_param("Timeout"),
            "validation_configurations": self._get_param("ValidationConfigurations"),
        }
        return json.dumps(self.databrew_backend.update_profile_job(**kwargs).as_dict())

    def create_recipe_job(self) -> str:
        # https://docs.aws.amazon.com/databrew/latest/dg/API_CreateRecipeJob.html
        kwargs = {
            "name": self._get_param("Name"),
            "role_arn": self._get_param("RoleArn"),
            "database_outputs": self._get_param("DatabaseOutputs"),
            "data_catalog_outputs": self._get_param("DataCatalogOutputs"),
            "dataset_name": self._get_param("DatasetName"),
            "encryption_key_arn": self._get_param("EncryptionKeyArn"),
            "encryption_mode": self._get_param("EncryptionMode"),
            "log_subscription": self._get_param("LogSubscription"),
            "max_capacity": self._get_int_param("MaxCapacity"),
            "max_retries": self._get_int_param("MaxRetries"),
            "outputs": self._get_param("Outputs"),
            "project_name": self._get_param("ProjectName"),
            "recipe_reference": self._get_param("RecipeReference"),
            "tags": self._get_param("Tags"),
            "timeout": self._get_int_param("Timeout"),
        }
        return json.dumps(self.databrew_backend.create_recipe_job(**kwargs).as_dict())

    def update_recipe_job(self) -> str:
        name = self._get_path().rstrip("/").rsplit("/", 1)[1]
        # https://docs.aws.amazon.com/databrew/latest/dg/API_UpdateRecipeJob.html
        kwargs = {
            "name": name,
            "role_arn": self._get_param("RoleArn"),
            "database_outputs": self._get_param("DatabaseOutputs"),
            "data_catalog_outputs": self._get_param("DataCatalogOutputs"),
            "encryption_key_arn": self._get_param("EncryptionKeyArn"),
            "encryption_mode": self._get_param("EncryptionMode"),
            "log_subscription": self._get_param("LogSubscription"),
            "max_capacity": self._get_int_param("MaxCapacity"),
            "max_retries": self._get_int_param("MaxRetries"),
            "outputs": self._get_param("Outputs"),
            "timeout": self._get_int_param("Timeout"),
        }
        return json.dumps(self.databrew_backend.update_recipe_job(**kwargs).as_dict())

    # endregion
