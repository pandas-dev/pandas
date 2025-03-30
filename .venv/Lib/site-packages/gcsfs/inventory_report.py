from datetime import datetime


class InventoryReport:
    """
    A utility class for fetching and processing inventory reports from GCS.

    The 'InventoryReport' class provides logic to support logic to fetch
    inventory reports, and process their content to obtain a final snapshot
    of objects in the latest inventory reports.

    High-Level Functionality:
    ------------------------
    1. Fetching Inventory Reports:
       - The class offers methods to fetch inventory report configurations and
         metadata from GCS.
       - It validates the inventory report information provided by the user.
       - Inventory report configurations include options for parsing CSV format
         and specifying the bucket and destination path.

    2. Parsing and Processing Inventory Report Content:
       - The class processes the raw content of inventory reports to extract
         object details such as name, size, etc.
       - It supports listing objects using a snapshot option or filtering
         based on a user-defined prefix.
       - The class handles CSV parsing, removes header (if specified), and
         fetches required object metadata.

    3. Constructing the Final Snapshot:
       - If the user wishes to use the snapshot to do listing directly, the
         snapshot will contain the relevant object details and subdirectory
         prefixes, filtered by the prefix.

       - If the user wishes to use the snapshot as a starting point for async
         listing, the snapshot will only contain a list of object names,
         filtered by the prefix.

    Note:
    -----
    - The class should only be internally used in the 'GCSFileSystem' as an
      optional configuration during listing.

    Example Usage:
    --------------
    # Should already be instanted in 'core.py'
    gcs_file_system = GCSFileSystem(...)

    # User defines inventory report information
    inventory_report_info = {
        "use_snapshot_listing": True,
        "location": "us-east1",
        "id": "inventory_report_id"
    }

    # User defines a prefix for filtering objects
    prefix = "prefix/"

    # Fetch the snapshot based on inventory reports
    items, prefixes = await InventoryReport.fetch_snapshot(
    gcs_file_system, inventory_report_info, prefix)
    """

    # HTTP endpoint of the Storage Insights Service.
    BASE_URL = "https://storageinsights.googleapis.com/v1"

    @classmethod
    async def fetch_snapshot(cls, gcs_file_system, inventory_report_info, prefix):
        """
        Main entry point of the 'InventoryReport' class.
        Fetches the latest snapshot of objects based on inventory report configuration.

        Parameters:
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            inventory_report_info (dict): A client-configured dictionary
            containing inventory report information.
            prefix (str): Listing prefix specified by the client.

        Returns:
            tuple: A tuple containing two lists: the 'items' list representing
            object details for the snapshot, and the 'prefixes' list containing
            subdirectory prefixes.

            Note: when 'use_snapshot_listing' in 'inventory_report_info' is set
            to False, the 'prefixes' list will be empty, and the 'items' list
            will contain only the object names.
        """
        # Validate the inventory report info that the user passes in.
        cls._validate_inventory_report_info(inventory_report_info)

        # Parse the inventory report info.
        use_snapshot_listing = inventory_report_info.get("use_snapshot_listing")
        inventory_report_location = inventory_report_info.get("location")
        inventory_report_id = inventory_report_info.get("id")

        # Fetch the inventory report configuration.
        raw_inventory_report_config = await cls._fetch_raw_inventory_report_config(
            gcs_file_system=gcs_file_system,
            location=inventory_report_location,
            id=inventory_report_id,
        )

        # Parse the inventory report configuration.
        inventory_report_config = cls._parse_raw_inventory_report_config(
            raw_inventory_report_config=raw_inventory_report_config,
            use_snapshot_listing=use_snapshot_listing,
        )

        # Use the config to fetch all inventory report metadata.
        unsorted_inventory_report_metadata = await cls._fetch_inventory_report_metadata(
            gcs_file_system=gcs_file_system,
            inventory_report_config=inventory_report_config,
        )

        # Sort the metadata based on reverse created time order.
        inventory_report_metadata = cls._sort_inventory_report_metadata(
            unsorted_inventory_report_metadata=unsorted_inventory_report_metadata
        )

        # Download the most recent inventory reports in raw form.
        bucket = inventory_report_config.bucket
        inventory_report_content = await cls._download_inventory_report_content(
            gcs_file_system=gcs_file_system,
            inventory_report_metadata=inventory_report_metadata,
            bucket=bucket,
        )

        # Parse the raw inventory reports into snapshot objects.
        objects = cls._parse_inventory_report_content(
            gcs_file_system=gcs_file_system,
            inventory_report_content=inventory_report_content,
            inventory_report_config=inventory_report_config,
            use_snapshot_listing=use_snapshot_listing,
            bucket=bucket,
        )

        # Construct the final snapshot based on the fetched objects.
        snapshot = cls._construct_final_snapshot(
            objects=objects, prefix=prefix, use_snapshot_listing=use_snapshot_listing
        )

        # Return the final snapshot.
        return snapshot

    def _validate_inventory_report_info(inventory_report_info):
        """
        Validates the inventory report information dictionary that user
        passes in.

        Parameters:
            inventory_report_info (dict): A dictionary containing the inventory
            report information with the following keys:
                - "use_snapshot_listing" (bool): A flag indicating whether
                  to use snapshot listing in the inventory report.
                - "location" (str): The location of the inventory report in GCS.
                - "id" (str): The ID of the inventory report in GCS.

        Raises:
            ValueError: If any required key (use_snapshot_listing, location, id)
            is missing from the inventory_report_info dictionary.
        """
        if "use_snapshot_listing" not in inventory_report_info:
            raise ValueError("Use snapshot listing is not configured.")
        if "location" not in inventory_report_info:
            raise ValueError("Inventory report location is not configured.")
        if "id" not in inventory_report_info:
            raise ValueError("Inventory report id is not configured.")

    async def _fetch_raw_inventory_report_config(gcs_file_system, location, id):
        """
        Fetches the raw inventory report configuration from GCS based on the
        specified location and ID.

        Parameters:
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            location (str): The location of the inventory report in GCS.
            id (str): The ID of the inventory report in GCS.

        Returns:
            dict: A dictionary containing the raw inventory report
            configuration retrieved from GCS.

        Raises:
            Exception: If there is an error while fetching the inventory
            report configuration.
        """
        project = gcs_file_system.project
        url = "{}/projects/{}/locations/{}/reportConfigs/{}"
        url = url.format(InventoryReport.BASE_URL, project, location, id)
        try:
            raw_inventory_report_config = await gcs_file_system._call(
                "GET", url, json_out=True
            )
            return raw_inventory_report_config
        except Exception as e:
            raise ValueError(
                f"Error encountered when fetching inventory report config: {e}."
            )

    def _parse_raw_inventory_report_config(
        raw_inventory_report_config, use_snapshot_listing
    ):
        """
        Parses the raw inventory report configuration and validates its properties.

        Parameters:
            raw_inventory_report_config (dict): A dictionary containing the raw
            inventory report configuration retrieved from GCS.
            use_snapshot_listing (bool): A flag indicating whether to use snapshot
            listing in the inventory report.

        Returns:
            InventoryReportConfig: An instance of the InventoryReportConfig
            class representing the parsed inventory report configuration.

        Raises:
            ValueError: If the current date is outside the start and
            end range specified in the inventory report config.
            ValueError: If the "name" field is not present in the metadata
            fields of the report config.
            ValueError: If "size" field is not present in the metadata
            fields and use_snapshot_listing is True.
        """
        # Parse the report config.
        frequency_options = raw_inventory_report_config.get("frequencyOptions")
        start_date = InventoryReport._convert_obj_to_date(
            frequency_options.get("startDate")
        )
        end_date = InventoryReport._convert_obj_to_date(
            frequency_options.get("endDate")
        )
        object_metadata_report_options = raw_inventory_report_config.get(
            "objectMetadataReportOptions"
        )
        storage_destination_options = object_metadata_report_options.get(
            "storageDestinationOptions"
        )

        # Save relevant report config properties.
        csv_options = raw_inventory_report_config.get("csvOptions")
        bucket = storage_destination_options.get("bucket")
        destination_path = storage_destination_options.get("destinationPath")
        metadata_fields = object_metadata_report_options.get("metadataFields")

        # Validate date, making sure the current date is within the start and end range.
        today = datetime.now()
        if today < start_date or today > end_date:
            raise ValueError(
                f"Current date {today} is outside the range \
                {start_date} and {end_date} specified by the inventory report config."
            )

        # Validate object name exists in the metadata fields.
        # Note that the size field is mandated to be included in the
        # config when the client sets up the inventory report.
        obj_name_idx = metadata_fields.index("name")

        # If the user wants to do listing based on the snapshot, also
        # validate the report contains size metadata for each object.
        if use_snapshot_listing:
            try:
                metadata_fields.index("size")
            except ValueError:
                raise ValueError(
                    "If you want to use the snapshot for listing, the object size \
                        metadata has to be included in the inventory report."
                )

        # Finally, construct and return the inventory report config.
        inventory_report_config = InventoryReportConfig(
            csv_options=csv_options,
            bucket=bucket,
            destination_path=destination_path,
            metadata_fields=metadata_fields,
            obj_name_idx=obj_name_idx,
        )

        return inventory_report_config

    async def _fetch_inventory_report_metadata(
        gcs_file_system, inventory_report_config
    ):
        """
        Fetches all inventory report metadata from GCS based on the specified
        inventory report config.

        Parameters:
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            inventory_report_config (InventoryReportConfig): An instance of
            the InventoryReportConfig class representing the inventory report
            configuration.

        Returns:
            list: A list containing dictionaries representing the metadata of
            objects from the inventory reports.

        Raises:
            ValueError: If the fetched inventory reports are empty.
        """
        # There might be multiple inventory reports in the bucket.
        inventory_report_metadata = []

        # Extract out bucket and destination path of the inventory reports.
        bucket = inventory_report_config.bucket
        destination_path = inventory_report_config.destination_path

        # Fetch the first page.
        page = await gcs_file_system._call(
            "GET", "b/{}/o", bucket, prefix=destination_path, json_out=True
        )

        inventory_report_metadata.extend(page.get("items", []))
        next_page_token = page.get("nextPageToken", None)

        # Keep fetching new pages as long as next page token exists.
        # Note that the iteration in the while loop should most likely
        # be minimal. For reference, a million objects is split up into
        # two reports, and if the report is generated daily, then in a year,
        # there will be roughly ~700 reports generated, which will still be
        # fetched in a single page.
        while next_page_token is not None:
            page = await gcs_file_system._call(
                "GET",
                "b/{}/o",
                bucket,
                prefix=destination_path,
                json_out=True,
                pageToken=next_page_token,
            )

            inventory_report_metadata.extend(page.get("items", []))
            next_page_token = page.get("nextPageToken", None)

        # If no reports are fetched, indicates there is an error.
        if len(inventory_report_metadata) == 0:
            raise ValueError(
                "No inventory reports to fetch. Check if \
                your inventory report is set up correctly."
            )

        return inventory_report_metadata

    def _sort_inventory_report_metadata(unsorted_inventory_report_metadata):
        """
        Sorts the inventory report metadata based on the 'timeCreated' field
        in reverse chronological order.

        Parameters:
            unsorted_inventory_report_metadata (list): A list of dictionaries
            representing the metadata of objects from the inventory reports.

        Returns:
            list: A sorted list of dictionaries representing the inventory
            report metadata, sorted in reverse chronological order based
            on 'timeCreated'.
        """
        return sorted(
            unsorted_inventory_report_metadata,
            key=lambda ir: InventoryReport._convert_str_to_datetime(
                ir.get("timeCreated")
            ),
            reverse=True,
        )

    async def _download_inventory_report_content(
        gcs_file_system, inventory_report_metadata, bucket
    ):
        """
        Downloads the most recent inventory report content from GCS based on
        the inventory report metadata.

        Parameters:
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            inventory_report_metadata (list): A list of dictionaries
            representing the metadata of objects from the inventory reports.
            bucket (str): The name of the GCS bucket containing
            the inventory reports.

        Returns:
            list: A list containing the content of the most recent inventory
            report as strings.
        """
        # Get the most recent inventory report date.
        most_recent_inventory_report = inventory_report_metadata[0]
        most_recent_date = InventoryReport._convert_str_to_datetime(
            most_recent_inventory_report.get("timeCreated")
        ).date()

        inventory_report_content = []

        # Run a for loop here, since there might be multiple inventory reports
        # generated on the same day. For reference, 1 million objects will be
        # split into only 2 inventory reports, so it is very rare that there
        # will be many inventory reports on the same day. But including this
        # logic for robustness.
        for metadata in inventory_report_metadata:
            inventory_report_date = InventoryReport._convert_str_to_datetime(
                metadata["timeCreated"]
            ).date()

            if inventory_report_date == most_recent_date:
                # Download the raw inventory report if the date matches.
                # Header is not needed, we only need to process and store
                # the content.
                _header, encoded_content = await gcs_file_system._call(
                    "GET", "b/{}/o/{}", bucket, metadata.get("name"), alt="media"
                )

                # Decode the binary content into string for the content.
                decoded_content = encoded_content.decode()

                inventory_report_content.append(decoded_content)

        return inventory_report_content

    def _parse_inventory_report_content(
        gcs_file_system,
        inventory_report_content,
        inventory_report_config,
        use_snapshot_listing,
        bucket,
    ):
        """
        Parses the raw inventory report content and extracts object details.

        Parameters:
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            inventory_report_content (list): A list of strings containing the
            raw content of the inventory report.
            inventory_report_config (InventoryReportConfig): An instance of the
            InventoryReportConfig class representing the inventory report
            configuration.
            use_snapshot_listing (bool): A flag indicating whether to use snapshot
            listing in the inventory report.
            bucket (str): The name of the GCS bucket containing the inventory
            reports.

        Returns:
            list: A list of dictionaries representing object details parsed
            from the inventory report content.
        """
        # Get the csv configuration for each inventory report.
        csv_options = inventory_report_config.csv_options
        record_separator = csv_options.get("recordSeparator", "\n")
        delimiter = csv_options.get("delimiter", ",")
        header_required = csv_options.get("headerRequired", False)

        objects = []

        for content in inventory_report_content:
            # Split the content into lines based on the specified separator.
            lines = content.split(record_separator)

            # Remove the header, if present.
            if header_required:
                lines = lines[1:]

            # Parse each line of the inventory report.
            for line in lines:
                obj = InventoryReport._parse_inventory_report_line(
                    inventory_report_line=line,
                    use_snapshot_listing=use_snapshot_listing,
                    gcs_file_system=gcs_file_system,
                    inventory_report_config=inventory_report_config,
                    delimiter=delimiter,
                    bucket=bucket,
                )

                objects.append(obj)

        return objects

    def _parse_inventory_report_line(
        inventory_report_line,
        use_snapshot_listing,
        gcs_file_system,
        inventory_report_config,
        delimiter,
        bucket,
    ):
        """
        Parses a single line of the inventory report and extracts object details.

        Parameters:
            inventory_report_line (str): A string representing a single line of
            the raw content from the inventory report.
            use_snapshot_listing (bool): A flag indicating whether to use snapshot
            listing in the inventory report.
            gcs_file_system (GCSFileSystem): An instance of the 'GCSFileSystem'
            class (see 'core.py').
            inventory_report_config (InventoryReportConfig): An instance of the
            InventoryReportConfig class representing the inventory report
            configuration.
            delimiter (str): The delimiter used in the inventory report content
            to separate fields.
            bucket (str): The name of the GCS bucket containing the inventory
            reports.

        Returns:
            dict: A dictionary representing object details parsed from the
            inventory report line.
        """
        obj_name_idx = inventory_report_config.obj_name_idx
        metadata_fields = inventory_report_config.metadata_fields

        # If the client wants to do listing from the snapshot, we need
        # to fetch all the metadata for each object. Otherwise, we only
        # need to fetch the name.
        if use_snapshot_listing is True:
            obj = gcs_file_system._process_object(
                {
                    key: value
                    for key, value in zip(
                        metadata_fields, inventory_report_line.strip().split(delimiter)
                    )
                },
                bucket,
            )
        else:
            obj = {"name": inventory_report_line.strip().split(delimiter)[obj_name_idx]}

        return obj

    def _construct_final_snapshot(objects, prefix, use_snapshot_listing):
        """
        Constructs the final snapshot based on the retrieved objects and prefix.

        Parameters:
            objects (list): A list of dictionaries representing object details
            from the inventory report.
            prefix (str): A prefix used to filter objects in the snapshot based
            on their names.
            use_snapshot_listing (bool): A flag indicating whether to use snapshot
            listing in the inventory report.

        Returns:
            tuple: A tuple containing two lists: the 'items' list representing
            object details for the snapshot, and the 'prefixes' list containing
            subdirectory prefixes. If 'use_snapshot_listing' is set to False,
            'prefix' will also be empty, and 'items' will contains the object
            names in the snapshot.
        """
        if prefix is None:
            prefix = ""

        # Filter the prefix and returns the list if the user does not want to use
        # the snapshot for listing.
        if use_snapshot_listing is False:
            return [obj for obj in objects if obj.get("name").startswith(prefix)], []

        else:
            # If the user wants to use the snapshot, generate both the items and
            # prefixes manually.
            items = []
            prefixes = set()

            for obj in objects:
                # Fetch the name of the object.
                obj_name = obj.get("name")

                # If the object name doesn't start with the prefix, continue.
                # In the case where prefix is empty, it will always return
                # true (which is the expected behavior).
                if not obj_name.startswith(prefix):
                    continue

                # Remove the prefix.
                object_name_no_prefix = obj_name[len(prefix) :]

                # Determine whether the object name is a directory.
                first_delimiter_idx = object_name_no_prefix.find("/")

                # If not, then append it to items.
                if first_delimiter_idx == -1:
                    items.append(obj)
                    continue

                # If it is, recompose the directory and add to the prefix set.
                dir = object_name_no_prefix[:first_delimiter_idx]
                obj_prefix = (
                    prefix.rstrip("/")
                    + ("" if prefix == "" else "/")
                    + dir
                    + ("" if dir == "" else "/")
                )
                prefixes.add(obj_prefix)

        return items, list(prefixes)

    @staticmethod
    def _convert_obj_to_date(obj):
        """
        Converts a dictionary representing a date object to a datetime object.

        Parameters:
            obj (dict): A dictionary representing a date object with keys "day",
            "month", and "year".

        Returns:
            datetime: A datetime object representing the converted date.
        """
        day = obj["day"]
        month = obj["month"]
        year = obj["year"]
        return datetime(year, month, day)

    @staticmethod
    def _convert_str_to_datetime(str):
        """
        Converts an ISO-formatted date string to a datetime object.

        Parameters:
            date_string (str): An ISO-formatted date string with or without
            timezone information (Z).

        Returns:
            datetime: A datetime object representing the converted date and time.
        """
        return datetime.fromisoformat(str.replace("Z", "+00:00"))


class InventoryReportConfig(object):
    """
    Represents the configuration for fetching inventory reports.

    Attributes:
        csv_options (dict): A dictionary containing options for parsing CSV
        format in the inventory reports.
        bucket (str): The name of the GCS bucket from which to fetch the
        inventory reports.
        destination_path (str): The path within the GCS bucket where the
        inventory reports are stored.
        metadata_fields (list): A list of strings representing metadata
        fields to be extracted from the inventory reports.
        obj_name_idx (int): The index of the "name" field in the 'metadata_fields'
        list, used to identify object names.
    """

    def __init__(
        self, csv_options, bucket, destination_path, metadata_fields, obj_name_idx
    ):
        self.csv_options = csv_options
        self.bucket = bucket
        self.destination_path = destination_path
        self.metadata_fields = metadata_fields
        self.obj_name_idx = obj_name_idx
