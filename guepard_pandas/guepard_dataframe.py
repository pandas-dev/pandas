import pandas as pd
import requests

class GuepardDataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.api_url = "https://api.guepard.com"
        self.dataset_id = kwargs.get('dataset_id', 'default')
    
    def commit(self, message=""):
        version_id = self._generate_version_id()
        data = self.to_parquet()
        response = requests.post(f"{self.api_url}/datasets/{self.dataset_id}/versions",
                                 files={"data": data},
                                 data={"message": message, "version_id": version_id})
        response.raise_for_status()
        return version_id
    
    def list_versions(self):
        response = requests.get(f"{self.api_url}/datasets/{self.dataset_id}/versions")
        response.raise_for_status()
        return response.json()
    
    def rollback(self, version_id):
        response = requests.get(f"{self.api_url}/datasets/{self.dataset_id}/versions/{version_id}")
        response.raise_for_status()
        data = response.content
        df = pd.read_parquet(data)
        self.__init__(df)
    
    def next_version(self):
        return self.commit()

    def _generate_version_id(self):
        from datetime import datetime
        return datetime.now().strftime("%Y%m%d_%H%M%S")