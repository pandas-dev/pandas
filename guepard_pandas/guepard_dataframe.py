import pandas as pd
import os
import pickle
from datetime import datetime

class GuepardDataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.version_dir = kwargs.pop('version_dir', './versions')
        if not os.path.exists(self.version_dir):
            os.makedirs(self.version_dir)
    
    def commit(self, message=""):
        version_id = self._generate_version_id()
        version_path = os.path.join(self.version_dir, f"{version_id}.pkl")
        with open(version_path, 'wb') as f:
            pickle.dump(self, f)
        return version_id
    
    def list_versions(self):
        versions = []
        for filename in os.listdir(self.version_dir):
            if filename.endswith(".pkl"):
                version_id = filename.split('.')[0]
                versions.append(version_id)
        return versions
    
    def rollback(self, version_id):
        version_path = os.path.join(self.version_dir, f"{version_id}.pkl")
        if not os.path.exists(version_path):
            raise ValueError("Version ID not found")
        with open(version_path, 'rb') as f:
            df = pickle.load(f)
        self.__init__(df)
    
    def next_version(self):
        return self.commit()

    def _generate_version_id(self):
        return datetime.now().strftime("%Y%m%d_%H%M%S")

# Example usage:
# df = GuepardDataFrame(pd.read_csv("data.csv"), version_dir="path/to/versions")
# df["new_col"] = df["existing_col"] * 2
# df.commit("Added new column")
# print(df.list_versions())
# df.rollback(version_id="20240326_123456")