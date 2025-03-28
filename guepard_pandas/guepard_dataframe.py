import pandas as pd
import os
import pickle
from datetime import datetime

class GuepardDataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        version_dir = kwargs.pop('version_dir', './versions')
        super().__init__(*args, **kwargs)
        self.current_version_path = os.path.join(version_dir, 'current_version.pkl')
        self.version_dir = version_dir
        self.versions_meta_file = os.path.join(version_dir, 'versions_meta.pkl')
        if not os.path.exists(self.version_dir):
            os.makedirs(self.version_dir)
        self._load_current_version()
    
    def _load_current_version(self):
        if os.path.exists(self.current_version_path):
            with open(self.current_version_path, 'rb') as f:
                df = pickle.load(f)
            super().__init__(df)
    
    def commit(self, message=""):
        version_id = self._generate_version_id()
        self._save_current_version()
        self._store_version_meta(version_id, message)
        return version_id
    
    def _save_current_version(self):
        with open(self.current_version_path, 'wb') as f:
            pickle.dump(self, f)
    
    def _store_version_meta(self, version_id, message):
        versions_meta = self._load_versions_meta()
        versions_meta.append({'version_id': version_id, 'message': message, 'timestamp': datetime.now()})
        with open(self.versions_meta_file, 'wb') as f:
            pickle.dump(versions_meta, f)
    
    def _load_versions_meta(self):
        if os.path.exists(self.versions_meta_file):
            with open(self.versions_meta_file, 'rb') as f:
                return pickle.load(f)
        return []
    
    def list_versions(self):
        versions_meta = self._load_versions_meta()
        return [{'version_id': meta['version_id'], 'message': meta['message'], 'timestamp': meta['timestamp']} for meta in versions_meta]
    
    def rollback(self, version_id):
        version_path = os.path.join(self.version_dir, f"{version_id}.pkl")
        if not os.path.exists(version_path):
            raise ValueError("Version ID not found")
        with open(version_path, 'rb') as f:
            df = pickle.load(f)
        self.__init__(df)
        self._save_current_version()
    
    def save_version(self, version_id):
        version_path = os.path.join(self.version_dir, f"{version_id}.pkl")
        with open(version_path, 'wb') as f:
            pickle.dump(self, f)
    
    def _generate_version_id(self):
        return datetime.now().strftime("%Y%m%d_%H%M%S")

