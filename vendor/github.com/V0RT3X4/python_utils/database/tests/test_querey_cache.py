# cd database
import logging

from vortexa_utils.database.default_factories import DevFactory
from vortexa_utils.database.query_cache import QueryCache

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.DEBUG)

# factory = DevFactory()
# engine = factory.engine()
# qc = QueryCache()

# %time df = qc.read_sql("clarksons", engine)


def test_filename():
    qc = QueryCache()
    assert qc.filename("some random query") == "qAdzxvMgeSc=.parquet.snappy"
    assert qc.filename("banned_words") == "LoRkfDuNmuA=.parquet.snappy"
