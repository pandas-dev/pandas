from zipfile import ZipFile

from openpyxl.chart._chart import ChartBase
from openpyxl.drawing.image import Image

def find_images(archive: ZipFile, path: str) -> tuple[list[ChartBase], list[Image]]: ...
