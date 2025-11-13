Loading Data into Pandas in Google Colab
Mount Google Drive to Access Files:
This connects your Google Drive storage to the Colab environment so files can be read as if local.

python
from google.colab import drive
drive.mount('/content/drive')
After executing, follow the prompt to authorize access. Your Drive files will appear under /content/drive/.

Locate Your File Path in Drive:
Use the file explorer pane in Colab to navigate to your file, right-click it, and select "Copy path." The path will be something like /content/drive/MyDrive/path_to_file.csv.

Load the File into a Pandas DataFrame:
Using pandas, read the data with:

python
import pandas as pd
df = pd.read_csv('/content/drive/MyDrive/path_to_file.csv')
df.head()
This works similarly for Excel files (pd.read_excel) or JSON files (pd.read_json).

Alternative: Upload Files Directly from Local to Colab:
If you don't want to use Drive, upload files interactively:

python
from google.colab import files
uploaded = files.upload()
import io
df = pd.read_csv(io.BytesIO(uploaded['filename.csv']))
Replace 'filename.csv' with the uploaded file’s name.

Suggested Addition to Pandas Documentation for Google Colab
A dedicated section on the IO page explaining this process would be useful for users running pandas in Colab, covering Drive mounting, direct uploading, and basic file reading commands.
