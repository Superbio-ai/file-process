Superbio.ai

Package that can validate and return preview for different files.

Currenty supported extensions:

- .h5, .h5ad
- .csv

Usage:

```
pip install file_process
```

```
from file_process import preview_file
target_names, var_preview, obs_preview = preview file(filename, file, model_metadata_file)
```

where:

- file: an object of io.BytesIo or FileStorage which will be validated and previewed
- filename: name of the validated file (only it's extention will be used)
- model_metadata_file (optional parameter): file with metadata of a model that will be used for validation. If this file is provided, the code will check that this file has the same set of columns as the validated file. If this file is not provided, no validation will be applied.

The code returns a list of targets (columns from the validated file), var preview and obs preview.
If some is not applicable for the file, None will be returned).
