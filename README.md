Superbio.ai

Package that can validate and return preview for different files.

Currenty supported extensions:

- .h5ad
- .csv
- .tsv
- .txt
- .pdb
- .fasta

Usage:

```
pip install file_process
```

```
from file_process import preview_file
try:
    processor = FileProcessFactory.get(file_data.filename, file_data.file)
except WrongExtension as exc:
    raise CustomException from exc
processor.validate(model_metadata_file)
var_target_names, obs_target_names, obs_preview, var_preview, text_preview  = processor.get_preview()
```

where:

- file: an object of io.BytesIo or FileStorage which will be validated and previewed
- filename: name of the validated file (only it's extention will be used)
- model_metadata_file (optional parameter): file with metadata of a model that will be used for validation. If this file is provided, the code will check that this file has the same set of columns as the validated file. If this file is not provided, no validation will be applied.

The code returns a list of var targets and obs targets (columns from the validated file), var preview and obs preview and text preview.
If some is not applicable for the file, None will be returned).

How to release:

1. Make your changes
2. Change the version in setup.py
3. Merge changes into the main branch; it will be automatically released to [pypi](https://pypi.org/project/file-process/#history)
