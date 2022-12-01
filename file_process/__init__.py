from typing import Optional, List

from file_process.file_processor_factory import FileProcessFactory


# pylint: disable=too-many-function-args
def preview_file(filename: str, file, *, model_metadata_file=None, **kwargs) \
        -> (List[str], Optional[List[dict]], Optional[List[dict]]):
    processor = FileProcessFactory.get(filename)
    target_names, var_preview, obs_preview = processor.process(file, model_metadata_file, **kwargs)
    return target_names, processor.create_tabular_response(var_preview), processor.create_tabular_response(obs_preview)
