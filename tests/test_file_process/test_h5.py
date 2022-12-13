import pytest
from anndata import AnnData

from werkzeug.datastructures import FileStorage

from file_process.exceptions import ModelFileValidationError, NoColumnsError
from file_process.h5 import H5ADFileProcessor
from tests.test_file_process import H5AD_INPUT_FILES_PATH, get_remote_file_obj


class TestH5ADFileProcessor:
    path = f'{H5AD_INPUT_FILES_PATH}/heart_sample.h5ad'
    valid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_model.csv'
    invalid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_invalid.csv'

    def test_read_remote_file(self):
        file_obj = get_remote_file_obj(self.path)
        res = H5ADFileProcessor(file_obj)
        assert res

    def test_read_local_file(self):
        with open(self.path, 'rb') as file:
            file_obj = FileStorage(file)
            res = H5ADFileProcessor(file_obj)
            assert isinstance(res.data, AnnData)

    def test_get_preview(self):
        file_obj = get_remote_file_obj(self.path)
        file_processor = H5ADFileProcessor(file_obj)
        target_names, obs_preview, var_preview = file_processor.get_preview()
        assert target_names == ['NRP', 'age_group', 'cell_source', 'cell_type', 'donor', 'gender', 'n_counts',
                                'n_genes', 'percent_mito', 'percent_ribo', 'region', 'sample',  'scrublet_score',
                                'source', 'type', 'version', 'cell_states', 'Used', '_scvi_batch', '_scvi_labels']
        assert obs_preview == [
            {'Feature Name': 'AACTCCCCACGAGAGT-1-HCAHeart7844001', 'NRP': 'Yes', 'age_group': '65-70', 'cell_source': 'Sanger-CD45', 'cell_type': 'Myeloid', 'donor': 'D6', 'gender': 'Male', 'n_counts': 1420.0, 'n_genes': 738, 'percent_mito': 0.05000000074505806, 'percent_ribo': 0.05999999865889549, 'region': 'LA', 'sample': 'HCAHeart7844001', 'scrublet_score': 0.11, 'source': 'CD45+', 'type': 'DCD', 'version': 'V2', 'cell_states': 'LYVE1+MØ1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'ATAACGCAGAGCTGGT-1-HCAHeart7829979', 'NRP': 'No', 'age_group': '70-75', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'D4', 'gender': 'Female', 'n_counts': 844.0, 'n_genes': 505, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'RV', 'sample': 'HCAHeart7829979', 'scrublet_score': 0.09, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'vCM1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'GTCAAGTCATGCCACG-1-HCAHeart7702879', 'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Fibroblast', 'donor': 'D2', 'gender': 'Male', 'n_counts': 1491.0, 'n_genes': 862, 'percent_mito': 0.0, 'percent_ribo': 0.009999999776482582, 'region': 'RA', 'sample': 'HCAHeart7702879', 'scrublet_score': 0.2, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'FB2', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'GGTGATTCAAATGAGT-1-HCAHeart8102858', 'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-CD45', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 2167.0, 'n_genes': 1115, 'percent_mito': 0.05999999865889549, 'percent_ribo': 0.029999999329447746, 'region': 'LA', 'sample': 'HCAHeart8102858', 'scrublet_score': 0.11, 'source': 'CD45+', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC10_CMC-like', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'AGAGAATTCTTAGCAG-1-HCAHeart8102863', 'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Cells', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 7334.0, 'n_genes': 2505, 'percent_mito': 0.09000000357627869, 'percent_ribo': 0.03999999910593033, 'region': 'RA', 'sample': 'HCAHeart8102863', 'scrublet_score': 0.13, 'source': 'Cells', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC5_art', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'AGGTTACCAGATAAAC-1-H0037_LV', 'NRP': 'No', 'age_group': '55-60', 'cell_source': 'Harvard-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'H4', 'gender': 'Male', 'n_counts': 3397.0, 'n_genes': 1597, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'LV', 'sample': 'H0037_LV', 'scrublet_score': 0.11, 'source': 'Nuclei', 'type': 'DBD', 'version': 'V3', 'cell_states': 'vCM2', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'AGCCTAATCTACCAGA-1-HCAHeart7835148', 'NRP': 'No', 'age_group': '65-70', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Myeloid', 'donor': 'D5', 'gender': 'Female', 'n_counts': 1273.0, 'n_genes': 797, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'LV', 'sample': 'HCAHeart7835148', 'scrublet_score': 0.05, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'LYVE1+MØ1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'AGCGCTGAGGCTTAGG-1-HCAHeart8287128', 'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 1347.0, 'n_genes': 830, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'AX', 'sample': 'HCAHeart8287128', 'scrublet_score': 0.09, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC5_art', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'GAGTCATTCTCCGTGT-1-HCAHeart8287128', 'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Fibroblast', 'donor': 'D11', 'gender': 'Female', 'n_counts': 1785.0, 'n_genes': 1098, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'AX', 'sample': 'HCAHeart8287128', 'scrublet_score': 0.24, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V3', 'cell_states': 'FB3', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0},
            {'Feature Name': 'TGGGAAGAGTCGAGTG-1-HCAHeart7829978', 'NRP': 'No', 'age_group': '70-75', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'D4', 'gender': 'Female', 'n_counts': 4269.0, 'n_genes': 1955, 'percent_mito': 0.0, 'percent_ribo': 0.0, 'region': 'LV', 'sample': 'HCAHeart7829978', 'scrublet_score': 0.16, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'vCM1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0}
        ]
        assert var_preview == [
            {'Feature Name': 'AL627309.1', 'gene_ids-Harvard-Nuclei': 'ENSG00000238009', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000238009', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000238009', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000238009', 'feature_types-Sanger-CD45': 0, 'n_counts': 249.0},
            {'Feature Name': 'AC114498.1', 'gene_ids-Harvard-Nuclei': 'ENSG00000235146', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000235146', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000235146', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000235146', 'feature_types-Sanger-CD45': 0, 'n_counts': 28.0},
            {'Feature Name': 'AL669831.2', 'gene_ids-Harvard-Nuclei': 'ENSG00000229905', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000229905', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000229905', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000229905', 'feature_types-Sanger-CD45': 0, 'n_counts': 3.0},
            {'Feature Name': 'AL669831.5', 'gene_ids-Harvard-Nuclei': 'ENSG00000237491', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000237491', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000237491', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000237491', 'feature_types-Sanger-CD45': 0, 'n_counts': 1342.0},
            {'Feature Name': 'FAM87B', 'gene_ids-Harvard-Nuclei': 'ENSG00000177757', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000177757', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000177757', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000177757', 'feature_types-Sanger-CD45': 0, 'n_counts': 15.0},
            {'Feature Name': 'LINC00115', 'gene_ids-Harvard-Nuclei': 'ENSG00000225880', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000225880', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000225880', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000225880', 'feature_types-Sanger-CD45': 0, 'n_counts': 330.0},
            {'Feature Name': 'FAM41C', 'gene_ids-Harvard-Nuclei': 'ENSG00000230368', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000230368', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000230368', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000230368', 'feature_types-Sanger-CD45': 0, 'n_counts': 47.0},
            {'Feature Name': 'AL645608.7', 'gene_ids-Harvard-Nuclei': 'ENSG00000272438', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000272438', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000272438', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000272438', 'feature_types-Sanger-CD45': 0, 'n_counts': 11.0},
            {'Feature Name': 'AL645608.1', 'gene_ids-Harvard-Nuclei': 'ENSG00000223764', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000223764', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000223764', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000223764', 'feature_types-Sanger-CD45': 0, 'n_counts': 43.0},
            {'Feature Name': 'SAMD11', 'gene_ids-Harvard-Nuclei': 'ENSG00000187634', 'feature_types-Harvard-Nuclei': 'Gene Expression', 'gene_ids-Sanger-Nuclei': 'ENSG00000187634', 'feature_types-Sanger-Nuclei': 0, 'gene_ids-Sanger-Cells': 'ENSG00000187634', 'feature_types-Sanger-Cells': 0, 'gene_ids-Sanger-CD45': 'ENSG00000187634', 'feature_types-Sanger-CD45': 0, 'n_counts': 171.0}
        ]

    def test_model_file_validation_with_h5(self):
        model_file_obj = get_remote_file_obj(self.valid_model_path)
        test_file_obj = get_remote_file_obj(self.path)
        file_processor = H5ADFileProcessor(test_file_obj)
        file_processor.model_file_validation(model_file_obj)

    def test_model_file_validation_with_h5_invalid_model(self):
        model_file_obj = get_remote_file_obj(self.invalid_model_path)
        test_file_obj = get_remote_file_obj(self.path)
        file_processor = H5ADFileProcessor(test_file_obj)
        with pytest.raises(ModelFileValidationError):
            file_processor.model_file_validation(model_file_obj)

    def test_validate_no_columns(self):
        test_file_obj = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/pbmc3k_raw.h5ad')
        file_processor = H5ADFileProcessor(test_file_obj)
        with pytest.raises(NoColumnsError):
            file_processor.validate()
