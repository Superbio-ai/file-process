import pytest
from numpy import nan

from file_process.exceptions import ModelFileValidationError, NoColumnsError
from file_process.h5 import H5ADFileProcessor
from tests.test_file_process import H5AD_INPUT_FILES_PATH, get_remote_file_obj


class TestH5ADFileProcessor:
    path = f'{H5AD_INPUT_FILES_PATH}/heart_sample.h5ad'
    valid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_model.csv'
    invalid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_invalid.csv'

    def test_read_file(self):
        file_bytes_io = get_remote_file_obj(self.path)
        res = H5ADFileProcessor(file_bytes_io)
        assert res

    def test_get_preview(self):
        file_bytes_io = get_remote_file_obj(self.path)
        target_names, obs_preview, var_preview = H5ADFileProcessor(file_bytes_io).get_preview()
        assert target_names == ['NRP', 'age_group', 'cell_source', 'cell_type', 'donor', 'gender', 'n_counts',
                                'n_genes', 'percent_mito', 'percent_ribo', 'region', 'sample',  'scrublet_score',
                                'source', 'type', 'version', 'cell_states', 'Used', '_scvi_batch', '_scvi_labels']
        assert obs_preview == [
            {'NRP': 'Yes', 'age_group': '65-70', 'cell_source': 'Sanger-CD45', 'cell_type': 'Myeloid', 'donor': 'D6', 'gender': 'Male', 'n_counts': 1420.0, 'n_genes': 738, 'percent_mito': 0.05492957681417465, 'percent_ribo': 0.06478872895240784, 'region': 'LA', 'sample': 'HCAHeart7844001', 'scrublet_score': 0.11347517730496454, 'source': 'CD45+', 'type': 'DCD', 'version': 'V2', 'cell_states': 'LYVE1+MØ1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'AACTCCCCACGAGAGT-1-HCAHeart7844001'},
            {'NRP': 'No', 'age_group': '70-75', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'D4', 'gender': 'Female', 'n_counts': 844.0, 'n_genes': 505, 'percent_mito': 0.0011848341673612595, 'percent_ribo': 0.0011848341673612595, 'region': 'RV', 'sample': 'HCAHeart7829979', 'scrublet_score': 0.0855457227138643, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'vCM1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'ATAACGCAGAGCTGGT-1-HCAHeart7829979'},
            {'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Fibroblast', 'donor': 'D2', 'gender': 'Male', 'n_counts': 1491.0, 'n_genes': 862, 'percent_mito': 0.0, 'percent_ribo': 0.005365526303648949, 'region': 'RA', 'sample': 'HCAHeart7702879', 'scrublet_score': 0.19786096256684493, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'FB2', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'GTCAAGTCATGCCACG-1-HCAHeart7702879'},
            {'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-CD45', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 2167.0, 'n_genes': 1115, 'percent_mito': 0.06414397805929184, 'percent_ribo': 0.027226580306887627, 'region': 'LA', 'sample': 'HCAHeart8102858', 'scrublet_score': 0.11347517730496454, 'source': 'CD45+', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC10_CMC-like', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'GGTGATTCAAATGAGT-1-HCAHeart8102858'},
            {'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Cells', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 7334.0, 'n_genes': 2505, 'percent_mito': 0.09353695064783096, 'percent_ribo': 0.04049631953239441, 'region': 'RA', 'sample': 'HCAHeart8102863', 'scrublet_score': 0.13214990138067062, 'source': 'Cells', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC5_art', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'AGAGAATTCTTAGCAG-1-HCAHeart8102863'},
            {'NRP': 'No', 'age_group': '55-60', 'cell_source': 'Harvard-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'H4', 'gender': 'Male', 'n_counts': 3397.0, 'n_genes': 1597, 'percent_mito': 0.0014718869933858514, 'percent_ribo': 0.0005887547740712762, 'region': 'LV', 'sample': 'H0037_LV', 'scrublet_score': 0.10806174957118353, 'source': 'Nuclei', 'type': 'DBD', 'version': 'V3', 'cell_states': 'vCM2', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'AGGTTACCAGATAAAC-1-H0037_LV'},
            {'NRP': 'No', 'age_group': '65-70', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Myeloid', 'donor': 'D5', 'gender': 'Female', 'n_counts': 1273.0, 'n_genes': 797, 'percent_mito': 0.00235663796775043, 'percent_ribo': 0.003142183879390359, 'region': 'LV', 'sample': 'HCAHeart7835148', 'scrublet_score': 0.046610169491525424, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'LYVE1+MØ1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'AGCCTAATCTACCAGA-1-HCAHeart7835148'},
            {'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Endothelial', 'donor': 'D11', 'gender': 'Female', 'n_counts': 1347.0, 'n_genes': 830, 'percent_mito': 0.0029695620760321617, 'percent_ribo': 0.0029695620760321617, 'region': 'AX', 'sample': 'HCAHeart8287128', 'scrublet_score': 0.09375, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V3', 'cell_states': 'EC5_art', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'AGCGCTGAGGCTTAGG-1-HCAHeart8287128'},
            {'NRP': 'Yes', 'age_group': '60-65', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Fibroblast', 'donor': 'D11', 'gender': 'Female', 'n_counts': 1785.0, 'n_genes': 1098, 'percent_mito': 0.0005602241144515574, 'percent_ribo': 0.0033613445702940226, 'region': 'AX', 'sample': 'HCAHeart8287128', 'scrublet_score': 0.24290220820189273, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V3', 'cell_states': 'FB3', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'GAGTCATTCTCCGTGT-1-HCAHeart8287128'},
            {'NRP': 'No', 'age_group': '70-75', 'cell_source': 'Sanger-Nuclei', 'cell_type': 'Ventricular_Cardiomyocyte', 'donor': 'D4', 'gender': 'Female', 'n_counts': 4269.0, 'n_genes': 1955, 'percent_mito': 0.0021082221064716578, 'percent_ribo': 0.0018739751540124416, 'region': 'LV', 'sample': 'HCAHeart7829978', 'scrublet_score': 0.16473317865429235, 'source': 'Nuclei', 'type': 'DCD', 'version': 'V2', 'cell_states': 'vCM1', 'Used': 'Yes', '_scvi_batch': 0, '_scvi_labels': 0, 'Feature Name': 'TGGGAAGAGTCGAGTG-1-HCAHeart7829978'}
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

    def test_get_preview_file_with_nans(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/follicular_sample.h5ad')
        _, obs_preview, _ = H5ADFileProcessor(file_bytes_io).get_preview()
        for item in obs_preview:
            for value in item.values():
                assert value is not nan

    def test_get_preview_file_no_rows(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/liver_sample.h5ad')
        target_names, obs_preview, var_preview = H5ADFileProcessor(file_bytes_io).get_preview()
        assert obs_preview != []
        assert var_preview == []

    def test_model_file_validation_with_h5(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_file_bytes_io = get_remote_file_obj(self.valid_model_path)
        _ = H5ADFileProcessor(file_bytes_io).model_file_validation(model_file_bytes_io)

    def test_model_file_validation_with_h5_invalid_model(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_file_bytes_io = get_remote_file_obj(self.invalid_model_path)
        with pytest.raises(ModelFileValidationError):
            _ = H5ADFileProcessor(file_bytes_io).model_file_validation(model_file_bytes_io)

    def test_validate_no_columns(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/pbmc3k_raw.h5ad')
        with pytest.raises(NoColumnsError):
            H5ADFileProcessor(file_bytes_io).validate()
