import pytest

class Test_SampleAnnotation:

    @pytest.fixture(scope="class")
    def sampleAnnotation(self, dropConfig):
        return dropConfig.sampleAnnotation

    def test_columns(self, sampleAnnotation):
        parsed_cols = set(list(sampleAnnotation.annotationTable))
        def_cols = set(sampleAnnotation.SAMPLE_ANNOTATION_COLUMNS)
        assert def_cols <= parsed_cols

    def test_mapping(self, sampleAnnotation):
        # ID mappings/groups
        assert sampleAnnotation.idMapping.shape == (14, 2)
        assert sampleAnnotation.sampleFileMapping.shape == (24, 4)
        true_mapping = {'batch_0': 1, 'batch_1': 2, 'mae': 3, 'outrider_external': 8, 'outrider': 10, 
                        'fraser': 10, 'fraser_external': 8}
        assert true_mapping == {k: len(v) for k, v in sampleAnnotation.rnaIDs.items()}
        assert true_mapping == {k: len(v) for k, v in sampleAnnotation.dnaIDs.items()}

    @pytest.mark.parametrize(
        "sample_id,file_type,file_name",
        [
            ("HG00096", "RNA_BAM_FILE", "Data/rna_bam/HG00096_ncbi.bam"),
            ("HG00103", "RNA_BAM_FILE", "Data/rna_bam/HG00103.bam"),
            ("HG00096", "DNA_VCF_FILE", "Data/dna_vcf/demo_chr21_ncbi.vcf.gz"),
            ("HG00103", "DNA_VCF_FILE", "Data/dna_vcf/demo_chr21.vcf.gz"),
            ("HG00178", "GENE_COUNTS_FILE", "Data/external_count_data/geneCounts.tsv.gz"),
            ("HG00201", "SPLICE_COUNTS_DIR", "Data/external_count_data")
        ]
    )
    def test_filePaths(self, demo_dir, sampleAnnotation, sample_id, file_type, file_name):
        true_path = f"{demo_dir}/{file_name}"
        test_path = sampleAnnotation.getFilePath(sample_id, file_type)
        assert true_path == test_path

    @pytest.mark.parametrize(
        "annotation,group,files",
        [
            ("v29", "outrider_external", {'Data/external_count_data/geneCounts.tsv.gz'})
        ]
    )
    def test_import(self, demo_dir, sampleAnnotation, annotation, group, files):
        true_import_files = {f"{demo_dir}/{file}" for file in files}
        test_import_files = sampleAnnotation.getImportCountFiles(annotation, group)
        assert true_import_files == test_import_files

