from pathlib import Path
from snakemake.logging import logger


class Genome:

    def __init__(self, annotation, assembly, reference):
        self.annotation = annotation
        self.assembly = assembly

        # Allow for old drop config stylings, where the file was a string under MAE
        # -> force into dictionary
        self.reference = {reference: reference} if isinstance(reference, str) else reference

    def getGeneAnnotations(self):
        return self.annotation

    def getGeneVersions(self):
        return self.annotation.keys()

    def getGeneAnnotationFile(self, annotation):
        """
        :param annotation: config-defined annotation key
        :return: GTF file from config key 'geneAnnotations'
        """
        return self.annotation[annotation]

    def getFastaFiles(self):
        """
        :return: dictionary of genome name -> genome file
        """
        if isinstance(self.reference, str):
            return {self.reference: self.reference}
        else:
            return self.reference

    def getFastaList(self):
        return list(self.getFastaFiles().values())

    def getFastaDict(self, fasta_file):
        return Path(fasta_file).with_suffix(".dict")
    
    def getDictList(self):
        return [self.getFastaDict(i) for i in self.getFastaList()]

    def getBSGenomeName(self):
        assemblyID = self.assembly

        if assemblyID == 'hg19':
            return "BSgenome.Hsapiens.UCSC.hg19"
        if assemblyID == 'hs37d5':
            return "BSgenome.Hsapiens.1000genomes.hs37d5"
        if assemblyID == 'hg38':
            return "BSgenome.Hsapiens.UCSC.hg38"
        if assemblyID == 'GRCh38':
            return "BSgenome.Hsapiens.NCBI.GRCh38"

        raise ValueError("Provided genome assembly not known: " + assemblyID)

    def getBSGenomeVersion(self):
        assemblyID = self.assembly

        if assemblyID in ['hg19', 'hs37d5']:
            return 37
        if assemblyID in ['hg38', 'GRCh38']:
            return 38

        raise ValueError("Provided genome assembly not known: " + assemblyID)

    def getMafDbName(self):
        assemblyID = self.assembly

        if assemblyID in ['hg19', 'hs37d5']:
            return "MafDb.gnomAD.r2.1.hs37d5"
        if assemblyID in ['hg38', 'GRCh38']:
            return "MafDb.gnomAD.r2.1.GRCh38"

        raise ValueError("Provided genome assembly not known: " + assemblyID)
