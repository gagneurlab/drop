import os
import pandas as pd
import yaml
from pathlib import Path
from drop.config.SampleAnnotation import SampleAnnotation


def fixSampleAnnotation(input_file, output_file):
    """
    :param input_file: sample annotation file with relative paths
    :param output_file: location for sample annotation file with absolute paths
    """
    sa = pd.read_csv(input_file, sep='\t')
    for key in SampleAnnotation.FILE_TYPES:
        if key not in sa.columns:
            continue
        sa[key] = [str(Path(x).resolve()) if not pd.isna(x) else x for x in sa[key] if x]
    sa.to_csv(output_file, sep='\t', index=False)


def fixConfig(input_file, output_file):
    """
    :param input_file: config file with relative paths
    :param output_file: location for config file with absolute paths
    """
    with open(input_file, "r") as f:
        config = yaml.load(f, Loader=yaml.Loader)
    path_keys = {"root": None,
                 "htmlOutputPath": None,
                 "sampleAnnotation": None,
                 "v29": "geneAnnotation",
                 "ncbi": "genome", "ucsc":"genome",
                 "qcVcf": "mae",
                 "highQualityVCFs":"rnaVariantCalling",
                 "dbSNP":"rnaVariantCalling",
                 "repeat_mask":"rnaVariantCalling"
                 }

    for key, sub in path_keys.items():
        # iterate to key and entry
        dict_ = config
        if sub is not None:
            dict_ = dict_[sub]
            # set absolute path
            if type(dict_[key]) == str:
                dict_[key] = str(Path(dict_[key]).resolve())
            elif type(dict_[key]) == list:
                dict_[key] = [str(Path(x).resolve()) for x in dict_[key]]

        else:
            dict_[key] = str(Path(dict_[key]).resolve())

    with open(output_file, "w") as f:
        yaml.safe_dump(config.copy(), f, default_flow_style=False, sort_keys=False)
