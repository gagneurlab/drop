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
                 "v29": ["geneAnnotation"],
                 "genome": ["mae"], "qcVcf": ["mae"]}

    for key, sub in path_keys.items():
        # iterate to key and entry
        dict_ = config
        if sub is not None:
            for x in sub:
                dict_ = dict_[x]
        # set absolute path
        dict_[key] = str(Path(dict_[key]).resolve())

    with open(output_file, "w") as f:
        yaml.safe_dump(config.copy(), f, default_flow_style=False, sort_keys=False)