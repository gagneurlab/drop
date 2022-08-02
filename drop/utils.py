from pathlib import Path
from snakemake.logging import logger
import wbuild
import copy


def returnPath(path, str_=True):
    return str(path) if str_ else Path(path)


def checkFileExists(files):
    """
    :param files: single filename or iterable of filenames
    :return: list of existing files
    """
    files = [files] if isinstance(files, str) else files
    return [f for f in files if Path(f).exists()]


def createDir(directory):
    directory = Path(directory)
    if not directory.exists():
        logger.debug(f"creating {directory}")
        directory.mkdir(parents=True)


def checkKeys(dict_, keys=None, check_files=False):
    """
    :param dict_: config dictionary
    :param keys: keys that are expected to be in dict_
    :param file_values: if True, do file checks on values
    :return: dictionary with absolute file paths
    """
    keys = dict_.keys() if keys is None else keys
    for key in keys:
        if key not in dict_.keys():
            raise KeyError(f"{key} is mandatory but missing")

        if check_files:
            existing = checkFileExists(dict_[key])
            # get real path
            if type(dict_[key]) == str:
                if len(existing) != 1:
                    raise FileExistsError(dict_[key])
                dict_[key] = str(Path(dict_[key]).expanduser())
            elif type(dict_[key]) == list:
                if len(existing) < len(dict_[key]):
                    raise FileExistsError(set(dict_[key]) - set(existing))
                dict_[key] = [str(Path(f).expanduser()) for f in dict_[key]]
    return dict_


def setKey(dict_, sub, key, default):
    if sub is not None:
        if not isinstance(sub, list):
            raise TypeError(f"{sub} is not of type list")
        for x in sub:
            dict_ = dict_[x]
    if key not in dict_ or dict_[key] is None:
        logger.debug(f'{key} not in config{sub}, using default')
        dict_[key] = default
    return dict_[key]


def getWBuildPath(str_=True):
    return returnPath(Path(wbuild.__file__).parent, str_=str_)


def getWBuildSnakefile(str_=True):
    wb_path = getWBuildPath(str_=False)
    return returnPath(wb_path / "wBuild.snakefile", str_=str_)


def subsetBy(df, column, values):
    """
    Subset by one or more values of different columns from data frame
    :param df: data frame
    :param column: column to subset by
    :param values: values to subset by
    :return: df subset by values and column
    """
    if values is None:
        return df
    
    inner_regex = values
    if not isinstance(values, str) :
        inner_regex = "(" + "|".join(values) + ")"
    
    if df[column].isnull().all():
        return df[[False for i in range(df.shape[0])]]
    return  df[df[column].str.contains("(?:^|,)" + inner_regex + "(?:,|$)", na = False)]
    
def deep_merge_dict(dict1: dict, dict2: dict, inplace: bool = False):
    """
    Merges two dictionaries and all is children recursively
    
    :param dict1: dictionary to be merged into
    :param dict2: dictionary to be merged
    :param inplace: if False, default, a new dictionary will be returned als in-place merging is performed.
    """
    if not inplace:
        dict1 = copy.deepcopy(dict1)
        dict2 = copy.deepcopy(dict2)
    
    for k, v in dict2.items():
        if isinstance(dict1.get(k), dict) and isinstance(v, dict):
            dict1[k] = deep_merge_dict(dict1[k], v, inplace=inplace)
        elif k not in dict1:
            dict1[k] = v
        elif isinstance(dict1.get(k), list) and isinstance(v, list):
            dict1[k] = list(dict.fromkeys(dict1[k] + v))
        elif isinstance(dict1.get(k), str) and isinstance(v, str):
            dict1[k] = [dict1.get(k), v]
        else:
            raise TypeError(f"{k} has different types that can not be merged.")
        
    return dict1
