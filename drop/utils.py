from pathlib import Path
from snakemake.logging import logger
import wbuild


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
        filename = Path(dict_[key])
        if not filename.exists():
            raise FileExistsError(filename)
        # get real path
        dict_[key] = str(filename.expanduser())
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


def subsetBy(df, column, values, exact_match=True):
    """
    Subset by one or more values of different columns from data frame
    :param df: data frame
    :param column: column to subset by
    :param values: values to subset by
    :param exact_match: default True. when False match substrings. Important for subsetting drop groups
    :return: df subset by values and column
    """
    if values is None:
        return df
#    elif isinstance(values, str):
#        return df[df[column].str.contains("(^|,)" + values + "(,|$)")]
#    else:
#        if(values.__len__() > 1):
#            raise ValueError(f"Values too long! Please report this. Values are: {values}.")
    elif exact_match:
        if isinstance(values, str):
            return df[df[column] == values]
        else:
            return df[df[column].isin(values)]
    elif isinstance(values,str):
        return df[df[column].str.contains(values)]
    else:
        # this does not work. If you have a group drop and drop_2 
        # first group will match the second. Also if you have drop,drop_2 
        # it will not detect on or the other. 
        return df[df[column].str.contains("|".join(values))]
