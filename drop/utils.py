from pathlib import Path
from snakemake.logging import logger
import wbuild

def returnPath(path, str_=True):
    return str(path) if str_ else Path(path)

def createIfMissing(directory):
    directory = Path(directory)
    if not directory.exists():
        logger.debug(f"creating {directory}")
        directory.mkdir(parents=True)

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

def getRuleFromPath(path, prefix=False):
    path = str(path)
    if not path.startswith("Scripts"):
        raise ValueError(f"{path} is invalid for wBuild rule")
    rule = path.replace("/", "_")
    if prefix:
        return rule.split(".")[0]
    else:
        return rule.replace(".", "_")

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
    elif isinstance(values, str):
        return df[df[column] == values]
    else:
        return df[df[column].isin(values)]
