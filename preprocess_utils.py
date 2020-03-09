"""
Preprocessing utilities
"""

def preprocess(windows, **config):
    """
    Apply preprocessing to the given windows
    by using the specified configuration
    :param windows:
    :param config:
    :return:
    """

    if "remove_outliers" in config.keys():
        remove_outliers(windows=windows, criteria=config["remove_outliers_criteria"])

    if "normalize_windows_data" in config.keys():
        remove_outliers(windows=windows, criteria=config["normalize_windows_criteria"])

    if "normality_check" in config.keys():
        for name in config["normality_check"]:
            normality_check(windows=windows, test_name=name)
    return windows

def remove_outliers(windows, criteria):
    """
    Remove the outlier windows from the list of the
    windows by using the given criteria
    :return:
    """
    pass


def normalize_windows_data(windows, criteria):
    """
    normalize the data in the given windows list
    """
    pass


def normality_check(windows, test_name):
    """
    Perform a normality test, specified by the test_name variable, for the given windows
    :param windows: The list of windows to use
    :param test_name: The name of the statistical test to use
    :return:
    """
    pass