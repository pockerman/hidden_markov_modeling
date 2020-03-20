"""
Preprocessing utilities
"""

from pomegranate import*

from exceptions import Error

VALID_DISTS = ['gaussian', 'uniform', 'poisson']


def fit_distribution(data, dist_name="gaussian", **kwargs):

    """
    Fits a distribution within the given dataset
    :param data:
    :param dist_name:
    :param kwargs:
    :return: appropriate distribution object
    """

    if dist_name not in VALID_DISTS:
        raise Error("Invalid distribution name. Name {0} not in {1}".format(dist_name, VALID_DISTS))

    if dist_name == 'gaussian':
        dist = NormalDistribution.from_samples(data)
        return dist
    elif dist_name == 'uniform':
        dist = UniformDistribution.from_samples(data)
        return dist
    elif dist_name == 'poisson':
        dist = PoissonDistribution.from_samples(data)
        return dist


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