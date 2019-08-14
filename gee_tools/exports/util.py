"""
Author: Anthony Perez

Utility functions for loading data.
"""
import datetime

import ee

from gee_tools.exports import constants

def check_empty_bands(scene):
    """Raise a ValueError if bands are empty (len(scene.bandNames().getInfo()) == 0)"""
    if len(scene.bandNames().getInfo()) == 0:
        raise ValueError("No bands available.")


def stack_feature_colls(feature_colls):
    """Return a feature collection with elements from all feature collections in the list feature_colls."""
    if len(feature_colls) == 0:
        raise ValueError("feature_colls was empty.")
    output = feature_colls.pop(0)
    for fc in feature_colls:
        output = output.merge(fc)
    return output


def start_date_before_end(start_date, end_date):
    """Raise a ValueError if the start date string is after the end date string. YYYY-MM-DD format."""
    s_date = parse_date(start_date)
    e_date = parse_date(end_date)
    if s_date > e_date:
        raise ValueError("Start date after end date (Start: {}, End: {}".format(s_date, e_date))


def parse_date(d_str):
    """Convert a date string of the form YYYY-MM-DD to a python date object. Single digit months and days work."""
    if isinstance(d_str, datetime.date) or isinstance(d_str, datetime.datetime):
        return d_str
    return datetime.datetime.strptime(d_str, constants.DATE_STR_FORMAT).date()


def date_to_str(date):
    """Convert a python date object to YYYY-MM-DD format."""
    if isinstance(date, str):
        return date
    if isinstance(date, datetime.date) or isinstance(date, datetime.datetime):
        return date.strptime(constants.DATE_STR_FORMAT)
    raise ValueError("Could not recognize date.\n{}".format(date))


def python_datetime_to_ee_date(py_datetime):
    """
    Convert a python datetime object into an ee.Date object.  Accurate to milliseconds.
    Args:
        py_datetime (datetime.datetime):  The date to convert.
    Returns:
        (ee.Date): The converted date.
    """
    date_str = py_datetime.isoformat(' ')  # yyyy-MM-dd HH:mm:ss
    seconds_fraction_digits = int(py_datetime.microsecond / (10.0**3)) % 10000
    date_str += '.' + str(seconds_fraction_digits)

    format_str = 'yyyy-MM-dd HH:mm:ss.SSSS'
    return ee.Date.parse(format=format_str, date=date_str)
