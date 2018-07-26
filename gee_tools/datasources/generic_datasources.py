"""
Author: Anthony Perez

A collection of ImageDatasource classes which allow loading of generic image collections by name.
"""
import ee
from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, DatasourceError

class GenericSingleImageDatasource(SingleImageDatasource):
    """Generic SingleImageDatasource for loading an image by name."""

    def build_img_coll(self, name=None):
        if name is None:
            raise ValueError('name must be provided, but was None')
        self.img = ee.Image(name)
        self.coll = ee.ImageCollection([self.img])

    def get_img_coll(self):
        return self.coll


class GenericGlobalImageDatasource(GlobalImageDatasource):
    """Generic GlobalImageDatasource for loading an image by name."""

    def build_img_coll(self, name=None):
        if name is None:
            raise ValueError('name must be provided, but was None')
        self.coll = ee.ImageCollection(name) \
                      .filterDate(self.start_date, self.end_date) \
                      .sort('system:time_start')

    def get_img_coll(self):
        return self.coll


class GenericMultiImageDatasource(MultiImageDatasource):
    """Generic GlobalImageDatasource for loading an image by name."""

    def build_img_coll(self, name=None):
        if name is None:
            raise ValueError('name must be provided, but was None')
        self.coll = ee.ImageCollection(name) \
                      .filterDate(self.start_date, self.end_date) \
                      .filterBounds(self.filterpoly) \
                      .sort('system:time_start')

    def get_img_coll(self):
        return self.coll
