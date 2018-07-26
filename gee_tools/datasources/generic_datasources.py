"""
Author: Anthony Perez

A collection of ImageDatasource classes which allow loading of generic image collections by name.
"""
import ee
from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, DatasourceError

class GenericSingleImageDatasource(SingleImageDatasource):
    """Generic SingleImageDatasource for loading an image by name or arguments to ee.Image."""

    def build_img_coll(self, image_args=None):
        if image_args is None:
            raise ValueError('image_args must be provided, but was None')
        self.img = ee.Image(image_args)
        self.coll = ee.ImageCollection([self.img])

    def get_img_coll(self):
        return self.coll


class GenericGlobalImageDatasource(GlobalImageDatasource):
    """Generic GlobalImageDatasource for loading an image by name or arguments to ee.ImageCollection."""

    def build_img_coll(self, coll_args=None):
        if coll_args is None:
            raise ValueError('coll_args must be provided, but was None')
        self.coll = ee.ImageCollection(coll_args) \
                      .filterDate(self.start_date, self.end_date) \
                      .sort('system:time_start')

    def get_img_coll(self):
        return self.coll


class GenericMultiImageDatasource(MultiImageDatasource):
    """Generic GlobalImageDatasource for loading an image by name or arguments to ee.ImageCollection."""

    def build_img_coll(self, coll_args=None):
        if coll_args is None:
            raise ValueError('coll_args must be provided, but was None')
        self.coll = ee.ImageCollection(coll_args) \
                      .filterDate(self.start_date, self.end_date) \
                      .filterBounds(self.filterpoly) \
                      .sort('system:time_start')

    def get_img_coll(self):
        return self.coll
