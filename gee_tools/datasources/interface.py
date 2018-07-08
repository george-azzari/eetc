"""Defines the interfaces for datasources"""
# Author: Anthony Perez

class DatasourceError(RuntimeError):
    pass

class ImageDatasource(object):
    """
    Defines the interface for optical datasources.

    Not all image datasources will be able to filter by geometry or date.
    See MultiImageDatasouerce GlobalImageDatasource and SingleImageDatasource
    """

    def build_img_coll(self, **kwargs):
        """Builds Image Collection(s)"""
        raise NotImplementedError

    def get_img_coll(self):
        """Returns one image collection"""
        raise NotImplementedError

class MultiImageDatasource(ImageDatasource):
    """
    Defines the datasource interface constructor for image collections where images
    can be filtered by both region and time.
    """

    def __init__(self, filterpoly, start_date, end_date, **kwargs):
        """
        filterpoly: ee.Geometry, used to filter the image collection
        start_date, end_date: strings of the form YYYY-MM-DD, used to filter the image collection
        kwargs are passed to self.build_img_coll
        """
        # TODO check arguments for errors
        self.filterpoly = filterpoly
        self.start_date = start_date
        self.end_date = end_date
        self.build_img_coll(**kwargs)

class GlobalImageDatasource(ImageDatasource):
    """
    Defines the datasource interface for image collections where images are global
    i.e. images collections where filtering by region does not make sense
    """

    def __init__(self, start_date, end_date, **kwargs):
        """
        start_date, end_date: strings of the form YYYY-MM-DD, used to filter the image collection
        kwargs are passed to self.build_img_coll
        """
        # TODO check arguments for errors
        self.start_date = start_date
        self.end_date = end_date
        self.build_img_coll(**kwargs)

class SingleImageDatasource(ImageDatasource):
    """
    Defines the datasource interface for image datasources corresponding to single images
    i.e. images collections where filtering by region or date does not make sense
    """

    def __init__(self, **kwargs):
        self.build_img_coll(**kwargs)
