"""
Author: Anthony Perez

Implements the logic necessary to take an arbitrary FeatureCollection with point geomerty and gather tiles around it for down stream tasks.
"""
import ee

from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, ImageDatasource, DatasourceError
from gee_tools import imgtools
from gee_tools.ai_io import ee_tf_exports as tfexp
from gee_tools.exports import constants, util


class ImageSpec(object):
    """
    A configuration object that outlines the components needed to create a scene and provides a helper
        function to create the corresponding scene.

    These components are:
    1) A set of images collections                             ee.ImageCollection
    2) Their compositing methods                               Function from ee.ImageCollection to ee.Image
    3) A compositing time period (i.e. start and end dates)    Pair of python datetime.date objects
    4) A compositing region (i.e. Africa bounding box)         ee.Geometry
    5) Projection and scale                                    string (projection id) and number (meters)
    """

    def __init__(self, start_date, end_date, filterpoly, scale, projection=constants.EPSG3857):
        """
        The config object should contain values for
        :param start_date:  Passed to MultiImageDatasource and GlobalImageDatasource 
            classes to filter by date
        :type string_date: string or date or datetime or ee.Date
        :param end_date:  See start_date
        :param filterpoly: Passed to MultiImageDatasource classes to filter by geometry / region
        :type filterpoly: ee.Geometry
        :param projection: A projection (crs) designation string.  All images are reprojected 
            using this projection. (default EPSG:3857)
        :type projection: string
        :param scale:  The scale parameter passed to the reproject method
        :type scale: int
        """

        self.start_date = start_date
        self.end_date = end_date
        self.region = filterpoly
        self.projection = projection
        self.scale = scale

        if not isinstance(self.region, ee.Geometry):
            raise ValueError("Expected filterpoly to be an ee.Geometry instance.  Got: {}".format(type(self.region)))

        if not (isinstance(self.start_date, ee.Date) or isinstance(self.end_date, ee.Date)):
            util.start_date_before_end(self.start_date, self.end_date)
        if not isinstance(self.start_date, ee.Date):
            self.start_date = util.date_to_str(self.start_date)
        if not isinstance(self.end_date, ee.Date):
            self.end_date = util.date_to_str(self.end_date)

        self.data_sources = []
        self._static_scenes = []
        self.scene = None
        self._scene_has_latlon = False  # Consider moving to constructor arguments.

    def add_datasource(self, datasource_class, composite_function, ds_kwargs=None):
        """
        Add a data source to this Image Specification
        :param datasource_class: A class that inherits from ImageDatasource.
            Note that this should not be an instance of the class, but rather the class itself.
            Must be one of MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource
        :param composite_function:  The desired compositing function used to transform the datasource
            into a singe image.  The function must take an ee.ImageCollection
            as an argument and return an ee.Image
        :param ds_kwargs:  Key word arguments passed to the datasource_class at construction.
        """
        if ds_kwargs is None:
            ds_kwargs = {}
        if not issubclass(datasource_class, ImageDatasource):
            raise ValueError("data_sources contained non ImageDatasource element: {}".format(datasource_class))
        if not isinstance(datasource_class, type):
            raise ValueError("data_sources elements should be classes, not instances: {}".format(datasource_class))
        self.data_sources.append((datasource_class, ds_kwargs, composite_function))
        self.scene = None  # Reset the scene so that it must be recomputed

    def add_static_scene(self, scene):
        """
        Adds an ee.Image to the scene computed by this ImageSpec.
        The added image will not by filtered by date or region, but
        will be reprojected according to the constructor arguments.

        Args:
            scene (ee.Image):  The scene to add.
        """
        self._static_scenes.append(scene)
        self.scene = None

    def get_scene(self, add_latlon=True):
        """
        Return the scene corresponding to this ImageSpec instance.
        """
        if len(self.data_sources) == 0:
            raise ValueError("Empty data_sources.  No data sources were specified.")

        if add_latlon != self._scene_has_latlon:
            self.scene = None
            self._scene_has_latlon = add_latlon

        if self.scene is None:
            self.scene = self._get_scene(self, add_latlon=add_latlon)
        return self.scene

    def set_specification(self, start_date=None, end_date=None,
                          filterpoly=None, scale=None,
                          projection=None):
        """
        Change the constructor arguments to this ImageSpec

        Args:
            start_date (Optional[Union[str, ee.String, ee.Date]]): Override the default start_date.
            end_date (Optional[Union[str, ee.String, ee.Date]]): Override the default end_Date.
            filterpoly (Optional[ee.Geometry]): Override the default filterpoly.
            scale (Optional[Union[int, float, ee.Number]]): Override the default scale.
            projection (Optional[Union[str, ee.String]]): Override the default projection/CRS.
        """
        self.start_date = self.start_date if start_date is None else start_date
        self.end_date = self.end_date if end_date is None else end_date
        self.region = self.region if filterpoly is None else filterpoly
        self.projection = self.projection if projection is None else projection
        self.scale = self.scale if scale is None else scale
        self.scene = None

    @staticmethod
    def _get_scene(image_spec, add_latlon=True, error_check=False):
        """Return the scene corresponding to an ImageSpec instance"""
        processed_imagery = []
        start, end = image_spec.start_date, image_spec.end_date
        for data_source in image_spec.data_sources:
            ds_class, kwargs, comp_fn = data_source

            if issubclass(ds_class, MultiImageDatasource):
                ds_args = (image_spec.region, start, end)
            elif issubclass(ds_class, GlobalImageDatasource):
                ds_args = (start, end)
            elif issubclass(ds_class, SingleImageDatasource):
                ds_args = ()
                if comp_fn is None:
                    comp_fn = lambda img_coll: ee.Image(img_coll.first())
            else:
                raise ValueError("Invalid image_spec.  {} was not a known datasource type".format(ds_class))

            ds = ds_class(*ds_args, **kwargs)
            img_coll = ds.get_img_coll()
            img = comp_fn(img_coll)
            img = img.reproject(image_spec.projection, None, image_spec.scale)
            processed_imagery.append(img)

            if error_check:
                util.check_empty_bands(img)

        for img in image_spec._static_scenes:

            img = img.reproject(image_spec.projection, None, image_spec.scale)
            processed_imagery.append(img)

            if error_check:
                util.check_empty_bands(img)

        if len(processed_imagery) == 0:
            # This should be unreachable since len(data_source) > 0
            raise ValueError("No Imagery Found")

        scene = processed_imagery.pop(0)
        for bands in processed_imagery:
            scene = scene.addBands(bands)
        if add_latlon:
            scene = ee.Algorithms.If(
                scene.bandNames().size().eq(0),
                scene,
                imgtools.add_latlon(scene)
            )
            scene = ee.Image(scene)
        return scene


def add_imagery_scene(ft, scene, scale, projection, output_size):
    """
    Take a featureCollection (tf) where each row has point geometry and 
    return a featureCollection with image bands added.

    :param ft: ee.FeatureCollection, each row must have point geometry
    :type ft: ee.FeatureCollection
    :param scale: The projection's scale
    :type scale: int
    :param projection: The projection's string designation
    :type string:
    :param output_size: The outputsize in pixels (final output is an (2 * output_size + 1) by (2 * output_size + 1) square)
    :type output_size: int

    :returns: a new feature collection with an output_size by output_size tile added to each row.
    The tile's bands are stored in separate columns.
    """
    if projection != constants.EPSG3857:
        # gee_tools.ai_io.ee_tf_exports.get_array_patches currently only samples with projection EPSG3857
        raise NotImplementedError("Projection must be EPSG:3857")

    ft = ee.Algorithms.If(
        scene.bandNames().size().eq(0),
        ft,
        tfexp.get_array_patches(
            scene, scale, output_size, ft,
            False, False, None, None,
            None, None, None
        )
    )
    ft = ee.FeatureCollection(ft)

    return ft


def add_imagery(ft, image_spec, output_size, add_latlon=True):
    """
    Take a featureCollection (tf) where each row has point geometry and
    return a featureCollection with image bands added.

    Args:
        ft (ee.FeatureCollection):  A feature collection, each feature must have point geometry.
        image_spec (ImageSpec):  The image spec to sample imagery from.
        output_size (int): The outputsize in pixels (final output is an (2 * output_size + 1) by (2 * output_size + 1) square)

    Returns:
        (ee.FeatureCollection): A new feature collection with an output_size by output_size tile added to each row.
        The tile's bands are stored in separate columns.
    """
    scene = image_spec.get_scene(add_latlon=add_latlon)
    return add_imagery_scene(ft, scene, image_spec.scale, image_spec.projection, output_size)
