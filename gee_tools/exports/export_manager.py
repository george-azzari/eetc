"""
Author: Anthony Perez

Defines the ExportManager class
"""
import logging
from enum import Enum as PythonEnum

import ee

from gee_tools.datasources.generic_datasources import GenericSingleImageDatasource
from gee_tools.datasources.interface import SingleImageDatasource
from gee_tools.exports.task_scheduler import TaskScheduler
from gee_tools.exports.image_spec import ImageSpec, add_imagery
from gee_tools.assetsmanager import asset_exists

logger = logging.getLogger(__name__)


class ExportManagerError(RuntimeError):
    pass


class ExportManager(object):
    """
    The purpose of this class is to consume one
    large configuration file describing a desired scene
    and then export that scene or sample tiles under
    different configurations.
    """

    def __init__(self, datasources_config):
        """
        Args:
            datasources_config (Dict[str, Dict[str, Any]]):
                A dictionary specifying a mapping from datasource names
                to datasource configurations.  Each data source configuration
                must have the following format:
                    {
                        "class": (Union[MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource]),
                        "args": (Optiona[Dict[str, Any]])  Will be passed to the class constructor as kwargs
                            after filterpoly, start_date, and end_date as appropriate,
                        "composite_fn": (Callable[[ee.ImageCollection], ee.Image])  A function to convert the image collection to an image.
                            If None and the class inherits from SingleImageDatasource then .first() will be used.,
                        "bands": (List[str])  The list of expected band names,
                        "tag": (Optional[Union[str, Enum, List[str], List[Enum]]])  An ID or set of IDs
                            used for filtering as described in public methods.  Defaults to [].,
                        "cache_asset_id": (Optional[str])  If this is present and the class inherits from SingleImageDatasource
                            then the output will be saved to the asset ID represented by this string before exports occur.  This
                            can be useful for jobs that would otherwise run out of memory.,
                    }

                {
                    "landsat": {
                        "class": optx.LandsatSR,
                        "args": {},
                        "composite_fn": cmp_fns.select_median(LANDSATSR_COMMON_BANDS),
                        "bands": LANDSATSR_COMMON_BANDS,
                        "tag": 'CNN_TILES',
                    },
                    "lm_worldpop": {
                        "class": GenericSingleImageDatasource,
                        "args": { 'image_args': ee.Image("projects/atlasaipbc/clients/world_bank_et/inputs/worldpop_v2_adj/AFR_PPP_2015_adj_v2").select([0], ['WORLD_POP']) },
                        "composite_fn": None,
                        "bands": ['WORLD_POP'],
                        "tag": 'LM_FEATURES',
                        # caching does not make sense here, but as an example if the composite function has intense featurization this may be useful.
                        "cache_asset_id": 'projects/atlasaipbc/clients/world_bank_et/linear_model_features/lm_worldpop',
                    },
                }
        """
        self.datasources_config = dict(datasources_config)

        for input_name in list(self.datasources_config.keys()):
            input_config = self.datasources_config[input_name]

            # Create shallow copy with defaults.
            input_config = dict(input_config)
            if 'tag' not in input_config:
                input_config['tag'] = []

            in_tag = input_config['tag']

            if isinstance(in_tag, list):
                input_config['tag'] = set(in_tag)
            elif isinstance(in_tag, str) or isinstance(in_tag, PythonEnum):
                input_config['tag'] = set([in_tag])
            else:
                raise ExportManagerError('Unrecognized tag type (must be str, Enum, or List): {}'.format(in_tag))

            self.datasources_config[input_name] = input_config


    def _get_datasources_by_tag(self, tags=None):
        """
        Args:
            tags (Optiona[Iterable[Union[str, Enum]]]):
        Returns:
            (Dict[str, Dict[str, Any]]): self.datasources_config filtered by tags.
                If tags is None, return all datasource configs.
        """
        data_sources = self.datasources_config
        if tags is None:
            return data_sources
        tags = set(tags)

        def matches_tags(input_config):
            in_tags = input_config['tag']
            return len(in_tags.intersection(tags)) > 0

        return {
            input_name: input_config
            for input_name, input_config in data_sources.items()
            if matches_tags(input_config)
        }


    @staticmethod
    def _convert_to_image_spec(image_spec):
        """
        Convert a dict or an ImageSpec into an ImageSpec
        """
        if isinstance(image_spec, ImageSpec):
            return image_spec
        else:
            return ImageSpec(**image_spec)


    @staticmethod
    def _populate_cache(image_spec, datasources):
        scheduler = TaskScheduler()
        export_region = None

        image_spec_kwargs = {
            'start_date': image_spec.start_date,
            'end_date': image_spec.end_date,
            'filterpoly': image_spec.region,
            'projection': image_spec.projection,
            'scale': image_spec.scale,
        }

        for input_name, input_config in datasources.items():

            output_asset_id = input_config.get('cache_asset_id', None)
            if output_asset_id is None:
                continue

            if asset_exists(output_asset_id):
                logger.warn('{} already exists.  Will not precompute.'.format(output_asset_id))
                continue

            if not issubclass(input_config['class'], SingleImageDatasource):
                raise ExportManagerError(
                    'Cannot cache {}.  The provided class '
                    'is not a SingleImageDatasource.'.format(input_name)
                )

            image_spec = ImageSpec(**image_spec_kwargs)
            image_spec.add_datasource(
                datasource_class=input_config['class'],
                composite_function=input_config['composite_fn'],
                ds_kwargs=input_config['args'],
            )

            scene = image_spec.get_scene(add_latlon=False)
            scene = scene.clip(image_spec_kwargs['filterpoly'])

            if export_region is None:
                export_region = image_spec.region.bounds().getInfo()['coordinates']

            task = ee.batch.Export.image.toAsset(**{
                'image': scene,
                'description': input_name,
                'assetId': output_asset_id,
                'region': export_region,
                'scale': image_spec_kwargs['scale'],
                'crs': image_spec_kwargs['projection'],
                'maxPixels': 1e13,
            })
            logger.info('Will precompute {}'.format(input_name))
            scheduler.add_task(task, output_asset_id)

        if len(scheduler) > 0:
            scheduler.run(verbose=999, error_on_fail=True)

    
    @staticmethod
    def _populate_image_spec(image_spec, datasources):
        output_bands = []
        for _, input_config in datasources.items():

            cache_asset_id = input_config.get('cache_asset_id', None)

            if cache_asset_id is None:
                image_spec.add_datasource(
                    datasource_class=input_config['class'],
                    composite_function=input_config['composite_fn'],
                    ds_kwargs=input_config['args'],
                )
            else:
                image_spec.add_datasource(
                    datasource_class=GenericSingleImageDatasource,
                    composite_function=None,
                    ds_kwargs={ 'image_args': cache_asset_id },
                )

            output_bands.extend(input_config['bands'])

        if 'LAT' in output_bands or 'LON' in output_bands:
            raise ValueError(
                'Output bands contained the key LAT or LON.  '
                'This is not allowed, LAT and LON will be added automatically.'
            )

        output_bands += [u'LAT', u'LON']
        return output_bands


    def _get_image_spec_helper(self, image_spec, tags=None):
        datasources = self._get_datasources_by_tag(tags=tags)
        image_spec = ExportManager._convert_to_image_spec(image_spec)
        ExportManager._populate_cache(image_spec, datasources)
        output_bands = ExportManager._populate_image_spec(image_spec, datasources)
        return image_spec, output_bands


    def get_scene(self, image_spec, tags=None):
        """
        Take a featureCollection (fc) where each row has point geometry and 
        return a featureCollection with image bands added.

        Args:
            image_spect (Union[ImageSpec, Dict[str, Any]]):  Either an ImageSpec instance or a dict formated as follows:
                {
                    'start_date': Union[ee.Date, str],
                    'end_date': Union[ee.Date, str],
                    'filterpoly': ee.Geoemtry,
                    'projection': str,  # CRS
                    'scale': Union[int, float],
                }
                If an ImageSpec is passed, datasources already in the ImageSpec instance will be included in the output.
            tags (Optional[Union[str, Enum, Collection[str, Enum]]]):  A collection of tags matching those passed in
                the 'datasources_config' constructor argument.  Only tags contained in the tags argument
                will be used when generating the return value.  If None, all datasources are included.
                Defaults to None.

        Returns:
            (Tuple[ee.Image, List[str]]):
                First element: The image represented by combining datasources according to
                the specifications in the image_spec argument.
                Second element:  The list of output bands.

            If a collection in the constructor argument is filtered in such a way that it becomes the empty
            collection, it's bands will be omitted from the output but will still be included in the second
            return element.
        """
        image_spec, output_bands = self._get_image_spec_helper(image_spec, tags)
        return image_spec.get_scene(), output_bands


    def sample_tiles(self, fc, image_spec, export_radius, tags=None):
        """
        Take a featureCollection (fc) where each row has point geometry and 
        return a featureCollection with image bands added.

        Args:
            fc (ee.FeatureCollection):  A feature collection, all features must have point geometries.
            image_spect (Union[ImageSpec, Dict[str, Any]]):  Either an ImageSpec instance or a dict formated as follows:
                {
                    'start_date': Union[ee.Date, str],
                    'end_date': Union[ee.Date, str],
                    'filterpoly': ee.Geoemtry,
                    'projection': str,  # CRS
                    'scale': Union[int, float],
                }
                If an ImageSpec is passed, datasources already in the ImageSpec instance will be included in the output.
            tags (Optional[Union[str, Enum, Collection[str, Enum]]]):  A collection of tags matching those passed in
                the 'datasources_config' constructor argument.  Only tags contained in the tags argument
                will be used when generating the return value.  If None, all datasources are included.
                Defaults to None.
            export_radius (int): The outputsize in pixels (final output is an (2 * output_size + 1) by (2 * output_size + 1) square)

        Returns:
            (Tuple[ee.FeatuerCollection, List[str]]): 
            First element: A new feature collection with an output_size by output_size tile added to each row.
            The tile's bands are stored in separate columns.
            Second element:  The list of output bands.

            If a collection in the constructor argument is filtered in such a way that it becomes the empty
            collection, it's bands will be omitted from the output but will still be included in the second
            return element.
        """
        image_spec, output_bands = self._get_image_spec_helper(image_spec, tags)
        fc = add_imagery(fc, image_spec, output_size=export_radius)
        return fc, output_bands
