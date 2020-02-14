from calendar import monthrange

import ee


def get_closest_to_date(img_coll, date):

    def closest(current_scene, closest_scene):
        current_scene = ee.Image(current_scene)
        closest_scene = ee.Image(closest_scene)
        current_date = ee.Date(current_scene.date())
        closest_date = ee.Date(closest_scene.date())
        curr_diff = ee.Number(date.difference(current_date, 'day')).abs()
        closest_diff = ee.Number(date.difference(closest_date, 'day')).abs()
        return ee.Algorithms.If(closest_diff.gte(curr_diff),
                                current_scene,
                                closest_scene)

    return ee.Image(img_coll.iterate(closest, img_coll.first()))


def get_center_date(start_date, end_date):
    start_date = ee.Date(start_date)
    end_date = ee.Date(end_date)
    diff = ee.Number(end_date.difference(start_date, 'day'))
    advance = ee.Number(diff).divide(ee.Number(2.0))
    return start_date.advance(advance, 'day')


def compute_monthly_seasonal_annual_and_long_running(base_img_coll, long_running_years, season_to_int_months, process_img_coll):
    """
    Args:
        base_img_coll (ee.ImageCollection):
        long_running_years (List[int]):  The list of years to include in the long running average.
            When computing the overall long running average, only the min and max are used.
            The full list is used when computing monthly or seasonal averages.
        season_to_int_months (Dict[str, List[int]]):  A dict from season names to the months in the season (as integers)
        process_img_coll (Callable[ [ee.ImageCollection, str], Tuple[ee.Image, List[str]]]):
            A function mapping an image collection and a string suffix to
            an output image and a set of bands in that image.
            This will be called for each time slice with the appropriate suffix
            e.g. monthly, seasonally, annually, long_running_monthly, long_running_seasonally, long_running_annually
    Returns:
        (Tuple[ee.Image, List[str]]):  Final image and bands
    """
    bands = []
    imgs = []

    def _process_img_coll(img_coll, suffix):
        img, _bands = process_img_coll(img_coll, suffix)
        imgs.append(img)
        bands.extend(_bands)

    for month in range(1, 12 + 1):
        long_running_img_coll = None

        for year in long_running_years:
            img_coll = base_img_coll.filterDate(
                '{}-{:02d}-01'.format(year, month),
                '{}-{:02d}-{}'.format(year, month, monthrange(year, month)[1])
            )

            if long_running_img_coll is None:
                long_running_img_coll = img_coll
            else:
                long_running_img_coll = long_running_img_coll.merge(img_coll)

        _process_img_coll(long_running_img_coll, 'month_{}_long_running_average'.format(month))

    for season, months in season_to_int_months.items():
        long_running_img_coll = None

        for year in long_running_years:
            img_coll = base_img_coll.filterDate(
                '{}-{:02d}-01'.format(year, min(months)),
                '{}-{:02d}-{}'.format(year, max(months), monthrange(year, max(months))[1])
            )

            if long_running_img_coll is None:
                long_running_img_coll = img_coll
            else:
                long_running_img_coll = long_running_img_coll.merge(img_coll)

        _process_img_coll(long_running_img_coll, '{}_long_running_average'.format(season))

    long_running_img_coll = base_img_coll.filterDate(
        '{}-01-01'.format(min(long_running_years)),
        '{}-12-31'.format(max(long_running_years)),
    )
    _process_img_coll(long_running_img_coll, 'long_running_average')

    # Final computation
    final_img = imgs.pop()
    for img in imgs:
        final_img = final_img.addBands(img)

    return final_img, bands
