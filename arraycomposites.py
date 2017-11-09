import ee


def show(desc, image, lon, lat):
    """
    Auxiliary function to inspect array-images at given pixel
    :param desc:
    :param image:
    :param lon:
    :param lat:
    :return:
    """
    print('==> ' + desc + ' <==')
    print(ee.ImageCollection([image]).getRegion(ee.Geometry.Point(lon, lat), 30))


def sorted_array(array, sortindex):
    """
    //Sort an array-image along the image axis (0) by values of bandname.
    :param array:
    :param sortindex:
    :return:
    """
    imageAxis = 0
    bandAxis = 1
    quality = array.arraySlice(bandAxis, sortindex, sortindex.add(1))
    values = array.arraySlice(bandAxis, 0)
    sorted = values.arraySort(quality)
    return sorted


def sorted_collection(collection, bandname):
    """
    Return a sorted array-image from a collection.
    Sorting is done along the image axis (0) by band called 'bandname'.
    :param collection:
    :param bandname:
    :return:
    """
    bandindex = ee.Image(collection.first()).bandNames().indexOf(bandname)
    array = collection.toArray()
    return sorted_array(array, bandindex)


def nthquantile_range(collection, bandname, quantile):
    """
    Return a sorted array-image from a collection.
    Sorting is done along the image axis (0) by band called 'bandname',
    but only values included between the minimum (first array element)
    and the given quantile (last element) are included.
    :param collection:
    :param bandname:
    :param quantile:
    :return:
    """
    imageAxis = 0
    sorted = sorted_collection(collection, bandname)
    imageCount = sorted.arrayLengths().arrayGet(imageAxis)
    percIndex = imageCount.multiply(quantile).toInt()
    lowestperc = sorted.arraySlice(imageAxis, 0, percIndex.add(1))
    return lowestperc


def nthquantile_mosaic(collection, bandname, quantile):
    """
    Return a mosaic built by sorting each stack of pixels by the input band
    in ascending order, and taking the image corresponding to the specified
    quantile.
    :param collection:
    :param bandname:
    :param quantile:
    :return:
    """
    lowestperc = nthquantile_range(collection, bandname, quantile)
    imageAxis = 0
    bandAxis = 1
    nthperc = lowestperc.arraySlice(imageAxis, -1)  # get last element in quantile range
    bandNames = collection.min().bandNames().slice(0)
    return nthperc.arrayProject([bandAxis]).arrayFlatten([bandNames])


def reducedrange_mosaic(collection, bandname, quantile, reducer):
    """
    Return a mosaic built by reducing the min:nth-quantile range of band 'bandname' with given reducer.
    :param collection:
    :param bandname:
    :param quantile:
    :param reducer:
    :return:
    """
    imageAxis = 0
    bandAxis = 1
    lowestperc = nthquantile_range(collection, bandname, quantile)
    reduced = lowestperc.arrayReduce(reducer, [imageAxis])
    bandNames = collection.min().bandNames().slice(0)
    return reduced.arrayProject([bandAxis]).arrayFlatten([bandNames])

