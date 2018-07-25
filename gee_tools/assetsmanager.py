import ee
import json
import re
from googlecloud import list_objects as lobjs


def upload_asset_core(gsfilepath, eefolderpath, nodata=-32768):
    """
    Upload and ingest asset in EE starting from its path in Google Cloud Storage.
    :param gsfilepath: complete file path to file, e.g. "gs://mybucket/mybucketdir/myfile.tif"
    :param eefolderpath: has to end with `/`, e.g. "users/georgeazzari/eedir/"
    :param nodata: no data value for masking during ingestion
    :return: ingestion task
    """
    fname = gsfilepath.split('/')[-1]
    eeasset_name = fname.split('.')[0]  # getting rid of file extension
    if eefolderpath[-1] == '/':
        eeasset_path = eefolderpath+eeasset_name
    else:
        eeasset_path = eefolderpath + '/' + eeasset_name

    print 'Uploading from GCP source: ' + gsfilepath
    print 'Ingesting to EE asset: ' + eeasset_path

    request = {"id": eeasset_path,
               "tilesets": [
                   {"sources": [
                       {"primaryPath": gsfilepath,
                        "additionalPaths": []}
                   ]}
               ],
               "bands": [],
               "pyramidingPolicy": "MEAN",
               "missingData": {"value": nodata}}

    taskid = ee.data.newTaskId(1)[0]
    t = ee.data.startIngestion(taskid, request)

    return t


def auto_upload(gsbucket, gsprefix, eefolderpath, nodata=-32768):
    """
    Upload and ingest all assets contained in Google Cloud Storage bucket.
    :param gsbucket:
    :param gsprefix:
    :param eefolderpath:
    :param nodata:
    :return:
    """

    fx = lambda f: (f.split('.')[-1] == 'tiff') and (gsprefix in f)
    bucketlist = lobjs.list_bucket(gsbucket)
    filenames = [b['name'] for b in bucketlist if fx(b['name'])]
    tasks = []

    for f in filenames:
        gsfilepath = "gs://{0}/{1}".format(gsbucket, f)
        t = upload_asset_core(gsfilepath, eefolderpath, nodata=nodata)
        tasks.append(t)

    return tasks


def autoupdate_assets(eefolderpath, propdict, sep):
    """
    Update the properties of all assets in 'eefolderpath' based on their file name. Mapping between file name
    chunks (separated by 'sep') and asset properties must be provided in 'propdict'.
    :param eefolderpath:
    :param propdict:
    :param sep:
    :return:
    """

    assets_dicts = ee.data.getList(dict(id=eefolderpath))

    for ad in assets_dicts:

        info = ad['id'].split('/')[-1].split(sep)
        print('Update properties of asset ' + ad['id'])

        for k in propdict.keys():

            info_k = info[int(k)]
            prop = propdict[k]
            ee.data.setAssetProperties(ad['id'], {prop: info_k})


def autoupdate_wrsassets(eefolderpath, yearindx=-2, wrsindx=2, readers=['dblobell@gmail.com']):
    """
    Note: eefolderpath has to be without final '/', e.g. 'users/georgeazzari/mydir'
    """
    assets_dicts = ee.data.getList(dict(id=eefolderpath))
    for ad in assets_dicts:
        info = ad['id'].split('/')[-1].split('_')
        wrs = re.split('(\d+)', info[wrsindx])
        ee.data.setAssetProperties(ad['id'], dict(year=int(info[yearindx]), WRS_PATH=int(wrs[1]), WRS_ROW=int(wrs[3]),
                                                  WRSLAB="P{0}R{1}".format(str(wrs[1]), str(wrs[3]))))
        if readers is not None:
            d = dict(writers=[], readers=readers)
            j = json.dumps(d)
            ee.data.setAssetAcl(ad['id'], j)


def copy_to(sourcepath, destpath, delete_originals=False):
    """
    Example: move_tocoll('users/georgeazzari/scym_usa_seamless_v2b',
                         'users/georgeazzari/scym_usa_seamless_alongpaths_v2b_coll')
    :param sourcepath:
    :param destpath:
    :param delete_originals:
    :return:
    """
    assets_dicts = ee.data.getList(dict(id=sourcepath))
    for ad in assets_dicts:
        ee.data.copyAsset(ad['id'], destpath + '/' + ad['id'].split('/')[-1])
    if delete_originals:
        for ad in assets_dicts:
            ee.data.deleteAsset(ad['id'])