import ee
import json
import googlecloud
import re


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
    eeasset_path = eefolderpath+eeasset_name
    request = {"id": eeasset_path,
               "tilesets": [
                   {"sources": [
                       {"primaryPath": gsfilepath,
                        "additionalPaths": []}
                   ]}
               ],
               "bands": [],
               "reductionPolicy": "MEAN",
               "missingData": {"value": nodata}}
    taskid = ee.data.newTaskId(1)[0]
    t = ee.data.startIngestion(taskid, request)
    return t


def auto_upload(gsbucket, eefolderpath, nodata=-32768):
    """
    Upload and ingest all assets contained in Google Cloud Storage bucket.
    :param gsbucket:
    :param eefolderpath:
    :param nodata:
    :return:
    """
    bucketlist = googlecloud.list_bucket(gsbucket)
    filenames = [b['name'] for b in bucketlist]
    for f in filenames:
        gsfilepath = "gs://{0}/{1}".format(gsbucket, f)
        upload_asset_core(gsfilepath, eefolderpath, nodata=nodata)


def autoupdate_assets(eefolderpath, yearindx=-1, stateindx=0, readers=['dblobell@gmail.com']):
    """
    This is based on this type of asset id:
    'users/georgeazzari/scym_usa_v0/illinois_maize_yield_2000'
    """
    assets_dicts = ee.data.getList(dict(id=eefolderpath))
    for ad in assets_dicts:
        info = ad['id'].split('/')[-1].split('_')
        if yearindx is not None:
            ee.data.setAssetProperties(ad['id'], dict(year=int(info[yearindx])))
        if stateindx is not None:
            ee.data.setAssetProperties(ad['id'], dict(state=info[stateindx]))
        if readers is not None:
            d = dict(writers=[], readers=readers)
            j = json.dumps(d)
            ee.data.setAssetAcl(ad['id'], j)


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