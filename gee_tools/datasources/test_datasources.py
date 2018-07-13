"""nonoptical_datasources unit test"""
from nonoptical_datasources import DMSPCalV4, VIIRSMonthlyStrCorr, SRTMElevation, Palsar, DMSPCalVIIRSJoined 
import ee

class constants:
    EPSG3857 = 'EPSG:3857'

def test1(start_date, end_date, use_viirs, test_palsar=True):
    scale = 30
    if use_viirs:
        nl = VIIRSMonthlyStrCorr(start_date, end_date).viirs
    else:
        nl = DMSPCalV4(start_date, end_date).dmsp
    nl = nl.median().reproject(constants.EPSG3857, None, scale)

    topo = SRTMElevation().topo.reproject(constants.EPSG3857, None, scale)
    palsar = Palsar(start_date, end_date).palsar.mean().reproject(constants.EPSG3857, None, scale)
    if len(nl.bandNames().getInfo()) != 1:
        raise ValueError("Invalid nightlights bands")
    if len(topo.bandNames().getInfo()) != 3:
        raise ValueError("Invalue elevation bands")
    palsar_bands = palsar.bandNames().getInfo()
    if test_palsar and len(palsar_bands) != 5:
        raise ValueError("Invalid pulsar bands\n{}".format(palsar_bands))

def test2():
    nl = DMSPCalVIIRSJoined('2009-1-1', '2011-12-31').get_img_coll()

    size = nl.size().getInfo()
    if size != 2:
        raise ValueError("DMSP VIIRS Joined should have 3 image for the '2009-1-1', '2011-12-31' period.  Had: {}".format(size))

    nl = DMSPCalVIIRSJoined('2012-1-1', '2014-12-31').get_img_coll()

    if nl.size().getInfo() != 12:
        raise ValueError("DMSP VIIRS Joined should have 1 image for the '2012-1-1', '2014-12-31' period")

    nl = DMSPCalVIIRSJoined('2015-1-1', '2017-12-31').get_img_coll()

    if nl.size().getInfo() != 36:
        raise ValueError("DMSP VIIRS Joined should have 3 images for the '2015-1-1', '2017-12-31' period")

    nl = DMSPCalVIIRSJoined('2012-1-1', '2013-12-31').get_img_coll()
    if nl.size().getInfo() != 0:
        raise ValueError("No imagery should exist for the '2012-1-1', '2013-12-31' period for DMSP VIIRS Joined")

    nl = DMSPCalVIIRSJoined('2015-1-1', '2015-12-31').get_img_coll().median().getInfo()


if __name__ == "__main__":
    ee.Initialize()
    test1('2009-1-1', '2011-12-31', False)
    test1('2012-1-1', '2014-12-31', True, test_palsar=False)
    test1('2015-1-1', '2017-12-31', True)
    test2()
