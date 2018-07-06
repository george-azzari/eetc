"""
Author: George Azzari
"""
import subprocess
from osgeo import gdal, ogr, osr
import os
import sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import copy
import matplotlib.mlab as mlab
from matplotlib import rc
gdal.AllRegister()


class GeoProps:
    def __init__(self):
        self.eDT = None
        self.Proj = None
        self.GeoTransf = None
        self.Driver = None
        self.Flag = False
        #self.pycuda_flag = False
        
    def SetGeoProps(self, eDT, Proj, GeoTransf, Driver):
        self.eDT = eDT
        self.Driver = Driver
        self.Proj = Proj
        self.GeoTransf = GeoTransf
        self.Flag = True
        
    def copy_geoprops(self, component):
        component.open_gdal()
        self.import_geogdal(component.gdal_dataset)
        #Component.CloseGdal()

    def import_geogdal(self, gdal_dataset):
        self.eDT = gdal_dataset.GetRasterBand(1).DataType
        self.Proj = gdal_dataset.GetProjection()
        self.GeoTransf = gdal_dataset.GetGeoTransform()
        self.Driver = gdal.GetDriverByName("GTiff")
        self.xOrigin = self.GeoTransf[0]
        self.yOrigin = self.GeoTransf[3]
        self.pixelWidth = self.GeoTransf[1] 
        self.pixelHeight = self.GeoTransf[5]
        self.srs = osr.SpatialReference()
        self.srs.ImportFromWkt(self.Proj)
        self.srsLatLon = self.srs.CloneGeogCS()
        self.Flag = True

    def get_affinecoord(self, geolon, geolat):
        """Returns coordinates in meters (affine) from degrees coordinates (georeferenced)"""
        ct = osr.CoordinateTransformation(self.srsLatLon, self.srs)
        tr = ct.TransformPoint(geolon, geolat)
        xlin = tr[0]
        ylin = tr[1]
        return xlin, ylin

    def get_georefcoord(self, xlin, ylin):
        """Returns coordinates in degrees (georeferenced) from coordinates in meters (affine)"""
        ct = osr.CoordinateTransformation(self.srs, self.srsLatLon)
        tr = ct.TransformPoint(xlin, ylin)
        geolon = tr[0]
        geolat = tr[1]
        return geolon, geolat

    def lonlat2colrow(self, lon, lat):
        """ Returns the (col, row) of a pixel given its coordinates (in meters)"""
        col = int((lon - self.xOrigin) / self.pixelWidth)
        row = int((lat - self.yOrigin) / self.pixelHeight)
        #print "(long,lat) = (",GeoX, ",", GeoY,") --> (col,row) = (",xOffset,",",yOffset,")"
        #NOTE: watch out! if you're using this to read a 2D np.array, remember 
        #that xOffset = col, yOffset = row
        return [col, row]
 
    def colrow2lonlat(self, col, row):
        """ Returns the (lon, lat) of a pixel given its (col, row)"""
        lon = col * self.pixelWidth + self.xOrigin
        lat = row * self.pixelHeight + self.yOrigin
        return [lon, lat]

    def get_center_coord(self, raster_array_shape, affine=False):
        """ Input: raster_array_shape is the output of gdalobject.np_array.shape, which is (#rows, #cols)
            Returns: coordinate (lon, lat) of the center of the raster."""
        s = raster_array_shape
        ul = self.colrow2lonlat(0, 0)
        lr = self.colrow2lonlat(s[1], s[0])
        lon_ext_m = lr[0] - ul[0]
        lat_ext_m = ul[1] - lr[1]
        lon_cntr = ul[0] + lon_ext_m/2
        lat_cntr = lr[1] + lat_ext_m/2
        if affine:
            return lon_cntr, lat_cntr
        if not affine:
            return self.get_georefcoord(lon_cntr, lat_cntr)

    def get_raster_extent(self, raster_array_shape):
        """ Input: raster_array_shape is the output of gdalobject.np_array.shape, which is (#rows, #cols)
            Returns: extent (in meters) of the raster."""
        s = raster_array_shape
        ul = self.colrow2lonlat(0, 0)
        lr = self.colrow2lonlat(s[1], s[0])
        lon_ext_m = lr[0] - ul[0]
        lat_ext_m = ul[1] - lr[1]
        return lon_ext_m, lat_ext_m

    def get_small_pxlwin(self, GeoX, GeoY, dpx, dpy):
        cen_col, cen_row = self.lonlat2colrow(GeoX, GeoY)
        rows = range(cen_row-dpx, cen_row+dpx+1, 1)
        columns = range(cen_col-dpy, cen_col+dpy+1, 1)
        row_indx = []
        col_indx = []
        for i in rows:
            for j in columns:
                row_indx.append(i)
                col_indx.append(j)
        return np.array(row_indx), np.array(col_indx)

    def get_polygon_corners(self, affine=True, **polyargs):
        """ Returns the extent of the polygon in georeferenced or affine coordinates.
            Georeferenced coordinates are directly read from the polygon file.
            Affine coordinates are calculated using the transformation set with the current instance of the class."""
        d = ogr.GetDriverByName(polyargs['driver'])
        poly = d.Open(polyargs['fpath'])
        layer = poly.GetLayerByName(polyargs['layer_name'])
        ext = layer.GetExtent()
        ulx, lrx, lry, uly = ext[0], ext[1], ext[2], ext[3]
        if not affine:
            return ulx, lrx, lry, uly
        if affine:
            lin_ul = self.get_affinecoord(ulx, uly)
            lin_lr = self.get_affinecoord(lrx, lry)
            top, bottom, left, right = lin_ul[1], lin_lr[1], lin_ul[0], lin_lr[0]
            return left, bottom, right, top  # this is the order for warping

    def get_polygon_center(self, affine=False, **polyargs):
        left, bottom, right, top = self.get_polygon_corners(affine=True, **polyargs)  # affine=True is right, think!
        lon_extent = right - left
        lat_extent = top - bottom
        aff_cntr_lon = left + lon_extent/2
        aff_cntr_lat = bottom + lat_extent/2
        if affine:
            return aff_cntr_lon, aff_cntr_lat
        if not affine:
            return self.get_georefcoord(aff_cntr_lon, aff_cntr_lat)


class GeoComponent:
    """This class should be maintained as much "pure" as possible. It is the very core for all rasters.
    The class is supposed to have only attributes/methods shared by ALL possible rasters/components in the library.
    The philosophy is: whatever is general enough to go here (could be ran by any component) SHOULD go here.
    It shouldn't be locked to only Landsat related stuff, but it should be open for future integration of more sensors"""
    def __init__(self, filepath):
        self.dirpath = "None"
        self.filename = "None"
        self.filepath = filepath
        self.label = "General GeoComponent"
        self.gdal_dataset = []
        self.geoprops = GeoProps()  # Initial GeoProps are empty
        self.array_flag = False
        self.gpuarray_flag = False
        self.gdal_dataset_flag = False
        self.init_flag = False
        self.write_flag = True  # By default writes component on disk
        self.wrp_win_flag = False
        self.plot_color = "Black"

    def _close_gdal(self):
        if self.gdal_dataset_flag:
            self.gdal_dataset = None
            self.gdal_dataset_flag = False

    def open_gdal(self):
        if not self.gdal_dataset_flag:
            if self.check_gdalfile():
                print "        + Reading {0} from last file...".format(self.label)
                self.gdal_dataset = gdal.Open(self.filepath)
                self.gdal_dataset_flag = True
                self.geoprops.import_geogdal(self.gdal_dataset)
            else:
                sys.exit("GeoComponent: couldn't find any {0} file in the requested path".format(self.label))

    def check_gdalfile(self):
        if not os.access(self.filepath, os.F_OK):  # Check if the destination folder exists already
            return False
        else:
            return True

    def get_data_array(self):
        self.open_gdal()
        array = self.gdal_dataset.GetRasterBand(1).ReadAsArray()
        self._close_gdal()
        return array

    def get_masked_array(self, **mask_args):
        """
        The input mask_args must include the following arguments:
        Mask arguments: "mask" (array), and "mask_value" (float)
        """
        masked = np.ma.masked_where(mask_args['mask'] == mask_args['mask_value'], self.get_data_array())
        return masked

    def get_nonmasked_values(self, **mask_args):
        """
        Returns 1D array of non masked values
        The input mask_args must include the following arguments:
        Mask arguments: "mask" (array), and "mask_value" (float)
        """
        masked = self.get_masked_array(**mask_args)
        nonmasked = np.ma.MaskedArray.compressed(masked)
        return nonmasked

    def get_threshold_mask(self, threshold):
        """"
        Returns a 2D array (mask) corresponding to pixels greater than a given threshold.
        """
        mask = np.zeros(self.get_data_array().shape)
        indx = np.where(self.get_data_array() > threshold)
        mask[indx] = 1
        return mask

    def get_threshold_maskargs(self, threshold):
        """
        Returns a mask dictionary resulting from selecting pixels greater than a given threshold.
        The mask dictionary includes the following:
        "masked" (bool), "mask" (array),"mask_value" (float), and "mask_id" (only if masked=True)
        """
        mask_args = dict(masked=True)
        mask_args.update(mask=self.get_threshold_mask(threshold))
        mask_args.update(mask_value=0)
        mask_args.update(mask_id="ThresholdMaksk_{0}".format(threshold))
        return mask_args

    def get_flat_array(self, **mask_args):
        """
        The input mask_args must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)
        """
        if mask_args["masked"]:
            data = self.get_nonmasked_values(**mask_args)
        elif not mask_args["masked"]:
            data = self.get_data_array().flat[:]
        return data

    def get_data_minmax(self, **mask_args):
        """ The input kwargs must include the following arguments:
            Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)"""
        data = self.get_flat_array(**mask_args)
        return np.min(data), np.max(data)

    def get_kde(self, **mask_args):
        """
        Gaussian Kernel data density (pdf) estimation (normalized?).
        Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)
        """
        data = self.get_flat_array(**mask_args)
        try:
            kde = stats.gaussian_kde(data)  # Calculate Kernel Density
            kdeisok = True
        except np.linalg.linalg.LinAlgError as err:
            if 'singular matrix' in err.message:
                kde = np.nan
                kdeisok = False
            else:
                sys.exit("Something went wrong when estimating KDE")
        kde_dict = dict(kde=kde, valid_kde=kdeisok)
        return kde_dict

    def get_kde_bins(self, **kwargs):
        """
        The input kwargs must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)
        Bins arguments: "default_bins", "min_bin", "max_bin" (only if default_bins=False).
        """
        if kwargs["default_bins"]:
            data = self.get_flat_array(**kwargs)
            m, M = self.get_data_minmax(**kwargs)
        else:
            m = kwargs["min_bin"]
            M = kwargs["max_bin"]
        kde_bins = np.arange(m, M, (M-m)/250.)  # support
        return kde_bins

    def evaluate_kde(self, **kwargs):
        """
        The input kwargs must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)
        Bins arguments: "default_bins", "min_bin", "max_bin" (only if default_bins=False).
        """
        kde_bins = self.get_kde_bins(**kwargs)
        kde_dict = self.get_kde(**kwargs)
        if kde_dict["valid_kde"]:
            kde = kde_dict["kde"].evaluate(kde_bins)
        else:
            kde = np.ones(kde_bins.shape) * np.nan
        kde_eval_dict = dict(kde_eval=kde, bins=kde_bins, valid_kde=kde_dict["valid_kde"])
        return kde_eval_dict

    def get_max_kde(self, **kwargs):
        """
        The input kwargs must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array), and "mask_value" (float) (only if masked=True)
        Bins arguments: "default_bins", "min_bin", "max_bin" (only if default_bins=False).
        """
        kde_eval_dict = self.evaluate_kde(**kwargs)
        if kde_eval_dict["valid_kde"]:
            max_indx = np.argmax(kde_eval_dict["kde_eval"])
            max_kde = kde_eval_dict["bins"][max_indx]
        else:
            max_kde = np.nan
        return max_kde

    def get_kde_plot(self, ax, **kwargs):
        """
        The input kwargs must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array),"mask_value" (float), and "mask_id" (only if masked=True)
        Bins arguments: "default_bins", "min_bin", "max_bin" (only if default_bins=False).
        """
        kde_eval_dict = self.evaluate_kde(**kwargs)
        line_args = dict(label="{0} KDE".format(self.label), color=self.plot_color, lw=1.5)
        kde_pl = ax.plot(kde_eval_dict["bins"], kde_eval_dict["kde_eval"], **line_args)
        kde_fill = fill_args = dict(color=self.plot_color, facecolor=self.plot_color, alpha=0.2)
        ax.fill_between(kde_eval_dict["bins"], 0.0, kde_eval_dict["kde_eval"], **fill_args)
        #vline = vline_args = dict(ymin=0, ymax=kde_max/self.limits.ymax, color=colors[j], ls=":", lw=1.5)
        #ax.axvline(kde.pdf_peak_bin, **vline_args)
        return kde_pl, kde_fill,  # vline

    def plot_kde(self, figw=12, figh=5, **kwargs):
        """
        The input kwargs must include the following arguments:
        Mask arguments: "masked" (bool), "mask" (array),"mask_value" (float), and "mask_id" (only if masked=True)
        Bins arguments: "default_bins", "min_bin", "max_bin" (only if default_bins=False).
        """
        fig, ax = plt.subplots(1, figsize=(figw, figh))
        pl = self.get_kde_plot(ax, **kwargs)
        ax.grid(b="on",  color="Grey")
        ax.set_xlabel("Pixel Value")
        ax.set_ylabel("Normalized Gaussian Kernel")
        ax.set_title("PDF Estimation (mask = {0})".format(kwargs["mask_id"]))
        
    def get_histo(self, bins_num=150, normed_sel=True, **mask_args):
        """
        Mask arguments: "masked" (bool), "mask" (array),"mask_value" (float), and "mask_id" (only if masked=True)
        :param bins_num:
        :param normed_sel:
        :param mask_args:
        :return:
        """
        #Histogram data density (pdf) estimation
        if normed_sel:
            hist_label = "Normed Histogram"
        elif not normed_sel:
            hist_label = "Histogram"
        histo, bins = np.histogram(self.get_flat_array(**mask_args), bins_num,  normed=normed_sel)
        bincntrs = 0.5*(bins[1:]+bins[:-1])
        return bincntrs, histo

    def plot_histo(self, **mask_args):
        figw = 12
        figh = 5
        fig, ax = plt.subplots(1, figsize=(figw, figh))
        bincntrs, histo = self.get_histo(**mask_args)
        h = ax.plot(bincntrs, histo, 'k', linewidth=1, label=self.label)
        ax.grid(b="on",  color ="Grey")
        ax.set_xlabel("Pixel Value")
        ax.set_ylabel(self.label)
        ax.set_title("PDF Estimation for {0}".format(self.label))

    def pdf_fits(self, bins, **mask_args):
        #Normal density fit
        mu, std = stats.norm.fit(self.get_nonmasked_values(**mask_args))  # So this is an actual fit of the data distribution
        normfit = stats.norm.pdf(bins, mu, std)
        return mu, std, normfit


