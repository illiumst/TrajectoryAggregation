
from osgeo import ogr, osr
import os
import pyproj
import numpy as np


# Geographic Services
esri_driver = ogr.GetDriverByName("ESRI Shapefile")
osmdriver = ogr.GetDriverByName("OSM")
WGS84 = osr.SpatialReference()
WorldMercator = osr.SpatialReference()
WGS84.ImportFromEPSG(4326)
WorldMercator.ImportFromEPSG(3395)
trans_to_m = osr.CoordinateTransformation(WGS84, WorldMercator)
trans_to_degree = osr.CoordinateTransformation(WorldMercator, WGS84)


def showLayerInfo(Layer):
    layerDefinition = Layer.GetLayerDefn()

    print "Name  -  Type  Width  Precision"
    for i in range(layerDefinition.GetFieldCount()):
        fieldName = layerDefinition.GetFieldDefn(i).GetName()
        fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
        fieldType = layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
        fieldWidth = layerDefinition.GetFieldDefn(i).GetWidth()
        GetPrecision = layerDefinition.GetFieldDefn(i).GetPrecision()

        print fieldName + " - " + fieldType + " " + str(fieldWidth) + " " + str(GetPrecision)


def CreateShapefile(ShapeName, FeatureType, Geometry_OR_Layer):
    """

    Parameters
    ----------
    ShapeName: basestring as Name for the new Shapefile
    FeatureType: osgeo.ogr.wkbtype e.g. ogr.wkbPoint, wkbLineString,
                 wkbPolygon, wkbMultiPoint, wkbMultiLineString, wkbMultiPolygon,
                 wkbGeometryCollection
    Geometry_OR_Layer: osgeo.ogr.Geometry / osgeo.ogr.Layer class object

    Returns
    -------
    None but creates an ShapeName.shp Shapefile in the current directory.

    """
    ShapeName = ShapeName if ShapeName.endswith(".shp") else "".join([ShapeName, ".shp"])
    if os.path.exists(ShapeName):
        esri_driver.DeleteDataSource(ShapeName)
    datasource = esri_driver.CreateDataSource(ShapeName)
    layer = datasource.CreateLayer("layer", WGS84, FeatureType, options=['ENCODING=UTF-8'])
    if type(Geometry_OR_Layer) == ogr.Geometry:
        if Geometry_OR_Layer.GetGeometryName() == "GEOMETRYCOLLECTION":
            for idx in range(Geometry_OR_Layer.GetGeometryCount()):
                geometry = Geometry_OR_Layer.GetGeometryRef(idx)
                feature = ogr.Feature(layer.GetLayerDefn())
                feature.SetGeometry(geometry)
                layer.CreateFeature(feature)

        else:
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(Geometry_OR_Layer)
            layer.CreateFeature(feature)

    elif type(Geometry_OR_Layer) == ogr.Layer:
        Geometry_OR_Layer.ResetReading()
        OLDDefn = Geometry_OR_Layer.GetLayerDefn()
        for i in range(OLDDefn.GetFieldCount()):
            layer.CreateField(OLDDefn.GetFieldDefn(i))
        featurecount = 0
        for feature in Geometry_OR_Layer:
            featurecount += 1
            if type(feature) == ogr.Feature:
                layer.CreateFeature(feature)
            else:
                print type(feature)
                continue
        Geometry_OR_Layer.ResetReading()
        print featurecount, "features written!"
    # elif type(Geometry_OR_Layer) ==  ogr.Geometr
    else:
        raise TypeError("This is not a geometry or a layer.")


def near_search_n(QueryPoint, PointList, n):
    """
    find N nearest neighbours of QueryPoint among PointList

    Parameters
    ----------
    QueryPoint: single row 2d-nparray with x / y coordinates
    PointList: multi row 2d-nparray with x / y coordinates
    n: int number of points to be returned

    Returns
    -------
    list of indices in order of ascending distance

    """

    ndata = PointList.shape[1]
    n = n if n < ndata else ndata
    # euclidean distances from the other points
    sqd = np.sqrt(((PointList - QueryPoint[:, :ndata]) ** 2).sum(axis=0))
    idx = np.argsort(sqd)  # sorting
    # return the indexes of K nearest neighbours
    return idx[:n].tolist()


def distance((x1, y1), (x2, y2)):
    geod = pyproj.Geod(ellps="WGS84")
    heading1, heading2, measured_distance = geod.inv(x1, y1, x2, y2)
    return measured_distance

