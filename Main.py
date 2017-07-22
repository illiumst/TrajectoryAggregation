# -*- coding: utf-8 -*-

import numpy as np
import utilities as ut
from math import cos
from osgeo import ogr, osr
import pytess, operator
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt


# gdal.SetConfigOption('OGR_ENABLE_PARTIAL_REPROJECTION', 'YES')
##############################################################################


# Boot Up Phase
print 'initialising geographical toolchain'
WGS84 = osr.SpatialReference()
WorldMercator = osr.SpatialReference()
WGS84.ImportFromEPSG(4326)
WorldMercator.ImportFromEPSG(3395)
trans_to_m = osr.CoordinateTransformation(WGS84, WorldMercator)
trans_to_degree = osr.CoordinateTransformation(WorldMercator, WGS84)
esriDriver = ogr.GetDriverByName('ESRI Shapefile')

print 'create an temporal datasource in memory'
memDriver = ogr.GetDriverByName('MEMORY')
mem_source = memDriver.CreateDataSource('memdrive')


def load_from_textfile(filename):
    memoryLayer = mem_source.CreateLayer('textfile', WGS84, ogr.wkbPoint)
    memoryLayer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('userID', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('lon', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('lat', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('info1', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('info2', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('dtime_str', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('status', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('info3', ogr.OFTString))
    memoryLayer.CreateField(ogr.FieldDefn('info4', ogr.OFTString))
    with open(filename) as f:
        for line in f:
            if f is not None or f != '':
                attrList = line.split(',')
                textFeature = ogr.Feature(memoryLayer.GetLayerDefn())
                geometry = ogr.Geometry(ogr.wkbPoint)
                geometry.AddPoint_2D(float(attrList[2]),float(attrList[3]))
                textFeature.SetGeometry(geometry)
                for i, value in enumerate(attrList):
                    textFeature.SetField(i, value)
                memoryLayer.CreateFeature(textFeature)
    return memoryLayer


def load_to_memory(ShapeName):
    print 'loading local files'
    ShapeName = ShapeName if ShapeName.endswith('.shp') else '%s%s' % (ShapeName, '.shp')
    shapefile = esriDriver.Open(ShapeName, 0)
    layer = mem_source.CreateLayer('layer', WGS84, ogr.wkbPoint)
    inputLayer = shapefile.GetLayer()
    OLDDefn = inputLayer.GetLayerDefn()
    print 'copying to memory'
    for i in range(OLDDefn.GetFieldCount()):
        layer.CreateField(OLDDefn.GetFieldDefn(i))
    for feature in inputLayer:
        layer.CreateFeature(feature)
    print 'finished loading - closing Input - cleanup'
    layer.ResetReading()
    return layer


def transform_to_array(layer):
    print 'transfering shapefile to numpy Array'
    array = np.array([[feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()] for feature in layer])
    layer.ResetReading()
    return array


def apply_dbscan(layer, epsFactor, minpts=0, eps=0, show=False, plot=None):
    def plot_results():
        # Plot result
        # Black removed and is used for noise instead.
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True

        unique_labels = set(labels)
        colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'

            class_member_mask = (labels == k)

            xy = array[class_member_mask & core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=14)

            xy = array[class_member_mask & ~core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=0.2)

        plt.title('Estimated number of clusters: %d' % n_clusters_)
        if show:
            plt.show()
        if plot is not None:
            plt.savefig(plot, bbox_inches='tight')

    if minpts <= 0:
        print 'calculate minPts'
        geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
        for feature in layer:
            geomcol.AddGeometry(feature.GetGeometryRef())

        layer.ResetReading()
        convexhull = geomcol.ConvexHull()
        convexhull.Transform(trans_to_m)
        minpts = (convexhull.GetArea() / layer.GetFeatureCount()) / 4000  # per 5km**2
        print 'automatic minPts at', minpts
    else:
        print 'determined minPts =', minpts

    if eps <= 0:
        print 'calculate eps'
        minX, maxX, minY, maxY = layer.GetExtent()
        eps = epsFactor * cos(minY + (maxY - minY)/2)  # Change the factor for a variation in search Radius
        print 'eps at', eps
    else:
        print 'determined eps =', eps

    # Compute DBSCAN
    array = transform_to_array(layer)
    print 'starting DBSCAN'
    db = DBSCAN(eps=eps, min_samples=minpts).fit(array)
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    if plot or show:
        plot_results()
    return labels


def build_cluster_layer(labels, array, shapename=None):
    print 'creating cluster layer'
    xValues, yValues = [], []
    # Create an new Layer with Cluster centroids
    layer = mem_source.CreateLayer('ClusterLayer', WGS84, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn('lon', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('lat', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('cluster', ogr.OFTInteger64))

    print 'sorting points to clusters'
    geomColDict = dict()
    for num, label in enumerate(labels):
        if label != -1:
            if geomColDict.get(label, None) is None:
                geomColDict[label] = ogr.Geometry(ogr.wkbGeometryCollection)
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(array[num, 0], array[num, 1])
            geomColDict[label].AddGeometry(point)

    print 'create clusters centroid features'
    for label in set(labels):
        if label != -1:
            feature = ogr.Feature(layer.GetLayerDefn())
            centre = geomColDict[label].Centroid()
            feature.SetGeometry(centre)
            feature.SetField('lon', centre.GetX())
            xValues.append(centre.GetX())
            feature.SetField('lat', centre.GetY())
            yValues.append(centre.GetY())
            feature.SetField('cluster', label)
            layer.SetFeature(feature)

    print 'Cluster Centroids Created'
    if shapename is not None:
        ut.CreateShapefile(shapename, ogr.wkbPoint, layer)
    return [layer, xValues, yValues]


def build_voronois(xValues, yValues, shapename=None):
    """
    Create an new Layer with Voronoi Polygons from clusters
    :return:
    """
    print 'building Voronoi Polygons from Cluster Centroids'
    layer = mem_source.CreateLayer('VoronoiLayer', WGS84, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn('lon', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('lat', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('cluster', ogr.OFTInteger64))
    voronois = pytess.voronoi(zip(xValues, yValues))

    print 'wirte features in Voronoi Layer'
    for label, (centre, coordinates) in enumerate(voronois):
        wkt = 'Polygon (('
        coordinates.append(coordinates[0])
        for x, y in coordinates:
            wkt = '%s %f %f,' % (wkt, x, y)
        wkt = '%s))' % wkt[:-1]

        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkt(wkt))
        feature.SetField('lon', centre[0] if centre is not None else 0)
        feature.SetField('lat', centre[1] if centre is not None else 0)
        feature.SetField('cluster', label)
        layer.CreateFeature(feature)

    if shapename is not None:
        ut.CreateShapefile(shapename, ogr.wkbPolygon, layer)
    print 'Voronois Done'
    return layer


def build_trajectorys(inputLayer, inputSource, voronoiLayer, shapename=None, initial=False, generalized=False, final=False):
    print 'Preparing User Dictionaries and setting Layers'
    layerDefn = inputLayer.GetLayerDefn()
    memInitialTrajectories = mem_source.CreateLayer('InitialTrajecLayer', WGS84, ogr.wkbPoint)
    for i in range(layerDefn.GetFieldCount()):
        memInitialTrajectories.CreateField(layerDefn.GetFieldDefn(i))
    memTrajectoryLayer = mem_source.CreateLayer('TrajectoryLayer', WGS84, ogr.wkbPoint)
    for i in range(layerDefn.GetFieldCount()):
        memTrajectoryLayer.CreateField(layerDefn.GetFieldDefn(i))

    trajectoryDict, featureDict = dict(), dict()

    print 'Creating All initial Trajectories'
    for i in range(0, inputLayer.GetFeatureCount()):
        feature = inputLayer.GetFeature(i)
        featureGeom = feature.GetGeometryRef()
        pointX = pointY = None
        for v in range(0, voronoiLayer.GetFeatureCount()):
            areaFeature = voronoiLayer.GetFeature(v)
            areaGeom = areaFeature.GetGeometryRef()
            if featureGeom.Intersects(areaGeom):
                pointX, pointY = areaFeature.GetField('lon'), areaFeature.GetField('lat')
                # featureGeom = areaGeom.Centroid()      # If u want to use the geometric center of a Vornonoi Shape
                break
        if pointX and pointY is not None:
            userID = feature.GetFieldAsDouble('userID')
            if trajectoryDict.get(userID, None) is None:
                trajectoryDict[userID] = ogr.Geometry(ogr.wkbLineString)
            trajectoryDict[userID].AddPoint(pointX, pointY)
            if featureDict.get(userID, None) is None:
                featureDict[userID] = feature

    fieldID = 0
    for trajectory, userID in zip(trajectoryDict.values(), trajectoryDict):
        if trajectory.GetPointCount() >= 2 and trajectory.Length() >= 0.0:
            trajectoryFeature = featureDict[userID]
            trajectoryFeature.SetGeometry(trajectory)
            trajectoryFeature.SetFID(fieldID)
            fieldID += 1
            memInitialTrajectories.CreateFeature(trajectoryFeature)

    if generalized or final:
        percentageBorder, initalCount, similarityDict, FID, deletedFeature = \
            0, memInitialTrajectories.GetFeatureCount(), dict(), 0, 0
        print 'mapping Trajectories to Cluster centroids'
        # This is just for visualization of progress
        for i in range(0, initalCount):
            trajectoryFeature = memInitialTrajectories.GetFeature(i)
            if i >= percentageBorder:
                print 'processing %i of %i' % (i, initalCount)
                print str(i / float(initalCount) * 100)[:2], '%'
                percentageBorder += (initalCount/20.0)

            newTrajectoryGeom = ogr.Geometry(ogr.wkbLineString)
            trajectoryGeom = trajectoryFeature.GetGeometryRef()

            pointList = [trajectoryGeom.GetPoint(x) for x in range(0, trajectoryGeom.GetPointCount())]
            for ((x1, y1, z1), (x2, y2, z2)) in zip(pointList[:-1], pointList[1:]):
                wkt = 'LINESTRING (%f %f, %f %f)' % (x1, y1, x2, y2)
                segmentGeom = ogr.CreateGeometryFromWkt(wkt)
                pointAndDistList = []
                start = ogr.Geometry(ogr.wkbPoint)
                start.AddPoint(x1, y1)

                for v in range(0, voronoiLayer.GetFeatureCount()):
                    areaFeature = voronoiLayer.GetFeature(v)
                    areaGeom = areaFeature.GetGeometryRef()
                    if segmentGeom.Intersects(areaGeom):
                        point = ogr.Geometry(ogr.wkbPoint)
                        point.AddPoint(areaFeature.GetField('lon'), areaFeature.GetField('lat'))
                        segmentCentre = point
                        # segmentCentre = areaGeom.Centroid()
                        pointAndDistList.append((segmentCentre, start.Distance(segmentCentre),
                                                 areaFeature.GetField('lon'), areaFeature.GetField('lat')))

                pointAndDistList.sort(key=operator.itemgetter(1))

                for point, distance, lon, lat in pointAndDistList:
                    newTrajectoryGeom.AddPoint(lon, lat)

            if newTrajectoryGeom.Length() == 0.0:
                memInitialTrajectories.DeleteFeature(i)
                deletedFeature += 1
                continue  # TODO: Find a better solution for just deleting the specific Feature, e.g. Count this occurence

            else:
                pointList = [newTrajectoryGeom.GetPoint(i) for i in range(0, newTrajectoryGeom.GetPointCount())]
                for ((x1, y1, z1), (x2, y2, z2)) in zip(pointList[:-1], pointList[1:]):
                    if x1 != x2 and y1 != y2:
                        wkt = 'LINESTRING (%s %s, %s %s)' % (x1, y1, x2, y2)
                        splitLine = ogr.CreateGeometryFromWkt(wkt)
                        trajectoryFeature.SetGeometry(splitLine)
                        trajectoryFeature.SetFID(FID)
                        FID += 1
                        memTrajectoryLayer.CreateFeature(trajectoryFeature)
                        if similarityDict.get(wkt, None) is None:
                            similarityDict[wkt] = [trajectoryFeature.Clone()]
                        else:
                            similarityDict[wkt].append(trajectoryFeature.Clone())
    print deletedFeature, 'Features were deleted due to inappropiated Lengths'
    if final:
        userfollow, userfriend, userstatus = 0, 0, 0
        for i in range(0, memTrajectoryLayer.GetFeatureCount()):
            feature = memTrajectoryLayer.GetFeature(i)
            userfollow += feature.GetFieldAsDouble('userfollow')
            userfriend += feature.GetFieldAsDouble('userfriend')
            userstatus += feature.GetFieldAsDouble('userstatus')
            pass

        userfollow /= memTrajectoryLayer.GetFeatureCount()
        userfriend /= memTrajectoryLayer.GetFeatureCount()
        userstatus /= memTrajectoryLayer.GetFeatureCount()
        finalTrajectories = mem_source.CreateLayer('FinalTrajectoriLayer', WGS84, ogr.wkbPoint)
        finalTrajectories.CreateField(ogr.FieldDefn('MaxFollow', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MaxFriend', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MaxStatus', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MinFollow', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MinFriend', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MinStatus', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MeanFollow', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MeanFriend', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('MeanStatus', ogr.OFTReal))
        finalTrajectories.CreateField(ogr.FieldDefn('Group', ogr.OFTInteger64))
        finalTrajectories.CreateField(ogr.FieldDefn('Size', ogr.OFTInteger64))
        groupDict = dict()
        for num, key in enumerate([(True, True, True), (True, True, False), (True, False, True), (True, False, False),
                                   (False, True, True), (False, True, False), (False, False, True), (False, False, False)],
                                  start=1):
            groupDict[key] = num

        for featureList in similarityDict.values():
            summary = np.array(zip(*[
                (feat.GetFieldAsDouble('userfollow'),
                 feat.GetFieldAsDouble('userfriend'),
                 feat.GetFieldAsDouble('userstatus')) for feat in featureList]))
            Max, Min, Mean = [max(summary[0, ]), max(summary[1, ]), max(summary[2, ])], \
                             [min(summary[0, ]), min(summary[1, ]), min(summary[2, ])],\
                             [sum(summary[0, ]) / float(len(summary[0, ])),
                              sum(summary[1, ]) / float(len(summary[1, ])),
                              sum(summary[2, ]) / float(len(summary[2, ]))]
            key = (Mean[0] >= userfollow, Mean[1] >= userfriend, Mean[2] >= userstatus)

            group = groupDict[key]
            allLists = []
            [allLists.extend(el) for el in [Max, Min, Mean, [group], [len(featureList)]]]

            finalFeature = ogr.Feature(finalTrajectories.GetLayerDefn())
            finalFeature.SetGeometry(featureList[0].GetGeometryRef())
            for i, value in enumerate(allLists):
                finalFeature.SetField(i, value)
            finalTrajectories.CreateFeature(finalFeature)

    if shapename is not None:
        if initial:
            ut.CreateShapefile('%s_initial' % shapename, ogr.wkbLineString, memInitialTrajectories)
        if generalized:
            ut.CreateShapefile('%s_general' % shapename, ogr.wkbLineString, memTrajectoryLayer)
        if final:
            ut.CreateShapefile('%s_final' % shapename, ogr.wkbLineString, finalTrajectories)
    print len(similarityDict), 'was the length of the dict'


if __name__ == '__main__':
    shanghai_all = load_from_textfile('sjtu_taxigps20070201.txt')
    #weibo_all = load_to_memory('weibo_all.shp')
    weibo_array = transform_to_array(shanghai_all)
    cluster_list = apply_dbscan(shanghai_all, 0.003, show=False)
    #cluster_list = apply_dbscan(weibo_all, 0.010, minpts=20, plot='smaller_cluster_minpts_20_%f.png' % 0.010)
    cluster_layer, xs, ys = build_cluster_layer(cluster_list, weibo_array)
    voronoi_layer = build_voronois(xs, ys)

    # Test Algorithm Mode
    #shapefileName = 'weibo_test.shp'
    #shapefile = esriDriver.Open(shapefileName, 0)
    #testlayer = shapefile.GetLayer()
    #build_trajectorys(testlayer, shapefile, voronoi_layer,
    #                  shapename='test', initial=False, generalized=False, final=True)

    #Full Mode
    build_trajectorys(shanghai_all, mem_source, voronoi_layer, 'trajec', generalized=True, final=False)


