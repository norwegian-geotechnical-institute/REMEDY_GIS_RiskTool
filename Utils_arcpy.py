import arcpy
import os
import sys
import json
import importlib
import BegrensSkadeLib
importlib.reload(BegrensSkadeLib)

'''
The library contains common utility functions for the BegrensSkade GIBV ArcGIS wrappers.
All functions are dependent on arcpy and needs an ArcGIS lisence.
'''


def extractBuildingsFromFL(buildings, features, buildingsOut, logger):
    logger.debug("TIME - Deleting old buildings: {}".format(buildingsOut))
    arcpy.Delete_management(buildingsOut)
    arcpy.SelectLayerByLocation_management(buildings, "INTERSECT", features, "", "NEW_SELECTION")
    logger.debug("Copy features")
    arcpy.CopyFeatures_management(buildings, buildingsOut)
    logger.debug("Builings extracted")

def addLayer(pMap, dataPath, lyrPath, name=""):
    try:
        newLyr = arcpy.mp.LayerFile(lyrPath)
        newLay = pMap.addLayer(newLyr)[0]
        conProp = newLay.connectionProperties
        dataFolder = os.path.dirname(dataPath)
        filename = os.path.basename(dataPath)
        conProp["connection_info"]["database"] = dataFolder
        conProp["dataset"] = filename
        newLay.updateConnectionProperties(newLay.connectionProperties, conProp, False, True)
        if name != "":
            newLay.name = name
        return newLay
    except:
        arcpy.AddError("Error in addLayer()")
        sys.exit()
        raise


def addLayerToGroup(pMap, dataPath, lyrPath, targetGroupLayer, name=""):
    try:
        newLyr = arcpy.mp.LayerFile(lyrPath)
        newLay = pMap.addLayer(newLyr)[0]
        conProp = newLay.connectionProperties
        dataFolder = os.path.dirname(dataPath)
        filename = os.path.basename(dataPath)
        conProp['connection_info']['database'] = dataFolder
        conProp['dataset'] = filename
        newLay.updateConnectionProperties(newLay.connectionProperties, conProp, False, True)
        if name != "":
            newLay.name = name
        pMap.addLayerToGroup(targetGroupLayer, newLay)
        pMap.removeLayer(newLay)
    except:
        arcpy.AddError("Error in addLayerToGroup()")
        sys.exit()
        raise

def getProjCodeFromFC(fc):
    return arcpy.Describe(fc).spatialReference.PCSCode

def getConstructionAsJson(uploaded_construction_outline):
    if arcpy.ProductInfo() == 'ArcServer':
        bounding_as_json = uploaded_construction_outline.JSON
        return json.loads(bounding_as_json)
    else:
        with arcpy.da.SearchCursor(uploaded_construction_outline, ['SHAPE@']) as cursor:
            for row in cursor:
                geom = row[0]
                bounding_as_json = geom.JSON
                return json.loads(bounding_as_json)

def getBuildingsClipExtentFromConstruction(constructionJson, margin, working_proj, logger):
    logger.debug("constructionJson: {}".format(constructionJson))
    construction_area_corners = BegrensSkadeLib.get_construction_corners_from_ArcGIS_json(constructionJson, 2, logger)
    sr = arcpy.SpatialReference(working_proj)
    minx = 99999999
    miny = 99999999
    maxx = -99999999
    maxy = -99999999

    for corner in construction_area_corners:
        if corner.x > maxx:
            maxx = corner.x
        if corner.x < minx:
            minx = corner.x
        if corner.y > maxy:
            maxy = corner.y
        if corner.y < miny:
            miny = corner.y

    ll = arcpy.Point(minx - margin, miny - margin)
    ul = arcpy.Point(minx - margin, maxy + margin)
    ur = arcpy.Point(maxx + margin, maxy + margin)
    lr = arcpy.Point(maxx + margin, miny - margin)

    rectangle = arcpy.Polygon(arcpy.Array([ll, ul, ur, lr, ll]), sr)

    logger.debug("Rectangle: {}".format(rectangle.extent.JSON))
    return rectangle

def CheckOutLicense(ext, extName):
    if arcpy.CheckExtension(ext) == "Available":
        arcpy.AddMessage("Checking out a " + extName + " Extension...")
        arcpy.CheckOutExtension(ext)
        return True
    else:
        return False