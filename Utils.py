import sys
import os
import numpy as np
import json
import math
import glob
from typing import Dict, Any

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import zipfile


try:
    from gdalconst import GA_ReadOnly
except:
    from osgeo.gdalconst import GA_ReadOnly

"""
Biblioteket inneholder basisfunksjonalitet for håndtering av gis-data med osgeo/gdal og json/requests. 
Ikke avhengig av arcpy eller ArcGIS lisens.
"""

epsgMapping = {
    "ETRS_1989_UTM_Zone_32N": "25832",
    "ETRS_1989_UTM_Zone_33N": "25833",
    "ETRS_1989_NTM_Zone_5": "5105",
    "ETRS_1989_NTM_Zone_6": "5106",
    "ETRS_1989_NTM_Zone_7": "5107",
    "ETRS_1989_NTM_Zone_8": "5108",
    "ETRS_1989_NTM_Zone_9": "5109",
    "ETRS_1989_NTM_Zone_10": "5110",
    "ETRS_1989_NTM_Zone_11": "5111",
    "ETRS_1989_NTM_Zone_12": "5112",
    "ETRS_1989_NTM_Zone_13": "5113",
    "ETRS_1989_NTM_Zone_14": "5114",
    "ETRS_1989_NTM_Zone_15": "5115",
    "ETRS_1989_NTM_Zone_16": "5116",
    "ETRS_1989_NTM_Zone_17": "5117",
    "ETRS_1989_NTM_Zone_18": "5118",
    "ETRS_1989_NTM_Zone_19": "5119",
    "ETRS_1989_NTM_Zone_20": "5120",
}


"""
Beskrivelse:
    Konverterer boolske variable lagret som tekst til python-boolske variable
Forutsetninger:
    .
Parametre:
    boolStr:     Boolsk variabel lagret som tekst
Return:
    Python-boolsk variabel
"""


def setBooleanParameter(boolStr):
    if str(boolStr).lower() == "true":
        return True
    else:
        return False


"""
Beskrivelse:
    Generelle vektorfunksjonalitet
Forutsetninger:
    .
Parametre:
    x1,y1 og x2,y2:     Vektorer
    xref, yref:         Referansevektor
Return:
    Vinkel mellom vektorer i riktig kvadrant

"""


def get_angle(x1, y1, x2, y2):

    if x1 == x2 and y1 == y2:
        return "Fail"

    if x2 == x1:
        if y2 > y1:
            return np.pi / 2, 4
        elif y2 < y1:
            return -3 * np.pi / 2, 2

    if y2 == y1:
        if x2 > x1:
            return 0, 1
        elif x2 < x1:
            return np.pi, 3

    angle = np.arctan((y2 - y1) / (x2 - x1))
    quadrant = 1

    if x2 < x1 and y2 > y1:
        angle += np.pi
        quadrant = 2
    elif x2 < x1 and y2 < y1:
        angle += np.pi
        quadrant = 3
    elif x2 > x1 and y2 < y1:
        angle += 2 * np.pi
        quadrant = 4

    return angle, quadrant


def getAngleFromDir(x, y, xref=0, yref=0):
    angle_quadr = abs(np.arctan((y - yref) / (x - xref)))
    if x < xref and y > yref:
        angle = np.pi - angle_quadr
    elif x < xref and y < yref:
        angle = np.pi + angle_quadr
    elif x > xref and y < yref:
        angle = 2 * np.pi - angle_quadr
    else:
        angle = angle_quadr

    return 180 * angle / np.pi


def getEPSGFromShape(sfn):
    # sf = gpd.read_file(sfn)
    # epsgStr = sf.crs.to_epsg()

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(sfn, 0)
    layer = dataSource.GetLayer()

    spatialRef = layer.GetSpatialRef()

    projcs = spatialRef.GetAttrValue("PROJCS", 0)
    print("projcs:{} ".format(projcs))

    epsg = epsgMapping[projcs]

    # epsg = layer.GetSpatialRef().GetAttrValue('AUTHORITY', 1)

    print("epsg:{} ".format(epsg))
    return epsg


def getEPSGFromRaster(rasterFile):
    ds = gdal.Open(rasterFile)
    prj = ds.GetProjection()
    print("prj:{} ".format(prj))
    # print(prj)

    srs = osr.SpatialReference(wkt=prj)
    if srs.IsProjected:
        projcs = srs.GetAttrValue("projcs")
        print("projcs:{} ".format(projcs))
        epsg = epsgMapping[projcs]
        return epsg
    else:
        # print srs.GetAttrValue('geogcs')
        raise Exception(
            "Not able to get projection info.  Supported projections EUREF89 UTM32 and 33, and NTM"
        )


"""
Beskrivelse:
    Zipper shapefil til zipfil med samme navn.
Forutsetninger:
    Zipper bare shp, shx, dbf og prj.
Parametre:
    shapeFileFolder:    Katalog der shapefil er.
    zipFileFolder:      Katalog der zipfilen havner.
    baseName:           navn på shapefil (uten .shp) og navnet på zipfilen (uten.zip)
Return:

"""


def zipShapefile(shapeFileFolder, zipFileFolder, baseName):

    # shapeFNFull = basePath + "\\" + shapeFileBaseName + ".shp"
    zipFileFN = zipFileFolder + "\\" + baseName + ".zip"

    zipObj = zipfile.ZipFile(zipFileFN, "w")

    # Add multiple files to the zip
    zipObj.write(shapeFileFolder + "\\" + baseName + ".shp", baseName + ".shp")
    zipObj.write(shapeFileFolder + "\\" + baseName + ".dbf", baseName + ".dbf")
    zipObj.write(shapeFileFolder + "\\" + baseName + ".prj", baseName + ".prj")
    zipObj.write(shapeFileFolder + "\\" + baseName + ".shx", baseName + ".shx")

    # close the Zip File
    zipObj.close()
    return zipFileFN


"""
Beskrivelse:
    Lager polygoner basert på koordinater i ring.
Forutsetninger:
    Stotter ikke hull etc.
Parametre:
    coords:    liste med liste med koordinatsett
"""


def createPolygon(coords):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[0], coord[1])

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


"""
Beskrivelse:
    Lager linje.
Forutsetninger:
    .
Parametre:
    coords:    liste med liste med koordinatsett
"""


def createLine(coords, logger=None):
    # if not logger is None:
    #     logger.debug("createLine, logger sent, coords: {}".format(
    #         coords))
    line = ogr.Geometry(type=ogr.wkbLineString)
    for coord in coords:
        line.AddPoint_2D(coord[0], coord[1])
    return line


"""
Beskrivelse:
    Lager punkter.
Forutsetninger:
    .
Parametre:
    x:    x koordinat for punkt
    y:    y koordinat for punkt
"""


def createPoint(x, y):
    wkt = "POINT(%f %f)" % (float(x), float(y))
    point = ogr.CreateGeometryFromWkt(wkt)
    return point


"""
Beskrivelse:
    Lager en feiltdefinisjon.
Forutsetninger:
    .
Parametre:
    fieldName:  Feltnavn
    fieldType:  Felttype
    fieldWidth: lengde på felt, brukes pr idag kun på stringfelt
    fieldPrec:  presisjon på felt, brukes pr idag ikke
Return:

"""


def createOGRFieldDef(filedName, fieldType, fieldWidth=100, fieldPrec=0):
    ft = None
    if fieldType.lower() == "string" or fieldType.lower() == "text":
        ft = ogr.OFTString
    elif fieldType.lower() == "float":
        ft = ogr.OFTReal
    elif fieldType.lower() == "double":
        ft = ogr.OFTReal
    elif fieldType.lower() == "int":
        ft = ogr.OFTInteger
    elif fieldType.lower() == "long":
        ft = ogr.OFTInteger64
    elif fieldType.lower() == "date":
        ft = ogr.OFTDate
    else:
        raise Exception(
            "Field definition %s not of legal type.  Available types: string, float, double, int, long, date"
            % fieldType
        )

    field = ogr.FieldDefn(filedName, ft)
    if ft == ogr.OFTString:
        field.SetWidth(fieldWidth)
    # elif ft == ogr.OFTReal:
    #     field.SetWidth(fieldWidth)
    #     field.SetPrecision(fieldPrec)
    return field


"""
Beskrivelse:
    Lager en shapefil med angitt felt, geometritype.  Støtter point og polygon.
Forutsetninger:
    .
Parametre:
    fileName:    Sti til shapefil som skal lages.
    geomType:    Geometri type (point eller polygon)
    projection:  Geometritype som tall ESPG nummer.  25832 UTM32, 25833 UTM33
    fieldNames:  Liste med feltnavn
    fieldTypes:  Liste med felttyper
Return:

"""


def createShapefile(fileName, geomType, projection, fields):
    if os.path.isfile(fileName):
        raise Exception("Shapefile %s already exist." % fileName)
    # set up the shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")

    # create the data source
    data_source = driver.CreateDataSource(fileName)

    # create the spatial reference
    srs = osr.SpatialReference()
    # srs.ImportFromEPSG(projection)
    res = srs.ImportFromEPSG(projection)  # projection is a the code as int
    if res != 0:
        raise Exception(repr(res) + ": could not import Spatial reference from EPSG")
        sys.exit()

    # create the layer
    ln = fileName.split(".")[0]
    layer = None
    if geomType.lower() == "point" or geomType.lower() == "points":
        layer = data_source.CreateLayer(ln, srs, ogr.wkbPoint)
    elif (
        geomType.lower() == "line"
        or geomType.lower() == "lines"
        or geomType.lower() == "polyline"
    ):
        layer = data_source.CreateLayer(ln, srs, ogr.wkbLineString)
    elif geomType.lower() == "polygon" or geomType.lower() == "polygons":
        layer = data_source.CreateLayer(ln, srs, ogr.wkbPolygon)

    # Add the fields we're interested in
    countFields = len(fields)
    for field in fields:
        ftn = field[1]
        fn = field[0]
        field = createOGRFieldDef(fn, ftn)
        layer.CreateField(field)
    data_source = None


def fieldExistsInShapefile(fileName, fieldName):
    if not os.path.isfile(fileName):
        raise Exception("Shapefile %s does not exist." % fileName)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    source = driver.Open(fileName)
    layer = source.GetLayer()
    ldefn = layer.GetLayerDefn()
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        print(fdefn.name)
        if fieldName.lower() == fdefn.name.lower():
            return True
    return False


"""
Beskrivelse:
    Legger til nytt felt for eksisterende shapefil.
Forutsetninger:
    Shapefilen må finnes.
Parametre:
    fileName:    Sti til shapefil som skal lages.
    fieldType:  Felttype
    fieldWidth: lengde på felt, brukes pr idag kun på stringfelt
    fieldPrec:  presisjon på felt, brukes pr idag ikke
Return:
"""


def addFieldToShapefile(fileName, fieldName, fieldType, fieldWidth=100, fieldPrec=0):
    if not os.path.isfile(fileName):
        raise Exception("Shapefile %s does not exist." % fileName)
    if fieldExistsInShapefile(fileName, fieldName):
        raise Exception("Shapefile allready has fieldname %s ." % fieldName)
    infile = ogr.Open(fileName, 1)
    inlyr = infile.GetLayerByIndex(0)
    fieldDef = createOGRFieldDef(fieldName, fieldType, fieldWidth, fieldPrec)
    inlyr.CreateField(fieldDef)


def getFieldDefinition(layer, fieldName, logger=None):
    layerDefinition = layer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        # logger.debug("FieldName: {}, fieldDefName: {}".format(
        #     fieldName, layerDefinition.GetFieldDefn(i).GetName()))
        if fieldName.lower() == layerDefinition.GetFieldDefn(i).GetName().lower():
            return layerDefinition.GetFieldDefn(i)
    return None


def addValueToField(feature, fieldNameTmp, value, layer, logger=None):
    fieldName = fieldNameTmp[:10]
    # ft = ""
    # if (not logger is None):
    # logger.debug("Logger sent, field name: {}, value: {}".format(
    #     fieldName, value))
    # if (not layer is None):
    #     if (logger.level == logging.DEBUG):
    fieldDef = getFieldDefinition(layer, fieldName, logger)
    ft = fieldDef.GetFieldTypeName(fieldDef.GetType())
    # logger.debug("Field name: {}, value: {}, fieldDefinition: {}".format(
    #     fieldName, value, ft))
    # else:
    #     logger.debug("Field name: {}, value: {}".format(fieldName, value))
    if not value is None:
        try:
            if ft == "Real" or ft == "Integer64" or ft == "Integer":
                x = float(value)
                isNumber = not math.isnan(x)
                if isNumber:
                    if ft == "Real":
                        feature.SetField(fieldName, float(value))
                    else:
                        feature.SetField(fieldName, int(value))
            else:
                feature.SetField(fieldName, value)
        except Exception as e:
            if not logger is None:
                logger.error(type(e))
                logger.error(e.args)
                logger.error(e)
                logger.debug(
                    "addValueToField: Setting field value: field name: {}, value: {}, fieldDefinition: {}".format(
                        fieldName, value, ft
                    )
                )


def writeOneFeatureToShapefile(outSHPfn, geom, fieldNames, fieldValues):
    # Open the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    outDataSource = shpDriver.Open(outSHPfn, 1)
    outLayer = outDataSource.GetLayer()

    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(geom)

    countFields = len(fieldNames)

    for i in range(countFields):
        outFeature.SetField(fieldNames[i], fieldValues[i])
    outLayer.CreateFeature(outFeature)

    outFeature = None
    # Save and close the data source
    outDataSource = None


"""
Beskrivelse:
    Skriver til en shapefil fra array med egenskaper og geometrier.
Forutsetninger:
    geometrien i geometriListen er laget med utilityfunksjonen greatePoint eller createPolygon.
    Antall elementer i geomList og valuesList må være lik.  fieldNames inneholder bare en rad.
Parametre:
    outSHPfn:    Sti til shapefil som skal skrives.
    fieldNames:  liste med feltnavn
    geomList:    Liste med geometrier laget i createPoint eller createPolygon
    valuesList:  Liste med verdier for hver feature.
Return:

"""


def writeToShapefileFromArray(outSHPfn, fieldNames, geomList, valuesList, logger=None):
    # Open the output shapefile
    if not logger is None:
        logger.debug("Logger sent, Shapefile: {}".format(outSHPfn))
        logger.debug("Logger sent, fieldNames: {}".format(fieldNames))
        logger.debug("Logger sent, geomList: {}".format(geomList))
        logger.debug("Logger sent, valuesList: {}".format(valuesList))
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    outDataSource = shpDriver.Open(outSHPfn, 1)
    outLayer = outDataSource.GetLayer()

    countFields = len(fieldNames)
    listLen = len(geomList)
    featureDefn = outLayer.GetLayerDefn()

    for fNum in range(listLen):
        geom = geomList[fNum]
        values = valuesList[fNum]
        # Create the feature and set values
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geom)

        for i in range(countFields):
            outFeature.SetField(fieldNames[i], values[i])
        outLayer.CreateFeature(outFeature)

    outFeature = None
    # Save and close the data source
    outDataSource = None


"""
Beskrivelse:
    Konverterer features i shapefil til json ihht ArcGIS format.
Forutsetninger:

Parametre:
    shapeFN:     Sti til shapefil som skal konverteres.
    geomType:    Geometritype for shapefil (point eller polygon)
    token:          Token hentet fra AGOL basert på client_id og client_secret
Return:
    Liste med json objekter med features.  En feature pr element
"""


def getShapefileAsJson(shapeFN, logger=None):
    if not logger is None:
        logger.debug("Logger sent, shapeFN name: {}".format(shapeFN))

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFN, 0)
    layer = dataSource.GetLayer()
    lyrDefn = layer.GetLayerDefn()
    geomType = lyrDefn.GetGeomType()
    # ogr.wkbLineString

    epsgNum = 25833

    jsonList = []
    layerDefinition = layer.GetLayerDefn()

    fieldNames = []

    for i in range(layerDefinition.GetFieldCount()):
        fieldNames.append(layerDefinition.GetFieldDefn(i).GetName())

    for feature in layer:
        gjson = feature.ExportToJson()
        pj = json.loads(gjson)

        x = ""
        if geomType == ogr.wkbPolygon:
            x = [
                {
                    "geometry": {
                        "rings": pj["geometry"]["coordinates"],
                        "spatialReference": {"wkid": epsgNum},
                    },
                    "attributes": {fn: pj["properties"][fn] for fn in fieldNames},
                }
            ]
        elif geomType == ogr.wkbPoint:
            x = [
                {
                    "geometry": {
                        "x": pj["geometry"]["coordinates"][0],
                        "y": pj["geometry"]["coordinates"][1],
                        "spatialReference": {"wkid": epsgNum},
                    },
                    "attributes": {fn: pj["properties"][fn] for fn in fieldNames},
                }
            ]

        xjs = json.dumps(x)
        jsonList.append(xjs)
    return jsonList


def deleteShapefilAndZip(folder, basename, logger=None):

    # get a recursive list of file paths that matches pattern including sub directories
    fileList = glob.glob(folder + "/" + basename + ".*", recursive=False)

    # Iterate over the list of filepaths & remove each file.
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            if not logger is None:
                logger.debug("Error while deleting file: {}".format(filePath))


def getRasterExtent(rasterFile):
    data = gdal.Open(rasterFile, GA_ReadOnly)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    return minx, miny, maxx, maxy


def getRasterMinMax(rasterFile):
    # open raster and choose band to find min, max
    # raster = r'C:\path\to\your\geotiff.tif'
    gtif = gdal.Open(rasterFile)
    srcband = gtif.GetRasterBand(1)

    # Get raster statistics
    stats = srcband.GetStatistics(True, True)
    return stats[0], stats[1]


def appendZValuesFromRaster(shapefile, raster, fieldName):
    # shp_filename = shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")

    src_ds = gdal.Open(raster)
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)

    ds = driver.Open(shapefile, 1)
    lyr = ds.GetLayer()
    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx, my = geom.GetX(), geom.GetY()  # coord in map units

        # Convert from map to pixel coordinates.
        # Only works for geotransforms with no rotation.
        px = int((mx - gt[0]) / gt[1])  # x pixel
        py = int((my - gt[3]) / gt[5])  # y pixel

        intval = rb.ReadAsArray(px, py, 1, 1)
        feat.SetField(fieldName, str(intval[0][0]))
        lyr.SetFeature(feat)
    ds.Destroy()


def updateFieldValues(layer, fieldIndex, value):
    for feat in layer:
        feat.SetField(fieldIndex, value)
        layer.SetFeature(feat)


"""
Beskrivelse:
    Projiserer raster til angitt projeksjon. Resampling nearest neighbour.
    Det skrives ikke world fil, projeksjonsinformasjon legges i rasterert
Forutsetninger:

Parametre:
    fromRaster: raster som blir projesert.  Må ha definert projeksjon.
    toRaster:   raster som blir skrevet.
    toESPG:     projeksjonsdefinisjon.  For eksempel "25832"

Div info:
   https://gdal.org/python/osgeo.gdal-module.html#WarpOptions


    Flere options på denne måten (må implementeres i funksjon):
    opt = gdal.WarpOptions(srcSRS='EPSG:25833',dstSRS='EPSG:25832')
    ds = gdal.Warp(rfnu, rfn, options=opt)
    del ds

    resmapling mode kan være {nearest (default),bilinear,cubic,cubicspline,lanczos,average,mode}
    for eksempel
    resampleAlg="bilinear"
"""


def projectRaster(fromRaster, toRaster, toESPG):
    if not os.path.isfile(fromRaster):
        raise Exception(
            "Sourcefile %s does not exist, not able to project raster." % fromRaster
        )
    if os.path.isfile(toRaster):
        raise Exception(
            "Toraster %s allready exist, not able to project raster" % fromRaster
        )
    toESPGStr = "EPSG:" + str(toESPG)
    opt = gdal.WarpOptions(dstSRS=toESPGStr)
    ds = gdal.Warp(toRaster, fromRaster, options=opt)
    del ds


def getRasterExentV1(rasterFile):
    ras = gdal.Open("/foo/bar/filename")
    upx, xres, xskew, upy, yskew, yres = ras.GetGeoTransform()
    cols = ras.RasterXSize
    rows = ras.RasterYSize

    ulx = upx + 0 * xres + 0 * xskew
    uly = upy + 0 * yskew + 0 * yres

    llx = upx + 0 * xres + rows * xskew
    lly = upy + 0 * yskew + rows * yres

    lrx = upx + cols * xres + rows * xskew
    lry = upy + cols * yskew + rows * yres

    urx = upx + cols * xres + 0 * xskew
    ury = upy + cols * yskew + 0 * yres


def getProjections():
    projections: Dict[Any, Any] = {}

    script_dir = os.path.dirname(__file__)
    file_path = os.path.join(script_dir, "projections.json")

    with open(file_path, "r") as json_file:
        data = json.load(json_file)
        projections = data
    return projections


def projectLayer(
    fromShape,
    outputShapefile,
    fromESPG,
    toESPG,
    geomType,
    logger=None,
    usingFiddler=False,
):
    doVerify = True
    if usingFiddler:
        doVerify = False

    driver = ogr.GetDriverByName("ESRI Shapefile")

    projections = getProjections()

    from_def_desc = projections.get(str(fromESPG))
    if not from_def_desc:
        raise Exception(f"SRID: {fromESPG} is not supported")
    from_definition = from_def_desc["definition"]["data"]

    to_def_desc = projections.get(str(toESPG))
    if not from_def_desc:
        raise Exception(f"SRID: {toESPG} is not supported")
    to_definition = to_def_desc["definition"]["data"]

    # input SpatialReference
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromProj4(from_definition)

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4(to_definition)

    if not logger is None:
        logger.debug("In proj {0}, Out proj {1}".format(inSpatialRef, outSpatialRef))

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # get the input layer
    inDataSet = driver.Open(fromShape)
    inLayer = inDataSet.GetLayer()

    # create the output layer
    # outputShapefile = toShape
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)

    outLayer = None
    if geomType.lower() == "point" or geomType.lower() == "points":
        outLayer = outDataSet.CreateLayer(
            "basemap_" + toESPG, outSpatialRef, ogr.wkbPoint
        )
    elif (
        geomType.lower() == "line"
        or geomType.lower() == "lines"
        or geomType.lower() == "polyline"
    ):
        outLayer = outDataSet.CreateLayer(
            "basemap_" + toESPG, outSpatialRef, ogr.wkbLineString
        )
    elif geomType.lower() == "polygon" or geomType.lower() == "polygons":
        outLayer = outDataSet.CreateLayer(
            "basemap_" + toESPG, outSpatialRef, ogr.wkbPolygon
        )
    # outLayer = outDataSet.CreateLayer("basemap_" + toESPG, geom_type=ogr.wkbMultiPolygon)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(
                outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i)
            )
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()

    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None
