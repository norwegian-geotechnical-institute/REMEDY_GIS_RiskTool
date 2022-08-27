coding_guide = 0 #avoids some sort of coding interpretation bugs
# Prepared for open source release August 2022

log_path = r'C:\Users\AOL\Documents\ArcGIS\BegrensSkadeCode\log'
lyr_path = r'C:\Users\AOL\Documents\ArcGIS\BegrensSkadeCode\lyr'

import arcpy
import sys
import os
import traceback
import logging.handlers
sys.path.append(log_path)
import importlib
import Utils
import Utils_arcpy
import BegrensSkade
import BegrensSkadeLib
importlib.reload(Utils)
importlib.reload(Utils_arcpy)
importlib.reload(BegrensSkade)
importlib.reload(BegrensSkadeLib)

CALCULATION_RANGE = 380

##############  SETUP LOGGERS ##############################
maxLoggerFileSize = 2 * 1024 * 1024
logger = logging.getLogger("BegrensSkade_EXCAVATION")
if not len(logger.handlers):
    logFile = log_path + "//BegrensSkadeII_ArcGISPro_EXCAVATION.log"
    hdlr = logging.handlers.RotatingFileHandler(logFile, "a", maxLoggerFileSize, 20)
    formatter = logging.Formatter("%(asctime)s %(levelname)s Thread %(thread)d %(message)s ")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
############################################################


##############  READ PARAMETERS ############################
building_polys_fl = arcpy.GetParameter(0)
excavation_polys_fl = arcpy.GetParameter(1)
output_folder = arcpy.GetParameterAsText(2)
feature_name = arcpy.GetParameterAsText(3)
coord_syst = arcpy.GetParameterAsText(4)

sr = arcpy.SpatialReference()
sr.loadFromString(coord_syst)
output_proj = sr.PCSCode

corner_name = feature_name + "_CORNER"
wall_name = feature_name + "_WALL"
building_name = feature_name + "_BUILDING"

bShortterm = arcpy.GetParameter(5)
if bShortterm:
    excavation_depth = arcpy.GetParameter(6)
    short_term_curve = arcpy.GetParameterAsText(7)
else:
    excavation_depth = None
    short_term_curve = None

bLongterm = arcpy.GetParameter(8)
if bShortterm == False and bLongterm == False:
    arcpy.AddError("Please choose Short term or Long term settlements, or both")
    sys.exit()

if bLongterm:
    dtb_raster = arcpy.GetParameter(9)
    pw_reduction_curve = arcpy.GetParameterAsText(10)
    porewp_red = arcpy.GetParameter(11)
    dry_crust_thk = arcpy.GetParameter(12)
    dep_groundwater = arcpy.GetParameter(13)
    density_sat = arcpy.GetParameter(14)
    OCR = arcpy.GetParameter(15)
    janbu_ref_stress = arcpy.GetParameter(16)
    janbu_const = arcpy.GetParameter(17)
    janbu_m = arcpy.GetParameter(18)
    consolidation_time = arcpy.GetParameter(19)
else:
    dtb_raster = None
    pw_reduction_curve = None
    dry_crust_thk = None
    dep_groundwater = None
    density_sat = None
    OCR = None
    porewp_red = None
    janbu_ref_stress = None
    janbu_const = None
    janbu_m = None
    consolidation_time = None

bVulnerability = arcpy.GetParameter(20)
if bVulnerability:
    fields = arcpy.ListFields(building_polys_fl)
    field_map = {}
    for idx, field in enumerate(fields):
        field_map[field.name] = idx - 2
    vuln_idx_count = 0

    try:
        foundation_field = field_map[arcpy.GetParameterAsText(21)]
        vuln_idx_count += 1
    except:
        foundation_field = None
    try:
        structure_field = field_map[arcpy.GetParameterAsText(22)]
        vuln_idx_count += 1
    except:
        structure_field = None
    try:
        status_field = field_map[arcpy.GetParameterAsText(23)]
        vuln_idx_count += 1
    except:
        status_field = None

    if vuln_idx_count > 0:
        arcpy.AddMessage("Vulnerability enabled")
    else:
        arcpy.AddWarning("No valid vulnerability input - disabling vulnerability!")
        bVulnerability = False
else:
    arcpy.AddMessage("Vulnerability disabled")
    foundation_field = None
    structure_field = None
    status_field = None

##############  SET PROJECTION ###########################
working_proj = Utils_arcpy.getProjCodeFromFC(building_polys_fl)
excavation_proj = Utils_arcpy.getProjCodeFromFC(excavation_polys_fl)

########  PROJECT EXCAVATION POLY IF NECESSARY ###########
excavation_polys_projected = False
if excavation_proj != working_proj:
    arcpy.AddMessage("Projecting excavation polygon..")
    excavation_polys_projected = output_folder + os.sep + "exc_proj.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_projected, working_proj)
    excavation_polys_fl = excavation_polys_projected

################ GET EXCAVATION INFO #####################
excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)
buildingsClipExtent = Utils_arcpy.getBuildingsClipExtentFromConstruction(excavation_outline_as_json, CALCULATION_RANGE, working_proj, logger)

################ EXTRACTING BUILDINGS ##################
buildings_clip = output_folder + os.sep + "buildings_clip.shp"
logger.debug("TIME - Starting extraction of buildings")
Utils_arcpy.extractBuildingsFromFL(building_polys_fl, buildingsClipExtent, buildings_clip, logger)
logger.info("TIME - Done extraction of buildings.")

############### HANDELING OF INPUT RASTER ################

if bLongterm:
    # If necessary, projects raster to the working projection
    raster_desc = arcpy.Describe(dtb_raster)
    dtb_raster_proj = raster_desc.SpatialReference.PCSCode
    if (str(dtb_raster_proj) != str(working_proj)):
        logger.info("START raster projection")
        arcpy.AddMessage("Projecting raster...")
        dtb_proj_raster = "temp_raster"
        if os.path.exists(dtb_proj_raster):
            os.remove(dtb_proj_raster)
        arcpy.ProjectRaster_management(dtb_raster, dtb_proj_raster, working_proj)
        dtb_raster = dtb_proj_raster
    logger.info("DONE raster projection")

    # Create a tif file from the raster. Necessary for input to GDAL.
    raster_desc = arcpy.Describe(dtb_raster)
    if raster_desc.extension != ".tif":
        logger.info("START raster to TIFF conversion")
        arcpy.AddMessage("Converting raster...")
        dtb_raster_tiff = output_folder + os.sep + raster_desc.name + ".tif"
        #Delete existing rasters with the same name
        if os.path.exists(dtb_raster_tiff):
            os.remove(dtb_raster_tiff)
        arcpy.RasterToOtherFormat_conversion(raster_desc.name, output_folder, "TIFF")
        dtb_raster_str = dtb_raster_tiff
    logger.info("DONE raster to TIFF conversion")
else:
    dtb_raster_str = None

############  RUN BEGRENS SKADE CORE FUNCTIONS   ##############
arcpy.AddMessage("Running mainBegrensSkade_Excavation...")
try:
    outputFiles = BegrensSkade.mainBegrensSkade_Excavation(
        logger,
        buildings_clip,
        excavation_outline_as_json,
        output_folder,
        feature_name,
        working_proj,
        output_proj,
        bShortterm,
        excavation_depth=excavation_depth,
        short_term_curve=short_term_curve,
        bLongterm=bLongterm,
        dtb_raster=dtb_raster_str,
        pw_reduction_curve=pw_reduction_curve,
        dry_crust_thk=dry_crust_thk,
        dep_groundwater=dep_groundwater,
        density_sat=density_sat,
        OCR=OCR,
        porewp_red=porewp_red,
        janbu_ref_stress=janbu_ref_stress,
        janbu_const=janbu_const,
        janbu_m=janbu_m,
        consolidation_time=consolidation_time,
        bVulnerability=bVulnerability,
        fieldNameFoundation=foundation_field,
        fieldNameStructure=structure_field,
        fieldNameStatus=status_field,
        )
except Exception:
    # Print original traceback info
    arcpy.AddError("UNEXPECTED ERROR:\n" + traceback.format_exc())
    arcpy.AddError(sys.exc_info()[1])
    sys.exit()

if excavation_polys_projected:
    arcpy.Delete_management(excavation_polys_projected)

############################### HANDLE THE RESULT ############################################
buildings_Shapefile_result = outputFiles[0]
walls_Shapefile_result = outputFiles[1]
corners_Shapefile_result = outputFiles[2]

arcpy.SelectLayerByAttribute_management(building_polys_fl, "CLEAR_SELECTION")

arcpy.AddMessage("Adding symbology layer to map...")
p = arcpy.mp.ArcGISProject("CURRENT")
pMap = p.activeMap

if bVulnerability:
    addRiskAngle = Utils.setBooleanParameter(arcpy.GetParameter(24))
    addRiskSettl = Utils.setBooleanParameter(arcpy.GetParameter(25))
addImpactAngle = Utils.setBooleanParameter(arcpy.GetParameter(26))
addImpactSettl = Utils.setBooleanParameter(arcpy.GetParameter(27))
addWalls = Utils.setBooleanParameter(arcpy.GetParameter(28))
addCorners = Utils.setBooleanParameter(arcpy.GetParameter(29))

lyr_corners = lyr_path + os.sep + "CORNER_SV_mm.lyrx"
lyr_walls = lyr_path + os.sep + "WALL_ANGLE.lyrx"
lyr_building_sv_max = lyr_path + os.sep + "BUILDING_TOTAL_SV_MAX_mm.lyrx"
lyr_building_a_max = lyr_path + os.sep + "BUILDING_TOTAL_ANGLE_MAX.lyrx"
lyr_building_risk_sv = lyr_path + os.sep + "BUILDING_RISK_SV_gdal.lyrx"
lyr_building_risk_a = lyr_path + os.sep + "BUILDING_RISK_ANGLE_gdal.lyrx"
lyr_group = lyr_path + os.sep + "GIBV_RUN_.lyrx"

lyr_group = pMap.addLayer(arcpy.mp.LayerFile(lyr_group), "TOP")[0]
lyr_group.name = feature_name

if addCorners:
    Utils_arcpy.addLayerToGroup(pMap, corners_Shapefile_result, lyr_corners, lyr_group)
if addWalls:
    Utils_arcpy.addLayerToGroup(pMap, walls_Shapefile_result, lyr_walls, lyr_group)
if addImpactSettl:
    Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_sv_max, lyr_group)
if addImpactAngle:
    Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_a_max, lyr_group)
if bVulnerability:
    if addRiskSettl:
        Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_sv, lyr_group)
    if addRiskAngle:
        Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_a,lyr_group)

logger.info("------------------------------DONE-------------------------------")



