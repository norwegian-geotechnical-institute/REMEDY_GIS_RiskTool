coding_guide = 0 #avoids some sort of coding interpretation bugs
# Prepared for open source release July 2022

lyr_path = r'C:\Users\AOL\Documents\ArcGIS\BegrensSkadeCode\lyr'

import arcpy
import sys
import os
import traceback
import logging.handlers

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

##############  READ PARAMETERS ############################
building_polys_fl = arcpy.GetParameter(0)
tunnel_poly_fl = arcpy.GetParameter(1)
output_folder = arcpy.GetParameterAsText(2)
feature_name = arcpy.GetParameterAsText(3)
coord_syst = arcpy.GetParameterAsText(4)

##############  SETUP LOGGERS ##############################
log_path = output_folder
sys.path.append(log_path)
maxLoggerFileSize = 2 * 1024 * 1024
logger = logging.getLogger("BegrensSkade_TUNNEL")
if not len(logger.handlers):
    logFile = log_path + "//BegrensSkadeII_ArcGISPro_TUNNEL.log"
    hdlr = logging.handlers.RotatingFileHandler(logFile, "a", maxLoggerFileSize, 20)
    formatter = logging.Formatter("%(asctime)s %(levelname)s Thread %(thread)d %(message)s ")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
############################################################(

output_spatial_ref = arcpy.SpatialReference()
output_spatial_ref.loadFromString(coord_syst)
output_proj = output_spatial_ref.PCSCode

corner_name = feature_name + "_CORNER"
wall_name = feature_name + "_WALL"
building_name = feature_name + "_BUILDING"

bShortterm = arcpy.GetParameter(5)
if bShortterm:
    tunnel_depth = arcpy.GetParameter(6)  # z_0
    tunnel_diameter = arcpy.GetParameter(7)  # D
    volume_loss = arcpy.GetParameter(8)
    trough_width = arcpy.GetParameter(9)  # K
else:
    tunnel_depth = None
    tunnel_diameter = None
    volume_loss = None
    trough_width = None

bLongterm = arcpy.GetParameter(10)
if bShortterm == False and bLongterm == False:
    arcpy.AddError("Please choose Short term or Long term settlements, or both")
    sys.exit()

if bLongterm:
    porewp_calc_type = arcpy.GetParameterAsText(11)
    tunnel_leakage = arcpy.GetParameter(12)  # Q
    porewp_red_at_site_m = arcpy.GetParameter(13)
    dtb_raster = arcpy.GetParameter(14)
    dry_crust_thk = arcpy.GetParameter(15)
    dep_groundwater = arcpy.GetParameter(16)
    density_sat = arcpy.GetParameter(17)
    OCR = arcpy.GetParameter(18)
    janbu_ref_stress = arcpy.GetParameter(19)
    janbu_const = arcpy.GetParameter(20)
    janbu_m = arcpy.GetParameter(21)
    consolidation_time = arcpy.GetParameter(22)
else:
    porewp_calc_type = None
    porewp_red_at_site_m = None
    tunnel_leakage = None
    dtb_raster = None
    dry_crust_thk = None
    dep_groundwater = None
    density_sat = None
    OCR = None
    janbu_ref_stress = None
    janbu_const = None
    janbu_m = None
    consolidation_time = None

bVulnerability = arcpy.GetParameter(23)
if bVulnerability:
    fields = arcpy.ListFields(building_polys_fl)
    field_map = {}
    for idx, field in enumerate(fields):
        field_map[field.name] = idx - 2
    vuln_idx_count = 0

    try:
        foundation_field = field_map[arcpy.GetParameterAsText(24)]
        vuln_idx_count += 1
    except:
        foundation_field = None
    try:
        structure_field = field_map[arcpy.GetParameterAsText(25)]
        vuln_idx_count += 1
    except:
        structure_field = None
    try:
        status_field = field_map[arcpy.GetParameterAsText(26)]
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

##############  GET INPUT PROJECTIONS ###########################
building_spatial_ref = arcpy.Describe(building_polys_fl).spatialReference
tunnel_spatial_ref = arcpy.Describe(tunnel_poly_fl).spatialReference

########  GET EXCAVATION AND BUILDINGS ON SAME PROJECTION ###########
tunnel_poly_matched = False
if tunnel_spatial_ref != building_spatial_ref:
    arcpy.AddMessage("Matching input projections before clip..")
    tunnel_poly_matched = output_folder + os.sep + "tun_match.shp"
    arcpy.Project_management(tunnel_poly_fl, tunnel_poly_matched, building_spatial_ref)
    tunnel_poly_fl = tunnel_poly_matched

################ GET TUNNEL INFO #####################
tunnel_outline_as_json = Utils_arcpy.getConstructionAsJson(tunnel_poly_fl)
buildingsClipExtent = Utils_arcpy.getBuildingsClipExtentFromConstruction(tunnel_outline_as_json, CALCULATION_RANGE, building_spatial_ref, logger)

################ EXTRACTING BUILDINGS ##################
buildings_clip = output_folder + os.sep + "buildings_clip.shp"
logger.debug("TIME - Starting extraction of buildings")
Utils_arcpy.extractBuildingsFromFL(building_polys_fl, buildingsClipExtent, buildings_clip, logger)
logger.info("TIME - Done extraction of buildings.")

#arcpy.SelectLayerByAttribute_management(building_polys_fl, "CLEAR_SELECTION")

######  PROJECT BUILDING AND EXCAVATION TO OUTPUT  ########
building_polys_projected = False
tunnel_poly_projected = False
buildings_clip_projected = False

if building_spatial_ref != output_spatial_ref:

    arcpy.AddMessage("Projecting bulidings polygon..")
    buildings_clip_projected = output_folder + os.sep + "buil_proj.shp"
    arcpy.Project_management(buildings_clip, buildings_clip_projected, output_proj)
    building_polys_fl = buildings_clip_projected

    arcpy.AddMessage("Projecting tunnel polygon..")
    tunnel_poly_projected = output_folder + os.sep + "tun_proj.shp"
    arcpy.Project_management(tunnel_poly_fl, tunnel_poly_projected, output_proj)
    tunnel_poly_fl = tunnel_poly_projected

    tunnel_outline_as_json = Utils_arcpy.getConstructionAsJson(tunnel_poly_fl)

else:
    building_polys_fl = buildings_clip


############### HANDELING OF INPUT RASTER ################

dtb_proj_raster = False
if bLongterm:
    # If necessary, projects raster to the working projection
    raster_desc = arcpy.Describe(dtb_raster)
    dtb_raster_proj = raster_desc.SpatialReference.PCSCode
    if (str(dtb_raster_proj) != str(output_proj)):
        logger.info("START raster projection")
        arcpy.AddMessage("Projecting bedrock raster...")
        dtb_proj_raster = "temp_raster"
        if os.path.exists(dtb_proj_raster):
            os.remove(dtb_proj_raster)
        arcpy.ProjectRaster_management(dtb_raster, dtb_proj_raster, output_proj)
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
arcpy.AddMessage("Running mainBegrensSkade_Tunnel...")
try:
    outputFiles = BegrensSkade.mainBegrensSkade_Tunnel(
        logger,
        building_polys_fl,
        tunnel_outline_as_json,
        output_folder,
        feature_name,
        output_proj,
        bShortterm,
        tunnel_depth=tunnel_depth,
        tunnel_diameter=tunnel_diameter,
        volume_loss=volume_loss,
        trough_width=trough_width,
        bLongterm=bLongterm,
        porewp_calc_type=porewp_calc_type,
        tunnel_leakage=tunnel_leakage,
        porewp_red_at_site_m=porewp_red_at_site_m,
        dtb_raster=dtb_raster_str,
        dry_crust_thk=dry_crust_thk,
        dep_groundwater=dep_groundwater,
        density_sat=density_sat,
        OCR=OCR,
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

if dtb_proj_raster:
    try:
        arcpy.Delete_management(output_folder + os.sep + dtb_proj_raster + ".tif")
    except:
        arcpy.AddWarning("Failed to delete temporary bedrock raster!")

arcpy.Delete_management(buildings_clip)
if tunnel_poly_matched:
    arcpy.Delete_management(tunnel_poly_matched)
if building_polys_projected:
    arcpy.Delete_management(building_polys_projected)
if buildings_clip_projected:
    arcpy.Delete_management(buildings_clip_projected)
if tunnel_poly_projected:
    arcpy.Delete_management(tunnel_poly_projected)

############################### HANDLE THE RESULT ############################################
buildings_Shapefile_result = outputFiles[0]
walls_Shapefile_result = outputFiles[1]
corners_Shapefile_result = outputFiles[2]


arcpy.AddMessage("Adding symbology layer to map...")
p = arcpy.mp.ArcGISProject("CURRENT")
pMap = p.activeMap

if bVulnerability:
    addRiskAngle = Utils.setBooleanParameter(arcpy.GetParameter(27))
    addRiskSettl = Utils.setBooleanParameter(arcpy.GetParameter(28))
addImpactAngle = Utils.setBooleanParameter(arcpy.GetParameter(29))
addImpactSettl = Utils.setBooleanParameter(arcpy.GetParameter(30))
addWalls = Utils.setBooleanParameter(arcpy.GetParameter(31))
addCorners = Utils.setBooleanParameter(arcpy.GetParameter(32))

lyr_corners = lyr_path + os.sep + "CORNER_SV.lyrx"
lyr_walls = lyr_path + os.sep + "WALL_ANGLE.lyrx"
lyr_building_sv_max = lyr_path + os.sep + "BUILDING_TOTAL_SV_MAX_mm.lyrx"
lyr_building_a_max = lyr_path + os.sep + "BUILDING_TOTAL_ANGLE_MAX.lyrx"
lyr_building_risk_sv = lyr_path + os.sep + "BUILDING_RISK_SV_gdal.lyrx"
lyr_building_risk_a = lyr_path + os.sep + "BUILDING_RISK_ANGLE_gdal.lyrx"
lyr_group = lyr_path + os.sep + "GIBV_RUN.lyrx"

lyr_group = pMap.addLayer(arcpy.mp.LayerFile(lyr_group), "TOP")[0]
lyr_group.name = feature_name

if addImpactSettl:
    Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_sv_max, lyr_group,
                                name="Maximum settlement each building experience")
if addImpactAngle:
    Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_a_max, lyr_group,
                                name="Maximum tilt each building experience")
if bVulnerability:
    if addRiskSettl:
        Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_sv, lyr_group,
                                    name="Building risk category based on settlement damage")
    if addRiskAngle:
        Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_a, lyr_group,
                                    name="Building risk category based on tilt damage")
if addWalls:
    Utils_arcpy.addLayerToGroup(pMap, walls_Shapefile_result, lyr_walls, lyr_group,
                                name="Tilt on building walls")
if addCorners:
    Utils_arcpy.addLayerToGroup(pMap, corners_Shapefile_result, lyr_corners, lyr_group,
                                name="Settlement on building corners")

#arcpy.SelectLayerByAttribute_management(buildings_Shapefile_result, "CLEAR_SELECTION")


logger.info("------------------------------DONE-------------------------------")

