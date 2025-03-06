import sys
import os
import csv
import numpy as np
from datetime import datetime
from osgeo import gdal
from osgeo.gdalconst import *
import importlib
import Utils
import BegrensSkadeLib
importlib.reload(Utils)
importlib.reload(BegrensSkadeLib)

def mainBegrensSkade_Excavation(
    logger, #logging.handlers object "logging.getLogger("name")" for debugging
    buildingsFN, #Shapefile containing buildings
    excavationJson, #Excavation outline as json
    output_ws, #Output folder for result storage
    feature_name, #Name of the result shapefiles
    output_proj, #PCSC Code of the Coordinate system for result shapefiles, refer to Utils.epsgMapping{}
    bShortterm, #Boolean flag for short term settlement calculations
    excavation_depth=10, #Depth of excavation, necessary for short term calculation
    short_term_curve="1 % av byggegropdybde", #Short term settlement curve, refer to manual for options
    bLongterm=False, #Bolean flag for long term settlement calculations
    dtb_raster="", #Depth to bedrock tiff raster for long term settlement (path string)
    porewp_red_m=10, #Reduction of porewater pressure close to excavation in meters
    dry_crust_thk=5, #Thickness of overburden not affected by porewater drawdown
    dep_groundwater=3, #Depht to groundwater table
    density_sat=18.5, #Soil saturation density
    OCR=1.2, #Over consolidation ratio
    janbu_ref_stress=0, #Janbu reference stress
    janbu_const=4, #Janbu constant, refer to manual
    janbu_m=15, #Janbu compression modulus
    consolidation_time=1000, #Consolidation time in years
    bVulnerability=False, #Boleean flag for building vulnerability
    fieldNameFoundation=None, #Shapefile field for building foundation, refer to manual
    fieldNameStructure=None, #Shapefile field for building structure, refer to manual
    fieldNameStatus=None, #Shapefile field for building status, refer to manual
) -> []:


    parameter_log = open(output_ws + os.sep + "GIBV_Excavation_params_" + feature_name + ".txt", "w")
    parameter_log.write("PARAMETER LOG FOR GIBV Excavation \n")
    parameter_log.write(feature_name+ "\n")
    parameter_log.write("Time of run: " + str(datetime.today()) + "\n")
    parameter_log.write("-----------------------------\n")
    parameter_log.write("Output coordinate system code: " + str(output_proj) + "\n")
    parameter_log.write("Short term enabled: " + str(bShortterm)+ "\n")
    if bShortterm:
        parameter_log.write("Excavation depth: " + str(excavation_depth) + "\n")
        parameter_log.write("Short term curve: "+ str(short_term_curve) + "\n")
    parameter_log.write("Long term enabled: " + str(bLongterm) + "\n")
    if bLongterm:
        parameter_log.write("Dtb raster: " + str(dtb_raster) + "\n")
        parameter_log.write("Depth unaffected pore water pressure (\"dry crust\"): " + str(dry_crust_thk) + "\n")
        parameter_log.write("Depth groundwater table: " + str(dep_groundwater)+"\n")
        parameter_log.write("Density: " + str(density_sat) +"\n")
        parameter_log.write("OCR: " +str(OCR)+ "\n")
        parameter_log.write("Porewater pressure reduction (m): " + str(porewp_red_m)+"\n")
        parameter_log.write("Janbu reference stress: " +str(janbu_ref_stress)+ "\n")
        parameter_log.write("Janbu constant: " + str(janbu_const)+ "\n")
        parameter_log.write("Janbu m: " + str(janbu_m)+  "\n")
        parameter_log.write("Consolidation time: " + str(consolidation_time)+"\n")
    parameter_log.write("Vulnerability enabled: " +str(bVulnerability)+ "\n")
    if bVulnerability:
        parameter_log.write("Foundation field column: " + str(fieldNameFoundation)+ "\n")
        parameter_log.write("Structure field column: " + str(fieldNameStructure)+"\n")
        parameter_log.write("Status field column: " + str(fieldNameStatus)+"\n")
    parameter_log.close()


    # A - Til berg - Direktefundamentert, peler
    # B - På løsmasser - Hel plate (betong, såle)
    # C - På løsmasser - Stripefundament (heller)
    # D - På løsmasser - Punkt- og trefundamenter (banketter)

    FOUNDATION_A = [
        "To bedrock",
        "Peler",
        "A - Til berg - Direktefundamentert, peler",
    ]  # 0
    FOUNDATION_B = [
        "Raft",
        "Betong",
        "B - På løsmasser - Hel plate (betong, såle)",
    ]  # 5
    FOUNDATION_C = [
        "Strip",
        "Grunnmur",
        "C - På løsmasser - Stripefundament (heller)",
    ]  # 20
    FOUNDATION_D = [
        "Wooden piles",
        "Trepeler",
        "D - På løsmasser - Punkt- og trefundamenter (banketter)",
    ]  # 50

    STRUCTURE_A = ["Steel", "A - Stål"]
    STRUCTURE_B = ["Reinforced concrete", "B - Armert betong"]
    STRUCTURE_C = ["Mixed", "C - Tre eller varierende"]
    STRUCTURE_D = ["Masonry", "D - Murstein eller spesiell type"]

    STATUS_A = ["Excellent", "A - Meget god tilstand"]
    STATUS_B = ["Good", "B - God tilstand"]
    STATUS_C = ["Medium", "C - Brukbar tilstand"]
    STATUS_D = ["Bad", "D - Dårlig"]

    VULN_WEIGHT_LENGTH = 0.75
    VULN_WEIGHT_SHAPE = 0.75
    VULN_WEIGHT_STRUCTURE = 1.0
    VULN_WEIGHT_FOUNDATION = 1.5
    VULN_WEIGHT_STATUS = 1.0

    WALL_CORNER_ANGLE_THRESHOLD = 5  # degrees
    CONSTR_RESAMPLE_LEN = 2.0

    # "Loading construction zone and building polygons to memory..."

    logger.info("-------------INPUT PARAMETERS:-----------")
    logger.info("janbu_const: " + str(janbu_const))
    logger.info("consolidation_time: " + str(consolidation_time))
    logger.info("fieldname_foundation: " + str(fieldNameFoundation))
    logger.info("fieldname_structure: " + str(fieldNameStructure))
    logger.info("fieldname_Status: " + str(fieldNameStatus))
    logger.info("---------------------------------------")

    # _array(corner)
    logger.debug( "Calling get_construction_corners_from_ArcGIS_json with json: {}".format(excavationJson))
    construction_area_corners = BegrensSkadeLib.get_construction_corners_from_ArcGIS_json(excavationJson, CONSTR_RESAMPLE_LEN, logger)
    # array(buildings)
    logger.info("TIME - getting buildings from shapefile")

    if bLongterm:
        buildings = BegrensSkadeLib.get_buildings_with_dtb(
            buildingsFN,
            dtb_raster,
            fieldNameFoundation=fieldNameFoundation,
            fieldNameStructure=fieldNameStructure,
            fieldNameStatus=fieldNameStatus,
            logger=logger,
        )
    else:
        buildings = BegrensSkadeLib.get_buildings(
            buildingsFN,
            fieldNameFoundation=fieldNameFoundation,
            fieldNameStructure=fieldNameStructure,
            fieldNameStatus=fieldNameStatus,
            logger=logger,
        )

    logger.info("TIME - gotten " + str(len(buildings)) + " buildings from shapefile")

    count_adj = 0
    logger.info("TIME - calculate settlements")

    buildingsWithBedrock = []

    for building in buildings:
        # filter corner duplicates and straight-wall corners

        building.filter_duplicates()
        building.filter_straights(WALL_CORNER_ANGLE_THRESHOLD)

        # if bLongterm, add bedrock depth information
        if bLongterm:
            building.corners = BegrensSkadeLib.appendZValuesFromRaster(building.corners, dtb_raster, logger)
            if building.corners is None:
                # logger.debug("Skipping building")
                continue
            else:
                buildingsWithBedrock.append(building)

        max_sv_short = 0.0
        max_sh_short = 0.0
        max_sv_total = 0.0
        near_closest = 9999.0
        near_furthest = 0.0

        for corner in building.corners:

            # Evaluating distance from every building corner to construction zone...
            near_dist, near_angle = BegrensSkadeLib.near_analysis(corner.x, corner.y, construction_area_corners)

            corner.near_angle = near_angle
            corner.near_dist = near_dist

            near_closest = min(near_dist, near_closest)
            near_furthest = max(near_dist, near_furthest)

            # Evaluating settlements at corners

            if bShortterm:
                if short_term_curve in ["Wall not to bedrock - regression",
                                        "0,5 % av byggegropdybde",
                                        "Spunt installert til berg med høy sikkerhet",
                                        "Norm_setning_0.5"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_a(near_dist, excavation_depth)
                elif short_term_curve in ["Wall not to bedrock - discrete",
                                          "1 % av byggegropdybde",
                                          "Spunt installert til berg med lav sikkerhet",
                                          "Norm_setning_1"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_b(near_dist, excavation_depth)
                elif short_term_curve in ["Tie-back anchors - regression",
                                          "2 % av byggegropdybde",
                                          "Svevespunt høy sikkerhet",
                                          "Norm_setning_2"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_c(near_dist, excavation_depth)
                elif short_term_curve in ["Tie-back anchors - regression",
                                          "3 % av byggegropdybde",
                                          "Svevespunt lav sikkerhet",
                                          "Norm_setning_3"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_d(near_dist, excavation_depth)
                else:
                    raise Exception("Not a valid regression curve: " + str(short_term_curve))

            else:
                sv_short = 0

            porewp_red_atdist = -1

            if bLongterm:

                #Oct 2024: New method for porewater pressure reduction calculation (DOI 10.1061/JGGEFK/GTENG-12490)

                porewp_red_atdist = 10*max((porewp_red_m - near_dist/50), 0) #pore water pressure reduction in kPa

                sv_long, red_adj = BegrensSkadeLib.get_sv_long_janbu(
                    corner.dtb,
                    dry_crust_thk,
                    dep_groundwater,
                    density_sat,
                    OCR,
                    porewp_red_atdist, #pore water pressure reduction in kPa
                    janbu_ref_stress,
                    janbu_const,
                    janbu_m,
                    consolidation_time,
                )
                if red_adj:
                    count_adj += 1
            else:
                sv_long = 0

            corner.sv_short = sv_short
            corner.sh_short = 0
            corner.sv_long = sv_long
            corner.sv_tot = sv_short + sv_long
            corner.porewp_red = porewp_red_atdist

            max_sv_short = max(max_sv_short, sv_short)
            max_sh_short = 0
            max_sv_total = max(max_sv_total, (sv_short + sv_long))

        building.create_walls()

        max_angle = 0.0
        max_strain = 0.0
        max_principal_strain = 0.0

        # logger.debug("Number of walls:{}".format(len(building.walls)))
        for wall in building.walls:
            max_angle = max(max_angle, abs(wall.slope_angle))
            max_strain = max(max_strain, abs(wall.hor_tensile_strain))
            max_principal_strain = max(max_principal_strain, abs(wall.principal_tensile_strain))

        # logger.debug(str(building.bid) + ", " + str(building.area) + ", "+str(building.circumf) + ", " + str(max_sv_short) + ", "+str(
        #    max_sh_short) + ", "+str(max_sv_total) + ", " + str(max_angle) + ", " + str(max_strain) + ", " + str(max_principal_strain))
        building.max_sv_short = max_sv_short
        building.max_sh_short = max_sh_short
        building.max_sv_total = max_sv_total
        building.max_angle = max_angle
        building.max_strain = max_strain
        building.max_principal_strain = max_principal_strain

        # Evaluate vulnerability
        if bVulnerability:

            if (fieldNameFoundation):
                foundation = building.foundation

                if foundation in FOUNDATION_A:
                    foundation_cvi = 0
                elif foundation in FOUNDATION_B:
                    foundation_cvi = 5
                elif foundation in FOUNDATION_C:
                    foundation_cvi = 20
                elif foundation in FOUNDATION_D:
                    foundation_cvi = 50
                else:
                    foundation_cvi = 50  # default
                building.foundation = foundation_cvi


            if fieldNameStructure:
                structure = building.structure

                if structure in STRUCTURE_A:
                    structure_cvi = 0
                elif structure in STRUCTURE_B:
                    structure_cvi = 5
                elif structure in STRUCTURE_C:
                    structure_cvi = 20
                elif structure in STRUCTURE_D:
                    structure_cvi = 50
                else:
                    structure_cvi = 50  # default

                building.structure = structure_cvi

            if fieldNameStatus:
                status = building.status

                if status in STATUS_A:
                    status_cvi = 0
                elif status in STATUS_B:
                    status_cvi = 5
                elif status in STATUS_C:
                    status_cvi = 20
                elif status in STATUS_D:
                    status_cvi = 50
                else:
                    status_cvi = 50  # default

                building.status = status_cvi

            # logger.debug("Got status structure and foundation cvi for building: {}".format(building.bid))
            isosquare = 16 * building.area / (building.circumf ** 2)
            radial_buil_len = near_furthest - near_closest
            # logger.debug("building.area: " + str(building.area) + ", building.circumf: " + str(building.circumf) + ", isosquare " + str(isosquare))
            # logger.debug("near_furthest: " + str(near_furthest) + ", near_closest: " + str(near_closest) + ", radial_buil_len: " + str(radial_buil_len))
            # logger.debug("building.foundation: " + str(building.foundation) + ", building.structure: " + str(building.structure) + ", building.status: " + str(building.status))

            length_cvi = BegrensSkadeLib.get_buil_len_cvi(radial_buil_len)
            shape_cvi = BegrensSkadeLib.get_buil_shape_cvi(isosquare)

            try:
                vulnerability = (VULN_WEIGHT_LENGTH * length_cvi + VULN_WEIGHT_SHAPE * shape_cvi)
                maxsum = 50.0 * (VULN_WEIGHT_LENGTH + VULN_WEIGHT_SHAPE)

                if fieldNameStructure:
                    vulnerability += structure_cvi * VULN_WEIGHT_STRUCTURE
                    maxsum += 50.0 * VULN_WEIGHT_STRUCTURE
                if fieldNameFoundation:
                    vulnerability += foundation_cvi * VULN_WEIGHT_FOUNDATION
                    maxsum += 50.0 * VULN_WEIGHT_FOUNDATION
                if fieldNameStatus:
                    vulnerability += status_cvi * VULN_WEIGHT_STATUS
                    maxsum += 50.0 * VULN_WEIGHT_STATUS

                vulnerability = vulnerability / maxsum
                building.vulnerability = vulnerability

                vuln_cvi = BegrensSkadeLib.get_buil_vuln_cvi(vulnerability)

                impact_angle_cvi = BegrensSkadeLib.get_buil_impact_angle_cvi(max_angle)
                impact_totset_cvi = BegrensSkadeLib.get_buil_impact_totset_cvi(max_sv_total)

                # logger.debug( "Calculating risk_totset and risk_angle for building: {}".format(building.bid))

                building.risk_totset = BegrensSkadeLib.get_risk_cvi(vuln_cvi, impact_totset_cvi)
                building.risk_angle = BegrensSkadeLib.get_risk_cvi(vuln_cvi, impact_angle_cvi)
            except Exception as e:
                logger.error("Error, type; {0}, error args: {1}, error: {2}".format(type(e), e.args, e))
                # logger.error(e.args)
                # logger.error(e)

    now = datetime.now()  # current date and time
    date_time_str = now.strftime("_%Y%m%d_%H%M%S")

    corner_name = feature_name + date_time_str + "_C"
    wall_name = feature_name + date_time_str + "_W"
    building_name = feature_name + date_time_str + "_B"
    if bShortterm and bLongterm:
        corner_name += "_SL"
        wall_name += "_SL"
        building_name += "_SL"
    elif bLongterm:
        corner_name += "_L"
        wall_name += "_L"
        building_name += "_L"
    elif bShortterm:
        corner_name += "_S"
        wall_name += "_S"
        building_name += "_S"


    logger.info("TIME - done calculating settlements")

    filterValue = feature_name

    # logger.debug("Delete threads returned")
    logger.debug( "Length buildings {0}, length buildings with bedrock: {1}".format(len(buildings), len(buildingsWithBedrock)))

    if bLongterm:
        buildings = buildingsWithBedrock

    logger.info("TIME - writing results to shape")
    building_shapefile = output_ws + os.sep + building_name + ".shp"
    building_shapefile_prj = output_ws + os.sep + building_name + "_prj.shp"
    BegrensSkadeLib.writeBuildingsToShape(building_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(building_shapefile,building_shapefile_prj,str(working_proj), str(output_proj), "polygon")

    wall_shapefile = output_ws + os.sep + wall_name + ".shp"
    wall_shapefile_prj = output_ws + os.sep + wall_name + "_prj.shp"
    BegrensSkadeLib.writeWallsToShape(wall_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(wall_shapefile,wall_shapefile_prj,str(working_proj), str(output_proj),"line")

    corner_shapefile = output_ws + os.sep + corner_name + ".shp"
    corner_shapefile_prj = output_ws + os.sep + corner_name + "_prj.shp"
    BegrensSkadeLib.writeCornersToShape(corner_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(corner_shapefile,corner_shapefile_prj,str(working_proj), str(output_proj),"point")
    logger.info("TIME - written all shapefiles")

    return [building_shapefile, wall_shapefile, corner_shapefile]


def mainBegrensSkade_ImpactMap(
    logger, #logging.handlers object "logging.getLogger("name")" for debugging
    excavationJson, #Excavation outline as json
    output_ws, #Output folder for result storage
    output_name, #Name of the result raster
    CALCULATION_RANGE, #Range of longterm settlement calculation
    output_proj, #PCSC Code of the Coordinate system for result shapefiles, refer to Utils.epsgMapping{}
    dtb_raster="", #Depth to bedrock tiff raster for long term settlement (path string)
    porewp_red_m=10,  # Reduction of porewater pressure close to excavation in meters
    dry_crust_thk=5, #Thickness of overburden not affected by porewater drawdown
    dep_groundwater=3, #Depht to groundwater table
    density_sat=18.5, #Soil saturation density
    OCR=1.2, #Over consolidation ratio
    porewp_red=50, #Reduction of porewater pressure close to excavation
    janbu_ref_stress=0, #Janbu reference stress
    janbu_const=4, #Janbu constant, refer to manual
    janbu_m=15, #Janbu compression modulus
    consolidation_time=1000, #Consolidation time in years
    bShortterm = False, #Boolean flag for shortterm settlement calculation
    excavation_depth = None, #Excavation depth for short-term calculations
    short_term_curve = "1 % av byggegropdybde", #Short term settlement curve, refer to manual for options

) -> []:

    parameter_log = open(output_ws + os.sep + "GIBV_ImpactMap_params_" + output_name + ".txt", "w")
    parameter_log.write("PARAMETER LOG FOR GIBV ImpactMap \n")
    parameter_log.write(output_name+ "\n")
    parameter_log.write("Time of run: " + str(datetime.today()) + "\n")
    parameter_log.write("-----------------------------\n")
    parameter_log.write("Output coordinate system code: " + str(output_proj) + "\n")
    parameter_log.write("Short term enabled: " + str(bShortterm)+ "\n")
    if bShortterm:
        parameter_log.write("Excavation depth: " + str(excavation_depth) + "\n")
        parameter_log.write("Short term curve: "+ str(short_term_curve) + "\n")
    parameter_log.write("Long term enabled: True (default)\n")
    parameter_log.write("Dtb raster: " + str(dtb_raster) + "\n")
    parameter_log.write("Depth unaffected pore water pressure (\"dry crust\"): " + str(dry_crust_thk) + "\n")
    parameter_log.write("Depth groundwater table: " + str(dep_groundwater)+"\n")
    parameter_log.write("Density: " + str(density_sat) +"\n")
    parameter_log.write("OCR: " +str(OCR)+ "\n")
    parameter_log.write("Porewater pressure reduction (m): " + str(porewp_red_m)+"\n")
    parameter_log.write("Janbu reference stress: " +str(janbu_ref_stress)+ "\n")
    parameter_log.write("Janbu constant: " + str(janbu_const)+ "\n")
    parameter_log.write("Janbu m: " + str(janbu_m)+  "\n")
    parameter_log.write("Consolidation time: " + str(consolidation_time)+"\n")
    parameter_log.close()

    CONSTR_RESAMPLE_LEN = 2.0

    # "Loading construction zone to memory..."

    logger.debug("Get construction corners from json")
    #logger.debug("Calling get_construction_corners_from_ArcGIS_json with json: {}".format(excavationJson))
    construction_area_corners = BegrensSkadeLib.get_construction_corners_from_ArcGIS_json(excavationJson, CONSTR_RESAMPLE_LEN, logger)

    logger.info("START calculate settlements")

    dataset = gdal.Open(str(dtb_raster))
    band = dataset.GetRasterBand(1)
    transform = dataset.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    inData = band.ReadAsArray(0, 0, cols, rows)
    outData = np.zeros((rows, cols), np.float32)

    count_crust = 0
    count_clay = 0
    min_near_dist = 9999
    max_sv_long = 0

    tot_px = rows*cols
    px = 0
    progress = 0

    #logger.info("xOrigin, yOrigin, rows, cols, pixelWidht, pixelHeight: "+ str(xOrigin)+ ", "+ str(yOrigin)+ ", "+ str(rows)+ ", "+ str(cols)+ ", "+ str(pixelWidth)+ ", "+ str(pixelHeight))

    for row in range(0, rows):
        for col in range(0, cols):
            x = xOrigin + col * pixelWidth
            y = yOrigin - row * pixelHeight
            near_dist_sqr = BegrensSkadeLib.near_analysis_sqr(x, y, construction_area_corners)
            if near_dist_sqr > CALCULATION_RANGE**2:
                px += 1
                new_progress = int(100 * px / tot_px)
                if new_progress > progress:
                    progress = new_progress
                    logger.info("Progress: " + str(progress) + " %. Outside.")
                continue
            dtb = inData[row][col]

            if dtb > 500 or dtb < -5:
                logger.info("DTB outside range")
                continue
            near_dist = np.sqrt(near_dist_sqr)
            min_near_dist = min(min_near_dist, near_dist)

            if dtb <= dry_crust_thk:
                sv_long = 0.0
                count_crust += 1

            else:
                # Evaluating Janbu long term settlements

                #Oct 2024: New method for porewater pressure reduction calculation (DOI 10.1061/JGGEFK/GTENG-12490)

                porewp_red_atdist = 10*max((porewp_red_m - near_dist/50), 0) #pore water pressure reduction in meters

                if porewp_red_atdist > 0:
                    sv_long, red_adj = BegrensSkadeLib.get_sv_long_janbu(
                        dtb,
                        dry_crust_thk,
                        dep_groundwater,
                        density_sat,
                        OCR,
                        porewp_red_atdist,
                        janbu_ref_stress,
                        janbu_const,
                        janbu_m,
                        consolidation_time,
                    )
                else:
                    sv_long = 0

                #if count_clay < 50:
                #    logger.info( "row, col, SV_LONG, near_dist, dtb: "+ str(row)+ ","+ str(col)+ ","+ str(sv_long)+ ", "+ str(near_dist)+ ", "+ str(dtb))
                count_clay += 1
                max_sv_long = max(max_sv_long, sv_long)

            if bShortterm:
                # Evaluating short term settlements
                if short_term_curve in ["Wall not to bedrock - regression",
                                        "0,5 % av byggegropdybde",
                                        "Spunt installert til berg med høy sikkerhet",
                                        "Norm_setning_0.5"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_a(near_dist, excavation_depth)
                elif short_term_curve in ["Wall not to bedrock - discrete",
                                          "1 % av byggegropdybde",
                                          "Spunt installert til berg med lav sikkerhet",
                                          "Norm_setning_1"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_b(near_dist, excavation_depth)
                elif short_term_curve in ["Tie-back anchors - regression",
                                          "2 % av byggegropdybde",
                                          "Svevespunt høy sikkerhet",
                                          "Norm_setning_2"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_c(near_dist, excavation_depth)
                elif short_term_curve in ["Tie-back anchors - regression",
                                          "3 % av byggegropdybde",
                                          "Svevespunt lav sikkerhet",
                                          "Norm_setning_3"]:
                    sv_short, W = BegrensSkadeLib.get_sv_short_d(near_dist, excavation_depth)
                else:
                    raise Exception("Not a valid regression curve: " + str(short_term_curve))

            else:
                sv_short = 0.0

            outData[row, col] = sv_short + sv_long
            px += 1
            new_progress = int(100*px/tot_px)
            if  new_progress > progress:
                progress = new_progress
                logger.info("Progress: " + str(progress) + " %")

    logger.info("Px, tot px: " + str(px) + ", " + str(tot_px))
    logger.info("Count crust: " + str(count_crust))
    logger.info("Count clay: " + str(count_clay))
    logger.info("Min near dist: " + str(min_near_dist))
    logger.info("Max near dist: " + str(max_sv_long))
    logger.info("STOP calculating settlements")
    logger.info("START writing results")

    # An alternative to store data
    driver = dataset.GetDriver()
    outFile = output_ws + os.sep + output_name + ".tif"
    logger.info(outFile)
    outDs = driver.Create(outFile, cols, rows, 1, gdal.GDT_Float32)
    if outDs is None:
        print("Could not create " + output_name)
        sys.exit(1)
    outBand = outDs.GetRasterBand(1)

    # write the data
    outBand.WriteArray(outData, 0, 0)
    outBand.SetNoDataValue(-99)

    # Set the geotransform and projection from the input raster because the input
    # is already in the correct projection
    outDs.SetGeoTransform(dataset.GetGeoTransform())
    outDs.SetProjection(dataset.GetProjection())
    outDs.FlushCache()

    del outData

    logger.info("DONE writing results")

    return outFile


def mainBegrensSkade_Tunnel(
    logger, #logging.handlers object "logging.getLogger("name")" for debugging
    buildingsFN, #Shapefile containing buildings
    tunnelJson, #Tunnel outline as json
    output_ws, #Output folder for result storage
    feature_name, #Name of the result shapefiles
    output_proj, #PCSC Code of the Coordinate system for result shapefiles, refer to Utils.epsgMapping{}
    bShortterm, #Boolean flag for short term settlement calculations
    tunnel_depth="", #Depth to tunnel
    tunnel_diameter="", #Tunnel diameter
    volume_loss="", #Net volume loss (tunnel cross section), refer to manual
    trough_width="", #Estimated trough width, refer to manual
    bLongterm=False, #Boolean flag for long term settlements
    tunnel_leakage="", #Leakage rate for long term settlement calculation
    porewp_calc_type="Mean", #Type of porewater reduction calculation, refer to manual
    porewp_red_at_site_m=0, #Porewater reduction at tunnel site, refer to manual
    dtb_raster="", #Depth to bedrock tiff raster for long term settlement (path string)
    dry_crust_thk=5, #Thickness of overburden not affected by porewater drawdown
    dep_groundwater=3, #Depht to groundwater table
    density_sat=18.5, #Soil saturation density
    OCR=1.2, #Over consolidation ratio
    janbu_ref_stress=0, #Janbu reference stress
    janbu_const=4, #Janbu constant, refer to manual
    janbu_m=15, #Janbu compression modulus
    consolidation_time=1000, #Consolidation time in years
    bVulnerability=False, #Boleean flag for building vulnerability
    fieldNameFoundation=None, #Shapefile field for building foundation, refer to manual
    fieldNameStructure=None, #Shapefile field for building structure, refer to manual
    fieldNameStatus=None, #Shapefile field for building status, refer to manual
) -> []:

    parameter_log = open(output_ws + os.sep + "GIBV_Tunnel_params_" + feature_name + ".txt", "w")
    parameter_log.write("PARAMETER LOG FOR GIBV Tunnel \n")
    parameter_log.write(feature_name+ "\n")
    parameter_log.write("Time of run: " + str(datetime.today()) + "\n")
    parameter_log.write("-----------------------------\n")
    parameter_log.write("Output coordinate system code: " + str(output_proj) + "\n")
    parameter_log.write("Short term enabled: " + str(bShortterm)+ "\n")
    if bShortterm:
        parameter_log.write("Tunnel depth: " + str(tunnel_depth) + "\n")
        parameter_log.write("Tunnel diameter: "+ str(tunnel_diameter) + "\n")
        parameter_log.write("Volume loss: " + str(volume_loss) + "\n")
        parameter_log.write("Trough width: " + str(trough_width) + "\n")
    parameter_log.write("Long term enabled: " + str(bLongterm) + "\n")
    if bLongterm:
        parameter_log.write("Dtb raster: " + str(dtb_raster) + "\n")
        parameter_log.write("Tunnel leakage: " + str(tunnel_leakage) + "\n")
        parameter_log.write("Depth unaffected pore water pressure (\"dry crust\"): " + str(dry_crust_thk) + "\n")
        parameter_log.write("Depth groundwater table: " + str(dep_groundwater)+"\n")
        parameter_log.write("Density: " + str(density_sat) +"\n")
        parameter_log.write("OCR: " +str(OCR)+ "\n")
        parameter_log.write("Porewater reduction calculation type: " + str(porewp_calc_type)+"\n")
        parameter_log.write("Porewater pressure reduction at site (m): " + str(porewp_red_at_site_m) + "\n")
        parameter_log.write("Janbu reference stress: " +str(janbu_ref_stress)+ "\n")
        parameter_log.write("Janbu constant: " + str(janbu_const)+ "\n")
        parameter_log.write("Janbu m: " + str(janbu_m)+  "\n")
        parameter_log.write("Consolidation time: " + str(consolidation_time)+"\n")
    parameter_log.write("Vulnerability enabled: " +str(bVulnerability)+ "\n")
    if bVulnerability:
        parameter_log.write("Foundation field column: " + str(fieldNameFoundation)+ "\n")
        parameter_log.write("Structure field column: " + str(fieldNameStructure)+"\n")
        parameter_log.write("Status field column: " + str(fieldNameStatus)+"\n")
    parameter_log.close()

    # A - Til berg - Direktefundamentert, peler
    # B - På løsmasser - Hel plate (betong, såle)
    # C - På løsmasser - Stripefundament (heller)
    # D - På løsmasser - Punkt- og trefundamenter (banketter)

    FOUNDATION_A = [
        "To bedrock",
        "Peler",
        "A - Til berg - Direktefundamentert, peler",
    ]  # 0
    FOUNDATION_B = [
        "Raft",
        "Betong",
        "B - På løsmasser - Hel plate (betong, såle)",
    ]  # 5
    FOUNDATION_C = [
        "Strip",
        "Grunnmur",
        "C - På løsmasser - Stripefundament (heller)",
    ]  # 20
    FOUNDATION_D = [
        "Wooden piles",
        "Trepeler",
        "D - På løsmasser - Punkt- og trefundamenter (banketter)",
    ]  # 50

    STRUCTURE_A = ["Steel", "A - Stål"]
    STRUCTURE_B = ["Reinforced concrete", "B - Armert betong"]
    STRUCTURE_C = ["Mixed", "C - Tre eller varierende"]
    STRUCTURE_D = ["Masonry", "D - Murstein eller spesiell type"]

    STATUS_A = ["Excellent", "A - Meget god tilstand"]
    STATUS_B = ["Good", "B - God tilstand"]
    STATUS_C = ["Medium", "C - Brukbar tilstand"]
    STATUS_D = ["Bad", "D - Dårlig"]

    VULN_WEIGHT_LENGTH = 0.75
    VULN_WEIGHT_SHAPE = 0.75
    VULN_WEIGHT_STRUCTURE = 1.0
    VULN_WEIGHT_FOUNDATION = 1.5
    VULN_WEIGHT_STATUS = 1.0

    WALL_CORNER_ANGLE_THRESHOLD = 5  # degrees
    CONSTR_RESAMPLE_LEN = 2.0

    logger.debug(
        f"### mainBegrenSkade called with: \nbuildingsFN: {buildingsFN}, tunnelJson: {tunnelJson},output_ws: {output_ws}, "
        f"feature_name: {feature_name}, coord_syst: {output_proj}, bShortterm: {bShortterm}, "
        f"tunnel_depth: {tunnel_depth}, tunnel_diameter: {tunnel_diameter}, volume_loss: {volume_loss}, trough_width: {trough_width},  bLongterm: {bLongterm},"
        f"tunnel_leakage: {tunnel_leakage}, porewp_calc_type: {porewp_calc_type}, porewp_red_at_site_m: {porewp_red_at_site_m}, dtb_raster: {dtb_raster}, "
        f"dry_crust_thk: {dry_crust_thk}, density_sat: {density_sat},OCR: {OCR}, janbu_ref_stress: {janbu_ref_stress}, janbu_const: {janbu_const}, "
        f"janbu_m: {janbu_m}, consolidation_time: {consolidation_time}, bVulnerability: {bVulnerability}, fieldNameFoundation: {fieldNameFoundation}, "
        f"fieldNameStructure: {fieldNameStructure}, fieldNameStatus: {fieldNameStatus}"
    )

    # _array(corner)
    logger.debug("Calling get_construction_corners_from_ArcGIS_json with json: {}".format(tunnelJson))
    construction_area_corners = BegrensSkadeLib.get_construction_corners_from_ArcGIS_json(tunnelJson, CONSTR_RESAMPLE_LEN, logger)

    # array(buildings)
    logger.info("TIME - getting buildings from shapefile")

    if bLongterm:
        buildings = BegrensSkadeLib.get_buildings_with_dtb(
            buildingsFN,
            dtb_raster,
            fieldNameFoundation=fieldNameFoundation,
            fieldNameStructure=fieldNameStructure,
            fieldNameStatus=fieldNameStatus,
            logger=logger,
        )
    else:
        buildings = BegrensSkadeLib.get_buildings(
            buildingsFN,
            fieldNameFoundation=fieldNameFoundation,
            fieldNameStructure=fieldNameStructure,
            fieldNameStatus=fieldNameStatus,
            logger=logger,
        )

    logger.info("Gotten buildings from shapefile")

    count_adj = 0
    logger.info("Calculate settlements")

    buildingsWithBedrock = []

    for building in buildings:

        # filter corner duplicates and straight-wall corners
        building.filter_duplicates()
        building.filter_straights(WALL_CORNER_ANGLE_THRESHOLD)

        # if bLongterm, add bedrock depth information
        if bLongterm:
            building.corners = BegrensSkadeLib.appendZValuesFromRaster(building.corners, dtb_raster, logger)
            if building.corners is None:
                logger.debug("Skipping building")
                continue
            else:
                buildingsWithBedrock.append(building)

        max_sv_short = 0.0
        max_sh_short = 0.0
        max_sv_total = 0.0
        near_closest = 9999.0
        near_furthest = 0.0

        for corner in building.corners:

            # Evaluating distance from every building corner to construction zone...
            near_dist, near_angle = BegrensSkadeLib.near_analysis(corner.x, corner.y, construction_area_corners)

            corner.near_angle = near_angle
            corner.near_dist = near_dist

            near_closest = min(near_dist, near_closest)
            near_furthest = max(near_dist, near_furthest)

            porewp_red_atdist = -1

            # Evaluating settlements at corners

            if bShortterm:
                sv_short = BegrensSkadeLib.get_sv_short_Peck(near_dist, tunnel_depth, tunnel_diameter, volume_loss, trough_width)
            else:
                sv_short = 0

            if bLongterm:

                if porewp_calc_type == "Øvre":
                    porewp_red_atdist = 10 * ((float(tunnel_leakage) - 4.17) / 4.54 - 0.02 * near_dist)
                elif porewp_calc_type == "Typisk":
                    porewp_red_atdist = 10 * ( (float(tunnel_leakage) - 2) / 3.11 - 0.02 * near_dist)
                elif porewp_calc_type == "Nedre":
                    porewp_red_atdist = 10 * ((float(tunnel_leakage) - 1.61) / 1.59 - 0.02 * near_dist)
                elif porewp_calc_type == "Manuell":
                    porewp_red_atdist = 10 * max((porewp_red_at_site_m - (0.02 * near_dist)), 0) # pore water pressure reduction in kPa
                else:
                    raise Exception(
                        "Not valid parameter",
                        "Not valid porewater pressure calculation type",
                    )

                sv_long, red_adj = BegrensSkadeLib.get_sv_long_janbu(
                    corner.dtb,
                    dry_crust_thk,
                    dep_groundwater,
                    density_sat,
                    OCR,
                    porewp_red_atdist,
                    janbu_ref_stress,
                    janbu_const,
                    janbu_m,
                    consolidation_time,
                )
                if red_adj:
                    count_adj += 1
            else:
                sv_long = 0

            corner.sv_short = sv_short
            corner.sh_short = 0
            corner.sv_long = sv_long
            corner.sv_tot = sv_short + sv_long
            corner.porewp_red = porewp_red_atdist

            max_sv_short = max(max_sv_short, sv_short)
            max_sh_short = 0
            max_sv_total = max(max_sv_total, (sv_short + sv_long))

        building.create_walls()

        max_angle = 0.0
        max_strain = 0.0
        max_principal_strain = 0.0

        # logger.debug("Number of walls:{}".format(len(building.walls)))
        for wall in building.walls:
            max_angle = max(max_angle, abs(wall.slope_angle))
            max_strain = max(max_strain, abs(wall.hor_tensile_strain))
            max_principal_strain = max(
                max_principal_strain, abs(wall.principal_tensile_strain)
            )

        building.max_sv_short = max_sv_short
        building.max_sh_short = max_sh_short
        building.max_sv_total = max_sv_total
        building.max_angle = max_angle
        building.max_strain = max_strain
        building.max_principal_strain = max_principal_strain

        # Evaluate vulnerability
        if bVulnerability:

            logger.debug("Calculating vulnerability for building: {}".format(building.bid))
            if (fieldNameFoundation):

                foundation = building.foundation

                if foundation in FOUNDATION_A:
                    foundation_cvi = 0
                elif foundation in FOUNDATION_B:
                    foundation_cvi = 5
                elif foundation in FOUNDATION_C:
                    foundation_cvi = 20
                elif foundation in FOUNDATION_D:
                    foundation_cvi = 50
                else:
                    foundation_cvi = 50  # default
                building.foundation = foundation_cvi


            if fieldNameStructure:
                structure = building.structure

                if structure in STRUCTURE_A:
                    structure_cvi = 0
                elif structure in STRUCTURE_B:
                    structure_cvi = 5
                elif structure in STRUCTURE_C:
                    structure_cvi = 20
                elif structure in STRUCTURE_D:
                    structure_cvi = 50
                else:
                    structure_cvi = 50  # default

                building.structure = structure_cvi


            if fieldNameStatus:
                status = building.status

                if status in STATUS_A:
                    status_cvi = 0
                elif status in STATUS_B:
                    status_cvi = 5
                elif status in STATUS_C:
                    status_cvi = 20
                elif status in STATUS_D:
                    status_cvi = 50
                else:
                    status_cvi = 50  # default

                building.status = status_cvi

            isosquare = 16 * building.area / (building.circumf ** 2)
            radial_buil_len = near_furthest - near_closest
            logger.debug(
                "building.area: "
                + str(building.area)
                + ", building.circumf: "
                + str(building.circumf)
                + ", isosquare "
                + str(isosquare)
            )
            logger.debug(
                "near_furthest: "
                + str(near_furthest)
                + ", near_closest: "
                + str(near_closest)
                + ", radial_buil_len: "
                + str(radial_buil_len)
            )

            length_cvi = BegrensSkadeLib.get_buil_len_cvi(radial_buil_len)
            shape_cvi = BegrensSkadeLib.get_buil_shape_cvi(isosquare)

            try:
                vulnerability = ( VULN_WEIGHT_LENGTH * length_cvi + VULN_WEIGHT_SHAPE * shape_cvi)
                maxsum = 50.0 * (VULN_WEIGHT_LENGTH + VULN_WEIGHT_SHAPE)

                if fieldNameStructure:
                    vulnerability += structure_cvi * VULN_WEIGHT_STRUCTURE
                    maxsum += 50.0 * VULN_WEIGHT_STRUCTURE
                if fieldNameFoundation:
                    vulnerability += foundation_cvi * VULN_WEIGHT_FOUNDATION
                    maxsum += 50.0 * VULN_WEIGHT_FOUNDATION
                if fieldNameStatus:
                    vulnerability += status_cvi * VULN_WEIGHT_STATUS
                    maxsum += 50.0 * VULN_WEIGHT_STATUS

                vulnerability = vulnerability / maxsum
                building.vulnerability = vulnerability
                vuln_cvi = BegrensSkadeLib.get_buil_vuln_cvi(vulnerability)

                impact_angle_cvi = BegrensSkadeLib.get_buil_impact_angle_cvi(max_angle)
                impact_totset_cvi = BegrensSkadeLib.get_buil_impact_totset_cvi(max_sv_total)

                # logger.debug("Calculating risk_totset and risk_angle for building: {}".format(building.bid))

                building.risk_totset = BegrensSkadeLib.get_risk_cvi(vuln_cvi, impact_totset_cvi)
                building.risk_angle = BegrensSkadeLib.get_risk_cvi(vuln_cvi, impact_angle_cvi)
            except Exception as e:
                logger.debug("Error, type; {0}, error args: {1}, error: {2}".format(type(e), e.args, e))
                # logger.error(e.args)
                # logger.error(e)

    now = datetime.now()  # current date and time
    date_time_str = now.strftime("_%Y%m%d_%H%M%S")

    corner_name = feature_name + date_time_str + "_C_"
    wall_name = feature_name + date_time_str + "_W_"
    building_name = feature_name + date_time_str + "_B_"
    if bShortterm:
        corner_name += "S"
        wall_name += "S"
        building_name += "S"
    if bLongterm:
        corner_name += "L"
        wall_name += "L"
        building_name += "L"

    logger.info("TIME - done calculating settlements")

    filterValue = feature_name

    # logger.debug("Delete threads returned")
    logger.debug(f"Number of buildings {len(buildings)}, Number of buildings with bedrock: {len(buildingsWithBedrock)}")

    if bLongterm:
        buildings = buildingsWithBedrock

    logger.info("TIME - writing results to shape")
    building_shapefile = output_ws + os.sep + building_name + ".shp"
    #building_shapefile_prj = output_ws + os.sep + building_name + "_prj.shp"
    BegrensSkadeLib.writeBuildingsToShape(building_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(building_shapefile,building_shapefile_prj,str(working_proj), str(output_proj), "polygon")

    wall_shapefile = output_ws + os.sep + wall_name + ".shp"
    #wall_shapefile_prj = output_ws + os.sep + wall_name + "_prj.shp"
    logger.debug("wall_shapefile:{}".format(wall_shapefile))
    BegrensSkadeLib.writeWallsToShape(wall_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(wall_shapefile,wall_shapefile_prj,str(working_proj), str(output_proj), "line")

    corner_shapefile = output_ws + os.sep + corner_name + ".shp"
    #corner_shapefile_prj = output_ws + os.sep + corner_name + "_prj.shp"
    logger.debug("corner_shapefile:{}".format(corner_shapefile))
    BegrensSkadeLib.writeCornersToShape(corner_shapefile, buildings, output_proj, filterValue, logger)
    #Utils.projectLayer(corner_shapefile,corner_shapefile_prj,str(working_proj), str(output_proj),"point")
    logger.info("TIME - written all shapefiles")

    return [building_shapefile, wall_shapefile, corner_shapefile]

