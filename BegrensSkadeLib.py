import numpy as np
from osgeo import gdal
from osgeo import ogr
import logging.handlers
# import Logger
import Utils

'''
This library contains classes and related functions that are specific to the BegrensSkade GIBV programs,
empirical curves for short and long term settlements, 
geotechnical calculators for settlements in time and space,
and vulnerability index functions.
'''

# For debugging or quality control of stepwise, detailed geotechnical calulations,
# the user has the possibility to enable output of a csv log containing such details:
write_csv = False

################ CLASSES FOR CORNER, WALLS AND BUILDINGS ##################

class Corner:
    def __init__(self, cid, x, y):
        self.cid = cid  # corner id to keep track on the line direction
        self.x = x
        self.y = y
        self.dtb = None
        self.sv_short = None
        self.sv_long = None
        self.sh_short = None
        self.near_dist = None
        self.near_angle = None
        self.porewp_red = None

class Wall:
    def __init__(self, wid, corner_start, corner_end, slope_angle=None, hor_tensile_strain=None, principal_tensile_strain=None):
        self.wid = wid
        self.corner_start = corner_start
        self.corner_end = corner_end
        self.slope_angle = slope_angle
        self.hor_tensile_strain = hor_tensile_strain
        self.principal_tensile_strain = principal_tensile_strain

class Building:
    def __init__(self, bid, corners, area=None, circumf=None, foundation=None, structure=None, status=None, logger=None):
        self.bid = bid
        self.corners = corners  # array[corner]
        self.walls = []  # array[wall]

        self.area = area
        self.circumf = circumf
        self.foundation = foundation
        self.structure = structure
        self.status = status

        self.max_sh_sho = None
        self.max_sv_sho = None
        self.max_sv_tot = None
        self.max_angle = None
        self.max_strain = None
        self.max_pstrai = None

        self.vulnerability = None
        self.risk_totset = None
        self.risk_angle = None
        self.logger = logger

    def filter_duplicates(self):

        security = 0
        i = 0
        while i < len(self.corners):
            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y
            dist = np.sqrt((y - y_prev) ** 2 + (x - x_prev) ** 2)
            if dist < 0.1:
                # Removing duplicate
                del self.corners[i - 1]
            else:
                i += 1

            if security > 10000:
                self.logger.error("While overrun in filter duplicates")
                raise Exception("While overrun in filter duplicates")
            security += 1

    def filter_straights(self, WALL_CORNER_ANGLE_THRESHOLD):
        # remove points that are along a straigth wall
        security = 0
        # i starts at 0, that means, first prev point has index -1, in this way first and last point are compared.
        i = 0
        while i < (len(self.corners) - 1):
            oid = self.corners[i].cid
            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y
            x_next = self.corners[i + 1].x
            y_next = self.corners[i + 1].y
            try:
                prev_dir, prev_quadr = Utils.get_angle(x_prev, y_prev, x, y)
                next_dir, next_quadr = Utils.get_angle(x, y, x_next, y_next)
            except Exception("Angle calculation failed, skipping building " + str(self.bid)):
                break

            prev_dir_deg = 180 * prev_dir / np.pi
            next_dir_deg = 180 * next_dir / np.pi

            if abs(prev_dir_deg - next_dir_deg) < WALL_CORNER_ANGLE_THRESHOLD:  # degree
                # Removing straight pnt
                del self.corners[i]
            else:
                i += 1

            if security > 10000:
                self.logger.error("While overrun in filter straights")
                raise Exception("While overrun in filter straights")
            security += 1

    def create_walls(self):
        cid = 0

        for i in range(0, len(self.corners)):

            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y

            xy_dist = np.sqrt((x - x_prev) ** 2 + (y - y_prev) ** 2)

            sv_short_prev = self.corners[i - 1].sv_short
            sv_short = self.corners[i].sv_short

            sh_short_prev = self.corners[i - 1].sh_short
            sh_short = self.corners[i].sh_short

            sv_long_prev = self.corners[i - 1].sv_long
            sv_long = self.corners[i].sv_long

            sv_short_term_diff = abs(sv_short_prev - sv_short)
            sh_short_term_diff = abs(sh_short_prev - sh_short)
            sv_total_diff = abs(
                (sv_short_prev + sv_long_prev) - (sv_short + sv_long))

            hor_tensile_strain = sh_short_term_diff / xy_dist

            if xy_dist < 1:
                slope_angle = 0
            else:
                slope_angle = sv_total_diff / xy_dist

            try:
                thetamax = 0.5 * np.arctan(slope_angle / hor_tensile_strain)
            except:
                thetamax = None

            try:
                principal_tensile_strain = hor_tensile_strain * (np.cos(thetamax)) ** 2 + slope_angle * np.sin(thetamax) * np.cos(thetamax)
            except:
                principal_tensile_strain = None

            self.walls.append(Wall(cid,
                                   self.corners[i-1], self.corners[i], slope_angle, hor_tensile_strain, principal_tensile_strain))

            cid += 1


####################### SHORT TERM SETTLEMENT CURVES #######################

def get_sv_short_a(near_dist, excavation_depth):
    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 3

    if x < 0.3:
        return [(0.15 + ((0.5-0.15)/0.3)*x)*excavation_depth*0.01, W]

    elif x >= 0.3 and x < 1:
        return [(0.5 + ((0.1-0.5)/(1-0.3))*(x-0.3))*excavation_depth*0.01 , W]

    elif x >= 1 and x < 1.5:
        return [(0.1 + ((-0.1)/(1.5-1))*(x-1)) *excavation_depth*0.01, W]

    else:
        return [0.0, W]

def get_sv_short_b(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 3

    if x < 0.5:
        return [(0.3 + ((1-0.3)/0.5)*x)*excavation_depth*0.01, W]

    elif x >= 0.5 and x < 2:
        return [(1 + ((0.2-1)/(2-0.5))*(x-0.5))*excavation_depth*0.01 , W]

    elif x >= 2 and x < 3:
        return [(0.2 + ((-0.2)/(3-2))*(x-2))*excavation_depth*0.01,W]

    else:
        return [0.0, W]

def get_sv_short_c(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 5

    if x < 0.7:
        return [(1 + ((2-1)/0.7)*x)*excavation_depth*0.01, W]

    elif x >= 0.7 and x < 2.5:
        return [(2 + ((0.5-2)/(2.5-0.7))*(x-0.7))*excavation_depth*0.01 , W]

    elif x >= 2.5 and x < 4:
        return [(0.5 + ((-0.5)/(4-2.5))*(x-2.5)) *excavation_depth*0.01, W]

    else:
        return [0.0, W]

def get_sv_short_d(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 5

    if x < 1:
        return [(1.5 + ((3 - 1.5) / 1) * x)*excavation_depth*0.01, W]

    elif x >= 1 and x < 3:
        return [(3 + ((0.75 - 3) / (3 - 1)) * (x - 1))*excavation_depth*0.01, W]

    elif x >= 3 and x < 5:
        return [(0.75 + ((-0.75) / (5 - 3)) * (x - 3))*excavation_depth*0.01, W]

    else:
        return [0.0, W]


################## PECK SETTLEMENT CURVE FOR SOIL TUNNELS ##################

def get_sv_short_Peck(near_dist, tunnel_depth, tunnel_diameter, volume_loss, trough_width):
    D = tunnel_diameter
    i = trough_width * tunnel_depth
    s0 = volume_loss/100*((np.pi*D**2)/4)/(np.sqrt(2*np.pi)*i) 
    s = s0 * np.exp(-(near_dist**2)/(2*i**2)) 
    return s 


################## JANBU LONG TERM SETTLEMENT CALCULATION #################

def get_sv_long_janbu(dtb, dry_crust_thk, dep_groundwater, density_sat, OCR, porewp_red, p_ref, janbu_const, janbu_m, consolidation_time):
    density_water = 10  # kN/m3
    permeability = 1e-9  # (m/s)
    adj = False
    dep_to_clay = max(dry_crust_thk, dep_groundwater)
    clay_thk = dtb - dep_to_clay
    sv_acc = 0  # accumulated horizontal settlement in clay layer

    janbu_M_avgs = []  #depth average of janbu_M to be used in terzagi time equation
    METHOD = "LYR"
    #METHOD = "AVG"

    for clay_dep_i in range(int(round(dep_to_clay, 0)), int(round(dtb, 0))):

        clay_dep_i += 0.5
        pz0 = density_sat * clay_dep_i
        uz0 = density_water * (clay_dep_i - dep_groundwater)
        pz0_eff = pz0 - uz0
        pz_pre = OCR * pz0_eff

        if pz_pre < p_ref:
            #This can happen right below the surface where pz0_eff is low, implies negative stiffness and does not make sense
            raise Exception("Negative stiffness in Janbu model. Try to decrease the reference pressure or increase the OCR")

        #porewp_red = porewp_red * longterm_porewr

        uz0_tot = density_water * (dtb - dep_groundwater)
        if porewp_red > uz0_tot:
            # arcpy.AddMessage("porewp_red_atdist: " + str(porewp_red_atdist))
            # arcpy.AddMessage("density_water * (dtb- dep_groundwater): " + str(density_water * (dtb - dep_groundwater)))
            adj = True
            porewp_red = density_water * (dtb - dep_groundwater)

        clay_dep_lyr = (clay_dep_i - dep_to_clay)

        du = (porewp_red/clay_thk)*clay_dep_lyr
        janbu_M = janbu_const * janbu_m * pz_pre
        janbu_M_avgs.append(janbu_M)

        if pz0_eff + du < pz_pre:
            # implicitly multiplied by 1 m clay thickness (iteration step)
            sv = du / janbu_M
            # regime = "flat"

        else:
            # implicitly multiplied by 1 m clay thickness (iteration step)
            sv = (pz_pre - pz0_eff)/janbu_M + (1/janbu_m) * \
                np.log((pz0_eff + du - p_ref)/(pz_pre - p_ref))
            # regime = "linear"

        # Calculation of time-dependent consolidation. Default is maxiumn consolidation (1000 years - gives 99.9 % consolidation for 50 m clay)
        if METHOD == "LYR":
            t = 60 * 60 * 24 * 365 * consolidation_time  # sec
            c = janbu_M * permeability / density_water
            T = c * t / (clay_thk) ** 2
            sv_inf = sv
            sv = sv * U_draintop_b(T, 10)

        sv_acc += sv




    #Calculation of time-dependent consolidation. Default is maxiumn consolidation (1000 years - gives 99.9 % consolidation for 50 m clay)
    if METHOD == "AVG" and sv_acc > 0:
        janbu_M_avg = sum(janbu_M_avgs)/len(janbu_M_avgs)
        t = 60 * 60 * 24 * 365 * consolidation_time  # sec
        c = janbu_M_avg * permeability / density_water
        T = c * t / (clay_thk) ** 2
        sv_acc = sv_acc * U_draintop_b(T, 10)

    return sv_acc, adj


################## LONG TERM POREWATER REDUCTION CURVES ###################

def get_longterm_porewr_min(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 220:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 25
        K = 0.62
        L = 0.6
        M = 0.96

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0

def get_longterm_porewr_max(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 360:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 26
        K = 1
        L = 0.625
        M = 0.985

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0

def get_longterm_porewr_mean(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 340:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 30
        K = 0.8
        L = 0.68
        M = 0.985

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0


#################### TERZAGI CONSOLIDATION TIME CURVES ####################

def U_drainboth(T, m_max):

    sum = 0

    for m in range(0, m_max):
        M = 0.5 * np.pi * (2*m + 1)

        sum += (2/(M**2))*np.exp(-T*M**2)

    return 1 - sum

def U_draintop_a(T, m_max): #Trekant med spiss oppad

    sum = 0

    for m in range(0, m_max):

        M = 0.5 * np.pi * (2*m - 1)

        sum += 2*((-1)**(m+1))*np.exp(-T*M**2)/M**3

    return 1 - sum

def U_draintop_b(T, m_max): #Trekant med spiss nedad

    return 2*U_drainboth(T, m_max) - U_draintop_a(T, m_max)


######################## VULNERABILITY INDEX FUNCTIONS ####################

def get_buil_len_cvi(buil_len):
    if buil_len < 10:
        return 0
    elif buil_len < 15:
        return 5
    elif buil_len < 30:
        return 20
    else:
        return 50

def get_buil_shape_cvi(squareness):
    if squareness < 0.35:
        return 50
    elif squareness < 0.5:
        return 20
    elif squareness < 0.75:
        return 5
    else:
        return 0

def get_buil_impact_totset_cvi(tot_set):
    if tot_set < 0.010:
        return 1
    elif tot_set < 0.050:
        return 2
    elif tot_set < 0.075:
        return 3
    else:
        return 4

def get_buil_impact_angle_cvi(angle):
    if angle < (1/500):
        return 1
    elif angle < (1/200):
        return 2
    elif angle < (1/50):
        return 3
    else:
        return 4

def get_buil_vuln_cvi(vuln):
    if vuln < 0.25:
        return 1
    elif vuln < 0.5:
        return 2
    elif vuln < 0.75:
        return 3
    else:
        return 4

def get_risk_cvi(vuln_cvi, impact_cvi):
    if vuln_cvi * impact_cvi < 3:
        return 1
    elif vuln_cvi * impact_cvi < 5 and vuln_cvi < 4 and impact_cvi < 4:
        return 2
    elif vuln_cvi * impact_cvi < 8:
        return 3
    elif vuln_cvi * impact_cvi < 12:
        return 4
    else:
        return 5


################ GIBV SPECIFIC SHAPEFILE I/O FUNCTIONALITY #################

def get_construction_corners(shapeFile, CONSTR_RESAMPLE_LEN, logger):
    construction_area_corner_pnts = []
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile, 0)
    layer = dataSource.GetLayer()
    cid = 0

    for feature in layer:
        # gjson = feature.ExportToJson()
        # pj = json.loads(gjson)
        # geom = feature.geometry()

        cid = 1

        for part in feature.geometry():
            # prev_X = part[0].X
            # prev_Y = part[0].Y
            first = True
            for p in part.GetPoints():
                X = p[0]
                Y = p[1]
                if first:
                    prev_X = X
                    prev_Y = Y
                    first = False
                    continue

                # print("X: {}  Y:{}".format(x,y))

                wall_len = np.sqrt((X-prev_X)**2 + (Y-prev_Y)**2)
                if wall_len > CONSTR_RESAMPLE_LEN:
                    n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
                    seg_len = wall_len / n_segments
                    for seg in range(0, n_segments):
                        X_sampl = prev_X + seg * seg_len * \
                            (X - prev_X) / wall_len
                        Y_sampl = prev_Y + seg * seg_len * \
                            (Y - prev_Y) / wall_len

                        construction_area_corner_pnts.append(
                            Corner(cid, X_sampl, Y_sampl))
                else:
                    construction_area_corner_pnts.append(Corner(cid, X, Y))
                cid += 1
                prev_X = X
                prev_Y = Y
    return construction_area_corner_pnts

def get_buildings(features, fieldNameFoundation=None, fieldNameStructure=None, fieldNameStatus=None, logger=None) -> []:
    input_buildings = []
    bid = 1
    numBuildings = 0
    numSkippedBuildings = 0

    if (not logger is None):
        logger.debug("Opening shapefile: {}".format(features))

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = shp_driver.Open(features, 1)
    layer = dataSource.GetLayer()
    foundation = None
    structure = None
    status = None

    for buildingFeature in layer:
        numBuildings +=1
        cid = 0

        g = buildingFeature.geometry().GetGeometryRef(0)

        if (g is None):
            # logger.debug("Geometry is None, jumping to next")
            continue
        if (g.GetPoints() is None):
            # logger.debug("g.GetPoints() is None, jumping to next")
            continue

        corners = []

        for pnt in g.GetPoints():
            corners.append(Corner(cid, pnt[0], pnt[1]))
            cid += 1

        if (not fieldNameFoundation is None):
            foundation = buildingFeature.GetField(fieldNameFoundation)
        if (not fieldNameStructure is None):
            structure = buildingFeature.GetField(fieldNameStructure)
        if (not fieldNameStatus is None):
            status = buildingFeature.GetField(fieldNameStatus)

        input_buildings.append(
            Building(bid, corners, g.GetArea(), g.Length(),
                     foundation=foundation, structure=structure, status=status, logger=logger))
        #logger.debug(f"Building {bid}, num corners: {len(corners)}, area: {g.GetArea()}, length: {g.Length()} - get_buildings")
        bid += 1
    return input_buildings

def get_buildings_with_dtb(features, dtb_filename, fieldNameFoundation=None, fieldNameStructure=None, fieldNameStatus=None, logger=None) -> []:
    input_buildings = []
    bid = 1
    numBuildings = len(input_buildings)
    numSkippedBuildings = 0

    if (not logger is None):
        logger.debug("Opening shapefile: {}".format(features))

    src_ds = gdal.Open(str(dtb_filename))
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)
    no_data_value = rb.GetNoDataValue()
    logger.debug(f"No data value for raster = {no_data_value}")

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = shp_driver.Open(features, 1)
    layer = dataSource.GetLayer()
    foundation = None
    structure = None
    status = None
    
    min_x = gt[0]
    max_y = gt[3]
    max_x = min_x + (gt[1] * src_ds.RasterXSize)
    min_y = max_y + (gt[5] * src_ds.RasterYSize)  # gt[5] is usually negative
    
    # Counter for points outside raster bounds
    numPointsOutsideRaster = 0


    for buildingFeature in layer:

        cid = 0

        g = buildingFeature.geometry().GetGeometryRef(0)

        if (g is None):
            # logger.debug("Geometry is None, jumping to next")
            continue
        if (g.GetPoints() is None):
            # logger.debug("g.GetPoints() is None, jumping to next")
            continue

        corners = []
        foundNullVal = False


        for pnt in g.GetPoints():
            if pnt:
                x, y = pnt[0], pnt[1]

                # Skip points outside the raster bounds
                if x < min_x or x > max_x or y < min_y or y > max_y:
                    numPointsOutsideRaster += 1
                    continue
                
                px = int((pnt[0] - gt[0]) / gt[1])  # x pixel
                py = int((pnt[1] - gt[3]) / gt[5])  # y pixel
                dtb = rb.ReadAsArray(px, py, 1, 1)                
                if dtb is None or dtb == no_data_value:
                    numSkippedBuildings += 1
                    foundNullVal = True
                    #if one of the corners has no dbt value, skip to next building
                    #logger.debug("Skipping buildingNumber of skipped buildings: {}".format(numSkippedBuildings))
                    break
                corners.append(Corner(cid, pnt[0], pnt[1]))              
                cid += 1

        if foundNullVal:
            continue
        if (not fieldNameFoundation is None):
            foundation = buildingFeature.GetField(fieldNameFoundation)
        if (not fieldNameStructure is None):
            structure = buildingFeature.GetField(fieldNameStructure)
        if (not fieldNameStatus is None):
            status = buildingFeature.GetField(fieldNameStatus)

        input_buildings.append(
            Building(bid, corners, g.GetArea(), g.Length(),
                     foundation=foundation, structure=structure, status=status, logger=logger))
        bid += 1
    # Log the total number of skipped points **only once**
    if numPointsOutsideRaster > 0:
        logger.debug(f"Total points skipped due to falling outside raster bounds: {numPointsOutsideRaster}")

    logger.debug(f"Number of skipped buildings: {numSkippedBuildings} of {numBuildings} - get_buildings_with_dtb")
    return input_buildings

def get_construction_corners_from_ArcGIS_json(jsonBody, CONSTR_RESAMPLE_LEN, logger):
    construction_area_corner_pnts = []

    # logger.debug("json: {}".format(jsonBody))
    if "features" in jsonBody:
        logger.debug("FINNES")
        jsonBody = jsonBody["features"]
        if "geometry" in jsonBody[0]:
            logger.debug("features FINNES")
            jsonBody = jsonBody[0]["geometry"]
        else:
            logger.debug("features FINNES IKKE")
    else:
        logger.debug("Cant fint fatures in json body, getting rings directly")

    # logger.debug("json etter avskrelling: {}".format(jsonBody))
    # features = jsonBody["features"]
    # geometry = features[0]["geometry"]
    rings = jsonBody["rings"]
    # logger.debug("rings: {}".format(rings))
    for ring in rings:
        # gjson = feature.ExportToJson()
        # pj = json.loads(gjson)
        # geom = feature.geometry()
        logger.debug("Getting corners in excavation")

        cid = 1
        first = True
        for p in ring:
            # logger.debug("point: {}".format(p))
            # for p in part.GetPoints():
            X = p[0]
            Y = p[1]
            if first:
                prev_X = X
                prev_Y = Y
                first = False
                continue

            # print("X: {}  Y:{}".format(x,y))

            wall_len = np.sqrt((X-prev_X)**2 + (Y-prev_Y)**2)
            if wall_len > CONSTR_RESAMPLE_LEN:
                n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
                seg_len = wall_len / n_segments
                for seg in range(0, n_segments):
                    X_sampl = prev_X + seg * seg_len * \
                        (X - prev_X) / wall_len
                    Y_sampl = prev_Y + seg * seg_len * \
                        (Y - prev_Y) / wall_len

                    construction_area_corner_pnts.append(
                        Corner(cid, X_sampl, Y_sampl))
            else:
                construction_area_corner_pnts.append(Corner(cid, X, Y))
            cid += 1
            prev_X = X
            prev_Y = Y
    # logger.debug("construction_area_corner_pnts: {}".format(
    #    construction_area_corner_pnts))
    return construction_area_corner_pnts

def createBuildingCornersDict(shapeFile, fieldNames, bLongterm, short_term_curve, byggegrop_depth, hor_vert_ratio, dry_crust_thk, dep_groundwater,
                              density_sat, OCR, porewp_red, janbu_ref_stress, janbu_m, pw_reduction_curve, logger) -> {}:

    bid_prev = "start"
    pnt_array = []
    result_building_corner_dict = {}

    first_time = True

    # Loop through the newly created point feature that contains the building corners.
    # Sh and sv and sv_long are calculated for every points.
    # Store corner information in dictionary building_corners_dict.
    # The key is the building_id (normally object ID), the value is a table containing the corner points for that particular building.
    # The purpose of this is to have a faster accessible data structure for the corner points, that can be edited/fitered for duplicates and straigt-wall points.
    count_all = 0
    count_adj = 0

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile, 0)
    layer = dataSource.GetLayer()

    for pnt in layer:
        oid = pnt.GetField(0)
        bid = pnt.GetField(1)
        x = pnt.GetField(2)
        y = pnt.GetField(3)
        near_dist = pnt.GetField(4)

        if bLongterm:
            dtb = pnt.GetField(6)  # 5??

        # arcpy.AddMessage(str(oid) + ", " + str(bid) + ", " + str(x) + ", " + str(y) + ", " + str(near_dist) + ", " + str(angle))


        if short_term_curve in ["Wall not to bedrock - regression", "0,5 % av byggegropdybde", "Spunt installert til berg med høy sikkerhet", "Norm_setning_0.5"]:
            sv_short, W = get_sv_short_a(near_dist, byggegrop_depth)
        elif short_term_curve in ["Wall not to bedrock - discrete", "1 % av byggegropdybde", "Spunt installert til berg med lav sikkerhet", "Norm_setning_1"]:
            sv_short, W = get_sv_short_b(near_dist, byggegrop_depth)
        elif short_term_curve in ["Tie-back anchors - regression", "2 % av byggegropdybde", "Svevespunt høy sikkerhet", "Norm_setning_2"]:
            sv_short, W = get_sv_short_c(near_dist, byggegrop_depth)
        elif short_term_curve in ["Tie-back anchors - regression", "3 % av byggegropdybde", "Svevespunt lav sikkerhet", "Norm_setning_3"]:
            sv_short, W = get_sv_short_d(near_dist, byggegrop_depth)
        else:
            raise Exception("Not a valid regression curve: " + str(short_term_curve) )

        norm_dist = near_dist / byggegrop_depth  # normalized distance from byggegrop

        sh_short = -hor_vert_ratio * (1 + 2*norm_dist/float(W))*sv_short

        if bLongterm:
            sv_long, red_adj = get_sv_long_janbu(
                 dtb, dry_crust_thk, dep_groundwater, density_sat, OCR, porewp_red, janbu_ref_stress, janbu_m)
            if red_adj:
                count_adj += 1
        else:
            sv_long = 0

        if bid_prev == "start":
            bid_prev = bid

        if bid != bid_prev:
            # New building. Store line data.
            # arcpy.AddMessage("Adding new building: " + str(bid_prev) + ", " + str(len(pnt_array)))
            result_building_corner_dict[bid_prev] = pnt_array
            pnt_array = []
            # arcpy.AddMessage([oid, x, y, near_dist, sv_short, sh_short, sv_long])
            pnt_array.append(
                [oid, x, y, near_dist, sv_short, sh_short, sv_long])

        else:
            pnt_array.append(
                [oid, x, y, near_dist, sv_short, sh_short, sv_long])

        bid_prev = bid
        count_all += 1

    if bLongterm:
        logger.info("Percentage adjustment of porewp red. due to shallow bedrock : " +
                    str(round(100*(count_adj/count_all), 1)))

    # write the last building
    result_building_corner_dict[bid_prev] = pnt_array
    return result_building_corner_dict

def near_analysis_sqr(xref, yref, construction_area_corners):

    near_dist_sqr = 999999

    for corner in construction_area_corners:
        x = corner.x
        y = corner.y

        dist_sqr = (xref-x)**2+(yref-y)**2

        if dist_sqr < near_dist_sqr:
            near_dist_sqr = dist_sqr

    return near_dist_sqr

def near_analysis(xref, yref, construction_area_corners):

    near_dist = 999999
    near_angle = 0  # unit vector

    for corner in construction_area_corners:
        x = corner.x
        y = corner.y

        dist = np.sqrt((xref-x)**2+(yref-y)**2)
        angle = Utils.getAngleFromDir(x, y, xref, yref)

        if dist < near_dist:
            near_dist = dist
            near_angle = angle

    return near_dist, near_angle

def writeBuildingsToShape(buildingsFN, buildings, projection, filterValue, logger):
    fields = [["bid", "long"], ["circumf", "float"], ["foundation", "float"], ["structure", "float"], ["status", "float"],
              ["max_sv_shr", "float"], ["max_sh_shr", "float"], [
                  "max_sv_tot", "float"], ["max_angle", "float"],
              ["max_strain", "float"], ["max_p_stra", "float"], ["vulnerab", "float"], ["risk_tots", "float"], ["risk_angle", "float"], ["filter_val", "string"]]

    try:
        numBuildings = 0
        Utils.createShapefile(buildingsFN, "polygon", projection, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(buildingsFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()

        for building in buildings:
            corners = []
            for corner in building.corners:
                c = [corner.x, corner.y]
                corners.append(c)
            geom = Utils.createPolygon(corners)
            # logger.debug(
            #     "writeBuildingsToShape, polygon created:{}".format(geom))
            outFeature = ogr.Feature(featureDefn)

            # logger.debug(
            #     "writeBuildingsToShape, feature created:{}".format(outFeature))
            outFeature.SetGeometry(geom)
            # logger.debug(
            #     "writeBuildingsToShape, feature created and geometry sat:{}".format(outFeature))
            Utils.addValueToField(
                outFeature, "bid", building.bid, outLayer, logger)
            Utils.addValueToField(outFeature, "circumf",
                                  building.circumf, outLayer, logger)
            Utils.addValueToField(
                outFeature, "foundation", building.foundation, outLayer, logger)
            Utils.addValueToField(outFeature, "structure",
                                  building.structure, outLayer, logger)
            Utils.addValueToField(outFeature, "status",
                                  building.status, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sv_shr", building.max_sv_short, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sh_shr", building.max_sh_short, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sv_tot", building.max_sv_total, outLayer, logger)
            Utils.addValueToField(outFeature, "max_angle",
                                  building.max_angle, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_strain", building.max_strain, outLayer, logger)
            Utils.addValueToField(outFeature, "max_p_stra",
                                  building.max_principal_strain, outLayer, logger)
            Utils.addValueToField(
                outFeature, "vulnerab", building.vulnerability, outLayer, logger)
            Utils.addValueToField(
                outFeature, "risk_tots", building.risk_totset, outLayer, logger)
            Utils.addValueToField(
                outFeature, "risk_angle", building.risk_angle, outLayer, logger)
            Utils.addValueToField(outFeature, "filter_val",
                                  filterValue, outLayer, logger)

            # logger.debug(
            #     "writeBuildingsToShape, feature created and all attributes added:{}".format(outFeature))
            outLayer.CreateFeature(outFeature)
            # logger.debug(
            #     "writeBuildingsToShape, feature created and added to outLayer:{}".format(outFeature))
            # outLayer = None
            # logger.debug("Building nr: {0}: geom {1}".format(numBuildings, geom.ExportToWkt()))
            outFeature = None
            numBuildings += 1
        # del outFeature
        # Save and close the data source
        outDataSource = None
        logger.debug(
            "Number of buildings written to shape: {}".format(numBuildings))

    except Exception as e:
        logger.error("Error writing buildings: {}".format(type(e)))
        logger.error("Error writing buildings: {}".format(e.args))
        logger.error("Error writing buildings: {}".format(e))
        raise e



def writeWallsToShape(wallsFN, buildings, working_proj, filterValue, logger):

    fields = [["wid", "long"], ["slope_ang", "float"], [
        "h_te_stra", "float"], ["p_te_stra", "float"], ["filter_val", "string"]]

    try:
        Utils.createShapefile(wallsFN, "line", working_proj, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(wallsFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()
        numWalls = 0
        for building in buildings:
            for wall in building.walls:

                coord = [[wall.corner_start.x, wall.corner_start.y],
                         [wall.corner_end.x, wall.corner_end.y]]

                # logger.debug(
                #     "writeWallsToShape, coord created:{}".format(coord))
                geom = Utils.createLine(coord, logger)
                # logger.debug(
                #     "writeWallsToShape, geom created:{}".format(geom))

                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(geom)
                Utils.addValueToField(
                    outFeature, "wid", wall.wid, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "filter_val", filterValue, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "slope_ang", wall.slope_angle, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "h_te_stra", wall.hor_tensile_strain, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "p_te_stra", wall.principal_tensile_strain, outLayer, logger)
                outLayer.CreateFeature(outFeature)
                numWalls += 1
                outFeature = None
            # del outFeature
        # Save and close the data source
        outDataSource = None
        logger.debug("Number of walls written to shape: {}".format(numWalls))

    except Exception as e:
        logger.error(type(e))
        logger.error(e.args)
        logger.error(e)
        raise e

def writeCornersToShape(cornersFN, buildings, projection, filterValue, logger):

    #logger.debug(f"Writing corners to shape, Number of buildings:{len(buildings)}, shapefile: {cornersFN}")

    fields = [["cid", "long"], ["dtb", "float"], ["sv_short", "float"], ["sv_long", "float"], ["sv_tot", "float"], [
        "sh_short", "float"], ["near_dist", "float"], ["near_angle", "float"], ["filter_val", "string"], ["porewp_red", "float"]]

    try:
        Utils.createShapefile(cornersFN, "point", projection, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(cornersFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()

        numCorners = 0
        for building in buildings:
            for corner in building.corners:

                geom = Utils.createPoint(
                    corner.x, corner.y)
                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(geom)
                Utils.addValueToField(
                    outFeature, "cid", corner.cid, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "dtb", corner.dtb, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "filter_val", filterValue, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_short", corner.sv_short, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_long", corner.sv_long, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sh_short", corner.sh_short, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "near_dist", corner.near_dist, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "near_angle", corner.near_angle, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_tot", corner.sv_tot, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "porewp_red", corner.porewp_red, outLayer, logger)
                outLayer.CreateFeature(outFeature)
                outFeature = None
                numCorners += 1
            # del outFeature
        # Save and close the data source
        # del outDataSource
        outDataSource = None
        logger.debug(
            "Number of corners written to shape: {}".format(numCorners))
    except Exception as e:
        logger.error(type(e))
        logger.error(e.args)
        logger.error(e)
        raise e

def appendZValuesFromRaster(corner_list, rasterFN, logger=None):
    # if (not logger is None):
    # logger.debug("Logger sent, corner_list: {}".format(corner_list))
    # logger.debug("Logger sent, rasterFN: {}".format(rasterFN))

    src_ds = gdal.Open(str(rasterFN))
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)

    result_corner_list = []

    # numWithoutBr = 0
    for corner_i in corner_list:

        # Convert from map to pixel coordinates.
        # Only works for geotransforms with no rotation.
        px = int((corner_i.x - gt[0]) / gt[1])  # x pixel
        py = int((corner_i.y - gt[3]) / gt[5])  # y pixel

        dtb = rb.ReadAsArray(px, py, 1, 1)

        if (dtb is None):
            # if (not logger is None):
            #     logger.debug("Returning None")
            return None
            # numWithoutBr += 1
            # continue

        if (dtb[0][0] < -50):
            # if (not logger is None):
                # logger.debug("DTB: {}".format(dtb[0][0]))
            return None
        corner_i.dtb = dtb[0][0]
        result_corner_list.append(corner_i)

        # if (numWithoutBr > 0):
        #      logger.debug("Skipping corner in point with no detph to bedrock for {} corners".format(numWithoutBr))
    return result_corner_list

