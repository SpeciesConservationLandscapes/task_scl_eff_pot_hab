import argparse
import ee
from datetime import datetime, timezone
from task_base import SCLTask


import ee
from datetime import date
ee.Initialize()

str_hab_collection = ee.ImageCollection("projects/SCL/v1/Panthera_tigris/geographies/Sumatra/hab/structural_habitat")
hist_range = ee.Image("projects/SCL/v1/Panthera_tigris/source/Inputs_2006/hist")
elev = ee.Image("CGIAR/SRTM90_V4")
water_mask = ee.Image("projects/HII/v1/source/phys/watermask_jrc70_cciocean")
extirpated_image = ee.Image("projects/SCL/v1/Panthera_tigris/source/Inputs_2006/extirp_fin")
hii_collection = ee.ImageCollection("projects/HII/v1/sumatra_poc/final/weighted_hii")
tcl_2006 = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/scl_poly/2006-01-01/scl_species")
eco_country_image = ee.Image("projects/SCL/v1/Panthera_tigris/source/Inputs_2006/eco_cntry")
density_image = ee.Image("projects/SCL/v1/Panthera_tigris/source/Inputs_2006/density")
sumatra_aoi = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/sumatra_poc_aoi")
tcl_survey = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/scl_poly/2006-01-01/scl_survey")
tcl_restoration = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/scl_poly/2006-01-01/scl_restoration")
tcl_frag = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/scl_poly/2006-01-01/scl_fragment")

# Get current date and set as ee date
currentDate = ee.Date.parse("YYYYMMDD", date.today().strftime("%Y%m%d"))

# Define threshold variables
threshold_lc = 50
threshold_elev = 3350
threshold_lc_patch = 5
thershold_hii = 30
thershold_prob = 0.4
thershold_dispersal = 2000

# Remap lists for tiger density estimates to minium core size and minimum stepping stone size
densityList = [10, 20, 30, 40, 50, 60, 70, 80, 90]
minPatchSizeList = [250, 250, 625, 250, 140, 70, 45, 30, 250]
minStpStnSizeList = [25, 25, 63, 25, 14, 7, 5, 3, 25]

# Define key functions

def f_gt_threshold_mask(image, threshold):
    return image.gt(threshold).selfMask()

def f_lt_threshold_mask (image, threshold):
    return image.lt(threshold).selfMask()

def f_coll_filter(collection, date, maxAge):
    image = ee.Image(collection.filterDate(date.advance(-maxAge, 'year'), date) \
    .sort('system:time_start', False).first())
    return image

def f_eff_pot_hab(elev_image, str_hab_image, water_mask_image, extirpated_image, historic_image, hii_image, elev_threshold, str_hab_threshold, min_patch_size, hii_threshold):
    elev_mask = f_lt_threshold_mask(elev_image, elev_threshold)
    str_hab_mask = f_gt_threshold_mask(str_hab_image, str_hab_threshold)
    str_hab = str_hab_mask.updateMask(elev_mask) \
        .updateMask(water_mask)
    con_str_hab = f_gt_threshold_mask(str_hab
        .connectedPixelCount(6, True), min_patch_size) \
        .selfMask() \
        .reproject(crs = "EPSG:4326", scale = 1000)
    current = extirpated_image.updateMask(historic_image.eq(0)) \
        .eq(1) \
        .selfMask()
    con_str_hab = con_str_hab.updateMask(current)
    pot_hab = con_str_hab.updateMask(hii_image.lt(hii_threshold)) \
        .selfMask()
    excl_hab = con_str_hab.updateMask(hii_image.gte(hii_threshold)) \
        .selfMask()
    return pot_hab.addBands(excl_hab) \
        .rename(['eff_pot_hab', 'excl_hab'])


def f_prob(elev_image, hii_image, prob_threshold):
    probablity = elev_image.expression(
        '1 - (2.71828 ** (-(2.71828 ** (- 0.069 - 0.044 * ((dem - 639) / 592) - 0.12 * ((hii - 7.72) / 5.47)))))', {
        'dem' : elev_image,
        'hii' : hii_image
        })
    high = probablity.gt(prob_threshold).multiply(2)
    low = probablity.lte(prob_threshold).multiply(1)
    return high.add(low).selfMask()

def f_hii_sum(eff_pot_hab_img):
    return eff_pot_hab_img.select(0) \
        .unmask(0) \
        .add(eff_pot_hab_img \
        .select(1) \
        .unmask(0) \
        .multiply(2)) \
        .selfMask() \
        .remap([1,2],[2,3])

polygonReducer = ee.Reducer.max().combine(ee.Reducer.mode().unweighted(), 'hii_') \
    .combine(ee.Reducer.sum().unweighted(), 'area_') \
    .combine(ee.Reducer.mode().unweighted(), 'patch_') \
    .combine(ee.Reducer.mode().unweighted(), 'stpStn_') \

def f_scl_species(eff_pot_hab_img, prob_img, eco_country_img, density_img, density_list, minimum_patch_size, minimum_stepping_stone_size, reducer):
    hii_sum = f_hii_sum(eff_pot_hab_img)
    hii_hab = eco_country_img.add(hii_sum)
    minPatchSize = density_img.remap(density_list, minimum_patch_size) \
        .updateMask(hii_hab)
    minStpStnSize = density_img.remap(density_list, minimum_stepping_stone_size) \
        .updateMask(hii_hab)
    hii_grp = hii_hab.addBands([prob_img, hii_sum, ee.Image.pixelArea().multiply(0.000001), minPatchSize, minStpStnSize]) \
        .reproject(crs = 'EPSG:4326', scale = 1000)

    hii_grp_poly = hii_grp.reduceToVectors(
        reducer = polygonReducer,
        geometry = sumatra_aoi.geometry().bounds(),
        crs = 'EPSG:4326',
        scale = 1000,
        maxPixels = 1e10,
        geometryInNativeProjection = True
        ) \
        .select(['area_sum', 'hii_mode', 'label', 'max', 'patch_mode', 'stpStn_mode'], \
            ['patch_area', 'hii', 'label', 'prob', 'min_patch_area', 'min_stp_stn_area'])
    return(hii_grp_poly)

def connectivityBuffer(ft):
    return ft.buffer(2000)

str_hab_img = f_coll_filter(str_hab_collection, currentDate, 10)
hii_img = ee.Image(hii_collection.first())

eff_pot_hab = f_eff_pot_hab(elev, str_hab_img, water_mask, extirpated_image, hist_range, hii_img, threshold_elev, threshold_lc, threshold_lc_patch, thershold_hii) \
    .set({'system:time_start': currentDate.millis()})

prob_image = f_prob(elev, hii_img, thershold_prob)

scl_species = f_scl_species(eff_pot_hab, prob_image, eco_country_image, density_image, densityList, minPatchSizeList, minStpStnSizeList, polygonReducer)

core = scl_species.filter(ee.Filter.And(ee.Filter.lessThanOrEquals(rightField = 'patch_area', leftField = 'min_patch_area'), \
    ee.Filter.eq('prob', 2), \
    ee.Filter.eq('hii', 2)))

survey = scl_species.filter(ee.Filter.And(ee.Filter.lessThanOrEquals(rightField = 'patch_area', leftField = 'min_patch_area'),
    ee.Filter.eq('prob', 1)))

restoration = scl_species.filter(ee.Filter.And(ee.Filter.lessThanOrEquals(rightField = 'patch_area', leftField = 'min_patch_area'),
    ee.Filter.eq('prob', 2),
    ee.Filter.eq('hii', 3)))

stepping_stones = scl_species.filter(ee.Filter.And(ee.Filter.greaterThan(rightField = 'patch_area', leftField = 'min_patch_area'),
    ee.Filter.lessThanOrEquals(rightField = 'patch_area', leftField = 'min_stp_stn_area'),
    ee.Filter.eq('prob', 2),
    ee.Filter.eq('hii', 2)))

coreConnectivity = core.distance(4000).unmask(0)

survey = coreConnectivity.reduceRegions(
  collection = survey,
  reducer = ee.Reducer.max(),
  scale = 1000,
  crs = 'EPSG:4326'
)

restoration = coreConnectivity.reduceRegions(
  collection = restoration,
  reducer = ee.Reducer.max(),
  scale = 1000,
  crs = 'EPSG:4326'
)

stepping_stones = coreConnectivity.reduceRegions(
  collection = stepping_stones,
  reducer = ee.Reducer.max(),
  scale = 1000,
  crs = 'EPSG:4326'
)

connectedSurvey = survey.filter(ee.Filter.greaterThan('max', 0))
nonConnectedSurvey = survey.filter(ee.Filter.equals('max', 0))
connectedRestoration = restoration.filter(ee.Filter.greaterThan('max', 0))
nonConnectedRestoration = restoration.filter(ee.Filter.equals('max', 0))
connectedSteppingStones = stepping_stones.filter(ee.Filter.greaterThan('max', 0))
nonConnectedSteppingStones = stepping_stones.filter(ee.Filter.equals('max', 0))

coreBuffered = core.map(connectivityBuffer)

scl_species = coreBuffered.merge(connectedSteppingStones).union(1)

# task = ee.batch.Export.image.toAsset(
#     image = eff_pot_hab,
# 	description = 'eff_pot_hab',
# 	assetId = 'projects/SCL/v1/Panthera_tigris/source/Outputs_2020_test/eff_pot_hab_1',
# 	crs = 'EPSG:4326',
# 	region = sumatra_aoi.geometry().bounds(),
# 	scale = 1000)
# task.start()

task = ee.batch.Export.table.toAsset(
    collection = scl_species,
	description = 'scl_species',
	assetId = 'projects/SCL/v1/Panthera_tigris/source/Outputs_2020_test/scl_species_1')
task.start()


# class SCLeffPotHab(EETask):
#     ee_rootdir = "projects/SCL/v1"
#
#     inputs = {
#         "extent": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "Panthera_tigris/sumatra_poc_aoi",
#             "maxage": 5  # years
#         },
#         "structural_hab": {
#             "ee_type": EETask.IMAGECOLLECTION,
#             "ee_path": "Panthera_tigris/geographies/Sumatra/hab/structural_habitat",
#             "maxage": 1  # years
#         },
#         "historic_range": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "Panthera_tigris/source/Inputs_2006/hist",
#             "maxage": 10  # years
#         },
#         "elevation": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "CGIAR/SRTM90_V4",
#             "maxage": 10  # years
#         },
#         "water_mask": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
#             "maxage": 10  # years
#         },
#         "hii": {
#             "ee_type": EETask.IMAGECOLLECTION,
#             "ee_path": "projects/HII/v1/sumatra_poc/final/weighted_hii",
#             "maxage": 1  # years
#         },
#         "extirpated_range": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "Panthera_tigris/source/Inputs_2006/extirp_fin",
#             "maxage": 10  # years
#         },
#         "density": {
#             "ee_type": EETask.IMAGE,
#             "ee_path": "Panthera_tigris/source/Inputs_2006/density",
#             "maxage": 10  # years
#         },
#         "country_aoi": {
#             "ee_type": SCLTask.FEATURECOLLECTION,
#             "ee_path": "Panthera_tigris/sumatra_poc_aoi",
#             "maxage": 10,
#         }
#         }
#
#         thresholds = {
#             'structural_lc': 50,
#             'elevation': 3350,
#             'structural_patch': 5,
#             'hii': 35,
#             'probability': 0.4
#         }
#
#
#
#
#
# #     gpw_cadence = 5
# #
# #
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))
#
#     def calc(self):
#
#
#         watermask = ee.Image(self.inputs['watermask']['ee_path'])
#         print(watermask)
#         print('taskdate is {}'.format(self.taskdate))
# #
# #         ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
# #         gpw_prior = gpw.filterDate(ee_taskdate.advance(-self.gpw_cadence, 'year'), ee_taskdate).first()
# #         gpw_later = gpw.filterDate(ee_taskdate, ee_taskdate.advance(self.gpw_cadence, 'year')).first()
# #         gpw_diff = gpw_later.subtract(gpw_prior)
# #         numerator = ee_taskdate.difference(gpw_prior.date(), 'day')
# #         gpw_diff_fraction = gpw_diff.multiply(numerator.divide(self.gpw_cadence * 365))
# #         gpw_taskdate = gpw_prior.add(gpw_diff_fraction)
# #         gpw_taskdate_300m = gpw_taskdate.resample().reproject(crs=self.crs, scale=self.scale)
# #
# #
# #
# #         hii_popdens_driver = gpw_taskdate_300m.add(ee.Image(1))\
# #             .log()\
# #             .multiply(ee.Image(3.333))\
# #             .updateMask(watermask)
# #
# #         self.export_image_ee(hii_popdens_driver, '{}/{}'.format(self.ee_driverdir, 'hii_popdens_driver'))
# #
# #     def check_inputs(self):
# #         super().check_inputs()
# #         # add any task-specific checks here, and set self.status = self.FAILED if any fail
# #
# #
# # if __name__ == "__main__":
# #     parser = argparse.ArgumentParser()
# #     parser.add_argument('-d', '--taskdate', default=datetime.now(timezone.utc).date())
# #     options = parser.parse_args()
# #     popdens_task = HIIPopulationDensity(**vars(options))
# #     popdens_task.run()
