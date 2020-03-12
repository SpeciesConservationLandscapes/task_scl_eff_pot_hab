import argparse
import ee
from datetime import datetime, timezone
from task_base import SCLTask

class SCLeffPotHab(EETask):
    ee_rootdir = "projects/SCL/v1"

    inputs = {
        "extent": {
            "ee_type": EETask.IMAGE,
            "ee_path": "Panthera_tigris/sumatra_poc_aoi",
            "maxage": 5  # years
        },
        "structural_hab": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": "Panthera_tigris/geographies/Sumatra/hab/structural_habitat",
            "maxage": 1  # years
        },
        "historic_range": {
            "ee_type": EETask.IMAGE,
            "ee_path": "Panthera_tigris/source/Inputs_2006/hist",
            "maxage": 10  # years
        },
        "elevation": {
            "ee_type": EETask.IMAGE,
            "ee_path": "CGIAR/SRTM90_V4",
            "maxage": 10  # years
        },
        "water_mask": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "maxage": 10  # years
        },
        "hii": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/sumatra_poc/final/weighted_hii",
            "maxage": 1  # years
        },
        "extirpated_range": {
            "ee_type": EETask.IMAGE,
            "ee_path": "Panthera_tigris/source/Inputs_2006/extirp_fin",
            "maxage": 10  # years
        },
        "density": {
            "ee_type": EETask.IMAGE,
            "ee_path": "Panthera_tigris/source/Inputs_2006/density",
            "maxage": 10  # years
        },
        "country_aoi": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "Panthera_tigris/sumatra_poc_aoi",
            "maxage": 10,
        }
        }

        thresholds = {
            'structural_lc': 50,
            'elevation': 3350,
            'structural_patch': 5,
            'hii': 35,
            'probability': 0.4
        }





#     gpw_cadence = 5
#
#
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    def calc(self):


        watermask = ee.Image(self.inputs['watermask']['ee_path'])
        print(watermask)
        print('taskdate is {}'.format(self.taskdate))
#
#         ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
#         gpw_prior = gpw.filterDate(ee_taskdate.advance(-self.gpw_cadence, 'year'), ee_taskdate).first()
#         gpw_later = gpw.filterDate(ee_taskdate, ee_taskdate.advance(self.gpw_cadence, 'year')).first()
#         gpw_diff = gpw_later.subtract(gpw_prior)
#         numerator = ee_taskdate.difference(gpw_prior.date(), 'day')
#         gpw_diff_fraction = gpw_diff.multiply(numerator.divide(self.gpw_cadence * 365))
#         gpw_taskdate = gpw_prior.add(gpw_diff_fraction)
#         gpw_taskdate_300m = gpw_taskdate.resample().reproject(crs=self.crs, scale=self.scale)
#
#
#
#         hii_popdens_driver = gpw_taskdate_300m.add(ee.Image(1))\
#             .log()\
#             .multiply(ee.Image(3.333))\
#             .updateMask(watermask)
#
#         self.export_image_ee(hii_popdens_driver, '{}/{}'.format(self.ee_driverdir, 'hii_popdens_driver'))
#
#     def check_inputs(self):
#         super().check_inputs()
#         # add any task-specific checks here, and set self.status = self.FAILED if any fail
#
#
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-d', '--taskdate', default=datetime.now(timezone.utc).date())
#     options = parser.parse_args()
#     popdens_task = HIIPopulationDensity(**vars(options))
#     popdens_task.run()
