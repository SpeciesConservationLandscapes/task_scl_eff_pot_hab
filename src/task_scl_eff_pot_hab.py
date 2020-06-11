import argparse
import ee
from datetime import datetime, timezone
from task_base import SCLTask

import os
import json
service_account_key = os.environ.get("SERVICE_ACCOUNT_KEY")
service_account_name = json.loads(service_account_key)["client_email"]
credentials = ee.ServiceAccountCredentials(
    service_account_name, key_data=service_account_key
)
ee.Initialize(credentials)

# sumatra_aoi = ee.FeatureCollection("projects/SCL/v1/Panthera_tigris/sumatra_poc_aoi")

class SCLPolygons(SCLTask):
    ee_rootdir = "projects/SCL/v1"
    inputs = {
        "structural_habitat": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/geographies/Sumatra/hab/structural_habitat",
            "maxage": 1,  # years
        },
        "hii": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/sumatra_poc/hii",
            "maxage": 1,
        },
        "historic_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/hist",
            "maxage": 1 / 365,
        },
        "extirpated": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/extirp_fin",
            "maxage": 1 / 365,
        },
        "elevation": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "CGIAR/SRTM90_V4",
            "maxage": 10,
        },
        "water": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "maxage": 1,
        },
        "ecoregion_country": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/eco_cntry",
            "maxage": 10,
        },
        "density": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/density",
            "maxage": 1,
        },
        "aoi": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/sumatra_poc_aoi",
            "maxage": 1,
        }
    }
    thresholds = {
        "str_lc": 50,
        "elev": 3350,
        "patch_size": 5,
        "hii": 10,
        "probability": 0.4
    }
    remap_lists = {
        "density": [10, 20, 30, 40, 50, 60, 70, 80, 90],
        "min_patch_size": [250, 250, 625, 250, 140, 70, 45, 30, 250],
        "min_stepping_stone_size": [25, 25, 63, 25, 14, 7, 5, 3, 25]
    }

    def _scl_path(self, scltype):
        if scltype is None or scltype not in self.inputs:
            raise TypeError("Missing or incorrect scltype for setting scl path")
        return "{}/{}/scl_poly/{}/{}".format(
            self.ee_rootdir, self.species, self.taskdate, scltype
        )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.error_margin = ee.ErrorMargin(1)
    #     self.area_proj = "EPSG:5070"
    #     self.structural_habitat = ee.Image(ee.ImageCollection(self.inputs["structural_habitat"]["ee_path"]).first())
    #     self.hii = ee.Image(ee.ImageCollection(self.inputs["hii"]["ee_path"]).first())
    #     self.historic_range = ee.Image(self.inputs["historic_range"]["ee_path"])
    #     self.extirpated = ee.Image(self.inputs["extirpated"]["ee_path"])
    #     self.elevation = ee.Image(self.inputs["elevation"]["ee_path"])
    #     self.ecoregion_country = ee.Image(self.inputs["ecoregion_country"]["ee_path"])
    #     self.density = ee.Image(self.inputs["density"]["ee_path"])
    #     self.water = ee.Image(self.inputs["water"]["ee_path"])

    def delineate_polygons(self):
        structural_habitat = ee.Image(ee.ImageCollection(self.inputs["structural_habitat"]["ee_path"]).first())
        hii = ee.Image(ee.ImageCollection(self.inputs["hii"]["ee_path"]).first())
        historic_range = ee.Image(self.inputs["historic_range"]["ee_path"])
        extirpated = ee.Image(self.inputs["extirpated"]["ee_path"])
        elevation = ee.Image(self.inputs["elevation"]["ee_path"])
        ecoregion_country = ee.Image(self.inputs["ecoregion_country"]["ee_path"])
        density = ee.Image(self.inputs["density"]["ee_path"])
        water = ee.Image(self.inputs["water"]["ee_path"])
        aoi = ee.FeatureCollection(self.inputs["aoi"]["ee_path"])

        elev_mask = elevation.lt(self.thresholds["elev"]).selfMask()
        str_hab_mask = structural_habitat.gt(self.thresholds["str_lc"]).selfMask()
        str_hab = str_hab_mask.updateMask(elev_mask).updateMask(water)
        str_hab_patch = str_hab.connectedPixelCount(16, True).gte(self.thresholds["patch_size"]).selfMask().reproject(crs = "EPSG:4326", scale = 300)
        current_range = extirpated.updateMask(historic_range.eq(0)).eq(1).selfMask()
        connected_structural_habitat = str_hab_patch.updateMask(current_range)
        potential_habitat = connected_structural_habitat.updateMask(hii.lt(self.thresholds["hii"])).selfMask()
        excluded_habitat = connected_structural_habitat.updateMask(hii.gte(self.thresholds["hii"])).selfMask()
        effective_potential_habitat = potential_habitat.addBands(excluded_habitat).rename(['eff_pot_hab', 'excl_hab'])

        probability_calc = elevation.expression(
            '1 - (2.71828 ** (-(2.71828 ** (- 0.069 - 0.044 * ((dem - 639) / 592) - 0.12 * ((hii - 7.72) / 5.47)))))', {
            'dem' : elevation,
            'hii' : hii
            })
        high_probability = probability_calc.gt(self.thresholds["probability"]).multiply(2)
        low_probability = probability_calc.lte(self.thresholds["probability"]).multiply(1)
        probability = high_probability.add(low_probability).selfMask()

        hii_sum = effective_potential_habitat.select(0).unmask(0) \
            .add(effective_potential_habitat.select(1).unmask(0).multiply(2)) \
            .selfMask() \
            .remap([1,2],[2,3])

        hii_ecoregion = ecoregion_country.add(hii_sum)

        minPatchSize = density.remap(self.remap_lists["density"], self.remap_lists["min_patch_size"]) \
            .updateMask(hii_ecoregion)
        minStpStnSize = density.remap(self.remap_lists["density"], self.remap_lists["min_stepping_stone_size"]) \
            .updateMask(hii_ecoregion)

        hii_grp = hii_ecoregion.addBands([probability, hii_sum, ee.Image.pixelArea().multiply(0.000001), minPatchSize, minStpStnSize]) \
            .reproject(crs = 'EPSG:4326', scale = 300)

        polygonReducer = ee.Reducer.max().combine(ee.Reducer.mode().unweighted(), 'hii_') \
            .combine(ee.Reducer.sum().unweighted(), 'area_') \
            .combine(ee.Reducer.mode().unweighted(), 'patch_') \
            .combine(ee.Reducer.mode().unweighted(), 'stpStn_') \

        hii_grp_poly = hii_grp.reduceToVectors(
            reducer = polygonReducer,
            geometry = aoi.geometry().bounds(),
            crs = 'EPSG:4326',
            scale = 300,
            maxPixels = 1e13,
            geometryInNativeProjection = True
            ) \
            .select(['area_sum', 'hii_mode', 'label', 'max', 'patch_mode', 'stpStn_mode'], \
                ['patch_area', 'hii', 'label', 'prob', 'min_patch_area', 'min_stp_stn_area'])

        return hii_grp_poly.first()
        # blob = "ls_stats/{}/{}/{}".format(self.species, self.taskdate, landscape_key)
        # self.export_fc_ee(ls_countries_biomes_pas, "scl-pipeline", blob)

foo = SCLPolygons.delineate_polygons(SCLPolygons)
print(foo.getInfo())




    def calc(self):
        self.delineate_polygons()

    def check_inputs(self):
        super().check_inputs()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    parser.add_argument("-s", "--species", default="Panthera_tigris")
    options = parser.parse_args()
    sclstats_task = SCLPolygons(**vars(options))
    sclstats_task.run()
