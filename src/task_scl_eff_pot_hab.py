import argparse
import ee
from datetime import datetime, timezone
from task_base import SCLTask

class SCLPolygons(SCLTask):
    ee_rootdir = "projects/SCL/v1"
    scale = 300
    inputs = {
        "structural_habitat": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/geographies/Sumatra/hab/structural_habitat",
            "maxage": 1,
        },
        "hii": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/sumatra_poc/hii",
            "maxage": 1,
        },
        "probability": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/geographies/Sumatra/hab/probability",
            "maxage": 1,
        },
        "historic_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/hist",
            "static": True
        },
        "extirpated": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/extirp_fin",
            "static": True
        },
        "elevation": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "CGIAR/SRTM90_V4",
            "static": True
        },
        "water": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True
        },
        "ecoregion_country": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/eco_cntry",
            "static": True
        },
        "ecoregions": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "RESOLVE/ECOREGIONS/2017",
            "static": True
        },
        "density": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/density/tiger_obs_density",
            "static": True
        },
    }
    thresholds = {
        "str_lc": 50,
        "elev": 3350,
        "patch_size": 5,
        "hii": 6,
        "probability": 1,
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.structural_habitat, _ = self.get_most_recent_image(
             ee.ImageCollection(self.inputs["structural_habitat"]["ee_path"])
        )
        self.probability, _ = self.get_most_recent_image(
             ee.ImageCollection(self.inputs["probability"]["ee_path"])
        )
        self.hii, _ = self.get_most_recent_image(
             ee.ImageCollection(self.inputs["hii"]["ee_path"])
        )
        self.historic_range = ee.Image(self.inputs["historic_range"]["ee_path"])
        self.extirpated = ee.Image(self.inputs["extirpated"]["ee_path"])
        self.elevation = ee.Image(self.inputs["elevation"]["ee_path"])
        self.ecoregion_country = ee.Image(self.inputs["ecoregion_country"]["ee_path"])
        self.ecoregions = ee.FeatureCollection(self.inputs["ecoregions"]["ee_path"])
        self.density = ee.FeatureCollection(self.inputs["density"]["ee_path"])
        self.water = ee.Image(self.inputs["water"]["ee_path"])
        self.set_aoi_from_ee(
            "{}/{}/sumatra_poc_aoi".format(self.ee_rootdir, self.species)
        )

    def calc(self):
        elev_mask = self.elevation.lt(self.thresholds["elev"]).selfMask()
        str_hab_mask = self.structural_habitat.gt(self.thresholds["str_lc"]).selfMask()
        str_hab = str_hab_mask.updateMask(elev_mask).updateMask(self.water)
        str_hab_patch = (
            str_hab.connectedPixelCount(16, True)
            .gte(self.thresholds["patch_size"])
            .selfMask()
            .reproject(crs="EPSG:4326", scale=self.scale)
        )
        current_range = (
            self.extirpated.updateMask(self.historic_range.eq(0)).eq(1).selfMask()
        )

        connected_structural_habitat = str_hab_patch.updateMask(current_range)
        potential_habitat = connected_structural_habitat.updateMask(
            self.hii.lt(self.thresholds["hii"])
        ).selfMask()
        excluded_habitat = connected_structural_habitat.updateMask(
            self.hii.gte(self.thresholds["hii"])
        ).selfMask()
        effective_potential_habitat = potential_habitat.addBands(
            excluded_habitat
        ).rename(["eff_pot_hab", "excl_hab"])

        high_probability = self.probability.gte(self.thresholds["probability"]).multiply(
            2
        )

        low_probability = self.probability.lt(self.thresholds["probability"])

        probability_calc = high_probability.add(low_probability).selfMask()

        hii_sum = (
            effective_potential_habitat.select(0)
            .unmask(0)
            .add(effective_potential_habitat.select(1).unmask(0).multiply(2))
            .selfMask()
            .remap([1, 2], [2, 3])
        )

        hii_ecoregion = self.ecoregion_country.add(hii_sum)

        def density_to_patch_size(ft):
            med_density_eco = ee.Number(ft.get('MED_DENSITY_ECO'))
            med_density_biome = ee.Number(ft.get('MED_DENSITY_ECO'))
            med_density_eco_biome = ee.Number(ft.get('MED_DENSITY_ECO'))
            density_val = ee.Algorithms.If(med_density_eco.gt(0),
                med_density_eco,
                ee.Algorithms.If(med_density_biome.gt(0),
                    med_density_biome,
                    ee.Algorithms.If(med_density_eco_biome.gt(0),
                        med_density_eco_biome,
                        1)))
            min_core_size = ee.Number(500).divide(density_val).int()
            return ee.Feature(None, {'ECO_ID': ft.get('ECO_ID'), 'min_core_size': min_core_size})

        density = self.density.map(density_to_patch_size)

        ecoregion_image = self.ecoregions.reduceToImage(
            properties=['ECO_ID'],
            reducer=ee.Reducer.mode()
        )

        ecoregion_id = density.aggregate_array('ECO_ID')
        core_size = density.aggregate_array('min_core_size')

        min_patch_size = ecoregion_image.remap(
            ecoregion_id, core_size).int().clamp(30, 625).updateMask(hii_ecoregion)
        min_stepping_stone_size = min_patch_size.divide(10).int().updateMask(hii_ecoregion)

        hii_grp = hii_ecoregion.addBands(
            [
                probability_calc,
                hii_sum,
                ee.Image.pixelArea().multiply(0.000001),
                min_patch_size,
                min_stepping_stone_size,
            ]
        ).reproject(crs="EPSG:4326", scale=self.scale)

        polygonReducer = (
            ee.Reducer.mode()
            .combine(ee.Reducer.mode().unweighted(), "hii_")
            .combine(ee.Reducer.sum().unweighted(), "area_")
            .combine(ee.Reducer.mode().unweighted(), "patch_")
            .combine(ee.Reducer.mode().unweighted(), "stp_stn_")
        )
        hii_grp_poly = hii_grp.reduceToVectors(
            reducer=polygonReducer,
            geometry=ee.Geometry.MultiPolygon(self.aoi).bounds(),
            crs=self.crs,
            scale=self.scale,
            maxPixels=self.ee_max_pixels,
            geometryInNativeProjection=True,
        ).select(
            [
                "area_sum",
                "hii_mode",
                "label",
                "mode",
                "patch_mode",
                "stp_stn_mode"
            ],
            [
                "patch_area",
                "hii",
                "label",
                "prob",
                "min_patch_area",
                "min_stp_stn_area",
            ],
        )

        scl_eff_pot_hab = "{}/geographies/Sumatra/scl_poly/{}/scl_eff_pot_hab_test".format(
            self.species, self.taskdate
        )

        self.export_fc_ee(hii_grp_poly, scl_eff_pot_hab)

    def check_inputs(self):
        super().check_inputs()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    parser.add_argument("-s", "--species", default="Panthera_tigris")
    options = parser.parse_args()
    sclstats_task = SCLPolygons(**vars(options))
    sclstats_task.run()
