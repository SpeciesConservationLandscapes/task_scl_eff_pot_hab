import argparse
import ee
import time
from datetime import datetime, timezone
from task_base import SCLTask


class SCLPolygons(SCLTask):
    scale = 1000
    inputs = {
        "structural_habitat": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/structural_habitat",
            "maxage": 10,  # until full image collection of structural habitat
        },
        "hii": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/hii",
            "maxage": 10,
        },
        "probability": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "scenario_probability",
            "static": True,
        },
        "historic_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/historical_range_img_200914",
            "static": True,
        },
        "extirpated": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/extirp_fin",
            "static": True,
        },
        "elevation": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "CGIAR/SRTM90_V4",
            "static": True,
        },
        "water": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
        "ecoregions": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "RESOLVE/ECOREGIONS/2017",
            "static": True,
        },
        "density": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/density/tiger_obs_density",
            "static": True,
        },
    }
    thresholds = {
        "structural_habitat": 4,
        "elev": 3350,
        "patch_size": 5,  # sq km
        "upper_limit_core_size": 625,  # sq km or number of pixels at 1000m scale
        "lower_limit_fragment_size": 3,  # sq km or number of pixels at 1000m scale
        "hii": 12,
        "probability": 1,
        "connectity_distance": 2,
    }
    density_values = {
        "n_core_animals": 5,
        "core_to_step_ratio": 0.1,
        "core_size_limits": {"min": 30, "max": 625},
        "step_size_limits": {"min": 3, "max": 63},
    }

    def scenario_probability(self):
        return f"{self.ee_rootdir}/probability"

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
        self.ecoregions = ee.FeatureCollection(self.inputs["ecoregions"]["ee_path"])
        self.density = ee.FeatureCollection(self.inputs["density"]["ee_path"])
        self.water = ee.Image(self.inputs["water"]["ee_path"])
        self.set_aoi_from_ee(self.inputs["probability"]["ee_path"])

    def calc(self):
        def density_to_patch_size(ft):
            med_density_eco = ee.Number(ft.get("MED_DENSITY_ECO"))
            med_density_biome = ee.Number(ft.get("MED_DENSITY_BIOME"))
            med_density_eco_biome = ee.Number(ft.get("MED_DENSITY_ECO_BIOME"))
            density_val = ee.Algorithms.If(
                med_density_eco.gt(0),
                med_density_eco,
                ee.Algorithms.If(
                    med_density_biome.gt(0),
                    med_density_biome,
                    ee.Algorithms.If(
                        med_density_eco_biome.gt(0), med_density_eco_biome, 1
                    ),
                ),
            )
            min_core_size = (
                ee.Number(self.density_values["n_core_animals"])
                .divide(density_val)
                .multiply(100)
                .int()
            )
            return ee.Feature(
                None, {"ECO_ID": ft.get("ECO_ID"), "min_core_size": min_core_size}
            )

        def dilate(image, distance):
            dialated_image = image.fastDistanceTransform(distance).sqrt()
            return dialated_image.lte(distance).selfMask()

        elev_mask = self.elevation.lt(self.thresholds["elev"]).selfMask()
        str_hab_mask = self.structural_habitat.gt(
            self.thresholds["structural_habitat"]
        ).selfMask()
        current_range = self.historic_range.updateMask(self.extirpated).selfMask()
        low_hii_mask = self.hii.lte(self.thresholds["hii"]).selfMask()
        str_hab = (
            str_hab_mask.updateMask(elev_mask)
            .updateMask(low_hii_mask)
            .updateMask(self.water)
            .updateMask(current_range)
        )

        connected_potential_habitat = str_hab.connectedPixelCount(
            self.density_val["core_size_limits"]["max"], True
        ).reproject(crs="EPSG:4326", scale=self.scale)

        potential_habitat = connected_potential_habitat.updateMask(
            connected_potential_habitat.gte(
                self.density_values["step_size_limits"]["min"]
            )
        ).selfMask()

        density = self.density.map(density_to_patch_size)

        ecoregion_image = self.ecoregions.reduceToImage(
            properties=["ECO_ID"], reducer=ee.Reducer.mode()
        )

        ecoregion_id = density.aggregate_array("ECO_ID")
        core_size = density.aggregate_array("min_core_size")

        min_patch_size = (
            ecoregion_image.remap(ecoregion_id, core_size)
            .int()
            .clamp(
                self.density_values["core_size_limits"]["min"],
                self.density_values["core_size_limits"]["max"],
            )
            .updateMask(potential_habitat)
        )

        min_stepping_stone_size = (
            min_patch_size.multiply(self.density_values["core_to_step_ratio"])
            .int()
            .clamp(
                self.density_values["step_size_limits"]["min"],
                self.density_values["step_size_limits"]["max"],
            )
            .updateMask(potential_habitat)
        )

        potential_core = (
            ee.Image(0).where(potential_habitat.gte(min_patch_size), 1).selfMask()
        )

        potential_stepping_stone = (
            ee.Image(0)
            .where(
                potential_habitat.lt(min_patch_size).And(
                    potential_habitat.gte(min_stepping_stone_size)
                ),
                1,
            )
            .selfMask()
        )

        potential_core = (
            dilate(potential_core, self.thresholds["connectity_distance"])
            .reproject(crs="EPSG:4326", scale=self.scale)
            .multiply(3)
            .unmask(0)
            .updateMask(self.water)
        )
        potential_stepping_stone = (
            dilate(potential_stepping_stone, self.thresholds["connectity_distance"])
            .reproject(crs="EPSG:4326", scale=self.scale)
            .unmask(0)
            .updateMask(self.water)
        )

        scl_polys = potential_core.add(potential_stepping_stone).selfMask()

        # Temp geometry for RFE scenario
        # geometry = ee.Geometry.Polygon(
        #     [
        #         [
        #             [129.28241561583914, 51.490498638741165],
        #             [129.28241561583914, 42.20681310459823],
        #             [140.35663436583914, 42.20681310459823],
        #             [140.35663436583914, 51.490498638741165],
        #         ]
        #     ],
        #     None,
        #     False,
        # )

        scl_polys = (
            ee.Image(1)
            .updateMask(scl_polys)
            .addBands(scl_polys)
            .addBands(self.probability)
            .rename(["constant", "scl_poly", "probability"])
            .reduceToVectors(
                reducer=ee.Reducer.max(),
                geometry=self.extent,
                scale=self.scale,
                crs="EPSG:4326",
                maxPixels=self.ee_max_pixels,
            )
        )

        export_path = "scl_poly/{}/scl_polys".format(self.taskdate)

        self.export_fc_ee(scl_polys, export_path)

    def check_inputs(self):
        super().check_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    parser.add_argument("-s", "--species", default="Panthera_tigris")
    parser.add_argument("--scenario", default=SCLTask.CANONICAL)
    options = parser.parse_args()
    sclstats_task = SCLPolygons(**vars(options))
    sclstats_task.run()
