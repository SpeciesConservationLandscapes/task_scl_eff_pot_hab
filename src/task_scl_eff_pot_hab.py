import argparse
import ee
from datetime import datetime, timezone
from task_base import SCLTask


class SCLPolygons(SCLTask):
    scale = 1000
    ZONES_LABEL = "Biome_zone"
    inputs = {
        "structural_habitat": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/canonical/structural_habitat",
            # TODO: change maxage once final image collection is incorporated
            "maxage": 10,
        },
        "hii": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/hii",
            # TODO: change maxage once final HII image collection is used
            "maxage": 10,
        },
        "probability": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "scenario_probability",
            # TODO: change static to maxage once collection is produced
            "static": True,
        },
        # TODO: add effort input when ready
        "historic_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/historical_range_img_200914",
            "static": True,
        },
        "extirpated_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/source/Inputs_2006/extirp_fin",
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
            "ee_path": "projects/SCL/v1/Panthera_tigris/biome_density",
            "static": True,
        },
        "zones": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "projects/SCL/v1/Panthera_tigris/zones",
            "static": True,
        },
        "zones_image": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "projects/SCL/v1/Panthera_tigris/zones_img",
            "static": True,
        },
    }
    thresholds = {
        "structural_habitat": 0.5,
        "reduce_res_input_pixels": 0.5,
        "structural_habitat_patch_size": 5,  # sq km
        "hii": 18,  # TODO: use this as a default value, if dynamic thresholding fails
        "probability": 1,  # TODO: probability threshold may change based on final probability format
        "connectivity_distance": 2,  # km (1/2 of actual dispersal distance)
        "landscape_size": 3,
        "current_range": 2,
        "landscape_probability": 1,
        "landscape_survey_effort": 1,
    }
    density_values = {
        "n_core_animals": 5,
        "core_to_step_ratio": 0.1,
        "core_size_limits": {"min": 30, "max": 625},  # size min/max in sqkm
        "step_size_limits": {"min": 3, "max": 63},  # size min/max in sqkm
    }

    def histogram_reducer(self, min, max, number):
        return (
            ee.Reducer.count()
            .combine(ee.Reducer.fixedHistogram(min, max, number), None, True)
            .repeat(2)
            .group(2, "zone")
        )

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
        self.extirpated_range = ee.Image(self.inputs["extirpated_range"]["ee_path"])
        self.ecoregions = ee.FeatureCollection(self.inputs["ecoregions"]["ee_path"])
        self.density, density_date = self.get_most_recent_featurecollection(
            self.inputs["density"]["ee_path"]
        )
        self.water = ee.Image(self.inputs["water"]["ee_path"])
        self.zones = ee.FeatureCollection(self.inputs["zones"]["ee_path"])
        self.zones_image = ee.Image(
            self.inputs["zones_image"]["ee_path"]
        )  # TODO: systemize creation/storage of zones_image

        self.scl_poly_filters = {
            "scl_species": ee.Filter.And(
                ee.Filter.gte("size", self.thresholds["landscape_size"]),
                ee.Filter.eq("range", self.thresholds["current_range"]),
                ee.Filter.gte("probability", self.thresholds["landscape_probability"]),
                ee.Filter.gte("effort", self.thresholds["landscape_survey_effort"]),
            ),
            "scl_restoration": ee.Filter.Or(
                ee.Filter.And(
                    ee.Filter.gte("size", self.thresholds["landscape_size"]),
                    ee.Filter.eq("range", self.thresholds["current_range"]),
                    ee.Filter.lt(
                        "probability", self.thresholds["landscape_probability"]
                    ),
                    ee.Filter.gte("effort", self.thresholds["landscape_survey_effort"]),
                ),
                ee.Filter.And(
                    ee.Filter.gte("size", self.thresholds["landscape_size"]),
                    ee.Filter.lt("range", self.thresholds["current_range"]),
                ),
            ),
            "scl_survey": ee.Filter.And(
                ee.Filter.gte("size", self.thresholds["landscape_size"]),
                ee.Filter.eq("range", self.thresholds["current_range"]),
                ee.Filter.gte("probability", self.thresholds["landscape_probability"]),
                ee.Filter.lt("effort", self.thresholds["landscape_survey_effort"]),
            ),
            "scl_fragment": ee.Filter.And(
                ee.Filter.lt("size", self.thresholds["landscape_size"]),
                ee.Filter.eq("range", self.thresholds["current_range"]),
                ee.Filter.gte("probability", self.thresholds["landscape_probability"]),
                ee.Filter.gte("effort", self.thresholds["landscape_survey_effort"]),
            ),
        }

    def distance_km_to_pixels(self, distance_km, image_resolution):
        image_resolution = ee.Number(image_resolution)
        scale_km = image_resolution.divide(1000)
        pixel_distance = ee.Number(distance_km).divide(scale_km).int()
        return pixel_distance

    def area_km_to_pixels(self, area_km, image_resolution):
        image_resolution = ee.Number(image_resolution)
        resolution_km = image_resolution.divide(1000)
        scale_km = resolution_km.pow(2)
        n_pixels = ee.Number(area_km).divide(scale_km).int()
        return n_pixels

    def density_to_patch_size(self, ft):
        med_density_eco = ee.Number(ft.get("MED_DENSITY_ECO"))
        med_density_biome = ee.Number(ft.get("MED_DENSITY_BIOME"))
        med_density_eco_biome = ee.Number(ft.get("MED_DENSITY_ECO_BIOME"))
        density_val = ee.Algorithms.If(
            med_density_eco.gt(0),
            med_density_eco,
            ee.Algorithms.If(
                med_density_biome.gt(0),
                med_density_biome,
                ee.Algorithms.If(med_density_eco_biome.gt(0), med_density_eco_biome, 1),
            ),
        )
        min_core_size = (
            ee.Number(self.density_values["n_core_animals"])
            .divide(density_val)
            .multiply(100)
            .int()
        )
        min_core_size_pixels = self.area_km_to_pixels(min_core_size, self.scale)
        return ee.Feature(
            None,
            {"ECO_ID": ft.get("ECO_ID"), "min_core_size": min_core_size_pixels},
        )

    def dilate(self, image, distance):
        pixel_distance = self.distance_km_to_pixels(distance, self.scale)
        dialated_image = image.fastDistanceTransform(pixel_distance).sqrt()
        return dialated_image.lte(pixel_distance).selfMask()

    def poly_export(self, polys, scl_name):
        size_test = polys.size().gt(0).getInfo()
        path = "pothab/" + scl_name
        if size_test:
            self.export_fc_ee(polys, path)
        else:
            print("no " + scl_name + " polygons delineated")

    def build_threshold_calc_image(
        self, threshold_variable_image, probability_image, variable_name
    ):
        variable_high_probability = threshold_variable_image.updateMask(
            probability_image.selfMask()
        )
        return (
            threshold_variable_image.addBands(
                [variable_high_probability, self.zones_image]
            )
            .updateMask(self.zones_image)
            .rename(
                [
                    f"{variable_name}_all_vals",
                    f"{variable_name}_high_prob_vals",
                    "zone",
                ]
            )
        )

    def histogram_calc(self, threshold_calc_image, hist_reducer):
        return ee.List(
            ee.Dictionary(
                threshold_calc_image.reduceRegion(
                    reducer=hist_reducer,
                    geometry=ee.Geometry.Polygon(self.extent),
                    scale=self.scale,
                    crs=self.crs,
                    maxPixels=self.ee_max_pixels,
                )
            ).get("groups")
        )

    def zone_threshold_calc(self, zone_histogram_objects):
        # histogram_format is mapped over zone_histogram_objects
        # returns a list of lists [[zone_number,...], [threshold,...]]
        def histogram_format(zone_object):
            zone_dictionary = ee.Dictionary(zone_object)
            zone_number = zone_dictionary.get("zone")
            pixel_counts = ee.List(zone_dictionary.get("count"))
            pixel_count_all = ee.Number(pixel_counts.get(0))
            pixel_count_high = ee.Number(pixel_counts.get(1))
            histograms = ee.List(zone_dictionary.get("histogram"))

            bin = ee.Array(histograms.get(0)).slice(1, 0, 1).toList().flatten()
            histogram_all_pixels = (
                ee.Array(histograms.get(0)).slice(1, 1, 2).toList().flatten()
            )
            histogram_high_probability_pixels = (
                ee.Array(histograms.get(1)).slice(1, 1, 2).toList().flatten()
            )
            histogram_combine = bin.zip(
                histogram_all_pixels.zip(histogram_high_probability_pixels)
            )
            # threshold_calc is mapped over hist_combine
            # each histogram object is a list of lists: [bin_value, frequency_all_pixels, frequency_high_probablity_pixels]
            # returns list of bin values where propotional difference is gte 0
            def threshold_calc(histogram_object):
                histogram_item_list = ee.List(histogram_object)
                bin_value = ee.Number(histogram_item_list.get(0))
                frequency_values = ee.List(histogram_item_list.get(1))
                frequency_all = ee.Number(frequency_values.get(0))
                frequency_high = ee.Number(frequency_values.get(1))
                proportion_all = frequency_all.divide(pixel_count_all)
                proportion_high = frequency_high.divide(pixel_count_high)
                difference = proportion_high.subtract(proportion_all)
                # if difference is gte to 0, returns bin value; if difference is < 0, returns 0
                return ee.Number(bin_value.multiply(difference.gte(0))).int()

            # returns the maximum bin value with a positive proportional difference
            threshold = histogram_combine.map(threshold_calc).reduce(ee.Reducer.max())
            return [zone_number, threshold]

        # converts [[zone_number, threshold],...] to [[zone_number,...], [threshold,...]]
        zone_threshold_list = (
            ee.Array(zone_histogram_objects.map(histogram_format)).transpose().toList()
        )
        # TODO: 1. identify and handle edge cases; 2. return default value if thresholding fails
        return zone_threshold_list

    def zone_threshold_image(self, zone_threshold_lists, variable_name):
        zone_numbers = zone_threshold_lists.get(0)
        threshold_vals = zone_threshold_lists.get(1)
        return self.zones_image.remap(zone_numbers, threshold_vals).rename(
            f"{variable_name}_thresholds"
        )

    def calc(self):
        str_hab_mask = self.structural_habitat.gte(
            self.thresholds["structural_habitat"]
        ).selfMask()

        str_hab_resolution = self.structural_habitat.projection().nominalScale()
        str_hab_connected_pixels = self.area_km_to_pixels(
            self.thresholds["structural_habitat_patch_size"], str_hab_resolution
        )

        connected_str_hab = (
            str_hab_mask.connectedPixelCount(str_hab_connected_pixels, True)
            .gte(str_hab_connected_pixels)
            .selfMask()
            .reproject(self.crs, None, str_hab_resolution)
        )

        hii_thresholding_image = self.build_threshold_calc_image(
            self.hii,
            str_hab_mask,  # TODO: replace str_hab_mask with probability when available
            "HII",
        )

        hii_thresholding_histograms = self.histogram_calc(
            hii_thresholding_image,
            self.histogram_reducer(0, 100, 100),
        )

        hii_threshold_values = self.zone_threshold_calc(hii_thresholding_histograms)

        hii_threshold_image = self.zone_threshold_image(hii_threshold_values, "HII")

        low_hii_mask = self.hii.lte(hii_threshold_image).selfMask()

        eff_pot_hab = (
            str_hab_mask.updateMask(low_hii_mask)
            .rename("eff_pot_hab")
            .unmask(0)
            .reduceResolution(ee.Reducer.mean())
            .gte(self.thresholds["reduce_res_input_pixels"])
            .selfMask()
            .rename("eff_pot_hab")
        )

        eff_pot_hab_export = (
            connected_str_hab.updateMask(low_hii_mask)
            .unmask(0)
            .reduceResolution(ee.Reducer.mean())
            .gte(self.thresholds["reduce_res_input_pixels"])
            .selfMask()
            .rename("eff_pot_hab")
        )

        max_core_pixels = self.area_km_to_pixels(
            self.density_values["core_size_limits"]["max"], self.scale
        )

        min_core_pixels = self.area_km_to_pixels(
            self.density_values["core_size_limits"]["min"], self.scale
        )

        max_step_pixels = self.area_km_to_pixels(
            self.density_values["step_size_limits"]["max"], self.scale
        )

        min_step_pixels = self.area_km_to_pixels(
            self.density_values["step_size_limits"]["min"], self.scale
        )

        connected_potential_habitat = eff_pot_hab.connectedPixelCount(
            max_core_pixels, True
        ).reproject(crs=self.crs, scale=self.scale)

        potential_habitat = connected_potential_habitat.updateMask(
            connected_potential_habitat.gte(min_step_pixels)
        ).selfMask()

        density = self.density.map(self.density_to_patch_size)

        ecoregion_image = self.ecoregions.reduceToImage(
            properties=["ECO_ID"], reducer=ee.Reducer.mode()
        )

        ecoregion_id = density.aggregate_array("ECO_ID")

        core_size = density.aggregate_array("min_core_size")

        min_patch_size = (
            ecoregion_image.remap(ecoregion_id, core_size)
            .int()
            .clamp(
                min_core_pixels,
                max_core_pixels,
            )
            .updateMask(potential_habitat)
        )

        min_stepping_stone_size = (
            min_patch_size.multiply(self.density_values["core_to_step_ratio"])
            .int()
            .clamp(
                min_step_pixels,
                max_step_pixels,
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
            self.dilate(potential_core, self.thresholds["connectivity_distance"])
            .reproject(crs=self.crs, scale=self.scale)
            .multiply(3)
            .unmask(0)
            .updateMask(self.water)
        )
        potential_stepping_stone = (
            self.dilate(
                potential_stepping_stone, self.thresholds["connectivity_distance"]
            )
            .reproject(crs=self.crs, scale=self.scale)
            .unmask(0)
            .updateMask(self.water)
        )

        scl_polys = potential_core.add(potential_stepping_stone).selfMask()

        # TODO: the probability threshold process will change with final version of probability
        probability_thresh = (
            self.probability.add(1).gte(self.thresholds["probability"]).selfMask()
        )

        historic = self.historic_range.unmask(0)
        extirpated = self.extirpated_range.eq(0)

        range_binary = (
            ee.Image(0)
            .where(historic.eq(1), ee.Image(2))
            .where(extirpated.eq(1), ee.Image(1))
        ).selfMask()

        # TODO: incorporate actual survey effort input
        survey_effort = ee.Image.constant(1)

        scl_polys = (
            ee.Image(1)
            .updateMask(scl_polys)
            .addBands([scl_polys, range_binary, probability_thresh, survey_effort])
            .rename(["scl_poly", "size", "range", "probability", "effort"])
            .reduceToVectors(
                reducer=ee.Reducer.max(),  # TODO: may need to consider a unique reducer for each band to delineate polygons
                geometry=ee.Geometry.Polygon(self.extent),
                scale=self.scale,
                crs=self.crs,
                maxPixels=self.ee_max_pixels,
            )
        )

        scl_species = scl_polys.filter(self.scl_poly_filters["scl_species"])

        scl_survey = scl_polys.filter(self.scl_poly_filters["scl_survey"])

        scl_restoration = scl_polys.filter(self.scl_poly_filters["scl_restoration"])

        scl_fragment = scl_polys.filter(self.scl_poly_filters["scl_fragment"])

        self.poly_export(scl_species, "scl_species")
        self.poly_export(scl_survey, "scl_survey")
        self.poly_export(scl_restoration, "scl_restoration")
        self.poly_export(scl_fragment, "scl_fragment")

        self.export_image_ee(eff_pot_hab_export, "pothab/potential_habitat")

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
