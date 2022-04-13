import argparse
import ee
from task_base import SCLTask


class SCLEffectivePotentialHabitat(SCLTask):
    BIOME_ZONE_LABEL = "Zone"
    inputs = {
        "structural_habitat": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "structural_habitat_path",
            "maxage": 1,
        },
        "extirpated_range": {
            "ee_type": SCLTask.IMAGE,
            "ee_path": "extirpated_range_path",
            "static": True,
        },
        "density": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "density_path",
            "static": True,
        },
        "zones": {
            "ee_type": SCLTask.FEATURECOLLECTION,
            "ee_path": "zones_path",
            "static": True,
        },
        "hii": {
            "ee_type": SCLTask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/hii",
            "maxage": 1,
        },
    }
    thresholds = {
        "structural_habitat": 0.5,
        "reduce_res_input_pixels": 0.5,
        "structural_habitat_patch_size": 5,  # sq km
        "hii": {"zone_1": 14.4, "zone_2": 7.2, "zone_3": 4.9, "zone_4": 4.9},
        "dispersal_distance": 4,
    }
    density_values = {
        "n_core_animals": 5,
        "core_to_step_ratio": 0.1,
        "core_size_limits": {"min": 30, "max": 625},  # size min/max in sqkm
        "step_size_limits": {"min": 3, "max": 63},  # size min/max in sqkm
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.structural_habitat, _ = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["structural_habitat"]["ee_path"])
        )
        self.hii, _ = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["hii"]["ee_path"])
        )
        self.extirpated_range = (
            ee.FeatureCollection(self.inputs["extirpated_range"]["ee_path"])
            .filter(
                ee.Filter.And(
                    ee.Filter.lte("ext_year", self.taskdate.year),
                    ee.Filter.gt("ext_revert", self.taskdate.year),
                )
            )
            .reduceToImage(["diss"], ee.Reducer.first())
            .unmask(0)
        )
        self.density = ee.FeatureCollection(self.inputs["density"]["ee_path"])
        self.zones = ee.FeatureCollection(self.inputs["zones"]["ee_path"])

    def structural_habitat_path(self):
        return f"{self.ee_rootdir}/structural_habitat"

    def extirpated_range_path(self):
        return f"{self.speciesdir}/extirpated_range"

    def density_path(self):
        return f"{self.speciesdir}/biome_density/biome_density_2021-10-26"

    def zones_path(self):
        return f"{self.speciesdir}/zones"

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
        density_val = ee.Algorithms.If(
            med_density_eco.gt(0),
            med_density_eco,
            ee.Algorithms.If(med_density_biome.gt(0), med_density_biome, 1),
        )
        min_core_size = (
            ee.Number(self.density_values["n_core_animals"])
            .divide(density_val)
            .multiply(100)
            .int()
        )
        # min_core_size_pixels = self.area_km_to_pixels(min_core_size, self.scale)
        # return km2 required to be a patch
        return ee.Feature(
            None,
            {"ECO_ID": ft.get("ECO_ID"), "min_core_size": min_core_size},
        )

    def dilate(self, image, distance):
        pixel_distance = self.distance_km_to_pixels(distance, self.scale)
        dialated_image = image.fastDistanceTransform(pixel_distance).sqrt()
        return dialated_image.lte(pixel_distance).selfMask()

    def calc(self):
        structural_habitat_mask = self.structural_habitat.gte(
            self.thresholds["structural_habitat"]
        ).selfMask()

        str_hab_resolution = self.structural_habitat.projection().nominalScale()
        str_hab_connected_pixels = self.area_km_to_pixels(
            self.thresholds["structural_habitat_patch_size"], str_hab_resolution
        )
        connected_structural_habitat = (
            structural_habitat_mask.connectedPixelCount(str_hab_connected_pixels, True)
            .gte(str_hab_connected_pixels)
            .selfMask()
            .reproject(self.crs, None, str_hab_resolution)
        )

        zone_image = self.zones.reduceToImage(
            [self.BIOME_ZONE_LABEL], ee.Reducer.first()
        )

        hii_threshold_image = zone_image.remap(
            [1, 2, 3, 4],
            [
                self.thresholds["hii"]["zone_1"],
                self.thresholds["hii"]["zone_2"],
                self.thresholds["hii"]["zone_3"],
                self.thresholds["hii"]["zone_4"],
            ],
        )
        low_hii_mask = self.hii.divide(100).lte(hii_threshold_image).selfMask()

        eff_pot_hab = (
            structural_habitat_mask.updateMask(low_hii_mask)
            .unmask(0)
            .reduceResolution(ee.Reducer.mean())
            .gte(self.thresholds["reduce_res_input_pixels"])
            .selfMask()
            .rename("eff_pot_hab")
        )
        eff_pot_hab_export = (
            connected_structural_habitat.updateMask(low_hii_mask)
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

        country_image = self.countries.reduceToImage(
            properties=["ISONUMERIC"], reducer=ee.Reducer.mode()
        )
        ecoregion_image = self.ecoregions.reduceToImage(
            properties=["ECO_ID"], reducer=ee.Reducer.mode()
        )
        biome_image = self.ecoregions.reduceToImage(
            properties=["BIOME_NUM"], reducer=ee.Reducer.mode()
        )
        eco_country = country_image.multiply(1000).add(ecoregion_image)
        density = self.density.map(self.density_to_patch_size)
        ecoregion_id = density.aggregate_array("ECO_ID")
        core_size = density.aggregate_array("min_core_size")

        min_patch_size = (
            ecoregion_image.remap(ecoregion_id, core_size)
            .int()
            .clamp(
                min_core_pixels,
                max_core_pixels,
            )
        )
        min_stepping_stone_size = (
            min_patch_size.multiply(self.density_values["core_to_step_ratio"])
            .int()
            .clamp(
                min_step_pixels,
                max_step_pixels,
            )
        )

        connectivity_distance = self.thresholds["dispersal_distance"] / 2
        potential_core = (
            self.dilate(
                ee.Image(0).where(potential_habitat.gte(min_patch_size), 1).selfMask(),
                connectivity_distance,
            )
            .reproject(crs=self.crs, scale=self.scale)
            .multiply(3)
            .unmask(0)
            .updateMask(self.watermask)
        )
        potential_stepping_stone = (
            self.dilate(
                ee.Image(0)
                .where(
                    potential_habitat.lt(min_patch_size).And(
                        potential_habitat.gte(min_stepping_stone_size)
                    ),
                    1,
                )
                .selfMask(),
                connectivity_distance,
            )
            .reproject(crs=self.crs, scale=self.scale)
            .unmask(0)
            .updateMask(self.watermask)
        )

        allpotential = potential_core.add(potential_stepping_stone).selfMask()

        # range reclass:
        # 0: neither historic or extirpated range
        # 1: extirpated
        # 2: historical
        range_class = (
            ee.Image(0)
            .where(self.historical_range.eq(1), ee.Image(2))
            .where(self.extirpated_range.eq(1), ee.Image(1))
        ).selfMask()

        paimage = self.pas.reduceToImage(["WDPAID"], ee.Reducer.first())

        area = ee.Image.pixelArea().divide(1000000).updateMask(self.watermask)
        eff_pot_hab_area = area.updateMask(eff_pot_hab)
        potential_habitat_area = area.updateMask(potential_habitat)
        pa_area = area.updateMask(paimage)

        mode_bands = [
            "range",
            "country",
            "ecoregion",
            "biome",
            "min_patch_size",
            "min_stepping_stone_size",
        ]
        sum_bands = [
            "polygon_area",
            "eff_pot_hab_area",
            "connected_eff_pot_hab_area",
            "pa_area",
        ]

        scl_image = (
            eco_country.updateMask(allpotential)
            .addBands(
                [
                    range_class,
                    country_image,
                    ecoregion_image,
                    biome_image,
                    min_patch_size,
                    min_stepping_stone_size,
                    area,
                    eff_pot_hab_area,
                    potential_habitat_area,
                    pa_area,
                ]
            )
            .rename(["scl_poly"] + mode_bands + sum_bands)
        )
        scl_polys = scl_image.reduceToVectors(
            reducer=ee.Reducer.mode()
            .forEach(mode_bands)
            .combine(ee.Reducer.sum().forEach(sum_bands)),
            geometry=ee.Geometry.Polygon(self.extent),
            scale=self.scale,
            crs=self.crs,
            maxPixels=self.ee_max_pixels,
            tileScale=8,
        ).filter(ee.Filter.gt("eff_pot_hab_area", 0))

        def _attribute(item):
            item = ee.List(item)
            feature = ee.Feature(item.get(0))
            poly_id = ee.Number(item.get(1)).int()

            pa_proportion = ee.Number(feature.get("pa_area")).divide(
                ee.Number(feature.get("polygon_area"))
            )

            feature = feature.set({"poly_id": poly_id, "pa_proportion": pa_proportion})
            return_properties = feature.propertyNames().filter(
                ee.Filter.neq("item", "pa_area")
            )
            return feature.select(return_properties)

        ids = ee.List.sequence(1, scl_polys.size())
        scl_poly_list = ee.List(scl_polys.toList(scl_polys.size()))
        scl_polys_assigned = ee.FeatureCollection(
            scl_poly_list.zip(ids).map(_attribute)
        )

        self.export_image_ee(eff_pot_hab_export, "pothab/potential_habitat")
        self.export_fc_ee(scl_polys_assigned, "pothab/scl_polys")
        self.export_image_ee(scl_image, "pothab/scl_image")

    def check_inputs(self):
        super().check_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate")
    parser.add_argument("-s", "--species")
    parser.add_argument("--scenario")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing outputs instead of incrementing",
    )
    options = parser.parse_args()
    effpothab_task = SCLEffectivePotentialHabitat(**vars(options))
    effpothab_task.run()
