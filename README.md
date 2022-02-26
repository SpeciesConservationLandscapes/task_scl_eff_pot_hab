Species Conservation Landscapes effective potential habitat task
----------------------------------------------------------------

Task for calculating effective potential habitat for a species on a date. Structural habitat, produced by 
[task_scl_structural_habitat](https://github.com/SpeciesConservationLandscapes/task_scl_structural_habitat)
and defined as areas that are habitat from a remote sensing point of view, is masked by the 
[Human Impact Index](https://github.com/SpeciesConservationLandscapes/task_hii_weightedsum)
using thresholds that are determined by year and by species "zone". Then a filter for patch size is applied, using
a "stepping stone" logic that grows habitat areas according to species dispersal distance. Finally, the task
outputs an image for effective potential habitat itself, an image for species conservation landscapes with 
many of the bands necessary for classification by
[task_scl_classification](https://github.com/SpeciesConservationLandscapes/task_scl_classification), 
and a vectorized version of the polygons ready to be classified.

## Usage

*All parameters may be specified in the environment as well as the command line.*

```
/app # python task.py --help
usage: task.py [-h] [-d TASKDATE] [-s SPECIES] [--scenario SCENARIO] [--overwrite]

optional arguments:
  -h, --help            show this help message and exit
  -d TASKDATE, --taskdate TASKDATE
  -s SPECIES, --species SPECIES
  --scenario SCENARIO
  --overwrite           overwrite existing outputs instead of incrementing
```

### License
Copyright (C) 2022 Wildlife Conservation Society
The files in this repository  are part of the task framework for calculating 
Human Impact Index and Species Conservation Landscapes (https://github.com/SpeciesConservationLandscapes) 
and are released under the GPL license:
https://www.gnu.org/licenses/#GPL
See [LICENSE](./LICENSE) for details.
