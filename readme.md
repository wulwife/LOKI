# LOKI

LOKI (**LO**cation of seismic events through traveltime sta**KI**ng)
is a code that performs earthquake detection and location
using waveform coherence analysis (waveform stacking).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. You should have received a copy of
the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.


### Installation
It is always good habit to create (or work) in environments.
For more information please see [conda](https://docs.conda.io/en/latest/) or [pipenv](https://docs.pipenv.org/)


For example:

```bash
$ conda create -n loki python=3.7
$ conda activate loki
```

To use LOKI you will need:

- `GCC6` or later versions or any C compiler that supports OpenMP calls
- Python3 libraries, `numpy` (needed also at installation stage) and `obspy`
- Python3 library: `LatLonUTMconverter`, provided here but developed by third parts. (soon this dependence will be removed by adding an internal function)
 From `v1.0.0` the main module has been incorporated in loki (`from loki import LatLonUTMconverter`)


The software now comes with an installer!  simply digit:

```bash
$ # activate the envionment (optional)
$ pip install numpy  # v1.18 or higher (needed for C-compilation)
$ cd WHERE_LOKI_IS_STORED
$ pip install .
```

_NB: before running the `pip install` command, make sure to change accordingly the `os.environ["CC"]` variable
with a compiler that supports OpenMP calls._


### Testing
If installation didn't throw any error, you should be fine.
If you still want to double-check, the ultimiate installation, you could find a `tests` dir with an executable script. Run it anc compare the output results with the expected_results dir.
Results will be equals up to significant float digits.


### Usage
To use LOKI you could use the carbon-print `RunLoki.py` script (in `bin` folder).
You could copy it anywhere in the system, just make sure to make it executable.
Open the file to change the configuration's options at your preference.


### Citing
Please cite the following articles in documents showing
outputs of LOKI:

For _location_:

- Grigoli et al. 2013,
Automated seismic event location via traveltime stacking:
An application to mining induced seismicity,
Seismological Research Letters,
`https://doi.org/10.1785/0220120191`

- Grigoli et al 2014,
Automated seismic event location by waveform coherence analysis,
Geophysical Journal international,
`https://doi.org/10.1093/gji/ggt477`

- Grigoli et al. 2016,
Automated microseismic event location using Master-Event Waveform Stacking,
Scientific Reports,
`https://doi.org/10.1038/srep25744`


For _detection_:

- Grigoli et al. 2018,
Pick-and waveform-based techniques for real-time detection of induced seismicity,
Geophysical Journal International,
`https://doi.org/10.1093/gji/ggy019`

#### Author contact details:
Copyright(C) 2013 Francesco Grigoli

Author: Francesco Grigoli

Mail: <fsco.grigoli@gmail.com> or <francesco.grigoli@sed.ethz.ch>


This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
