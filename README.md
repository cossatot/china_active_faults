# China Active Faults

This is a repository containing the data and script necessary to run a joint 
geodetic/geologic block inversion to get slip rates for *HimaTibetMap v.2.0*, 
and the results (as of the most recent commit).

## Running the inversion

The block inversion uses [`Oiler`](github.com/cossatot/Oiler). Please see that 
repository for installation instructions.

Once `Oiler` has been installed, go into the `scripts/` directory, start a 
`julia` interpreter in the terminal, and type
```julia
include chn_blocks.jl
```

First, the code will compile, then the inversion will run.


## Accessing the results

The results are stored in the 
[`results`](https://github.com/cossatot/china_active_faults/tree/master/results) 
directory. There are three files:

- `block_vels.csv`, which is a CSV file that contains the predicted block 
  velocities relative to stable Eurasia. The fields are `fid` (the block 
  index), `lon` and `lat` (the coordinates of the block centroid), and the 
  velocity components with 1-sigma uncertainties, `ve`, `vn`, `ee`, `en`, and 
  the covariance `cen`.

- `chn_gnss_results.csv`, which is a CSV file that contains the predicted GNSS 
  velocities (which are the sum of the block motions and earthquake cycle 
  effects). The fields are `lon` and `lat` (the geographic coordinates of each 
  GNSS station), `name`, `fix` (the index of the reference frame, 1111 for all, 
  indicating stable Eurasia, `mov` (the `fid` of the block), and a number of 
  fields with self-explanatory names denoting the observed, predicted and 
  residual GNSS velocities with 1-sigma uncertainties.

- `chn_faults_out.geojson`, which is a LineString (polyline) GIS file with the 
  fault traces, additional geometric information, estimated slip rates, and 
  some other metadata.


The `chn_faults_out.geojson` file has the following attributes:

- `fid`: the ID of the feature (the fault trace).
- `dip`: the average dip of the fault. Note that 89 is the max, so that `hw` 
  and `fw` are defined.
- `dip_dir`: the cardinal direction the fault dips in.
- `usd`: the upper seismogenic locking depth (km) in the inverison.
- `lsd`: the lower seismogenic locking depth (km) in the inversion.
- `hw`: the `fid` of the block in the hanging wall.
- `fw`: the `fid` of the block in the footwall.
- `name`: the name of the fault, if known (may be empty).
- `dextral_rate`: the dextral slip rate of the fault, in mm/yr. Sinistral 
  negative.
- `dextral_err`: the 1-s.d. uncertainty of the dextral slip rate, in mm/yr.
- `extension_rate`: the extensional slip rate of the fault, in mm/yr. 
  Contraction is negative.  Note that this is also the horizontal extension 
  rate across the fault.
- `extension_err`: the 1-s.d. uncertainty of the extensional slip rate, in 
  mm/yr.
- `cde`: the covariance of the `dextral_err` and `extension_err`.
- `net_slip_rate`: the vector magnitude of the slip rate, in mm/yr.
- `net_slip_rate_err`: the 1-s.d. uncertainty of the `net_slip_rate`.

## Visualizing the block motions with the interactive web viewer

The best way to understand the block results is to use an interactive viewer 
that is provided. This uses a web browser and requires the user to start a web 
server; if the user has Python installed, this should work without further 
configuration.

In a terminal, go into the `web_viewer/` and type

```bash
python -m http.server
```

and a web server should start from this directory, serving pages at 
`localhost:8000`. Then, type [`localhost:8000`](localhost:8000) in the browser 
to see the page.

A blank globe (latitude and longitude lines) should pop up, and there should be 
a button that says 'Draw Blocks', and a slider that says 'Ma'. Click on 'Draw 
Blocks', and after a little while (a minute or two, for the first time), the 
block model should appear on the globe. You can click on it and move it around, 
zoom in and out with the mouse, etc. Perhaps more importantly, you can move the 
'Ma' slider back and forth to simulate *instantaneous* block motions projected 
forward or backward 5 million years in time. *Please note this is not meant to 
be a real reconstruction or project, but simply a way to visualize the current 
block motions better!*
