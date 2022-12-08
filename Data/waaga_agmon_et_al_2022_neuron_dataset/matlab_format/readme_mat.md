# Dataset

## Dataset description:
Long Evans rats were implanted in the medial entorhinal cortex with neuropixels probes and recorded during a random foraging open field protocol in darkness (task 1) then in light (task 2).

## Dataset file names:
- `moserlab_waaga_25843_2019-09-13_22-54-22_v1.mat`
- `moserlab_waaga_26018_2019-12-10_15-25-47_v1.mat`
- `moserlab_waaga_26018_2019-12-14_16-03-44_v1.mat`
- `moserlab_waaga_26718_2020-09-16_17-23-51_v1.mat`
- `moserlab_waaga_26820_2020-11-05_11-03-13_v1.mat`

## Dataset file format:
`.mat (Matlab .MAT v.7)`

Each file contains tracked rat behaviour and all spiketimes from cells recognized as grid cells (see methods of paper).

---

## data loading:
The dataset for each rat is stored in one single matlab `.mat` file.

```matlab
% Load data from rat #26018
dataset = load("moserlab_waaga_26018_2019-12-10_15-25-47_v1.mat")
```

## dataset organisation:
```
dataset = struct with fields:

              name: <string>       name of dataset
              task: <1×n struct>   data for each task
       description: <string>       description of dataset
           unit_id: <n×1 double>   cluster id of each cell
         module_id: <n×1 double>   grid module id of each cell


dataset.task(i) = struct with fields:

              name: <string>      name of task
   spike_timestamp: <n×1 double>  all spiketimes for all cells
  spike_cluster_id: <n×1 double>  cluster index of each spike
          tracking: <1×1 struct>  positions of the rat
    tracking_units: <string>      tracking unit of measurement


dataset.task(i).tracking = struct with fields:

         timestamp: <n×1 double>  tracking timestamps
                 x: <n×1 double>  x-position
                 y: <n×1 double>  y-position
                 z: <n×1 double>  z-position
    head_direction: <n×1 double>  azimuth head direction (radians)
```
Each entry in the `spike_timestamp` array has a corresponding entry in the `spike_cluster_id` array which maps the spiketime to a specific cell id. i.e. the spikes times of cell `42` are located at the index of `spike_cluster_id` with value `42`.

```matlab
% get spiketimes of cluster 42
darkness = 1

st = dataset.task(darkness).spike_timestamp;
cid = dataset.task(darkness).spike_cluster_id;

cell_42_spiketimes = st(cid == 42);
```

---

Torgeir Waaga,
Moser lab,
Kavli Institute for Systems Neuroscience,
NTNU, Trondheim, Norway.