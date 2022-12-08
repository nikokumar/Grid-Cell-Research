# Dataset

## Dataset description:
Long Evans rats were implanted in the medial entorhinal cortex with neuropixels probes and recorded during a random foraging open field protocol in darkness (task 1) then in light (task 2).

## Dataset file names:
- `moserlab_waaga_25843_2019-09-13_22-54-22_v1.npy`
- `moserlab_waaga_26018_2019-12-10_15-25-47_v1.npy`
- `moserlab_waaga_26018_2019-12-14_16-03-44_v1.npy`
- `moserlab_waaga_26718_2020-09-16_17-23-51_v1.npy`
- `moserlab_waaga_26820_2020-11-05_11-03-13_v1.npy`

## Dataset file format:
`.npy (python numpy)`

Each file contains tracked rat behaviour and all spiketimes from cells recognized as grid cells (see methods of paper).

---

## data loading:
The dataset for each rat is stored in one single python numpy `.npy` file.

```python
# load data from rat #26018
import numpy
dataset = np.load("moserlab_waaga_26018_2019-12-10_15-25-47_v1.npy", allow_pickle=True).item()
```

## dataset organisation:

```
dataset = dict with keys:

              name: <string>      name of dataset
              task: <1×n dict>    data for each task
       description: <string>      description of dataset
           unit_id: <n×1 int>     cluster id of each cell
         module_id: <n×1 int>     grid module id of each cell


dataset["task"][0] = dict with keys:

              name: <string>      name of task
   spike_timestamp: <n×1 float>   all spiketimes for all cells
  spike_cluster_id: <n×1 int>     cluster index of each spike
          tracking: <1×1 dict>    positions of the rat
    tracking_units: <string>      tracking unit of measurement


dataset["task"][0]["tracking"] = dict with keys:

                 t: <n×1 float>   tracking timestamps
                 x: <n×1 float>   x-position
                 y: <n×1 float>   y-position
                 z: <n×1 float>   z-position
                hd: <n×1 float>   azimuth head direction (radians)
```
Each entry in the `spike_timestamp` array has a corresponding entry in the `spike_cluster_id` array which maps the spiketime to a specific cell id. i.e. the spikes times of cell `42` are located at the index of `spike_cluster_id` with value `42`.

```python
# get spiketimes of cluster 42
darkness = 0
st = dataset["task"][darkness]["spike_timestamp"]
cid = dataset["task"][darkness]["spike_cluster_id"]

cell_42_spiketimes = st[cid == 42]
```

---

Torgeir Waaga,
Moser lab,
Kavli Institute for Systems Neuroscience,
NTNU, Trondheim, Norway.