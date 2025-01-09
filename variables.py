""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ AVASPEC                                                                  │
  └──────────────────────────────────────────────────────────────────────────┘
"""

dev_handle = 0
pixels = 4096
wavelength = [0.0] * 4096
spectraldata = [0.0] * 4096
NrScanned = 0
lasttime = 0
wavelength_sel = 0
cont = 0
inttime = 0.1
averages = 10
fiber = 400e-6

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ MMS                                                                      │
  └──────────────────────────────────────────────────────────────────────────┘
"""

mms_serial = None
move_speed = 10000
travel_speed = 10000

x0abs, y0abs, z0abs = 0,0,0

WAIT = False

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ PROBHOLDER                                                               │
  ├──────────────────────────────────────────────────────────────────────────┤
  │   ╭────────────────────────────────────────────╮                         │
  │   │ ╭───────╮  ╭───────╮  ╭───────╮  ╭───────╮ │                         │
  │   │ │ [1,2] │  │ [2,2] │  │ [3,2] │  │ [4,2] │ │                         │
  │   │ │       │  │       │  │       │  │       │ │                         │
  │   │ ╰───────╯  ╰───────╯  ╰───────╯  ╰───────╯ │                         │
  │   │                                            │                         │
  │   │ ╭───────╮  ╭───────╮  ╭───────╮  ╭───────╮ │                         │
  │   │ │ [1,1] │  │ [2,1] │  │ [3,1] │  │ [4,1] │ │                         │
  │   │ │       │  │       │  │       │  │       │ │                         │
  │   │ ╰───────╯  ╰───────╯  ╰───────╯  ╰───────╯ │                         │
  │   ╰────────────────────────────────────────────╯                         │
  │                                                                          │
  │ Default array for all positions:                                         │
  │                                                                          │
  │ probholder_location = [[1,1],[2,1],[3,1],[4,1],[1,2],[2,2],[3,2],[4,2]]  │
  └──────────────────────────────────────────────────────────────────────────┘
"""

probholder_location = [[1,1]]

probholder_xdelta = 30
probholder_ydelta = 40

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ MEASUREMENT                                                              │
  └──────────────────────────────────────────────────────────────────────────┘
"""

steps = None
xprob, yprob = 22.5, 22.5
scan_axis = "yx"
scan_grid, move_grid = 0.05, 0.25

motion_time = -1
start_time = 0
log_time = 1

BACK_TO_ORIGIN = True
PARK_HEAD = False
LOG = False
STOPWATCH = False
SPECTRUM = False

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ EVALUATION                                                               │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ NEEDED: h5py (pip install h5py)                                          │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ Current structure of hdf5 file (Add new elements into the list):         │
  │                                                                          │
  │ data_{file_time}.h5     (hdf5 file)                                      │
  │ ├─── position           (hdf5 group)                                     │
  │ │    ├─── X             (1D np.array)                                    │
  │ │    ├─── Y             (1D np.array)                                    │
  │ │    ├─── Z             (1D np.array)                                    │
  │ │    └─── meas_time     (1D np.array)                                    │
  │ └─── spectrum           (hdf5 group)                                     │
  │      ├─── wavelength    (1D np.array)                                    │
  │      ├─── timestamp     (1D np.array)                                    │
  │      └─── intensity     (2D np.array)                                    │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ For example, to access position data:                                    │
  │                                                                          │
  │ import h5py                                                              │
  │ f = h5py.File(data_{file_time}.h5, "r")                                  │
  │                                                                          │
  │ x_coord = f["position/X"][:]                                             │
  │                                                                          │
  │ # Close hdf5 file !IMPORTANT! otherwise data gets currupted              │
  │ f.close()                                                                │
  ├──────────────────────────────────────────────────────────────────────────┤
  │ To save new data (dict) change following parameter accordingly:          │
  │                                                                          │
  │ Add these elements                      │                                │
  │                                         V                                │
  │ data_group = ["position", "spectrum", "YOUR_DATA"]                       │
  │ data       = [None,       None,       None       ]                       │
  │                                                                          │
  │ And save data into variable: variables.data[index of YOUR_DATA]          │
  └──────────────────────────────────────────────────────────────────────────┘
"""

path = r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\mai\Filme\2024-05-24\11-13-03\data_11-13-03.h5' #ändern

pic_file = ["intensity_map", "intensity_xy"]

file_data = 0
file_time = 0

data_file  = path
data_group = ["position", "spectrum"]
data       = [None,        None     ]

interpolation = 'quadric'
wavelength_eval = 590

EVAL = True
SHOW_PLOT = True
SAVE_PLOT = False
INTENSITY_PLOT = False