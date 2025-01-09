'''
AUTHOR: Manuel Lippert (GitHub: ManeLippert)
LAST EDIT: 08.06.2024

DESCRIPTION:
This python file contains variables to set up initially parameter for MMS and change + store them dynamically.
Changes should be made at own risk.
'''



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
switch_speed = 1000

x0abs, y0abs, z0abs = 0,0,0

WAIT = False

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ CARRIER                                                                  │
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
  │ prob_location = [[1,1],[2,1],[3,1],[4,1],[1,2],[2,2],[3,2],[4,2]]        │
  └──────────────────────────────────────────────────────────────────────────┘
"""

prob_location   = [[1,1],[2,1],[3,1],[4,1],[1,2],[2,2],[3,2],[4,2]]
prob_input_file = [["./input_sample.json" for i in range(8)]]

carrier_xdelta = 30
carrier_ydelta = 40
carrier_zdelta = 4

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ CARRIERBOX                                                               │
  ├──────────────────────────────────────────────────────────────────────────┤
  │   View from above of MMS                                                 │  
  │                                                                          │  
  │   ╭────────────────╮               ╭────────────────╮                    │
  │   │╭──╮╭──╮╭──╮╭──╮│               │╭──╮╭──╮╭──╮╭──╮│                    │
  │   │╰──╯╰──╯╰──╯╰──╯│               │╰──╯╰──╯╰──╯╰──╯│                    │
  │   │╭──╮╭──╮╭──╮╭──╮│               │╭──╮╭──╮╭──╮╭──╮│                    │
  │   │╰──╯╰──╯╰──╯╰──╯│               │╰──╯╰──╯╰──╯╰──╯│                    │
  │   ╰────────────────╯               ╰────────────────╯                    │
  │                                                                          │  
  │   Measurement Stack (Left)         Store Stack (Right)                   │
  └──────────────────────────────────────────────────────────────────────────┘
"""

carrierbox_max_num  = 5

carrierbox_meas_num  = 3

carrierbox_meas_xabs = 84.5
carrierbox_meas_yabs = 390.2
carrierbox_meas_zabs = 2

carrierbox_store_num = 0

carrierbox_store_xabs = 354
carrierbox_store_yabs = 390.8
carrierbox_store_zabs = 2

carrierbox_save_zabs = 50

# Possible states: "open", "tight", "close"
SWITCH_STATE = "open"

erel_close    = 9.4
erel_tight    = 2.9

""" 
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ MEASUREMENT                                                              │
  └──────────────────────────────────────────────────────────────────────────┘
"""

prob_name = "prob"
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
  │ variables.data.update(YOUR_DATA)                                         │
  │                                                                          │
  │ Dictionary with YOUR_DATA should have the following form                 │
  │                                                                          │
  │ dic = {                                                                  │
  │     "group_name": {                                                      │
  │         "data1" : [],                                                    │
  │         "data2" : []                                                     │
  │     }                                                                    │
  │ }                                                                        │
  └──────────────────────────────────────────────────────────────────────────┘
"""

path = r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Juli\2024-06-19\M17\data_M17_13-14-37.h5'

pic_file = ["intensity_map", "intensity_xy"]

file_date = 0
file_time = 0

data_file  = "data"
data       = {}

interpolation = 'quadric'
wavelength_eval = 590

EVAL = True
SHOW_PLOT = True
SAVE_PLOT = True
INTENSITY_PLOT = False
AUTOMATIC_PATH = False