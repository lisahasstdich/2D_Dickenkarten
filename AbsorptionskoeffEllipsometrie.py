import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.optimize import minimize_scalar
import pandas as pd
import variables, h5py
import variablesRef, h5py

measured_thickness = 58.21 # in nm angeben
x_pos =  17.55 #in mm angeben
y_pos = 11.5 #in mm angeben


class Absorption:
    def Initializing(self):
        self.intensity_array = None
        self.intensityRef_array = None
        
    def intensities(self,h5file = None, interpolation = None):
       h5file        = h5file        or f'{variablesRef.path}'
       interpolation = interpolation or variables.interpolation
       
       g = h5py.File(h5file, "r")
       wavelength_min = np.min(g['spectrum/wavelength']) #Wellenlängen, die man ausgewertet haben möchte
       wavelength_max = np.max(g['spectrum/wavelength']) #mit dieser Einstellung hat man ALLE Wellenlängen, kann aber wg teilweise positivem Alpha zu RuntimeWarnings kommen
       
       diff_wavelength = wavelength_max-wavelength_min
       wavelengths_array = g['spectrum/wavelength'][:]
       
       g.close()
       
       def Int0(self, h5file = None, scan_axis = None,
                             scan_grid = None, move_grid = None,
                             interpolation = None):
             def minimize_shift(x,n_move,n_scan,timestamps,t_array):
                 num_stripes=10
                 interv=int(round(n_move/(num_stripes+1)))
                 
                 diff=0
                 for j in range(1,num_stripes):
                     dist=int(round(n_move/num_stripes))
                     for i in range(n_scan):
                         m=np.searchsorted(timestamps-x*0.01,t_array[j*dist,i])
                         n=np.searchsorted(timestamps-x*0.01,t_array[j*dist+1,i])
                         
                         diff=diff+abs(intensityRef[m,0]-intensityRef[n,0])
                 
                 return diff
             # Initialize varibales from variablesRef.py
             h5file        = h5file        or f'{variablesRef.path}'
             g = h5py.File(h5file, 'r')
             try:
                 scan_grid     = scan_grid     or g['measurement/scan/scan_grid'][()]                
                 move_grid     = move_grid     or g['measurement/scan/move_grid'][()]                
                 scan_axis     = scan_axis     or g['measurement/scan/scan_axis'][()].decode('utf-8')
             except KeyError:
                 scan_grid     = scan_grid     or variablesRef.scan_grid
                 move_grid     = move_grid     or variablesRef.move_grid
                 scan_axis     = scan_axis     or variablesRef.scan_axis
             
             # Initialize data for intensityRef
             timestamps = g['spectrum/timestamp'][:]
             wavelength_array = g['spectrum/wavelength']
             
             wavelength_length = np.zeros_like(wavelength_array)

             # Add wavelength_length to wavelength_array
             wavelength_array += wavelength_length
             wavelength_values = g['spectrum/wavelength'][:]
            
             
             sorted_indices = np.argsort(np.abs(wavelength_values))
             # index_min= np.searchsorted(wavelength_values, wavelength_min)
             # index_max= np.searchsorted(wavelength_values, wavelength_max)
             # wavelength_array=wavelength_array[index_min:index_max]
             
             indices_wave = np.searchsorted(wavelength_values, wavelength_array)
             
             
             
             intensityRef = g['spectrum/intensity'][:,indices_wave]
             
             
             x, y, z, t = g['position/X'][:], g['position/Y'][:], g['position/Z'][:], g['position/meas_time'][:]
             xmin, xmax, ymin, ymax, zmin, zmax = x.min(), x.max(), y.min(), y.max(), z.min(), z.max()                           
             
             
              # Switch evaluation according the move direction
             if scan_axis == 'yx':
                  scan_dir, scan_min, scan_max, scan_label = y, ymin, ymax, 'y'   
                  move_dir, move_min, move_max, move_label = x, xmin, xmax, 'x'
          
             if scan_axis == 'xy':
                  scan_dir, scan_min, scan_max, scan_label = x, xmin, xmax, 'x'
                  move_dir, move_min, move_max, move_label = y, ymin, ymax, 'y'
             
              # Amount of steps for given grid size
             n_move = int(round((move_max-move_min)/move_grid+1)-1)
             n_scan = int(round((scan_max-scan_min)/scan_grid+1))
             
             scan_ax=np.linspace(scan_min,scan_max,n_scan)
             move_ax=np.linspace(move_min,move_max,n_move)
         
             scan_indices=np.searchsorted(scan_ax,scan_dir)
             move_indices=np.searchsorted(move_ax,move_dir)
             
             single=np.full(len(scan_indices),True)
             
             for i in range(1,len(scan_indices)):
                  if scan_indices[i]==scan_indices[i-1] and move_indices[i]==move_indices[i-1]: 
                      single[i]=False
                      single[i-1]=False
         
             filled=np.full((n_move,n_scan),-1.0)
             for i in range(len(scan_indices)):
                  filled[move_indices[i],scan_indices[i]]=t[i]
         
             indices=np.empty((n_move,n_scan))
             n=0
             
             for j in range(n_move):
                  try:
                      line=filled[j]
                      sel=line>-1
                      line_sel=line[sel]
                      if line_sel[1]>line_sel[0]:
                          for i in range(n_scan):
                              indices[j,i]=n
                              n=n+1
                      else:
                          for i in reversed(range(n_scan)):
                              indices[j,i]=n
                              n=n+1
                  except IndexError:
                      pass
         
             move_array, scan_array = np.indices((n_move,n_scan))
             scan_flat=np.ravel(scan_array)
             move_flat=np.ravel(move_array)
             filledflat=np.ravel(filled)
             indflat=np.ravel(indices)
             sort=np.argsort(indflat)
         
             scan_sorted=scan_flat[sort]
             move_sorted=move_flat[sort]
             indsorted=indflat[sort]
             filledsorted=filledflat[sort]
             sel=filledsorted>-1
             t_sel=filledsorted[sel]
             ind_sel=indsorted[sel]
             t_all=np.interp(indsorted,ind_sel,t_sel)
             t_array=np.empty((n_move,n_scan))
             for j in range(len(t_all)):
                  t_array[move_sorted[j],scan_sorted[j]]=t_all[j]                 
             # Calculate shift between different timestamps
             
             
             res = minimize_scalar(minimize_shift,args=(n_move,n_scan,timestamps,t_array))
             shift=res.x*0.01
             
             print("Shift: " + str(shift))
             # nur test
             #shift = -0.0145837
             intensityRef_array=np.empty((n_scan,n_move,len(indices_wave)))
             indices_wave=indices_wave[np.newaxis, np.newaxis, :]
             intensityRef_array= intensityRef_array + indices_wave
             
             for j in range(n_move):
                    for i in range(n_scan):
                            m=np.searchsorted(timestamps-shift,t_array[j,i])
                            try:
                                intensityRef_array[i,j,:]=intensityRef[m,:]
                            except IndexError:
                                    pass
             
             
             
               # Create 2D image of measurement
             
             
             g.close()
             return intensityRef_array
       def Int1(self, h5file = None, scan_axis =None,
                             scan_grid = None, move_grid = None,
                             interpolation = None):
             def minimize_shift(x,n_move,n_scan,timestamps,t_array):
                 num_stripes=10
                 interv=int(round(n_move/(num_stripes+1)))
                 
                 diff=0
                 for j in range(1,num_stripes):
                     dist=int(round(n_move/num_stripes))
                     for i in range(n_scan):
                         m=np.searchsorted(timestamps-x*0.01,t_array[j*dist,i])
                         n=np.searchsorted(timestamps-x*0.01,t_array[j*dist+1,i])
                         
                         diff=diff+abs(intensity[m,0]-intensity[n,0])
                 
                 return diff
             # Initialize varibales from variablesRef.py
             h5file        = h5file        or f'{variables.path}'
             
             g = h5py.File(h5file, 'r')
             try:
                 scan_grid     = scan_grid     or g['measurement/scan/scan_grid'][()]                
                 move_grid     = move_grid     or g['measurement/scan/move_grid'][()]                
                 scan_axis     = scan_axis     or g['measurement/scan/scan_axis'][()].decode('utf-8')
             except KeyError:
                 scan_grid     = scan_grid     or variables.scan_grid
                 move_grid     = move_grid     or variables.move_grid
                 scan_axis     = scan_axis     or variables.scan_axis
             
            
             # Initialize data for intensityRef
             timestamps = g['spectrum/timestamp'][:]
             wavelength_array = g['spectrum/wavelength']
             
             wavelength_length = np.zeros_like(wavelength_array)

             # Add wavelength_length to wavelength_array
             wavelength_array += wavelength_length
             
             wavelength_values = (g['spectrum/wavelength'][:])
             
             sorted_indices = np.argsort(np.abs(wavelength_values))
             # index_min= np.searchsorted(wavelength_values, wavelength_min)
             # index_max= np.searchsorted(wavelength_values, wavelength_max)
             # wavelength_array=wavelength_array[index_min:index_max]
             indices_wave = np.searchsorted(wavelength_values, wavelength_array)
             
            
             
             intensity =  g['spectrum/intensity'][:,indices_wave]
             
             
             x, y, z, t = g['position/X'][:], g['position/Y'][:], g['position/Z'][:], g['position/meas_time'][:]
             xmin, xmax, ymin, ymax, zmin, zmax = x.min(), x.max(), y.min(), y.max(), z.min(), z.max()                            
            
            
             # Switch evaluation according the move direction
             if scan_axis == 'yx':
                scan_dir, scan_min, scan_max, scan_label = y, ymin, ymax, 'y'   
                move_dir, move_min, move_max, move_label = x, xmin, xmax, 'x'
        
             if scan_axis == 'xy':
                scan_dir, scan_min, scan_max, scan_label = x, xmin, xmax, 'x'
                move_dir, move_min, move_max, move_label = y, ymin, ymax, 'y'
            
             # Amount of steps for given grid size
             n_move = int(round((move_max-move_min)/move_grid+1)-1)
             n_scan = int(round((scan_max-scan_min)/scan_grid+1))
            
             scan_ax=np.linspace(scan_min,scan_max,n_scan)
             move_ax=np.linspace(move_min,move_max,n_move)
        
             scan_indices=np.searchsorted(scan_ax,scan_dir)
             move_indices=np.searchsorted(move_ax,move_dir)
            
             single=np.full(len(scan_indices),True)
            
             for i in range(1,len(scan_indices)):
                 if scan_indices[i]==scan_indices[i-1] and move_indices[i]==move_indices[i-1]: 
                     single[i]=False
                     single[i-1]=False
        
             filled=np.full((n_move,n_scan),-1.0)
             for i in range(len(scan_indices)):
                 filled[move_indices[i],scan_indices[i]]=t[i]
        
             indices=np.empty((n_move,n_scan))
             n=0
            
             for j in range(n_move):
                 try:
                     line=filled[j]
                     sel=line>-1
                     line_sel=line[sel]
                     if line_sel[1]>line_sel[0]:
                         for i in range(n_scan):
                             indices[j,i]=n
                             n=n+1
                     else:
                         for i in reversed(range(n_scan)):
                             indices[j,i]=n
                             n=n+1
                 except IndexError:
                     pass
         
             move_array, scan_array = np.indices((n_move,n_scan))
             scan_flat=np.ravel(scan_array)
             move_flat=np.ravel(move_array)
             filledflat=np.ravel(filled)
             indflat=np.ravel(indices)
             sort=np.argsort(indflat)
         
             scan_sorted=scan_flat[sort]
             move_sorted=move_flat[sort]
             indsorted=indflat[sort]
             filledsorted=filledflat[sort]
             sel=filledsorted>-1
             t_sel=filledsorted[sel]
             ind_sel=indsorted[sel]
             t_all=np.interp(indsorted,ind_sel,t_sel)
             t_array=np.empty((n_move,n_scan))
             for j in range(len(t_all)):
                  t_array[move_sorted[j],scan_sorted[j]]=t_all[j]                 
             # Calculate shift between different timestamps
             
             
             res = minimize_scalar(minimize_shift,args=(n_move,n_scan,timestamps,t_array))
             shift=res.x*0.01
             
             print("Shift: " + str(shift))
             # nur test
             #shift = -0.0145837
             intensity_array=np.empty((n_scan,n_move,len(indices_wave)))
             indices_wave=indices_wave[np.newaxis, np.newaxis, :]
             intensity_array= intensity_array + indices_wave
             
             for j in range(n_move):
                    for i in range(n_scan):
                            m=np.searchsorted(timestamps-shift,t_array[j,i])
                            try:
                                intensity_array[i,j,:]=intensity[m,:]
                            except IndexError:
                                    pass
            
             
              # print('intensityRef_array', intensityRef_array)
               # Create 2D image of measurement
             
            
             
             g.close()
             
             return intensity_array, scan_ax, move_ax
                              
       Int0=Int0(r"C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Juli\2024-06-26\Ref\data_Ref_15-26-33.h5")
       Int1, scan_ax, move_ax=Int1( r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Juli\2024-06-26\E5\data_E5_15-23-33.h5')
       int_ratio = Int0/Int1
       x_index = np.searchsorted(move_ax, x_pos)
       y_index = np.searchsorted(scan_ax, y_pos)
       alpha = np.log10(int_ratio[y_index, x_index, :])/measured_thickness
       np.savetxt(r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\alpha ellipsometrie\0810\1010\E5_log10_2.txt', np.column_stack([wavelengths_array, alpha]), header='wavelengths,alpha', delimiter=' ')
       plt.plot(wavelengths_array, alpha)
       plt.grid(True)
       plt.show()
       g.close()
evaluation = Absorption()
evaluation.intensities()          