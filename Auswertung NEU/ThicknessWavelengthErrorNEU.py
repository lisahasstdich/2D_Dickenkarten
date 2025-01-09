import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.optimize import minimize_scalar
import pandas as pd
import variables, h5py
import variablesRef, h5py
import locale
# Setze die Lokalisierung auf Deutsch
locale.setlocale(locale.LC_NUMERIC, 'de_DE.UTF-8')
unten = 0
oben = 450 #zuschneideparameter aus Dickekarte klauen
#file = r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Oktober\Absorptionskoeffizienten\Alpha_E5_dick.txt' #file für die berechneten Absorptionskoeffizienten
file = r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Oktober\Absorptionskoeffizienten\Mean_alpha_CFUNDCFCB.txt'
plt_show = False # einzelne Plots der dicke für die verschiedenen Wellenlängen, wenn man zb nur ein paar einzelne Wellenlängen anschauen will
calculate_thickness_difference = True #hier ist der Plot für die Differenz mitgemeint
save_difference_plot = True
plot_spektrum = False #Spektrum an drei Stellen des Films
plt.rcParams.update({
    'font.size': 13,        # Allgemeine Schriftgröße
    'axes.titlesize': 13,   # Titel der Achsen
    'axes.labelsize': 13,   # Achsenbeschriftung
    'xtick.labelsize': 13,  # x-Achsen-Ticks
    'ytick.labelsize': 13,  # y-Achsen-Ticks
    'legend.fontsize': 13,  # Legende
    'figure.titlesize': 16.5,  # Titel der gesamten Abbildung
    'axes.formatter.use_locale' : True
})
class Absorption:
    def Initializing(self):
        self.intensity_array = None
        self.intensityRef_array = None
        
    def intensities(self,h5file = None, interpolation = None):
       h5file        = h5file        or f'{variablesRef.path}'
       interpolation = interpolation or variables.interpolation
       
       g = h5py.File(h5file, "r")
       # wavelength_min = np.min(g['spectrum/wavelength']) #Wellenlängen, die man ausgewertet haben möchte
       # wavelength_max = np.max(g['spectrum/wavelength']) #mit dieser Einstellung hat man ALLE Wellenlängen, kann aber wg teilweise positivem Alpha zu RuntimeWarnings kommen
       wavelength_min = 570
       wavelength_max = 620
       diff_wavelength = wavelength_max-wavelength_min
       wavelengths_array = g['spectrum/wavelength'][:]
       index_min_1= np.searchsorted(wavelengths_array , wavelength_min)
       index_max_1= np.searchsorted(wavelengths_array, wavelength_max)
       df = pd.read_csv(file, delimiter=' ')  # delimiter für Leerzeichen
       alphas = df.iloc[:, 1]
      # alpha_errors = df.iloc[:,2]
       alpha_array = np.array(alphas)
      # alpha_error =np.array(alpha_errors) # importieren aus excel datei oder sonstwo!!! sollte das gemittelte aus verschiedenen berechneten von alphavaluesarraylike sein!
       alpha_array=alpha_array[index_min_1:index_max_1]
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
             index_min= np.searchsorted(wavelength_values, wavelength_min)
             index_max= np.searchsorted(wavelength_values, wavelength_max)
             wavelength_array=wavelength_array[index_min:index_max]
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
             
             if plt_show == True:
                 for i in range(0,len(indices_wave)):
                   fig, ax = plt.subplots()
                   
                   
                   ax.set_xlabel('x (mm)')
                   ax.set_ylabel('y (mm)')
                   ax.set_title(f'Intensity für Wellenlänge{i+wavelength_min}')
                   # intensity_arry could have to be transpose to get the right orientation of the picture
                   if variables.scan_axis == "xy":
                       intensityRef_array[:,:,i] = np.rot90(intensityRef_array)
                     
                       im=ax.imshow(intensityRef_array,aspect='equal',
                                interpolation=interpolation,
                                extent=[scan_min,scan_max,move_min,move_max],
                                origin='upper')
             
                   if variables.scan_axis == "yx":
                      im=ax.imshow(intensityRef_array[:,:,i],aspect='equal',
                                interpolation=interpolation,
                                extent=[scan_min,scan_max,move_min,move_max],
                                origin='lower')
                   cbar = plt.colorbar(im)
                   plt.show()
                   plt.close()
             
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
             index_min= np.searchsorted(wavelength_values, wavelength_min)
             index_max= np.searchsorted(wavelength_values, wavelength_max)
             wavelength_array=wavelength_array[index_min:index_max]
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
             
             if plt_show == True:
                 for i in range(0,np.max(indices_wave)):
                   fig, ax = plt.subplots()
                   
                   
                   ax.set_xlabel('x (mm)')
                   ax.set_ylabel('y (mm)')
                   ax.set_title(f'Intensity für Wellenlänge{i+wavelength_min}')
                   # intensity_arry could have to be transpose to get the right orientation of the picture
                   if variables.scan_axis == "xy":
                       intensity_array[:,:,i] = np.rot90(intensity_array)
                     
                       im=ax.imshow(intensity_array,aspect='equal',
                                interpolation=interpolation,
                                extent=[scan_min,scan_max,move_min,move_max],
                                origin='upper')
             
                   if variables.scan_axis == "yx":
                      im=ax.imshow(intensity_array[:,:,i],aspect='equal',
                                interpolation=interpolation,
                                extent=[scan_min,scan_max,move_min,move_max],
                                origin='lower')
                   cbar = plt.colorbar(im)
                   plt.show()
                   plt.close()
             
             g.close()
             
             return intensity_array, scan_ax, move_ax
                              
       Int0=Int0(r"C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Juli\2024-06-26\Ref\data_Ref_15-26-33.h5")
       Int1, scan_ax, move_ax=Int1( r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Juli\2024-06-26\E7\data_E7_15-22-03.h5')
       int_ratio = Int0/Int1
       thickness = (np.log10(int_ratio)/alpha_array)
       thickness = thickness[unten:oben, :,:]
       thickness = np.clip(thickness, 0, None)
       
       scan_ax = scan_ax[unten:oben] #Dicke ist schon so zugeschnitten!!
       
       if plt_show == True : #hab keine ahnung was hier abgeht
           for i in range(diff_wavelength):
               fig, ax = plt.subplots()
               
               ax.set_title(f'Dicke(nm) für Wellenlänge{i+wavelength_min}')
               ax.set_xlabel('x (mm)')
               ax.set_ylabel('y (mm)')
               plt.imshow(thickness[:,:,i], cmap=plt.cm.get_cmap('Blues', 15), aspect='equal',
                          interpolation=interpolation,
                          extent=[np.min(move_ax), np.max(move_ax), np.min(scan_ax), np.max(scan_ax)],
                          origin='lower')
               plt.colorbar(extend='max')
               plt.clim(0, 5)
               plt.show()
               plt.close()
       if calculate_thickness_difference == True:
           
           max_thickness = thickness.max(axis=2)
           min_thickness = thickness.min(axis=2)
           diff = max_thickness - min_thickness
           fig, ax = plt.subplots(figsize=(4,3))
           ax.set_title(f'Maximale Dickendifferenz (nm) im \n Wellenlängenbereich {wavelength_min}-{wavelength_max}nm')
           ax.set_xlabel('x (mm)')
           ax.set_ylabel('y (mm)')
           plt.imshow( diff, cmap='viridis', aspect='equal',
                       interpolation=interpolation,
                       extent=[np.min(move_ax), np.max(move_ax), np.min(scan_ax), np.max(scan_ax)],
                       origin='lower')
           plt.colorbar(extend='max')
           plt.clim(0,20)
           if save_difference_plot == True:
               fig.savefig(r'C:\Users\Lisa\Desktop\DickenkartenSVG\Fehlerkarten\E5_fehler_570620_dünn.svg', format='svg')
           plt.show()
           plt.close()
       if plot_spektrum == True:
             
             x_1 = 5
             x_2 = 11
             x_3 = 17
             x_1_index = np. searchsorted(move_ax, x_1)
             x_2_index = np. searchsorted(move_ax, x_2)
             x_3_index = np. searchsorted(move_ax, x_3)
             y_mitte = 11.25
             y_mitte_index = np.searchsorted(scan_ax, y_mitte)
             alpha_spek_1 = np.log10(Int0[y_mitte_index, x_1,:]/Int1[y_mitte_index, x_1_index, :])
             alpha_spek_2 = np.log10(Int0[y_mitte_index, x_2,:]/Int1[y_mitte_index, x_2_index, :]) 
             alpha_spek_3 = np.log10(Int0[y_mitte_index, x_3,:]/Int1[y_mitte_index, x_3_index, :])
             plt.plot(wavelengths_array[index_min_1:index_max_1], alpha_spek_1, color='red', linewidth =2, linestyle='-', label = 'x=5mm')
             plt.plot(wavelengths_array[index_min_1:index_max_1], alpha_spek_2, color='blue', linewidth =2, linestyle='-', label = 'x=11mm')
             plt.plot(wavelengths_array[index_min_1:index_max_1], alpha_spek_3, color='orange', linewidth =2, linestyle='-',label = 'x=17mm')
             plt.xlabel('wavelength (nm)')
             plt.ylabel('absorption')
             plt.title('Spectrum for E7 at y=11.25mm')
             plt.legend()
             plt.show()
             plt.close() 
       g.close()
evaluation = Absorption()
evaluation.intensities()          
