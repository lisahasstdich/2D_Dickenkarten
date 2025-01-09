import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import pandas as pd
from scipy.optimize import minimize_scalar
import variables, h5py
import variablesRef, h5py
import locale
# Setze die Lokalisierung auf Deutsch
locale.setlocale(locale.LC_NUMERIC, 'de_DE.UTF-8')
sample = 'M34'
used_alpha ='dünn'
savepath = r'C:\Users\Lisa\Desktop\DickenkartenSVG\Profile\\'
file = r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\Oktober\Absorptionskoeffizienten\Mean_alpha_CFUNDCFCB.txt'
#file = r'C:/Users/Lisa/Desktop/Uni/Bachelorarbeit/Oktober/Absorptionskoeffizienten/Alpha_E5_dick.txt'
#file = r'C:/Users/Lisa/Desktop/Uni/Bachelorarbeit/Oktober/Absorptionskoeffizienten/Alpha_E5_mittel.txt'
unten = 10
oben = 395 #Indices bei denen abgeschnitten werden soll
links = 5
rechts = 70
wavelength_evaluate = 590
plot_show_default = False
plot_show_zugeschnitten = True
plot_einzelprofil_y = True
plot_profiltiefe_y = False
plot_profiltiefe_x = True
plt.rcParams.update({
    'font.size': 11.5,        # Allgemeine Schriftgröße
    'axes.titlesize': 11.5,   # Titel der Achsen
    'axes.labelsize': 11.5,   # Achsenbeschriftung
    'xtick.labelsize': 11.5,  # x-Achsen-Ticks
    'ytick.labelsize': 11.5,  # y-Achsen-Ticks
    'legend.fontsize': 10.5,  # Legende
    'figure.titlesize': 17,  # Titel der gesamten Abbildung
    'axes.formatter.use_locale' : True
})
class Absorption:
    def Initializing(self):
        self.intensity_array = None
        self.intensityRef_array = None
        
    def intensities(self,h5file = None,interpolation = None):
        h5file        = h5file        or f"{variablesRef.path}"
        interpolation = interpolation or variables.interpolation
        
        g = h5py.File(h5file, "r")
        wavelength_min = np.min(g[f"{variablesRef.data_group[1]}/wavelength"]) #Wellenlängen, die man ausgewertet haben möchte
        wavelength_max = np.max(g[f"{variablesRef.data_group[1]}/wavelength"])
        diff_wavelength = wavelength_max-wavelength_min
        wavelengths_array = g[f"{variablesRef.data_group[1]}/wavelength"][:]
        index_wavelength = np.searchsorted(wavelengths_array, wavelength_evaluate)
        
        g.close()
        
        def Int0(self, h5file = None,
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
                          
                          diff=diff+abs(intensityRef[m]-intensityRef[n])
                  
                  return diff
              # Initialize varibales from variablesRef.py
              scan_grid     = scan_grid     or variablesRef.scan_grid
              move_grid     = move_grid     or variablesRef.move_grid
              interpolation = interpolation or variablesRef.interpolation
              h5file        = h5file        or f"{variablesRef.path}"
              
              g = h5py.File(h5file, "r")
              # Initialize data for intensityRef
              timestamps = g[f"{variablesRef.data_group[1]}/timestamp"][:]
              intensityRef = g[f'{variablesRef.data_group[1]}/intensity'][:,index_wavelength]
              
              
              x, y, z, t = g[f"{variablesRef.data_group[0]}/X"][:], g[f"{variablesRef.data_group[0]}/Y"][:], g[f"{variablesRef.data_group[0]}/Z"][:], g[f"{variablesRef.data_group[0]}/meas_time"][:]
              xmin, xmax, ymin, ymax, zmin, zmax = x.min(), x.max(), y.min(), y.max(), z.min(), z.max()                           
              
              
               # Switch evaluation according the move direction
              if variablesRef.scan_axis == "yx":
                   scan_dir, scan_min, scan_max, scan_label = y, ymin, ymax, "y"   
                   move_dir, move_min, move_max, move_label = x, xmin, xmax, "x"
          
              if variablesRef.scan_axis == "xy":
                   scan_dir, scan_min, scan_max, scan_label = x, xmin, xmax, "x"
                   move_dir, move_min, move_max, move_label = y, ymin, ymax, "y"
              
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
              intensityRef_array=np.empty((n_scan,n_move))
              
              
              for j in range(n_move):
                     for i in range(n_scan):
                             m=np.searchsorted(timestamps-shift,t_array[j,i])
                             try:
                                 intensityRef_array[i,j]=intensityRef[m]
                             except IndexError:
                                     pass
              
              
              
                # Create 2D image of measurement
              g.close()
              return intensityRef_array
        def Int1(self, h5file = None,
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
                          
                          diff=diff+abs(intensity[m]-intensity[n])
                  
                  return diff
              # Initialize varibales from variablesRef.py
              scan_grid     = scan_grid     or variables.scan_grid
              move_grid     = move_grid     or variables.move_grid
              interpolation = interpolation or variables.interpolation
              h5file        = h5file        or f"{variables.path}"
              
              g = h5py.File(h5file, "r")
              # Initialize data for intensityRef
              timestamps = g[f"{variables.data_group[1]}/timestamp"][:]
              intensity = g[f'{variables.data_group[1]}/intensity'][:,index_wavelength]
              
              
              x, y, z, t = g[f"{variables.data_group[0]}/X"][:], g[f"{variables.data_group[0]}/Y"][:], g[f"{variables.data_group[0]}/Z"][:], g[f"{variables.data_group[0]}/meas_time"][:]
              xmin, xmax, ymin, ymax, zmin, zmax = x.min(), x.max(), y.min(), y.max(), z.min(), z.max()                           
              
              
               # Switch evaluation according the move direction
              if variables.scan_axis == "yx":
                   scan_dir, scan_min, scan_max, scan_label = y, ymin, ymax, "y"   
                   move_dir, move_min, move_max, move_label = x, xmin, xmax, "x"
          
              if variables.scan_axis == "xy":
                   scan_dir, scan_min, scan_max, scan_label = x, xmin, xmax, "x"
                   move_dir, move_min, move_max, move_label = y, ymin, ymax, "y"
              
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
              intensity_array=np.empty((n_scan,n_move))
              for j in range(n_move):
                     for i in range(n_scan):
                             m=np.searchsorted(timestamps-shift,t_array[j,i])
                             try:
                                 intensity_array[i,j]=intensity[m]
                             except IndexError:
                                     pass
             
              g.close()
              
              return intensity_array, scan_ax, move_ax
        Int0=Int0( r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\mai\Filme\2024-05-23\15-59-45\data_15-59-45.h5')
        Int1, scan_ax, move_ax=Int1(r'C:\Users\Lisa\Desktop\Uni\Bachelorarbeit\mai\Filme\2024-05-23\15-56-56\data_15-56-56.h5')
        int_ratio = Int0/Int1 
        df = pd.read_csv(file, delimiter=' ')  # delimiter für Leerzeichen
        alphas = df.iloc[:, 1]
        #alpha_errors = df.iloc[:,2]
        alpha = np.array(alphas)
       # alpha_error =np.array(alpha_errors)
        alpha = alpha[np.newaxis, np.newaxis, index_wavelength]
        print('alpha', alpha)
        dicke = np.log10(int_ratio) / alpha
        
        dicke = np.clip(dicke, 0, None)
        if plot_show_default == True:
            fig, ax = plt.subplots(figsize=(3,3))
            
            ax.set_title('Dicke (nm)')
            ax.set_xlabel('x (mm)')
            ax.set_ylabel('y (mm)')
            plt.imshow(dicke, cmap=plt.cm.get_cmap('Blues', 15), aspect='equal',
                       interpolation=interpolation,extent=[np.min(move_ax), np.max(move_ax), np.min(scan_ax), np.max(scan_ax)],
                       origin='lower')
            plt.colorbar(extend='max')
            plt.clim(0,100)
            fig.savefig(savepath + f'{sample}_dicke_{used_alpha}.svg', format='svg')
            plt.show()
            plt.close()
            
        if plot_show_zugeschnitten ==True:
            move_ax=move_ax[links:rechts]
            scan_ax=scan_ax[unten:oben]
            dicke = dicke[unten:oben,links:rechts]
            fig, ax = plt.subplots()
            
            ax.set_title('Dicke (nm)')
            ax.set_xlabel('x (mm)')
            ax.set_ylabel('y (mm)')
            plt.imshow(dicke, cmap=plt.cm.get_cmap('Blues', 15), aspect='equal',
                       interpolation=interpolation,extent=[np.min(move_ax), np.max(move_ax), np.min(scan_ax), np.max(scan_ax)],
                       origin='lower')
            plt.colorbar(extend='max')
            plt.clim(0, 70)
           # fig.savefig(savepath + f'{sample}_dicke_passend_zugeschnitten.png', bbox_inches='tight', dpi=300)
            plt.show() 
            plt.close()
            if plot_einzelprofil_y == True:
                y_einzel = [2, 6,10,14,18]
                y_index = np.searchsorted(scan_ax, y_einzel)
                # Erstellen der Abbildung und der Achsen
                fig, ax = plt.subplots(figsize=(4.5,2.5))
                
                # Für jeden Y-Wert das Profil in X-Richtung plotten
                for y in y_index:
                    ax.plot(move_ax, dicke[y,:], linestyle='-', label=f'{scan_ax[y]:.1f}'.replace('.', ','))
                # Plot-Anpassungen
                ax.set_title('Dickenprofile in x-Richtung')
                ax.set_xlabel('x (mm)')
                ax.set_ylabel('Dicke (nm)')
                ax.legend(title='y (mm)',loc='upper right', ncol=1)
                ax.grid(True)
                
                # ax.xaxis.set_major_locator(MultipleLocator(2.5))
                ax.set_xticks(np.arange(2.5, 17.5, 2.5))
                # Anzeigen des Plots
                fig.savefig(savepath + f'{sample}_y_einzelprofil_{used_alpha}.svg', format='svg')
                plt.show()
                plt.close()
        if plot_profiltiefe_y == True:
            move_ax=move_ax[links:rechts]
            scan_ax=scan_ax[unten:oben]
            dicke = dicke[unten:oben,links:rechts]
            dickemax = np.max(dicke, axis=None)
            dickemin = np.min(dicke, axis=None)
            
            desired_y_values = [2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14,15,16,17,18,19]           
            desired_y_index = np.searchsorted(scan_ax, desired_y_values)
            dicke_neu_y = dicke[desired_y_index,:]
            mean_dicke_neu_y = np.mean(dicke_neu_y, axis=0)
            mean_dicke_diff_y = np.max(mean_dicke_neu_y) - np.min(mean_dicke_neu_y)
            plt.plot(move_ax, mean_dicke_neu_y, color="red", linewidth=2., linestyle="-")
            plt.xlabel('x (mm)')
            plt.ylabel('Mittlere Dicke in x-Richtung (nm)')
            plt.xlim(np.min(move_ax), np.max(move_ax))
           # plt.savefig(savepath + f'{sample}_Profil_y_gemittelt.png', bbox_inches='tight', dpi=300)
            np.savetxt(savepath+f'{sample}_Profil_y_gemittelt_{used_alpha}.txt', np.column_stack([move_ax,mean_dicke_neu_y]), header='x, mean_dicke_y', footer=f'gemittelt über y-Werte: {desired_y_values}', delimiter=' ')
            plt.show()           
            plt.close()
        if plot_profiltiefe_x == True:
            move_ax=move_ax[links:rechts]
            scan_ax=scan_ax[unten:oben]
            dicke = dicke[unten:oben,links:rechts]
            dickemax = np.max(dicke, axis=None)
            dickemin = np.min(dicke, axis=None)
            desired_x_values = [3,4,5,6,7,8,9, 10, 11, 12, 13, 14,15,16]           
            desired_x_index = np.searchsorted(move_ax, desired_x_values)
            dicke_neu_x = dicke[:,desired_x_index]
            mean_dicke_neu_x = np.mean(dicke_neu_x, axis=1)
            mean_dicke_diff_x = np.max(mean_dicke_neu_x) - np.min(mean_dicke_neu_x)
            plt.figure(figsize=(3,2.5))
            plt.plot(scan_ax, mean_dicke_neu_x, color="red", linewidth=2., linestyle="-")
            plt.xlabel('y (mm)')
            plt.title('Dickenprofil in y-Richtung')
            plt.ylabel('Mittlere Dicke (nm)')
            plt.xlim(np.min(scan_ax), np.max(scan_ax))
            plt.grid(True)
            plt.savefig(savepath + f'{sample}_{used_alpha}_Profil_x_gemittelt.svg', format='svg')
            np.savetxt(savepath+f'{sample}_{used_alpha}_Profil_x_gemittelt.txt', np.column_stack([scan_ax,mean_dicke_neu_x]), header='y, mean_dicke_x', footer=f'gemittelt über x-Werte: {desired_x_values}', delimiter=' ')
            plt.show()           
            plt.close()
evaluation = Absorption()
evaluation.intensities()          