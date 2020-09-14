"""
    PROJECT NAME: CONTRAST CURVE
    
    DESCRIPTION: Turns contrast/radius curves to mass/radius

    INSTRUCTION: 
    
    
    VERSION	3.0
    DATE 11/09/2020

    Note: 

    1. Beware of the age of the model


    To-do: 
"""

import glob
import os
import numpy as np
import math
import csv
from scipy import interpolate
from scipy.interpolate import griddata

"""
    DEFINED VARIABLES
"""
INPUT_NAMES = 'star_names.txt'           #file in which to find names of stars to analysis and age, distance, mag data
INPUT_BARAFFE_DATA = 'baraffe_new.txt'  #data from baraffe
seperation_column = 3                      #column of input data containg seperation in arcsec
contrast_column = 0                       #cloumn of input data containing contrast
ONE_PARSEC = 206268.0                    #AU scale to use on axis of output
JUP_MASS = 1047.34
AGE_STAR_UNIT = 1.0e6
AGE_MODEL_UNIT = 1.0e9
PLOT_PRECISION = 100               #defines grid side for interpolation
SAVEFIG = True                    # Change to 'False' to show plots of mass/radius while program is running (does not save fig)
SHOW_SURFACE = False #Change to 'True' if plot of surface of baraffe data is required. requires mpl_toolkits to be available

"""
    FUNCTION DEFINITIONS
"""

#returns file conaining the contrast data
def get_contrast_data_filename(name):
    os.chdir(os.getcwd())
    for file in glob.glob("*%s*" %(name)):
        if os.path.isfile(file):
            print(file)
            return file
    else:
        print("\nError: Could not find file")
        exit()

#loads contrast data from file into array and returns orbital seperation and mag
def load_contrast_data(filename):
    data = np.loadtxt(filename,skiprows=1)

    contrast_curve = np.zeros( (len(data[:,0]), 2) )    
    contrast_curve[:,0] = data[:,seperation_column]
    contrast_curve[:,1] = data[:,contrast_column]
    
    return contrast_curve


#returns absolute magnitude
def absolute_mag(distance, mag):
    M = mag - ((5 * math.log10(distance)) - 5.0)
    return M

#returns the magnitude of the companion as array
def companion_mag_curve(contrast_data, abs_mag):
    return (np.absolute(2.5 * np.log10(contrast_data)) + abs_mag)

#turns angular seperation into physical speration and returns
def orbital_sep(angular_sep, distance):
    rad_coef = math.pi /180.0 #convert from degrees to radians
    arc_sec_coef = 1 / 3600.0 #convert from arcseconds to degrees
    VIP_scaling_coef = 0.00953 #factor of 100 in demoninator due to unknown scaling of data files from VIP
    angular_sep = angular_sep * rad_coef * arc_sec_coef * VIP_scaling_coef
    return distance * np.tan(angular_sep) * ONE_PARSEC

#plots or saves plot of mass againts orbital seperation sensitivity curve
def plot_rad_mass(rad_mass,file_dest, name):
    import matplotlib.pyplot as plt
    plt.clf()

    from matplotlib.ticker import MaxNLocator
    ax = plt.figure().gca()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    best, = plt.plot(rad_mass[:,0], rad_mass[:,1],'r', label='Best age estimate')
    oldest, = plt.plot(rad_mass[:,0],rad_mass[:,2],'b', label='Oldest age estimate')
    youngest, = plt.plot(rad_mass[:,0],rad_mass[:,3],'g', label='Youngest age estimate')

    plt.fill_between(rad_mass[:,0], rad_mass[:,3], rad_mass[:,2],color='grey')
    plt.yscale('linear')
    plt.legend(loc='upper right')
    plt.ylabel(r'Mass ($\mathrm{M_{Jup}}$)')
    plt.xlabel('Orbital separation (AU)')
    plt.legend(frameon=False)
    #plt.title('Mass - Orbital Separation Sensitivity Plot (%s)' %(name))
    
    if SAVEFIG:
        plt.savefig(file_dest,dpi=300,transparent=True)
    else:
        plt.show()


#finds nan's in interpolated grid and changes them to nearest neighbur value (flatterns area) but stops interpolation errors
def check_nan(grid):
    lum=0.0
    for i in range(0,PLOT_PRECISION):
        for j in range(0,PLOT_PRECISION):
            if np.isnan(grid[i,j]):
                grid[i,j]=lum
            else:
                lum = grid[i,j]
    return grid

#does chi suqered test of data to determine how good a fit function is to baraffe data
def chi_squared(points, lum, func):
    chi = 0.0
    for i in range(0,len(points)):
        chi += ((lum[i]-func(points[i,0], points[i,1]))**2)/func(points[i,0], points[i,1])
    return chi

#saves figures and output data to seperate files (creatng file if not present)
def save_data(name_data, rad_mass):
    current_dir = os.getcwd()
    final_dir = os.path.join(current_dir, r'%s_data' %(str(int(name_data[i,0]))))
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
    
    data_file = os.path.join(final_dir, r'rad_mass_data_%s.txt' %(str(int(name_data[i,0]))))
    fig_file = os.path.join(final_dir, r'rad_mass_%s.png' %(str(int(name_data[i,0]))))

    plot_rad_mass(rad_mass, fig_file, str(int(name_data[i,0])))
    np.savetxt(data_file, rad_mass,header="Orbital seperation (AU) Mass_Limit(Best) Mass_Limit(Oldest) Mass_Limit(Youngest)")

def plot_surface_baraffe(ai_grid, mi_grid, lumi):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig=plt.figure()
    ax=fig.gca(projection='3d')
    ax.plot_surface(ai_grid, mi_grid,lumi)
    plt.show()

"""
    CODE IMPLENTATION
"""
#load baraffe data and star data
name_data = np.loadtxt(INPUT_NAMES)             
print (name_data)
baraffe_data = np.loadtxt(INPUT_BARAFFE_DATA)


#seperates data points from baraffe into points(age,mass) and magnitude list, lum
points = np.zeros((len(baraffe_data),2))        
points[:,0] = baraffe_data[:,0]
points[:,1] = baraffe_data[:,1]
lum = np.zeros((len(baraffe_data),1))
lum = baraffe_data[:,2]

#Sets up mesh of values to interpolate surface from the scattered data points from baraffe
ai = np.linspace(np.amin(points[:,0]),np.amax(points[:,0]),num = PLOT_PRECISION, endpoint=True) 

#the precision of the grid i set by PLOT_PRECISION (ai=age, mi=mass)
mi = np.linspace(np.amin(points[:,1]),np.amax(points[:,1]),num = PLOT_PRECISION, endpoint=True) 
ai_grid, mi_grid = np.meshgrid(ai, mi)


#interpolate from baraffe data points to surface grid lumi
lumi = griddata(points, lum, (ai_grid, mi_grid), method='cubic')
#remove any nan values form grid  
lumi = check_nan(lumi)
#create function from grid                                          
func_lum = interpolate.interp2d(ai, mi, lumi,kind='quintic')  

#print results of chi squared test    
print('Chi Squared test of function verse known data points =')   
print(chi_squared(points, lum, func_lum))

if SHOW_SURFACE:
    plot_surface_baraffe(ai_grid, mi_grid, lumi)

#loop over all stars in file outputing data for each star
for i in range(0,len(name_data)):                   
    #load correct contrast data star
    contrast_curve_file = get_contrast_data_filename(str(int(name_data[i,0]))) 
    contrast_curve = load_contrast_data(contrast_curve_file)
    
    #creates array of orbital seperation and companion magnitude
    mag_curve  = np.zeros((len(contrast_curve),2))      
    mag_curve[:,0] = orbital_sep(contrast_curve[:,0], name_data[i,4])
    mag_curve[:,1] = companion_mag_curve(contrast_curve[:,1], absolute_mag(name_data[i,4], name_data[i,5]))
    
    star_age = np.zeros((3))
    #creates new array to contain mass and orbital seperation
    rad_mass = np.zeros((len(mag_curve),4)) 
    rad_mass[:,0] = mag_curve[:,0]
    
    for j in range(0,3):
        # Convert the age of the star to the same unit as the age of the model. 
        star_age[j] = np.multiply(AGE_STAR_UNIT,name_data[i,j+1])/AGE_MODEL_UNIT 
        massi = np.linspace(np.amin(points[:,1]),np.amax(points[:,1]),num=PLOT_PRECISION, endpoint=True) #creates new linespace of mass to create function of mass in terms of magnitude
        
        #takes slice of baraffe data surface defined by star age
        lumi = func_lum(star_age[j], massi)        
        check = False
        
        if lumi[0] > lumi[-1]:
            lumi = lumi[::-1]
            massi = massi[::-1]
            check=False

        #interpolates function mass(magnitude) from baraffe data slice and new mass linespace
        func_mass = interpolate.interp1d(np.ravel(lumi), massi, kind='linear', fill_value='extrapolate') 
        mag = mag_curve[:,1]
   
        if np.amin(mag) < np.amin(lumi):
            print('Warning data below Baraffe model bounds....Extrapolating')
                    
        if np.amax(mag) > np.amax(mag):
            print('Warning data above Baraffe model bounds....Extrapolating')

        #uses function mass(magnitude) to output mass for a give magnitude for each point of the load ontrast curves from VIP
        rad_mass[:,j+1]=np.multiply(JUP_MASS, func_mass(mag))
    
    save_data(name_data, rad_mass)  #save relevant data

