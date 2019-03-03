"""
FORS DATA REDUCTION SCRIPT - v1.0 (2019)

N.B. the script structure is similar to a similar script used for MUSE data reduction

To run this script call you need firt to first copy the science files into the main folder
After that, you need to launch the following command:

> python FORS_pipeluca_v1.py

The script will check for the grism(s) for the data into the folder
and then it will ask which grism to reduce first

REQUISITES: astropy, CPL libraries, shutil
"""

from astropy.io import fits
import glob,sys,os,shutil
from subprocess import call

# Ask the user for the directory where the calibration data is located
science_data_directory = os.getcwd()
# this is an input example: /Users/lucaizzo/Documents/ENGRAVE/FORS_data_170817
data_files = glob.glob(str(science_data_directory + '/FORS2.*fits'))

#determine the number of grisms in the dataset
grisms = []
for i in range(len(data_files)):
    temp = fits.open(data_files[i])
    try:
        hierarch = temp[0].header['HIERARCH ESO INS GRIS1 NAME']
    except:
        pass
    if hierarch in grisms:
        pass
    else:
        grisms.append(hierarch)

#create output folders for each grism
for item in grisms:
    try:
        os.makedirs('output_'+item)
    except:
        pass

#show which grisms have been observed
print('')
if len(grisms) == 1:
	print("There is only the "+str(grisms[0])+" dataset")
	grism = grisms[0]
elif len(grisms) == 0:
	print("Something went wrong with the static calibration files ! Check better.")
else:
    print("There are "+str(len(grisms))+" grisms in the dataset.")
    print('')
    print("List of datasets/grisms : ")
    for i in range(len(grisms)):
        print('')
        print(grisms[i])
        print('')
    grism = input('Which grism do you want to reduce ? ')
    if grism in grisms:
        print('')
    else:
        print('This grism is not in the dataset ! Sorry, you need to start from the beinning :(')
        sys.exit[0]


# Prints out each line of a list line by line
# This was used during testing to ensure that the lists contain the correct file names.
def print_lines(input_list):
	for i in range(len(input_list)):
     		print(input_list[i])



# Define variables to be used for input files
# These files must be supplied.
MASTER_LINECAT = 'master_linecat.fits'
EXTINCT_TABLE = 'extinct_table.fits'
STD_FLUX_TABLE = 'std_flux_table.fits'
TELLURIC_CONTAMINATION = 'telluric_contamination.fits'
GRISM_TABLE = 'grism_table.fits'
GLOBAL_DISTORTION_TABLE_1 = 'global.fits'

# Default names for files that can be created.
MASTER_BIAS = 'master_bias.fits'
DISP_COEFF_LSS = 'disp_coeff_lss.fits'
FLAT_SED_LSS = 'flat_sed_lss.fits'
MASTER_NORM_FLAT_LSS = 'master_norm_flat_lss.fits'
SLIT_LOCATION_LSS = 'slit_location_lss.fits'
SPECPHOT_TABLE = 'specphot_table.fits'

# This for loop checks the header for each data file and determines what type it is.
# Now, we define functions that will help in creating sof files
# Creates a text file from a given input list.
def write_list(type_list, file_name):
    name = file_name
    try:
        f = open(name, "w")
        f.write("\n".join(map(lambda x: str(x), type_list)))
        f.close()
    except:
        print('Something went wrong! Can\'t write to the file: ' + name)
        sys.exit(0) # quit Python


# The following are the definitions on how to construct the set of frame (sof) files.
# Each definition also checks to ensure that the required files already exist in the directory.
#this task must be launched for each single observation !!!
def create_science_sof():
    if True:
        for i in range(len(data_files)):
            science_list = []
            science_list.append(data_files[i] + ' ' + 'SCIENCE_LSS')
            temp = fits.open(data_files[i])
            try:
                grism = temp[0].header['HIERARCH ESO INS GRIS1 NAME']
            except:
                print('There is some problem in the input file - check with the WG-IS')
                sys.exit(0)
            science_list.append('./'+grism+'/'+MASTER_BIAS + ' ' + 'MASTER_BIAS')
            science_list.append('./'+grism+'/'+GRISM_TABLE + ' ' + 'GRISM_TABLE')
            science_list.append('./'+grism+'/'+MASTER_NORM_FLAT_LSS + ' ' + 'MASTER_NORM_FLAT_LSS')
            science_list.append('./'+grism+'/'+DISP_COEFF_LSS + ' ' + 'DISP_COEFF_LSS')
            science_list.append('./'+grism+'/'+SLIT_LOCATION_LSS + ' ' + 'SLIT_LOCATION_LSS')
            science_list.append('./'+grism+'/'+FLAT_SED_LSS + ' ' + 'FLAT_SED_LSS')
            science_list.append('./'+grism+'/'+TELLURIC_CONTAMINATION + ' ' + 'TELLURIC_CONTAMINATION')
            science_list.append('./'+grism+'/'+EXTINCT_TABLE + ' ' + 'EXTINCT_TABLE')
            science_list.append('./'+grism+'/'+SPECPHOT_TABLE + ' ' + 'SPECPHOT_TABLE')
            write_list(science_list, 'science'+str(i)+'.sof')
            print('')
            print('Created science'+str(i)+'.sof')
        #sys.exit(0)



# Here start the data reduction !!!

if True:
    try:
        create_science_sof()
        print('---Running FORS_science for target(s)---', shell = True)
        for i in range(len(data_files)):
            call('esorex --log-file=science.log fors_science --skyalign=-1 --skyglobal=TRUE --skylocal=FALSE --cosmics=TRUE science'+str(i)+'.sof', shell = True)
            try:
                grism = temp[0].header['HIERARCH ESO INS GRIS1 NAME']
            except:
                break
            os.rename('mapped_all_sci_lss.fits', 'output_'+grism+'/mapped_all_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('unmapped_sky_sci_lss.fits', 'output_'+grism+'/unmapped_sky_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('unmapped_sci_lss.fits', 'output_'+grism+'/unmapped_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('mapped_sky_sci_lss.fits', 'output_'+grism+'/mapped_sky_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('object_table_sci_lss.fits', 'output_'+grism+'/object_table_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('reduced_flux_sci_lss.fits', 'output_'+grism+'/reduced_flux_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('reduced_flux_error_sci_lss.fits', 'output_'+grism+'/reduced_flux_error_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('reduced_sci_lss.fits', 'output_'+grism+'/reduced_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('reduced_sky_sci_lss.fits', 'output_'+grism+'/reduced_sky_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('reduced_error_sci_lss.fits', 'output_'+grism+'/reduced_error_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('mapped_flux_sci_lss.fits', 'output_'+grism+'/mapped_flux_sci_'+str(i)+'_'+str(grism)+'.fits')
            os.rename('mapped_sci_lss.fits', 'output_'+grism+'/mapped_sci_'+str(i)+'_'+str(grism)+'.fits')
            end_call = True
    except:
        print('ERROR: Could not run FORS_science for target(s)')
        sys.exit(0)


#removing files
for i in range(len(data_files)):
    os.remove('science'+str(i)+'.sof')

try:
    os.remove('disp_coeff_sci_lss.fits')
    os.remove('sky_shifts_slit_sci_lss.fits')
    os.remove('global_sky_spectrum_lss.fits')
    os.remove('wavelength_map_sci_lss.fits')
    os.remove('science.log')
except:
    pass

print('You have successfully reduced the '+str(len(data_files))+' spectra of you target !!!')
print('')
print('2D flux-calibrated spectra are saved in the corresponding output_GRISM folder')
print('')
print('Flux-calibrated 2D spectra are stored in the files starting with "mapped_flux_sci"')
print('')
print('Good work and enjoy !!!')
