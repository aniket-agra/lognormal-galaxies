#!/usr/bin/env python3
#============================================================================
# This is the Python script of calculating the power spectrum from the
# mock galaxy catalogue generated in /genmock folder.
#
# First parts calculate the number of galaxyes N_g, and the number of 
# random samples N_s inside of a grid.
# [1] calc_Fourier_box::
#		calculate the box size, and 
# [2] galaxy_to_grid::
#		put galaxies into the Fourier grid
# [3] calc_Luminosity_fn::
#		code to estimate n(L) from given galaxy sample.
#
# the resulting files are as follow:
# 1) size of the Fourier box
# ../Grids/fourier_dimension_[selection name]_[survey name].dat 
# 2) galaxy numbers in the grid
# ../Grids/gpicc_[survey name]_[selection name]_nmesn=[nmesh].bin 
# 3) galaxy luminosity function n(L)
# ../Grids/gpicc_[survey name]_nmesn=[nmesh].bin 
# 
# The second part generates the random sample inside of the box, 
# applies exactly the same selection as LAEs, and assign the numbers to
# the Fourier grid.
# [1] rsample_to_grid::
#		Select random samples for given selection scheme
#		selection=0	::	no selection
#		selection=1	::	geometry
#		selection=2	::	shot only
#		selection=3	::	shot+mask
#		selection=4	::	radial only
#		selection=5	::	shot+mask+radial
#
# the resulting files are:
# 1) number of random smaples in the grid
# ../Grids/random_[survey name]_[selection name]_rspace.bin
#
# In the third part, we estimate the power spectrum of LAEs from the 
# number of LAEs/random samples in grids.
#
# 24 July 2011
# by Donghui Jeong
#============================================================================

import os
import random
class executable:
    """A class for running executables"""
    def __init__(self,exename):
        self.exename=exename
        return
    
    def __call__(self,params):
        #print(params,len(params))
        nparams = len(params)
        paramfname = params[nparams-1]
        if os.path.exists(paramfname):
            os.remove(paramfname)
        pfile=open(paramfname,'x')
        for iparam in range(0,nparams-1):
            try:
                pfile.write(params[iparam]+'\n')
            except:
                if type(params[iparam])==int:
                    pfile.write('%d\n' % params[iparam])
                else:
                    pfile.write('%g\n' % params[iparam])
        pfile.close()
        cmd = 'time ./'+self.exename+' <'+paramfname
        print(cmd)
        os.system(cmd)
        return

#-----------------------------------------------------------------------------
# function to read number of lines in the file
#-----------------------------------------------------------------------------
def file_len(fname):
	f=open(fname)
	for i, l in enumerate(f):
		pass
	f.close()
	return i + 1

#-----------------------------------------------------------------------------
# how many realisation? initial seed?
#-----------------------------------------------------------------------------
##nrealisation = int(raw_input('N_realisation: '))
##iseed = int(raw_input('iseed: '))
#tmp = raw_input('N_realisation & iseed [separated by space]: ')
#nrealisation = int(tmp.split()[0])
#iseed = int(tmp.split()[1])

#print 'the number of realisation = ',nrealisation
#print 'the initial seed = ',iseed

nrealisation = 1
iseed = 100
print('the number of realisation = ',nrealisation)
random.seed(iseed)
print('the initial seed = ',iseed)

#-----------------------------------------------------------------------------
# set default values 
#-----------------------------------------------------------------------------
# cosmological parameters
Omegam=0.272		# Omega matter for computing the velocity field
Omegade=1-Omegam
zz=2.2			    # redshift for computing the velocity field

def fgrowth(zred):
    Omz = Omegam*(1+zred)**3/(Omegam*(1+zred)**3 + 1-Omegam)
    return Omz**0.545

fgrowth = fgrowth(zz)		# f=dlnD/dlna
aHz=100.*pow(Omegam*pow(1.+zz,3)+Omegade,0.5)/(1.+zz)

# parameters for lognormal realization
Ngalaxies = 8345000	# number of galaxy in integer
bias = 1.5	 # linear bias for computing the matter fluctuation and then velocity field (need to be consistent with input P(k))
mbias = 1	 # enter value used in log-transformed matter power spectrum file name					#modified by Aniket

Lx = 2880	 # box size in x [h^-1 Mpc]
Ly =  480	 # box size in y [h^-1 Mpc]
Lz =  632	 # box size in z [h^-1 Mpc]
#Pnmax = 2500	# mesh number of the longest dimension for lognormal realization
Pnmax = 1440

# file names
inputpkfname = 'pkG_rmax10000_b'+str(bias)+'.dat'	# input log-transformed power spectrum (format from calc_pkG)
num_data_inputpk = file_len(inputpkfname)
print('number of data in '+inputpkfname+' is '+str(num_data_inputpk))
inputmpkfname = 'pkG_rmax10000_b'+str(mbias)+'.dat'  # input log-transformed matter power spectrum (format from calc_pkG) 
num_data_inputmpk = file_len(inputmpkfname)												
print('number of data in '+inputmpkfname+' is '+str(num_data_inputmpk))								  
inputffname = 'wavenumber_fnu.txt'													
num_data_inputf = file_len(inputffname)													
print('number of data in '+inputffname+' is '+str(num_data_inputf))

DIR_RESULT = '/Users/shsaito/Desktop/tmp/'
do_weight_LAE = 2
numbin_log10L = 100
log10Lmin = 41.00
log10Lmax = 44.00
use_cpkG = 0
output_gal = 1
output_matter = 0

#random.seed(iseed)
#print 'iseed=',iseed
for i in range(0,nrealisation):
	seed1 = random.randint(1,100000)
	seed2 = random.randint(1,100000)
	seed3 = random.randint(1,100000)
	p1fname = DIR_RESULT+'params.gen_Poisson_mock_LogNormal_iseed'+str(iseed)+'_rlz'+str(i)
	Poissonfname = DIR_RESULT+'galaxies_iseed'+str(iseed)+'_rlz'+str(i)+'_Pnmax'+str(Pnmax)+'_b'+str(bias)+'.bin'
    #output lognormal galaxy catalog in (x,y,z,vx,vy,vz) [h^-1 Mpc] and [km/s]
	Densityfname = DIR_RESULT+'density_lognormal_is'+str(iseed)+'_rlz'+str(i)+'_Pnmax'+str(Pnmax)+'_b'+str(bias)+'.bin'
	executable('generate_Poisson/gen_Poisson_mock_LogNormal')([inputpkfname, inputmpkfname, use_cpkG, inputmpkfname, 
                                                               Lx,Ly,Lz,Pnmax,
                                                               do_weight_LAE, numbin_log10L, log10Lmin, log10Lmax, 
                                                               aHz,inputffname, bias, seed1, seed2, seed3,
                                                               output_gal, Poissonfname,
                                                               output_matter, Densityfname,
                                                               p1fname])
