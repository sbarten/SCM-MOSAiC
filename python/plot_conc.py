# --------------------------------------------------------------------------
# This python script reads the outfile '1d_chem.out', written by the 
# subroutine budget.f of the 1-D version of the coupled chemistry RACMO model.
# The outputfile can be found in the directory ../echam/case/output and shows 
# the concentrations of the different species calculated by the chemistry 
# routine versus time (x-axis) and altitude (y-axis), 
#
# written by: Laurens Ganzeveld, 27-11-98, updated for python 2019
# --------------------------------------------------------------------------

#!/usr/bin/env python
import os,sys,re,glob
import numpy
from subprocess import call
from numpy import *
from pylab import *
from matplotlib import colors
import time
import datetime
import matplotlib.dates as mdates

def makeZeroColormap(zero_point, base_map, num_segments=256):
        # base_map should be an existing colormap, such as cm.bwr
        # zero_value, between 0 and 1, should be where 0.5 in the old colormap should be mapped to
        output_cmap = dict.fromkeys(['red', 'green', 'blue'])
        red_list = []
        blue_list = []
        green_list = []
        lower_half = linspace(0.,zero_point,num_segments)
        for i, init_point in enumerate(linspace(0.,0.5,num_segments)):
            r,g,b,alpha = base_map(init_point)
            red_list.append((lower_half[i], r,r))
            green_list.append((lower_half[i],g,g))
            blue_list.append((lower_half[i],b,b))
        upper_half = linspace(zero_point, 1., num_segments)
        for i, init_point in enumerate(linspace(0.5,1.0,num_segments)):
            if i>0:
                r,g,b,alpha = base_map(init_point)
                red_list.append((upper_half[i], r,r))
                green_list.append((upper_half[i],g,g))
                blue_list.append((upper_half[i],b,b))
        output_cmap['red'] = tuple(red_list)
        output_cmap['green'] = tuple(green_list)
        output_cmap['blue'] = tuple(blue_list)
        return output_cmap

# declaration of some parameters
nstep0=0
nlevpr=0

# January 2004: define number of days of each month

ndaysmonth=zeros((13))
ndaysmonth[1]=31
ndaysmonth[2]=28
ndaysmonth[3]=31
ndaysmonth[4]=30
ndaysmonth[5]=31
ndaysmonth[6]=30
ndaysmonth[7]=31
ndaysmonth[8]=31
ndaysmonth[9]=30
ndaysmonth[10]=31
ndaysmonth[11]=30
ndaysmonth[12]=31

# determining if data should be read

readfile=-1
while int(readfile) < 0 or int(readfile) > 1: 
	readfile = raw_input('Reading of input data? (0=no,1=yes)  ')
	type(readfile)

if int(readfile) == 1:

# determining if absolute concentrations or the differences between
# two concentration fields are interpreted

  rfield=-1
  while int(rfield) < 0 or int(rfield) > 1:
  	rfield = raw_input('Do you want to see one field (0)'\
    	' or the difference between two fields (1)? ')
	type(rfield)

  fname1 = raw_input('Which file ? (default 1d_chem.out in directory: ../SCM/racmo/echam/case/output)')
  fname2='1d_chem.out'
  type(fname1)
  if fname1 == '':
    fname1=fname2

  if int(rfield) == 1:
    if int(readfile) == 1:
    	fname2 = raw_input('Which second file (1d_chem.out in other directory: ../SCM/racmo/echam/case/output)')
    	type(fname2)
    	if fname2 == '':
       		fname2='1d_chem.out'

  print(fname1)

# opening of data file #1: definition of default directory
  dir = '../case/output'
  print ' Reading the input file: ' ,dir+"/"+fname1

  n = 0
  i = 0
  nstart= 10
  nend=10
  npar=0

  with open(dir+"/"+fname1, "r") as my_file:
	  for line in my_file:
		  #print(line)
		  if n == 0:  # first line with some general parameters
      			  inp = line.split()
			  #print(inp)
			  nstep0=int(inp[0])
			  nlevatm1=int(inp[1])
			  nlevel1=int(inp[2])
			  iprval1=int(inp[3])
			  dtime1=float(inp[4])
			  #print(nstep0,nlevatm1,nlevel1,dtime1)
		  if n == 1: # second line with some general parameters
      			  inp = line.split()
			  ichemtype1=int(inp[0])
			  ntrac=int(inp[1])
		  if n == 2: # third line with more general parameters
      			  inp = line.split()
			  nlev1=int(inp[0])
			  nlev2=int(inp[1])
			  nstep=int(nstep0/iprval1+1)
			  nday=int(nstep0*dtime1/(3600.*24.))
			  nday=nday > 1
			  nstep_day=int((nstep-1)/nday)

			  print 'no. of timesteps',nstep0/iprval1,' and levels: ',nlevel1

                          # initialization of height, pressure and rhoa
			  height = zeros((nstep,nlevel1))
			  press = zeros((nstep,nlevel1))
			  rhoa = zeros((nstep,nlevel1))

			  jday = zeros((nstep))
			  gmttime = zeros((nstep))
			  loctime = zeros((nstep))  
			  hc = zeros((nstep))
			  trphgt = zeros((nstep))
			  pblhgt = zeros((nstep))

			  field =zeros((nstep,nlevel1,ntrac))
			  timestep =zeros((nstep))
                          ldateltime = ["" for x in range(nstep)]
                          ldate = ["" for x in range(nstep)]
			  ltime = ["" for x in range(nstep)]

                          nstart=n+1
			  nend=int(n+1+nlevel1/3)
			  #print(nstart,nend,npar)
			  nn=0

	   	  time.sleep(0.1)
   		  #print(inp[0])
		  print(n)
		  n = n + 1
		  if n == 3:
			  break

  # start reading the rest of the file doing the timestepping
  nlines = 0
  with open(dir+"/"+fname1, "r") as my_file:
          for _ in range(3):             # skipping the first three lines
        	  next(my_file)

	  for line in my_file:
		  # start the timestepping

		  #print(i,n,nstart,nend)
		  if (npar == 0 and n >= nstart and n <= nend):   # reading the heights of the levels: NOTE that the number oflines read is determined by the 
		                                                		  # previous nend statement. Since the heights/press/rhoa are written out in three columns
										  # is is nlevel1/3
			  inp = line.split()

			  for ii in range(0,3):
				  height[i,ii+nn*3]=float(inp[ii])
				  #print(i,ii+nn*3,height[i,ii+nn*3])

			  nn=int(nn+1)

			  if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+nlevel1/3)
				  npar=int(npar+1)
				  nn=0
				  #print(nstart,nend)

		  if (npar == 1 and n >= nstart and n <= nend):   # reading the pressure of the levels

			  inp = line.split()
			  for ii in range(0,3):
				  press[i,ii+nn*3]=float(inp[ii])
				  #print(ii+nn*3,press[i,ii+nn*3])

			  nn=int(nn+1)

			  if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+nlevel1/3)
				  npar=int(npar+1)
				  nn=0
				  #print(nstart,nend)

		  if (npar == 2 and n >= nstart and n <= nend):   # reading the rhoa of the levels

			  inp = line.split()
			  for ii in range(0,3):
				  rhoa[i,ii+nn*3]=float(inp[ii])
				  #print(ii+nn*3,rhoa[i,ii+nn*3])

          		  nn=int(nn+1)

                          if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+1)
				  npar=int(npar+1)
				  nn=0
				  #print(nstart,nend)

		  if (npar == 3 and n >= nstart and n <= nend):   # read the ldateltime

			  inp = line.split()  # reading the date/time information
			  ldate[i]=inp[0]
			  ltime[i]=inp[1]
                          ldateltime[i]=line[0:14]
			  timestep[i] = float(i)
	
                  	  nn=int(nn+1)

			  if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+1)
				  npar=int(npar+1)
				  nn=0

		  if (npar == 4 and n >= nstart and n <= nend):   # read the next line with time parameters and canopy height

			  inp = line.split() 

			  jday[i]=float(inp[1])
			  gmttime[i]=float(inp[2])
			  loctime[i]=float(inp[3])  
			  hc[i]=float(inp[4])

			  nn=int(nn+1)

                          if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+1)
				  npar=int(npar+1)
				  nn=0

		  if (npar == 5 and n >= nstart and n <= nend):   # read the next line with the tropopause and Bl height

			  inp = line.split() 

			  trphgt[i]=float(inp[0])
			  pblhgt[i]=float(inp[1])

			  nn=int(nn+1)

                          if (n == nend-1):
				  nstart=n+1
				  nend=int(n+1+ntrac)
				  npar=int(npar+1)
				  nn=0

		  if (npar == 6 and n >= nstart and n <= nend):   # and now really reading the tracer fields

			  inp = line.split()  
			  for j in range(0,nlevel1):
				  field[i,j,nn]=float(inp[j+2])
				  #print(j,nn,conc[i,j,nn])

                 	  nn=int(nn+1)

                          if (n == nend-1):
				  nstart=n+1
                                  nend=int(n+1+nlevel1/3)  # jumping back to read all the height levels
				  npar=0                   # resetting the parameter index, starting again reading the height, etc.
				  nn=0

# to check carefully the proper reading of the file: activate the following timer
		  #time.sleep(0.1)
   		  #print(inp[0])

		  print 'reading line #',n
		  n = n + 1
		  nlines = nlines + 1
		  if nlines == int(3*(nlevel1/3)+3+ntrac):   # going to next timestep when all the parameters are read
			  nlines = 0
			  i = i + 1
                  if i == nstep+1:
			  break


# opening of data file #2: definition of default directory
  if rfield == 1:
    dir = '../case/output'
    print ' Reading the input file: ' ,dir+"/"+fname2

    n = 0
    i = 0
    nstart= 10
    nend=10
    npar=0

    with open(dir+"/"+fname2, "r") as my_file:
	    for line in my_file:
		    #print(line)
		    if n == 0:  # first line with some general parameters
      			    inp = line.split()
			    #print(inp)
			    nstep02=int(inp[0])
			    nlevatm2=int(inp[1])
			    nlevel2=int(inp[2])
			    iprval2=int(inp[3])
			    dtime2=float(inp[4])
			    #print(nstep0,nlevatm1,nlevel1,dtime1)
		    if n == 1: # second line with some general parameters
      			    inp = line.split()
			    ichemtype2=int(inp[0])
			    ntrac2=int(inp[1])
		    if n == 2: # third line with more general parameters
      			    inp = line.split()
			    nlev12=int(inp[0])
			    nlev22=int(inp[1])
			    nstep2=int(nstep0/iprval1+1)
			    nday2=int(nstep0*dtime1/(3600.*24.))
			    nday2=nday > 1
			    nstep_day2=int((nstep-1)/nday)

			    print 'no. of timesteps',nstep0/iprval1,' and levels: ',nlevel1

                            # initialization of height, pressure and rhoa
			    height2 = zeros((nstep2,nlevel2))
			    press2 = zeros((nstep2,nlevel2))
			    rhoa2 = zeros((nstep2,nlevel2))

			    jday2 = zeros((nstep2))
			    gmttime2 = zeros((nstep2))
			    loctime2 = zeros((nstep2))  
			    hc2 = zeros((nstep2))
			    trphgt2 = zeros((nstep2))
			    pblhgt2 = zeros((nstep2))

			    field2 = zeros((nstep2,nlevel2,ntrac2))
			    timestep2 = zeros((nstep2))
                            ldateltime2 = ["" for x in range(nstep)]
                            ldate2 = ["" for x in range(nstep2)]
			    ltime2 = ["" for x in range(nstep2)]

                            nstart=n+1
			    nend=int(n+1+nlevel2/3)
			    #print(nstart,nend,npar)
			    nn=0

	   	    time.sleep(0.1)
   		    #print(inp[0])
		    print(n)
		    n = n + 1
		    if n == 3:
			    break

    # for two output files with different frequency of writing output 
    # other input files need to be read

    if iprval1*dtime1 != iprval2*dtime2:
      print 'Input files with results written away at a different time resolution'
      print 'Other input files must be read'
      stop

    # for two output files with results of different chemical schemes
    # other input files need to be read

    if ichemtype1 != ichemtype2:
      print 'Input files with results of different chemical schemes'
      print 'Other input files must be read'
      stop

    ip2=(nlev2+1)-nlev1

    # declaring new height with the dimension of the maximum number of levels

    ip=ip1
    if ip2 > ip1:
      ip=ip2
      nlev=ip1

    # start reading the rest of the file doing the timestepping
    nlines = 0
    with open(dir+"/"+fname2, "r") as my_file:
            for _ in range(3):             # skipping the first three lines
        	    next(my_file)

	    for line in my_file:
		    # start the timestepping

		    #print(i,n,nstart,nend)
		    if (npar == 0 and n >= nstart and n <= nend):   # reading the heights of the levels: NOTE that the number oflines read is determined by the 
		                                                		    # previous nend statement. Since the heights/press/rhoa are written out in three columns
										    # is is nlevel1/3
			    inp = line.split()

			    for ii in range(0,3):
				    height2[i,ii+nn*3]=float(inp[ii])
				    #print(i,ii+nn*3,height[i,ii+nn*3])

			    nn=int(nn+1)

			    if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+nlevel2/3)
				    npar=int(npar+1)
				    nn=0
				    #print(nstart,nend)

		    if (npar == 1 and n >= nstart and n <= nend):   # reading the pressure of the levels

			    inp = line.split()
			    for ii in range(0,3):
				    press2[i,ii+nn*3]=float(inp[ii])
				    #print(ii+nn*3,press[i,ii+nn*3])

			    nn=int(nn+1)

			    if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+nlevel2/3)
				    npar=int(npar+1)
				    nn=0
				    #print(nstart,nend)

		    if (npar == 2 and n >= nstart and n <= nend):   # reading the rhoa of the levels

			    inp = line.split()
			    for ii in range(0,3):
				    rhoa2[i,ii+nn*3]=float(inp[ii])
				    #print(ii+nn*3,rhoa[i,ii+nn*3])

          		    nn=int(nn+1)

                            if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+1)
				    npar=int(npar+1)
				    nn=0
				    #print(nstart,nend)

		    if (npar == 3 and n >= nstart and n <= nend):   # read the ldateltime

			    inp = line.split()  # reading the date/time information
			    ldate2[i]=inp[0]
			    ltime2[i]=inp[1]
                            ldateltime2[i]=line[0:14]
			    timestep2[i]=float(i)

                  	    nn=int(nn+1)

			    if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+1)
				    npar=int(npar+1)
				    nn=0

		    if (npar == 4 and n >= nstart and n <= nend):   # read the next line with time parameters and canopy height

			    inp = line.split() 

			    jday2[i]=float(inp[1])
			    gmttime2[i]=float(inp[2])
			    loctime2[i]=float(inp[3])  
			    hc2[i]=float(inp[4])

			    nn=int(nn+1)

                            if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+1)
				    npar=int(npar+1)
				    nn=0

		    if (npar == 5 and n >= nstart and n <= nend):   # read the next line with the tropopause and Bl height

			    inp = line.split() 

			    trphgt2[i]=float(inp[0])
			    pblhgt2[i]=float(inp[1])

			    nn=int(nn+1)

                            if (n == nend-1):
				    nstart=n+1
				    nend=int(n+1+ntrac)
				    npar=int(npar+1)
				    nn=0

		    if (npar == 6 and n >= nstart and n <= nend):   # and now really reading the tracer fields

			    inp = line.split()  
			    for j in range(0,nlevel2):
				    field2[i,j,nn]=float(inp[j+2])
				    #print(j,nn,conc[i,j,nn])

                 	    nn=int(nn+1)

                            if (n == nend-1):
				    nstart=n+1
                                    nend=int(n+1+nlevel2/3)  # jumping back to read all the height levels
				    npar=0                   # resetting the parameter index, starting again reading the height, etc.
				    nn=0

   # to check carefully the proper reading of the file: activate the following timer
		    #time.sleep(0.1)
   		    #print(inp[0])

		    print 'reading line #',n
		    n = n + 1
		    nlines = nlines + 1
		    if nlines == int(3*(nlevel2/3)+3+ntrac2):   # going to next timestep when all the parameters are read
			    nlines = 0
			    i = i + 1
                    if i == nstep+1:
			    break


  # end reading of input file(s)      
  ################################################################
  # definition of the specifics of the chemistry scheme being used

  ichemtype=ichemtype1

  # definition of extra number of parameters which are being calculated
  # from the parameters in the file

  dparam=3 
  iccn=0
  isoa=0

  # ESS_lg_20100728+ definition of number of aerosol species and some other properties
  ntrac_aer=0
  if ((ichemtype == 2 & ntrac == 70) or (ichemtype == 5 & ntrac == 193)):
     ntrac_aer=4
     print('The number of aerosol scheme tracers is: ',ntrac_aer)
     print('ENTER to continue')
     read,a

  ntrac_gas=ntrac+ntrac_aer

  # ESS_lg_20100728-
  # definition of no. of species and assigning the names and numbers

  if (ichemtype == 1):
    nspec=ntrac+dparam
    ntract=ntrac+dparam
    spec = ['']*nspec
    spec[0:nspec-1]=\
    ['NX-FAM ','O3','CH4','CO','HNO3','H2O2',
     'CH3O2H','DMS','SO2','SO4','RADON','NO','NO2','NO3',
     'N2O5','HNO4','OH','HO2','O3S','O1D','CH3O',
     'CH3O2','CH2O','H+','NOX','NOY','NOx/NOX']

    inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
    idms=7;iso2=8;iso4=9;irad=10;ino=11;ino2=12;ino3=13;in2o5=14;
    ihno4=15;ioh=16;iho2=17;io3s=18;io1d=19;ich3o=20;
    ich3o2=21;ich2o=22;ihplus=23;inox=24;inoy=25;inoxinoy=26

  if (ichemtype == 2):
    nspec=ntrac+dparam
    ntract=ntrac+dparam
    ntrac_gas=ntract-ntrac_aer
    spec = ['']*nspec
    spec = ['']*nspec

  # new scheme, including the CO2

    if ntrac == 50:
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','NO','NO2','NO3','N2O5','HNO4','OH',
      'HO2','CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3',
      'H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      ino=34;ino2=35;ino3=36;in2o5=37;ihno4=38;ioh=39;iho2=40;
      ich3o2=41;ic2o3=42;ixo2=43;iror=44;ixo2n=45;irxpar=46;
      ibxo2n=47;imc3o3=48;ihplus=49;inox=50;inoy=51;inoxinoy=52
  # not in this chemistry scheme
      io1d=53

  # new scheme, including the monoterpenes

    if ntrac == 52:
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','NO','NO2','NO3','N2O5',
      'HNO4','OH','HO2','CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N',
      'MC3O3','H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;ino=36;ino2=37;ino3=38;in2o5=39;ihno4=40;
      ioh=41;iho2=42;ich3o2=43;ic2o3=44;ixo2=45;iror=46;ixo2n=47;irxpar=48;
      ibxo2n=49;imc3o3=50;ihplus=51;inox=52;inoy=53;inoxinoy=54
  # not in this chemistry scheme
      io1d=55

  # new scheme, including CH3CL and CHCL3

    if ntrac == 54 :
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','CH3CL','CHCL3',
      'NO','NO2','NO3','N2O5','HNO4','OH','HO2','CH3O2','C2O3','XO2',
      'ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;ich3cl=36;ichcl3=37;ino=38;ino2=39;ino3=40;
      in2o5=41;ihno4=42;ioh=43;iho2=44;ich3o2=45;ic2o3=46;ixo2=47;
      iror=48;ixo2n=49;irxpar=50;ibxo2n=51;imc3o3=52;ihplus=53;inox=54;
      inoy=55;inoxinoy=56
  # not in this chemistry scheme
      io1d=57

  # new scheme, including sesquiterpenes

    if ntrac == 55 :
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'NO','NO2','NO3','N2O5','HNO4','OH','HO2','CH3O2','C2O3','XO2',
      'ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ino=39;ino2=40;
      ino3=41;in2o5=42;ihno4=43;ioh=44;iho2=45;ich3o2=46;ic2o3=47;ixo2=48;
      iror=49;ixo2n=50;irxpar=51;ibxo2n=52;imc3o3=53;ihplus=54;inox=55;
      inoy=56;inoxinoy=57
  # not in this chemistry scheme
      io1d=58

  # new scheme, including HONO

    if ntrac == 56 :
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','NO','NO2','NO3','N2O5','HNO4','OH','HO2','CH3O2','C2O3','XO2',
      'ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ino=40;ino2=41;ino3=42;in2o5=43;ihno4=44;ioh=45;iho2=46;ich3o2=47;
      ic2o3=48;ixo2=49;iror=50;ixo2n=51;irxpar=52;ibxo2n=53;imc3o3=54;
      ihplus=55;inox=56;inoy=57;inoxinoy=58
  # not in this chemistry scheme
      io1d=59

  # new scheme, including hydroxy-methyl/alkylhydroxyperoxide (HMHP/HAHP)

    if ntrac == 58 :
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;
      ino=42;ino2=43;ino3=44;in2o5=45;ihno4=46;ioh=47;iho2=48;ich3o2=49;
      ic2o3=50;ixo2=51;iror=52;ixo2n=53;irxpar=54;ibxo2n=55;imc3o3=56;
      ihplus=57;inox=58;inoy=59;inoxinoy=60
  # not in this chemistry scheme
      io1d=61

  # new scheme, including CH3OH and CH3CN

    if ntrac == 60 :
      dparam=9 # mz_lg_20060102+ added MVK+METHAC, organic peroxides and some others
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','CH3OH','CH3CN','NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOX',
      'MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;ich3oh=42;ich3cn=43;
      ino=44;ino2=45;ino3=46;in2o5=47;ihno4=48;ioh=49;iho2=50;ich3o2=51;
      ic2o3=52;ixo2=53;iror=54;ixo2n=55;irxpar=56;ibxo2n=57;imc3o3=58;
      ihplus=59;inox=60;inoy=61;inoxinoy=62;imvkmethac=63;
      iisopprodsum=64;imvkmethacrat=65;iisopprodrat=65;
      iprodisoprat=67;iorgperox=68;

  # not in this chemistry scheme
      io1d=69

  # new scheme, including 2- and 1- double bond sesquiterpenes and other
  # monoterpenes

    if ntrac == 63 :
      dparam=9 # mz_lg_20060102+ added MVK+METHAC, organic peroxides and some others
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','CH3OH','CH3CN','SQTERP2B','SQTERP1B',
      'MTTERP','NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOy',
      'MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;ich3oh=42;ich3cn=43;isqterp2b=44;isqterp1b=45;
      imtterp=46;ino=47;ino2=48;ino3=49;in2o5=50;ihno4=51;ioh=52;iho2=53;ich3o2=54;
      ic2o3=55;ixo2=56;iror=57;ixo2n=58;irxpar=59;ibxo2n=60;imc3o3=61;
      ihplus=62;inox=63;inoy=64;inoxinoy=65;imvkmethac=66;
      iisopprodsum=67;imvkmethacrat=68;iisopprodrat=69;
      iprodisoprat=70;iorgperox=71;

  # not in this chemistry scheme
      io1d=72

  # ESS_lg_20100628+ new scheme, including a selection of anthropogenic alkenes

    if ntrac == 66 :
      dparam=9 # mz_lg_20060102+ added MVK+METHAC, organic peroxides and some others
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','CH3OH','CH3CN','SQTERP2B','SQTERP1B',
      'MTTERP','HEXANE','BUTADIENE','TMBENZENE','NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOy',
      'MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;ich3oh=42;ich3cn=43;isqterp2b=44;isqterp1b=45;
      imtterp=46;ihexane=47;ibutadiene=48;itmbenzene=49;ino=50;ino2=51;ino3=52;in2o5=53;
      ihno4=54;ioh=55;iho2=56;ich3o2=57;ic2o3=58;ixo2=59;iror=60;ixo2n=61;irxpar=62;
      ibxo2n=63;imc3o3=64;ihplus=65;inox=66;inoy=67;inoxinoy=68;imvkmethac=69;
      iisopprodsum=70;imvkmethacrat=71;iisopprodrat=72;iprodisoprat=73;iorgperox=74;

  # not in this chemistry scheme
      io1d=75

  # MAQ_lg_20180717+ new scheme, including COS

    if ntrac == 67 :
      dparam=9 # mz_lg_20060102+ added MVK+METHAC, organic peroxides and some others
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','CH3OH','CH3CN','SQTERP2B','SQTERP1B',
      'MTTERP','HEXANE','BUTADIENE','TMBENZENE','COS','NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+','NOx','NOY','NOx/NOy',
      'MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;ich3oh=42;ich3cn=43;isqterp2b=44;isqterp1b=45;
      imtterp=46;ihexane=47;ibutadiene=48;itmbenzene=49;icos=50;ino=51;ino2=52;ino3=53;in2o5=54;
      ihno4=55;ioh=56;iho2=57;ich3o2=58;ic2o3=59;ixo2=60;iror=61;ixo2n=62;irxpar=63;
      ibxo2n=64;imc3o3=65;ihplus=66;inox=67;inoy=68;inoxinoy=69;imvkmethac=70;
      iisopprodsum=71;imvkmethacrat=72;iisopprodrat=73;iprodisoprat=74;iorgperox=75;

  # not in this chemistry scheme
      io1d=76

  # ESS_lg_20120117+ new scheme, including a selection of anthropogenic alkenes and SOA species 

    if ntrac == 70 :
      isoa=1
      dparam=9 # mz_lg_20060102+ added MVK+METHAC, organic peroxides and some others
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','CO2','MATERP','MBTERP','SQTERP','CH3CL','CHCL3',
      'HONO','CH2OHO2H','RCHOHO2H','CH3OH','CH3CN','SQTERP2B','SQTERP1B',
      'MTTERP','HEXANE','BUTADIENE','TMBENZENE',
      'APINp1g','APINp1a','APINp2g','APIN2p2a',    # ESS_lg_20120118+
      'NO','NO2','NO3','N2O5','HNO4','OH','HO2',
      'CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3','H+',
      'NOx','NOY','NOx/NOy','MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;ico2=33;
      imaterp=34;imbterp=35;isqterp=36;ich3cl=37;ichcl3=38;ihono=39;
      ich2oho2h=40;irchoho2h=41;ich3oh=42;ich3cn=43;isqterp2b=44;isqterp1b=45;
      imtterp=46;ihexane=47;ibutadiene=48;itmbenzene=49;
      iapinp1g=50;iapinp1a=51;iapinp2g=52;iapinp2a=53;
      ino=54;ino2=55;ino3=56;in2o5=57;
      ihno4=58;ioh=59;iho2=60;ich3o2=61;ic2o3=62;ixo2=63;iror=64;ixo2n=65;irxpar=66;
      ibxo2n=67;imc3o3=68;ihplus=69;
      inox=70;inoy=71;inoxinoy=72;imvkmethac=73;
      iisopprodsum=74;imvkmethacrat=75;iisopprodrat=76;iprodisoprat=77;iorgperox=78;

  # not in this chemistry scheme
      io1d=79

  # ESS_lg_20100728-
      AVOug=1.e12/6.023e23
      moleweight=fltarr(ntrac_gas)
      moleweight[iapinp1g:iapinp2a]=[170.,170.,170.,170.] 
      molec2ug=fltarr(ntrac_gas)
      molec2ug=moleweight*AVOug
  # ESS_lg_20100728-

  # old scheme, without the CO2

    if ntrac == 49 :
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','CH2O','O3S','ALD2','PAR','OLE','ETH',
      'PAN','ACET','ISOP','MGLY','ISOPRD','METHAC','MVK','MEK','MPAN',
      'NTR','DMS','SO2','SO4','RADON','ISONTR','HCOOH','CH3CO2H',
      'NH2','NH3','NH4','NO','NO2','NO3','N2O5','HNO4','OH',
      'HO2','CH3O2','C2O3','XO2','ROR','XO2N','RXPAR','BXO2N','MC3O3',
      'H+','NOx','NOY','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      ich2o=7;io3s=8;iald2=9;ipar=10;iole=11;ieth=12;
      ipan=13;iacet=14;iisop=15;imgly=16;iisoprd=17;
      imethac=18;imvk=19;imek=20;impan=21;intr=22;
      idms=23;iso2=24;iso4=25;irad=26;iisontr=27;
      ihcooh=28;ich3co2h=29;inh2=30;inh3=31;inh4=32;
      ino=33;ino2=34;ino3=35;in2o5=36;ihno4=37;ioh=38;iho2=39;
      ich3o2=40;ic2o3=41;ixo2=42;iror=43;ixo2n=44;irxpar=45;
      ibxo2n=46;imc3o3=47;ihplus=48;inox=49;inoy=50;inoxinoy=51
  # not in this chemistry scheme
      io1d=52;ico2=53

  if (ichemtype == 3) :
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['Nx-FAM','O3','CH4','CO','HNO3','H2O2',
      'CH3O2H','DMS','SO2','SO4','ROR','XO2',
      'C2O3','PAN','XO2N','ALD2','MGLY','PAR','ETH',
      'OLE','ISOP','ROOH','RXPAR','ORGNTR','RADON','NO','NO2','NO3',
      'N2O5','HNO4','OH','HO2','O3S','O1D','CH3O','CH3O2','CH2O',
      'H+','NOx','NOy','NOx/NOX']
      inxfam=0;io3=1;ich4=2;ico=3;ihno3=4;ih2o2=5;ich3o2h=6;
      idms=7;iso2=8;iso4=9;iror=10;ixo2=11;ic2o3=12;
      ipan=13;ixo2n=14;iald2=15;imgly=16;ipar=17;
      ieth=18;iole=19;iisop=20;irooh=21;irxpar=22;iorgntr=23;irad=24;
      ino=25;ino2=26;ino3=27;in2o5=28;ihno4=29;ioh=30;
      iho2=31;io3s=32;io1d=33;ich3o=34;ich3o2=35;ich2o=36;
      ihplus=37;inox=38;inoy=39;inoxinoy=40
  # not in this chemistry scheme
      io1d=50

  if (ichemtype == 4) :
      nspec=ntrac+dparam
      ntract=ntrac+dparam
      spec = ['']*nspec
      spec[0:nspec-1]=\
      ['dummy1','dummy2','CH4','ACET','N2O5','CH3OH',
      'HCOOH','H2O2','CH3COOH','HNO4','CH3CO3H','PAN',
      'CH3CO2H','HNO3','MPAN','ACETOOH','ISOOH','ACETO2',
      'MACROOH','CO','ONIT','CH2O','HACET','ISOP','MACR','MGLY',
      'C2O3','NO2','CH3O2','OH','IISO2','MACRO2',
      'NO3','HO2','O3','NO','Nx-FAM','DMS','SO2','SO4',
      'RADON','H+','NOx','NOy','NOx/NOX']
      idummy1=0;idummy2=1;ich4=2;iacet=3;in2o5=4;ich3oh=5;ihcooh=6;
      ih2o2=7;ich3cooh=8;ihno4=9;ich3co3h=10;ipan=11;ich3co2h=12;
      ihno3=13;impan=14;iacetooh=15;iisooh=16;iaceto2=17;
      imacrooh=18;ico=19;ionit=20;ich2o=21;ihacet=22;iisop=23;
      imacr=24;imgly=25;ic2o3=26;ino2=27;ich3o2=28;ioh=29;
      iiso2=30;imacro2=31;ino3=32;iho2=33;io3=34;ino=35;
      inxfam=36;idms=37;iso2=38;iso4=39;irad=40;ihplus=42;inox=43;inoy=44;
      inoxinoy=45
  # not in this chemistry scheme
      io1d=50

  # uncomment to test reading/assignment of parameters for MECCA chemistry code
  #ichemtype=5
  #ntrac=266

  # ESS_lg_20100228+ added MECCA's chemistry scheme
  if (ichemtype == 5) :
     dparam=10 # ESS_lg_20100228+ added extra diagnostic tracers, see below
     nspec=ntrac+dparam
     ntract=ntrac+dparam
     spec = ['']*nspec

  #  opening and reading of mecca species file

     speciesindx= ['']*ntrac
     speciesname= ['']*ntrac
     spec= ['']*ntrac

  # the next file should resemble the one found in ../work/mecca1/v** but additional tracers can be added that are
  # not explicitly considered in MECCA1

     if (ntrac == 74):  
	 fname='mecca1/species_mecca1v19c_NTR74.inc'
     if (ntrac == 189):
	 fname='mecca1/species_mecca1v19c_NTR189.inc'
     if (ntrac == 251):
	 fname='mecca1/species_mecca1v19c_NTR251.inc'
     if (ntrac == 266):
	 fname='mecca1/species_mecca1v19c_NTR266.inc'

     jt=0
     with open(fname, "r") as my_file:
       for _ in range(1):             # skipping the first line
	 next(my_file)

       for line in my_file:
	 inp = line.split()
	 speciesindx[jt]=inp[0]
	 speciesname[jt]=inp[2]

  #   commands to reduce read-in tracer name to usable format and assign index 
  #   thanks to Joost en de Brugh (meteorologie, March 2010)

	 speciesname_length = len(speciesname[jt])
	 spec[jt] =speciesname[jt][1:-1]  # to remove the accent at the beginning and the end of the speciesname
	 #print(spec[jt])

	 #specind='i'+string(format='(1a)',spec[jt])
	 #void = Execute("i"+spec[jt]+' = '+string[jt])

       jt=jt+1

       spec[ntrac:ntract-1]=\
       ['NOx','NOY','NOx/NOy',
       'MVK+METHAC','ISOP+MVK+METHAC','MVK/METHAC','ISOP/(MVK+METHAC)','(MVK+METHAC)/ISOP','CO2H','Radon']

       inox=ntrac+1
       inoy=ntrac+2
       inoxinoy=ntrac+3
       imvkmethac=ntrac+4
       iisopprodsum=ntrac+5
       imvkmethacrat=ntrac+6
       iisopprodrat=ntrac+7
       iprodisoprat=ntrac+8
       iorgperox=ntrac+9
       irad=ntrac+10

  ##  commands to produce the new tracer listing file including plot options; this 
  ##  file is read in later on in this script.
  #
  #   inew_tracer_file=0
  #   if inew_tracer_file == 1 :
  #   
  ##    writing print parameters
  #
  #     openw,2,'../python/parameter_files/1d_chem_mecca1.par'
  #
  #     noyes=intarr[ntrac]
  #     loga=intarr(ntrac)
  #     minval=fltarr(ntrac)
  #     maxval=fltarr(ntrac)
  #     xtic=fltarr(ntrac)
  #     minmax=fltarr(ntrac)
  #     title=strarr(ntrac)
  #     unit=strarr(ntrac)
  #
  #     irow=1
  #     ikol=1
  #     xps=1
  #     eps=0
  #     pola=0
  #
  #     printf,2,irow,ikol
  #     printf,2,xps,eps
  #     printf,2,pola 
  #     for jt in range(0,ntrac-1):
  #       noyes(jt)=1
  #       loga(jt)=0
  #       minval(jt)=0
  #       maxval(jt)=1
  #       xtic(jt)=4
  #       minmax(jt)=1
  #       title(jt)=spec(jt)
  #       unit(jt)='[ppbv]'
  #       if (jt ge ntrac_gas) then unit(jt)='[ugm3]'
  #
  #       printf,2,format='(i1,1x,i1,1x,f4.1,1x,f4.1,1x,i1,1x,i1,1x,1a,1x,1a)',$
  #         noyes(jt),loga(jt),minval(jt),maxval(jt),xtic(jt),minmax(jt),$
  #         title(jt),unit(jt)
  #     end
  #     close,2
  #   end
  #
  # end

  # end reading of chemistry scheme specifics  
  ##########################################################################
  # definition of some general parameters and recalculatio of concentrations

  # air mass, avogadro
  zma=28.9644
  avo=6.022e23
  recalc=1e9*zma/avo
  halflife_radon=5.5*86400 # halflife expressed in seconds
  Rdecay_radon=0.69/halflife_radon # 0.69=ln2

  # recalculation of the read in concentrations (molecules cm-3) to the units for plotting, 
  # generally this is mixing ratio

  conc =zeros((nstep,nlevel1,ntract))

  for jt in range (0,ntrac-1): 
     if jt != ioh and jt != iho2 and jt != io1d and jt != irad and jt != ico2:
       conc[:,:,jt]=recalc*field[:,:,jt]
     if jt == ioh or jt == iho2 or jt == io1d or jt == irad:
       conc[:,:,jt]=field[:,:,jt]

  # ESS_lg_20121216+ radon fieldentration in mBq m-3
     if jt == irad: 
       conc[:,:,jt]=field[:,:,jt]*1e3*1e6*Rdecay_radon

  # ESS_lg_20100412+ added CCN
     if ichemtype == 5 and jt == iccn:
       conc[:,:,jt]=field[:,:,jt]

  # ESS_lg_20100728+ added the SOA aerosols
     if isoa == 1 :
       if (jt == iapinp1g or jt == iapinp1a or jt == iapinp2g or jt == iapinp2a):
	 conc[:,:,jt]=field[:,:,jt]*molec2ug(jt)

  # recalculation of OH, O1D and HO2 concentrations to plottable 
  # order of magnitude

  # ESS_lg_20100626+ plotting OH, HO2 default in molecules cm-3
  # if not then set the switch to 0

  plot_OHHO2_molec=1
  if plot_OHHO2_molec == 1:
     print 'The OH and HO2 concentrations are plotted in molec cm-3: '
     print 'Modify the value of plot_OHHO2_molec' \
           ' (0 for mixing ratio, 1 for molec.)'

  # ESS_lg_20100626+ modified
  if plot_OHHO2_molec == 1:
      conc[:,:,ioh]=1e-6*conc[:,:,ioh]
      conc[:,:,iho2]=1e-8*conc[:,:,iho2]
  if plot_OHHO2_molec != 1:
      conc[:,:,ioh]=recalc*conc[:,:,ioh]/rhoa[:,:]
      conc[:,:,iho2]=recalc*conc[:,:,iho2]/rhoa[:,:]

  if (ichemtype1 == 1):
      conc[:,:,io1d]=conc[:,:,io1d]

  # calculating the NOx concentration 

  conc[:,:,inox]=conc[:,:,ino]+conc[:,:,ino2]

  # calculating the NOy concentration

  if (ichemtype == 1):
      conc[:,:,inoy]=conc[:,:,inxfam]+conc[:,:,ihno3]
  if (ichemtype == 2):
      conc[:,:,inoy]=conc[:,:,inxfam]+conc[:,:,ihno3]+conc[:,:,ipan]
  if (ichemtype == 3):
      conc[:,:,inoy]=conc[:,:,inxfam]+conc[:,:,ihno3]+conc[:,:,ipan]+conc[:,:,iorgntr]

  # calculating the NOx to NOX ratio

  conc[:,:,inoxinoy]=conc[:,:,inox]/conc[:,:,inoy]

  # mz_lg_20050102+ adding MVK+METHAC

  if (ichemtype1 == 2 and ntrac == 60):
      conc[:,:,imvkmethac]=conc[:,:,imvk]+conc[:,:,imethac]
      conc[:,:,iisopprodsum]=conc[:,:,iisop]+conc[:,:,imvk]+conc[:,:,imethac]
      conc[:,:,iisopprodrat]=conc[:,:,iisop]/(conc[:,:,imvk]+conc[:,:,imethac])
      conc[:,:,iorgperox]=conc[:,:,ich3o2h]+conc[:,:,ich2oho2h]+conc[:,:,irchoho2h]

  if (ichemtype1 == 2 and ntrac == 63):
      conc[:,:,imvkmethac]=conc[:,:,imvk]+conc[:,:,imethac]
      conc[:,:,iisopprodsum]=conc[:,:,iisop]+conc[:,:,imvk]+conc[:,:,imethac]
      conc[:,:,imvkmethacrat]=conc[:,:,imvk]/conc[:,:,imethac]
      conc[:,:,iisopprodrat]=conc[:,:,iisop]/(conc[:,:,imvk]+conc[:,:,imethac])
      conc[:,:,iprodisoprat]=(conc[:,:,imvk]+conc[:,:,imethac])/conc[:,:,iisop]
      conc[:,:,iorgperox]=conc[:,:,ich3o2h]+conc[:,:,ich2oho2h]+conc[:,:,irchoho2h]

################################################################
# start of the actual plotting

# March 2002, assigning the local actual date, day, month and year, and 
# time for use in defining plotting domain

lday= ['']*nstep
lmonth= ['']*nstep
lyear= ['']*nstep
lhr= ['']*nstep
lmin= ['']*nstep
time=zeros((nstep))
sec=zeros((nstep))

for i in range (0,nstep):
  lday[i]=ldateltime[i][0:2]
  lmonth[i]=ldateltime[i][3:5]
  lyear[i]=ldateltime[i][6:9]
  lhr[i]=ldateltime[i][9:11]
  lmin[i]=ldateltime[i][12:14]

## ESS_lg_12112009+ correction to plot in GMT time
#  if itimeGMT eq 1 then begin
#    lhr(ic)=lhr(ic)+dtimeGMT
#    if (lhr(ic) ge 24) then begin
#      lhr(ic)=lhr(ic)-24
#      lday(ic)=lday(ic)+1
#      jday(ic)=jday(ic)+1
#     endif
#  endif
## ESS_lg_12112009-
#
#  lmin(ic)=strmid(ldatltime(ic),12,2)

# calculation of total number of days of simulation, the hour and the minutes

  if i == 0:
     nday0=jday[i]
     nyear0=lyear[i] # added 01102004
  tday=jday[i]-nday0

  if i == 0:
    tdayold=tday

  hr=int(loctime[i]/100)
  minute=loctime[i]-hr*100
  time[i]=hr+(minute/60)

# calculation of time

# January 2004, 01102004 further modified adding year: 
  sec[i]=(float(lyear[i])-float(nyear0))*365*24*3600+ \
    float(jday[i])*24*3600+ \
    float(lhr[i])*3600.+    \
    float(lmin[i])*60.  

# correction for proper time calcuation in case of timestep < 60 s

  if i > 0 and dtime1 < 60. and sec[i] == 0. and hr == 0 and tdayold == tday:
    tday=tday+1

#   January 2004, 01102004 further modified adding year: 
    sec[i]=(float(lyear[i])-float(nyear0))*365*24*3600+ \
      float(jday[i]+1.)*24*3600+ \
      float(lhr[i])*3600.+    \
      float(lmin[i])*60.  

# definition of the time interval for plotting

# March 2002, modified
tintval=-1
while int(tintval) < 0 or int(tintval) > 2:
   tintval = raw_input('Whole integration period (0), a specific timeframe (1) or do you want to see an average diurnal cycle (2) ')
   type(tintval)

ts1=0
ts2=nstep-1

if int(tintval) == 0:
  starttime=' '
  starttime=ldateltime[0]
  endtime=' '
  endtime=ldateltime[nstep-1]

# ESS_lg_20111231+ added tintval eq 2 
if int(tintval) == 1 or int(tintval) == 2:
  print 'Which period of the simulation starting at dd:mm:yr hr:min ',\
        lday[0],':',lmonth[0],':',lyear[0],' ',lhr[0],':',lmin[0]
  print '                             end ending at dd:mm:yr hr:min ',\
        lday[nstep-1],':',lmonth[nstep-1],':',lyear[nstep-1],' ',lhr[nstep-1],':',lmin[nstep-1]

  print 'Current starttime is: %s ' %(starttime) 
  print 'Current endtime is  : %s ' %(endtime) 
  enter = raw_input('Type 1 to change this: ')
  type(enter)
  if (int(enter) == 1):
     starttime=raw_input('Give begin time (dd:mm:yr hr:mm): ')
     type(starttime)
     endtime=raw_input('Give end time (dd:mm:yr hr:mm): ')
     type(endtime)

startlday=starttime[0:2]
startlmonth=starttime[3:5]
startlyear=starttime[6:9]
startlhr=starttime[9:11]
startlmin=starttime[12:14]

# January 2004
nmonth=0
secmonth=0
nmonth=int(startlmonth)

for i in range(1,nmonth):
   secmonth=secmonth+ndaysmonth[i]*24.*3600.+ \
     (float(startlyear)-float(nyear0))*365*24.*3600. 

startsec=secmonth+float(startlday)*24*3600+ \
      float(startlhr)*3600.+float(startlmin)*60. 

endlday=endtime[0:2]
endlmonth=endtime[3:5]
endlyear=endtime[6:9]
endlhr=endtime[9:11]
endlmin=endtime[12:14]

# January 2004
nmonth=0
secmonth=0
nmonth=int(endlmonth)
for i in range(1,nmonth):
   secmonth=secmonth+ndaysmonth[i]*24.*3600.+ \
     (float(endlyear)-float(nyear0))*365*24.*3600. 

endsec=secmonth+float(endlday)*24*3600+ \
      float(endlhr)*3600.+float(endlmin)*60. 

if int(tintval) == 2 and tday < 1:
  print 'too short integration period to construct diurnal cycle: select longer period !!'

if int(tintval) == 0:
  istep = nstep
if int(tintval) == 1:
  istep = nstep/nday+1

if int(tintval) == 0:
  ts1=0
  ts2=istep-1

# ESS_lg_20111231+ added tintval eq 2
if int(tintval) == 1 or int(tintval) == 2:
  for ii in range(1,nstep):

     print(ii-1)
# ESS_lg_20120103+ modified and added ic1 and ic2
     if sec[ii] > startsec and sec[ii-1] <= startsec:
       ts1=ii
     if sec[ii] > endsec and sec[ii-1] <= endsec:
       ts2=ii-1
# ESS_lg_20120103-

istep=int(ts2-ts1)

# MAQ_lg_20190812+ not yet introduced: see plot_conc.pro to also have the following option!
# May 2008# introduced the option to show distance axis instead


# asking which species to plot

jjt=0
trc=-1
titla=['']*int(ntract)

while int(trc) < 0 or int(trc) > jjt:
  for jt in range (0,ntract-1):
     jjt=jjt+1
     print 'no.: ',jjt,': ',spec[jt]

     #  defining plot title
     titla[jjt-1]=spec[jt]

  trc = raw_input('Give tracer number: ')
  type(trc)

# assigning jt 
jt=int(trc)-1

# asking for which plot types the user wants to see

choice=-1
while int(choice) < 0 or int(choice) > 2:
  print 'Would you like to see the time series at each specific height (0) or'
  print 'Would you like to see vertical profiles for a specific time (1) or'
  choice=raw_input('Would you like to see diurnal contours (2) ')
  type(choice)

# definition of the height interval
nlev=nlevel1
if int(choice) == 0:
   for i in range (nlev-1,0,-1):
     ii = int(i-1)
     print 'level: %i alt.: %8.1f' % (i,height[0,ii])

   print'Which levels ? (type lowest level and highest level, max=8 for L19, max=12 for L60)'
   levh1=raw_input('lowest level: ')
   type(levh1)
   levh1=int(levh1)
   levh2=raw_input('highest level: ')
   type(levh2)
   levh2=int(levh2)

   if nlev < 21 and int(levh2)-int(levh1) > 8:
     levh2=int(levh1)+8
   if nlev < 62 and int(levh2)-int(levh1) > 12:
     levh2=int(levh1)+12

   nlevpr=int(levh2)-int(levh1)+1

# plotting vertical profiles, not yet updated in python
#
# if choice eq 1 then begin
#
#   ; February 2006
#   nprofmax=1
#   j12: print,'Do you want to see one profile (0) ',$
#              'all profiles for a specific timeframe (1) (with a maximum of 12) ',$
#              'or the timeframe average profiles (2)?'
#   read,iprof_time
#   if iprof_time lt 0 or iprof_time gt 2 then goto, j12
#
#   if iprof_time eq 0 then begin
#     j10:print,'Which time of the simulation starting at dd:mm:yr hr:min ',$
#          lday(0),':',lmonth(0),':',lyear(0),' ',lhr(0),':',lmin(0)
#        print,'                             end ending at dd:mm:yr hr:min ',$
#          lday(ic-1),':',lmonth(ic-1),':',lyear(ic-1),' ',lhr(ic-1),':',lmin(ic-1)
#     print,'Give time in (dd:mm:yr hr:mm): '
#     proftime=' ' ; December 2005
#     read,proftime
#     proflday=strmid(proftime,0,2)
#     proflmonth=strmid(proftime,3,2)
#     proflyear=strmid(proftime,6,2)
#     proflhr=strmid(proftime,9,2)
#     proflmin=strmid(proftime,12,2)
#
#;    December 2005, modified calculation of prof sec, also due to the diff. no.
#;    of days in months
#     ndays=0
#     for i=1,lmonth(0) do begin
#       ndays=ndays+ndaysmonth(i)
#     endfor
#     profsec=ndays*86400+ $
#           float(proflday)*24*3600+ $
#           float(proflhr)*3600.+float(proflmin)*60. 
#
#;    if time lt 1 or time gt nstep then goto,j10
#
#     for i=1,ic-1 do begin
#      if sec(i) ge profsec and sec(i-1) lt profsec then begin
#        tsprof=i 
#        print,i,sec(i),profsec,tsprof
#	read,a
#      endif
#      tsprof_end=tsprof  ; February 2006
#     endfor
# 
#   end else begin
#
#     ; February 2006
#     nprofmax=nstep
#     nprofmax_plot=12
#     j14:print,'Which time of the simulation starting at dd:mm:yr hr:min ',$
#           lday(0),':',lmonth(0),':',lyear(0),' ',lhr(0),':',lmin(0)
#        print,'                             end ending at dd:mm:yr hr:min ',$
#           lday(ic-1),':',lmonth(ic-1),':',lyear(ic-1),' ',lhr(ic-1),':',lmin(ic-1)
#     print,'Give start time in (dd:mm:yr hr:mm): '
#     profstarttime=' ' ; December 2005
#     read,profstarttime
#     print,'Give end time in (dd:mm:yr hr:mm): '
#     profendtime=' ' ; December 2005
#     read,profendtime
#
#     proflstartday=strmid(profstarttime,0,2)
#     proflstartmonth=strmid(profstarttime,3,2)
#     proflstartyear=strmid(profstarttime,6,2)
#     proflstarthr=strmid(profstarttime,9,2)
#     proflstartmin=strmid(profstarttime,12,2)
#;    December 2005, modified calculation of prof sec, also due to the diff. no.
#;    of days in months
#     ndays=0
#     for i=1,lmonth(0) do begin
#       ndays=ndays+ndaysmonth(i)
#     endfor
#     profstartsec=ndays*86400+ $
#           float(proflstartday)*24*3600+ $
#           float(proflstarthr)*3600.+float(proflstartmin)*60. 
#
#     proflendday=strmid(profendtime,0,2)
#     proflendmonth=strmid(profendtime,3,2)
#     proflendyear=strmid(profendtime,6,2)
#     proflendhr=strmid(profendtime,9,2)
#     proflendmin=strmid(profendtime,12,2)
#     profendsec=ndays*86400+ $
#           float(proflendday)*24*3600+ $
#           float(proflendhr)*3600.+float(proflendmin)*60. 
#
#;    if time lt 1 or time gt nstep then goto,j13
#
#     proftime_hr=fltarr(nstep)
#     proftime_min=fltarr(nstep)
#
#     for i=0,ic-1 do begin ; March 2006, note the change, 0,ic-1 for calc. of time
#
#      if i gt 1 then begin
#       print,i,sec(i),profstartsec,profendsec,tsprof
#       read,a
#       if sec(i) ge profstartsec and sec(i-1) lt profstartsec then begin
#         tsprof=i
#	 read,a
#       endif 
#       if sec(i) ge profendsec and sec(i-1) le profendsec then tsprof_end=i-1
#
#      end
#
#      proftime_hr(i)=(sec(i)-    $
#          (ndays*86400+float(lday(i))*24*3600+float(lmin(i))*60.))/3600. 
#      proftime_min(i)=(sec(i)-    $
#          (ndays*86400+float(lday(i))*24*3600+float(lhr(i))*3600.))/60. 
#
#     endfor
#
#   end
#
#
# end
#
#; definition of countour domain
#

# some general settings, to be updated for modified code
nlev=nlevel1
laxis=1
dfield=1

if int(choice) == 2:
   height_cont=0.
   while (float(height_cont) < height[1,0]): 
      print 'Contours up to which altitude ?(does not have to resemble alt. of levels)'
      print 'Altitude in [m] and at least higher than: ',height[0,0]
      for ii in range(nlev-1,0,-1):
         print 'level: %i alt.: %8.1f' % (ii,height[0,ii])

      height_cont=raw_input('Give altitude in [m]: ')
      type(height_cont)

   levh=-1
   while int(levh) < 1 or int(levh) > nlev: 
     for ii in range(0,nlev-1):
       if float(height_cont) < height[0,ii]:
          levh=int(ii)
	  break
  
#; fill contours, just contour lines or both
#
#  j111:
#  print,'Do you want filled contours (without contours) (0) ', $
#   'or just contour lines (1) ?'
#  read,cont
#  if cont lt 0 or cont gt 1 then goto,j111
#
# end

# Generate list with timestamps
dhr=(loctime[1]-loctime[0])/60.
base = datetime.datetime(year=int(lyear[ts1]),month=int(lmonth[ts1]),day=int(lday[ts1]),hour=int(lhr[ts1]))
date_list = [base + datetime.timedelta(hours=dhr*x) for x in range(istep)]

# line plots
if int(choice) == 0:

  xa =zeros((istep))
  param =zeros((nlevpr+1,istep))

  # assigning the parameter to be plotted
  for i in range (0,istep):
        ii=int(i+ts1)
        xa = date_list
	for j in range (0,int(nlevpr)):
            param[j,i] = conc[ii,j,jt]

# definition of the units for the different species
  unit = ['']*nspec
  unit_res = ['']*nspec
  minmax = [0]*nspec
  minval = [0]*nspec
  maxval = [0]*nspec

  unit[jt]=' ppbv] '
  unit_res=' [ppbv] '

  if (jt == ioh and plot_OHHO2_molec == 1) :
    unit_res=' [1e6 molec cm-3] '
  if (jt == iho2 and plot_OHHO2_molec == 1) :
    unit_res=' [1e8 molec cm-3] '
  if (jt == io1d) :
    unit_res=' [molec cm-3] '
  if (jt == irad) :
    unit_res=' [atoms cm-3] '

# ESS_lg_20121216+ added the radon concentration in m B==uerel m-3
  if (jt == irad) :
    unit_res=' [mBq m-3] '
# ESS_lg_20121216-

  if (jt == inoxinoy) :
    unit_res=' [-] '
  if (jt == ico2) :
    unit_res=' [ppmv] '
  if (jt == imvkmethacrat or jt == iisopprodrat or jt == iprodisoprat) :
    unit_res=' [-] '
  if (ichemtype == 5 and jt == iccn) :
    unit_res=' [# cm-3] '
# ESS_lg_20100728+ added the SOA aerosols
  if isoa == 1 :
    if (jt == iapinp1g or jt == iapinp1a or jt == iapinp2g or jt == iapinp2a) :
      unit_res=' [ug m-3]'

# calculating the minimum and maximum value for plotting, only if minmax
# is set to 1 in 1d_hcchem.par

  minmax[jt]=1
  if minmax[jt] == 1 :
    ymin=amin(param)
    ymax=amax(param)
  if minmax[jt] != 1 :    
    ymin=minval[jt]
    ymax=maxval[jt]

# plot title, x- and y-title

  if laxis == 1:
     titlx=' Time  '
  if laxis == 2:
     titlx=' Distance from ref. point [km]'

  titly=spec[jt]+' mixing ratio'+unit_res
  if (jt == ioh and plot_OHHO2_molec == 1) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == iho2 and plot_OHHO2_molec == 1) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == io1d) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == irad) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == inoxinoy) :
    titly=spec[jt]+' ratio'+unit_res
  if (jt == imvkmethacrat or jt == iisopprodrat or jt == iprodisoprat) :
    titly=spec[jt]+' ratio'+unit_res
  if (ichemtype == 5 and jt == iccn) :
    titly=spec[jt]+unit_res
# ESS_lg_20100728+ added the aerosols
  if (jt > ntrac_gas) :
    titly=spec[jt]+' conc. '+unit_res

  if dfield == 2 :
    if type == 1:
       titla[jt]='C2-C1 for '+spec[jt]
    if type == 2 :
       titla[jt]='100*(C2-C1)/C2 for '+spec[jt]
       unit[jt]=' [%]'
       unit_res=' [%]'
       ymax=min(100,ymax)
       ymin=max(-100.,ymin) 
    if type != 1 and type != 2 :
       titla[jt]=spec[jt]

# reseting the units for small values

  if (jt != ioh and plot_OHHO2_molec != 1) or (jt != iho2 and plot_OHHO2_molec != 1) or \
    jt != io1d or jt != irad or jt != inoxinoy or \
    jt != imvkmethacrat or jt != iisopprodrat or jt != iprodisoprat: 
    print 'resetting the units'
    if (ichemtype != 5 and jt != iccn) or ymax == 1.e-30:
      min_val=1
      ymult=1

      while abs(ymax-ymin) < 1 and abs(ymax-ymin) > 1e-3:
        param[:,:]=param[:,:]*1e3
        ymin=ymin*1e3
        ymax=ymax*1e3
        ymult=ymult*1e-3

        #unit_res='  ['+string(format='(1e6.0)',ymult*min_val)+unit[jt]

  iset_minmax=0

  print 'The minimum value is: %4.1f  %4s' % (ymin,unit_res)
  print 'The maximum value is: %4.1f  %4s' % (ymax,unit_res)

  print 'Type 1 to reset these values'
  iset_minmax= raw_input()
  type(iset_minmax)

  iset_minmax=int(iset_minmax)
  if iset_minmax == 1 :
    ymin= raw_input('Give the minimal value: ')
    type(ymin)   
    ymax= raw_input('Give the maximum value: ')
    type(ymax)

  print 'The plot title will be: ',titla[jt]+unit_res
  iadd_title=0
  add_title=''
  print 'To add something type 1 ?'
  iadd_title = raw_input()
  type(iadd_title)
  if int(iadd_title) == 1 :
    print 'Type additional title text '
    add_title = raw_input()
    type(add_title)
    plot_title=''# titla[jt]+unit_res+add_title

  titly=spec[jt]+' mixing ratio'+unit_res
  if dfield == 2:
     titly='d'+spec[jt]+' mixing ratio'+unit_res

  if (jt == ioh and plot_OHHO2_molec == 1) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == iho2 and plot_OHHO2_molec == 1) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == io1d) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == irad) :
    titly=spec[jt]+' conc. '+unit_res
  if (jt == inoxinoy) :
    titly=spec[jt]+' ratio'+unit_res
  if (jt == imvkmethacrat or jt == iisopprodrat or jt == iprodisoprat) :
    titly=spec[jt]+' ratio'+unit_res
  if (ichemtype == 5 and jt == iccn) :
    titly=spec[jt]+unit_res
# ESS_lg_20100728+ added the aerosols
  if (jt > ntrac_gas) :
    titly=spec[jt]+' conc. '+unit_res

  dy=ymax-ymin

# end assigning plotting parameters for choice == 0

  f,axs = subplots()

  # Options for plotting
  #axs.xaxis.set_major_locator(mdates.HourLocator(interval=2))
  #axs.xaxis.set_major_formatter(mdates.DateFormatter('%d%h'))

  axs.xaxis.set_major_locator(mdates.DayLocator(interval=1))
  axs.xaxis.set_major_formatter(mdates.DateFormatter('%d%h'))

  for j in range(0,nlevpr): 
     plot(xa,param[j,:])

  axs.set_xlabel('LT [hr]')
  axs.set_ylabel(spec[jt]+unit_res)
  #axs.set_title(spec[jt])

  show()

# contour plots
if int(choice) == 2:
  xa = date_list
  ya = height[0,0:levh]
  param =zeros((levh,istep,))

  # assigning the parameter to be plotted
  for i in range (0,istep):
      ii=int(i+ts1)
      for j in range (0,int(levh)):
          param[j,i] = conc[ii,j,jt]

  f,axs = subplots()

  # Options for plotting
  axs.xaxis.set_major_locator(mdates.DayLocator(interval=1))
  axs.xaxis.set_major_formatter(mdates.DateFormatter('%d%h'))

  vmax = numpy.nanmax(param)
  vmin = numpy.nanmin(param)
  zero_point = vmin/(vmin-vmax)
  levels=numpy.linspace(vmin,vmax,50)

  cmap = colors.LinearSegmentedColormap('my_colormap', makeZeroColormap(zero_point, cm.bwr))

  #cax = axs.pcolor(xa,ya,param)
  #cax = axs.contourf(xa,ya,param)

  #xa, ya  = meshgrid(xa, ya)
  #cax = axs.contourf(xa, ya, param, levels=levels, vmin=vmin, vmax=vmax, cmap=cm.bwr)
  #using the colormap Spectral
  #cax = axs.contourf(xa, ya, param, levels=levels, vmin=vmin, vmax=vmax, cmap=cm.Spectral)
  # reversed colorscale
  cax = axs.contourf(xa, ya, param, levels=levels, vmin=vmin, vmax=vmax, cmap=cm.Spectral_r)

  axs.set_xlabel('Time')
  axs.set_ylabel('Height [m]')
  axs.set_title(spec[jt]+unit_res)
  cbar = f.colorbar(cax)
  show()

# introducing the option to write plotted data to output file for
# processing in other graphical programs

# writing the unmanipulated dataset

if int(tintval) == 0 or int(tintval) == 1:

       # open file for writing:
       f = open('inputmod.out','w')
       print 'Writing data to output file "inputmod.out"'
       f.write("{:14s}".format("Time")) 
       for i in range(levh1-1,levh2):
          f.write( "{:10.3f}".format(height[0,i])) 
       f.write("\n")

       for it in range(ts1,ts2):
          iit = int(it-ts1)
          #        writing to output file
          f.write("{:14s}".format(ldateltime[it])) 
          for i in range(levh1-1,levh2):
             f.write(""+"{:10.3f}".format(param[i,iit]))
          f.write("\n")

       f.close()

##     writing the average diurnal cycle to file
#if int(tintval) == 2:
#
#      ts1=0
#      ts2=no_intervals-1
#
#; ESS_lg_20111231+ added the median
#
#      y = fltarr(no_intervals,3,npar1)
#      xx = fltarr(no_intervals)
#
#      names2=strarr(3,npar1)
#      names2(0,0:npar1-1)=height(levh1-1:levh2-1,ts1)
#      names2(1,0:npar1-1)='median.'
#      names2(2,0:npar1-1)='std.'
#
#      openw,1,'inputmod_diurnavg.out'
#      print,'Writing data to output file "inputmod_diurnavg.out"'
#      printf,1,format='(80a15)','Time',names2(0:2,0:npar1-1)
#
#      for it = ts1,ts2 do begin
#
#        y(it,0,0:npar1-1) = avgparam(it,0:npar1-1)
#	y(it,1,0:npar1-1) = medianparam(it,0:npar1-1)
#        y(it,2,0:npar1-1) = stdparam(it,0:npar1-1)
#
#;       since the model output is written away such that the averaged values
#;       of a specific interval reflect the average of the interval ending with
#;       the assigned time in the output files which is the number of printing
#;       intervals + half the timestep 
#
#        xx(it) = etime(it)
#
#;       writing to output file
#
#       printf,1,format='(f15.2,80(5x,e10.3))',xx(it),y(it,0:2,0:npar1-1)
#
#; ESS_lg_20111231-
#
#      endfor
# 
#     close,1
#    end

#;   end writing data
