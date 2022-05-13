# SCM-MOSAiC

In this directory you will find a .zip file that contains the model code for the Single-Column atmospheric chemistry and meteorological Model (SCM) set up for 1-year simulations for the Multidisciplinary drifting Observatory for the Study of Arctic Climate (MOSAiC) campaign.
Feel free to use the model code and output. If you have any questions, please contact Sjoerd Barten (sjoerd.barten@wur.nl) or Laurens Ganzeveld (laurens.ganzeveld@wur.nl).

The structure of the main directories is as follows:
- 'case/': This is where the study cases are set up and run.
- 'case/input/': This is where the main input namelists are defined.
- 'case/output/': This is where you can find the model output files. (This is the place to go if you want to use the model output)
- 'work/': This is where all the routines and parameterizations of the model are located.

Note that not all simulated variables are by default in the model output. This would cause extremely large data streams and slowing down of model code. If additional parameters are requested we can run these model simulation again and include these parameters.
Furthermore, the SCM will not run unless input files from ECMWF ERA5 and CAMS reanalysis are provided. Due to the very large size of these files and limited storage space these are not included in the directories. These files can be downloaded through the Climate Data Store (ERA5) and Atmosphere Data Store (CAMS) or requested via email.

Sjoerd Barten

sjoerd.barten@wur.nl; Wageningen University; Meteorology & Air Quality Section
