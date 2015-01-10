findIonsat.py: read all the files and plot Id-Vg for Vd=0.86 V, in linear and log scale
findIonlin.py: read all the files and plot Id-Vg for Vd=0.05 V, in linear and log scale
findIonIoff.py:plot Ion vs Ioff plot
extractresults.py:extract results script from I-V results
plotelectricalextraction.py: plot data of simulations such as Vth, SS, Ion vs Ioff etc.
extractresultscorner.py: extract results script from I-V results of corner simulations, additional code had to be added since some name of files were wrongly written (missing Nfin parameter for 1e18 body doping)
plotelectricalextractioncorner.py: plot data of simulations such as Vth, SS, Ion vs Ioff etc.
extracDIBL.py: extract DIBL and write down all electrical measurements of each device in a file (sat and lin)
plotelectricalextractioncornersatlin.py: uses results from extracDIBL.py to plot different electrical quantities
