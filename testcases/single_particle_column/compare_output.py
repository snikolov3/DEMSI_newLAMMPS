from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
import string

varname = "iceVolumeCategory"

filenames = sorted(glob.glob("./output/*"))

minval = 1e30
maxval = -1e30

for filenameDEMSI in filenames:

    year  = int(filenameDEMSI[23:27])
    month = int(filenameDEMSI[28:30])
    day   = int(filenameDEMSI[31:33])
    hour  = int(filenameDEMSI[34:36])

    fileIn = Dataset(filenameDEMSI,"r")

    try:
        varDEMSI = fileIn.variables[varname][0,0,0]
    except:
        try:
            varDEMSI = fileIn.variables[varname][0,0]
        except:
            try:
                varDEMSI = fileIn.variables[varname][0]
            except:
                print(filenameDEMSI)

    fileIn.close()

    filenameMPAS = "/Users/akt/Work/DEMSI/Icepack/MPAS/rundir_sc/output/output.2000-%2.2i-%2.2i_%2.2i.00.00.nc" %(month,day,hour)

    try:
        fileIn = Dataset(filenameMPAS,"r")

        try:
            varMPAS = fileIn.variables[varname][0,0,0,0]
        except:
            try:
                varMPAS = fileIn.variables[varname][0,0,0]
            except:
                varMPAS = fileIn.variables[varname][0,0]

        fileIn.close()

        print(filenameDEMSI, varDEMSI-varMPAS, varDEMSI, varMPAS)
        minval = min(minval,varDEMSI-varMPAS)
        maxval = max(maxval,varDEMSI-varMPAS)

    except:
        pass

print(varname, minval, maxval)
