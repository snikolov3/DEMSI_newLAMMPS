import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt

filePrefix = "particles_out"

oceanVariablesToAdd = []
#oceanVariablesToAdd.append("")
oceanVariablesToAdd.append("atmosReferenceTemperature2mOcean")
oceanVariablesToAdd.append("atmosReferenceHumidity2mOcean")
oceanVariablesToAdd.append("albedoVisibleDirectOcean")
oceanVariablesToAdd.append("albedoVisibleDiffuseOcean")
oceanVariablesToAdd.append("albedoIRDirectOcean")
oceanVariablesToAdd.append("albedoIRDiffuseOcean")
oceanVariablesToAdd.append("longwaveUpOcean")
oceanVariablesToAdd.append("sensibleHeatFluxOcean")
oceanVariablesToAdd.append("latentHeatFluxOcean")
oceanVariablesToAdd.append("evaporativeWaterFluxOcean")
oceanVariablesToAdd.append("seaSurfaceTiltU")
oceanVariablesToAdd.append("seaSurfaceTiltV")
oceanVariablesToAdd.append("freezingMeltingPotential")
oceanVariablesToAdd.append("airOceanDragCoefficientRatio")

filein = open("varnames.txt","r")
varnames = filein.readlines()
filein.close()

#filenameRemap = "./output_remap/merged.nc"
#filenameRemapOcean = "./output_remap/merged_ocean.nc"
#filenameNoRemap = "./output_noremap/merged.nc"

filenameRemap = "./output_remap/merged.nc"
filenameRemapOcean = "./output_remap/merged_ocean.nc"
filenameNoRemap = "./output_noremap/merged.nc"

fileRemap = Dataset(filenameRemap,"r")
fileRemapOcean = Dataset(filenameRemapOcean,"r")
fileNoRemap = Dataset(filenameNoRemap,"r")

nRecords = len(fileRemap.dimensions["record"])

for varname in varnames:

    varname = varname.strip()
    print(varname)

    if (varname in fileRemap.variables):

        if (len(fileRemap.variables[varname].dimensions) == 2):

            VariableRemap   = fileRemap.variables[varname][:,0]
            VariableNoRemap = fileNoRemap.variables[varname][:,0]
    
        elif (len(fileRemap.variables[varname].dimensions) == 3):
                
            VariableRemap   = fileRemap.variables[varname][:,0,0]
            VariableNoRemap = fileNoRemap.variables[varname][:,0,0]

        if (varname in oceanVariablesToAdd):

            if (len(fileRemap.variables[varname].dimensions) == 2):
                VariableRemap = fileRemapOcean.variables[varname][:,0]
    
            elif (len(fileRemap.variables[varname].dimensions) == 3):
                VariableRemap = fileRemapOcean.variables[varname][:,0,0]

        filenameOut = "./figures/compare_%s.png" %(varname)

        plt.plot(VariableRemap,label="remap",linewidth=0.7)
        plt.plot(VariableNoRemap,label="noremap",linewidth=0.5)
        #plt.plot(VariableRemap[2350:2400],label="remap",linewidth=0.7)
        #plt.plot(VariableNoRemap[2350:2400],label="noremap",linewidth=0.5)
        #plt.plot(VariableNoRemap[0:50]-VariableRemap[0:50],linewidth=0.5)
        #plt.plot(VariableNoRemap[:]-VariableRemap[:],linewidth=0.5)
        plt.legend()
        #plt.xlim((0,50))
        plt.savefig(filenameOut,dpi=300)
        plt.cla()
        plt.close()

fileRemap.close()
fileRemapOcean.close()
fileNoRemap.close()

