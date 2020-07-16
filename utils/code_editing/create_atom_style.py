from __future__ import print_function, absolute_import, division
import numpy as np
import sys

#Dictionary containing per-atom property names to be added
# Keys are names of properties
# Values are: C++ data type, dimensionality, forward comm, reverse comm
#  Set forward comm True if quantity is being used by LAMMPS and changes at every step, False otherwise
#  Set reverse comm True if quantity is being modified by LAMMPS at every step, False otherwise
#  We expect going forward, DEMSI will not need addtional per-atom properties that are forward/reverse communicated,
#   meaning most entries will just have 'False, False' settings
new_props = {"forcing": ("double", 2, False, False), 
             "mean_thickness": ("double", 1, False, False), 
             "min_thickness": ("double" , 1, False, False),
             "ice_area": ("double", 1, False, False),
             "coriolis": ("double", 1, False, False),
             "ocean_vel": ("double", 2, False, False),
             "bvector": ("double", 2, False, False)
             }

#Location of DEMSI/LAMMPS src directory:
src_location = "../../LAMMPS/src/"

#New atom style name; the files that are generated will be named atom_vec_<name>.cpp/h (e.g. atom_vec_demsi.cpp)
name = "demsi"
name = name.lower()

#Convert entries to objects for ease of use interally in this script
class Property:
    pass
props = []
for k,v in new_props.items():
    prop = Property()
    prop.name = k
    prop.dtype = v[0]
    prop.ndim = v[1]
    prop.forward_comm = v[2]
    prop.reverse_comm = v[3]    
    props.append(prop)

#Basic scheme:
#For forward comm:
# buf[m++] = x[j][2] not followed by ubuf(tag[j]).d --> only add properties if forward comm
    # if followed by ubuf, add all properties
# x[i][2] = buf[m++] not followed by tag[i] = --> only add properties if forward comm
    # if followed by tag[i], add all properties
    
##Exception is pack_exchange and unpack_exchange, in which case all properties are added
# For pack_exchange, indexing is 'i'; for unpack_exchange, it's 'nlocal'
    
#For reverse comm (unlikely to ever need this):
# buf[m++] = torque[i][2] --> only add properties if reverse comm    
# torque[j][2[] += buf[m++] --> only add properties if forward comm


bondstring = "int **nspecial;\ntagint **special;\nint *num_bond;\nint **bond_type;\ntagint **bond_atom\n;\n"

try:
    hfile = open(src_location+"/atom_vec_sphere.h", "r")
except:
    print("Could not open ",src_location+hfile)
    sys.exit()

try:
    cppfile = open(src_location+"/atom_vec_sphere.cpp", "r")
except:
    print("Could not open ",src_location+cppfile)
    sys.exit()

try:
    hfileout = open("atom_vec_"+name+".h", "w")
    cppfileout = open("atom_vec_"+name+".cpp", "w")
except:
    print("Could not open output files")
    sys.exit()

for line in hfile:
    if "AtomStyle" in line:
        hfileout.write("AtomStyle("+name+", AtomVec"+name.capitalize()+")\n")
    elif "#ifndef" in line:
        hfileout.write("#ifndef LMP_ATOM_VEC_"+name.upper()+"_H\n")
    elif "#define" in line:
        hfileout.write("#define LMP_ATOM_VEC_"+name.upper()+"_H\n")
    elif "AtomVecSphere" in line:
        hfileout.write(line.replace("AtomVecSphere", "AtomVecDemsi"))
    elif "**omega,**torque" in line:
        hfileout.write(line)
        for p in props:
            hfileout.write(p.dtype+" "+"*"*p.ndim+" "+p.name+";\n")
    elif "int radvary" in line:
        hfileout.write(line)
        hfileout.write(bondstring)
    else:
        hfileout.write(line)
hfileout.close()
hfile.close()

for line in cppfile:
    if "atom_vec_sphere.h" in line:
        cppfileout.write(line.replace("atom_vec_sphere", "atom_vec_"+name))
    elif "cmath" in line:
        cppfileout.write(line.replace("cmath", "math.h"))
    elif "cstdlib" in line:
        cppfileout.write(line.replace("cstdlib", "stdlib.h"))
    elif "cstring" in line:
        cppfileout.write(line.replace("cstring", "string.h"))
    elif "AtomVecSphere" in line:
        cppfileout.write(line.replace("AtomVecSphere", "AtomVec"+name.capitalize()))
    elif "molecular" in line:
        cppfileout.write("molecular = 1;\n")
    elif "size_forward" in line:
        words = line.split()        
        cppfileout.write("size_forward = "+str(int(3+np.sum([p.forward_comm for p in props])))+";\n")
    elif "size_border" in line:
        cppfileout.write("size_border = "+str(int(8+np.sum([p.ndim for p in props])))+";\n")
    elif "radvary = 0" in line:
        cppfileout.write(line)
        cppfileout.write("comm_x_only = 0;\n")
        cppfile.readline()
    elif "size_reverse" in line:
        cppfileout.write("size_reverse = "+str(int(6+np.sum([p.reverse_comm for p in props])))+";\n")
    elif "atom->sphere_flag" in line:
        cppfileout.write(line)
        cppfileout.write("atom->demsi_flag = 1;\n")
    elif "atom->torque_flag" in line:
        cppfileout.write(line)
        cppfileout.write("bonds_allow = 1;\n")
    elif "torque = memory->grow" in line:
        cppfileout.write(line+"\n")
        for p in props:
            if p.ndim == 1:
                numstr = ""
            else:
                numstr = str(p.ndim)+","
            cppfileout.write("  "+p.name+" = memory->grow(atom->"+p.name+",nmax,"+numstr+"\"atom:"+p.name+"\");\n")
        
        cppfileout.write("\n  nspecial = memory->grow(atom->nspecial,nmax,3,\"atom:nspecial\");\n")
        cppfileout.write("  special = memory->grow(atom->special,nmax,atom->maxspecial,\"atom:special\");\n")
        cppfileout.write("  num_bond = memory->grow(atom->num_bond,nmax,\"atom:num_bond\");\n")
        cppfileout.write("  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,\"atom:bond_type\");\n")
        cppfileout.write("  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,\"atom:bond_atom\");\n")

    elif "omega = atom->omega;" in line:
        cppfileout.write(line+"\n")
        for p in props:
            cppfileout.write("  "+p.name+" = atom->"+p.name+";\n")
        cppfileout.write("\n  nspecial = atom->nspecial; special = atom->special;\n")
        cppfileout.write("  num_bond = atom->num_bond; bond_type = atom->bond_type;\n")
        cppfileout.write("  bond_atom = atom->bond_atom;\n")
        
    elif "omega[j][2] = omega[i][2];" in line:
        cppfileout.write(line+"\n")
        for p in props:
            if p.ndim > 1:
                for index in range(0,p.ndim):
                    cppfileout.write("  "+p.name+"[j]["+str(index)+"] = "+p.name+"[i]["+str(index)+"];\n");
            else:
                cppfileout.write("  "+p.name+"[j] = "+p.name+"[i];\n")
        cppfileout.write("  nspecial[j][0] = nspecial[i][0];\n")
        cppfileout.write("  nspecial[j][1] = nspecial[i][1];\n")
        cppfileout.write("  nspecial[j][2] = nspecial[i][2];\n")
        cppfileout.write("  for (int k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];\n")
        cppfileout.write("    num_bond[j] = num_bond[i];\n")
        cppfileout.write("  for (int k = 0; k < num_bond[j]; k++) {\n")
        cppfileout.write("    bond_type[j][k] = bond_type[i][k];\n")
        cppfileout.write("    bond_atom[j][k] = bond_atom[i][k];\n")
        cppfileout.write("}\n")
    
    #Forward packing
    elif "buf[m++] = x[j][2]" in line or "buf[m++] = x[i][2]" in line:
        cppfileout.write(line)
        nextline = cppfile.readline()        
        for p in props:
            if p.forward_comm or "ubuf(tag" in nextline or "buf[m++] = x[i][2]" in line:
                if "buf[m++] = x[i][2]" in line:
                    index0 = "i"
                else:
                    index0 = "j"
                if p.ndim > 1:
                    for index1 in range(0, p.ndim):
                        cppfileout.write("  buf[m++] = "+p.name+"["+index0+"]["+str(index1)+"];\n")
                else:
                    cppfileout.write("  buf[m++] = "+p.name+"["+index0+"];\n")
        cppfileout.write(nextline)
         
    #Forward unpacking
    elif "x[i][2] = buf[m++]" in line or "x[nlocal][2] = buf[m++]" in line:
        cppfileout.write(line)
        nextline = cppfile.readline()        
        for p in props:
            if p.forward_comm or "ubuf" in nextline or "x[nlocal][2] = buf[m++]" in line:
                if "x[nlocal][2] = buf[m++]" in line:
                    index0 = "nlocal"
                else:
                    index0 = "i"
                if p.ndim > 1:
                    for index1 in range(0, p.ndim):
                        cppfileout.write(p.name+"["+index0+"]["+str(index1)+"] = buf[m++];\n")
                else:
                    cppfileout.write(p.name+"["+index0+"] = buf[m++];\n")
        cppfileout.write(nextline)
    
    #Reverse packing
    elif "buf[m++] = torque[i][2];" in line:
        cppfileout.write(line)
        for p in props:
            if p.reverse_comm:
                if p.ndim > 1:
                    for index in range(0, p.ndim):
                        cppfileout.write("buf[m++] = "+p.name+"[j]["+str(index)+"];\n")
                else:
                    cppfileout.write("buf[m++] = "+p.name+"[j]\n;")          
    
    
    #Reverse unpacking    
    elif "torque[j][2] += buf[m++];" in line:
        cppfileout.write(line)
        for p in props:
            if p.reverse_comm:
                if p.ndim > 1:
                    for index in range(0, p.ndim):
                        cppfileout.write(p.name+"[j]["+str(index)+"] += buf[m++];\n")
                else:
                    cppfileout.write(p.name+"[j] += buf[m++]\n;")
                    
    #Bond info in pack_exchange
    elif "buf[m++] = omega[i][2];" in line:
        cppfileout.write(line+"\n")
        cppfileout.write("\n  buf[m++] = ubuf(num_bond[i]).d;\n")
        cppfileout.write("  for (int k = 0; k < num_bond[i]; k++) {\n")
        cppfileout.write("    buf[m++] = ubuf(bond_type[i][k]).d;\n")
        cppfileout.write("    buf[m++] = ubuf(bond_atom[i][k]).d;\n")
        cppfileout.write("  }\n")
        cppfileout.write("  buf[m++] = ubuf(nspecial[i][0]).d;\n")
        cppfileout.write("  buf[m++] = ubuf(nspecial[i][1]).d;\n")
        cppfileout.write("  buf[m++] = ubuf(nspecial[i][2]).d;\n")
        cppfileout.write("  for (int k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;\n")
    
    #Bond info in unpack_exchange
    elif "omega[nlocal][2] = buf[m++]" in line:
        cppfileout.write(line+"\n")
        cppfileout.write("num_bond[nlocal] = (int) ubuf(buf[m++]).i;\n")
        cppfileout.write("for (int k = 0; k < num_bond[nlocal]; k++) {\n")
        cppfileout.write("  bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;\n")
        cppfileout.write("  bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;\n")
        cppfileout.write("}\n")    
        cppfileout.write("nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;\n")
        cppfileout.write("nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;\n")
        cppfileout.write("nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;\n")
        cppfileout.write("for (int k = 0; k < nspecial[nlocal][2]; k++)\n")
        cppfileout.write("  special[nlocal][k] = (tagint) ubuf(buf[m++]).i;\n")
        
    elif "4.0*MY_PI/3.0" in line:
        if "rmass[nlocal]" in line:
            cppfileout.write("rmass[nlocal] = MY_PI*radius[nlocal]*radius[nlocal];\n")
            cppfile.readline()
        elif "buf" in line:
            cppfileout.write(line.replace("4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i]", "MY_PI*radius[i]*radius[i]"))
        
    elif "omega[nlocal][2] = 0.0" in line:
        cppfileout.write(line+"\n")
        for p in props:
            if "int" in p.dtype:
                str_zero = "0"
            elif "double" in p.dtype or "float" in p.dtype:
                str_zero = "0.0"
            if p.ndim > 1:
                for index in range(0, p.ndim):
                    cppfileout.write("  "+p.name+"[nlocal]["+str(index)+"] = "+str_zero+";\n")
            else:
                cppfileout.write("  "+p.name+"[nlocal] = "+str_zero+";\n")
        cppfileout.write("\nnum_bond[nlocal] = 0;\n")
        cppfileout.write("nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;\n")
    
    elif "bytes += memory->usage(torque" in line:
        cppfileout.write(line+"\n")
        for p in props:
            if p.ndim > 1:
                dim_str = ","+str(p.ndim)
            else:
                dim_str = ""
            cppfileout.write("  if (atom->memcheck(\""+p.name+"\")) bytes += memory->usage("+p.name+",nmax"+dim_str+");\n")
        cppfileout.write("  if (atom->memcheck(\"num_bond\")) bytes += memory->usage(num_bond,nmax);\n")
        cppfileout.write("  if (atom->memcheck(\"bond_type\"))\n")
        cppfileout.write("    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);\n")
        cppfileout.write("  if (atom->memcheck(\"bond_atom\"))\n")
        cppfileout.write("    bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);\n")
        
        cppfileout.write("  if (atom->memcheck(\"nspecial\")) bytes += memory->usage(nspecial,nmax,3);\n")
        cppfileout.write("  if (atom->memcheck(\"special\"))\n")
        cppfileout.write("    bytes += memory->usage(special,nmax,atom->maxspecial);\n")

    else:
        cppfileout.write(line)
        
cppfile.close()
cppfileout.close()

#Also process atom.cpp/atom.h
try:
    atomhin = open(src_location+"/atom.h", "r")
except:
    print("Could not open ",src_location+atomhin)
    sys.exit()

try:
    atomcin = open(src_location+"/atom.cpp", "r")
except:
    print("Could not open ",src_location+cppfile)
    sys.exit()

try:
    atomhout = open("atom.h", "w")
except:
    print("Could not open atom.h for writing\n")
    sys.exit()

try:
    atomcout = open("atom.cpp", "w")
except:
    print("Could not open atom.cpp for writing\n")
    sys.exit()
    

for line in atomhin:
    #Discard previous DEMSI edits
    if "DEMSI" in line:
        templine = atomhin.readline()    
        while templine.strip():
            templine = atomhin.readline()
    elif "int damage_flag" in line:
        atomhout.write(line+"\n")
        atomhout.write("  //USER-DEMSI package\n  int demsi_flag;\n\n")
    elif "int cc_species" in line:
        atomhout.write(line+"\n")
        atomhout.write("  //USER-DEMSI package\n")
        for p in props:
            atomhout.write("  "+p.dtype+" "+"*"*p.ndim+p.name+";\n")
        atomhout.write("\n")
    else:
         atomhout.write(line)

atomhout.close()
atomhin.close()

for line in atomcin:
    #Discard previous DEMSI edits
    if "DEMSI" in line:
        templine = atomcin.readline()    
        while templine.strip():
            templine = atomcin.readline()
    elif "edpd_temp = edpd_flux" in line:
        atomcout.write(line+"\n")
        atomcout.write("  //USER-DEMSI\n")
        for p in props:
            atomcout.write("  "+p.name+" = NULL;\n")
    else:
        atomcout.write(line)

atomcout.close()
atomcin.close()
    
            
