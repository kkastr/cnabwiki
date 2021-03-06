import numpy as np 
import sys 
import subprocess as sp
import os
import glob

#TODO add option to replace files with backups and re-do update

def insert_modifications(filename,identifier_string_list,list_string_to_insert):
    lines_from_read_file = open(filename,"r").readlines()
    
    
    isInserted = []

    for s in list_string_to_insert:
        line_already_exists = 0
        for line in lines_from_read_file:
            if line.startswith(s):
                line_already_exists = 1
                # print(f"|{s}| has already been added to {filename}\n")

        isInserted.append(line_already_exists)
    if np.sum(isInserted)<len(isInserted):
        file_to_write = open(filename,"w")
        for line in lines_from_read_file:
            file_to_write.write(line)
            for s in list_string_to_insert:
                idx = list_string_to_insert.index(s)
                if isInserted[idx] ==0:
                    if line.startswith(identifier_string_list[idx]):
                        new_line = list_string_to_insert[idx] +'\n'
                        file_to_write.write(new_line)
                    else:
                        pass
        file_to_write.close()


def make_backups(hoomd_path,current_dir,files_to_modify):
    bak_flag = 0
    create_bak_dir = 'mkdir -p backups'
    sp.run(create_bak_dir.split())
    files_in_dir = [os.path.basename(idx) for idx in glob.glob('backups/*')]
    
    for fidx in files_to_modify:
        if not fidx in files_in_dir:
            command = f'cp {hoomd_path}{fidx} {current_dir}/backups/'
            sp.run(command.split())
        else:
            bak_flag = 1
    
    if bak_flag==1:
        print('Backup of original files already exists.')
    else:
        print('Created backup of original files.')

    


if len(sys.argv) >1:
    path_to_hoomd = sys.argv[1]
else:
    print("No path to hoomd was given!\n")
    print("Please provide the full path to hoomd-blue/hoomd/md/.\n")
    print("Example path: /home/username/hoomd-blue/hoomd/md/")
    sys.exit()

cwd = sp.check_output(['pwd'],universal_newlines=True).rstrip()

print(f'Currently in {cwd}')

files_to_modify = ['AllExternalPotentials.h','PotentialExternalGPU.cu','module-md.cc','external.py']

make_backups(path_to_hoomd,cwd,files_to_modify)


evaluator_names =['NanoporeEfield','NanoporeGeometry','ThickNanoporeEfield','NanoporeFriction']
shortnames = ['np_efield','nanopore','thick_np_efield','pore_friction']

for en in evaluator_names:

    copy_evaluators = f'cp {cwd}/hoomd-evaluators/Evaluator{en}.h {path_to_hoomd}'
    os.system(copy_evaluators)

os.chdir(path_to_hoomd)

print('Should be inside the hoomd/md directory\n')

#make backups

# modify AllExternalPotentials.h

print('Modifying "AllExternalPotentials.h"...\n')

for en in evaluator_names:
    fn ='AllExternalPotentials.h'
    ids = ['#include "PotentialExternal.h"','typedef PotentialExternal<EvaluatorExternalElectricField> PotentialExternalElectricField;','typedef PotentialExternalGPU<EvaluatorExternalElectricField> PotentialExternalElectricFieldGPU;']

    stin = [f'#include "Evaluator{en}.h"',f'typedef PotentialExternal<Evaluator{en}> Potential{en};',f'typedef PotentialExternalGPU<Evaluator{en}> Potential{en}GPU;']



    insert_modifications(fn,ids,stin)

print('Modifying "AllExternalPotentials.h"... Done.\n')


#modify PotentialExternalGPU.cu

print('Modifying "PotentialExternalGPU.cu"...\n')

for en in evaluator_names:
    fn ='PotentialExternalGPU.cu'
    ids = ['#include "EvaluatorExternalElectricField.h"','template cudaError_t gpu_cpef<EvaluatorExternalElectricField>(const external_potential_args_t& external_potential_args, const typename EvaluatorExternalElectricField::param_type *d_params, const typename EvaluatorExternalElectricField::field_type *d_field);']

    stin = [f'#include "Evaluator{en}.h"',f'template cudaError_t gpu_cpef<Evaluator{en}>(const external_potential_args_t& external_potential_args, const typename Evaluator{en}::param_type *d_params,const typename Evaluator{en}::field_type *d_field);']


    insert_modifications(fn,ids,stin)

print('Modifying "PotentialExternalGPU.cu"... Done.\n')

#modify module-md.cc
print('Modifying "module-md.cc"...\n')
for en in evaluator_names:

    fn = 'module-md.cc'
    ids = ['    export_PotentialExternal<PotentialExternalElectricField>(m, "PotentialExternalElectricField");','    export_PotentialExternalGPU<PotentialExternalElectricFieldGPU, PotentialExternalElectricField>(m, "PotentialExternalElectricFieldGPU");']
    stin = [f'    export_PotentialExternal<Potential{en}>(m, "Potential{en}");',f'    export_PotentialExternalGPU<Potential{en}GPU, Potential{en}>(m, "Potential{en}GPU");']

   
    insert_modifications(fn,ids,stin)

print('Modifying "module-md.cc"... Done.\n')

#modify external.py


#electric field


en = evaluator_names[0]
fn = 'external.py'
idx = evaluator_names.index(en)
class_name = shortnames[idx]
string_block = f"""
class {class_name}(_external_force):

    def __init__(self, field, name=""):
        hoomd.util.print_status_line();

        # initialize the base class
        _external_force.__init__(self, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.Potential{en}(hoomd.context.current.system_definition,self.name);
        else:
            self.cpp_force = _md.Potential{en}GPU(hoomd.context.current.system_definition,self.name);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['V0', 'a', 'c']

        self.field_coeff = tuple(field)

    def process_coeff(self, coeff):
        V0 = coeff['V0']
        a = coeff['a']
        c = coeff['c']

        return _hoomd.make_scalar3(V0, a, c)
    def process_field_coeff(self, field):
        return _hoomd.make_scalar3(field[0],field[1],field[2])
"""

print("Appending the python class for the electric field to external.py...")

lines_from_read_file = open(fn ,"r").readlines()

class_already_exists = 0
for line in lines_from_read_file:
    if line.startswith(f'class {class_name}(_external_force)'):
        class_already_exists = 1
    else:
        pass
if class_already_exists==0:
    file_to_append = open(fn,"a")
    file_to_append.write(string_block)
else:
    pass
print("Appending the python class for the electric field to external.py... Done.")

# nanopore geometry

en = evaluator_names[1]
fn = 'external.py'
idx = evaluator_names.index(en)
class_name = shortnames[idx]
string_block = f"""
class {class_name}(_external_force):

    def __init__(self, field, name=""):
        hoomd.util.print_status_line();

        # initialize the base class
        _external_force.__init__(self, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.Potential{en}(hoomd.context.current.system_definition,self.name);
        else:
            self.cpp_force = _md.Potential{en}GPU(hoomd.context.current.system_definition,self.name);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['epsilon', 'sigma', 'rpore', 'tpore']

        self.field_coeff = tuple(field)

    def process_coeff(self, coeff):
        epsilon = coeff['epsilon']
        sigma = coeff['sigma']
        rpore = coeff['rpore']
        tpore = coeff['tpore']
        

        # can't pass these anymore as I need sigma to be variable for corrugated chains
        # lj1 = 4.0 * epsilon * math.pow(sigma, 12.0)
        # lj2 = 4.0 * epsilon * math.pow(sigma, 6.0)
        return _hoomd.make_scalar4(epsilon, sigma, rpore, tpore)

    #the field parameters are the center of the nanopore in x,y,z
    def process_field_coeff(self, field):
        return _hoomd.make_scalar3(field[0],field[1],field[2])
"""

print("Appending the python class for the nanopore to external.py...")

lines_from_read_file = open(fn ,"r").readlines()

class_already_exists = 0
for line in lines_from_read_file:
    if line.startswith(f'class {class_name}(_external_force)'):
        class_already_exists = 1
    else:
        pass
if class_already_exists==0:
    file_to_append = open(fn,"a")
    file_to_append.write(string_block)
else:
    pass
print("Appending the python class for the nanopore to external.py... Done.")


# nanopore electric field for pores with non vanishing thickness

en = evaluator_names[2]
fn = 'external.py'
idx = evaluator_names.index(en)
class_name = shortnames[idx]
string_block = f"""
class {class_name}(_external_force):

    def __init__(self, field, name=""):
        hoomd.util.print_status_line();

        # initialize the base class
        _external_force.__init__(self, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.Potential{en}(hoomd.context.current.system_definition,self.name);
        else:
            self.cpp_force = _md.Potential{en}GPU(hoomd.context.current.system_definition,self.name);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['V0', 'a', 'c', 'tpore']

        self.field_coeff = tuple(field)

    def process_coeff(self, coeff):
        V0 = coeff['V0']
        a = coeff['a']
        c = coeff['c']
        tpore = coeff['tpore']

        return _hoomd.make_scalar4(V0, a, c, tpore)

    #the field parameters are the center of the nanopore in x,y,z
    def process_field_coeff(self, field):
        return _hoomd.make_scalar3(field[0],field[1],field[2])
"""

print("Appending the python class for the thick nanopore electric field to external.py...")

lines_from_read_file = open(fn ,"r").readlines()

class_already_exists = 0
for line in lines_from_read_file:
    if line.startswith(f'class {class_name}(_external_force)'):
        class_already_exists = 1
    else:
        pass
if class_already_exists==0:
    file_to_append = open(fn,"a")
    file_to_append.write(string_block)
else:
    pass
print("Appending the python class for the thick nanopore electric field to external.py... Done.")



en = evaluator_names[3]
fn = 'external.py'
idx = evaluator_names.index(en)
class_name = shortnames[idx]
string_block = f"""
class {class_name}(_external_force):

    def __init__(self, field, name=""):
        hoomd.util.print_status_line();

        # initialize the base class
        _external_force.__init__(self, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.Potential{en}(hoomd.context.current.system_definition,self.name);
        else:
            self.cpp_force = _md.Potential{en}GPU(hoomd.context.current.system_definition,self.name);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['epsilon', 'alpha']

        self.field_coeff = tuple(field)

    def process_coeff(self, coeff):
        epsilon = coeff['epsilon']
        alpha = coeff['alpha']

        return _hoomd.make_scalar2(epsilon, alpha)

    #the field parameters are the center of the nanopore in x,y,z
    def process_field_coeff(self, field):
        return _hoomd.make_scalar3(field[0],field[1],field[2])
"""

print("Appending the python class for pore friction to external.py...")

lines_from_read_file = open(fn ,"r").readlines()

class_already_exists = 0
for line in lines_from_read_file:
    if line.startswith(f'class {class_name}(_external_force)'):
        class_already_exists = 1
    else:
        pass
if class_already_exists==0:
    file_to_append = open(fn,"a")
    file_to_append.write(string_block)
else:
    pass
print("Appending the python class for pore friction to external.py... Done.")
