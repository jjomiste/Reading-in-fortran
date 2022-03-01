import h5py
import numpy as np
import sys

radios=['1.', '1.05', '1.1','1.15', '1.2','1.25', '1.3','1.35', '1.4','1.45', '1.5','1.55', '1.6','1.65', '1.7', '1.75', '1.8','1.85', '1.9']

if len(sys.argv[:]) != 3:
    print("Tienes que incluir 2 argumentos: la ruta y la raiz del archivo")
    print("Deteniendo el código...")
    exit()

ruta =sys.argv[1]
raiz=sys.argv[2]
pes_salida='PES_'+raiz+'.dat'

f=open(pes_salida, 'w')

for radio in radios:
    archivo_in=ruta+'/'+raiz
    archivo_in=archivo_in+str(radio)+'.h5'
    f.write(radio)
    #Now we read the hdf5 file
    fras=h5py.File(archivo_in,'r')
    for energy in fras['ROOT_ENERGIES']:
        f.write('          ')
        f.write(str(energy))
    
    f.write('\n')


