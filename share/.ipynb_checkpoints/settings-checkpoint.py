import os
# Root directory of the QGPU modules
ROOT = os.path.dirname(os.path.realpath(__file__))
ROOT = ROOT.split('src')[0]

SLURM =  {
'KEBNE' : {'--nodes'            : '1',
           '--ntasks-per-node'  : '10',
           '--time'             : '0-12:00:00',  # d-hh:mm:ss
           '--account'          : 'SNIC2019-2-1'
          }
    
}

MODULES = {
'KEBNE' : ['CUDA/10.1.243',
           'iccifort/2019.5.281',
           'impi/2018.5.288',
           'PyCUDA/2019.1.2-Python-3.7.4',
          ]
}