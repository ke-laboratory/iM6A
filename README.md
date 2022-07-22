## iM6A: modeling m6A site-specific deposition by deep learning
#### Setup
##### Create environment
```bash
cd /pod/2/ke-lab/LUOZ

# Create writable containers
srun --job-name=bash --ntasks=1 --cpus-per-task=4 --time=8:00:00 --gres-flags=enforce-binding --mem=40000 --pty bash
module load singularity
singularity build --sandbox tensorflow1901/ docker://nvcr.io/nvidia/tensorflow:19.01-py2

# Open Jupyter notebook
singularity shell --writable-tmpfs --bind /pod/2/ke-lab/LUOZ:/mnt tensorflow1901/
PORT=$(shuf -i10000-11999 -n1)
echo "executing jupyter on http://$(hostname):$PORT"
jupyter notebook --no-browser --port=$PORT --ip=`hostname -i`

# Install packages into an external path:
conda create -n iM6A python=2.7
conda activate iM6A
cd /pod/2/ke-lab/LUOZ
pip install --target=/Singularity/iM6A biopython==1.76
pip install --target=/Singularity/iM6A scikit-learn==0.20.3
pip install --target=/Singularity/iM6A matplotlib==2.2.4
pip install --target=/Singularity/iM6A keras==2.0.5
```
##### Sbatch job.slurm
```bash
#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --mail-type=BEGIN,END
#SBATCH -q training
#SBATCH --partition=gpu
#SBATCH -t 120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem=40000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err       
module load singularity
PORT=$(shuf -i10000-11999 -n1)
echo "executing jupyter on http://$(hostname):$PORT"
singularity run --writable-tmpfs --bind /pod/2/ke-lab/LUOZ:/mnt --nv tensorflow1901/ jupyter notebook --no-browser --port=$PORT --ip=`hostname -i`
```

#### m6A Prediction
```python
import sys
sys.path.append("Singularity/iM6A")

from keras.models import load_model
from pkg_resources import resource_filename
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import keras.backend as kb

def one_hot_encode(seq):
    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]

def categorical_crossentropy_2d(y_true, y_pred):
    return - kb.mean(y_true[:, :, 0]*kb.log(y_pred[:, :, 0]+1e-10)
                   + y_true[:, :, 1]*kb.log(y_pred[:, :, 1]+1e-10))

paths = ('Model/mouseRRACH10000_c{}.h5'.format(x) for x in range(1, 6))
models = [load_model(y, custom_objects={'categorical_crossentropy_2d': categorical_crossentropy_2d}) for y in paths]
context = 10000

## Replace the input with your custom sequence
input_sequence = "CCAGCAATGAGCAGTCTTCTTTCTCCTGTCTCCTAGGCCTACTTAGGAGCCTACGTTTTTATCATCTTCGCTGCCTTCCTCATCTTCTTCCTAATCTTCACCTTCTTCAAAGTCCCGGAGACCAAAGGCAGGACTTTCGAGGACATTGCCCGGGCCTTCGAGGGGCAGGCGCACTCTGGAAAAGGCCCTGCCGGTGTGGAGTTGAACAGCATGCAGCCGGTCAAGGAGACCCCTGGCAACGCCTGAGCCGGGAGCACCTCCTTCACCTCCCTCCACTGTGGAAAGCCACCTCCCCTTAAGTGGGGAGACCTCATCAGGATGAACCAGGACTGCTTCTGAGTGCTGCTATTCTCACCCCACACACTCCATGAAAACTCAGCTGCACCCAATGCTGGGGTTGACCAGATCGCCAATGACTTTTAAGTGTTTGATTTGGGGGATACTTCCCTTGTAATCAGGAAAGACCAAGGAAGCCTACCTTTATATTGGGAGGGAAGGGGCCGCAGCTCCCCTTAGGTTCTAAAACCCGCTAACTAGGACAACTAGGGAGTAGGAGTAGGGTGCAGCACCGCCCCCACCACCACCTTAGTTTGCCAGGAAACAGATCTTCATACCCAGTGTGGAGGGCCCTGGGGGATTGAACAAAAGACCCCCTCCTCACTTGATACAGCTCTGCACAGCAAAGTAACTTGAGTTTTATTTATTTTATCTTCAGGTTGAATTGCATAAATATTTATTTTTTTAATTGTAATTTTACCAAATAATGAGACAGTAACGAAATTGAGGTTGGAGGGAGGTGTT"

x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
Probability = y[0, :, 1]
```

## Contact

Zhiyuan Luo: luozhiyuan0717@hotmail.com

Shengdong Ke: kelab018@gmail.com
