#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scvi_import_test.$JOB_ID.log
#$ -l s_vmem=32G
source /home/xzy0723/miniconda3/etc/profile.d/conda.sh
conda activate scvi
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
echo "host=$(hostname) vmem=$(ulimit -v)"
python -c "import scvi,torch,scanpy,jax; print('IMPORT OK | scvi',scvi.__version__,'| torch',torch.__version__,torch.version.cuda,'| jax',jax.__version__)"
