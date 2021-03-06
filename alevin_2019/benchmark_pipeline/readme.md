Benchmark pipeline
------------------

CGATPipelines pipeline to run and benchmark various tools.

Plan is to have self contained downloading of data.
Will assume CellRanger, Kallisto and umi-tools are installed.

Will generate as much of references and indexes as possible.

Plan is to have a directory `benchmark_pipeline` and then a src dir within that.
Pipeline will run in the outer directory (i.e. witin the alevin-paper repo). Outputs
can be added to the github.

How to run
------
**Note**: We are very actively developing the pipeline, if you face any problem while replicating the pipeline please raise an issue in this repo.

```
# CGAT-core requires these python libraries to be installed before installing it
# https://github.com/CGATOxford/CGATCore/blob/eceb5200efa2b81a00c60cd417887e565a90edd8/requires.txt

git clone git@github.com:CGATOxford/CGATCore.git
cd CGATCore; python setup.py develop; cd ..

git clone git@github.com:COMBINE-lab/alevin-paper-pipeline.git
python alevin-paper/benchmark_pipeline/pipeline_10xbenchmark.py -c1 -e -p2 -v5 make full --local
```

### Working Settings on macosx
```
conda create -n cgat -c anaconda python=3
source activate cgat
pip install pysam
pip install numpy
conda install -c anaconda cython
conda install -c anaconda networkx
pip install drmaa
conda install -c r rpy2
```
