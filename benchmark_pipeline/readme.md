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
* CGAT-core requires these python libraries to be installed before installing it
* https://github.com/CGATOxford/CGATCore/blob/eceb5200efa2b81a00c60cd417887e565a90edd8/requires.txt

git clone git@github.com:CGATOxford/CGATCore.git
cd CGATCore; python setup.py install; cd ..

git clone git@github.com:COMBINE-lab/alevin-paper-pipeline.git
python ../alevin-paper/benchmark_pipeline/pipeline_10xbenchmark.py -c1 -e -p2 -v5 show run_alevin_options --local
