Benchmark pipeline
------------------

CGATPipelines pipeline to run and benchmark various tools.

Plan is to have self contained downloading of data.
Will assume CellRanger, Kallisto and umi-tools are installed.

Will generate as much of references and indexes as possible.

Plan is to have a directory `benchmark_pipeline` and then a src dir within that.
Pipeline will run in the outer directory (i.e. witin the alevin-paper repo). Outputs
can be added to the github.


Issues
------

* How do we make everything go to nodes with same CPUs?
   - could name a specific node, but that might make the whole cluster submission thing
     a bit pointless.
     
* I thought we had agreed to use hgmm datasets. However, these are not compatible with
  kallisto I don't think. Or they might be but I can't figure out the Kallisto code
  enough to know.


ToDo:
------

- [X] Tasks to download data
- [X] Tasks to merge references, genesets etc.
- [X] Annotated genesets with tx2gene, mito_genes, rRNA genes
- [X] Implement alevin subpipe
   - [X] Test alevin tasks
- [X] Write timing file parser
- [X] Find set of job_options that will mean all tasks will run on machines with
     same processor
- [X] implement umi-tools subpipe
- [X] implement CellRanger subpipe
- [X] find either how to run v1 chemistry with alevin or v2 with kallisto
- [ ] implement kalisto
- [ ] run
- [ ] figure
