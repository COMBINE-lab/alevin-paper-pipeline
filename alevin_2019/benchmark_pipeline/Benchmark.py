import CGATCore.IOTools as IOTools

def makeCounterStatement(infile, outfile):
    '''Count number of reads in input files.
    This code is essentially a direct copy from
    CGATPipelines.PipelineMapping.Counter():
    https://github.com/CGATOxford/CGATPipelines
    '''

    statement = '''
        zcat %(infile)s
        | awk '{n+=1;} END {printf("nreads\\t%%%%i\\n",n/4);}'
        > %(outfile)s ''' % locals()
    
    return(statement)
