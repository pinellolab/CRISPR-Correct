###
### Read in sequences from a pre-parsed reporter file (the reporter file is deprecated, and instead replaced with the UMI-tools output parsing)
### 
import csv
import multiprocessing as mp
from collections import Counter
from functools import reduce
from typeguard import typechecked
from functools import partial

def gen_chunks(reader, chunksize=1000):
    """
    Chunk generator. Take a CSV `reader` and yield
    `chunksize` sized slices.
    Source: https://gist.github.com/miku/820490
    """
    chunk = []
    for index, line in enumerate(reader):
        if index % chunksize == 0 and index > 0:
            yield chunk
            chunk = []
        chunk.append(line)
    yield chunk

def process(tsv_chunk, include_surrogate = False):
    local_counter = Counter()
    
    for row in tsv_chunk.copy():
        if include_surrogate:
            local_counter[(row[1], row[2], row[3])] += 1
        else:
            local_counter[row[1]] += 1
    
    return local_counter


@typechecked
def map_sample_protospacers(parsing_demult_handler, include_surrogate = False, cores=1):
    tsv_reader = csv.reader(parsing_demult_handler, delimiter='\t')  # change delimiter for normal csv files
    header = next(tsv_reader)
      
    read_chunks = gen_chunks(tsv_reader, chunksize=10000)

    process_func = partial(process, include_surrogate=include_surrogate)

    combined_counter = None
    if cores > 1:
        with mp.Pool(processes=cores) as pool:
            chunk_counters = pool.map(process_func, read_chunks)
            combined_counter = reduce(lambda x, y: x + y, chunk_counters)
    else:
        chunk_counters = [process_func(read_chunk) for read_chunk in read_chunks]
        combined_counter = reduce(lambda x, y: x + y, chunk_counters)

    return combined_counter 


