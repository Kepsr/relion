'''
Module for reading STAR (Self-defining Text Archive and Retrieval) files
'''
from collections import OrderedDict
import re
import time


def recursivelydescend(d, keys):
    for k in keys:
        d = d[k]
    return d


def safe(loadmethod, max_tries=5, wait=10):
    '''
    Make a `loadmethod` safe.
    '''
    def safemethod(filename, expected=[]):
        for _ in range(max_tries):
            try:
                starfile = loadmethod(filename)
                # Ensure the expected keys are present
                # By descending through the dictionary
                entry = starfile
                recursivelydescend(entry, expected)
                return starfile
            except KeyError:
                print("Just tried (and failed) to read {}, expected key: {}".format(filename, expected))
                time.sleep(wait)
        raise Exception("Failed to read a star file: {}".format(filename))
    return safemethod


# @safe
# def load(filename):
#     '''
#     Load a STAR file, returning the datasets (of type `OrderedDict`).
#     '''

#     datasets = OrderedDict()
#     datablock, datanames = None, None

#     in_loop = 0
#     # 0: outside
#     # 1: reading colnames
#     # 2: reading data

#     with open(filename) as readfile:
#         datapattern = re.compile(r'data_(.*)')
#         looppattern = re.compile(r'loop_')
#         namepattern = re.compile(r'_(.*)')
#         for line in map(str.strip, readfile):

#             # Remove comments
#             comment_pos = line.find('#')
#             line = line[:comment_pos] if comment_pos > 0 else line

#             if not line:
#                 in_loop = 0 if in_loop == 2 else in_loop
#                 continue

#             datamatch = datapattern.match(line)
#             loopmatch = looppattern.match(line)
#             namematch = namepattern.match(line)

#             if datamatch:
#                 # Start parsing data block
#                 blockname = datamatch.group(1)
#                 datablock = OrderedDict()
#                 datasets[blockname] = datablock
#                 in_loop = 0

#             elif loopmatch:
#                 # Start parsing data loop
#                 datanames = []
#                 in_loop = 1

#             elif namematch:
#                 # Read data name
#                 dataitems = namematch.group(1).split()
#                 if len(dataitems) < 2:
#                     print(line)
#                 dataname, dataitem = dataitems[:2]
#                 if in_loop == 1:
#                     # We are inside a data loop
#                     datablock[dataname] = []
#                     datanames.append(dataname)
#                 else:
#                     # We are outside a data loop
#                     datablock[dataname] = dataitem
#                 if in_loop == 2:
#                     in_loop = 0

#             elif in_loop in (1, 2):
#                 in_loop = 2
#                 dataitems = line.split()
#                 assert len(dataitems) == len(datanames), (
#                     "Error in STAR file {}, number of elements in {} does not match number of column names {}".format(
#                         filename, dataitems, datanames
#                     )
#                 )
#                 for dataname, dataitem in zip(datanames, dataitems):
#                     datablock[dataname].append(dataitem)

#     return datasets


@safe
def load(filename):

    datasets = OrderedDict()
    current_data = None
    current_colnames = None

    in_loop = 0
    # 0: outside
    # 1: reading colnames
    # 2: reading data

    for line in open(filename):
        line = line.strip()
        # remove comments
        if '#' in line:
            line = line[:line.find('#')]

        if not line:
            if in_loop == 2:
                in_loop = 0
            continue

        if line.startswith("data_"):
            in_loop = 0

            data_name = line[5:]
            current_data = OrderedDict()
            datasets[data_name] = current_data

        elif line.startswith("loop_"):
            current_colnames = []
            in_loop = 1

        elif line.startswith("_"):
            if in_loop == 2:
                in_loop = 0

            elems = line[1:].split()
            key = elems[0]
            if in_loop == 1:
                current_colnames.append(key)
                current_data[key] = []
            else:
                value = elems[1]
                current_data[key] = value

        elif in_loop:
            in_loop = 2
            elems = line.split()
            assert len(elems) == len(current_colnames), (
                "Error in STAR file {}, number of elements in {} does not match number of column names {}".format(filename, elems, current_colnames)
            )
            for i, elem in enumerate(elems):
                current_data[current_colnames[i]].append(elem)

    return datasets
