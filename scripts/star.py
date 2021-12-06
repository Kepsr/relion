'''
Module for reading STAR (Self-defining Text Archive and Retrieval) files.
See:

https://www.iucr.org/__data/assets/file/0013/11416/star.5.html
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
#     Load a STAR file, returning the datablocks (of type `OrderedDict`).
#     '''

#     datablocks = OrderedDict()
#     datablock, datanames = None, None

#     loop_state = 0
#     # 0: outside
#     # 1: reading datanames
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
#                 loop_state = 0 if loop_state == 2 else loop_state
#                 continue

#             datamatch = datapattern.match(line)
#             loopmatch = looppattern.match(line)
#             namematch = namepattern.match(line)

#             if datamatch:
#                 # Start parsing data block
#                 blockcode = datamatch.group(1)
#                 datablock = OrderedDict()
#                 datablocks[blockcode] = datablock
#                 loop_state = 0

#             elif loopmatch:
#                 # Start parsing data loop
#                 datanames = []
#                 loop_state = 1

#             elif namematch:
#                 # Read data name
#                 dataitems = namematch.group(1).split()
#                 if len(dataitems) < 2:
#                     print(line)
#                 dataname, dataitem = dataitems[:2]
#                 if loop_state == 1:
#                     # We are inside a data loop
#                     datablock[dataname] = []
#                     datanames.append(dataname)
#                 else:
#                     # We are outside a data loop
#                     datablock[dataname] = dataitem
#                 if loop_state == 2:
#                     loop_state = 0

#             elif loop_state in (1, 2):
#                 loop_state = 2
#                 dataitems = line.split()
#                 assert len(dataitems) == len(datanames), (
#                     "Error in STAR file {}, number of elements in {} does not match number of column names {}".format(
#                         filename, dataitems, datanames
#                     )
#                 )
#                 for dataname, dataitem in zip(datanames, dataitems):
#                     datablock[dataname].append(dataitem)

#     return datablocks


def strip_comment(line):
    line = line.strip()
    if '#' in line:
        return line[:line.find('#')]
    return line


@safe
def load(filename):

    datablocks = OrderedDict()
    datablock = None
    datanames = None

    loop_state = 0
    # 0: outside
    # 1: reading datanames
    # 2: reading dataitems

    for line in map(strip_comment, open(filename)):

        if not line:
            if loop_state == 2:
                loop_state = 0
            continue

        if line.startswith("data_"):
            # Enter a data block
            loop_state = 0
            blockcode = line[5:]
            datablock = OrderedDict()
            datablocks[blockcode] = datablock

        elif line.startswith("loop_"):
            # Enter a data loop
            loop_state = 1
            datanames = []

        elif line.startswith("_"):
            # If inside a loop (reading data items)
            if loop_state == 2:
                # Exit loop
                loop_state = 0
            tokens = line[1:].split()

            # If inside a loop (reading data names)
            if loop_state == 1:
                dataname = tokens[0]
                datanames.append(dataname)
                datablock[dataname] = []
            # If outside a loop
            else:
                dataname, dataitem = tokens[:2]
                datablock[dataname] = dataitem

        # If inside a loop
        elif loop_state:
            # Read dataitems
            loop_state = 2
            dataitems = line.split()
            assert len(datanames) == len(dataitems), (
                "Error in STAR file {}, number of data items in {} does not match number of data names {}".format(filename, dataitems, datanames)
            )
            for dataname, dataitem in zip(datanames, dataitems):
                datablock[dataname].append(dataitem)
        
        # else:
        #     raise SyntaxError(f'Line could not be parsed: {line}')

    return datablocks
