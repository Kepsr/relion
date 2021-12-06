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


def strip_comment(line):
    line = line.strip()
    if '#' in line:
        return line[:line.find('#')]
    return line


@safe
def load(filename):

    datablocks = OrderedDict()

    loop_state = 0
    # 0: outside
    # 1: reading datanames
    # 2: reading dataitems

    blockcodepattern = re.compile(r'data_(.*)')
    datanamepattern = re.compile(r'_(.*)')

    with open(filename) as file:
        for line in map(strip_comment, file):

            if not line:
                if loop_state == 2:
                    loop_state = 0
                continue
            
            blockcodematch = blockcodepattern.match(line)
            datanamematch = datanamepattern.match(line)

            if blockcodematch:
                # Enter a data block
                loop_state = 0
                blockcode = blockcodematch.group(1)
                datablock = OrderedDict()
                datablocks[blockcode] = datablock

            elif line == "loop_":
                # Enter a data loop
                loop_state = 1
                datanames = []

            elif datanamematch:
                tokens = datanamematch.group(1).split()
                # If inside a loop (reading data names)
                if loop_state == 1:
                    dataname = tokens[0]
                    datanames.append(dataname)
                    # Begin new array
                    datablock[dataname] = []
                else:
                    loop_state = 0
                    # Make sure outside loop
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
