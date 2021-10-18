'''star.py
Module for reading STAR (Self-defining Text Archive and Retrieval) files.
'''
from collections import OrderedDict
import re


def load_star(filename):
    """
    Load a STAR (Self-defining Text Archive and Retrieval) file,
    returning the datasets (of type `OrderedDict`).
    """

    datasets = OrderedDict()
    datablock, datanames = None, None

    in_loop = 0
    # 0: outside
    # 1: reading colnames
    # 2: reading data

    with open(filename) as readfile:
        datapattern = re.compile(r'data_(.*)')
        looppattern = re.compile(r'loop_')
        namepattern = re.compile(r'_(.*)')
        for line in map(str.strip, readfile):

            # Remove comments
            comment_pos = line.find('#')
            line = line[:comment_pos] if comment_pos > 0 else line

            if not line:
                in_loop = 0 if in_loop == 2 else in_loop
                continue

            datamatch = datapattern.match(line)
            loopmatch = looppattern.match(line)
            namematch = namepattern.match(line)

            if datamatch:
                # Start parsing data block
                blockname = datamatch.group(1)
                datablock = OrderedDict()
                datasets[blockname] = datablock
                in_loop = 0

            elif loopmatch:
                # Start parsing data loop
                datanames = []
                in_loop = 1

            elif namematch:
                # Read data name
                dataitems = namematch.group(1).split()
                dataname, dataitem = dataitems[:2]
                if in_loop == 1:
                    # We are inside a data loop
                    datablock[dataname] = []
                    datanames.append(dataname)
                else:
                    # We are outside a data loop
                    datablock[dataname] = dataitem
                if in_loop == 2:
                    in_loop = 0

            elif in_loop in (1, 2):
                in_loop = 2
                dataitems = line.split()
                assert len(dataitems) == len(datanames), (
                    "Error in STAR file {}, number of elements in {} does not match number of column names {}".format(
                        filename, dataitems, datanames
                    )
                )
                for dataname, dataitem in zip(datanames, dataitems):
                    datablock[dataname].append(dataitem)

    return datasets

