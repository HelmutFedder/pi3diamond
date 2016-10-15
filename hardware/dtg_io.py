"""
    This file is part of pi3diamond, a toolkit for
    confocal scanning, anti-bunching, FLIM, pulsed ODMR / NMR,
    and more sophisticated quantum physics experiments,
    typically performed with NV centers in diamond,
    written in python using the enthought traits packages.

    pi3diamond is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pi3diamond is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diamond. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2009-2016 Helmut Fedder <helmut@fedder.net>
"""

import struct
import numpy
import binascii
import os

from dtg_recids import *

def read_tree(N, bytes):
    """
    read N bytes from a string 'bytes' and build up a tree
    """
    tree = []
    while N>0 and len(bytes) > 0:
        try:
            recid = struct.unpack('<H', bytes[:2])[0] 
        except:
            print "unpack impossible. offending string " + bytes[:2]
            return bytes, tree
        bytes = bytes[2:]
        N -= 2
        if (recid & internal_node):
            lint = struct.unpack('<I', bytes[:4])[0]
            bytes, subtree = read_tree(lint, bytes[4:])
            hr_recid = [key for key in recids.keys() if recids[key]==recid ]
            tree.append([hr_recid[0], lint, subtree])
            N -= lint + 4
        else: #leaf node
            data, lleaf, bytes = read_leaf_data(bytes)
            hr_recid = [key for key in recids.keys() if recids[key]==recid ]
            tree.append([hr_recid[0], lleaf-4, data])
            N -= lleaf
    return bytes, tree

def read_leaf_data(bytes):
    """
    read data from a leaf node, 
    bytes must start with an entry of type "length, data" (RECID must be chopped before)
    """
    len = struct.unpack('<I', bytes[:4])[0]
    if (len %2): 
        data = numpy.array(struct.unpack('<' + str(len) + 'B', bytes[4:4+len]), dtype=numpy.int8)
    else:
        data = numpy.array(struct.unpack('<' + str(len/2) + 'h', bytes[4:4+len]), dtype=numpy.short)        
    return [data, len + 4, bytes[4+len:]]

def dump_tree(tree, fil):
    """ 
    accept a DTG setup tree 'tree' and dump it into the file fil
    """
    for node in tree:
        recid_hr, len, data = node
        recid = recids[recid_hr]
        if recid & internal_node:
            fil.write(struct.pack('<H', recid))
            fil.write(struct.pack('<I', len))
            dump_tree(data, fil)
        else: #leaf node
            fil.write(struct.pack('<H', recid))
            fil.write(struct.pack('<I', len))
            fil.write(data.tostring())

def recalculate_space(tree):
    """
    given a DTG setup tree 'tree', update the size value of each node
    """
    total = 0
    for i, node in enumerate(tree):
        try:
            recid_hr, length, data = node
        except:
            print "could not unpack " + str(node) + '\n'
        
        recid = recids[recid_hr]
        if recid & internal_node:
            node[1], node[2] = recalculate_space(data)            
            tree[i] = node
        else: #leaf node
            node[1] = len(data.tostring())
            tree[i] = node
        total += node[1]+6 # add length of the treated node to total length
    return total, tree

def get_leaf(tree, leaf):
    """
    get_leaf: search a DTG tree 'tree' for a leaf with (human readable) recid 'leaf' and return it
    """
    for node in tree:
        recid_hr, dlen, data = node
        recid = recids[recid_hr]
        if recid & internal_node:
            res = get_leaf(data, leaf)
            if len(res):
                return res
        
        else:
            if recid_hr == leaf:
                return [dlen, data]
    
    return []

def change_leaf(leafs, leaf_id, newlen, newval):
    """
    change_leaf: descend the list of leafs 'leafs' and look for all leafs with (human readable string) id 'leaf_id'. 
    Once found, replace their content by 'newlen', 'newval'
    """
    for i, leaf in enumerate(leafs):
        if leaf[0] == leaf_id:
            leafs[i][1] = newlen
            leafs[i][2] = newval
    return leafs
    
def load_scaffold():
    """ 
    load_scaffold: load and return the DTG config scaffold
    """
    fil = open('scaffold.dtg', 'rb')
    bytes = fil.read()
    fil.close()
    return read_tree(12000, bytes)[1]

# global variable holding the DTG scaffold
scaffold = load_scaffold()

#test code to reorder entries in the file and dump it to disk
#view = tree[0][-1].pop()
#pattern = tree[0][-1].pop()

#tree[0][-1].append(view)
#tree[0][-1].append(pattern)

#pattern = get_leaf(tree, 'DTG_PATTERNDATA_RECID')[1]

#total, tree = recalculate_space(tree)

#fil = open('dtg.reversed.dat', 'wb')
#dump_tree(tree, fil)
#fil.close()
