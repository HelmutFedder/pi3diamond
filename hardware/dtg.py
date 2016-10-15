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

# ToDo: convert to class

import numpy

from visa import instrument
import struct
import dtg_io
            
DTG = instrument('GPIB0::1::INSTR')

#DTG.chunk_size=2**20
DTG.timeout=100
ChunkSize = 1000000

# block_granularity and min_block_size: Constants of the DTG. 
block_granularity = 4
min_block_size = 1000 #bytes

#DTG.write('TBAS:CRAN 15')
#DTG.write('TBAS:SMODe SOFTware')

#def WriteVector(Block, Channel, data, offset ):
#    s = '"'
#    for n in data:
#        s+= str(n)
#    s += '"'
#    DTG.write('BLOC:SEL "'+Block+'"')
#    DTG.write('SIGN:DATA "Group1['+str(Channel)+']", %i, %i, '%(offset, len(data)) + s )

ChannelMap = {'laser':0, 'mw':1, 'mw_a':1, 'mw_b':4, 'mw_y':1, 'mw_x':4, 'sequence':3, 'aom':2, 'ch0':0, 'ch1':1, 'ch2':2, 'ch3':3}

path_on_this_pc = 'Z:\setup.dtg'
path_on_dtg = 'C:\pulsefiles\setup.dtg'
#ChannelMapBin = {'laser':'\x00', 'mw':'\x01', 'trigger':'\x04', 'SequenceTrigger':'\x02'}

def ChannelsToInt(channels):
    bits = 0
    for channel in channels:
        bits |= 1 << ChannelMap[channel]
    return bits

def gen_block(length, id, name):
    """
    gen_block: generate a tree for a DTG block of given length, id and name
    'name' is assumed to be packed into a numpy array. 
    'length' and 'id' are integers
    """
    length_entry = numpy.array([int(length)])
    return [
                ['DTG_BLOCK_RECID',
                30,
                [
                ['DTG_ID_RECID', 2, numpy.array([id], dtype=numpy.int16)],
                ['DTG_NAME_RECID', name.nbytes, name],
                ['DTG_SIZE_RECID', length_entry.nbytes, length_entry]]]
            ]
def gen_pattern(blockid, pulses):
    """
    gen_pattern: generate a python tree for a pattern command to fill block 'blockid' 
    with the binary sequence 'pulses'
    """
    return [
                ['DTG_PATTERN_RECID',
                1,
                [['DTG_GROUPID_RECID', 2, numpy.array([0], dtype=numpy.int16)],
                ['DTG_BLOCKID_RECID', 2, numpy.array([blockid], dtype=numpy.int16)],
                ['DTG_PATTERNDATA_RECID',
                1,
                pulses]
                ]]                    
            ]

def gen_sequence(label, subname, Nrep, goto):
    """
    gen_sequence: generate a python tree for a sequence entry of a given label, 
    name of the subsequence, number of repetitions 'Nrep' and goto label 'goto'. 
    all strings are assumed to by numpy arrays. Nrep is an int
    """
    return  [
                ['DTG_MAINSEQUENCE_RECID',
                 51,
                 [
                 ['DTG_LABEL_RECID', label.nbytes, label],
                 ['DTG_WAITTRIGGER_RECID', 2, numpy.array([0], dtype=numpy.int16)],
                 ['DTG_SUBNAME_RECID', subname.nbytes, subname],
                 ['DTG_REPEATCOUNT_RECID', 4, numpy.array([Nrep], dtype=numpy.int32)], # was 0,0
                 ['DTG_JUMPTO_RECID', 1, numpy.array([0], dtype=numpy.int8)],
                 ['DTG_GOTO_RECID', goto.nbytes, goto]]
                ],                       
             ]      
        
def BlockToDtgTree(BLOCK):
    pulses = numpy.array([], dtype=numpy.int8)
    
    total = 0 # total length of last subsequence
    subseq = 2 #index to the current subsequence
    
    break_len = 1000 # length of the macrobreak
    break_entry = numpy.array([int(break_len)])
    start_label = numpy.array('START\x00', dtype='|S')

    #blocks, sequence, patterns: trees holding the sequence data. 
    #After initialization, they contain the break block only
    blocks = [
        ['DTG_BLOCK_RECID',
        30,
        [['DTG_ID_RECID', 2, numpy.array([0], dtype=numpy.int16)],
        ['DTG_NAME_RECID', 6, numpy.array('BREAK\x00', dtype='|S')],
        ['DTG_SIZE_RECID', break_entry.nbytes, break_entry]]],
    ]
    sequence = [
                ]
    
    patterns = [
        ['DTG_PATTERN_RECID',
        1,
        [['DTG_GROUPID_RECID', 2, numpy.array([0], dtype=numpy.int16)],
        ['DTG_BLOCKID_RECID', 2, numpy.array([0], dtype=numpy.int16)],
        ['DTG_PATTERNDATA_RECID',
        1,
        numpy.zeros((1,break_len), dtype = numpy.int8)]
        ]]
    ]
    
    complete = numpy.sum([length for channels,length in BLOCK])
    consumed = 0
    
    for channels, length in BLOCK:
        length = round(length)
        consumed += length
        #check whether the current step is a break. If so, split it up into a macrobreak (using the BREAK blcok)
        #and its remainder. Start a new sequence from the remainder
        #Do not do this for the last min_block_size part of the sequence. This will be appended as a single block
        if channels == [] and consumed < complete - min_block_size:
            #if the current sequence ('pulses' of length 'total') is not a multiple of block_granularity 
            #and not long enough yet, extend it by deducing a tax from the break
            tax = (block_granularity - total%block_granularity) % block_granularity
            tax += max(min_block_size - (total+tax), 0)
            tax = min(tax, length) # do not tax more than the whole break
            total += tax
            length -= tax
            pulses = numpy.append(pulses, ChannelsToInt(channels)*numpy.ones((1,tax), dtype=numpy.int8))
            
            Tmacro = int(length/break_len)
            if Tmacro > 0:

                # if we have accumulated a sequence, dump it into the tree
                if total > 0: 
                    blockname = numpy.array('BLOCK' + str(subseq) + '\x00', dtype='|S') 

                    blocks.extend(gen_block(int(total), subseq, blockname))
                    patterns.extend(gen_pattern(subseq, pulses))
                    sequence.extend(gen_sequence(numpy.array('\x00', dtype='|S'),
                                    blockname,
                                    numpy.array([1], dtype=numpy.int32),
                                    numpy.array('\x00', dtype='|S')
                                    ))
                    total = 0
                                    
                # append Tmacro repetitions of the macrobreak
                sequence.append(
                    ['DTG_MAINSEQUENCE_RECID',
                    51,
                    [
                    ['DTG_LABEL_RECID', 1, numpy.array([0], dtype=numpy.int8)],
                    ['DTG_WAITTRIGGER_RECID', 2, numpy.array([0], dtype=numpy.int16)],
                    ['DTG_SUBNAME_RECID', 6, numpy.array('BREAK\x00', dtype='|S')],
                    ['DTG_REPEATCOUNT_RECID', 4, numpy.array([Tmacro], dtype=numpy.int32)], # was 0,0
                    ['DTG_JUMPTO_RECID', 1, numpy.array([0], dtype=numpy.int8)],
                    ['DTG_GOTO_RECID', 1, numpy.array([0], dtype=numpy.int8)]
                    ]]
                    )

                # append the remainder to the binary sequence
                length = length - Tmacro*break_len
                total = length
                subseq +=1
                pulses = ChannelsToInt(channels)*numpy.ones((1,length), dtype=numpy.int8)
            else: # Tmacro=0, append the whole break to the sequence
                total += length
                pulses = numpy.append(pulses, ChannelsToInt(channels)*numpy.ones((1,length), dtype=numpy.int8))
        else: # mw or laser pulse. Append it to the sequence
            total += length
            pulses = numpy.append(pulses, ChannelsToInt(channels)*numpy.ones((1,length), dtype=numpy.int8))
    
    #append the accumulated pulse sequence
    if total > 0:
        blockname = numpy.array('BLOCK' + str(subseq) + '\x00', dtype='|S') 
        blocks.extend(gen_block(int(total), subseq, blockname))
        patterns.extend(gen_pattern(subseq, pulses))
        sequence.extend(gen_sequence(numpy.array('\x00', dtype='|S'),
                        blockname,
                        numpy.array([1], dtype=numpy.int32),
                        numpy.array('\x00', dtype='|S')
                        ))

    #finally, label the first subsequence as "start", point the goto of the last subsequence there
    sequence[0][2] = dtg_io.change_leaf(sequence[0][2], 'DTG_LABEL_RECID', 5, numpy.array('START\x00', dtype='|S'))
    sequence[-1][2] = dtg_io.change_leaf(sequence[-1][2], 'DTG_GOTO_RECID', 5, numpy.array('START\x00', dtype='|S'))
                    
    binary_pattern = []
    binary_pattern.extend(blocks)
    binary_pattern.extend(sequence)
    binary_pattern.extend(patterns)
    
    # construct a DTG tree by appending binary_pattern to the scaffold tree. 
    dtgfile = dtg_io.load_scaffold()
    view = dtgfile[-1][-1].pop()
    
    for node in binary_pattern:
        dtgfile[-1][-1].append(node)
    
    dtgfile[-1][-1].append(view)
    dtgfile = dtg_io.recalculate_space(dtgfile)[1]
    
    return dtgfile

    
def Sequence(BLOCK, loop=False):
    Stop()
    # set sequence length
    #if loop:
    #    DTG.write('SEQ:LENG inf')
    #else:
    length = sum([round(plen) for (mask, plen) in BLOCK])
    m = length % 4
    
    #if block size is not a multiple of 4, increase first laser wait to avoid block size granularity error
    if length % 4:  
        BLOCK[2] = (BLOCK[2][0], BLOCK[2][1] + 4-m)
    
    nlength = sum([round(plen) for (mask, plen) in BLOCK])
    #print "adjusted length from " + str(length) + " to " + str(nlength) + "\n"
    
    # generate a setup tree for the DTG containing the new pulse sequence
    dtgtree = BlockToDtgTree(BLOCK)
    
    # write this file onto the DTG (via Samba)
    fil = open(path_on_this_pc, 'wb')
    dtg_io.dump_tree(dtgtree, fil)
    fil.close()
    
#    fil = open('D:/Python/pi3diamond/current_pattern.dtg', 'wb')
#    dtg_io.dump_tree(dtgtree, fil)
#    fil.close()

    #load and execute the generated setup file
    DTG.write('MMEM:LOAD "'+path_on_dtg+'"')
    DTG.write('OUTP:STAT:ALL ON')
    while not int(DTG.ask('TBAS:RUN?')):
        DTG.write('TBAS:RUN ON')
    
    return nlength

def State(flag = None):
    if flag is not None:
        if (flag):
            DTG.write('TBAS:RUN ON')
        else:
            DTG.write('TBAS:RUN OFF')
    return DTG.ask('TBAS:RUN?')

def Run():
    DTG.write('OUTP:STAT:ALL ON')
    DTG.write('TBAS:RUN ON')
    return DTG.ask('TBAS:RUN?')

def Stop():
    DTG.write('TBAS:RUN OFF')
    DTG.write('OUTP:STAT:ALL OFF')
    return DTG.ask('TBAS:RUN?')

def Light():
    Stop()
    Sequence( [  (['laser'],  1000), ] )

def Night():
    Stop()
    Sequence( [  ([],  1000), ] )

def Open():
    Stop()
    Sequence( [  (['laser','mw'],  1000), ] )

def High(chans): #chans: sequence of channel id strings, like ['laser', 'mw']
    Stop()
    Sequence( [  (chans,  1000), ] )
