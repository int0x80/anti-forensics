#!/usr/bin/env python
# -----------------------------------------------------------
# helloyesthisisdog is a python script for embedding and 
# extracting parasite files into a host png file. several 
# options are available to protect your parasite file from
# detection or carving.
#
# when embedding a parasite file, the resulting file is named
# as 'infected.png' unless -o is used.
#
# the parasite file can first be 'encoded' with an xor using
# the -x option and specifying a key in hex representation. 
#
# your parasite file can also be split across two sections
# in the png so that carving tools will not extract a single,
# contiguous version of your parasite file. use -e for this.
#
# usage: helloyesthisisdog.py [options] <action> <action args>
# -----------------------------------------------------------
# 
# -----------------------------------------------------------
# actions
#   -m      <host file> <parasite file>
#
#           embeds parasite file into host file
#
#   -x      <host file> <parasite file>
#
#           extracts parasite file from host file
#
#   -p      <host file>
#
#           prints list of available ancillary png sections
#           these are sections to use for embedding files
# -----------------------------------------------------------
# options
#   -a
#
#           automatically find and use one available section.
#           this is for embedding a file contiguously.
#
#   -e
#
#           evil mode, split file across sections so
#           investigators don't get contiguous file extract.
#           this will randomly select a subset of sections.
#
#   -k key
#
#           xor the parasite file with key (hex).
#           e.g. helloyesthisisdog.py -k 41006100 embed host.png bad.exe
#
#   -o output_file
#
#           write results to output_file.
#           defaults to 'infected.png'.
#
#   -s section1,section2,...
#
#           embed/extract using specific section(s).
#           the order needs to remain consistent when 
#           embedding/extracting across multiple sections.
# -----------------------------------------------------------

from struct import *
from itertools import izip, cycle
from random import choice, seed, randint
import binascii, getopt, sys

png_sections = [
    'bKGD',
    'cHRM',
    'gAMA',
    'hIST',
    'iCCP',
    'iTXt',
    'pHYs',
    'sBIT',
    'sPLT',
    'sRGB',
    'sTER',
    'tEXt',
    'tIME',
    'tRNS',
    'zTXt'
]

# -----------------------------------------------------------
# embed contiguous file contents into a section
# starting at the back of the file 
# -----------------------------------------------------------
def embed_file(host_file_contents, chunk_type, parasite_file_contents, output_file):
    # -----------------------------------------------------------
    # find length and crc32 of parasite file
    # -----------------------------------------------------------
    file_len = len(parasite_file_contents)
    file_crc = binascii.crc32(parasite_file_contents) & 0xffffffff
    file_len_data = pack('>L', file_len)
    file_crc_data = pack('>L', file_crc)

    # -----------------------------------------------------------
    # parse the host file
    # -----------------------------------------------------------
    iend_offset = find_section_reverse(host_file_contents, 'IEND')
    host_head = host_file_contents[0:iend_offset]
    host_tail = host_file_contents[iend_offset:]
    
    # -----------------------------------------------------------
    # construct the new host file with parasite embedded
    # -----------------------------------------------------------
    infected_file = open(output_file, 'wb')
    infected_file.write(host_head)
    infected_file.write(file_len_data)
    infected_file.write(chunk_type)
    infected_file.write(parasite_file_contents)
    infected_file.write(file_crc_data)
    infected_file.write(host_tail)
    infected_file.close()

    return True


# -----------------------------------------------------------
# embed file into randomly selected available chunks
# -----------------------------------------------------------
def evil_embed(host_file_contents, host_free_sections, parasite_file_contents, output_file):
    # -----------------------------------------------------------
    # parse the host file
    # -----------------------------------------------------------
    iend_offset = find_section_reverse(host_file_contents, 'IEND')
    host_head = host_file_contents[0:iend_offset]
    host_tail = host_file_contents[iend_offset:]
    
    # -----------------------------------------------------------
    # construct the new host file with original start contents
    # -----------------------------------------------------------
    infected_file = open(output_file, 'wb')
    infected_file.write(host_head)
    
    # -----------------------------------------------------------
    # select random number of chunks, and random chunks to use
    # -----------------------------------------------------------
    num_chunks = randint(2, len(host_free_sections))
    if num_chunks == len(host_free_sections):
        evil_chunks = host_free_sections
    else:
        evil_chunks = {}
        while len(evil_chunks) < num_chunks:
            chunk = choice(host_free_sections)
            evil_chunks[chunk] = True
            host_free_sections.remove(chunk)
        
    # -----------------------------------------------------------
    # let's get it on
    # -----------------------------------------------------------
    chunks_left = num_chunks
    slice_count = 0
    offset = 0
    for chunk_type in evil_chunks.keys():
        # -----------------------------------------------------------
        # special case for last slice to write remainder of file
        # -----------------------------------------------------------
        if slice_count == (num_chunks - 1):
            file_slice = parasite_file_contents[offset:]
            offset = offset + len(file_slice)
        # -----------------------------------------------------------
        # else randomly determine chunk size and read that many bytes
        # -----------------------------------------------------------
        else:
            slice_size_max = (len(parasite_file_contents) - offset) / chunks_left
            slice_size = randint(3, slice_size_max)
            file_slice = parasite_file_contents[offset:(offset + slice_size)]
            chunks_left = chunks_left - 1
            offset = offset + slice_size
        
        # -----------------------------------------------------------
        # find length and crc32 of parasite file slice
        # -----------------------------------------------------------
        file_len = len(file_slice)
        file_crc = binascii.crc32(file_slice) & 0xffffffff
        # print 'Slice: %d \t Length: %d \t Offset: %04d \t CRC32: %08x' % (slice_count, file_len, offset, file_crc)
        file_len_data = pack('>L', file_len)
        file_crc_data = pack('>L', file_crc)

        # -----------------------------------------------------------
        # write the chunk out to file
        # -----------------------------------------------------------
        infected_file.write(file_len_data)
        infected_file.write(chunk_type)
        infected_file.write(file_slice)
        infected_file.write(file_crc_data)

        # -----------------------------------------------------------
        # on to the next slice
        # -----------------------------------------------------------
        slice_count = slice_count + 1    
    
    # -----------------------------------------------------------
    # close out the new host file
    # -----------------------------------------------------------
    infected_file.write(host_tail)
    infected_file.close()

    return evil_chunks


# -----------------------------------------------------------
# extract embedded contiguous file from infected host file
# -----------------------------------------------------------
def extract_file(host_file_contents, chunk_type, output_file):
    offset = find_section(host_file_contents, chunk_type)
    file_len = unpack('>L', host_file_contents[offset:(offset+4)])[0]
    parasite_file_contents = host_file_contents[(offset+8):(offset+8+file_len)]
    parasite_file = open(output_file, 'wb')
    parasite_file.write(parasite_file_contents)
    parasite_file.close()

    return True


# -----------------------------------------------------------
# retrieve, presumably, a portion of an embedded file
# -----------------------------------------------------------
def extract_section(host_file_contents, chunk_type):
    offset = find_section(host_file_contents, chunk_type)
    file_len = unpack('>L', host_file_contents[offset:4])
    chunk_data = host_file_contents[(offset+8):file_len]

    return chunk_data


# -----------------------------------------------------------
# find offset of section in png file
# return value of 0 means the section is not in use
# -----------------------------------------------------------
def find_section(host_file_contents, chunk_type):
    offset = 0
    while offset < (len(host_file_contents) - len(chunk_type)):
        if host_file_contents[offset:(offset+len(chunk_type))] == chunk_type:
            return offset - 4
        offset = offset + 1

    return 0


# -----------------------------------------------------------
# find offset of section in png file, starting at back
# return value of 0 means the section is not in use
# -----------------------------------------------------------
def find_section_reverse(host_file_contents, chunk_type):
    offset = len(host_file_contents) - 8
    while offset > 0:
        if host_file_contents[offset:(offset+len(chunk_type))] == chunk_type:
            return offset - 4
        offset = offset - 1

    return 0


# -----------------------------------------------------------
# print available ancillary png sections
# -----------------------------------------------------------
def free_sections(host_file_contents, verbose=True):
    if verbose == True:
        print "Sections available:"
    count = 0
    host_free_sections = []
    for section in png_sections:
        available = find_section(host_file_contents, section)
        if available == 0:
            host_free_sections.append(section)
            if verbose == True:
                print "\t" + section
            count = count + 1

    if verbose == True:
        print "%d free section(s)" % count
    return host_free_sections


# -----------------------------------------------------------
# auxiliary function to grab file contents
# -----------------------------------------------------------
def get_file(dat_file):
    dat_file_handle = open(dat_file, 'rb')
    dat_file_contents = dat_file_handle.read()
    dat_file_handle.close()

    return dat_file_contents


# -----------------------------------------------------------
# use it, abuse it, confuse it, or lose it
# hope this usage helps, derp
# -----------------------------------------------------------
def usage():
    print '''
usage: ''' + sys.argv[0] + ''' [options] <action> <action args>

actions
  -m      <host file> <parasite file>

          embeds parasite file into host file

  -x      <host file> <parasite file>

          extracts parasite file from host file

  -p      <host file>

          prints list of available ancillary png sections
          these are sections to use for embedding files

options
  -a

          automatically find and use one available section.
          this is for embedding a file contiguously.

  -e

          evil mode, split file across sections so
          investigators don't get contiguous file extract.
          this will randomly select a subset of sections.

  -k key

          xor the parasite file with key (hex).
          e.g. helloyesthisisdog.py -k 41006100 embed host.png bad.exe

  -o output_file

          write results to output_file.
          defaults to 'infected.png'.

  -s section1,section2,...

          embed/extract using specific section(s).
          the order needs to remain consistent when 
          embedding/extracting across multiple sections.'''

    return True


# -----------------------------------------------------------
# xor the file contents against the supplied key
# helps prevent detection by ids and carving tools
# -----------------------------------------------------------
def xor_file(parasite_file_contents, key):
    return ''.join(chr(ord(file_byte) ^ ord(key_byte)) for (file_byte, key_byte) in izip(parasite_file_contents, cycle(key)))


def main():
    # -----------------------------------------------------------
    # grab arguments passed in
    # -----------------------------------------------------------
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'aehmxpk:o:s:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)

    if len(sys.argv) == 1:
        usage()
        sys.exit(1)

    # -----------------------------------------------------------
    # initialize variables
    # -----------------------------------------------------------
    bln_auto    = False
    bln_evil    = False
    sections    = []
    xor_key     = 0
    output_file = 'infected.png'
    chunk_type  = ''

    # -----------------------------------------------------------
    # first pass: options
    # -----------------------------------------------------------
    for option, argument in optlist:
        if option == '-a':
            bln_auto = True

        elif option == '-e':
            bln_evil = True

        elif option == '-h':
            usage()
            sys.exit(1)

        elif option == '-k':
            xor_key = binascii.a2b_hex(argument)

        elif option == '-o':
            output_file = argument

        elif option == '-s':
            sections = argument.split(',')
            


    # -----------------------------------------------------------
    # second pass: action
    # -----------------------------------------------------------
    for option, argument in optlist:
        # -----------------------------------------------------------
        # embedding a parasite file
        # -----------------------------------------------------------
        if option == '-m':
            host_file = args[0]
            parasite_file = args[1]
            host_file_contents = get_file(host_file)
            parasite_file_contents = get_file(parasite_file)

            # -----------------------------------------------------------
            # xor the parasite file if -k was specified
            # -----------------------------------------------------------
            if xor_key != 0:
                parasite_file_contents = xor_file(parasite_file_contents, xor_key)

            # -----------------------------------------------------------
            # chunks enumerated for file embedding
            # -----------------------------------------------------------
            if len(sections) > 0:           # -s specified
                evil_sections = evil_embed(host_file_contents, sections, parasite_file_contents, output_file)
                print '[+] Embedded %s into chunks: %s.' % (parasite_file, ','.join(evil_sections))
                print '[+] Final file is %s' % output_file

            # -----------------------------------------------------------
            # hack all the things
            # -----------------------------------------------------------
            elif bln_evil == True:          # -e specified
                host_free_sections = free_sections(host_file_contents, False)
                evil_sections = evil_embed(host_file_contents, host_free_sections, parasite_file_contents, output_file)
                print '[+] Embedded %s into chunks: %s.' % (parasite_file, ','.join(evil_sections))
                print '[+] Final file is %s' % output_file
                
            # -----------------------------------------------------------
            # automatically select a chunk for contiguous file embed
            # -----------------------------------------------------------
            elif bln_auto == True:          # -a specified
                host_free_sections = free_sections(host_file_contents, False)
                chunk_type = choice(host_free_sections)
                embed_file(host_file_contents, chunk_type, parasite_file_contents, output_file)
                print '[+] Embedded %s into chunk %s.' % (parasite_file, chunk_type, output_file)
                print '[+] Final file is %s' % output_file

            print '[+] Remember the chunks and their order for extraction.'
                

        # -----------------------------------------------------------
        # extracting a parasite file
        # -----------------------------------------------------------
        elif option == '-x':
            host_file = args[0]
            parasite_file = args[1]
            host_file_contents = get_file(host_file)

            if len(sections) == 0:
                print '[x] You must specify at least one chunk for extraction.'

            elif len(sections) == 1:
                chunk_type = ''.join(sections)
                extract_file(host_file_contents, chunk_type, parasite_file)
                print '[+] Extracted file to %s' % parasite_file

            else:
                s = 1


        # -----------------------------------------------------------
        # print listing of available chunks in host file
        # -----------------------------------------------------------
        elif option == '-p':
            host_file = args[0]
            host_file_contents = get_file(host_file)
            host_free_sections = free_sections(host_file_contents)




if __name__ == "__main__":
    main()
