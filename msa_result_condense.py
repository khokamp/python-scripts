# If you have results from a multiple sequence aligner (MSA) that are too big to view, i.e. thousands of files
# then you can use this script to turn this alignment (in multi-fasta format) into a condensed format which can
# be loaded into a spreadsheet program.
# For each position in the alignment you will get a row in the table which indicates the most common residue
# at that position, the number of gaps and the alternative residues (in decreasing order).
# This script was initially developed for output from kalign but should work for any MSA results in multi-fasta format.
#
# Karsten Hokamp, March 2022 (kahokamp@tcd.ie)
#
# USAGE:
# python3 msa_result_condense.py multi_fasta_alignment_file
#
# output will be stored in multi_fasta_alignment_file.xls

import sys

# get name of input file from the command line
file = sys.argv[1]

# try to open it
fh = open(file, 'r')
print(f'Reading alignments from {file}...')

# use input file as root for output file
outfile = f'{file}.xls'
fout = open(outfile, 'w')


align = dict()
seqs = dict()
lens = dict()
total_seqs = 0

# read in all alignments
# and store in dict seqs based on header
headers = 0
for line in fh :
    if line.startswith('>') :
        header = line.strip().replace('>', '')
        headers += 1
        seqs[header] = ''
        total_seqs += 1
    else :
        seqs[header] += line.strip()

fh.close()
print(f'read {headers} sequences\nstarting the parse process...')

# go through each sequence and organise by position
for seqid in seqs :
    seq = seqs[seqid]
    seqlen = len(seq)

    # keep track of sequence lengths
    if seqlen in lens.keys() :
        lens[seqlen] += 1
    else :
        lens[seqlen] = 1

    # collect residues and their occurrences at each position of the alignment
    for i in range(0, seqlen) :
        residue = seq[i]
        if i in align.keys() :
            if residue in align[i].keys() :
                align[i][residue] += 1
            else :
                align[i][residue] = 1
        else :
            align[i] = dict()
            align[i][residue] = 1

# print summary info on alignment lengths
print('\nalignment lengths detected:\nlength\tsequences')
for alen in sorted(lens.keys()) :
    print(f'{alen}\t{lens[alen]}')

print(f"\nwriting output into '{outfile}'...")
# write parsed info into output file, one line per alignment position
fout.write('position\tmajor\tabsolute\tpercent\tgaps\tothers\n')
for i in sorted(align.keys()) :
    freq = dict()
    gaps = 0
    for residue in align[i].keys() :
        if residue == '-' :
            gaps = align[i][residue]
            continue
            
        num = align[i][residue]
        if not num in freq.keys() :
            freq[num] = list()
            
        freq[num].append(residue)

    # separate major residue from the rest
    num_order = sorted(freq.keys(), reverse=True)
    top_num = num_order.pop(0)
    top_residue = ','.join(freq[top_num])
    top_perc = top_num / total_seqs * 100
    others = list()

    # list alternative residues by decreasing occurrence
    for num in num_order :
        other_residue = ','.join(freq[num])
        others.append(f'{other_residue} ({num})')
        
    others_out = '; '.join(others)
    fout.write(f'{i}\t{top_residue}\t{top_num}\t{top_perc}\t{gaps}\t{others_out}\n')

fout.close()

print('Finished!')
