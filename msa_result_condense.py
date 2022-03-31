# parse kalign alignment results and turn into a tab-delimited file
# that can be loaded into a spreadsheet
# this might also work for output from other multiple-sequence aligner provided in multi-fasta format
#
# Karsten Hokamp, 2022 (kahokamp@tcd.ie)
#
# USAGE:
# kalign_results_parser.py multi_fasta_alignment_file [reference_multi_fasta_file]
#
# the optional file of reference sequences allows to split the output into smaller chunks
# for this the order of the reference sequences must be in the same order as file used for the alignment
# output files are named after the sequence id

import sys
import os

# when reference is provided an aligner needs to be called
# instead of the default (below) an aligner can be provided as third argument
aligner = 'kalign'

# get name of input file from the command line
args = sys.argv[1:]
if len(args) < 1:
    raise Exception('Please provide the name of a file containing MSA output!')
    
file = args.pop(0)

# if a second argument is given, read in set of Fasta sequences 
# that will be used  as reference for splitting the alignment
reference = ''
if args :
    reference = args.pop(0)

# an additional argument can be used to specify an aligner tool (kalign is the default)
if args :
    aligner = args.pop(0)
# make sure aligner is found
if reference :
    if not aligner.endswith('kalign') :
        print('Script is set up to run with kalign - for a different aligner, please adjust the job specification further down')
    if not os.path.isfile(aligner) :
        aligner_location = os.popen(f'which {aligner}').read().strip()
        if not os.path.isfile(aligner_location) :
            raise Exception(f'The specified aligner ({aligner}) does not seem to be available on your system.')
        else :
            print(f'aligner found in {aligner_location}')
                            
ref = dict()
ref_order = list()
ref_all = ''
if reference :
    with open(reference, 'r') as f:
        for line in f:
            if line.startswith('>') :
                header = line.strip().replace('>', '')
                ref_order.append(header)
                ref[header] = 0
            else :
                ref[header] += len(line.strip())
                ref_all += line.strip().upper()
                
# read in MSA:
with open(file, 'r') as fh:
    print(f'reading alignments from {file}...')

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
            seqs[header] += line.strip().upper()
        
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

print(f"\nwriting output into '{outfile}'...\n")
# write parsed info into output file, one line per alignment position
fout.write('position\tconsensus\tabsolute\tpercent\tgaps\tothers\n')
msa = ''
msa_out = dict()
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
    pos = i + 1
    fout.write(f'{pos}\t{top_residue}\t{top_num}\t{top_perc}\t{gaps}\t{others_out}\n')
    msa_out[i] = f'{pos}\t{top_residue}\t{top_num}\t{top_perc}\t{gaps}\t{others_out}\n'
    msa += top_residue

fout.close()

if reference :
    
    # split output into chunks according to sequences in reference file
    
    outfile = f'{file}.{reference}.all'
    with open(outfile, 'w') as fref :
        fref.write(f'>msa\n{msa}\n>ref\n{ref_all}')
    
    # run kalign on reference and msa consensus to find gaps
    msa_ref_out = outfile + '.msa'
    job = f'{aligner} -i "{outfile}" -o "{msa_ref_out}"'        
    print(f'running alignment job: {job}')
    msa_err = os.popen(job).read()
    if not os.path.isfile(msa_ref_out) :
        raise Exception('Something went wrong - the output of the alignment (msa consensus vs reference) could not be found in {msa_ref_out}')
        
    msa_kalign = dict()
    with open(msa_ref_out) as fkal :    
        for line in fkal :
            if line.startswith('>') :
                header = line.strip().replace('>', '')
                msa_kalign[header] = ''
            else :
                msa_kalign[header] += line.strip()
    
    if msa_kalign["msa"].find('-') > -1 :
        print(f'\nreference has extra residues (e.g. at {msa_kalign["msa"].find("-")})')
    if msa_kalign["ref"].find('-') > -1 :
        print(f'\nMSA has extra residues (e.g. at {msa_kalign["ref"].find("-")})')
    ref_ext = list(msa_kalign['ref'])
    msa_ext = list(msa_kalign['msa'])
    
    print('\nsplitting alignment...')
    
    # for each reference sequence, go through each residue in the alignment and adjust length for gaps
    i = 0
    start = 0
    for seqid in ref_order :
        msa_sub = ''
        ref_sub = ''
        reflen = ref[seqid]
        print(f'\n{i}: start of {seqid} ({reflen} bp)')
        outfile = f'{file}.{seqid}.xls'
        fout_part = open(outfile, 'w')
        fout_part.write('pos_ref\treference\tpos_msa\tconsensus\tabsolute\tpercent\tgaps\tothers\n')
        ref_sub_len = 0
        gaps_ref = 0
        msa_sub_len = 0
        gaps_msa = 0
        
        end = start + reflen
        pos = 0
        while pos < reflen :

            add_ref = ref_ext.pop(0)
            add_msa = msa_ext.pop(0)
            msa_sub += add_msa
            ref_sub += add_ref            
            
            while add_ref == '-':
#                pos += 1
                gaps_ref += 1
                fout_part.write(f'{pos}\t{add_ref}\t{msa_out[i]}')
                i += 1
                add_ref = ref_ext.pop(0)
                add_msa = msa_ext.pop(0)
                msa_sub += add_msa
                ref_sub += add_ref
            
            while add_msa == '-':
                pos += 1
                fout_part.write(f'{pos}\t{add_ref}\t-\n')
                gaps_msa += 1
#                msa_sub += '-'
                if pos >= reflen :
                    break
                add_ref = ref_ext.pop(0)
                add_msa = msa_ext.pop(0)
                msa_sub += add_msa
                ref_sub += add_ref

            if pos >= reflen :
                break

            pos += 1
            fout_part.write(f'{pos}\t{add_ref}\t{msa_out[i]}')
            i += 1
            

        fout_part.close()
        start = end
                
        if msa_sub != ref_sub :
            print(f'differences between ref and msa for {seqid} (gaps in ref: {gaps_ref}, msa: {gaps_msa})\nMSA: {msa_sub}\nREF: {ref_sub}')

print('\nFinished!')
