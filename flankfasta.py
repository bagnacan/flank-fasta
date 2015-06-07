# Quick and dirty script to transform FASTA files in FASTA files
# whose sequence alignment is optimized for performing Phylogenetic
# Analysis via calculation of Maximum Likelihood trees.
# This optimizer needs:
# - Tcoffee (https://github.com/cbcrg/tcoffee) configured with
#   Mafft, Probcons, Muscle;
# - Gblocks (http://molevol.cmima.csic.es/castresana/Gblocks.html)
#

import codecs
import getopt
import subprocess
import sys
import os
from os.path import isfile, join, abspath, dirname, basename

EMAIL   = 'email'
AMINO   = 'amino'
NUCLEIC = 'nucleic'
aln_ext = '.fasta_aln'
aln_ext__optimized = '_gb.txt'
aln_ext__flanked   = '_flanked.fasta'
script__gblocks_alignment_optimization = 'biocomp.gblocks'
GBLOCKS_FLANKS = "Flanks:"
GBLOCKS_FLANK_START = "["
GBLOCKS_FLANK_END   = "]"

# usage
def print_usage():
    '''
    Prints usage
    '''
    print 'u suc'



# read a FASTA file
#
def read_fasta(file_name):
    '''
    Reads the FASTA file at the given path and returns a dictionary
    representing the file's internal structure
    '''
    fasta = {}
    key = ''

    for line in open(file_name):
        line = line.rstrip()
        
        if line.startswith('>'):
            key = line.replace('>', '')
            fasta[key] = ''
        else:
            fasta[key] += line

    return fasta



# write a FASTA file
#
def write_fasta(fasta_dict, path):
    '''
    Writes a FASTA file from the given dictionary representation at
    the given path
    '''
    fasta = open(path, 'w')
    for key in list(fasta_dict.keys()):
        fasta.write('>' + key + '\n')
        fasta.write(fasta_dict[key] + '\n')
    fasta.close()



# return the flanks of the alignment
#
def get_flanks(aln_file):
    '''
    Returns the range positions (flanks) that optimize the alignment
    contained within the given file
    '''
    result = []

    for line in open(aln_file):
        line = line.rstrip()
        if line.startswith(GBLOCKS_FLANKS):
            line = line.replace(GBLOCKS_FLANKS, '')
            line = line.split()
            i = 0
            while i < len(line):
                line_range = line[i:(i+2)]
                min_range = int(line_range[0].replace(GBLOCKS_FLANK_START, ''))
                max_range = int(line_range[1].replace(GBLOCKS_FLANK_END, ''))
                flank = (min_range, max_range)
                result.append(flank)
                i += 2
    return result



# read a FASTA dictionary and return a flanked FASTA dictionary
#
def read_flanked_fasta(fasta, flanks):
    '''
    Returns the given FASTA dictionary that is comprised within the
    given flank regions
    '''
    result = {}

    for key in list(fasta.keys()):
        sequence = fasta[key]
        result[key] = ''
        for flank in flanks:
            result[key] += sequence[(flank[0]-1):flank[1]]

    return result



# get command line
#
def get_command_line(argv):
    '''
    Returns the configuraion from the given command line arguments:
    -a/--amino    directory whose files contain amino acid sequences
    -n/--nucleic  directory whose files contain nucleic acid sequences
    -e/--email    the user's e-mail
    '''
    email = ''
    nucleic_files = []
    amino_files   = []
    config = {}
    sopt_mail = 'e'
    sopt_amin = 'a'
    sopt_nucl = 'n'
    lopt_mail = 'email'
    lopt_amin = 'amino'
    lopt_nucl = 'nucleic'

    try:
        opts, args = getopt.getopt(
                argv,
                sopt_mail + ':' + sopt_amin + ':' + sopt_nucl + ':',
                [str('--' + lopt_mail), str('--' + lopt_amin), str('--' + lopt_nucl)])
    except getopt.GetoptError:
        print_usage()
        sys.exit()

    for opt, arg in opts:
        if opt in (str('-' + sopt_nucl), lopt_nucl):
            nucleic_files = [join(arg, i) for i in os.listdir(arg) if isfile(join(arg, i))]
        if opt in (str('-' + sopt_amin), lopt_amin):
            amino_files = [join(arg, i) for i in os.listdir(arg) if isfile(join(arg, i))]
        if opt in (str('-' + sopt_mail), lopt_mail):
            email = arg

    if ((len(amino_files) + len(nucleic_files)) == 0) or (email == ''):
        print_usage()
        sys.exit()

    config[NUCLEIC] = nucleic_files
    config[AMINO]   = amino_files
    config[EMAIL]   = email
    return config



# optimize the protein/nucleic-acid sequence alignment
#
def optimize_alignment(alignment_file, sequence_type):
    '''
    Writes a FASTA file containing the optimized version (Gblocks)
    of the given protein/nucleic-acid sequence alignment
    '''
    # launch Gblocks differentiating between protein sequences and
    # nucleic sequences
    if sequence_type == AMINO:
        gblocks_call = [script__gblocks_alignment_optimization,
                alignment_file, '-t=p', '-e=_gb', '-p=t', '-b4=5']
    else:
        gblocks_call = [script__gblocks_alignment_optimization,
                alignment_file, '-t=d', '-e=_gb', '-p=t', '-b4=5']
        
    fd = read_fasta(alignment_file)

    # Gblocks optimization
    subprocess.call(gblocks_call)

    # restrict the FASTA alignment with the Gblocks flanked positions
    optimized_alignment_file = alignment_file + aln_ext__optimized
    flanks = get_flanks(optimized_alignment_file)
    ofd = read_flanked_fasta(fd, flanks)

    # write the Gblocks optimized fasta alignment
    flanked_file = alignment_file + aln_ext__flanked
    write_fasta(ofd, flanked_file)


    
# Tcoffee sequence alignment
#
def sequence_alignment(in_file):
    '''
    Performs sequence alignment of the sequences in the given fasta
    file using Tcoffee sequence alignment methods.
    '''
    parent_dir = os.getcwd()
    os.chdir(dirname(abspath(in_file)))
    in_file = basename(in_file)

    tcoffee_seqmet  = '-in=S' + in_file + ',Mmafft_msa,Mmuscle_msa,Mprobcons_msa'
    tcoffee_cpus    = '2'
    tcoffee_outfmt  = '-output=fasta'
    tcoffee_outfile = '-outfile=' + str(in_file.replace('.fas', '.fasta_aln'))

    # Tcoffee alignment
    print '\n# Tcoffee sequence alignment on file \"' + in_file + '\"'
    print '#'
    subprocess.call(['t_coffee', tcoffee_seqmet, tcoffee_outfile,
        '-cpu', tcoffee_cpus, tcoffee_outfmt])


    # read the fasta alignment as a dictionary
    alignment_file = in_file.replace('.fas', aln_ext)
    optimize_alignment(alignment_file, NUCLEIC)

    os.chdir(parent_dir)



# Tcoffee structural alignment
#
def structural_alignment(in_file):
    '''
    Performs structural alignment of the sequences in the given fasta
    file using Tcoffee structural alignment utility.
    '''
    parent_dir = os.getcwd()
    os.chdir(dirname(abspath(in_file)))
    in_file = basename(in_file)

    # Tcoffee alignment
    print '\n# Tcoffee structural alignment on file \"' + in_file + '\"'
    print '#'
    
    # structural alignment step 1:
    #
    # Make a multiple profile alignment of all the sequences in the
    # given file using Tcoffee Psi-BLAST
    subprocess.call(['t_coffee', str('-in=S' + in_file),
        '-mode psicoffee', str('-email ' + config[EMAIL]),
        '-multi_core', '-output=fasta_aln'])

    # structural alignment step 2:
    #
    # Fetch candidate templates via BLAST against the PDB, and use
    # the candidate's structure to obtain the sequence alignment
    # out of the structural alignment.
    # In case of missing template is missing, proba_pair (the
    # default T-Coffee method for building libraries) is used
    # instead of SAP.
    # By default Expresso fetches only XRAY structures, but it is
    # possible to control this behavior via the flag -pdb_type=dn
    # where d for is for diffraction (X-Ray) and n for NMR.
    subprocess.call(['t_coffee', str('-in=S' + in_file),
        '-mode expresso', '-pdb_type=dn', str('-email ' + config[EMAIL]),
        '-multi_core', '-output=fasta_aln'])

    # structural alignment step 3:
    #
    # Refinement of the previous methods
    subprocess.call(['t_coffee', str('-in=S' + in_file),
        '-mode accurate', str('-email ' + config[EMAIL]),
        '-multi_core', '-output=fasta_aln'])


    # read the fasta alignment as a dictionary
    alignment_file = in_file.replace('.fas', aln_ext)
    optimize_alignment(alignment_file, AMINO)

    os.chdir(parent_dir)



if __name__ == '__main__':
    reload(sys)
    
    # get the list of nucleic/aminoacid sequence files
    cli = get_command_line(sys.argv[1:])

    # nucleic acid sequences from non-coding DNA strands are aligned
    # using Tcoffee and optimized using Gblocks
    for in_file in cli[NUCLEIC]:
        sequence_alignment(in_file)

    # amino acid sequences from coding DNA strand are aligned using
    # Tcoffee template-extrapolated structural alignment and
    # optimized using Gblocks
    for in_file in cli[AMINO]:
        structural_alignment(in_file)

