import os, sys, re, argparse, glob
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

complement_nucleotide = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def reverse_complement(sequence):
    return ''.join(complement_nucleotide.get(base, base) for base in reversed(sequence.upper()))
def reverse_complement_fasta(input_file):
    # Read the input FASTA file
    sequence = read_fasta(input_file)
    return reverse_complement(sequence)
def read_fasta(file_path):
    sequences = {}
    current_id = None
    current_seq = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id:
        sequences[current_id] = ''.join(current_seq)

    return sequences[current_id]

def read_primer3_out(path):
    # read table 
    if os.path.exists(path):
        table = pd.read_csv(path, delim_whitespace=True, skiprows=3, header=None)
        columns = ['ID','sequence','pi','length','N','GC','Tm','self_any_th','self_end_th','harpin','quality']
        table.columns = columns
        return table
class main():
    def __init__(self):
        parser = argparse.ArgumentParser(
        description='Script create the roi sequence wiht the alternative allele instead the reference',
        usage='''utils.py <command> [<args>]
        Available commands are:
        get_alt_sequence
        build_primer3_input
        build_primer3_check_primer
        ''')

        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help(file=sys.stderr)
            exit(1)
        getattr(self, args.command)(sys.argv[2:])

    def get_alt_sequence(self, argv):
        parser = argparse.ArgumentParser(description='''
        Given a folder files .for and .rev are iteratively
        readed. Those files are expecting a name format:
        [qtn_name]_[amplicon-size]_[direction].[for|rev]
        ''')
        
        parser.add_argument("ref_fasta", help="Fasta file of reference sequence")
        parser.add_argument("ref_allele", help="Reference allele")
        parser.add_argument("alt_allele", help="Alternative allele")
        parser.add_argument("window_len", help="Window length")
        parser.add_argument("out", help="filename of output fasta")
        args = parser.parse_args(argv)
        
        sequences = SeqIO.to_dict(SeqIO.parse(args.ref_fasta, "fasta"))
        
        for seq_id, sequence in sequences.items():
            seq = str(sequence.seq)
            nseq = seq[:int(args.window_len)] + args.alt_allele + seq[int(args.window_len)+len(args.ref_allele):]
            sequence.seq = Seq(nseq)

        SeqIO.write(sequences.values(), args.out, "fasta")
    def build_primer3_input(self, argv):
        parser = argparse.ArgumentParser(description='''
        Get primer3 input file using the input DNA sequence
        ''')
        parser.add_argument("name", help="name sequence")
        parser.add_argument("fasta", help="Fasta file of roi")
        parser.add_argument("target_pos", help="position of target allele")
        parser.add_argument("allele", help="target allele")
        parser.add_argument("orientation", help="convert the input to reverse complement [53|35]")
        parser.add_argument("sep", help="distance from the allele to search common")
        parser.add_argument("out", help="output file")
        args = parser.parse_args(argv)

        if args.orientation == '35':
            seq = reverse_complement_fasta(args.fasta)
            position = len(seq) - int(args.target_pos )
        else:
            seq = read_fasta(args.fasta)
            position = int(args.target_pos) + len(args.allele)
            
        with open(args.out, 'w') as file:
            file.write("SEQUENCE_ID={name}\n".format(name=args.name))
            file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
            file.write("SEQUENCE_TARGET={pos},{sep}\n".format(
                pos=position, sep=args.sep))
            file.write("P3_FILE_FLAG=1\n=".format(name=args.name))

    def build_primer3_check_primer(self, argv):
            parser = argparse.ArgumentParser(description='''
            Get primer3 input file using the input DNA sequence
            ''')
            parser.add_argument("name", help="name sequence")
            parser.add_argument("primer3_out", help="path of reference allele-specific primer3 out")
            parser.add_argument("fasta", help="Fasta file of roi")
            parser.add_argument("orientation", help="convert the input to reverse complement [53|35]")
            parser.add_argument("ref", help="ref allele")
            parser.add_argument("alt", help="alt allele")
            parser.add_argument("sep", help="distance from the allele to search common")
            parser.add_argument("out", help="output dir")
            args = parser.parse_args(argv)

            primers = read_primer3_out(args.primer3_out)

            for row_id, primer in primers.iterrows():
                name = "{base_name}_{primer_id}".format(base_name = args.name, primer_id =  primer.ID)
                
                if args.orientation == '35':
                    seq = reverse_complement_fasta(args.fasta)
                    name += '_35'
                    
                    complement_allele = ''.join([complement_nucleotide [nucleotide] for nucleotide in args.alt])
                    primer_seq = primer.sequence[:-len(args.ref)]
                    primer_seq += complement_allele[::-1]
                    
                else:
                    seq = read_fasta(args.fasta)
                    name += '_53'
                    primer_seq = primer.sequence[:-len(args.ref)]
                    primer_seq += args.alt

                
                if not os.path.exists(args.out):
                    os.makedirs(args.out)
                
                
                input_name = '{variant_id}_{allele}_{primer_id}.txt'.format(variant_id = args.name.split('_')[0],
                                                                    allele = "alt",
                                                                    primer_id = primer.ID)
                outpath = os.path.join(args.out, input_name)
                
                with open(outpath, 'w') as file:
                    file.write("SEQUENCE_ID={name}\n".format(name=name))
                    file.write("SEQUENCE_TEMPLATE={seq}\n".format(seq=seq))
                    file.write("SEQUENCE_PRIMER={primer_seq}\n".format(primer_seq=primer_seq))
                    file.write("P3_FILE_FLAG=1\n")
                    file.write("PRIMER_PICK_RIGHT_PRIMER=0\n")
                    file.write("PRIMER_PICK_LEFT_PRIMER=1\n=")



if __name__ == "__main__":
    main()