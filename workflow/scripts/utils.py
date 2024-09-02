import os, sys, re, argparse, glob
from Bio import SeqIO
from Bio.Seq import Seq



class main():
    def __init__(self):
        parser = argparse.ArgumentParser(
        description='Script create the roi sequence wiht the alternative allele instead the reference',
        usage='''get_alternative_roi.py <command> [<args>]
        Available commands are:
        get_sequence 
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

if __name__ == "__main__":
    main()