import os, sys, re, argparse, glob
import pandas as pd

# Util functions
def read_primer3_out(path):
    # read table 
    if os.path.exists(path):
        table = pd.read_csv(path, sep='\s+', skiprows=3, header=None)
        columns = ['ID','sequence','pi','length','N','GC','Tm','self_any_th','self_end_th','harpin','quality']
        table.columns = columns
        return table

def filter_primers_tm(row):
    # Filter allele-specific primers
    if row.direction == '53':
        if (row.type == 'for') and (58.3 <= row.Tm <= 63.5):
            return True
        elif (row.type == 'rev') and (62 <= row.Tm <= 75):
            return True
        else:
            return False
    else:
        if (row.type == 'rev') and (58.3 <= row.Tm <= 63.5):
            return True
        elif (row.type == 'for') and (62 <= row.Tm <= 75):
            return True
        else:
            return False
        
def get_suitable_combinations(df):
    suitable_combs = list()
    for direction, sdf in df.groupby('direction'):
        cross = pd.crosstab(sdf['amplicon_size'], sdf['type'])
        if cross.shape[1] == 2:
            for amplicon_l, row in cross.iterrows():
                if row['for'] == 1 and row['rev'] == 1:
                    suitable_combs.append([direction, amplicon_l])
    return suitable_combs

def output_suitable_combination(comb, df):
    output = dict()
    suitable_primers = df[(df['direction'] == comb[0]) & (df['amplicon_size'] == comb[1])]
    for primer_type, sdf in suitable_primers.groupby('type'):
        output[primer_type] = sdf
    return output
        
    
# Interface


class main():
    def __init__(self):
        parser = argparse.ArgumentParser(
        description='Script to select primer-sets according to Tm-shift ranges of Tm',
        usage='''select_primers_tm.py <command> [<args>]
        Available commands are:
        select      Select the primer sets (forwards and  reverses)
                    that meets the ranges of Want etal.,(2005):
                        * Tm allele-specific 54.3-63.5ºC
                        * Tm common primer 62-75ºC
                    Select iteratively as imput primer3 for and rev files.
        ''')

        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help(file=sys.stderr)
            exit(1)
        getattr(self, args.command)(sys.argv[2:])

    def select(self, argv):
        parser = argparse.ArgumentParser(description='''
        Given a folder files .for and .rev are iteratively
        readed. Those files are expecting a name format:
        [qtn_name]_[amplicon-size]_[direction].[for|rev]
        ''')
        
        parser.add_argument("dir", help="input dir of primer3 outputs")
        parser.add_argument("outdir", help="output directory")
        parser.add_argument("prefix", help="prefix of outputs")
        args = parser.parse_args(argv)
        
        dfs = list()
        for file in glob.glob(args.dir + '/*'):
            # Is primer3 output?
            if file[-3:] in ['rev', 'for']:
                name = file.split('/')[-1][:-4]
                params = name.split('_')
                
                if params == ['1', '53'] and file[-3:] == 'for':
                    df = read_primer3_out(file)
                    df['type'] = file[-3:]
                    df['amplicon_size'] = params[1]
                    df['direction'] = params[2]
                    
                
                

                dfs.append(df)
        compiled = pd.concat(dfs, ignore_index=True)
        
        compiled['Tm-shift-ranges'] = compiled.apply(filter_primers_tm, axis=1)
        
        counts_pass = compiled[compiled['Tm-shift-ranges']].groupby(['direction', 'type', 'amplicon_size'], as_index=False).count()
        
        pass_combs = get_suitable_combinations(counts_pass)
        
        rows = list()
        for comb in pass_combs:
            primers = output_suitable_combination(comb, compiled[compiled['Tm-shift-ranges']])
            primers['for'].to_csv('{outdir}/{prefix}_{ampsize}_{direction}_forward.csv'.format(
                outdir = args.outdir,
                prefix = args.prefix,
                ampsize = comb[1],
                direction = comb[0]), sep='\t', index=False)
            primers['rev'].to_csv('{outdir}/{prefix}_{ampsize}_{direction}_reverse.csv'.format(
                outdir = args.outdir,
                prefix = args.prefix,
                ampsize = comb[1],
                direction = comb[0]), sep='\t', index=False)
            
            
            row = [comb[0], comb[1], primers['for'].shape[0], primers['rev'].shape[0]]
            
            rows.append(row)
            
        summary = pd.DataFrame(rows, columns=['Orientation', 'Amplicon_size', 'No_A', 'No_B'])
        summary = summary.sort_values(['Orientation', 'Amplicon_size'])
        summary.to_csv('{outdir}/{prefix}_summary.csv'.format(outdir = args.outdir,
                    prefix = args.prefix), sep='\t', index=False)

        
if __name__ == "__main__":
    main()