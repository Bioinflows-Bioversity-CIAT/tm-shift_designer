import os, sys, re, argparse, glob
import pandas as pd
import re 
from collections import defaultdict

def read_primer3_out(path):
    # read table 
    if os.path.exists(path):
        table = pd.read_csv(path, sep=r'\s+', skiprows=3, header=None)
        columns = ['ID','sequence','pi','length','N','GC','Tm','self_any_th','self_end_th','harpin','quality']
        table.columns = columns
        return table



class Primer:
    def __init__(self, qtl_name, workdir):
        primer_id = None
        sequence = None
        start_pos = None
        length = None
        GC_content = None
        tm = None
        self_any = None
        self_end = None
        harpin = None
        quality = None
        primer_type = None
        orientation = None

        
    def load(self, row, primer_type, orientation):
        self.primer_id = row.ID
        self.sequence = row.sequence
        self.start_pos = row.pi
        self.length = row.length
        self.GC_content = row.GC
        self.tm = row.Tm
        self.self_any = row.self_any_th
        self.self_end = row.self_end_th
        self.harpin = row.harpin
        self.quality = row.quality
        self.primer_type = primer_type
        self.orientation = orientation
        
    def temperature_filter(self):
        
        if (self.primer_type == 'for') and (58.3 <= self.tm <= 63.5):
            return True
        elif (self.primer_type == 'rev') and (62 <= self.tm <= 75):
            return True
        else:
            return False
        
            
class Common(Primer):
    def __init__(self, qtl_name):
        self.primer_type = 'common'
        self.dist = None
        self.qtl_name = qtl_name
        
class Allele_specific(Primer):
    def __init__(self, qtl_name, workdir):
        self.primer_type = 'allele-specific'
        
        
        self.workdir = workdir
        self.qtl_name = qtl_name
        
        allele_ref= None
        allele_alt= None
    
    def get_alternative_path(self,):
        folder = '{workdir}/alt_allele/{qtl_name}_alt_1_{orient}_for_out/'.format(
            workdir = self.workdir,
            qtl_name=self.qtl_name,
            orient = self.orientation)
        
        file = '{qtl_name}_alt_1_for_{primer_id}_{orient}.for'.format(
            qtl_name=self.qtl_name,
            orient = self.orientation,
            primer_id=self.primer_id)
        
        path = os.path.join(folder, file)
        return path
    
    def read_alternative(self):
        path = self.get_alternative_path()
        
        try:
            primer_df = read_primer3_out(path)
        except FileNotFoundError:
            print('File not found!')
            
        for n, row in primer_df.iterrows():
            primer = Allele_specific(self.qtl_name, self.workdir)
            primer.load(row, self.primer_type, self.orientation)
            
        self.alternative = primer
        
class Allele_specific_space():
    def __init__(self, space):
        self.space = space
        self.tables = dict()
        self.columns = ['qtn_name','primer_id','sequence_ref','start_pos_allele',
                        'length_ref','GC_content_ref','tm_ref','self_any_ref',
                        'self_end_ref','harpin_ref','quality_ref',
                       'sequence_alt','length_alt',
                        'GC_content_alt','tm_alt',
                        'self_any_alt','self_end_alt',
                        'harpin_alt','quality_alt',
                       'pr_ref_name', 'pr_alt_name']
        self.create_space()
        
        
    def create_space(self,):
        primer_options = defaultdict(list)
        for orientation, values in self.space.items():
            count = 1
            for primer in values:
                primer.read_alternative()

                # Add the GC-rich tail according to temperature
                if primer.tm >= primer.alternative.tm:
                    primer.sequence = "gcgggcagggcggc" + primer.sequence
                    primer.alternative.sequence = "gcgggc" + primer.alternative.sequence
                else:
                    primer.sequence = "gcgggc" + primer.sequence
                    primer.alternative.sequence = "gcgggcagggcggc" + primer.alternative.sequence

                data = [primer.qtl_name, primer.primer_id,primer.sequence,primer.start_pos,
                        primer.length,primer.GC_content,float(primer.tm),
                        primer.self_any,primer.self_end,
                        primer.harpin,primer.quality,
                       primer.alternative.sequence,primer.alternative.length,
                        primer.alternative.GC_content,float(primer.alternative.tm),
                        primer.alternative.self_any,primer.alternative.self_end,
                        primer.alternative.harpin,primer.alternative.quality]
                
                zero_counter = str(count).zfill(2)

                if orientation == '53':
                    or_sfx = "F"
                else:
                    or_sfx = "R"
                    
                primer_name_ref = f"{primer.qtl_name}_{zero_counter}_{or_sfx}a"
                primer_name_alt = f"{primer.qtl_name}_{zero_counter}_{or_sfx}b"
                
                data.extend([primer_name_ref,primer_name_alt])
                primer_options[orientation].append(data)
                count += 1

        for orientation, values in primer_options.items():
            df = pd.DataFrame(values)
            df.columns = self.columns
            self.tables[orientation] = df
            
class Common_space():
    
    def __init__(self, space):
        self.space = space
        self.tables = dict()
        self.columns = ['qtn_name','primer_id','sequence_common','start_pos_common',
                        'length_common','GC_content_common','tm_common', 'distance_common','self_any_common',
                        'self_end_common','harpin_common','quality_common','orientation_common',
                        'pr_common_name']
        self.create_space()
        
        
    def create_space(self,):
        primer_options = defaultdict(list)
        for orientation, values in self.space.items():
            count = 1
            for primer in values:
                data = [primer.qtl_name, primer.primer_id,primer.sequence,primer.start_pos,
                        primer.length,primer.GC_content,float(primer.tm), int(primer.dist),
                        primer.self_any,primer.self_end,
                        primer.harpin,primer.quality, primer.orientation]
                
                
                zero_counter = str(count).zfill(2)

                if orientation == '53':
                    or_sfx = "R"
                else:
                    or_sfx = "F"
                    
                primer_name = f"{primer.qtl_name}_{or_sfx}{zero_counter}"

                
                data.extend([primer_name])
                primer_options[orientation].append(data)
                count += 1

        for orientation, values in primer_options.items():
            df = pd.DataFrame(values)
            df.columns = self.columns
            self.tables[orientation] = df
        
    
# Interface


class main():
    def __init__(self):
        parser = argparse.ArgumentParser(
        description='Script to select primer-sets according to Tm-shift ranges of Tm',
        usage='''select_primers_tm.py <command> [<args>]
        Available commands are:
        summarize      Select the primer sets (forwards and  reverses)
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

    def summarize(self, argv):
        parser = argparse.ArgumentParser(description='''
        Given a folder files .for and .rev are iteratively
        readed. Those files are expecting a name format:
        [qtn_name]_[amplicon-size]_[direction].[for|rev]
        ''')
        
        parser.add_argument("qtl_name", help="name of the qtl")
        parser.add_argument("basedir", help="path of working directory")
        parser.add_argument("outdir", help="path of output directory")
        args = parser.parse_args(argv)
        
        
        # Read allele-specific output read
        allele_53_path =os.path.join(args.basedir, f'ref_allele/{args.qtl_name}_ref_1_53.for')
        allele_35_path = os.path.join(args.basedir, f'ref_allele/{args.qtl_name}_ref_1_35.for')
        allele_paths = [allele_53_path, allele_35_path]
        allele_specific = defaultdict(list)
        allele_specific_fail = defaultdict(list)
        for output, orientation in zip(allele_paths, ['53', '35']):
            for n, row in read_primer3_out(output).iterrows():
                primer_type = output[-3:]
                primer = Allele_specific(args.qtl_name, args.basedir)
                primer.load(row, primer_type, orientation)
                if primer.temperature_filter():
                    allele_specific[orientation].append(primer)
                else:
                    allele_specific_fail[orientation].append(primer)
                    
                    
        # Common primer output read
        common = defaultdict(list)
        common_fail = defaultdict(list)
        for file in glob.glob(args.basedir + f"/ref_allele/{args.qtl_name}*"):
            if file[-3:] in ['rev']:
                primer_type = file[-3:]
                name = file.split('/')[-1][:-4]
                params = name.split('_')
            
                orientation = params[4]
                for n, row in read_primer3_out(file).iterrows():
                    primer = Common(args.qtl_name)
                    primer.load(row, primer_type, orientation)
                    primer.dist = params[3]
                    if primer.temperature_filter():
                        common[orientation].append(primer)
                    else:
                        common_fail[orientation].append(primer)
        
        # Remove orientations if incomplete
        orientations = list(set(allele_specific.keys()) & set(common.keys()))
        remove_orientation = [ori for ori in ["53","35"] if not ori in orientations]
        for ori in remove_orientation:
            if ori in allele_specific:
                del allele_specific[ori]
            if ori in common:
                del common[ori]
        
        
        options = Allele_specific_space(allele_specific)
        options_common = Common_space(common)
        merged_list = list()
        stacked_list = list()
        if len(list(options.tables.keys())) > 0:
            for orientation in options.tables.keys():
                merged_table = pd.merge(options.tables[orientation],
                                        options_common.tables[orientation],
                                        on = 'qtn_name', suffixes=('_allele', '_common'))
                merged_list.append(merged_table)



                stack_columns = ['pr_name', 'sequence', 'length', 'tm', 'GC_content']
                ref_columns = ['pr_ref_name', 'sequence_ref', 'length_ref', 'tm_ref', 'GC_content_ref']
                ref_df = options.tables[orientation][ref_columns]
                ref_df.columns = stack_columns 

                alt_columns = ['pr_alt_name', 'sequence_alt', 'length_alt', 'tm_alt', 'GC_content_alt']
                alt_df = options.tables[orientation][alt_columns]
                alt_df.columns = stack_columns 

                common_columns = ['pr_common_name', 'sequence_common', 'length_common', 'tm_common', 'GC_content_common']
                common_df = options_common.tables[orientation][common_columns]
                common_df.columns = stack_columns 

                stacked_list.extend([ref_df,alt_df,common_df])


            merged = pd.concat(merged_list, ignore_index=True)
            merged['amplicon_size_ref'] = (merged['start_pos_common'] + merged['length_common'] ) - merged['start_pos_allele']  
            merged['amplicon_size_alt'] = (merged['start_pos_common'] + merged['length_common'] ) - merged['start_pos_allele']  
            order_columns = ['qtn_name','orientation_common',
                             'pr_ref_name','pr_alt_name','pr_common_name',
                             'tm_ref','tm_alt', 'tm_common',
                             'GC_content_ref','GC_content_alt','GC_content_common',
                             'self_any_ref','self_any_alt','self_any_common',
                             'self_end_ref','self_end_alt','self_end_common',
                             'harpin_ref','harpin_alt','harpin_common',
                             'quality_ref','quality_alt','quality_common',
                             'sequence_ref','sequence_alt','sequence_common',
                             'distance_common','amplicon_size_ref','amplicon_size_alt']
            out_merged = merged[order_columns]
            out_merged['Tm_diff'] = abs(out_merged['tm_ref'] - out_merged['tm_alt'])
            out_merged = out_merged.sort_values(['Tm_diff', 'quality_ref', 'quality_common'])
            stacked_df = pd.concat(stacked_list, ignore_index=True)

            os.mkdir(args.outdir)
            merged_path = f'{args.outdir}/{args.qtl_name}_combinations.csv'
            out_merged.to_csv(merged_path, index = False)

            stacked_path = f'{args.outdir}/{args.qtl_name}_primers.csv'
            stacked_df.to_csv(stacked_path, index = False)
        else:
            os.mkdir(args.outdir)

        
if __name__ == "__main__":
    main()