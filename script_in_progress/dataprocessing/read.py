import os
import sys
import pandas as pd
import numpy as np
from .utility import empty_dataframe, pn, po

class Readfiles():
    def __init__(self,files):
        """
            Prepare list of input files.
        """

        self.files = []
        for f in files:
            if f.endswith('.txt'):
                self.files.append(f)
            elif os.path.isdir(f):
                self.files += self.parsedir(f)
        if len(self.files) == 0:
            sys.exit("error: %s" % ('No valid \'_info.txt\' files provided.'))

    def parsedir(self,rootdir):
        """
            Returns valid _info files from directory and it's subdirectories
        """
        for subdir, dirs, files in os.walk(rootdir):
            to_parse = []
            for f in files:
                if f.endswith(".txt"):
                    to_parse.append(subdir+ "/" + f)
        return to_parse

    def topandas(self):
        """
            Load info files to pandas dataframe
        """
        pn('Importing info files...')
        df =  empty_dataframe() #empty numpy dataframe with set datatypes and header labels
        
        #iterate over info files to parse into np
        for f in self.files:
            po("Importing: %s" % (f))
            df = pd.concat([df, self.parse_file(f)], ignore_index=True)
        
        #deal with any remaining PositionA_ datasets
        if 'PositionA_' in df:
            df = df.drop(['PositionA'], axis=1)
            df['PositionA_'] = np.uint64(df['PositionA_'])
            df.rename(columns={'PositionA_':'PositionA'}, inplace=True)

        pn("Imported %d files" % (len(self.files)))

        return df

    def parse_file(self,fn):
        """ 
            Import single _info.txt file 
        """
        try:
         INFILE = open(fn, "rt")
        except OSError:
            print("Could not open/read file:", fn)
            sys.exit()
        with INFILE:
            df = pd.read_table(
                INFILE,
                dtype={
                    "Genotype": "str",
                    "ChromA": "str",
                    "ChromB": "str",
                    "PositionA": "uint64",
                    "PositionB": "uint64",
                    "Genotype_call": "str",
                    "Call_Agenome": "str",
                    "Call_Bgenome": "str",
                    "Major allele": "str",
                    "Minor allele": "str",
                },
            )
            return df

