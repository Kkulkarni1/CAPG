import numpy as np
import pandas as pd

class Het():
    def __init__(self,outfile):
        self.het_name = outfile+'_het.txt'
        self.ho_name = outfile+'_Ho.txt'
        self.all_name = outfile+'_all.txt'

    def printHetFiles(self, df):
        #repace NaN with 0 or drop rows with NaN values
        df.fillna(0, inplace=True)
        #df.dropna(inplace=True)

       # df['PositionA'] = np.uint64(df['PositionA'])
       # df['PositionB'] = np.uint64(df['PositionB'])

        #print PP(hetA) PP(hetB) file
        df[['PositionA', 'PositionB','Genotype','PP(hetA)', 'PP(hetB)','CovA','CovB']].sort_values(by=['Genotype','PositionA','PositionB']).to_csv(self.het_name, sep='\t', line_terminator='\n', index=False, float_format='%e')

        #calculate PP(Ho) 
        #Print file with all calculations
        df[['PositionA', 'PositionB', 'Genotype','PP(homeo)','PP(homoA)','PP(homoB)','PP(hetA)', 'PP(hetB)']].sort_values(by=['Genotype','PositionA','PositionB']).to_csv(self.all_name, sep='\t', line_terminator='\n', index=False, float_format='%e')
        
        #Collapse by position to unique values only
        #print( df[['PositionA','PositionB','CovA','CovB']].to_csv(float_format='%e', index=False))
        #df2['CovA'] = df.groupby(['PositionA','PositionB'])['CovA'].agg('mean')
        df['CovA'] = df.groupby(by=['PositionA','PositionB'])['CovA'].transform(lambda x: x.agg('mean'))
        df['CovB'] = df.groupby(by=['PositionA','PositionB'])['CovB'].transform(lambda x: x.mean())
        df2 = df[['PositionA','PositionB', 'PP(homeo)', 'PP(homoA)','PP(homoB)', 'CovA','CovB']].drop_duplicates()
        #Print by position
        df2.to_csv(self.ho_name,sep='\t',index=False,line_terminator='\n', float_format='%e')
