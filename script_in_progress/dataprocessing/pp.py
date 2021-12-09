import pandas as pd
import numpy as np
from .utility import pn, po

class pp():
    def __init__(self, df):
        self.df = df

    def calc(self):
        pn("Calculating PP")
        self.df = self.df.groupby(by=['PositionA', 'PositionB', 'ChromA', 'ChromB']).apply(self.compute_pp)
        pn("PP calculated")
        return self.df


    def compute_pp(self, df):
        """
        Takes a dataframe grouped by 'Reference', 'PositionA_' and 'PositionB'
        and calculates posterior probability based on A/B genome calls
        """
        p = df.iloc[0]
        po("Calculating pp for: %s:%s\t%s:%s" % (p['ChromA'], p['PositionA'], p['ChromB'], p['PositionB']))
           
        unq_allele = set(''.join( df['Genotype_call'].unique().tolist()))
        unq_allele.discard("/")
        allele_a = set(''.join( df['Call_Agenome'].unique().tolist()))
        allele_b = set(''.join( df['Call_Bgenome'].unique().tolist()))
         
        """
        For all calculations, nan == discarded
        at position group, it is expected that pp_homo_a/b and pp_hetero is same for each genotype,
        but pp_hetero will (possibly) be variable per genotype
        """
        df['pp_het_a'] = self.computeHet_A(df)
        df['pp_het_b'] = self.computeHet_B(df)
        #calc a pp homo_a
        #allele_a is the number of unique alleles across all GT in this position in Call_Agenome
        df['pp_homo_a'] = self.computeHomo(df,'a',allele_a) #branching now done in computeHomo
        #df['pp_homo_a'] = 0
    
        #calc pp homo_b
        #allele_ is the number of unique alleles across all GT in this position in Call_Bgenome
        df['pp_homo_b'] = self.computeHomo(df,'b',allele_b) #branching now done in computeHomo
        #df['pp_homo_b'] = 0
        #calc pp homeo
        #unique_allele is the number of unique alleles across all GT in this position in Genotype_call (with / ignored)
        df['pp_homeo'] = self.computeHomeo(df, unq_allele)
        #df['pp_homeo'] = 0
    
        #calc pp hetero
        # No need to do anything. Comparing different files is not a requirement
        # PP_het_A = PPA(1)
        # PP_het_B = PPB(1)
        # for any given GT/[ps
    
        #returns the same df, with three new columns
        return df
    def computeHet_A(self,df):
        return(np.log(df["PA(1)"] / np.maximum(df["PA(0)"], df["PA(2)"])))
    
    def computeHet_B(self,df):
        return(np.log(df["PB(1)"] / np.maximum(df["PB(0)"], df["PB(2)"])))

    def chooseMajorAllele(self, x, alleles):
        
        #Counts occurence of individual alleles and returns [major,minor] list
        #if count N1 = count N2, defaults to whichever comes first in
        #lexographic order, otherwise major is allele with largest count
        
        alleles = sorted(alleles) #Put alleles in lexographic order
        n_g = []
        counts = []
        #setup major/minor allele for A
        # Peform allele count
        # allele with highest count = N1
        # other allele = N2
        for i, allele in enumerate(alleles):
            n_g.append(allele)
            counts.append(x.str.contains(allele).sum())
        if counts[0] < counts[1]: #if second allele is more common, flip array so n_a_2 is the first(major) allele
            n_g = n_g[::-1]
        return n_g
    
    def computeHomeo(self, x, alleles):
    #    """
    #        Calculate PP for homoeologous SNP
    #        PP_homoeologous = 0
    #        for N in N1, N2:			# consider genotype N1N1/N2N2, then N2N2/N1N1
    #        N' = {N1,N2}\N			# let N' be the other nucleotide
    #        PP = 1				# initialize product
    #        for k = 1 to n:			# for each individual
    #            PPi = 0			# default: current PP may not computed, in which case assume 0
    #            if Major allele[k] == N and Minor allele[k] == N':
    #                PPi = PP(0,2)[k]		# desired genotype is 0,2
    #            elif Minor allele[k] == N and Major allele[k] == N':
    #                PPi = PP(2,0)[k]		# desired genotype is 2,0
    #            PP *= PPi			# update product
    #            PP_homoeologous += PP		# PP is probability ALL genotypes are N1N1/N2N2 (1st loop) or N2N2/N1N1 (2nd loop)
    #    """
    
        #branches for too few/many unique alleles
        if len(alleles) == 1:
            return 0
        elif len(alleles) > 2:
            return np.nan
    
        n_g = self.chooseMajorAllele(x['Genotype_call'],alleles)
        def probN(x,n):
            
           # Every genotype in dataframe x has same major and minor allele, so we can calculate
           #PPgN for everything at once, saving calculation time
            
            row_0 = x.loc[x.index[0]]
            exp_alleles = [row_0['Major allele'], row_0['Minor allele']]
    
            x['n1'] = 0
            x['n2'] = 0
            x['n_max'] =0
            if n[0] == exp_alleles[0] and n[1] == exp_alleles[1]:
                #x['n1'] = x['PP(0,2)']
                x['n1'] = np.log(x['PP(0,2)'])
                x['n2'] = np.log(x['PP(2,0)'])
                x['n_max'] = np.log(np.vstack([x['PP(0,0)'], x['PP(0,1)'], x['PP(0,2)'], x['PP(1,0)'], x['PP(1,1)'], x['PP(1,2)'], x['PP(2,0)'], x['PP(2,1)'], x['PP(2,2)']]).max(axis=0))
            elif n[0] == exp_alleles[1] and n[1] == exp_alleles[0]:
                #x['n1'] = x['PP(2,0)']
                x['n1'] = np.log(x['PP(2,0)'])
                x['n2'] = np.log(x['PP(0,2)'])
                x['n_max'] = np.log(np.vstack([x['PP(0,0)'], x['PP(0,1)'], x['PP(0,2)'], x['PP(1,0)'], x['PP(1,1)'], x['PP(1,2)'], x['PP(2,0)'], x['PP(2,1)'], x['PP(2,2)']]).max(axis=0))
            return x
    
        x = x.groupby(by=['Major allele', 'Minor allele']).apply(probN,n_g)
        x['n1'] = x['n1'].replace(-np.inf, 0)
        #print(x['n2'])
        x['n2'] = x['n2'].replace(-np.inf, 0)
        #print(x['n2'])
        #x['n1'] = x['n1'].replace(-np.inf, 0, inplace=True)
        #x['n2'] = x['n2'].replace(-np.inf, 0, inplace=True)
        #Π[1{N1}Pr(0,2) + 1{N2}Pr(2,0)]
        #p_n1 = x['n1'].prod()
        p_n1 = x['n1'].sum()
        #Π[1{N2}Pr(0,2) + 1{N2}Pr(2,0)]
        #p_n2 = x['n2'].prod()
        p_n2 = x['n2'].sum()
        p_max = x['n_max'].sum()
        #print(p_n1, p_n2, p_max)
        #print(p_n1 + p_n2 - p_max)
        return p_n1 + p_n2 - p_max
    
    
    def computeHomo(self,x,g,alleles):
        """
        # Calculate PP for homologous SNP
        for g in A, B:
           PP_homo[g] = 1
           if g == A:			# because of the awkard naming; should use index instead of name
              PG = PA
           else:
              PG = PB
           for N in N1, N2:
               PPi = 1			# initialize product
               for k = 1 to n:
                   PPg = 0			# default, 0
                   if Major allele[k] == N:	# want PP of g = 0
                      PPg = PG(0)[k]
                   elif Minor allele[k] == N:# want PP of g = 2
                      PPg = PG(2)[k]
                   PPi *= PPA
                PP_homo[g] -= PPi		# PPi is probability ALL genotypes are NN
        """
        #TODO: Calculate A,B at same time removing need for g and shaving off a execution of this function
    
        #branches for too few/many unique alleles
        if len(alleles) == 1:
            return 0
        elif len(alleles) > 2:
            return np.nan 
    
        alleles = sorted(alleles) #Put alleles in lexographic order
        pp_0 = ''
        pp_1 = ''
        pp_2 = ''
        pp_max = ''
        pp_homo = 1.0
        pp_g = 0.0
        call = ''
        if g == 'a':
            pp_0 = 'PA(0)'
            pp_1 = 'PA(1)'
            pp_2 = 'PA(2)'
            call = 'Call_Agenome'
        elif g == 'b':
            pp_0 = 'PB(0)'
            pp_1 = 'PB(1)'
            pp_2 = 'PB(2)'
            call = 'Call_Bgenome'
        else:
            return np.nan
    
        n_g = self.chooseMajorAllele(x[call], alleles)
    
        def probN(x,n):
            
            #Every genotype in dataframe x has same major and minor allele, so we can calculate
            #PPgN for everything at once, saving calculation time
            
            row_0 = x.loc[x.index[0]]
            exp_alleles = [row_0['Major allele'], row_0['Minor allele']]
            x['n1'] = np.log(x[pp_0])
            x['n2'] = np.log(x[pp_2])
            x['n_max'] = np.log(np.vstack([x[pp_0], x[pp_1], x[pp_2]]).max(axis=0))

            if n[0] not in exp_alleles: #check N1 in {Major allele[k], Minor allele[k]}
                x['n1'] = 0
                x['n_max'] = 0
            elif n[0] == exp_alleles[1]: #check N1 == Minor allele[k]
                x['n1'] = np.log(x[pp_2])
                x['n_max'] = np.log(np.vstack([x[pp_0], x[pp_1], x[pp_2]]).max(axis=0))
    
            if n[1] not in exp_alleles: #check N2 in {Major allele[k], Minor allele[k]}
                x['n2'] = 0;
            elif n[1] == exp_alleles[0]: #check N2 == Major allele[k]
                x['n2'] = np.log(x[pp_0])
                x['n_max'] = np.log(np.vstack([x[pp_0], x[pp_1], x[pp_2]]).max(axis=0))
            return x
    
        #mass transform dataframe by logical groups of same Major and Minor allele at the current position
        x = x.groupby(by=['Major allele', 'Minor allele']).apply(probN,n_g)
        
        #PPg = 1 - {Π[1{N1}Pr(Pg0) + 1{N1}Pr(Pg2)] + Π[1{N2}Pr(Pg0) + 1{N2}Pr(Pg2)]
        #Π[1{N1}Pr(Pg0) + 1{N1}Pr(Pg2)]
        #p_n1 = x['n1'].prod()
        p_n1 = x['n1'].sum()
        #Π[1{N2}Pr(Pg0) + 1{N2}Pr(Pg2)]
        #p_n2 = x['n2'].prod()
        p_n2 = x['n2'].sum()
       
        p_max = x['n_max'].sum()
        #return 1 - p_n1 + p_n2
        #print(p_n1, p_n2, p_max)
        #print(p_n1 + p_n2 - p_max)
        return p_n1 + p_n2 - p_max 

