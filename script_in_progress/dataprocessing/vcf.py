import sys
import re
import datetime
import tempfile
from .utility import pn, po

class Vcf():
    def __init__(self,vcf_args):
        outfile_prefix = vcf_args[0][0]
        self.a_out = outfile_prefix + "_A.vcf"
        self.b_out = outfile_prefix + "_B.vcf"
        self.ab_out = outfile_prefix + "_AB.vcf"
        self.a_ref = False
        self.b_ref = False
        self.pseudo = False
        if vcf_args[1]:
            self.a_ref = vcf_args[1]
            self.a_data = self.load_fasta(vcf_args[1][0], vcf_args[4])
        if vcf_args[2]:
            self.b_data = self.load_fasta(vcf_args[2][0],vcf_args[4])
            self.b_ref = vcf_args[2]
        if vcf_args[3]:
            self.pseudo = True

    def print_vcf_from_dataframe(self,df):
        a_data = {}
        b_data = {}
        vcf_open = []
        datestamp = datetime.datetime.now().strftime("%Y%m%d")
        gtl = df["Genotype"].unique().tolist()
        gts = "\t".join(map(str, gtl))
    
        if self.a_ref:
            vcf_a = open(self.a_out, "wt")
            vcf_open.append(vcf_a)
            self.init_vcf(self.a_ref, vcf_a, 1, datestamp, gts)
        if self.b_ref:
            vcf_b = open(self.b_out, "wt")
            vcf_open.append(vcf_b)
            self.init_vcf(self.b_ref, vcf_b, 2, datestamp, gts)
        if self.pseudo:
            vcf_ab = open(self.ab_out, "wt")
            vcf_open.append(vcf_ab)
            self.init_vcf([self.a_ref[0], self.b_ref[0]], vcf_ab, 3, datestamp, gts)
    
        # control variable for printing files
        print_genomes = 0
        if self.pseudo:
            print_genomes = 4  # a homo, b homo, ab homeo
        elif self.a_ref and self.b_ref:
            print_genomes = 3  # a and b homo only
        elif self.a_ref:
            print_genomes = 1  # a-genome only
        else:
            print_genomes = 2  # b-genome only
        pn("Printing VCF files") 
        df.groupby(by=["ChromA", "ChromB", "PositionA", "PositionB"]).apply(
            self.print_vcf_line, vcf_open, print_genomes, gtl
        )
        pn("Finished printing VCF files")


    def load_fasta(self, fa_file, valid_keys):
        loaded_data = {}
        try:
            pn("Reading reference file %s" % (fa_file))
            f = open(fa_file, "rt")
        except OSError:
            sys.exit("error: %s\t%s\n" % ("cannot open/read file", fa_file))
    
        with f:
            current_key = ""
            load = False
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    load = False
                    current_key = re.split(">| - |\| ", line)[1]
                    if current_key in valid_keys:
                        po("Loading %s from fasta" % (current_key))
                        loaded_data[current_key] = {"data":[], "line_size":0}
                        load = True
                elif load:
                    loaded_data[current_key]["data"].append(line)
                    if loaded_data[current_key]["line_size"] == 0:
                        loaded_data[current_key]["line_size"] = len(line)
            
            pn("Loaded reference file %s" % (fa_file))
            return loaded_data


    def init_vcf(self, in_file, vcf, genome, datestamp, gts):
        lines = [
            "##fileformat=VCFv4.2\n",
            "##fileDate=%s\n" % (datestamp),
            "##source=info2pandas.py\n",
        ]
        vcf.writelines(lines)
        for ref in in_file:
            vcf.write("##reference=file://%s\n" % (ref))
        if genome == 1:
            vcf.write(
                '##INFO=<ID=PPHA,Number=A,Type=Float,Description="Posterior probability of being a homologous SNP in A genome">\n'
            )
        elif genome == 2:
            vcf.write(
                '##INFO=<ID=PPHA,Number=A,Type=Float,Description="Posterior probability of being a homologous SNP in B genome">\n'
            )
        else:
            vcf.write(
                '##INFO=<ID=PPHA,Number=A,Type=Float,Description="Posterior probability of being a homeologous SNP">\n'
            )
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        if genome == 1 or genome == 2:
            vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (gts))
        else:
            vcf.write(
                '##FORMAT=<ID=DPA,Number=1,Type=Integer,Description="Read Depth A">\n'
            )
            vcf.write(
                '##FORMAT=<ID=DPB,Number=1,Type=Integer,Description="Read Depth B">\n'
            )
            vcf.write("#CHROM\tPOS_A\tPOS_B\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (gts))


    def print_vcf_line(self, df, vcf, genome, gtl):
        pos_a = int(df["PositionA"].unique()[0])
        pos_b = int(df["PositionB"].unique()[0])
        chrom_a = df["ChromA"].unique()[0]
        chrom_b = df["ChromB"].unique()[0]
        po("Generating vcf for %s:%d %s:%d" % (chrom_a,pos_a,chrom_b,pos_b))
        vid = "."
        qual = "."
        vfilter = "."
        vformat = "GT:DP"
        cov_a = df["CovA"]
        cov_b = df["CovB"]
        gts = df["Genotype"].values
    
        def combine_series(ser_a, ser_b):
            return "%s:%s" % (ser_a, ser_b)
    
        def reformat_and_combine_series(ser_a, ser_b):
            # ref and alt are grabbed from the hosting scope
            ser_a = ser_a.replace("/", "")
            ser_a = ser_a.replace(ref, "0")
            for idx, alt in enumerate(alt_set, start=1):
                ser_a = ser_a.replace(alt, "%s" % (idx))
            ser_a = ser_a.replace("", "/")[1:-1]
            return combine_series(ser_a, ser_b)
    
        if genome != 2:  # print line in a
            if chrom_a in self.a_data:
                alt = "."
                line_size = self.a_data[chrom_a]["line_size"]
                idx_a = (pos_a-1) // line_size
                ref = self.a_data[chrom_a]["data"][idx_a][(pos_a - 1) % line_size]
                # filter out unique single bases for alt
                alt_set = set("".join(df["Call_Agenome"].unique().tolist()))
                alt_set.discard(ref)
                if len(alt_set) > 0:
                    alt = ",".join(alt_set)
                info = "PPHA=%e" % (df["PP(homoA)"].unique())
                format_data = []
                for gt in gtl:
                    if gt in gts:
                        dl = df.loc[df["Genotype"] == gt]
                        format_data.append(
                            "".join(map(
                                str,
                                dl["Call_Agenome"].combine(dl["CovA"], func=(reformat_and_combine_series)),
                            )
                        ))


                    else:
                        format_data.append("./.")
                format_data = "\t".join(format_data)
    
                # format_data = "\t".join(map(str, df["Call_Agenome"]))
                vcf[0].write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (chrom_a, pos_a, vid, ref, alt, qual, vfilter, info, vformat, format_data)
                )
            else:
                pn("Discarded line in a vcf: %s not found in %s\nPosA: %d" % (chrom_a, self.a_ref, pos_a))
    
        if genome != 1:  # print line in b
            if chrom_b in self.b_data:
                alt = "."
                line_size = self.b_data[chrom_b]["line_size"]
                idx_b = (pos_b-1) // line_size
                ref = self.b_data[chrom_b]["data"][idx_b][(pos_b - 1)%line_size]
                # filter out unique single bases for alt
                alt_set = set("".join(df["Call_Bgenome"].unique().tolist()))
                alt_set.discard(ref)
                if len(alt_set) > 0:
                    alt = ",".join(alt_set)
                info = "PPHA=%e" % (df["PP(homoB)"].unique())
                format_data = []
                for gt in gtl:
                    if gt in gts:
                        dl = df.loc[df["Genotype"] == gt]
                        format_data.append(
                            "".join(map(
                                str,
                                dl["Call_Bgenome"].combine(dl["CovB"], func=(reformat_and_combine_series)),
                            )
                        ))


                    else:
                        format_data.append("./.")
                format_data = "\t".join(format_data)
                # format_data = "\t".join(map(str, format_data))
                if genome == 2:
                    vcf[0].write(
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                        % (chrom_b, pos_b, vid, ref, alt, qual, vfilter, info, vformat, format_data)
                    )
                else:
                    vcf[1].write(
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                        % (chrom_b, pos_b, vid, ref, alt, qual, vfilter, info, vformat, format_data)
                    )
            else:
                pn("Discarded line in b vcf: %s not found in %s\nPosB: %d" % (chrom_b, self.b_ref,pos_b))
                
    
        if genome == 4:  # print line in ab
            if chrom_a in self.a_data and chrom_b in self.b_data:
                alt = "."
                chrom = chrom_a + "/" + chrom_b
                line_size_a = self.a_data[chrom_a]["line_size"]
                line_size_b = self.b_data[chrom_b]["line_size"]
                idx_a = (pos_a-1) // line_size_a
                idx_b = (pos_b-1) // line_size_b
                ref = self.a_data[chrom_a]["data"][idx_a][(pos_a - 1) % line_size_a] + self.b_data[chrom_b]["data"][idx_b][(pos_b - 1) % line_size_b]
                alt_set = set(
                    df["Call_Agenome"]
                    .append(df["Call_Bgenome"], ignore_index=True)
                    .unique()
                    .tolist()
                )
                alt_set.discard(ref)
                if len(alt_set) > 0:
                    alt = ",".join(alt_set)
                info = "PPH0=%e" % (df["PP(homeo)"].unique())
                format_data = df["Genotype_call"].combine(
                    cov_a, func=(reformat_and_combine_series)
                )
                format_data = "\t".join(
                    map(str, format_data.combine(cov_b, func=(combine_series)))
                )
                # format_data = "\t".join(map(str, format_data))
                vcf[2].write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (chrom, pos_a, pos_b, vid, ref, alt, qual, vfilter, info, vformat, format_data)
                )
            else:
                pn("Discarded line in pseudo vcf: %s or %s not found in fasta files\nPositionA: %d\t PositionB:%d" % (chrom_a, chrom_b, pos_a, pos_b))


