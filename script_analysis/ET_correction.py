import csv
import math
import argparse

def skip_first(seq, n):
    for i,item in enumerate(seq):
        if i >= n:
            yield item

p = 0.05
cal_p = -math.log(p,10)

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--in_file", required = True, nargs = '+', help = "input merged vcf file")
ap.add_argument("-o", "--out_file", required = True, help = "output text file")
args = ap.parse_args()

if args.out_file:
   OUTPUT = open(args.out_file, 'w')

for infile in args.in_file:
    with open (infile, 'r') as csvfile:
        data = csv.reader(csvfile, delimiter='\t')
        for row in skip_first(data, 16):
            for idx in range(len(row)):
                if idx < 9:
                   OUTPUT.write(row[idx] + '\t')
                else:
                   field = row[idx].split(':')
                   if (len(field)) > 3 and field[3] != ".":
                      #print(field[3])
                      if float(field[3]) >= cal_p:
                         #print("here")
                         result = "./."
                      else:
                         result = "0/1"

                      if result != field[0]:
                         field[0] = result
                         #print(result)
                      for field_elm in field:
                          print(field_elm)
                          OUTPUT.write(field_elm + ":")
                      OUTPUT.write('\t')
                   else:
                       OUTPUT.write(row[idx] + '\t')
            OUTPUT.write('\n')
OUTPUT.close()



