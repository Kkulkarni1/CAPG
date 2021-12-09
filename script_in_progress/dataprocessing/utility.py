import numpy as np
import pandas as pd
import sys

from decimal import Decimal
from io import StringIO

def empty_dataframe():
    headers = (
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
        % (
            "Genotype",
            "ChromA",
            "ChromB",
            "PositionA",
            "PositionB",
            "Genotype_call",
            "Call_Agenome",
            "Call_Bgenome",
            "Major allele",
            "Minor allele",
            "PP(0,0)",
            "PP(0,1)",
            "PP(0,2)",
            "PP(1,0)",
            "PP(1,1)",
            "PP(1,2)",
            "PP(2,0)",
            "PP(2,1)",
            "PP(2,2)",
            "PA(0)",
            "PA(1)",
            "PA(2)",
            "PB(0)",
            "PB(1)",
            "PB(2)",
            "CovA",
            "CovB",
        )
    )
    return pd.read_table(
        StringIO(headers),
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

def pn(msg):
    po(msg)
    sys.stderr.write('\n')
    sys.stderr.flush()

def po(msg):
    sys.stderr.write('\r')
    sys.stderr.write(msg.ljust(85,' '))
    sys.stderr.flush()
