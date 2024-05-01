import re, sys
import pprint as pp
import argparse
import time
from datetime import timedelta
from datetime import datetime

#sys.path.insert(0, './easyHLA-src')
import IPDtools as ipdt

# program interface

parser = argparse.ArgumentParser(
    prog = 'cat DAT | python3 easyHLA.py',
    description = 'A lightweight tool for extracting information \
    from IPD-HLA database',
    epilog = 'easyHLA takes one scan of the IPD-HLA database, extract\
    information for all alleles or sellected alleles by adding\
    option --alleleF or --alleleS based on the way of defining the name\
    list. Both HLA ID such as "HLA00001", and allele name such as\
    "HLA-A*01:01:01:01" are accepted. Four flag options are used for\
    defining output, only one flag is allowed for one run')

parser.add_argument('-A', '--alleleF', type=str, required=False,
                    help='[string] filename of allele names to be extracted')

parser.add_argument('-a', '--alleleS', type=str, required=False,
                    help='[string] allele names (sep=\',\') to be extracted')

parser.add_argument('-AN', '--alleleName', 
                    action='store_true',
                    help='[flag] extract all allele names')


parser.add_argument('-G', '--geneSeq', 
                    action='store_true',
                    help='[flag] extract gene sequence')

parser.add_argument('-P', '--protSeq', 
                    action='store_true',
                    help='[flag] extract protein sequence')

parser.add_argument('-C', '--cdsSeq', 
                    action='store_true',
                    help='[flag] extract CDS sequence')

parser.add_argument('-E', '--elements', 
                    action='store_true',
                    help='[flag] extract ranges for gene, CDS, exons, UTR, codon_start')

parser.add_argument('-V', '--verbose', default=True, 
                    action='store_true',
                   help='[flag] on/off running msg')  # on/off flag

sys.stderr.write("<---------------------------------------->\n")
sys.stderr.write("<---------------------------------------->\n")
startTime = time.monotonic()
sys.stderr.write("Start time: " + str(datetime.now()) + "\n")

args = parser.parse_args()
#if args.verbose : print(args)
# allele filtration to output

if args.alleleF and args.alleleS :
    sys.exit('only one option is allowed: --alleleF or --alleleS')

if args.alleleS :
    alleleSets = args.alleleS.split(',')
elif args.alleleF :
    try :
        alleleSets = open(args.alleleF, 'r').readlines()
        alleleSets = [s.strip('\n') for s in alleleSets]
    except :
        sys.exit('fail to open '+ args.alleleList)
else :
    alleleSets = None
    if args.verbose : sys.stderr.write(" <msg> all alleles are selected\n")

nout = int(args.geneSeq + args.protSeq + args.cdsSeq + args.elements + args.alleleName)

if nout == 0 :
    sys.exit('add one tag: --[geneSeq|protSeq|cdsSeq|elements|alleleName]')
elif nout > 1 :
    sys.exit('too many tags, use only one:\
             --[geneSeq|protSeq|cdsSeq|elements|alleleName]')
else :
    if args.verbose :
        if args.alleleName: 
            sys.stderr.write(" <msg> will extract all allele names\n")
        if args.geneSeq: 
            sys.stderr.write(" <msg> will extract gene sequence\n")
        if args.protSeq: 
            sys.stderr.write(" <msg> will extract protein sequence per gene\n")
        if args.cdsSeq: 
            sys.stderr.write(" <msg> will extract cds sequences per gene\n")
        if args.elements: 
            sys.stderr.write(" <msg> will extract elements per gene\n")

# main program

df = {}
buf = ''

firstLineTag = True

for line in sys.stdin:
    if 'Exit' == line.rstrip():
        break
    buf = buf + line
    if line[0:2] == '//':
        df.update(ipdt.readMETA(buf))
        #print(alleleSets)
        if alleleSets != None :
            if not ( 
                df['ID'] in alleleSets 
                or df['allele'] in alleleSets 
            ) : 
                df = {}
                buf = ''
                continue
            else :
                if df['ID'] in alleleSets : 
                    alleleSets.remove(df['ID'])
                if df['allele'] in alleleSets : 
                    alleleSets.remove(df['allele'])

        #print(df['ID'] + '-2')
        #if df['ID']== 'HLA06674': 
        df.update(ipdt.readGENE(buf))

        if args.alleleName:
            ipdt.writeSTDOUTalleleName(df)
        elif args.elements:
            df.update(ipdt.readGENESeq(buf))
            ipdt.writeSTDOUTgeneElement(df)
            #if firstLineTag :
            #    ipdt.writeSTDOUTgeneRecordsFirstline()
            #    firstLineTag = False
            #ipdt.writeSTDOUTgeneRecords(df)
        elif args.protSeq: 
            df.update(ipdt.readPROTEIN(buf))
            ipdt.writeSTDOUTprotein(df)
        elif args.geneSeq:
            df.update(ipdt.readGENESeq(buf))
            ipdt.writeSTDOUTgeneseq(df)
        elif args.cdsSeq:
            df.update(ipdt.readGENESeq(buf))
            ipdt.writeSTDOUTcds(df, args.verbose)
        #pp.pprint(df)
        df = {}
        buf = ''

endTime = time.monotonic()
#sys.stderr.write(timedelta(seconds=endTime - startTime))

if df == {} :
    if alleleSets :
        sys.stderr.write(" <msg> allele missed in the input data:\n")
        sys.stderr.write('\n'.join(alleleSets))
        sys.stderr.write('\n')

    sys.stderr.write("End time: " + str(datetime.now()) + "\n")
    sys.stderr.write("Duration: " 
                     + str(timedelta(seconds=endTime - startTime)) + "\n")
    sys.exit("<-------------------Done------------------>\n\n")
else :
    sys.exit("Warning: wrong input format, incomlpete entry for allele: " + df['gene'] + df['ID'] + df['allele'])
