import sys,re,gzip
import genebanktools as gbt


def writeFa(seq, width=80) :
    n = len(seq) // width
    ostr = ''
    for i in range(n + 1) :
        a = i*width
        b = a + width
        if b >= len(seq) :
            b = len(seq)
        elif b == a :
            break
        ostr += seq[a:b] + '\n'
    return(ostr)


inpfile = sys.argv[1]
outpref = sys.argv[2]

#output
outfile = outpref + ".exon.fa.gz"

buff = ''
with open(inpfile, 'r') as f:
    while True :
        line = f.readline()
        if not line : break
        buff += line
        if line[0:2] == "//" :
            exon = gbt.extractExon(buff)
            gene, codon_start, cds = gbt.extractCDS(buff)
            gseq = gbt.extractSeq(buff)
            buff = ''
if len(buff) > 10 :
    exon = gbt.extractExon(buff)
    gene, codon_start, cds = gbt.extractCDS(buff)
    gseq = gbt.extractSeq(buff)
    buff = ''

# output exon
exonarr = []
exon = exon.split(',')
minExon = 99999
for x in exon :
    a,bc = x.split(':')
    b, c = bc.split('..')
    a = int(a)
    b = int(b)
    c = int(c)
    if minExon > c - b + 1:
        minExon = c - b + 1
    exonarr.append([a, b, c])


cdsarr = []
cds = cds.split(',')
for x in cds :
    a, b = x.split('..')
    a = int(a)
    b = int(b)
    cdsarr.append([a, b])

if minExon < 20 :
    sys.stderr.write("<msg> Warning: min exon length = " + str(minExon) 
                     + ', too short for split alignment\n')
exonarr = sorted(exonarr, key = lambda x : x[1])
cdsarr = sorted(cdsarr, key = lambda x : x[1])

exonSeq = ''
for x in exonarr:
    gseqstr = gseq[x[1]-1:x[2]]
    if x[0] == 26 :
        AseqRegex = re.compile(r'(CCCTGT[ATCG]{6}TTAGAC)')
        BseqRegex = re.compile(r'(CTCTCT[ATCG]{6}ATACAT)')
        reta = AseqRegex.search(gseqstr)
        retb = BseqRegex.search(gseqstr)
        if reta :
            p0 = reta.start() + 1
            p1 = reta.end()
        elif retb :
            p0 = retb.start() + 1
            p1 = retb.end()
        else :
            sys.exit("Error: Cannot find out the key feature for A/B allele\n")
        p0 += len(exonSeq)
        p1 += len(exonSeq)
        AB_blockstr = str(p0) + '..' + str(p1)
    if x[0] == 10 :
        intron9_blockstr = str(len(exonSeq)) + '..' + str(len(exonSeq)+1)
    exonSeq += gseqstr
## new exon block
currentPos = 1
for i in range(len(exonarr)):
    x = exonarr[i]
    a = x[0]
    b = currentPos
    c = currentPos + x[2] - x[1]
    x.append(b)
    x.append(c)
    exonarr[i] = x
    currentPos = c + 1

## new CDS block based on exon seq
k = 0
l = len(exonarr)
for i in range(len(cdsarr)) :
    a, b = cdsarr[i]
    for j in range(k, l) :
        x = exonarr[j]
        c1 = x[1] <= a and a < x[2]
        c2 = x[1] < b and b <= x[2]
        if c1 and c2 :
            k = j + 1
            c = x[3] + a - x[1]
            d = x[3] + b - x[1]
            break
    cdsarr[i].append(c)
    cdsarr[i].append(d)

## search start and stop codon position
b0 = cdsarr[0][2] - 1
p0 = b0 + int(codon_start) - 1
if exonSeq[p0 : p0 +3] == "ATG" :
    start_codon_bloc = [p0+1, p0+3]
else :
    start_codon_bloc = []
    print("no start codon", file = sys.stderr)

p1 = cdsarr[-1][3] - 3
c0 = exonSeq[p1 : p1 +3] in ["TAG", "TGA", "TAA"]
c1 = p1 +6 <= len(exonSeq)
c2 = exonSeq[p1 +3 : p1 +6] in ["TAG", "TGA", "TAA"]
if c0 :
    stop_codon_bloc = [p1+1, p1+3]
elif c1 and c2 :
    stop_codon_bloc = [p1+4, p1+6]
else :
    stop_codon_bloc = []
    print("no stop codon", file = sys.stderr)

##split codon --> rare but possible
## to check whether the codon include an exon boundary
exonboundaries = []
for x in exonarr :
    exonboundaries.append(x[3])

## only one split is allowed
for x in exonboundaries :
    if stop_codon_bloc :
        if stop_codon_bloc[0] + 1 == x :
            stop_codon_bloc = [x-1, x-1, x, x+1]
        elif stop_codon_bloc[0] + 2 == x :
            stop_codon_bloc = [x-2, x-1, x, x]
    if start_codon_bloc :
        if start_codon_bloc[0] + 1 == x:
            start_codon_bloc = [x-1, x-1, x, x+1]
        elif start_codon_bloc[0] + 2 == x:
            start_codon_bloc = [x-2, x-1, x, x]

## start and stop codon string
startstr = ''
if start_codon_bloc :
    x = start_codon_bloc
    startstr = str(x[0]) + '..' + str(x[1])
    if len(x) == 4 :
        startstr += "," + str(x[2]) + '..' + str(x[3])
stopstr = ''
if stop_codon_bloc :
    x = stop_codon_bloc
    stopstr = str(x[0]) + '..' + str(x[1])
    if len(x) == 4 :
        stopstr += "," + str(x[2]) + '..' + str(x[3])

## exon block string
exonstr = []
for x in exonarr :
    exonstr.append(str(x[0]) + ':' + str(x[3]) + '..' + str(x[4]))
exonstr = ','.join(exonstr)
cdsstr = []
for x in cdsarr :
    cdsstr.append(str(x[2]) + '..' + str(x[3]))
cdsstr = ','.join(cdsstr)

fp = gzip.open(outfile, 'wt')
ostr = '>' + gene
ostr += " codon_start=" + codon_start
ostr += " start_codon=" + startstr
ostr += " stop_codon=" + stopstr
ostr += " ab_block=" + AB_blockstr
ostr += " intron9_block=" + intron9_blockstr
ostr += " cds_block=" + cdsstr
ostr += " exon_block=" + exonstr + '\n'
fp.write(ostr)
fp.write(writeFa(exonSeq))
fp.close()

