import sys,re,gzip
import pprint as pp
from datetime import date



def extractGenename(line) :
    regex = re.compile(r'gene_name\s"(.*)?";')
    ret = regex.search(line)
    if ret :
        genename = ret.group(1)
    else :
        genename = ''
    return(genename)

ipdfile = sys.argv[1]
c4file = sys.argv[2]
geneinfofile = sys.argv[3]
geneblacklist = sys.argv[4]


geneExclude = []
with open(geneblacklist, 'rt') as f:
    while True:
        line = f.readline()
        if not line:
            break
        geneExclude.append(line.replace('\n', ''))
 

geneinfo = {}
with gzip.open(geneinfofile, 'rt') as f:
    while True:
        line = f.readline()
        if not line:
            break
        line = line.replace('\n', '')
        llst = line.split(',')
        rec = {}
        rec['source'] = llst[1]
        rec['geneid'] = llst[2]
        rec['transid'] = llst[3]
        rec['copyindex'] = 0
        geneinfo[llst[0]] = rec

ctgs = {}
bodyout = []
with gzip.open(ipdfile, 'rt') as f :
    while True:
        line = f.readline()
        if not line:
            break
        bodyout.append(line)
        genename = extractGenename(line)
        llst = line.split('\t')
        if llst[1] not in ctgs.keys() :
            ctgs[llst[1]] = [llst[0]]
        else :
            if llst[0] not in ctgs[llst[1]] :
                ctgs[llst[1]].append(llst[0])
        if llst[2] == 'gene' :
            geneinfo[genename]['copyindex'] += 1

with gzip.open(c4file, 'rt') as f :
    while True:
        line = f.readline()
        if not line:
            break
        bodyout.append(line)
        genename = extractGenename(line)
        genename = genename[0:(len(genename)-1)]
        llst = line.split('\t')
        if llst[1] not in ctgs.keys() :
            ctgs[llst[1]] = [llst[0]]
        else :
            if llst[0] not in ctgs[llst[1]] :
                ctgs[llst[1]].append(llst[0])
        if llst[2] == 'gene' :
            geneinfo[genename]['copyindex'] += 1

bodyout = sorted(
    bodyout, key = lambda x : int(x.split('\t')[3]))

gene0 = []
gene1 = []
genem = []
headerout = []

#pp.pprint(geneinfo, stream=sys.stderr)

for g in geneinfo.keys() :
    cpnum = geneinfo[g]['copyindex']
    if g in geneExclude : continue
    if cpnum == 0 :
        gene0.append(g)
    elif cpnum == 1:
        gene1.append(g)
    else  :
        genem.append(g)

headerout.append("## format: gtf\n")
headerout.append("## date: " +  str(date.today()) + "\n")
if len(gene0) > 0 :
    headerout.append("## gene (copy num = 0): " + ",".join(gene0) + "\n")
if len(gene1) > 0 :
    headerout.append("## gene (copy num = 1): " + ",".join(gene1) + "\n")
if len(genem) > 0 :
    headerout.append("## gene (copy num > 1): " + ",".join(genem) + "\n")
for k in ctgs.keys() :
    headerout.append("## contigs for " + k + ": " + ",".join(ctgs[k]) + "\n")

sys.stdout.write(''.join(headerout))
sys.stdout.write(''.join(bodyout))
