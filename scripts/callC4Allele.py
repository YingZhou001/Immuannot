import sys,re
import gzip
import paftools as pt
import algntools as at

removeStopCodonFromCDS = True

def callC4allele(intron9seq, exon26keyseq) :
    L = len(intron9seq)
    if L > 5000 :
        size = 'L'
    else :
        size = 'S'
    # A or B allele
    allele = 'C4'
    AseqRegex = re.compile(r'(CCCTGT[ATCG]{6}TTAGAC)')
    BseqRegex = re.compile(r'(CTCTCT[ATCG]{6}ATACAT)')
    reta = AseqRegex.search(exon26keyseq)
    retb = BseqRegex.search(exon26keyseq)
    tag = ''
    if reta :
        allele = 'C4A'
        tag = reta.group(1)
    elif retb:
        allele = 'C4B'
        tag = retb.group(1)
    else :
        allele = 'C4X'
        tag = exon26keyseq
    if tag :
        tag = tag[0:6] + tag[12:18]
    
    return([allele, size, tag])
 
def gtfline(seqid, source, feature, fro, to, score, strand,
                 phase, attribute) :
        line = ''
        line += str(seqid)
        line += '\t' + str(source)
        line += '\t' + str(feature)
        line += '\t' + str(fro)
        line += '\t' + str(to)
        line += '\t' + str(score)
        line += '\t' + str(strand)
        line += '\t' + str(phase)
        line += '\t' + str(attribute)
        return(line)

exonPaffile = sys.argv[1] # exon based split alignment
exonCtgfile = sys.argv[2]
tarCtgfile = sys.argv[3]
geneinfofile = sys.argv[4]

removeStopCodonFromCDS = True


expafdata = []
with gzip.open(exonPaffile, 'rt') as f:
    while True :
        line = f.readline()
        if not line : break
        newdct = pt.parsingPaf(line)
        expafdata.append(newdct) ### add any filtrations?

expafdata = sorted(expafdata, key = lambda x: int(x['tfro']))
expafdata = sorted(expafdata, key = lambda x: int(x['match']))
expafdata = sorted(expafdata, key = lambda x: int(x['NM']))

ctgids = []
for x in expafdata :
    if not x['tstr'] in ctgids :
        ctgids.append(x['tstr'])

## select unique mappings
newexpafdata = []

for ctgid in ctgids :
    pafda = []
    selpafda = []
    for da in expafdata :
        if da['tstr'] == ctgid :
            pafda.append(da)
    for da in pafda :
        newtag = 1
        if selpafda :
            f = int(da['tfro'])
            t = int(da['tto'])
            for da0 in selpafda :
                f0 = int(da0['tfro'])
                t0 = int(da0['tto'])
                if (f >= f0 and f <= t0 ) or (t >= f0 and t <= t0 ) :
                    newtag = 0
                    continue
        if newtag == 1 :
            selpafda.append(da)
    for da in selpafda :
        newexpafdata.append(da)


exonSeq = ''
with gzip.open(exonCtgfile, 'rt') as f:
    while True :
        line = f.readline()
        if not line : break
        line = line.replace('\n', '')
        if line[0] == '>' :
            llst = line.split(' ')
            codonStart = llst[1].split('=')[1]
            startcodonblock = llst[2].split('=')[1]
            stopcodonblock = llst[3].split('=')[1]
            abblock = llst[4].split('=')[1]
            intron9block = llst[5].split('=')[1]
            cdsblock = llst[6].split('=')[1]
            exonblock = llst[7].split('=')[1]
        else :
            exonSeq += line

geneInfo = {}
with gzip.open(geneinfofile, 'rt') as f :
    while True :
        line = f.readline()
        if not line : break
        line = line.replace('\n', '')
        llst = line.split(',')
        #sys.stderr.write(str(llst)+'\n')
        newdict = {}
        newdict['source'] = llst[1]
        newdict['geneid'] = llst[2]
        newdict['transid'] = llst[3]
        newdict['copyindex'] = 0
        geneInfo[llst[0]] = newdict


for epaf in newexpafdata :
    ## determine gene range
    seqid = epaf['tstr']
    #seq = tarSeq[seqid]

    ##determint the gene type and size

    ## call A/B allele
    exon26keyseq = at.recoverTargetSeqFromQuery(
        exonSeq, epaf['cs'], abblock, epaf['strand'])
    intron9seq = at.recoverTargetSeqFromQuery(
        exonSeq, epaf['cs'], intron9block, epaf['strand'])

    genename, size, keytag = callC4allele(intron9seq, exon26keyseq)
    geneInfo[genename]['copyindex'] += 1
    source = geneInfo[genename]['source']
    copyindex = geneInfo[genename]['copyindex']
    #print(copyindex, file=sys.stderr)
    if copyindex > 1 :
        suffix = '.'+ str(copyindex) +'";' 
    else :
        suffix = '";' 
    attrpref = 'gene_id "' + geneInfo[genename]['geneid'] + suffix
    attrpref += ' transcript_id "' + geneInfo[genename]['transid']  + suffix
    attrpref += ' gene_name "' + genename + size + '";'


    bodyarr = []
    ## annotate exon
    exonArr = []
    #exonbodyarr = []
    for e in exonblock.split(',') :
        tag, rgs = e.split(":")
        exonwarning, rg = at.intervalMapQry2Ctg(
            rgs, epaf['cg'], epaf['strand'],
            int(epaf['qfro'])+1, epaf['qto'],
            int(epaf['tfro'])+1, epaf['tto'])
        if exonwarning == "imcomplete_coverage": 
            print(epaf)
            print(epaf['cg'])
            print(exonwarning, rgs, rg, epaf['strand'],
                  int(epaf['qfro'])+1, epaf['qto'],
                  int(epaf['tfro'])+1, epaf['tto'])
            break
        fro, to = rg.split("..")
        exonArr.append([int(fro), int(to)])
        attr = attrpref
        attr += 'exon_number "' + tag + '";'
        line = gtfline(seqid, source, 'exon', fro, to, '.', 
                       epaf['strand'], '.', attr)
        bodyarr.append(line + '\n')

    ## annotate start and stop codon
    for s in startcodonblock.split(',') :
        codonwarning, rg = at.intervalMapQry2Ctg(
            s, epaf['cg'], epaf['strand'],
            int(epaf['qfro'])+1, epaf['qto'],
            int(epaf['tfro'])+1, epaf['tto'])
        fro, to = rg.split('..')
        attr = attrpref
        line = gtfline(seqid, source, 'start_codon', fro, to, '.', 
                       epaf['strand'], '.', attr)
        bodyarr.append(line + "\n")

    stop_codon = []
    for s in stopcodonblock.split(',') :
        codonwarning, rg = at.intervalMapQry2Ctg(
            s, epaf['cg'], epaf['strand'],
            int(epaf['qfro'])+1, epaf['qto'],
            int(epaf['tfro'])+1, epaf['tto'])
        fro, to = rg.split('..')
        stop_codon.append(int(fro))
        stop_codon.append(int(to))
        attr = attrpref
        line = gtfline(seqid, source, 'stop_codon', fro, to, '.', 
                       epaf['strand'], '.', attr)
        bodyarr.append(line + "\n")

    ##transcript
    fro = epaf['tfro']
    to = epaf['tto']
    attr = attrpref
    attr += ' pipetide_key "' + at.translate(keytag) + '";'
    line = gtfline(seqid, source, 'transcript', str(int(fro) + 1), 
                   to, '.', epaf['strand'], '.', attr)
    transcriptbodystr = line + "\n"

    ## gene
    fro = epaf['tfro']
    to = epaf['tto']
    attr = 'gene_id "' + geneInfo[genename]['geneid']
    #print(copyindex, file=sys.stderr)
    if copyindex > 1 :
        attr += '.' + str(copyindex) + '";'
    else :
        attr += '";'
    attr += ' gene_name "' + genename + size + '";'
    line = gtfline(seqid, source, 'gene', str(int(fro) + 1), 
                   to, '.', epaf['strand'], '.', attr)
    genebodystr = line + "\n"

    ##CDS
    ### map new cds positions
    cdswarning, cdsrg = at.intervalMapQry2Ctg(
        cdsblock, epaf['cg'], epaf['strand'],
        int(epaf['qfro'])+1, epaf['qto'],
        int(epaf['tfro'])+1, epaf['tto'])
    cds = cdsrg.split(',')
    cdsdf = []
    cdsmin = 999999999
    cdsmax = 0
    for c in cds :
        line = [seqid, source]
        line.append('CDS')
        fro, to = c.split("..")
        #print(fro, to, file =sys.stderr)
        if stop_codon and removeStopCodonFromCDS :
            if epaf['strand'] == '+' and to == str(stop_codon[1]) :
                to = str(int(stop_codon[0]))
            if epaf['strand'] == '-' and fro == str(stop_codon[0]) :
                fro = str(int(stop_codon[1]))
        #print(stop_codon, file = sys.stderr)
        #print(fro, to, file = sys.stderr)
        line.append(fro)
        line.append(to)
        if cdsmin > int(fro) : cdsmin = int(fro)
        if cdsmax < int(to) : cdsmax = int(to)
        line.append('.')
        line.append(epaf['strand'])
        if epaf['strand'] == '+' and c == cds[0]:
            line.append(str(int(codonStart) -1))
        elif epaf['strand'] == '-' and c == cds[-1]:
            line.append(str(int(codonStart) -1))
        else :
            line.append('.')
        attr = attrpref
        line.append(attr)
        cdsdf.append(line)
        line = ''
    # update frame
    totalCDSlength = 0
    if epaf['strand'] == '+' :
        cdsdf = sorted(cdsdf, key = lambda x : int(x[3]))
    elif epaf['strand'] == '-' :
        cdsdf = sorted(cdsdf, key = lambda x : int(x[3]), reverse=True)
    ph = int(cdsdf[0][7])
    fro = int(cdsdf[0][3])
    to = int(cdsdf[0][4])
    phlast = (to - ph - fro + 1) % 3
    for i in range(1, len(cdsdf)) :
        ph = (3 - phlast) % 3
        cdsdf[i][7] = str(ph)
        fro = int(cdsdf[i][3])
        to = int(cdsdf[i][4])
        phlast = (to - ph - fro + 1) % 3

    ### write cds
    for x in cdsdf :
        # remove stop codon from CDS
        keepline = True
        if stop_codon and removeStopCodonFromCDS :
            a = min(stop_codon)
            b = max(stop_codon)
            if epaf['strand'] == '+' :
                if int(x[3]) < a and a <= int(x[4]) :
                    x[4] == str(a -1)
                elif int(x[3]) == a :
                    keepline = False
            if epaf['strand'] == '-' :
                if int(x[3]) <= b and b < int(x[4]) :
                    x[4] == str(b + 1)
                elif b == int(x[4]) :
                    keepline = False
        if keepline :
            line = "\t".join(x)
            bodyarr.append(line + "\n")
        # update cds range
        if cdsmin > int(x[3]) : cdsmin = int(x[3])
        if cdsmax < int(x[4]) : cdsmax = int(x[4])

    # define utr
    exonArr = sorted(exonArr, key = lambda x: x[0])
    i = 0
    utr = []
    for e in exonArr :
        if e[1] < cdsmin or e[0] > cdsmax:
            utr.append(e)
        elif e[0] < cdsmin and cdsmin < e[1] :
            utr.append([e[0], cdsmin-1])
        elif e[0] < cdsmax and cdsmax < e[1] :
            utr.append([cdsmax+1, e[1]])

    for u in utr :
        fro, to = u
        exonArr.append([int(fro), int(to)])
        attr = attrpref
        line = gtfline(seqid, source, 'UTR', fro, to, '.', 
                       epaf['strand'], '.', attr)

        bodyarr.append(line + "\n")

    bodystr = genebodystr + transcriptbodystr
    bodyarr = sorted(bodyarr, key=lambda x: int(x.split('\t')[4]))
    bodyarr = sorted(bodyarr, key=lambda x: int(x.split('\t')[3]))
    for line in bodyarr :
        bodystr += line
    sys.stdout.write(bodystr)
