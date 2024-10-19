import sys, re, gzip
import paftools as pt
import pprint as pp
import fastatools as ft
import algntools as at

#gene, geneID, transID = sys.argv[1].split(",")
paffile = sys.argv[1]
allelefile = sys.argv[2]
qctgfile = sys.argv[3]
outpref = sys.argv[4]
overlapcut = sys.argv[5]
diffcut = sys.argv[6]

outfile1 = outpref+"/tmp.gene.csv"
outfile2 = outpref+"/cds.fa.gz"
outfile3 = outpref+"/tmp.newgene.txt"
outfile4 = outpref+"/gene.filtered.paf"

file1 = open(outfile1, "w")
file2 = gzip.open(outfile2, "wt")
file3 = open(outfile3, "w")
file4 = open(outfile4, "w")

def readAllele(file):
    # read allele element
    alleleda = []
    alleles = []
    keybags = ['allele', 'codonShift', 'CDS']
    with gzip.open(file, 'rt') as f:
        while True:
            line = f.readline()
            if not line: break
            newdct = {}
            llst = line.split()
            for x in llst :
                k, v = x.split('=')
                if k in keybags :
                    newdct[k] = v
            alleles.append(newdct['allele'])
            alleleda.append(newdct)
    return([alleles, alleleda])


def checkGeneCopies(pafda) :
    rgs = []
    if not pafda :
        return(0)
    for da in pafda :
        qfro = int(da['qfro'])
        qto = int(da['qto'])
        if rgs :
            for fro,to in rgs :
                if qto < fro or qfro > qto :
                    rgs.append([qfro, qto])
        else :
            rgs.append([qfro, qto])
    return(len(rgs))

def checkCDScompleteness(cdsstr) :
    ret = []
    # check length
    L = len(cdsstr)
    if L % 3 != 0 :
        ret.append("partial_CDS")
    else :
        n = L // 3
        for i in range(n):
            codon = cdsstr[(i*3-1):i*3]
            if codon in ['TAG', 'TAA', 'TGA'] and i < n-1 :
                ret.append("inframe_stop")
                break
        if cdsstr[0:3] != 'ATG' :
            ret.append("no-start_codon")
        if not cdsstr[(L-3):L] in ['TAG', 'TAA', 'TGA'] :
            ret.append("no-stop_codon")
    return(",".join(ret))


def mergeInterval(a, b, cutoff) :
    # if intervals a and b overlap over the cutoff, a and b will be merged
    # into one interval, otherwise return a
    if not a :
        return([1, b])
    x0, y0 = a.split(',')
    x1, y1 = b.split(',')
    x0 = int(x0)
    y0 = int(y0)
    if y0 < x0 :
        tmp = x0
        x0 = y0
        y0 = tmp
    x1 = int(x1)
    y1 = int(y1)
    if y1 < x1 :
        tmp = x1
        x1 = y1
        y1 = tmp
    if y0 < x1 or x0 > y1 :
        return([1, b])
    else :
        l0 = y0 - x0
        l1 = y1 - x1
        l2 = abs(y1 - x0)
        minl = min(l0, l1)
        r = (l0+l1-l2)/minl/2
        if r > float(cutoff) :
            return([0, str(min(x0,x1)) + ',' + str(max(y0,y1))])
        else :
            return([0, a])


def detectSplitCtg(strand, qlen, qfro, qto, tlen, tfro, tto) :
    # assuming that target the is contig to detect
    ql = int(qlen)
    qf = int(qfro)
    qt = int(qto)
    tl = int(tlen)
    tf = int(tfro)
    tt = int(tto)

    if strand == '+' :
        c1 = qf > tf and tf < 3
        c2 = ((ql - qt) > (tl - tt) and (tl - tt) < 3)
    else :
        c1 = qf > (tl - tf) and (tl - tf) < 3
        c2 = ((ql - qt) > tt and tt < 3)
    if c1 or c2 :
        return(1)
    else :
        return(0)


def searchTemplatePGPC(pafDA, ctgname, qctgseqs) :
    ret = []
    rec = {}
    # search template per gene per contig
    pafda = []
    for da in pafDA :
        if da['tstr'] == ctgname :
            pafda.append(da)
    if not pafda :
        rec['attr'] = "no evidence"
        rec['ctgname'] = ctgname
        rec['gene'] = "NA"
        rec['template'] = "NA"
        rec['allele'] = "NA"
        rec['strand'] = "NA"
        ret.append(rec)
        return(ret)
    else :
        # decresing order
        pafda = sorted(pafda, key=lambda x: int(x['match']), reverse=True)
        # increasing order
        pafda = sorted(pafda, key=lambda x: int(x['adjNM']))

    # search mapping clusters when multiple gene mappings exist on the same
    ## contig
    clusters = [pafda[0]['tfro'] + ',' + pafda[0]['tto']]

    for da in pafda :
        newtag = 1
        testc = da['tfro'] + ',' + da['tto']
        clustMatch = []
        for c in clusters :
            tag, newc = mergeInterval(c, testc, 0.1)
            if tag == 0 :
                # record the matched cluster and update
                clustMatch.append([c, newc])
                if len(clustMatch) > 1 :
                    # multiple matches, abandon this mapping
                    newtag = 0
                    break
        if len(clustMatch) == 1 :
            c = clustMatch[0][0]
            newc = clustMatch[0][1]
            clusters[clusters.index(c)] = newc
            newtag = 0
        if newtag == 1:
            clusters.append(testc)
        #print(clusters)

    # load gene structure information
    alleles, alleleda = readAllele(allelefile)


    #print(clusters)
    for clust in clusters :
        subpafda = []
        for da in pafda :
            #subset paf record for each cluster
            testc = da['tfro'] + ',' + da['tto']
            tag, c = mergeInterval(clust, testc, 0.5)
            if tag == 0 :
                subpafda.append(da)

        rec = {}
        rec['ctgname'] = ctgname
        rec['allele'] = "NA"
        templateTag = 0
        # find allele for the perfect match:
        for da in subpafda :
            c1 = float(da['adjNM']) == 0
            c2 = int(da['qlen']) == int(da['qto']) - int(da['qfro'])
            """
            ## detect split contig
            c3 = 1 == detectSplitCtg(
                da['strand'], da['qlen'], da['qfro'], 
                da['qto'], da['tlen'], da['tfro'], da['tto'])
            """
            if c1 and c2 :
                rec['template'] = da['qstr']
                templateTag = 1
                rec['allele'] = da['qstr']
                rec['gene'] = da['qstr'].split('*')[0]
                rec['strand'] = da['strand']
                dct = alleleda[alleles.index(da['qstr'])]
                qcds = dct['CDS']
                frame = dct['codonShift']
                err, ccds = at.intervalMapQry2Ctg(
                    qcds, da['cg'], da['strand'], 
                    int(da['qfro'])+1, da['qto'], 
                    int(da['tfro'])+1, da['tto'])
                break

        # find template for imperfect match and extract cds string:
        if templateTag == 0 :
            # best case: find template with the equal CDS length
            for da in subpafda :
                err = ''
                dct = alleleda[alleles.index(da['qstr'])]
                qcds = dct['CDS']
                frame = dct['codonShift']
                err, ccds = at.intervalMapQry2Ctg(
                    qcds, da['cg'], da['strand'], 
                    int(da['qfro'])+1, da['qto'], 
                    int(da['tfro'])+1, da['tto'])
                #sys.stderr.write(">>" + da['qstr'] + ">>" + err + "\n")
                if not err :
                    templateTag = 1
                    break

        # moderate case: find template without frameshift mutation in CDS
        if templateTag == 0 :
            for da in subpafda :
                err = ''
                dct = alleleda[alleles.index(da['qstr'])]
                qcds = dct['CDS']
                frame = dct['codonShift']
                err, ccds = at.intervalMapQry2Ctg(
                    qcds, da['cg'], da['strand'], 
                    int(da['qfro'])+1, da['qto'], 
                    int(da['tfro'])+1, da['tto'])
                #sys.stderr.write(">>" + da['qstr'] + ">>" + err + "\n")
                if err == 'unequal_length':
                    templateTag = 1
                    break

        # worst case: find template that the CDS length may not equal
        if templateTag == 0 :
            err = ''
            da = subpafda[0]
            dct = alleleda[alleles.index(da['qstr'])]
            qcds = dct['CDS']
            frame = dct['codonShift']
            err, ccds = at.intervalMapQry2Ctg(
                qcds, da['cg'], da['strand'], 
                int(da['qfro'])+1, da['qto'], 
                int(da['tfro'])+1, da['tto'])
            #sys.stderr.write(">>>this category" + '\n')

        #construct target cds seq from query seq
        qseq = qctgseqs[da['qstr']]
        fro = int(da['qfro'])
        to = int(da['qto'])
        qseq = qseq[fro: to]
        csstr = da['cs']
        cdsstr = ""
        intronBound = ""
        """
        sys.stderr.write('*******\n')
        sys.stderr.write('templateTag:' + str(templateTag) +'\n')
        sys.stderr.write('query seq:' + da['qstr'] +'\n')
        sys.stderr.write(da['qfro'] + '..' + da['qto'] + '\n')
        sys.stderr.write('target seq:' + da['tstr'] +'\n')
        sys.stderr.write(da['tfro'] + '..' + da['tto'] + '\n')
        sys.stderr.write('qcds' + qcds +'\n')
        sys.stderr.write('ccds' + ccds +'\n')
        sys.stderr.write('*******\n')
        """
        for cds in qcds.split(",") :
            a,b = cds.split("..")
            a = int(a)
            b = int(b)
            if b > to or a < fro:
                sys.stderr.write("Warning, CDS seq may be incomplete\n"
                         + "query seq:" + da['qstr'] + '\n'
                         + "target seq:" + da['tstr'] + '\n')
                continue
            a = a - fro
            b = b - fro
            s0 = qseq[a-1:b]
            cdsrg = str(a) + '..' + str(b)
            # will recover positive or the same strand diretion as qseq
            s1 = at.recoverTargetSeqFromQuery(qseq, csstr, cdsrg, da['strand'])
            if s1 :
                #sys.stderr.write('*******\n')
                cdsstr += s1
                tmprg = str(a-2) + '..' + str(b + 2)
                s2 = at.recoverTargetSeqFromQuery(qseq, csstr, tmprg, da['strand'])
                """
                #sys.stderr.write('qstr ' + da['qstr'] + '_\n')
                #sys.stderr.write('tstr ' + da['tstr'] + '_\n')
                #sys.stderr.write('cdsrg_' + cds + '_\n')
                #sys.stderr.write('cdsrg_' + cdsrg + '_\n')
                #sys.stderr.write('cs_' + csstr + '_\n')
                #sys.stderr.write('s0_' + s0 + '_\n')
                #sys.stderr.write('s1_' + s1 + '_\n')
                #sys.stderr.write('s2_' + s2 + '_\n')
                """
                left = s2[0:2]
                right = s2[-2] + s2[-1]
                intronBound += left +  ":"  + right + "_"

        rec['template'] = da['qstr']
        rec['strand'] = da['strand']
        rec['gene'] = da['qstr'].split('*')[0]
        tmp = intronBound.split(":")[1:-1]
        intronBound = ":".join(tmp)
        rec['cdsstr'] = cdsstr
        rec['intronBound'] = intronBound
        attr = ''
        attr += "nm=" + da['adjNM'] + ';'
        attr += "match=" + da['match'] + ';'
        attr += "qlen=" + da['qlen'] + ";"
        attr += "qrg=" + da['qfro'] + ".." + da['qto'] + ";"
        attr += "trg=" + da['tfro'] + ".." + da['tto'] + ";"
        attr += "cg=" + da['cg'] + ";"
        cdserr = checkCDScompleteness(cdsstr)
        #sys.stderr.write(rec['gene'] + '\n')
        #sys.stderr.write(cdserr + '\n')
        if cdserr :
            attr += "warning=" + cdserr + ";"
        rec['attr'] = attr
        ret.append(rec)
    return(ret)




## main program

# read paf records
pafDA = []
keytags = ["qstr", 'qlen', 'qfro', 'qto', 'strand', 'tfro', 'tto', 
           'tlen', 'tstr', 'match', 'adjNM', 'cg', 'cs']
ctgids = []
qctgids = []
sys.stderr.write("<msg> loading paf records\n")
with gzip.open(paffile, 'rt') as f:
    while True:
        line = f.readline()
        if not line: break
        newdct = pt.parsingPaf(line)
        newdct['adjNM'] = pt.adjNM(newdct['cs'])
        if newdct:
            # query length should be shorter than the target length, otherwise
            # to be hardly annotated
            c0 = int(newdct['qlen']) <= int(newdct['tlen'])
            c1 = int(newdct['match']) / int(newdct['qlen']) > float(overlapcut)
            # to avoid that a mapping cover two gene regions
            qmlen = int(newdct['qto']) - int(newdct['qfro'])
            tmlen = int(newdct['tto']) - int(newdct['tfro'])
            c2 = qmlen / tmlen > 0.5
            c3 = int(newdct['adjNM']) <= 3
            c4 = int(newdct['adjNM']) / int(newdct['match']) < float(diffcut)
            if c0 and c1 and c2 and (c3 or c4) :
                x = { k: v for k, v in newdct.items() if k in keytags }
                pafDA.append(x)
                ctg = x['tstr']
                qctg = x['qstr']
                if not ctg in ctgids :
                    ctgids.append(ctg)
                if not qctg in qctgids :
                    qctgids.append(qctg)
            else :
                file4.write(line)
file4.close()

#sys.stderr.write(str(len(pafDA)) + '\n')
#load contigs
sys.stderr.write("<msg> loading contig seq\n")
qctgids = ",".join(qctgids)
qctgseqs = ft.loadFastagz(qctgfile, qctgids)

#pp.pprint(pafda)
sys.stderr.write("<msg> searching template: \n")
newgenes = []

# kir gene patches
## this is because kir genes share the same structures
## see https://www.ebi.ac.uk/ipd/kir/assets/images/about_figure_02-8e7817f508.png
kir_gene_patchs = {
    'KIR2DL1' : ['KIR2DL1', 'KIR2DL2', 'KIR2DL3'],
    'KIR2DL2' : ['KIR2DL1', 'KIR2DL2', 'KIR2DL3'],
    'KIR2DL3' : ['KIR2DL1', 'KIR2DL2', 'KIR2DL3'],
    'KIR2DS1' : ['KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5'],
    'KIR2DS2' : ['KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5'],
    'KIR2DS3' : ['KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5'],
    'KIR2DS4' : ['KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5'],
    'KIR2DS5' : ['KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5'],
    'KIR3DL1' : ['KIR3DL1', 'KIR3DL2'],
    'KIR3DL2' : ['KIR3DL1', 'KIR3DL2'],
    'KIR2DL4' : ['KIR2DL4', 'KIR2DL5A', 'KIR2DL5B'],
    'KIR2DL5A' : ['KIR2DL4', 'KIR2DL5A', 'KIR2DL5B'],
    'KIR2DL5B' : ['KIR2DL4', 'KIR2DL5A', 'KIR2DL5B'],
    'KIR3DS1' : ['KIR3DS1'],
    'KIR3DL3' : ['KIR3DL3'],
    'KIR3DP1' : ['KIR3DP1'],
    'KIR2DP1' : ['KIR2DP1'],
}


for ctgname in ctgids :
    #sys.stderr.write("<msg>    gene:" + gene + ", contig:" + ctgname + '\n')
    out = searchTemplatePGPC(pafDA, ctgname, qctgseqs)
    #write out
    ## CDS seq
    for i in range(len(out)) :
        x = out[i]
        gene = x['gene']
        if "cdsstr" in x.keys() :
            file2.write(">" + ctgname + '_' + gene + '_' + str(i+1) 
                    + " " + x['intronBound'] + "\n")
            file2.write(x['cdsstr'] + "\n")
            ###
            if not gene in newgenes :
                if gene in kir_gene_patchs.keys() :
                    for kir_gene in kir_gene_patchs[gene] :
                        newgenes.append(kir_gene)
                        file3.write(kir_gene + "\n")
                else :
                    newgenes.append(gene)
                    file3.write(gene + "\n")

        if x['template'] == "NA" : continue
        ostr = ctgname
        ostr += "\t" + gene
        ostr += "\t" + x['strand']
        ostr += "\t" + x['template']
        ostr += "\t" + x['allele']
        ostr += "\t" + x['attr']
        file1.write(ostr + "\n")

file1.close()
file2.close()
file3.close()
