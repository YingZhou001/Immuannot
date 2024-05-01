import sys,re,gzip
import pprint as pp
import algntools as at
from datetime import date


callfile = sys.argv[1]
allelefile = sys.argv[2]
geneinfofile = sys.argv[3]
removeStopCodonFromCDS = True

def geneStructMapping(templateStr, referenceRec, strand) :
    # extract information from templatestr
    tlst = templateStr.replace('"', '').split(";")
    tmp = {}
    for x in tlst :
        if x:
            k, v = x.split('=')
            tmp[k] = v
    cg = tmp['cg']
    qrg = tmp['qrg']
    trg = tmp['trg']
    # input from minimap2, which is zero based coordinate
    qfro, qto = qrg.split("..") 
    qfro = str(int(qfro) + 1)
    tfro, tto = trg.split("..")
    tfro = str(int(tfro) + 1)

    utr = referenceRec['UTR']
    cds = referenceRec['CDS']
    exon = referenceRec['exon']
    startCodon = referenceRec['startCodon']
    stopCodon = referenceRec['stopCodon']
    frame = int(referenceRec['codonShift']) - 1
    warnings = []

    ret = {}
    ret['gene'] = trg
    ret['strand'] = strand
    ret['frame'] = str(frame)
    if utr and utr != "NA" :
        msg, ret['utr'] = at.intervalMapQry2Ctg(
            utr, cg, strand, qfro, qto, tfro, tto)
        if msg : warnings.append("utr:" + msg)
    msg, ret['cds'] = at.intervalMapQry2Ctg(
        cds, cg, strand, qfro, qto, tfro, tto)
    if msg : warnings.append("cds:" + msg)

    newExon = []
    for x in exon.split(",") :
        ind, rg = x.split(":")
        msg, rg = at.intervalMapQry2Ctg(rg, cg, strand, qfro, qto, tfro, tto)
        if msg : warnings.append("exon" + ind + ":" + msg)
        newExon.append(ind + ":" + rg)
    ret['exon'] = ",".join(newExon)

    p,rg = startCodon.split(":")
    msg, rg = at.intervalMapQry2Ctg(rg, cg, strand, qfro, qto, tfro, tto)
    if msg : warnings.append("startCodon:" + msg)
    ret['startCodon'] = p + ":" + rg
    p,rg = stopCodon.split(":")
    msg, rg = at.intervalMapQry2Ctg(rg, cg, strand, qfro, qto, tfro, tto)
    if msg : warnings.append("stopCodon:" + msg)
    ret['stopCodon'] = p + ":" + rg

    ret['warnings'] = ",".join(warnings)
    return(ret)

geneInfo = {}
with gzip.open(geneinfofile, 'rt') as f :
    while True:
        line = f.readline()
        if not line : break
        line = line.replace('\n', '')
        llst = line.split(',')
        newdict = {}
        newdict['source'] = llst[1]
        newdict['geneid'] = llst[2]
        newdict['transid'] = llst[3]
        newdict['copyindex'] = 0
        geneInfo[llst[0]] = newdict


callData = []
tmpltAlleles = []
with open(callfile, 'r') as f :
    while True :
        line = f.readline()
        if not line: break
        line = line.replace('\n', '')
        llst = line.split("\t")
        x= {}
        #sys.stderr.write(str(llst)+'\n')
        x['sample'] = llst[0]
        x['gene'] = llst[1]
        x['strand'] = llst[2]
        x['template'] = llst[3]
        x['alleles'] = llst[4]
        #sys.stderr.write(line + '\n')
        k,v = llst[5].split(" ", 1)
        x[k] = v
        k,v = llst[6].split(" ", 1)
        x[k] = v
        callData.append(x)
        tmpltAlleles.append(x['template'])

#pp.pprint(callData, stream=sys.stderr)
#sys.exit()

alleleData = []
keybag = ['allele', 'codonShift', 'startCodon', 'stopCodon', 
          'CDS', 'UTR', 'exon', 'gene', 'type']
with gzip.open(allelefile, 'rt') as f :
    while len(tmpltAlleles) > 0 :
        line = f.readline()
        if not line: break
        tag = 0
        for x in tmpltAlleles :
            if x in line :
                tag = 1
                tmpltAlleles.remove(x)
                break
        if tag == 1 :
            llst = line.split("\t")
            y = {}
            for i in range(len(llst)) :
                k, v = llst[i].split("=")
                if k in keybag :
                    y[k] = v
            alleleData.append(y)

#pp.pprint(alleleData)

#headstr = ''
bodystr = ''

alleles = []
for x in alleleData :
    alleles.append(x['allele'])

for da in callData :
    tallele = da['template']
    rec = alleleData[alleles.index(tallele)]

    gstruct = geneStructMapping(da['template_mapping'], rec, da['strand']) 

    # output gtf format
    genename = da['gene']
    #sys.stderr.write('->>>>' + genename + '\n')
    pref = [da['sample'], geneInfo[genename]['source']]
    geneInfo[genename]['copyindex'] += 1
    copyindex = geneInfo[genename]['copyindex']
    if copyindex > 1 :
        suffix = '.'+ str(copyindex) +'";' 
    else :
        suffix = '";' 
    attrpref = 'gene_id "' + geneInfo[genename]['geneid']  + suffix
    attrpref += ' transcript_id "' + geneInfo[genename]['transid'] + suffix
    attrpref += ' ' + 'gene_name "' + genename + '";'

    ## gene
    line = pref.copy()
    line.append('gene')
    fro, to = gstruct['gene'].split("..")
    line.append(str(int(fro) + 1))
    line.append(to)
    line.append('.')
    line.append(gstruct['strand'])
    line.append('.')
    attr = 'gene_id "' + geneInfo[genename]['geneid'] + suffix
    attr += ' template_allele "' + tallele + '";'
    ts1, ts2=da['template_mapping'].replace('"', '').split(';', 1)
    ts1 = ts1.split('=')[1]
    attr += ' template_distance ' + ts1 + ';'
    #if gstruct['warnings'] :
    #    attr += ' ' + 'template_warning "' + gstruct['warnings'] + '";'
    attr += ' gene_name "' + genename + '";' 
    line.append(attr)
    line = "\t".join(line)

    bodystr += line + "\n"

    ##transcript
    line = pref.copy()
    line.append('transcript')
    fro, to = gstruct['gene'].split("..")
    line.append(str(int(fro) + 1))
    line.append(to) # zero based coordinate to one based coordinate
    line.append('.')
    line.append(gstruct['strand'])
    line.append('.')
    attr = attrpref

    ##only keep warning for ts2
    regex = re.compile(r'warning=(.*?);')
    ret = regex.search(ts2)
    if ret :
        ts2 = ret.group(1)
    else :
        ts2 = 'NA'

    if da['cds_mapping'] != '"NA"' : 
        cs0, cs1, cs2=da['cds_mapping'].replace('"', '').split(';', 2)
        cs0 = cs0.split('=')[1]
        cs1 = cs1.split('=')[1]
        cs2 = cs2.split('=')[1]
        attr += ' consensus "' + cs0 + '";'
        attr += ' alleles "' + da['alleles'] + '";'
        if ts2:
            attr += ' template_warning "' + ts2 + '";'
        attr += ' cds_distance ' + cs1 + ';'
        if int(cs1) > 0:
            attr += ' cds_mut "' + cs2 + '";'
    else :
        attr += ' consensus "' + da['alleles'] + '";'
        attr += ' alleles "' + da['alleles'] + '";'
        if ts2:
            attr += ' template_warning "' + ts2 + '";'
        #attr += ' cds_mapping ' + da['cds_mapping'] + ';'
    line.append(attr)
    line = "\t".join(line)

    bodystr += line + "\n"
    
    bodyarr = []
    ## utr
    if 'utr' in gstruct.keys() :
        utr = gstruct['utr']
        for u in utr.split(',') :
            line = pref.copy()
            line.append('UTR')
            fro, to = u.split("..")
            # remove elements that not well mapped
            if int(fro) < 0 or int(to) < 0 : continue
            line.append(fro)
            line.append(to)
            line.append('.')
            line.append(gstruct['strand'])
            line.append('.')
            attr = attrpref
            line.append(attr)
            line = "\t".join(line)
            bodyarr.append(line + "\n")

    if 'exon' in gstruct.keys() :
        exon = gstruct['exon']
        for c in exon.split(',') :
            line = pref.copy()
            line.append('exon')
            tag, rg = c.split(":")
            fro, to = rg.split("..")
            # remove elements that not well mapped
            if int(fro) < 0 or int(to) < 0 : continue
            line.append(fro)
            line.append(to)
            line.append('.')
            line.append(gstruct['strand'])
            line.append('.')
            attr = attrpref
            attr += ' exon_number "' + tag + '";'
            line.append(attr)
            line = "\t".join(line)
            bodyarr.append(line + "\n")

    if 'startCodon' in gstruct.keys() :
        sc = gstruct['startCodon']
        line = pref.copy()
        line.append('start_codon')
        tag, rg = sc.split(":")
        fro, to = rg.split("..")
        # remove elements that not well mapped
        if int(fro) >= 0 and int(to) >= 0 :
            line.append(fro)
            line.append(to)
            line.append('.')
            line.append(gstruct['strand'])
            line.append('0')
            attr = attrpref
            attr += ' codon "' + tag + '";'
            line.append(attr)
            line = "\t".join(line)
            if tag == "ATG" :
                bodyarr.append(line + "\n")


    stopCodon = []
    if 'stopCodon' in gstruct.keys() :
        sc = gstruct['stopCodon']
        #print(genename, sc, file=sys.stderr)
        line = pref.copy()
        line.append('stop_codon')
        tag, rg = sc.split(":")
        fro, to = rg.split("..")
        # remove elements that not well mapped
        if int(fro) >= 0 and int(to) >= 0 :
            line.append(fro)
            line.append(to)
            stopCodon = [fro, to]
            line.append('.')
            line.append(gstruct['strand'])
            line.append('0')
            attr = attrpref
            attr += ' codon "' + tag + '";'
            line.append(attr)
            line = "\t".join(line)
            if tag in ["TGA", "TAG", "TAA"] :
                bodyarr.append(line + "\n")


    if 'cds' in gstruct.keys() :
        cdsdf = []
        cds = gstruct['cds']
        cds = cds.split(',')
        for c in cds :
            line = pref.copy()
            line.append('CDS')
            fro, to = c.split("..")
            # remove elements that not well mapped
            if int(fro) < 0 or int(to) < 0 : continue

            # remove stop codon from CDS
            if stopCodon and removeStopCodonFromCDS :
                if gstruct['strand'] == '+' and to == stopCodon[1] :
                    to = str(int(stopCodon[0])-1)
                if gstruct['strand'] == '-' and fro == stopCodon[0] :
                    fro = str(int(stopCodon[1])+1)
            line.append(fro)
            line.append(to)
            line.append('.')
            line.append(gstruct['strand'])
            if gstruct['strand'] == '+' and c == cds[0]:
                line.append(gstruct['frame'])
            elif gstruct['strand'] == '-' and c == cds[-1]:
                line.append(gstruct['frame'])
            else :
                line.append('.')
            attr = attrpref
            line.append(attr)
            cdsdf.append(line)
            line = ''

        # update frame
        if gstruct['strand'] == '+' :
            cdsdf = sorted(cdsdf, key = lambda x : int(x[3]))
        elif gstruct['strand'] == '-' :
            cdsdf = sorted(cdsdf, key = lambda x : int(x[3]), reverse=True)

        #sys.stderr.write(str(cdsdf[0]) + '\n')
        if cdsdf[0][7] != '.' :
            ph = int(cdsdf[0][7])
            fro = int(cdsdf[0][3])
            to = int(cdsdf[0][4])
            phlast = (to - ph - fro + 1) % 3
            for i in range(1, len(cdsdf)) :
                ph = (3 - phlast) % 3
                fro = int(cdsdf[i][3])
                to = int(cdsdf[i][4])
                cdsdf[i][7] = str(ph)
                phlast = (to - ph - fro + 1) % 3

        for x in cdsdf :
            #print(x, file=sys.stderr)
            line = "\t".join(x)
            bodyarr.append(line + "\n")

    bodyarr = sorted(bodyarr, key=lambda x: int(x.split('\t')[4]))
    bodyarr = sorted(bodyarr, key=lambda x: int(x.split('\t')[3]))
    for line in bodyarr :
        bodystr += line

sys.stdout.write(bodystr)
