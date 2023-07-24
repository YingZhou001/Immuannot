import sys,re,gzip,os
import paftools as pt
import pprint as pp
import algntools as at

paffile = sys.argv[1]
cdsseqfile = sys.argv[2]
callfile = sys.argv[3]

#search mutation

def extractCDS(cdsfa, keys) :
    buf = ''
    ks = keys.split(',')
    ret = []
    buftag = 0 # empty : 0 otherwise is 1
    with gzip.open(cdsfa, 'rt') as f:
        while True :
            line = f.readline()
            if not line: break
            if line[0] == '>' :
                if buftag == 1:
                    ret.append({'allele' : k,
                                'frame' : frame,
                                'cdsseq': buf})
                    buf = ''
                    ks.remove(k)
                    if len(ks) == 0:
                        break
                for k in ks :
                    buftag = 0
                    if k in line :
                        regex = re.compile(r'frame=(\d*)')
                        frame = regex.search(line)
                        if frame :
                            frame = frame.group(1)
                        buftag = 1
                        break
            else :
                if buftag == 1 :
                    buf += line.replace('\n', '')
        
    if buf :
        ret.append({'allele' : k,
                    'frame' : frame,
                    'cdsseq': buf})
    return(ret)


def nameNewHlaAllele(refAllele, cdsdiff) :
    gene, fields = refAllele.split('*')
    fields = fields.split(':')
    if not cdsdiff or cdsdiff == "NA" :
        fields[-1] = 'new'
    else :
        regex = re.compile(r'(\w+)\(\w*\)<(\w+)\(\w*\)')
        ret = regex.findall(cdsdiff)
        a = ''
        b = ''
        frameshift = False
        if ret :
            for x in ret :
                if len(x[1]) != len(x[0]):
                    frameshift = True
                a += x[0]
                b += x[1]
        else :
            sys.stderr.write(refAllele)
            sys.stderr.write(cdsdiff + '\n')
            sys.exit("Error in calling new allele : nameNewAllele\n")
        if a == b :
            #synonymous mutation
            if len(fields) >= 3 :
                fields = fields[0:3]
                fields[2] = 'new'
            else :
                fields[-1] = 'new'
        else :
            #non-synonymous mutation
            if len(fields) >= 2 :
                fields[1] = 'new'
                fields = fields[0:2]
            elif len(fields) <= 1 :
                fields = ['new']
        if frameshift :
            if len(fields) >= 2 :
                fields = fields[0:2]
                fields[1] = 'new'
            elif len(fields) == 1 :
                fields = ['new']
    ostr = gene + "*" + ":".join(fields)
    return(ostr)


def nameNewKirAllele(refAllele, cdsdiff) :
    gene, fields = refAllele.split('*')
    L = len(fields)
    if L >= 7 :
        fields = [fields[0:3], fields[3:5], fields[5:]]
    elif L >= 5 :
        fields = [fields[0:3], fields[3:]]
    elif L >= 3 :
        fields = [fields[0:]]
    else :
        sys.exit('failed to parse allele:' + refAllele + '\n')

    if not cdsdiff or cdsdiff == "NA" :
        fields[-1] = 'new'
    else :
        regex = re.compile(r'(\w+)\(\w*\)<(\w+)\(\w*\)')
        ret = regex.findall(cdsdiff)
        a = ''
        b = ''
        frameshift = False
        if ret :
            for x in ret :
                if len(x[1]) != len(x[0]):
                    frameshift = True
                a += x[0]
                b += x[1]
        else :
            sys.exit("Error in calling new allele : nameNewAllele\n")
        if a == b :
            #synonymous mutation
            if len(fields) >= 2 :
                fields = fields[0:2]
                fields[1] = 'new'
            elif len(fields) == 1 :
                fields[-1] = 'new'
        else :
            #non-synonimous mutation
            fields=['new']
        if frameshift :
            fields=['new']
    ostr = gene + "*" + "".join(fields)
    #sys.stderr.write(refAllele + ' ' + cdsdiff + '\n')
    #sys.exit()
    return(ostr)


def nameNewAllele(refAllele, cdsdiff) :
    if ':' in refAllele :
        return(nameNewHlaAllele(refAllele, cdsdiff))
    else :
        return(nameNewKirAllele(refAllele, cdsdiff))

def consensusCall(alleles) :
    common = os.path.commonprefix(alleles)
    if not '*' in common : 
        return("undetermined")
    else :
        gene, fields = common.split('*')
    if gene[0:3] == "KIR" :
        L = len(fields)
        if L >= 7 :
            fields = fields[0:7]
        elif L >= 5 :
            fields = fields[0:5]
        elif L >= 3 :
            fields = fields[0:3]
        else :
            fields = ''
        return(gene + '*' + fields + "new")
    else :
        fields = fields.split(":")
        fields[-1] = "new"
        return(gene + '*' + ':'.join(fields))



def callNewAllele(cdspaf, cdsseqfile) :
    # features of the best match
    cdspaf = sorted(cdspaf, key=lambda x: int(x['match']), reverse=True)
    cdspaf = sorted(cdspaf, key=lambda x: int(x['NM']))

    nm = -1
    match = 0
    for x in cdspaf :
        refallele = x['qstr']
        if refallele[-1] != 'N' :
            nm = x['NM']
            match = x['match']
            break
    if nm == -1 :
        nm = cdspaf[0]['NM']
        match = cdspaf[0]['match']

    cands = []
    cs = []
    for x in cdspaf :
        if x['NM'] <= nm and x['match'] >= match :
            cands.append(x['qstr'])
            cs.append(x['cs'])
        else :
            break
    
    newallele = []
    mut = []
    nonsense = 0
    if int(nm) == 0 :
        for y in cands :
            tmp = nameNewAllele(y, "NA")
            tmptag = 0
            for x in newallele :
                if tmp.replace("new", '') in x :
                    tmptag = 1
                    break
            if tmptag == 0 :
                newallele.append(tmp)
            i = cands.index(y)
            mut.append("|".join([y, cs[i], '']))
    elif int(nm) > 0 :
        cdsseq=extractCDS(cdsseqfile, ",".join(cands))
        for y in cdsseq :
            i = cands.index(y['allele'])
            dfs = at.findCodingDiff(y['cdsseq'], cs[i], int(y['frame'])-1, '+')
            tmp = nameNewAllele(y['allele'], dfs)
            tmptag = 0
            for x in newallele :
                if tmp.replace("new", '') in x :
                    tmptag = 1
                    break
            if tmptag == 0 :
                newallele.append(tmp)
            mut.append("|".join([y['allele'], cs[i], dfs]))
                #pp.pprint(cdsseq)
                #cdsseq = ",".join(fd.loadfasta(qryfile, cands))

    if len(newallele) > 1 :
        newallele = [consensusCall(newallele)]
        #newallele=['"multiple calls"']
    return({'consensus' : ",".join(newallele),
            'NM' : nm,
            'call' : ",".join(cands),
            'mut' : ",".join(mut)})

# main program

# read cds paf records
pafda = []
keybags = ["qstr", 'qlen', 'qfro', 'qto', 'strand', 'tfro', 'tto', 
           'tlen', 'tstr', 'match', 'NM', 'cg', 'cs']
with gzip.open(paffile, 'rt') as f:
    while True:
        line = f.readline()
        newdct = {}
        if not line: break
        else : 
            line = line.replace('\n', '')
            newdct = pt.parsingPaf(line)
        if newdct:
            c1 = int(newdct['qto']) - int(newdct['qfro']) == int(newdct['qlen'])
            c2 = int(newdct['tto']) - int(newdct['tfro']) == int(newdct['tlen'])
            c3 = int(newdct['qlen']) / int(newdct['tlen']) > 0.95
            c4 = newdct['NM'] == '0'
            if (c1 or c2) and (c3 or ((not c3) and c4)):
                # for partial sample, we only allow nm = 0
                x = { k: v for k, v in newdct.items() if k in keybags }
                pafda.append(x)

with open(callfile, 'r') as f:
    while True:
        line = f.readline()
        line =line.replace('\n', '')
        if not line: break
        llst = line.split("\t")
        ctgid = llst[0]
        gene = llst[1]
        strand = llst[2]
        template = llst[3]
        allelecall = llst[4]
        template_mapping = llst[5]
        cds_mapping = 'NA'

        if allelecall != "NA" or template == "NA":
            ostr = [ctgid, gene, strand, template, allelecall,
                    "template_mapping " + '"' + template_mapping + '"',
                    "cds_mapping " + '"' + cds_mapping  + '"']
            ostr =  "\t".join(ostr)
            print(ostr)
            continue

        spaf = []
        for x in pafda :
            k =  ctgid + '_'+ gene
            if k == x['tstr'][0:len(k)] :
                spaf.append(x)

        if spaf :
            newcall = callNewAllele(spaf, cdsseqfile)
            #update gene name based on the allele calls
            ## for undetermined allele, use the best template match
            if '*' in newcall['consensus'] :
                gene = newcall['consensus'].split('*')[0]
            
            allelecalls = newcall['call']
            cds_mapping = "consensus=" + newcall['consensus'] + ";"
            cds_mapping += "nm=" + newcall['NM'] + ";"
            cds_mapping += "mut=" + newcall['mut']  + ';"'

        ostr = [ctgid, gene, strand, template, allelecalls,
                "template_mapping " + '"' + template_mapping + '"',
                "cds_mapping " + '"' + cds_mapping  + '"']
        ostr =  "\t".join(ostr)
        print(ostr)


