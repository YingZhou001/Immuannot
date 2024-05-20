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
    ret = {}
    buftag = 0 # empty : 0 otherwise is 1
    with gzip.open(cdsfa, 'rt') as f:
        while True :
            line = f.readline()
            if not line: break
            if line[0] == '>' :
                if buftag == 1:
                    ret[k] = {'frame' : frame, 'cdsseq': buf}
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
        ret[k] = {'frame' : frame, 'cdsseq': buf}
    return(ret)


def nameNewHlaAllele(refAllele, cdsdiff) :
    gene, fields = refAllele.split('*')
    fields = fields.split(':')
    # append '00' to the empty field
    while len(fields) < 4 : fields.append('00') 
    if not cdsdiff or cdsdiff == "NA" :
        fields[3] = 'new'
    else :
        regex = re.compile(r'(\w*)\(\w*\)<(\w*)\(\w*\)')
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
            fields = fields[0:3]
            fields[2] = 'new'
        else :
            #non-synonymous mutation
            fields = fields[0:2]
            fields[1] = 'new'
        if frameshift :
            fields = fields[0:2]
            fields[1] = 'new'
    ostr = gene + "*" + ":".join(fields)
    return(ostr)


def nameNewKirAllele(refAllele, cdsdiff) :
    gene, fields = refAllele.split('*')
    while len(fields) < 7 : fields += '0' 
    fields = [fields[0:3], fields[3:5], fields[5:]]

    if not cdsdiff or cdsdiff == "NA" :
        fields[2] = 'new'
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
            fields = fields[0:2]
            fields[1] = 'new'
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
    if 'KIR' in refAllele :
        return(nameNewKirAllele(refAllele, cdsdiff))
    else :
        return(nameNewHlaAllele(refAllele, cdsdiff))


def consensusCall(alleles) :

    # check gene names
    genes = []
    uniq_alleles = []
    for allele in alleles :
        gene = allele.split('*')[0]
        if gene not in genes :
            genes.append(gene)
        if allele not in uniq_alleles :
            uniq_alleles.append(allele)

    if len(genes) > 1 :
        return("undetermined")

    if len(uniq_alleles) == 1 :
        return(uniq_alleles[0])

    common = os.path.commonprefix(uniq_alleles)
    gene, fields = common.split('*')
    if gene[0:3] == "KIR" :
        L = len(fields)
        if L >= 5 :
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

    # determine the cutoff based on expression allele and non-expression allele
    nm_cutoff = -1
    match_cutoff = 0
    for x in cdspaf :
        refallele = x['qstr']
        if refallele[-1] != 'N' :
            nm_cutoff = x['NM']
            match_cutoff = x['match']
            break
    if nm_cutoff == -1 :
        nm_cutoff = cdspaf[0]['NM']
        match_cutoff = cdspaf[0]['match']

    cands = []
    cs = []
    qfros = []
    for x in cdspaf :
        if x['NM'] <= nm_cutoff and x['match'] >= match_cutoff :
            cands.append(x['qstr'])
            cs.append(x['cs'])
            qfros.append(x['qfro'])
        else :
            break

    newalleles = []
    muts = []

    cands_len = len(cands)
    if int(nm_cutoff) == 0 :
        for i in range(cands_len) :
            newalleles.append(nameNewAllele(cands[i], "NA"))
            muts.append("|".join([cands[i], cs[i], '']))
    elif int(nm_cutoff) > 0 :
        cdsseq=extractCDS(cdsseqfile, ",".join(cands))
        for i in range(cands_len) :
            cand_seq = cdsseq[cands[i]]['cdsseq']
            cand_seq_frame = cdsseq[cands[i]]['frame']
            dfs = at.findCodingDiff(cand_seq, cs[i], 
                                    int(cand_seq_frame)-1, '+',
                                    int(qfros[i]))
            newalleles.append(nameNewAllele(cands[i], dfs))
            muts.append("|".join([cands[i], cs[i], dfs]))


    # select allele with longest field length
    cands_field_len = []
    for i in range(cands_len) :
        cands_field_len.append(len(cands[i].split('*')[1]))
    max_cands_field_len = max(cands_field_len)

    sel_cands = []
    sel_muts = []
    sel_newalleles = []
    for i in range(cands_len) :
        if cands_field_len[i] == max_cands_field_len :
            sel_cands.append(cands[i])
            sel_newalleles.append(newalleles[i])
            sel_muts.append(muts[i])

    newallele_consensus = consensusCall(sel_newalleles)
    return({'consensus' : newallele_consensus,
            'NM' : nm_cutoff,
            'call' : ",".join(sel_cands),
            'mut' : ",".join(sel_muts)})

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


ctg_gene_index = {} # to make sure duplication genes to be annotated independently

with open(callfile, 'r') as f:
    while True:
        line = f.readline()
        if not line: break
        line =line.replace('\n', '')
        llst = line.split("\t")
        ctgid = llst[0]
        gene = llst[1]
        strand = llst[2]
        template = llst[3]
        allelecall = llst[4]
        template_mapping = llst[5]
        cds_mapping = 'NA'

        if ctgid not in ctg_gene_index :
            ctg_gene_index[ctgid] = 1
        else :
            ctg_gene_index[ctgid] += 1

        if allelecall != "NA" or template == "NA":
            ostr = [ctgid, gene, strand, template, allelecall,
                    "template_mapping " + '"' + template_mapping + '"',
                    "cds_mapping " + '"' + cds_mapping  + '"']
            ostr =  "\t".join(ostr)
            print(ostr)
            continue

        spaf = []
        for x in pafda :
            k =  ctgid + '_'+ gene + '_' + str(ctg_gene_index[ctgid])
            if x['tstr'] == k :
                spaf.append(x)

        if spaf :
            newcall = callNewAllele(spaf, cdsseqfile)
            #update gene name based on the allele calls
            ## for undetermined allele, use the best template match
            if '*' in newcall['consensus'] :
                gene = newcall['consensus'].split('*')[0]
            
            allelecall = newcall['call']
            cds_mapping = "consensus=" + newcall['consensus'] + ";"
            cds_mapping += "nm=" + newcall['NM'] + ";"
            cds_mapping += "mut=" + newcall['mut']  + ';"'


        ostr = [ctgid, gene, strand, template, allelecall,
                "template_mapping " + '"' + template_mapping + '"',
                "cds_mapping " + '"' + cds_mapping  + '"']
        ostr =  "\t".join(ostr)
        print(ostr)


