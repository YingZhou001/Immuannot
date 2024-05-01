# extract useful information from the dat file

import re
import sys

def readMETA(str) : 
    # return the seq id, gene seq length, allele type,
    # and description
    ret = {}
    regex = re.compile(r'ID\s{3}(.*); SV \d+; standard; DNA; HUM; (\d+) BP.')
    ret1 = regex.search(str)
    if ret1 :
        ret['ID'] = ret1.group(1)
        ret['genelength'] = ret1.group(2)
    else:
        ret['ID'] = "NA"
        ret['genelength'] = "NA"

    regex = re.compile(r'DE\s{3}(.*), (.*)\n')
    ret2 = regex.search(str)
    if ret2 :
        ret['allele'] = ret2.group(1)
        ret['description'] = ret2.group(2)
    else :
        ret['allele'] = "NA"
        ret['description'] = "NA"

    regex1 = re.compile(r'FT\s*UTR')
    regex2 = re.compile(r'FT\s*intron')
    ret1 = regex1.search(str)
    ret2 = regex2.search(str)
    if ret1 or ret2 :
        ret['seqtype'] = "gene"
    else :
        ret['seqtype'] = "CDS"

    regex = re.compile(r'FT\s*/translation="')
    ret2 = regex.search(str)
    if ret2 :
        ret['type'] = "gene"
    else :
        ret['type'] = "unkown"

    regex = re.compile(r'CC\s*(.*) Release Version (.*)\n')
    ret2 = regex.search(str)
    if ret2 :
        ret['version'] = ret2.group(1)+'-V'+ret2.group(2)
    else :
        ret['version'] = "NA"
    return ret

def readGENE(str) :
    # catch CDS intervals
    regex1 = re.compile(r'CDS\s*.*codon_start', re.DOTALL)
    ret1 = regex1.search(str)
    if ret1:
        regex2 = re.compile(r'(\d*\.\.\d*)')
        ret2 = regex2.findall(ret1.group())
    else :
        ret2 = None

	# catch codon starting position
    regex3 = re.compile(r'/codon_start=(\d*)\n')
    ret3 = regex3.search(str)
    if ret3 :
        ret3 = ret3.group(1)

    # catch UTR intervals
    regex4 = re.compile(r'UTR\s*(\d*\.\.\d*)\n')
    ret4 = regex4.findall(str)

    # catch exons
    regex5 = re.compile(
        r'exon\s*(\d*\.\.\d*)\nFT\s*/number="(\d*)"\n(FT\s*/pseudo)?')
    #regex5 = re.compile(r'exon\s*(\d*\.\.\d*)\nFT\s*/number="(\d*)"\n')
    ret5 = regex5.findall(str)
    #print(ret5, file=sys.stderr)

    #catch intron
    regex6 = re.compile(r'intron\s*(\d*\.\.\d*)\nFT\s*/number="(\d*)"\n')
    ret6 = regex6.findall(str)

    # catch gene name
    regex7 = re.compile(r'/gene="(.*)"\n')
    ret7 = regex7.search(str)
    if ret7 :
        ret7 = ret7.group(1)
    else :
        ret7 = "NA"

    ret = {}
    ret['CDS'] = ret2
    ret['codon_start'] = ret3
    ret['UTR'] = ret4
    ret['exon'] = ret5
    ret['intron'] = ret6
    ret['gene'] = ret7

    return ret

def readGENESeq(str) :
    ret = {}
    regex = re.compile(r'\nSQ.*\n((\s*.*\n)*)//')
    ret1 = regex.search(str)
    if ret1 :
        s = re.sub('\s|\d', '', ret1.group(1)).upper()
        ret['geneSeq'] = s
    else :
        ret['geneSeq'] = "NA"
    return ret

def readPROTEIN(str) :
    ret = {}
    regex = re.compile(r'translation=(.*\n(FT\s{19}.*\n)*)')
    ret1 = regex.search(str)
    if ret1 :
        s = re.sub('FT|"|\n|\s', '', ret1.group(1))
        ret['protein'] = s
    else :
        ret['protein'] = "NA"
    #print(s)
    return ret


def writeSTDOUTprotein(df) :
    s = df['protein']
    headline = ">" + df['allele'] + " " + df['ID'] + " " \
            + df['description'].replace(' ', '_') + " " \
            + str(len(s)) + '\n'
    sys.stdout.write(headline)
    #print('_' + s + '_')
    for i in range(len(s) // 60 + 1) :
        ostr = s[60 * i : 60 * (i + 1)] + '\n'
        sys.stdout.write(ostr)


def writeSTDOUTgeneRecordsFirstline() :
    firstline = 'gene' + "\t" + 'ID' + "\t" + 'allele' \
            + "\t" + 'description' + "\t" + 'codon_start' \
            + "\t" + "element" + "\t" + "index" \
            + "\t" + "from\tto" + '\n'
    sys.stdout.write(firstline)


def writeSTDOUTgeneRecords(df) :
    if df['gene'] and df['ID'] and df['allele'] \
       and df['description'] and df['codon_start'] :
        pref = df['gene'] + "\t" + df['ID'] + "\t" \
                + df['allele'] + "\t" + df['description'] \
                + "\t" + df['codon_start']
    else :
        return None
    # output gene
    ostr = pref + "\t" + "gene" + "\t" + "0" \
            + "\t" + "1\t" + df['genelength'] + '\n'
    sys.stdout.write(ostr)
    # output CDS
    cds = df['CDS']
    if cds:
        for i in range(len(cds)):
            ostr = pref + "\t" + "CDS" + "\t" + str(i+1) \
                    + "\t" + cds[i].replace('..', "\t") \
                    + '\n'
            sys.stdout.write(ostr)
    # output UTR
    utr = df['UTR']
    if utr:
        for i in range(len(utr)) :
            ostr = pref + "\t" \
                    + "UTR" + "\t" + str(i+1) + "\t" \
                    + utr[i].replace('..', "\t") + '\n'
            sys.stdout.write(ostr)

    # output exon
    exon = df['exon']
    if exon:
        for i in range(len(exon)) :
            ostr = pref + "\t"\
                    + "exon" + "\t" \
                    + list(exon[i])[1] + "\t" \
                    + list(exon[i])[0].replace('..', "\t") \
                    + '\n'
            sys.stdout.write(ostr)

def writeSTDOUTgeneElement(df) :
    if df['gene'] and df['ID'] and df['allele'] \
       and df['description'] and df['codon_start'] :
        pref = "gene=" + df['gene'] + "\t" +\
                "version=" + df['version'] + "\t" +\
                "ID=" + df['ID'] + "\t" +\
                "allele=" + df['allele'] + "\t" +\
                "description=" + df['description'].replace(' ', '_') + "\t" +\
                "type=" + df['type']
    else :
        return None

    ostr = ''
    # output gene
    ostr = pref + "\t" + "geneRange=" + "1.." + df['genelength']

    # output CDS
    cds = df['CDS']
    if cds:
        cdsstr = ",".join(cds)
    else :
        cdsstr = "NA"
    ## search for start codon and end codon
    cdsArr = []
    for x in cds :
        a, b = x.split("..")
        for y in range(int(a), int(b)+1) : cdsArr.append(y)
    cdsArr.sort()
    n = len(cdsArr)
    shift = int(df['codon_start'])
    nCodon = (n - shift + 1)// 3
    firstCodon = cdsArr[shift - 1 ]
    lastCodon = cdsArr[shift -1 + (nCodon - 1) * 3]

    #print([len(cdsArr), firstCodon, lastCodon])
    seq = df["geneSeq"]
    startCodon = seq[firstCodon-1:firstCodon+2]
    stopCodon = seq[lastCodon-1:lastCodon+2]

    # for the case that the stop codon is not included in the CDS
    if not stopCodon in ["TAA", "TAG", "TGA"] and lastCodon+5 <= len(seq) :
        tmp = seq[lastCodon+2:lastCodon+5]
        if len(tmp) == 3 and tmp in ["TAA", "TAG", "TGA"] :
            stopCodon = tmp
            lastCodon = lastCodon + 3

    # output UTR and exon
    utr = df['UTR']
    utrArr = []
    if utr:
        utrstr = ",".join(utr)
        for u in utr :
            a, b = u.split("..")
            utrArr.append([int(a), int(b)])
    else :
        utrstr = "NA"

    exon = df['exon']
    nexon = 0
    if exon:
        nexon = len(exon)
        exonstr = []
        for i in range(nexon) :
            lst = list(exon[i])
            if lst[2] != '' : lst[2] = "pseudo"
            if utr:
                # merge UTR into exon
                c, d = lst[0].split("..")
                c = int(c)
                d = int(d)
                rmv = []
                for x in utrArr :
                    a, b = x
                    #print(a,b, c, d)
                    if c == b + 1 : 
                        c = a
                        rmv.append(x)
                    elif c > b+1 :
                        break
                    if d == a -1 : 
                        d = b
                        rmv.append(x)
                    elif d < a-1:
                        break
                    #print(a,b, c, d)
                for x in rmv : utrArr.remove(x)
                lst[0] = str(c) + ".." + str(d)
            exonstr.append(lst[2] + lst[1] + ":" + lst[0])
        if len(utrArr) > 0 :
            # turn unmerged UTR into new exon
            for c, d in utrArr:
                exonstr.append("utr:" + str(c) + ".." + str(d))
        exonstr = ",".join(exonstr)


    intron = df['intron']
    if intron:
        intronstr = []
        for i in intron :
            a, b = i[0].split("..")
            a = int(a)
            b = int(b)
            l = seq[a - 1 : a + 1]
            r = seq[b - 2 : b ]
            x = l + "-" + r
            if not x in intronstr :
                intronstr.append(x)
        intronstr.sort()
        intronstr = ",".join(intronstr)
    else :
        intronstr = "NA"


    goodtemplate = False
    c1 = len(utrstr.split(",")) == 2 # UTR to UTR
    #c2 = startCodon == "ATG" and stopCodon in ["TAA", "TAG", "TGA"]
    c2 = True #startCodon == "ATG" and stopCodon in ["TAA", "TAG", "TGA"]
    if c1 and c2 :
        goodtemplate = True

    exonArr = []
    for x in exon :
        a, b = x[0].split("..")
        for y in range(int(a), int(b)+1) : exonArr.append(y)
    exonArr.sort()

    if len(exonArr) != len(cdsArr) and "pseudo" not in exonstr and False:
        sys.stderr.write(df['ID'] + ": CDS not euqals to Exon\n")
        sys.stderr.write("CDS = " + cdsstr + "\n")
        sys.stderr.write("Exon = " + exonstr + "\n")

    ostr = ostr + "\t" + "goodTemplate=" + str(goodtemplate)
    ostr = ostr + "\t" + "codonShift=" + df['codon_start']
    ostr = ostr + "\t" + "startCodon=" + startCodon + ":" 
    ostr = ostr + str(firstCodon) + ".." + str(firstCodon+2)
    ostr = ostr + "\t" + "stopCodon=" + stopCodon + ":"
    ostr = ostr + str(lastCodon) + ".." + str(lastCodon+2)
    ostr = ostr + "\t" + "CDS=" + cdsstr
    ostr = ostr + "\t" + "UTR=" + utrstr
    ostr = ostr + "\t" + "nexon=" + str(nexon)
    ostr = ostr + "\t" + "exon=" + exonstr
    ostr = ostr + "\t" + "intron=" + intronstr
    ostr = ostr + "\n"
    sys.stdout.write(ostr)


def writeSTDOUTalleleName(df) :
    if  df['gene'] != "NA" :
        headline = df['gene'] + " " + df['allele'] + " " + df['ID'] + " " \
                + df['description'].replace(' ', '_') + " " \
                + str(df['genelength']) + 'bp\n'
        sys.stdout.write(headline)

def writeSTDOUTgeneseq(df) :
    if df['seqtype'] == "gene" :
        s = df['geneSeq']
        headline = ">" + df['allele'] + " " + df['ID'] + " " \
                + df['description'].replace(' ', '_') + " " \
                + str(df['genelength']) + 'bp\n'
        sys.stdout.write(headline)
        for i in range(len(s) // 80 + 1) :
            ostr = s[80 * i : 80 * (i + 1)] + '\n'
            sys.stdout.write(ostr)

def writeSTDOUTcds(df, verbose) :
    seq = df['geneSeq']
    if df['CDS'] == None :
        if verbose: 
            sys.stderr.write(' <msg> CDS is missing for ' + df['ID'] + ", skip\n")
        return None

    cds = ''
    regex = re.compile(r'(\d*)..(\d*)')
    for i in df['CDS'] :
        rg = regex.search(i)
        fro = int(rg.group(1)) - 1
        to = int(rg.group(2))
        cds = cds + seq[fro:to]

    headline = ">" + df['allele'] + " " + df['ID'] + " " \
            + "frame=" + df['codon_start'] + " " + str(len(cds)) + 'bp\n'
    sys.stdout.write(headline)
    for i in range(len(cds) // 80 + 1) :
        ostr = cds[80 * i : 80 * (i + 1)] + '\n'
        sys.stdout.write(ostr)
