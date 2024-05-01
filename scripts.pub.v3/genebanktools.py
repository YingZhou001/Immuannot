import sys,re


def extractExon(buff) :
    regex = re.compile(
        r'exon\s*(\d*\.\.\d*).*?number=(\d*)\n', re.DOTALL
    )
    ret = regex.findall(buff)
    if ret :
        out = []
        for i in ret :
            x, y = i
            out.append(y + ':' + x)
        ret = ",".join(out)
    return(ret)


def extractCDS(buff) :
    regex = re.compile(
        r'CDS\s*join\((.*?)\)\n.*?gene="(.*?)".*?codon_start=(\d)', re.DOTALL
    )
    ret = regex.search(buff)
    if ret :
        s = ret.group(1)
        gene = ret.group(2)
        codon_start = ret.group(3)
        s = re.sub('[\n|\s]', '', s)
        return([gene, codon_start, s])

def extractSeq(buff) :
    regex = re.compile(
        r'ORIGIN(.*?)//', re.DOTALL
    )
    ret = regex.search(buff)
    if ret :
        s = ret.group(1)
        s = re.sub('[\n|\s|\d]', '', s).upper()
        return(s)

