import sys,re
import pprint as pp

def translate(codonseq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    prot = ''
    #s = ''
    l = len(codonseq)
    if l % 3 == 0 :
        n = l //3
        for i in range(n) :
            s = codonseq[3*i:3*i+3].upper()
            if s in table.keys() :
                prot += table[s]
            else :
                sys.stderr.write("##"+s + '\n')
    else :
        prot = codonseq
    return(prot)

def calOverlap(a, b) :
    x0, y0 = a
    x1, y1 = b
    if y0 <= x1 or x0 >= y1 :
        return([])
    else :
        return([max(x0, x1), min(y0, y1)])


def reversedComplimentStr(s) :
    map = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G',
          'a' : 't', 't' : 'a', 'g' : 'c', 'c' : 'g', 
           'N' : 'N'}
    ret = ''
    for i in list(reversed(s)) :
        ret += map[i]
    return(ret)

def recoverTargetSeqFromQuery(qseq, csstr, qrg, strand) :
    # recover the targeted seq based on the alignment
    # qrg is 1-based coordinates
    # always return the strand the same as qseq
    if strand == '-' :
        seq = reversedComplimentStr(qseq)
        l = len(qseq)
        a,b = qrg.split('..')
        fro = l - int(b)
        to = l - int(a) + 1
    elif strand == '+' :
        seq = qseq
        fro, to = qrg.split('..')
        fro = int(fro) - 1
        to = int(to)
    #sys.stderr.write('*****'+seq[1:400] + '\n')
    tseq = ''
    #print(qrg)
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)[1:]
    n = len(cslst) // 2
    cspairs = []
    for i in range(n) :
        cspairs.append([cslst[i*2], cslst[i*2+1]])
    curp = 0
    for i in range(n) :
        tag = cspairs[i][0]
        value = cspairs[i][1]
        #if fro-10 < curp and curp < to +10 :
            #print(curp, tag, value)
        if tag == ':' :
            v = int(value)
            newp = curp + v
            x = calOverlap([curp, newp], [fro, to])
            if x :
                #print(x, '***', [curp, newp], [fro, to])
                #print('R',tseq, seq[x[0]:x[1]].upper())
                #print('Q', seq[fro:to])
                tseq += seq[x[0]:x[1]]
        elif tag == '~' :
            v = int(re.sub('[a-z]', '', value))
            if fro < curp and curp < to :
                tseq += 'N'*v
        elif tag == '*' :
            newp = curp + 1
            if fro <= curp and curp < to :
                #print('R', tseq, value[0].upper())
                #print('Q', seq[fro:to])
                tseq += value[0].upper()
        elif tag == '+' :
            newp = curp + len(value)
        elif tag == '-' :
            newp = curp + 0
            if fro <= curp and curp < to :
                #print('R', tseq, value.upper())
                #print('Q', seq[fro:to])
                tseq += value.upper()
        #sys.stderr.write('>'+ tag + ' ' + value + ' '
        #                 + str([curp, newp]) + ' +' + str(newp-curp) 
        #                 + str([fro, to]) + '\n')
        curp = newp
    #print('R', tseq)
    #print('Q', seq[fro:to])
    
    if newp > len(seq) :
        #sys.stderr.write(qseq + str(len(qseq)) + '\n')
        #sys.stderr.write(seq + str(len(seq)) + '\n')
        #sys.stderr.write(csstr + '\n')
        #sys.stderr.write(qrg + '\n')
        #sys.stderr.write(str([newp, len(qseq)]) + '\n')
        sys.stderr.write("Warning : CS string length longer than the seq length\n")

    if strand == '-' :
        tseq = reversedComplimentStr(tseq)

    #if len(tseq) != abs(to - fro) :
    #    sys.stderr.write('#'+ strand + seq + ' ' + str(len(seq)) + '\n')
    #    sys.stderr.write('#'+csstr + ' ' + str(curp) + '\n')
    #    sys.stderr.write('#'+qrg + ' ' + str(fro) + ' ' + str(to) + '\n')
    #    sys.stderr.write('#'+tseq + '\n')

    return(tseq)


def calCodonPos0(pos, frame, strand, seqRange) :
    # all are all in 0-based coordinate
    # return 0, 1, and 2 represnet the first, second, and the last position in
    # a codon
    fro,to = seqRange
    if strand == '+' :
        phase = (pos - int(frame) - fro) % 3
        pfro = pos - phase
        pto = pos + (3 - phase)
    elif strand == '-' :
        phase = (to - 1 - pos - int(frame)) %3
        pfro = pos + phase + 1
        pto = pos - (2 - phase)
    return([phase, pfro, pto])

def findCodingDiff(qcdsseq, csstr, frame, strand, qstrand = '+') :
    # frame in zero based coordinate
    # always assume that the qseq is positive strand
    # strand is represent the difference between
    ret = []
    seql = len(qcdsseq)
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)[1:]
    n = len(cslst) // 2
    cspairs = []
    for i in range(n) :
        cspairs.append([cslst[i*2], cslst[i*2+1]])
    if strand == '-' :
        # reverse the cs str to make the position match qseq
        cspairs = list(reversed(cspairs))

    #locate target interval, 0 based coordinate
    out = [] 
    cur = 0
    for i in range(n) :
        tag = cspairs[i][0]
        value = cspairs[i][1]
        if tag == ':' :
            # equal
            l = int(value)
            cur += l
        elif tag == '*' : 
            #single base substitution
            l = 1
            p,b0,c0 = calCodonPos0(cur, frame, qstrand, [0,seql])
            p,b1,c1 = calCodonPos0(cur+l-1, frame, qstrand, [0,seql])
            out.append([min(b0, b1), max(c0, c1)])
            cur += l
        elif tag == '+' : 
            #insertion
            l = len(value)
            p,b0,c0 = calCodonPos0(cur, frame, qstrand, [0,seql])
            p,b1,c1 = calCodonPos0(cur+l-1, frame, qstrand, [0,seql])
            out.append([min(b0, b1), max(c0, c1)])
            cur += l
        elif tag == '-' : 
            # deletion
            l = 0
            p,b0,c0 = calCodonPos0(cur, frame, qstrand, [0,seql])
            p,b1,c1 = calCodonPos0(cur+l-1, frame, qstrand, [0,seql])
            out.append([min(b0, b1), max(c0, c1)])
            cur += l

    #pp.pprint(out)
    newout = []
    for x in out:
        if not x in newout :
            newout.append(x)
    out = newout

    # check CDS length change not times of 3
    # check early termination
    qslen = 0
    rslen = 0

    nonsense = 0

    for x in out :
        a,b = x
        qs = qcdsseq[a:b]
        qrg = str(a+1)+ '..' + str(b)
        rs = recoverTargetSeqFromQuery(qcdsseq, csstr, qrg, strand)
        if strand == '-' :
            qs = reversedComplimentStr(qs)
            rs = reversedComplimentStr(rs)
        tqs = translate(qs)
        trs = translate(rs)
        qslen += len(qs)
        rslen += len(rs)
        if tqs != trs and '_' in trs+tqs :
            nonsense = 1
        ret.append(translate(qs) + "(" + qs + ")" + "<" +
                   translate(rs) + "(" + rs + ")")
    if (qslen - rslen) %3 != 0 :
        nonsense = 1
    return(":".join(ret))


def cgMapQry2Ctg(x, cg, strand):
    # input x: 1 based coordinate in the query seq
    # return value: 1 based position in the mapped seq
    # a seq of length L is represented as [1, L]
    y = int(x)
    # parsing cigar string
    cglist = re.split(r'([A-Z])', cg)
    out = 0
    rgs = range(len(cglist) // 2)
    if strand == '-' : rgs = reversed(rgs)
    for i in rgs:
        k = i*2
        c = int(cglist[k])
        tag = cglist[k + 1]
        if tag == 'M' :
            out += min(y, c)
        elif tag == 'I' :
            out += 0
        elif tag == 'N' or tag == 'D':
            out += c
            y += c
        else :
            sys.exit("Wrong Cigar string:" + cg + "\n")
        if y <= c :
            y -= c
            break
        else :
            y -= c
    if y > 0 :
        out += y
        #sys.stderr.write("Warning: mapping position is out of cigar ranges\n")
        #sys.stderr.write(str(x) +' '+cg +'\n')
        #sys.exit()
    #    sys.exit("Error: mapping position longer than the query seq length\n")
    return(out)

def intervalMapQry2Ctg(qcds, cg, strand, qfro, qto, cfro, cto) :
    # q is a string of interval in the format as a..b,c..d,e..f
    # all input should be 1 based coordinate
    qarr = []
    warnings = []
    for x in qcds.split(",") :
        a, b = x.split("..")
        a = int(a)
        b = int(b)
        if min(a, b) < int(qfro) or max(a, b) > int(qto) :
            #check the cds region is fully covered or not
            wcode = "incomplete_mapping"
            if wcode not in warnings :
                warnings.append(wcode)
            if strand == '+' :
                qarr.append(str(-1) + ".." + str(-2))
            else :
                qarr.insert(0, str(-1) + ".." + str(-2))
            continue
        if strand == '+' :
            c = cgMapQry2Ctg(a, cg, strand) - int(qfro) + int(cfro)
            d = cgMapQry2Ctg(b, cg, strand) - int(qfro) + int(cfro)
            qarr.append(str(c) + ".." + str(d))
        elif strand == '-' :
            c = int(cto) - (cgMapQry2Ctg(a, cg, strand) - int(qfro))
            d = int(cto) - (cgMapQry2Ctg(b, cg, strand) - int(qfro))
            qarr.insert(0, str(d) + ".." + str(c))
        if abs(b - a) != abs(d -c) : 
            # check frameshift mutation
            wcode = "unequal_length"
            if wcode not in warnings :
                warnings.append(wcode)
    ret = [",".join(warnings), ",".join(qarr)]
    return(ret)
