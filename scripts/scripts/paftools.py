import re,sys
import pprint as pp

def parsingPaf(line) :
    # input one str of a row of minimap2
    # output a dict dormat
    ret = {}
    l = line.replace('\n', '')
    lst = l.split("\t")
    #First 12 field
    ret['qstr'] = lst[0]
    ret['qlen'] = lst[1]
    ret['qfro'] = lst[2]
    ret['qto'] = lst[3]
    ret['strand'] = lst[4]
    ret['tstr'] = lst[5]
    ret['tlen'] = lst[6]
    ret['tfro'] = lst[7]
    ret['tto'] = lst[8]
    ret['match'] = lst[9]
    ret['bpnum'] = lst[10]
    ret['mpqual'] = lst[11]
    # attribute columns
    if len(lst) > 12 :
        for s in lst[12:len(lst)]:
            key,tag,value = s.split(":", 2)
#            print(s,key,tag,value)
            ret[key] = value
    return(ret)


def adjNM(csstr) :
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)
    anm = 0
    L = len(cslst)
    for i in range(L-1) :
        if cslst[i] == '*' : anm += 1
        elif cslst[i] in ['+', '-'] :
            if len(cslst[i+1]) > 5 : anm += 2
            else : anm += 1

    return(str(anm))


def adjNM_bk(csstr) :
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)
    anm = 0
    for x in cslst :
        c0 = x == '*'
        c1 = x == '+'
        c2 = x == '-'
        if c0 or c1 or c2 :
            anm += 1
    return(str(anm))
