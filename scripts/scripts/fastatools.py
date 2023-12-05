import sys,re,gzip

def matchKeyward(line, keys) :
    retk = ''
    ks = keys.split(',')
    k = line.split(' ')[0][1:]
    if k in ks:
        retk = k
    return(retk)

def loadFastagz(fastafile, keys) :
    ret = {}
    buf = ''
    buftag = 0 # 0 is empty and 1 is full
    key = ''
    ks = keys.split(',')
    l0 = len(ks)

    with gzip.open(fastafile, 'rt') as f:
        while True:
            line = f.readline()
            if not line or len(ks) == 0: break

            if line[0] == '>' :
                # output the seq or not
                if key and key in ks:
                    ret[key] = buf
                    ks.remove(key)
                # clean the buffer
                buf = ''
                # find the key
                key = matchKeyward(line, keys)
            else:
                if key:
                    buf += line.replace('\n', '')
    if buf and key:
        ret[key] = buf
    return(ret)


def loadFasta(fastafile, keys) :
    ret = {}
    buf = ''
    buftag = 0 # 0 is empty and 1 is full
    key = ''

    with open(fastafile, 'rt') as f:
        while True:
            line = f.readline()
            if not line or len(ks) == 0: break

            if line[0] == '>' :
                # output the seq or not
                if key :
                    ret[key] = buf
                    ks.remove(key)
                # clean the buffer
                buf = ''
                # find the key
                key = matchKeyward(line, keys)
            else:
                buf += line.replace('\n', '')
    if buf and key:
        ret[key] = buf
    return(ret)
