import sys
import math

def get_tot_comp(f, seq_len):
    reader = open(f)
    lines = reader.readlines()
    tot_comp = 0
    for line in lines:
        tot_comp += int(float(line.split(" ")[1].strip()))

    dis = seq_len - tot_comp
    print ede(dis, seq_len)


def ede ( invdist, ngene ):
    ede_a = 0.595639
    ede_c = 0.457748
    ll=0.0
    tt=0.0
    kk=0.0
    pp=0.0
    dval=0.0
    newvalue=0

    kk = invdist / ( ngene + 0.0 )

    if kk >= 0.999999999999 :
        kk = 0.999999999999
    if kk <= 1 - ede_c :
        return invdist

    ll = ede_c * kk - ede_a
    tt = 4 * ede_a * ( 1 - kk ) * kk + ll * ll
    tt = ll + math.sqrt(tt)
    pp = tt / ( 2 * ( 1 - kk ) )
    pp *= ngene

    dval = pp
    newvalue = int(math.ceil(dval))
    return newvalue


get_tot_comp(sys.argv[1], int(sys.argv[2]))
