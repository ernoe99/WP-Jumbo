import datetime as dt
import pandas as pd

import numpy as np
import scipy as sp


def getval(sst, sdt):
    sstrange = rs[:, 1, 1:2]
    print(sstrange[0])
    sdtrange = rs[1, :, 2:3]
    print(sdtrange[0])
    if sst < sstrange[0] or sst > sstrange[-1] or sdt < sdtrange[0] or sdt > sdtrange[-1]:
        return 0

    for i in sstrange:
        if sdt > float(sstrange[i]):
            ir = [i - 1, i]

        for j in sdtrange:
            if sdt > float(sdtrange[j]):
                jr = [j - 1, j]

    print(ir)
    print(jr)

    res = [rs[jr[0], ir[0], :], rs[jr[1], ir[0], :], rs[jr[0], ir[1], :], rs[jr[1], ir[1], :]]
    print(res)

    return res


a = np.arange(6).reshape(2, 3)

print(a)

now = dt.time()

print(now)

time1 = dt.datetime.now()
print(time1)

time2 = dt.datetime(2021, 4, 16, 22)

print(time2)

delta = time2 - time1
print(delta)
print(type(delta))
print(delta.total_seconds())

syear = dt.timedelta(days=365.25)
print(syear.total_seconds())

print(delta.total_seconds() / syear.total_seconds())

cdTGH285 = pd.read_csv('..\\Datenfiles\\TGH285-noEcon.csv', header=9, sep=';')
indexx = pd.read_csv('..\\Datenfiles\\TGH285-noEcon.csv', header=9, sep=';', nrows=0)
print(indexx)

tgh285 = cdTGH285.to_numpy(copy=True)

print(tgh285[0, 0])

print(tgh285.ndim)
print(tgh285.shape)

rs = tgh285.reshape(17, 31, 39)

print(rs[:, 1, 1:2])
print(rs[1, :, 2:3])

rs_suc= rs[:, 1, 1:2]
rs_cond= rs[1, :, 2:3]
print(rs_suc.shape)
rsx = rs_suc.reshape(17)
print(rsx)

print(rs_cond.shape)

with open('TGH285_0P', 'wb') as wfile:
    np.save(wfile, np.rs)

print("File written")
with open('TGH285_0P', 'rb') as rfile:
    rss = np.load(rfile)

print(rss[1, 1])
print(shape.rss)



#l = getval(5.5, 20.1)
#print(l)

# print(cdTGH285.head(13))
# print(cdTGH285.tail(13))

# print(cdTGH285.ndim)
# print(cdTGH285.size)
# print(cdTGH285.shape)

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html#scipy.interpolate.RegularGridInterpolator
print("Hier")
print(rs[1, :, :])

tin = -14
tcond = 11

id1 = int((tin + 18) / 3)
print(id1)

jd1 = int((tcond+5)/2.5)

print(jd1)

slc = rs[id1:id1+2, jd1:jd1+2, :]

print("SLC")
print(slc)

res = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]  #40

xvec= np.empty(38)

print(xvec)

for i in range(1, 39, 1):

    # print(rs[id1, jd1, i])
    # print(rs[id1 + 1, jd1, i])
    w11 = rs[id1, jd1, i] + (rs[id1 + 1, jd1, i] - rs[id1, jd1, i]) / 3.0 * (tin - rs[id1,jd1,1])
#   print(tin - rs[id1, jd1+1, i])
#    print(w11)
    w12 = rs[id1, jd1 + 1, i] + (rs[id1 + 1, jd1 + 1, i] - rs[id1, jd1 + 1, i]) / 3.0 * (tin - rs[id1, jd1 + 1, 1])
#    print(w12)
    w23 = w11 + (w12 - w11)/2.5 * (tcond - rs[id1, jd1, 2])
#    print(w23)
    xvec[i-1] = w23


print(xvec)
print(xvec.shape)




