import fitsio
import numpy as np
import matplotlib.pyplot as plt

file = '/Volumes/Transcend/Isotropy/spectra/realspec/spec-3615-55445-0008.fits'
#file = '/Volumes/Transcend/Isotropy/spectra/realspec/spec-5964-56098-0162.fits'
data, header = fitsio.read(file,ext=0,header=True)
coadded = fitsio.read(file,ext=1,columns=['LOGLAM','FLUX'])

print header

exp1, head1 = fitsio.read(file,ext=4,columns=['LOGLAM','FLUX'],header=True)
exp2, head2 = fitsio.read(file,ext=5,columns=['LOGLAM','FLUX'],header=True)
exp3, head3 = fitsio.read(file,ext=6,columns=['LOGLAM','FLUX'],header=True)
exp4, head4 = fitsio.read(file,ext=7,columns=['LOGLAM','FLUX'],header=True)
exp5, head5 = fitsio.read(file,ext=8,columns=['LOGLAM','FLUX'],header=True)
exp6, head6 = fitsio.read(file,ext=9,columns=['LOGLAM','FLUX'],header=True)
exp7, head7 = fitsio.read(file,ext=10,columns=['LOGLAM','FLUX'],header=True)
exp8, head8 = fitsio.read(file,ext=11,columns=['LOGLAM','FLUX'],header=True)
exp9, head9 = fitsio.read(file,ext=12,columns=['LOGLAM','FLUX'],header=True)
exp10, head10 = fitsio.read(file,ext=13,columns=['LOGLAM','FLUX'],header=True)
exp11, head11 = fitsio.read(file,ext=14,columns=['LOGLAM','FLUX'],header=True)
exp12, head12 = fitsio.read(file,ext=15,columns=['LOGLAM','FLUX'],header=True)

s0 = np.size(coadded)
s1 = np.size(exp1)+np.size(exp7)
s2 = np.size(exp2)+np.size(exp8)
s3 = np.size(exp3)+np.size(exp9)
#s4 = np.size(exp4)+np.size(exp10)
#s5 = np.size(exp5)+np.size(exp11)
#s6 = np.size(exp6)+np.size(exp12)


wav0 = []
wav1 = []
wav2 = []
wav3 = []

flx0 = []
flx1 = []
flx2 = []
flx3 = []

for i in xrange(s0):
    wav0.append(10**coadded[i][1])
    flx0.append(coadded[i][0])
j=0
for i in xrange(s1):
    if i<np.size(exp1):
        wav1.append(10**exp1[j][1])
        flx1.append(exp1[j][0])
        j = j + 1
    else:
        t = i - j
        wav1.append(10**exp7[t][1])
        flx1.append(exp7[t][0])
j=0
for i in xrange(s2):
    if i<np.size(exp2):
        wav2.append(10**exp2[i][1])
        flx2.append(exp2[i][0])
        j = j + 1
    else:
        t = i - j
        wav2.append(10**exp8[t][1])
        flx2.append(exp8[t][0])
j=0
for i in xrange(s3):
    if i<np.size(exp3):
        wav3.append(10**exp3[i][1])
        flx3.append(exp3[i][0])
        j = j + 1
    else:
        t = i - j
        wav3.append(10**exp9[t][1])
        flx3.append(exp9[t][0])
    
pixperbin = 4
print s1, s2, s3
n0 = int(float(s0)/pixperbin)
n1 = int(float(s1)/pixperbin)
n2 = int(float(s2)/pixperbin)
n3 = int(float(s3)/pixperbin)

newwav0 = []
newwav1 = []
newwav2 = []
newwav3 = []

newflx0 = []
newflx1 = []
newflx2 = []
newflx3 = []

for i in xrange(0,s0,pixperbin):
    sumwavbin = 0
    sumflxbin = 0
    for j in xrange(i,i+pixperbin):
        if j==s0:
            break
        else:
            sumwavbin = sumwavbin + wav0[j]
            sumflxbin = sumflxbin + flx0[j]
    newwav0.append(float(sumwavbin)/(pixperbin))
    newflx0.append(float(sumflxbin)/pixperbin)
    
for i in xrange(0,s1,pixperbin):
    sumwavbin = 0
    sumflxbin = 0
    for j in xrange(i,i+pixperbin):
        if j==s1:
            break
        else:
            sumwavbin = sumwavbin + wav1[j]
            sumflxbin = sumflxbin + flx1[j]
    newwav1.append(float(sumwavbin)/pixperbin)
    newflx1.append(float(sumflxbin)/pixperbin)

for i in xrange(0,s2,pixperbin):
    sumwavbin = 0
    sumflxbin = 0
    for j in xrange(i,i+pixperbin):
        if j==s2:
            break
        else:
            sumwavbin = sumwavbin + wav2[j]
            sumflxbin = sumflxbin + flx2[j]
    newwav2.append(float(sumwavbin)/pixperbin)
    newflx2.append(float(sumflxbin)/pixperbin)
for i in xrange(0,s3,pixperbin):
    sumwavbin = 0
    sumflxbin = 0
    for j in xrange(i,i+pixperbin):
        if j==s3:
            print j,s3, 'break'
            break
        else:
            sumwavbin = sumwavbin + wav3[j]
            sumflxbin = sumflxbin + flx3[j]
    newwav3.append(float(sumwavbin)/pixperbin)
    newflx3.append(float(sumflxbin)/pixperbin)

tail = 200
del newwav0[-1]
del newwav1[-1]
del newwav2[-1]
del newwav3[-1]

del newflx0[-1]
del newflx1[-1]
del newflx2[-1]
del newflx3[-1]

plt.figure()
plt.plot(newwav0, newflx0, newwav1,newflx1, newwav2,newflx2, newwav3,newflx3)
plt.xlim(3579,5779)
plt.ylim(-2,50)
plt.show()
