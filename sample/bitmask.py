
flag = input("Enter flag ")

m = 65

for i in xrange(m):
    #print i, 2**i
    div, mod = divmod(flag,2)
    flag = div
    if mod==1:
        print "Bit {0:2d} TRUE".format(i)
