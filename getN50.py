import sys

lengths = [int(i.strip()) for i in open(sys.argv[2])]
lengths.sort(reverse=True)
total = sum(lengths)*1.0
count = 0
n = 0
nx = int(sys.argv[1])*1.0/100
for i in lengths:
    if count < total*nx:
        n+=1
        count += i
    else:
        print "N"+str(sys.argv[1])+": "+ str(i) + "\tn: "+str(n) + "\ttotal: "+ str(int(total))
        break
