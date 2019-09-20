import sys
import operator

data = open(sys.argv[1]).readlines()
cutoff = int(sys.argv[2])
start_aux = end_aux = 0
contig_aux = ""
idn_aux = ""
strand_aux = ""
idn_list = {}
prod = "merged"
key = "repeat_region"
source = ""
if len(sys.argv) == 4:
    source = sys.argv[3]

l_total = len(data)
for i in xrange(1, l_total):
    contig_pre = data[i-1].split()[0]
    contig_now = data[i].split()[0]
    if contig_pre == contig_now:
        start_pre = int(data[i-1].split()[3])
        end_pre = int(data[i-1].split()[4])
        start_now = int(data[i].split()[3])
        end_now = int(data[i].split()[4])
        idn_now = data[i].split()[8].split(";")[0].split("=")[1]
        idn_pre = data[i-1].split()[8].split(";")[0].split("=")[1]
        strand = data[i].split()[6]

        if (start_now >= end_pre and start_now - end_pre <= cutoff) or (start_now >= start_pre and start_now <= end_pre):
            if start_aux == 0:
                start_aux = start_pre
                end_aux = end_now
                idn_aux = idn_pre
                contig_aux = contig_pre
                strand_aux = strand
                idn_list[idn_now] = end_now - start_now
                idn_list[idn_aux] = end_pre - start_pre
            else:
                end_aux = end_now
                idn_list[idn_aux] = end_pre - start_pre
                idn_list[idn_now] = end_now - start_now
        else:
            if start_aux == 0:
                print data[i-1].strip()
            else:
                length = str(end_aux - start_aux)
                long_frag = max(idn_list.iteritems(), key=operator.itemgetter(1))
                print "\t".join(map(str, [contig_pre, source, key, start_aux, end_aux,".",strand_aux,".","ID="+long_frag[0]+";"]))
                start_aux = end_aux = 0
                contig_aux = ""
                strand_aux = ""
                idn_list = {}
    else:
        if start_aux == 0:
            print  data[i-1].strip()
        else:
            length = str(end_aux - start_aux)
            long_frag = max(idn_list.iteritems(), key=operator.itemgetter(1))
            print "\t".join(map(str, [contig_pre, source, key, start_aux, end_aux,".",strand_aux,".","id="+long_frag[0]+";"]))
            start_aux = end_aux = 0
            contig_aux = ""
            strand_aux = ""
            idn_list = {}
    if l_total == i+1:
        if start_aux == 0:
#            print  data[i-1].strip()+";"
            print  data[i].strip()+";"
        else:
            length = str(end_aux - start_aux)
            long_frag = max(idn_list.iteritems(), key=operator.itemgetter(1))
            print "\t".join(map(str, [contig_pre, source, key, start_aux, end_aux,".",strand_aux,".","id="+long_frag[0]+";"]))
            start_aux = end_aux = 0
            contig_aux = ""
            strand_aux = ""
            idn_list = {}

