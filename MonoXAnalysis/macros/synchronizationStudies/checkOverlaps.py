import sys

fi1N,fi2N = sys.argv[1],sys.argv[2]
outputName = sys.argv[3]

fi1 = open(fi1N)
fi2 = open(fi2N)

events1 = [line.split() for line in fi1.readlines()]
events2 = [line.split() for line in fi2.readlines()]

print "size before filter for list1 ",len(events1)," list2 ",len(events2);
events1 = filter(lambda x: len(x), events1)
events2 = filter(lambda x: len(x), events2)
print "size after filter for list1 ",len(events1)," list2 ",len(events2);

events_1_in_2 = filter(lambda x: x in events2,events1)
totalevents = events1+events2
totalevents = filter(lambda x: x not in events_1_in_2,totalevents)
totalevents += events_1_in_2

events_1_not_in_2 = filter(lambda x: x not in events2,events1)

numoverlaps = len(events_1_in_2)
numbernotoverlap = len(events_1_not_in_2)
numtotal    = len(totalevents)

print "produce output file "
output = open(outputName,"w");
for ev in events_1_not_in_2:
    line = ""
    for item in ev:
        line += item +" ";
    output.write(line+"\n")
output.close()
print "Total events = ", numtotal
print "Found %d overlapping events between %s and %s (%g %%) wrt total"%(numoverlaps, fi1N, fi2N,100*float(numoverlaps)/numtotal)
print "Found %d overlapping events between %s and %s (%g %%) wrt file1"%(numoverlaps, fi1N, fi2N,100*float(numoverlaps)/len(events1))
print "Found %d overlapping events between %s and %s (%g %%) wrt file2"%(numoverlaps, fi1N, fi2N,100*float(numoverlaps)/len(events2))
