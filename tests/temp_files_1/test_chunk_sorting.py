import pybedtools as bt

file = bt.BedTool("/home/iprada/mammalian_circles/faststorage/test_inigo/temp_files_63594/peaks.bed")


chunks = 3200

list = []
for i in range(0,chunks):
    list.append([])
print(list)

counter = 0
for interval in file:
    if counter == 3200:
        counter = 0

    if int(interval[2])-int(interval[1])>500:
        w_start = int(interval[1])
        while w_start < int(interval[2]):
            splitted = [interval.chrom,str(w_start),str(w_start+300)]
            w_start+=300
            list[counter].append(splitted)
    else:
        list[counter].append([interval.chrom,str(interval.start),str(interval.end)])
    counter +=1



print(list)
