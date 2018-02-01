
import pybedtools as bt
import pandas as pd
import time
import numpy as np
import sys

labels = ['chrom', 'start', 'end', 'read']
my_data = [['chr1', 784344, 800125, 'read1'],
           ['chr1', 784344, 800124, 'read2'],
           ['chr1', 784344, 800124, 'read3']]
my_data_pd = pd.DataFrame.from_records(my_data, columns=labels)





my_test_bed = bt.BedTool(my_data)





begin = time.time()
for i in range(0,100):
    sorted = my_test_bed.sort()
    sorted_grouped = sorted.groupby(g=[1,2,3],c=4,o='count_distinct')
end = time.time()

print((end-begin)/60,"time with pybedtools")


time.sleep(2)


#my try

begin = time.time()
for w in range(0,100):

    sortedd  = my_data_pd.groupby(['chrom','start','end'],sort=False).read.nunique()
    my_test = bt.BedTool(sortedd.index.tolist())

end = time.time()
print((end-begin)/60,"time with pandas read unique")

time.sleep(2)


begin = time.time()
for s in range(0,100):
    sortedd = my_data_pd.groupby(['chrom', 'start', 'end'], sort=False).read.nunique()
    L = [x for x in sortedd.reset_index().values.tolist()]
    a = bt.BedTool(L)

end = time.time()
print((end-begin)/60,"time with pandas list comprehension")
time.sleep(2)

begin = time.time()
for d in range(0,100):
    sortedd = my_data_pd.groupby(['chrom', 'start', 'end'], sort=False).read.nunique()
    idx = pd.MultiIndex.from_arrays(sortedd.reset_index().values.T)
    a = bt.BedTool(idx.tolist())

end = time.time()
print((end-begin)/60,"time with pandas Multiindex")
time.sleep(2)


begin = time.time()
for f in range(0,100):
    sortedd = my_data_pd.groupby(['chrom', 'start', 'end'], sort=False).read.nunique()
    a = bt.BedTool(sortedd.to_frame('val').set_index('val', append=True).index.tolist())


end = time.time()
print((end-begin)/60,"time with pandas set_index")





