################################################################
# libraries:
import re
import time
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

################################################################
# functions:
    
from sa import SuffixArray
# from PrefixDoubling import SuffixArray  # Cython converted version.
from sa import binary_search
from SEQsimulator import simulate_string
from SEQsimulator import get_exact_read

def SA_read_mapper(string, pattern):
    SA = SuffixArray(string)
    all_matches = binary_search(SA, string, pattern)
    return all_matches

def re_find(string, pattern):
    re_findings = [m.start() for m in re.finditer('(?={0})'.format(re.escape(pattern)), string)]
    return re_findings 

################################################################
# tests:

# 1 mio test:
string = simulate_string(100000)
start_time = time.time()
SuffixArray(string)
end_time = time.time()
print(end_time-start_time)


# Test suffixtree-algorithm vs naive-algorithm for same result:
for i in range(50000):
    print('Iteration nr: ', i+1)
    ref = simulate_string(random.randint(40,100))
    read = get_exact_read(ref, random.randint(1,30))
    if sorted(SA_read_mapper(ref, read)) != sorted(re_find(ref,read)):
        print('Algorithm mistake!')
        print(ref)
        print(read)
        break
    if i+1 == 50000:
        print('DONE')


# Runtimes for the array construction (varying ref lengths):
ref_lengths = [25000,50000,75000,100000,125000,150000,175000]
runtimes = []
for idx in range(7):
    print('Iteration nr: ', idx+1) 
    replicate = []
    for j in range(10):
        ref = simulate_string(ref_lengths[idx])
        start_time = time.time()
        SuffixArray(ref)
        end_time = time.time()
        replicate.append(end_time-start_time)
    runtimes.append(np.mean(replicate))
# plot running times:
# for i in range(len(runtimes)):
#     runtimes[i] = runtimes[i]/math.log10(ref_lengths[i])**2
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes, ax=ax)
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()


# Runtimes for mapping (varying read lengths):
ref_lengths = [2500,5000,7500,10000,12500]
read_lengths_10 = [100]*5
read_lengths_20 = [200]*5
read_lengths_30 = [300]*5
read_lengths_40 = [400]*5
read_lengths_50 = [500]*5
runtimes_10 = []
runtimes_20 = []
runtimes_30 = []
runtimes_40 = []
runtimes_50 = []
for idx in range(5):
    print('Iteration nr: ', idx+1)
    runtimes_10_replicate = []
    runtimes_20_replicate = []
    runtimes_30_replicate = []
    runtimes_40_replicate = []
    runtimes_50_replicate = []
    
    for i in range(100):
        ref = simulate_string(ref_lengths[idx])
        SA = SuffixArray(ref)
        
        read_10 = get_exact_read(ref, read_lengths_10[idx])
        read_20 = get_exact_read(ref, read_lengths_20[idx])
        read_30 = get_exact_read(ref, read_lengths_30[idx])
        read_40 = get_exact_read(ref, read_lengths_40[idx])
        read_50 = get_exact_read(ref, read_lengths_50[idx])
        
        # Dont know why it is nessesary to read run this chunk of code, but 
        # if i dont the first the first runtime is affected. Mayde something 
        # to do with loading the modules??
        binary_search(SA, ref, read_10)
        
        start_time = time.time()
        binary_search(SA, ref, read_10)
        end_time = time.time()
        runtimes_10_replicate.append(end_time-start_time)
        
        start_time = time.time()
        binary_search(SA, ref, read_20)
        end_time = time.time()
        runtimes_20_replicate.append(end_time-start_time)
        
        start_time = time.time()
        binary_search(SA, ref, read_30)
        end_time = time.time()
        runtimes_30_replicate.append(end_time-start_time)
        
        start_time = time.time()
        binary_search(SA, ref, read_40)
        end_time = time.time()
        runtimes_40_replicate.append(end_time-start_time)
        
        start_time = time.time()
        binary_search(SA, ref, read_50)
        end_time = time.time()
        runtimes_50_replicate.append(end_time-start_time)
        
    runtimes_10.append(np.mean(runtimes_10_replicate))
    runtimes_20.append(np.mean(runtimes_20_replicate))
    runtimes_30.append(np.mean(runtimes_30_replicate))
    runtimes_40.append(np.mean(runtimes_40_replicate))
    runtimes_50.append(np.mean(runtimes_50_replicate))

# plot running times:
runtimes = [runtimes_10,runtimes_20,runtimes_30,runtimes_40,runtimes_50]
for run in runtimes:
    for i in range(len(run)):
        run[i] = run[i]/math.log2(ref_lengths[i])
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes_50, ax=ax, label='read length = 50')
sns.lineplot(x=ref_lengths, y=runtimes_40, ax=ax, label='read length = 40')
sns.lineplot(x=ref_lengths, y=runtimes_30, ax=ax, label='read length = 30')
sns.lineplot(x=ref_lengths, y=runtimes_20, ax=ax, label='read length = 20')
sns.lineplot(x=ref_lengths, y=runtimes_10, ax=ax, label='read length = 10')
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()

################################################################
