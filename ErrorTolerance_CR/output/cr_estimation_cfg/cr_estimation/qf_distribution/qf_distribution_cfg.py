import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.figure as fig
import math

fname = ["astro","blast2_p", "bump", "dpot", "eddy", "fish", "sedov_p", "yf17_p", "yf17_t"]
error=["1e-11", "1e-10", "1e-9", "1e-8", "1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e-0"]

astro=[]
blast2_p = []
bump=[]
dpot=[]
eddy=[] 
fish=[] 
sedov_p=[] 
yf17_p=[] 
yf17_t=[]

for h in range(0,len(fname)):
    for i in range(0,len(error)):
        data = []
        f_name = fname[h]+"-1e-"+str(i)+"-qf.txt"
        with open(f_name) as inputfile:
            for line in inputfile:
                data.append(float(line.strip('\n')))

        count, bins, patches= plt.hist(data,color='b',bins=1000,density=True,alpha=1,label="histogram of real prediction error")
        plt.xlim((1000000,1100000))
        plt.savefig(fname[h]+str(i)+".png", dpi=150)
        plt.clf()


