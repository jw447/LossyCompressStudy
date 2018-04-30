from matplotlib import rc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
import math

#fname = ["astro", "blast2_p", "bump", "dpot", "eddy", "fish", "sedov_p", "yf17_p", "yf17_t"]
error = [1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0] # point-wise relative error bound
intvCapacity = 2097152

def hitratio_estimate(fname):

    hitratio_o = []

    with open(fname+"-hitratio.txt") as inputfile:
        for line in inputfile:
            hitratio_o.append(float(line.strip('\n')))

    hitratio_o=np.flip(hitratio_o,0)

    data=[]
    with open(fname+'-data.txt') as inputfile:
        for line in inputfile:
            data.append(float(line.strip('\n')))

    hitratio_esti = []

    for i in range(0,12):

        realPrecision = []
        for k in range(0,len(data)):
            if(k%32==0):
                if(data[k] !=0):
                    realPrecision.append(data[k]*error[i])
                elif data[k] == 0:
                    realPrecision.append(error[i])

        count = 0
        miss = 0

        prediction_esti = np.zeros(len(data))
        prediction_esti[0] = data[0]
        prediction_esti[1] = data[1]

        prederror_esti = np.zeros(len(data))

        h = 0
        realP = realPrecision[h]
        checkRadius = (intvCapacity-1)*realP
        interval = 2*realP
        h = h+1

        for j in range(2,len(data)):
            if(j % 32 == 0):

                realP = realPrecision[h]
                checkRadius = (intvCapacity-1)*realP
                interval = 2*realP
                h = h+1

            prediction_esti[j] = 2*data[j-1]-data[j-2]
            prederror_esti[j] = data[j] - prediction_esti[j]

            if (abs(prederror_esti[j]) < checkRadius): # Hit
                count = count+1
                state = (prederror_esti[j]/realP+1)/2

                if(data[j] >= prediction_esti[j]):
                    prediction_esti[j] = prediction_esti[j] + state*interval

                elif (data[j] < prediction_esti[j]):
                    prediction_esti[j] = prediction_esti[j] - state*interval

            else: # Miss
                prediction_esti[j]=(data[j])

        hitratio_esti.append(count/len(data))
    return hitratio_esti, hitratio_o
