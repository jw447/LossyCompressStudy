import numpy as np

error = [1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]

def nodecount_estimate(fname):

    nodecount = []
    qf0 = []
    with open(fname+"-1e-0-qf.txt") as inputfile:
        for line in inputfile:
            qf0.append(float(line.strip('\n')))

    s = set()
    for j in qf0:
        if int(j) != 0:
            if int(j) in s:
                continue
            else:
                s.add(int(j))
    nodecount.append(len(s))

    for i in range(0,len(error)-1):
        qf0=[]
        f_name0 = fname+"-1e-"+str(i)+"-qf.txt"

        with open(f_name0) as inputfile:
            for line in inputfile:
                qf0.append(float(line.strip('\n')))

        qf1=[]
        f_name1 = fname+"-1e-"+str(i+1)+"-qf.txt"

        with open(f_name1) as inputfile:
            for line in inputfile:
                qf1.append(float(line.strip('\n')))

        s = set()
        for j in qf1:
            if int(j) != 0:
                if int(j) in s:
                    continue
                else:
                    s.add(int(j))
        nodecount.append(len(s))


        tmp_data = []
        for j in range(0,len(qf1)-2):
            if(qf0[j]!=0 and qf1[j]==0):
                tmp_data.append(qf0[j])

        qf0 = np.array(qf0)
        qf0 = qf0[qf0.nonzero()]
        qf1 = np.array(qf1)
        qf1 = qf1[qf1.nonzero()]

        v_esti = (np.sum(np.power(qf0-np.mean(qf0),2))-np.sum(np.power(np.array(tmp_data)-np.mean(qf0),2)))/len(qf1)
        qf_esti = np.random.normal(loc=np.mean(qf0), scale=np.sqrt(v_esti), size=len(qf1))

    return nodecount
