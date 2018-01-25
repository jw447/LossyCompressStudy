import os
import numpy as np
import matplotlib.pyplot as plt

files = os.listdir()
# cr ci br bi pr pi pf o
# eddy_l = [s for s in files if "eddy" in s]
eddy_l = ["eddy_ci.dat.log",
		"eddy_bi.dat.log",
		"eddy_br.dat.log",
		"eddy_cr.dat.log",
		"eddy_pi.dat.log",
		"eddy_pf.dat.log",
		"eddy_pr.dat.log",
		"eddy_o.dat.log"
		]

# fish_l = [s for s in files if "fish" in s]
fish_l = ["fish_ci.dat.log",
		"fish_bi.dat.log",
		"fish_br.dat.log",
		"fish_cr.dat.log",
		"fish_pi.dat.log",
		"fish_pf.dat.log",
		"fish_pr.dat.log",
		"fish_o.dat.log"
		]

# yf17_l = [s for s in files if "yf17" in s]
yf17_l = ["yf17_p_ci.dat.log",
		"yf17_p_bi.dat.log",
		"yf17_p_br.dat.log",
		"yf17_p_cr.dat.log",
		"yf17_p_pi.dat.log",
		"yf17_p_pf.dat.log",
		"yf17_p_pr.dat.log",
		"yf17_p_o.dat.log"
		]



for i in range(1,4):
	for j in range(5,9):
		if (i==1): #fish
			if (j != 8):
				plt.subplot(4,3,(j-4)*3-2)
				quant = np.loadtxt(fish_l[j-1]);
				plt.hist(quant,200)
				plt.title('Fish Sample-'+fish_l[j-1][5:7])
				plt.xlim(0,10000)
				plt.ylim(0,10)
			elif j == 8:
				plt.subplot(4,3,(j-4)*3-2)
				quant = np.loadtxt(fish_l[j-1]);
				plt.hist(quant,500)
				plt.title('Fish Original')
				plt.xlim(0,10000)
				plt.ylim(0,200)
		if (i==2):#eddy
			if (j != 8):
				plt.subplot(4,3,(j-4)*3-1)
				quant = np.loadtxt(eddy_l[j-1]);
				plt.hist(quant,1000)
				plt.title('Eddy Sample-'+eddy_l[j-1][5:7])
				plt.xlim((0,65535))
				plt.ylim(0,20)
			elif j == 8:
				plt.subplot(4,3,(j-4)*3-1)
				quant = np.loadtxt(eddy_l[j-1]);
				plt.hist(quant,4000)
				plt.title('Eddy Original')
				plt.xlim((0,65535))
				plt.ylim(0,1000)
		if (i==3):#yf17
			if (j != 8):
				plt.subplot(4,3,(j-4)*3)
				quant = np.loadtxt(yf17_l[j-1]);
				plt.hist(quant,1000)
				plt.title('Yf17 Sample-'+yf17_l[j-1][5:7])
				plt.xlim(0,8192)
				plt.ylim(0,20)
			elif j == 8:
				plt.subplot(4,3,(j-4)*3)
				quant = np.loadtxt(yf17_l[j-1]);
				plt.hist(quant,4500)
				plt.title('Yf17 Original')
				plt.xlim(0,8192)
				plt.ylim(0,300)

# plt.subplot(2,2,2)
# quant = np.loadtxt(files[1]+".txt");
# plt.hist(quant,200)
# plt.title('YF17 Sample')
# plt.xlim(0,8192)
# plt.ylim(0,20)


# plt.subplot(2,2,3)
# quant = np.loadtxt(files[2]+".txt");
# plt.hist(quant,800)
# plt.title('Fish Origin')
# plt.xlim(0,10000)
# plt.ylim(0,50)


# plt.subplot(2,2,4)
# quant = np.loadtxt(files[3]+".txt");
# plt.hist(quant,800)
# plt.title('YF17 Origin')
# plt.xlim(0,8192)
# plt.ylim(0,400)


# # plt.subplot(2,3,5)
# # quant = np.loadtxt(files[4]+".txt");
# # plt.hist(quant,800)
# # plt.title('Eddy Origin')
# # plt.xlim(0,65535)
# # plt.ylim(0,300)


# # plt.subplot(2,3,6)
# # quant = np.loadtxt(files[5]+".txt");
# # plt.hist(quant,800)
# # plt.title('YF17 Origin')
# # plt.xlim(0,8192)
# # plt.ylim(0,400)

plt.show()

