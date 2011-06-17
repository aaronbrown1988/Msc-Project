#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

steps = [200,2000,2000]
colors= ['r','y','b','k']
fig = plt.figure();
ax = fig.add_subplot(111);
x = 0
for  i in  steps:
	filename = "/home/phrhbo/Model/output/jarzinski/sweeps/12x12/5.95/1e5/"+`i`+"/out-curves-0.075-5.95.tsv"
	print filename
	(pf, wpf, pr, wpr) = np.loadtxt(filename, unpack=True)
	ax.scatter(pr,wpr,c=colors[x])
	ax.scatter(pf,wpf,c=colors[x])
	x = x +1
	
plt.show()
