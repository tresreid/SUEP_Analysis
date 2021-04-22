import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

mind_dr =[]
print("opening")
with open('min_dR.txt') as f:
    for line in f.readlines():
      cols = line.rstrip().split(' ')
      mind_dr.append({"min_dR":float(cols[0])})

print("panda")
df = pd.DataFrame(mind_dr)

print("plotting")
fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Plot title")    
ax1.set_xlabel('x label')
ax1.set_ylabel('y label')
plt.yscale('log')
ax1.hist(df["min_dR"],bins=50,range=(0,0.3))

leg = ax1.legend()

plt.show()
print("done")
