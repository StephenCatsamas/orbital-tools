import matplotlib.pyplot as plt
from math import *
import csv

t = []
x = []

with open('vis.dat', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        t.append(float(row[0]))
        x.append(float(row[1]))


fig, ax = plt.subplots()
ax.plot(t, x)
xe = [exp(ta) for ta in t]
ax.plot(t, xe)

plt.show()