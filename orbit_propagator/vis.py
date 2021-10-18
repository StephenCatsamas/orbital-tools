import matplotlib.pyplot as plt
from math import *
import csv

t = []
x = []
y = []

with open('vis.dat', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        t.append(float(row[0]))
        x.append(float(row[1]))
        y.append(float(row[2]))


fig, ax = plt.subplots()
ax.plot(x,y)


plt.show()