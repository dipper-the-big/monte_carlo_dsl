import matplotlib.pyplot as plt
import prettyplotlib as ppl

f = open("lakimoca_data5.txt")
val = [x.split() for x in f]
val = val[:-1]
xVal = [float(x[1]) for x in val]
nI = [int(x[3]) for x in val]
nIC = [int(x[4]) for x in val]
nV = [int(x[5]) for x in val]
nVC = [int(x[6]) for x in val]
nIs = [x + y for (x, y) in zip(nI, nIC)]
nVs = [x + y for (x, y) in zip(nV, nVC)]

ppl.plot(xVal, nIs, '-', label="All interstitials", linewidth = 1.0)
ppl.plot(xVal, nIC, '-.', label="Interstitials in cluster", linewidth = 1.25)
ppl.plot(xVal, nVs, '--', label="All vacancies", linewidth = 1.75)
ppl.plot(xVal, nVC, ':', label="Vacancies in cluster", linewidth = 1.75)
plt.xscale('log')
plt.xlabel("Time (s)")
plt.ylabel("Number of defects")
ppl.legend(shadow=True, fancybox=True, ncol=2)
plt.savefig("lakimocaTrendsAutoSave.eps", format="eps", dpi=1200)
plt.show()
plt.close()
