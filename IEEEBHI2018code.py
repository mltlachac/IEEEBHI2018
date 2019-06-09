#Last upddated Jan 2018

import pandas as pd
import numpy as np
import math
from sklearn import linear_model
from sklearn.svm import SVR
import scipy.stats as stats
df = pd.read_csv('dataset.csv', header = 0)
df.head()


#for number of reports
antibiotic = []
bacteria = []
year1 = []
year2 = []
y1reports = []
y2reports = []
y1mean = []
y2mean = []
meanDiff = []
pchi = []
SEvalue = []
significant = []

notResistant1 = []
resistant1 = []
notResistant2 = []
resistant2 = []
totaly1 = []
totaly2 = []

pTtest = []
Ttest = []

for org in set(df.organism):
    for ab in set(df.component):
        y = 2014
        #for y in range(2006, 2014):
        tdf = df[(df.component == ab) & (df.organism == org)]
        ttdfy1 = tdf[(tdf['Report Year'] == y)]
        ttdfy2 = tdf[(tdf['Report Year'] == y+1)]
        if ((ttdfy1.shape[0]>=1) & (ttdfy2.shape[0]>=1)):
            ttdfy1=ttdfy1.set_index(np.arange(0,ttdfy1.shape[0]))
            ttdfy2=ttdfy2.set_index(np.arange(0,ttdfy2.shape[0]))
            antibiotic.append(ab)
            bacteria.append(org)
            year1.append(y)
            year2.append(y+1)
            y1 = 0
            y2 = 0
            for i in range(0, len(ttdfy1)):
                y1 = y1 + ttdfy1["Total Tests (by organism)"][i] * ttdfy1["Indicator Value (Pct)"][i]
            totaly1.append(sum(ttdfy1["Total Tests (by organism)"]))
            y1mean.append(y1/sum(ttdfy1["Total Tests (by organism)"]))

            for i in range(0, len(ttdfy2)):
                y2 = y2 + ttdfy2["Total Tests (by organism)"][i] * ttdfy2["Indicator Value (Pct)"][i]
            totaly2.append(sum(ttdfy2["Total Tests (by organism)"]))
            y2mean.append(y2/sum(ttdfy2["Total Tests (by organism)"]))

            notResistant1.append(int(sum(ttdfy1["Total Tests (by organism)"]) * ((y1/100)/sum(ttdfy1["Total Tests (by organism)"]))))
            resistant1.append(int(round(sum(ttdfy1["Total Tests (by organism)"]) * (1-((y1/100)/sum(ttdfy1["Total Tests (by organism)"]))))))
            notResistant2.append(int(round(sum(ttdfy1["Total Tests (by organism)"]) * ((y2/100)/sum(ttdfy2["Total Tests (by organism)"])))))
            resistant2.append(int(round(sum(ttdfy1["Total Tests (by organism)"]) * (1-((y2/100)/sum(ttdfy2["Total Tests (by organism)"]))))))

            y1reports.append(ttdfy1.shape[0]) 
            y2reports.append(ttdfy2.shape[0])
            
for n in y1reports:    
    se = 10/(n**(1/2))
    SEvalue.append(se)
for r in range(0, len(SEvalue)):
    diff = y2mean[r] - y1mean[r]
    meanDiff.append(diff)
    significant.append(-SEvalue[r] + abs(diff))


pchi = []
chi = []
for i in range(0, len(notResistant1)):
    dfTy1 = pd.DataFrame(["NotResistant"]*notResistant1[i] + ["Resistant"]*(resistant1[i]))
    dfTy2 = pd.DataFrame(["NotResistant"]*notResistant2[i] + ["Resistant"]*(resistant2[i]))
    expected = pd.crosstab(index=dfTy1[0], columns ="count")
    observed = pd.crosstab(index=dfTy2[0], columns ="count")
    chi_squared_stat = (((observed-expected)**2)/expected).sum()
    if resistant1[i] ==0:
        chi.append((((notResistant1[i]-notResistant2[i])**2)/(notResistant1[i]+0.0000000001)))
    elif notResistant1[i] ==0:
        chi.append((((resistant1[i]-resistant2[i])**2)/(resistant1[i]+0.0000000001)))
    else:
        chi.append((((resistant1[i]-resistant2[i])**2)/(resistant1[i]+0.0000000001)) + (((notResistant1[i]-notResistant2[i])**2)/(notResistant1[i]+0.0000000001)))
    newstat = stats.chi2.cdf(x=chi_squared_stat, df=1)[0]
    pchi.append(1-newstat)


seDF = pd.DataFrame()
seDF["antibiotic"] = antibiotic
seDF["bacteria"] = bacteria
seDF["year1"] = year1
seDF["year2"] = year2
seDF["y1reports"] = y1reports
seDF["y2reports"] = y2reports
seDF["y1samples"] = totaly1
seDF["y2samples"] = totaly2
seDF["y1mean"] = y1mean
seDF["y2mean"] = y2mean
seDF["meanDiff"] = meanDiff
seDF["pchi"] = pchi
seDF["chi"] = chi
seDF["SEvalue"] = SEvalue
seDF["significant"] = significant
seDF.head()


print("Standard Error")
print("Total: " + str(seDF.shape[0]))
print("sig: " + str(seDF[(seDF.significant>0)].shape[0]))
print("sig + neg:")
new = seDF[(seDF.significant>0) & (seDF.meanDiff<0)]
print(new.shape[0])


print("Chi-squared")
print("Total: " + str(seDF.shape[0]))
print("sig: " + str(seDF[(seDF.pchi)<=0.05].shape[0]))
print("sig + neg: ")
new = seDF[(seDF.pchi<=0.05)]
print(new.shape[0]) 
print(max(new.meanDiff))
print(min(new.meanDiff))


print("negative: " + str(seDF[(seDF.meanDiff<0)].shape[0]))
tagged = seDF[(seDF.pchi<=0.05) & (seDF.significant>0) & (seDF.meanDiff<0)]
print("tagged: " + str(tagged.shape[0]))


