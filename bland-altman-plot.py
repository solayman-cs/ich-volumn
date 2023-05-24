# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:07:00 2018
# Code for creating Bland Altman plot and Concordance Correlation Coefficient plots
@author: David
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
from numpy import polyfit

example_text = r'Example Usage: python bland_altman.py --fname=C:\Users\David\Desktop\Plaque_burden_data.csv \
--minlimit=0 --maxlimit=100 --var1="Computed Plaque Burden" --var2="Measured Plaque Burden"'

parser = argparse.ArgumentParser(description=\
'Create Bland-Altman and Concordance correlation coefficient plots', \
epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--fname', help='Enter filename')
parser.add_argument('--minlimit', type=int, help='Enter min limits')
parser.add_argument('--maxlimit', type=int, help='Enter max limits')
parser.add_argument('--var1', help='name of variable1')
parser.add_argument('--var2', help='name of variable2')


#optional arguments for x and y labels
args = parser.parse_args()
data = pd.read_csv(args.fname)

xlim = [args.minlimit, args.maxlimit]
ylim = [args.minlimit, args.maxlimit]

# drop an rows that contain nans
data = data.dropna(axis=0, how='any')

# determine mean and difference
mean = (data.iloc[:, 0] + data.iloc[:, 1])/2
diff = data.iloc[:, 0] - data.iloc[:, 1]

# get mean and standard deviation of the differences
mean_diff = diff.mean()
std_diff = diff.std()

# determine 95% limits of agreement
LoA_plus = mean_diff + 1.96*std_diff
LoA_neg = mean_diff - 1.96*std_diff

# plot figure
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.scatter(mean, diff, s=4)
ax1.set_xlim((xlim[0], xlim[1]))
xval = ax1.get_xlim()
ax1.plot(xval, (mean_diff, mean_diff))
ax1.plot(xval, (LoA_plus, LoA_plus), 'r--')
ax1.plot(xval, (LoA_neg, LoA_neg), 'r--')
# set text
ax1.text(xval[1] - (xval[1]/8), LoA_plus + xval[1]*0.015, '+1.96SD')
ax1.text(xval[1] - (xval[1]/8), LoA_neg + xval[1]*0.015, '-1.96SD')
ax1.text(xval[1] - (xval[1]/8), mean_diff + xval[1]*0.015, 'MEAN')
ax1.text(xval[1] - (xval[1]/8), LoA_plus - xval[1]*0.035, round(LoA_plus, 2))
ax1.text(xval[1] - (xval[1]/8), LoA_neg - xval[1]*0.035, round(LoA_neg, 2))
ax1.text(xval[1] - (xval[1]/8), mean_diff - xval[1]*0.035, round(mean_diff, 2))
if args.var1:
    var1 = args.var1
    var2 = args.var2
    ax1.set_xlabel('Mean of ' + var1 + ' and ' + var2)
    ax1.set_ylabel('Difference between ' + var1 + ' and ' + var2)

# concordance correlation coefficient
covar = data.cov()
ccc = 2*covar.iloc[0, 1]/(covar.iloc[0, 0] + covar.iloc[1, 1] + \
                  (data.iloc[:, 0].mean() - data.iloc[:, 1].mean())/2)

# determine fit
coeff = polyfit(data.iloc[:, 0], data.iloc[:, 1], 1)
print('Coefficients are {}'.format(coeff))

ax2 = fig.add_subplot(122)
ax2.scatter(data.iloc[:, 0], data.iloc[:, 1], s=4)
ax2.plot((xlim[0],xlim[1]), (coeff[0]*0 + coeff[1], coeff[0]*ylim[1] + coeff[1]))
# plot concordance line
ax2.plot((xlim[0],xlim[1]), (ylim[0], ylim[1]), 'k--')
ax2.text(xlim[0]+(xlim[1]/20), xlim[1]-(xlim[1]/10), 'CCC={}'.format(round(ccc, 2)))
ax2.set_xlim((xlim[0], xlim[1]))
ax2.set_ylim((ylim[0], ylim[1]))
if args.var1:
    var1 = args.var1
    var2 = args.var2
    ax2.set_xlabel(var1)
    ax2.set_ylabel(var2)
plt.show()