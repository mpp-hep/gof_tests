import json
import math
import matplotlib.pyplot as plt
import numpy as np

def doSmoothing(infile, qrs = ['01','10','30','50','70','90'], reduceRange=False):
  f = open(infile,)
  qr_unc = json.load(f)

  functions= {}
  for qr in qrs:
    qstring = 'q{}'.format(qr)
    x      = np.array([i[0] for i in qr_unc[qstring]])
    y      = np.array([i[1] for i in qr_unc[qstring]])
    y_down = np.fabs(y-np.array([i[3] for i in qr_unc[qstring]]))
    y_up   = np.fabs(y-np.array([i[4] for i in qr_unc[qstring]]))
    
    if reduceRange:
      x      = x     [:35] 
      y      = y     [:35] 
      y_down = y_down[:35] 
      y_up   = y_up  [:35] 

    fig, axs = plt.subplots(2)

    asymmetric_error = [y_down, y_up]

    axs[0].errorbar(x, y, yerr=asymmetric_error, fmt='-')
    axs[0].set_title(qstring)

    coefficients = np.polyfit(x, y, 5, w=1/(y_up+y_down))
    functions[qr]= np.poly1d(coefficients)
    new_y = poly(x)
    axs[0].plot(x, y, "o", x, new_y,color="red")
    ratio = y/new_y
    axs[1].plot(x, ratio, "o",color="blue")
    fig.savefig('smoothQR_{}_qr{}.pdf'.format(infile.split('/')[-1].replace('.json',''),qr))

  for qr in qrs:
    print("{} : np.poly1d({}),".format(1.-float(qr)/100,coeffs[qr]))
    
  return functions
