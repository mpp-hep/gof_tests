

#     xp = np.linspace(x[0], x[-1])

    fig, ax1 = plt.subplots(figsize=(10,10))

    asymmetric_error = [y_down, y_up]

    ax1.errorbar(x, y, yerr=asymmetric_error, fmt='-')
    ax1.set_title(qstring)
#     ax1.set_yscale('log')
    coefficients = np.polyfit(x, y, 5, w=1/(y_up+y_down))
    coeffs[qr] = coefficients
    print(qstring)
    print("np.poly1d({})".format(coefficients))
    poly = np.poly1d(coefficients)
    new_y = poly(x)
    plt.plot(x, y, "o", x, new_y,color="red")
    fig.savefig('fit_xs100_{}.pdf'.format(qr))
    fig2, ax2 = plt.subplots(figsize=(10,2))
    ratio = y/new_y
    plt.plot(x, ratio, "o",color="blue")
    plt.ylim([0.99,1.01])
    plt.show()
    fig2.savefig('ratio_xs100_{}.pdf'.format(qr))

for qr in qrs:
    print("{} : np.poly1d({}),".format(1.-float(qr)/100,coeffs[qr]))

f0 = open('/eos/project/d/dshep/TOPCLASS/DijetAnomaly/QR_models/envelope/cut_stats.json',)
f100 = open ('cut_stats_allQ_GtoWW35na.json',)
qr_unc0 = json.load(f0)
qr_unc100 = json.load(f100)

qrs = ['01','10','90']

for qr in qrs:
    qstring = 'q{}'.format(qr)
    xs0_x      = np.array([i[0] for i in qr_unc0[qstring]])
    xs0_y      = np.array([i[1] for i in qr_unc0[qstring]])
    xs0_y_down = np.fabs(xs0_y-np.array([i[3] for i in qr_unc0[qstring]]))
    xs0_y_up   = np.fabs(xs0_y-np.array([i[4] for i in qr_unc0[qstring]]))
    xs0_x      = xs0_x     [:35] 
    xs0_y      = xs0_y     [:35] 
    xs0_y_down = xs0_y_down[:35] 
    xs0_y_up   = xs0_y_up  [:35] 
    
    xs100_x      = np.array([i[0] for i in qr_unc100[qstring]])
    xs100_y      = np.array([i[1] for i in qr_unc100[qstring]])
    xs100_y_down = np.fabs(xs100_y-np.array([i[3] for i in qr_unc100[qstring]]))
    xs100_y_up   = np.fabs(xs100_y-np.array([i[4] for i in qr_unc100[qstring]]))
    xs100_x      = xs100_x     [:35] 
    xs100_y      = xs100_y     [:35] 
    xs100_y_down = xs100_y_down[:35] 
    xs100_y_up   = xs100_y_up  [:35] 
    

    fig, ax1 = plt.subplots(figsize=(10,10))

    xs0_asymmetric_error = [xs0_y_down, xs0_y_up]
    xs100_asymmetric_error = [xs100_y_down, xs100_y_up]

    ax1.errorbar(xs0_x, xs0_y, yerr=xs0_asymmetric_error, fmt='-',label='0 fb inj.')
    ax1.errorbar(xs100_x, xs100_y, yerr=xs100_asymmetric_error, fmt='-',label='100 fb inj.')
    ax1.legend(loc='upper left')
    ax1.set_title(qstring)
    fig.savefig('qrRatio{}.pdf'.format(qr))
    fig2, ax2 = plt.subplots(figsize=(10,2))
    ratio = xs0_y/xs100_y
    plt.plot(xs0_x, ratio, "o",color="blue")
    plt.ylim([0.99,1.01])
    plt.show()
    fig2.savefig('qrRatio{}.pdf'.format(qr))
  
