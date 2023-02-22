'''Running toggle switch simulations'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from eqs import *



pulse1 = []
pulse2 = []


t = np.linspace(0,1000,10000)
k = [
    0.17, #ktranscription
    1000, #Ks
    0.05, #Kd
    0.19, #mRNAdeg
    7.85, #ktranslation
    0.05, #repdeg
    2,   #N-HillFunction-REP
    1    #N-HillFunction-IND
]

y0 = [0, 0, 0, 0]
Y = [80, 80]
X1 = 50000
X2 = 0

indconc = 5e05

pulse_args1 = (200, 200, indconc, pulse1)
pulse_args2 = (600, 200, indconc, pulse2)


args = (k, Y, IND_pulse, pulse_args1, pulse_args2)

alt_args = (k, Y, X1, X2)


######################
### RUN ##
######################



timeseries = odeint(pulsed_toggle, y0, t, args=(args))


#timeseries = odeint(full_toggle, y0, t, args=(alt_args))


#plt.plot(t, timeseries)
#plt.show()



fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(5,5), constrained_layout=True)

ax1.plot(t, timeseries[:,0], color='r'    )
ax4.plot(t, timeseries[:,1], color='g'    )
ax3.plot(t, timeseries[:,2], color='b'    )
ax2.plot(t, timeseries[:,3], color='black')


ax1.set_title('mRNA - LacI', fontsize=6)
ax2.set_title('LacI'       , fontsize=6)
ax3.set_title('mRNA - TetR', fontsize=6)
ax4.set_title('TetR'       , fontsize=6)



plt.xlabel('Time', loc='center')
plt.savefig('/Users/kushaaldesai/Documents/Code/toggle/Toggle.png')