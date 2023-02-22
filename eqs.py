'''Modelling the action of the genetic toggle switch'''


#1 Modelling the half toggles
import numpy as np



def full_toggle(y, t, k, Y, X1, X2):

    RNAconc1 =y[0]
    REPconc1 =y[1]
    RNAconc2=y[2]
    REPconc2 =y[3]

    print(t, X1, X2)

    DNAconc1 = Y[0]
    DNAconc2 = Y[1]

    dMRNA1dt = (k[0]*DNAconc1/(1+(REPconc2/(1+((X1/k[1])**k[7])*k[2]))**k[6]))-(k[3]*RNAconc1)
    dREP1dt =  (k[4]*RNAconc2) - (k[5]*REPconc1)

    dMRNA2dt = (k[0]*DNAconc2/(1+(REPconc1/(1+((X2/k[1])**k[7])*k[2]))**k[6]))-(k[3]*RNAconc2)
    dREP2dt =  (k[4]*RNAconc1) - (k[5]*REPconc2)

    return [dMRNA1dt, dREP1dt, dMRNA2dt, dREP2dt]


def toggle_constmARi(y, t, k, Y, X1, X2):

    RNAconc1 =y[0]
    REPconc1 =y[1]
    RNAconc2=y[2]
    REPconc2 =y[3]

    print(t, X1, X2)

    DNAconc1 = Y[0]
    DNAconc2 = Y[1]

    dMRNA1dt = (k[0]*DNAconc1/(1+(REPconc2/(1+((X1/k[1])**k[7])*k[2]))**k[6]))-(k[3]*RNAconc1) - (k[8]*RNAconc1) 
    dREP1dt =  (k[4]*RNAconc2) - (k[5]*REPconc1)

    dMRNA2dt = (k[0]*DNAconc2/(1+(REPconc1/(1+((X2/k[1])**k[7])*k[2]))**k[6]))-(k[3]*RNAconc2) - (k[8]*RNAconc2) 
    dREP2dt =  (k[4]*RNAconc1) - (k[5]*REPconc2)

    return [dMRNA1dt, dREP1dt, dMRNA2dt, dREP2dt]



def IND_pulse(t, t_0, tau, x_0, pulse):
    
    pulse.append(np.logical_and(t >= t_0, t <= (t_0 + tau)) * x_0)

    return np.logical_and(t >= t_0, t <= (t_0 + tau)) * x_0


  


def pulsed_toggle(y, t, k, Y, x_fun, pulse_args1, pulse_args2):

    X1 = x_fun(t, *pulse_args1)

    X2 = x_fun(t, *pulse_args2)

    solution = full_toggle(y, t, k, Y, X1, X2)

    return solution


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import numpy as np

    pulse_args1 = (2, 2, 5e05)
    pulse_args2 = (6, 2, 5e05)


    t = np.linspace(0,10,1000)

    x = IND_pulse(t, *pulse_args1)
    x2 = IND_pulse(t, *pulse_args2)




    plt.plot(t, x)
    plt.plot(t, x2)
    plt.show()


