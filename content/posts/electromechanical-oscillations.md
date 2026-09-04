+++
title = 'Electromechanical Oscillations: Towards Subsynchronous Resonance'
date = 2026-09-04
draft = false
+++

## Introduction

Electromechanical oscillations are common phenomena in power systems that can lead to adverse consequences i.e. resonances. They occur between two oscillators, mechanical and electrical, which exchange energy. 
This paper demonstrates the physical nature behind a particular phenomenon typical of a synchronous generator: the mechanical part is a governor and the electrical part is a series-compensated transmission line. When the electrical resonance frequency is below the nominal system frequency, the difference frequency can interact with a mechanical torsional mode and create subsynchronous resonance. The paper also shows how such a system can be analyzed using a numerical modeling approach.   

## Electromechanical system overview

The figure below illustrates the system mentioned above: the governor is represented as a two-mass inertia model consisting of a turbine (inertia constant $H_t$ and speed $\omega_t$), a synchronous generator (inertia constant $H_g$ and speed $\omega_g$), and an elastic shaft connecting them. A two-mass system (a mechanical oscillator) is essential for the observation of resonance effects. Power system is represented as:

 - generator with internal voltage $E_g\angle\delta$ (an ideal model without losses for simplicity)
 - connected to an infinite bus with voltage $E_{ib}\angle0$
 - via a series-compensated transmission line of inductance $L$, resistance $R$ and capacitor size $C$. 

The capacitor is tuned to a resonance frequency, as explained below. Thus, the electrical system is also an oscillator. 

![Electromechanical system](/images/electromechanical-oscillations/ElectromechanicalSystem.jpg)  
*Figure 1: Electromechanical system used for modeling.*

## Mechanical system model

As mentioned above, the mechanical system is an oscillator that has a torsional mode - natural oscillations of the turbine and the generator masses against each other caused by the elasticity of the shaft. The following set of equations describes the system dynamics: 

$$
\frac{d\omega_t}{dt}=\frac{1}{2H_t}\Big(\frac{P_m}{\omega_t}-T_\text{shaft}\Big)\\\\
\text{ }\\\\
\frac{d\omega_g}{dt}=\frac{1}{2H_g}\Big(-\frac{P_e}{\omega_g}+T_\text{shaft}\Big)\\\\
\text{ }\\\\
T_\text{shaft}=D(\omega_t-\omega_g)+\tau\\\\
\text{ }\\\\
\frac{d\tau}{dt}=K(\omega_t-\omega_g)
$$

where $P_m$ and $P_e$ are the mechanical and electrical power, respectively, $T_\text{shaft}$ is the shaft torque, $D$ is teh damping factor and $K$ is the shaft stiffness. The last two equations are characterized by the difference of the speeds that is an effective indicator of torsional oscillations when opposite variation of the speeds takes place as it is demonstrated further. 

Let us simulate dynamic response of the system after a large disturbance for the following equilibrium point: $P_m=P_e=\tau=0$ p.u. and $\omega_t=\omega_g=1$ p.u. (in such case, all derivatives are zero). The Matlab code with all necessary model parameters realizing this is as follows:

```
function x_der=TwoMassInertia(x,x_der,Pm,Pe,Ht,Hg,K,D)#x(1)=wt,x(2)=wg,x(3)->shaft
  x_der(1)=(Pm/x(1)-(D*(x(1)-x(2))+x(3)))/(2*Ht);
  x_der(2)=(-Pe/x(2)+D*(x(1)-x(2))+x(3))/(2*Hg);
  x_der(3)=K*(x(1)-x(2));
end

dt=1e-4;
t=[0:dt:1.5];
Ht=0.7;
Hg=0.4;
D=5;
K=100;

wt=zeros(1,length(t));
wg=zeros(1,length(t));
Pm=ones(1,length(t));#step mechanical power input 0->1
Pe=Pm;
x=[1;1;0];
x_der=x*0;
for ii=1:length(t)
  wt(ii)=x(1);
  wg(ii)=x(2);
  x_der=TwoMassInertia(x,x_der,Pm(ii),Pe(ii),Ht,Hg,K,D);
  x=x+dt*x_der;#Euler method
endfor
```

The resulting dynamic response in terms of the speeds is shown in the figure below. As it was mentioned above, the turbine and generator speeds oscillate in opposite directions. It is also possible to determine natural oscillation frequency from the plot - it is about 2.1 Hz.

![Mechanical system step response](/images/electromechanical-oscillations/MechanicalSystemStepResponse.jpg)  
*Figure 2: Mechanical system step response.*

Note that explicit Euler integration is sufficient for demonstration purposes, although more advanced integration methods would generally be preferred for accurate power-system simulations.

### Small signal stability analysis

There is a direct equation for calculating the natural oscillation frequency of such system; however, the paper demonstrates how to use small signal stability analysis to determine it by applying numerical approach. The equations above are not a standard linear state-space model $\dot{x}=Ax+Bu$ because the states (the speeds) are also in denominators. In such case, linearization around an equilibrium point is required. For a large (and sometimes varying) system, this task can be cumbersome; therefore, the paper demonstrates the following approach to find $A$-matrix: 

$$
\frac{dx}{dt}=g=Ax+Bu\\\\
\frac{dg}{dx}=A
$$

Thus, the $A$-matrix can be evaluated numerically as:

$$
A\approx\frac{g(x+\epsilon)-g(x-\epsilon)}{2\epsilon}
$$

The Matlab code implementing this and calculating Eigenvalues is:

```
Pe=1;Pm=1;
x0=[1;1;1];#equilibrium point for given Pe and Pm after disturbance above
x0_der=[0;0;0];
eps=1e-6;
A=zeros(3,3);
for ii=1:3
  x0_rt=x0;x0_rt(ii)=x0_rt(ii)+eps;#states right step
  x0_lt=x0;x0_lt(ii)=x0_lt(ii)-eps;#states left step
  df1_rt=TwoMassInertia(x0_rt,x0_der,Pm,Pe,Ht,Hg,K,D);#state derivatives or g-function for right step
  df1_lt=TwoMassInertia(x0_lt,x0_der,Pm,Pe,Ht,Hg,K,D);#state derivatives or g-function for left step
  df1_dt=(df1_rt-df1_lt)/(2*eps);
  A(:,ii)=df1_dt;
endfor
lambda=eig(A)
```

The resulting Eigenvalues are:

$$
\lambda=
\begin{bmatrix}
0\\\\
-4.643+j13.19 \\\\
-4.643-j13.19
\end{bmatrix} 
$$

Thus, the natural oscillation frequency is $13.19/(2\pi)=2.1$ Hz which matches the value obtained from the simulation. The zero eigenvalue corresponds to the absence of a restoring torque for the absolute rotor angle, while the complex pair represents the torsional mode.

## Electrical system model

The electrical system is described by the following system of equations:

$$
u_g=\sqrt{2}E_g\sin{\delta}\\\\
\text{ }\\\\
u_{ib}=\sqrt{2}E_{ib}\sin{\omega_et}\\\\
\text{ }\\\\
\frac{di_L}{dt}=\frac{1}{L}\Big(u_g-u_{ib}-u_C-i_LR\Big)\\\\
\text{ }\\\\
\frac{du_C}{dt} = \frac{i_L}{C}\\\\
\text{ }\\\\
\frac{d\delta}{dt}=\omega_g\omega_e
$$

where $i_L$, $u_C$ and $\delta$ are states, $\omega_e$ is the base electrical angular frequency ($\omega_g$ is in per units here). The Matlab function implementing this model is:

```
function [x_der,Pe,Pf]=PowerSystem(x,x_der,Eg,R,L,C,we,uib,Pf)#x(4)=iL,x(5)=uC,x(6)=delta
  ug=sqrt(2)*Eg*sin(x(6));
  Pf(1:end-1)=Pf(2:end);Pf(end)=ug*x(4);
  Pe=mean(Pf);
  x_der(4)=(ug-uib-x(5)-x(4)*R)/L;
  x_der(5)=x(4)/C;
  x_der(6)=x(2)*we;
end
```

Note that one of the outputs of this function is electrical power $P_e$ that is used in the inertia model. Since input voltages are instantaneous, $P_e$ is calculated as a mean value over samples in one period. This function can be used to determine resonance frequency applying the approach shown above for Eigenvalues calculations and it coincides with the well-known equation $1/(2\pi\sqrt{LC})$.

### Electrical system parameters

The transmission line parameters are chosen such that with $E_g=E_{ib}=1$ p.u., electrical power of about 0.1 p.u. is achievable at relatively small $\delta$. Thus, inductance was chosen 10 mH which corresponds to a reactance of about 3 p.u. (with base impedance 1 Ohm).

**Important:** Capacitor should be tuned to electrical resonance frequency below the system nominal 50 Hz or 50 - 2.1 = 47.9 Hz to observe electromechanical resonance and 47.5 Hz was chosen (this is why it is called subsynchronous). The strongest response is observed when the subsynchronous frequency approaches the mechanical torsional mode frequency.

The resistance value should not be too high to avoid fast oscillation damping and quality factor 27 was chosen. 

The key idea behind subsynchronous resonance is that the electrical resonance frequency does not interact with the mechanical system directly. Instead, it produces a subsynchronous component whose frequency is equal to the difference between the nominal network frequency and the electrical resonance frequency. For example, a 47.5 Hz electrical resonance in a 50 Hz system creates a 2.5 Hz component, which is close to the 2.1 Hz torsional mode of the mechanical system. As these frequencies approach each other, energy exchange between the electrical and mechanical oscillators becomes more effective, resulting in amplification of torsional oscillations.

## Electromechanical system model and test results

We can now investigate dynamic responses of the whole system for different scenarios. The Matlab script below serves this purpose:

```
we=50*2*pi;#50 Hz system
t=[0:dt:10];
L=10e-3;
fr=47.5;#electrical resonance frequency
C=1/((fr*2*pi)^2*L);
R=fr*2*pi*L/27;
Eg=1;#generator
Eib=1;#infinite bus
uib=sqrt(2)*Eib*sin(we*t);#infinite bus-instantaneous

x=[1;1;0;0;0;0];#state vector
x_der=zeros(length(x),1);
wt=zeros(1,length(t));
wg=zeros(1,length(t));
Pm=ones(1,length(t))*0.1;#constant mechanical power 0.1pu
Pe=zeros(1,length(t));
Pf=zeros(1,round(2*pi/we/dt));
for ii=1:length(t)
  wt(ii)=x(1);
  wg(ii)=x(2);
  [x_der,Pe(ii),Pf]=PowerSystem(x,x_der,Eg,R,L,C,we,uib(ii),Pf);
  x_der=TwoMassInertia(x,x_der,Pm(ii),Pe(ii),Ht,Hg,K,D);
  x=x+dt*x_der;#Euler method
  x(end)=mod(x(end),2*pi);#avoid growth of delta
endfor
```

The figure below shows the resulting speed and power curves for the scenario when the capacitor is excluded from the system. It is possible to see that the transient process after applying mechanical power is oscillatory due to the inertia model, but it gradually attenuates without any presence of electromechanical resonance because the electrical oscillator was removed. Oscillation frequency is about 1 Hz and it is due to excitation of multiple dynamic modes: since there is no electrical resonance, the response is dominated by interactions among the remaining system modes rather than by a single torsional resonance.

![Electromechanical system step response without capacitor](/images/electromechanical-oscillations/ResponseWithoutC.jpg)  
*Figure 3: Electromechanical system step response without capacitor.*

The next figure demonstrates the case where the capacitor is included and it is possible to notice large oscillations. The right plot shows frequency spectrum for the signal $\omega_t-\omega_g$ that has a peak around 2.5 Hz - it corresponds to the coupled electromechanical resonance frequency. Note that the coupled electromechanical mode does not necessarily coincide exactly with the natural mechanical mode because interaction with the electrical network shifts the resulting resonance frequency.   

![Electromechanical system step response with capacitor](/images/electromechanical-oscillations/ResponseWithC.jpg)  
*Figure 4: Electromechanical system step response with capacitor.*

It is possible to prove that there are mutual oscillations between the mechanical and electrical systems by changing capacitor size or electrical resonance frequency used for its calculation. The table below demonstrates that the peak magnitude is growing and has maximum at frequency corresponding to the 3 Hz electromechanical resonance frequency. For a 50 Hz system, electrical resonance frequencies around 47 Hz produce subsynchronous components near 3 Hz, which are close to the torsional mode and therefore lead to stronger interaction.  

fr Hz |	Magnitude p.u.
--- | ---
46 |	0.00197
46.5 |	0.002
47 |	0.0054
47.5 |	0.0019
48 |	0.00185

## Conclusion

The paper demonstrates a model that reproduces the essential electromechanical resonance mechanism. The amplitude of the torsional modes (2.1 Hz) depends strongly on the electrical resonance frequency determined by the RLC network. When the electrical resonance approaches the expected subsynchronous frequency (47 Hz), the magnitude of the torsional mode increases significantly, demonstrating energy exchange between the electrical and mechanical systems. It is also worth noting that the peak disappears with high resistance. The presented model therefore captures the essential physics behind subsynchronous resonance while remaining sufficiently simple for numerical experimentation and educational purposes.

