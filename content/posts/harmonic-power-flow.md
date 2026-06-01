+++
title = 'Harmonic power flow: analysis of main approaches'
date = 2026-06-01
draft = false
+++

## Introduction

This paper shows algorithms of the two main approaches for harmonic power flow in detail. Such analysis is done in the frequency domain and often implemented in various RMS-type power system simulators. It also becomes of high interest due to the increasing share/adoption of inverter-based resources.  

The first approach is based on representation of harmonics as current sources with the following power flow analysis performed separately for each harmonic. In this procedure, impact of one harmonic on another is ignored. In this article, this approach will be denoted as **decoupled**. Its advantage is built around its reduced complexity and it is easy to implement in software, robust and has low computational expense. It is also sufficiently accurate. However, it does not account for harmonic voltage interactions that can dramatically reduce precision in weak grids, in the presence of strong voltage sensitivity or in heavily distorted waveforms. Resonance effects and interaction between harmonics can also be missed.  

These issues can be resolved (to some extent) by the second approach that will be referred to as **coupled**. Harmonics are coupled - each harmonic depends on the voltages of others. 

In further sections, detailed analysis of each approach is given.

## Test system

The figure below shows the test system used for the analysis. It is kept very simple to obtain concise equations facilitating analysis and to focus on understanding the algorithms. The system consists of the grid ideal voltage source $u_g$, a line with inductance $L$ and non-linear load with current $i_{ld}$ being a function of its voltage $u_{ld}$, i.e. $i_{ld}=f(u_{ld})$.  

![Test system](/images/harmonic-power-flow/elsys.jpg)  
*Figure 1: Test system.*

The load is essentially an inductor with non-linear flux linkage - current relationship. Its characteristics are presented in the figure below. It is possible to see that the current has dominant odd-harmonics. In this article, we will focus only on the 3<sup>rd</sup> and 5<sup>th</sup>.

<img src="/images/harmonic-power-flow/load_characteristics.jpg" width="1200">

*Figure 2: Non-linear load characteristics: left - current dependency on flux linkage, center - voltage and current waveforms, right - current harmonic spectrum.*

In this article, harmonics are calculated from instantaneous quantities $x(t)$ using Discrete Fourier Transform (DFT): 

$$
X_h=\frac{\sqrt{2}\exp(j\frac{\pi}{2})}{N}\sum_{n=1}^N x(n)\exp(-j2\pi\frac{n}{N}h)
$$

where $X_h$ is the $h$-order sin-based RMS harmonic phasor and $N$ is the number of samples in one period. 

### Time domain analysis

To obtain reference values for load voltage harmonics, time domain analysis is performed using the following equations: 

$$
u_g=u_L+u_{ld} \\\\
\lambda=\int u_{ld}dt \\\\
i_{ld}=f(\lambda) \\\\
u_L=L\frac{di_{ld}}{dt}
$$

where $u_L$ the voltage across the line and $\lambda$ is load flux linkage. 

Here is the Matlab script used for time domain analysis with all necessary system parameters:

```
w=2*pi*50;
dt=10e-6;
t=[0:dt:0.04];
u_g=sqrt(2)*sin(w*t+pi/2);
L=0.2/w;
u_L=0;
lambda=0;
u_ld=zeros(1,length(t));
i_ld=zeros(1,length(t));
for ii=2:length(t)
  u_ld(ii)=u_g(ii)-u_L;
  lambda=lambda+0.5*(u_ld(ii-1)+u_ld(ii))*dt;
  i_ld(ii)=sign(lambda)*(9.6539e+06*abs(lambda)^3+2.0879e+04*abs(lambda)^2);
  u_L=L*(i_ld(ii)-i_ld(ii-1))/dt;
endfor
```

The result is the distorted load voltage shown in the figure below. 

<img src="/images/harmonic-power-flow/volt_harm_ref.jpg" width="1200">

*Figure 3: Load voltage in time domain: left - load voltage and current waveforms, center - load voltage harmonic spectrum, right - 1<sup>st</sup>, 3<sup>rd</sup> and 5<sup>th</sup> harmonic in polar coordinates.*

**Note:** grid phase is set to $90^\circ$ to avoid core saturation due to the appearance of a DC offset in flux linkage. The figure below shows as an example of the load characteristics with the grid phase $0^\circ$. The core remains saturated due to absence of active losses in the system. 

![Saturated core](/images/harmonic-power-flow/saturated_core.jpg)  
*Figure 4: Load voltage and current waveforms with saturated core.*

## Decoupled approach

In the decoupled approach, each harmonic is treated independently. This means that harmonic currents are computed based only on the fundamental voltage, and harmonic voltages are then calculated separately without accounting for interaction between different harmonic orders. Thus, the algorithm is:
1. Find the load voltage 1<sup>st</sup> harmonic $U_{ld-h1}$.
2. Find load current 3<sup>rd</sup> and 5<sup>th</sup> harmonics $I_{ld-h3}$ and $I_{ld-h5}$ using $U_{ld-h1}$.
3. Find separately load voltage harmonics $U_{ld-h3}$ and $U_{ld-h5}$ using $I_{ld-h3}$ and $I_{ld-h5}$.

The first step is a standard power flow solution by Newton-Raphson method for the main harmonic with the following function to nullify: 

$$
G=(U_g-U_{ld-h1})Y-I_{ld-h1}
$$

with the grid voltage $U_g$ and the line admittance $Y=(j\omega L)^{-1}$. When we need to find load current 1<sup>st</sup> harmonic $I_{ld-h1}$ from voltage $U_{ld-h1}$, we will switch from frequency domain to time domain (DC offset is not present):

$$
u_{ld}=|U_{ld-h1}|\cdot\sin(\omega t+\angle U_{ld-h1})
$$

then, as it was done above, flux linkage is calculated $\lambda=\int u_{ld}$ and $i_{ld}=f(\lambda)$. From this time waveform, the 1<sup>st</sup> harmonic $I_{ld-h1}$ is obtained using DFT. The corresponding Matlab function is:

```
function I_ld_h1=nonlinear_load_current_main_harmonic(U_ld_h1,w,t,dt)
  N=length(t);
  u_ld=sqrt(2)*abs(U_ld_h1)*sin(w*t+angle(U_ld_h1));
  lambda=zeros(1,N);
  for ii=2:N
    lambda(ii)=lambda(ii-1)+0.5*(u_ld(ii-1)+u_ld(ii))*dt;
  endfor
  i_ld=sign(lambda).*(9.6539e+06*abs(lambda).^3+2.0879e+04*abs(lambda).^2);
  I_ld_h1=mean(i_ld.*exp(-1i*2*pi/N*[1:N]))*sqrt(2)*exp(1i*pi/2);
end
```

To find $U_{ld-h1}$, we use an iterative approach with the following linear equations:

$$
\begin{bmatrix}
\frac{dG_{re}}{dU_{ld-h1-re}} & \frac{dG_{re}}{dU_{ld-h1-im}} \\\\
\\\\
\frac{dG_{im}}{dU_{ld-h1-re}} & \frac{dG_{im}}{dU_{ld-h1-im}}
\end{bmatrix} 
\cdot
\begin{bmatrix}
\Delta U_{ld-h1-re} \\\\
\Delta U_{ld-h1-im}
\end{bmatrix}
\=
\begin{bmatrix}
G_{re} \\\\
G_{im}
\end{bmatrix}
$$

where complex values are split into real $re$ and imaginary $im$ parts. Jacobian (derivatives of $G$) is difficult to obtain analytically, therefore we will approximate it as: 

$$
\frac{dG}{dU}\approx\frac{G(U+d)-G(U-d)}{2d}
$$

with $d=0.001$ that seems reasonable and robust. The corresponding Matlab function is:

```
function Jac=NR_Jacobian(U_ld_h1,w,t,dt,U_g,Y)
  dm=0.001;
  I_ld_h1_re_rt=nonlinear_load_current_main_harmonic(U_ld_h1+dm,w,t,dt);
  I_ld_h1_re_lt=nonlinear_load_current_main_harmonic(U_ld_h1-dm,w,t,dt);
  I_ld_h1_im_rt=nonlinear_load_current_main_harmonic(U_ld_h1+1i*dm,w,t,dt);
  I_ld_h1_im_lt=nonlinear_load_current_main_harmonic(U_ld_h1-1i*dm,w,t,dt);

  g_re_rt=(U_g-(U_ld_h1+dm))*Y-I_ld_h1_re_rt;
  g_re_lt=(U_g-(U_ld_h1-dm))*Y-I_ld_h1_re_lt;
  g_im_rt=(U_g-(U_ld_h1+1i*dm))*Y-I_ld_h1_im_rt;
  g_im_lt=(U_g-(U_ld_h1-1i*dm))*Y-I_ld_h1_im_lt;

  Jac=zeros(2,2);
  Jac(1,1)=real(g_re_rt-g_re_lt)/(2*dm);
  Jac(1,2)=real(g_im_rt-g_im_lt)/(2*dm);
  Jac(2,1)=imag(g_re_rt-g_re_lt)/(2*dm);
  Jac(2,2)=imag(g_im_rt-g_im_lt)/(2*dm);
end
```

Finally, Matlab solver with 5 iterations and initial guess $U_{ld-h1}=j$ looks like this:

```
t=[0:dt:0.02]; ##only one period!
U_g=1i;
Y=1/(1i*w*L);
U_ld_h1=1i; ##initial guess
for iter=1:5
  Jac=NR_Jacobian(U_ld_h1,w,t,dt,U_g,Y);
  I_ld_h1=nonlinear_load_current_main_harmonic(U_ld_h1,w,t,dt);
  g=(U_g-U_ld_h1)*Y-I_ld_h1;
  delta_u=-Jac\[real(g);imag(g)];
  U_ld_h1=U_ld_h1+delta_u(1)+1i*delta_u(2);
endfor
```
It gives $|U_{ld-h1}|=0.89302$ that is very close to the reference $0.89502$. 

For the second step to find load current 3<sup>rd</sup> and 5<sup>th</sup> harmonics $I_{ld-h3}$ and $I_{ld-h5}$ using $U_{ld-h1}$, we use the same approach as above and go from frequency to time domain. Matlab script for this is:

```
N=length(t);
u_ld=sqrt(2)*abs(U_ld_h1)*sin(w*t+angle(U_ld_h1));
lambda=zeros(1,N);
for ii=2:N
  lambda(ii)=lambda(ii-1)+0.5*(u_ld(ii-1)+u_ld(ii))*dt;
endfor
i_ld=sign(lambda).*(9.6539e+06*abs(lambda).^3+2.0879e+04*abs(lambda).^2);

I_ld_h3=mean(i_ld.*exp(-1i*2*pi/N*[1:N]*3))*sqrt(2)*exp(1i*pi/2);
I_ld_h5=mean(i_ld.*exp(-1i*2*pi/N*[1:N]*5))*sqrt(2)*exp(1i*pi/2);
```
Finally, for the last step, we find load voltage harmonics $U_{ld-h3}$ and $U_{ld-h5}$ using $I_{ld-h3}$ and $I_{ld-h5}$ as (Matlab script):

```
U_ld_h3=-I_ld_h3/(Y/3);
U_ld_h5=-I_ld_h5/(Y/5);
```

The figure below demonstrates the results: magnitudes of the harmonics are similar but there are noticeable errors for 3<sup>rd</sup> and 5<sup>th</sup>. 

![Harmonic power flow A](/images/harmonic-power-flow/harmonic_pf_A.jpg)  
*Figure 5: Decoupled approach: comparison of load voltage harmonic spectrum for time and frequency domain solutions.*

## Coupled approach

In contrast to the decoupled method, the coupled approach solves all harmonic components simultaneously. This allows capturing the interaction between harmonics, since the non-linear load current at each harmonic depends on the full voltage waveform rather than on a single harmonic component. Thus, we keep all harmonic voltages in the equations from the beginning. If we introduce the function $G$ for each harmonic $h$ as:

$$
G_h=(U_{g-h}-U_{ld-h})\frac{Y}{h}-I_{ld-h}
$$

with $U_{g-h}=0$ when $h>1$, then the size of the Jacobian increases significantly. Here only 1<sup>st</sup> and 3<sup>rd</sup> harmonics are included for shorter demonstration:

$$
\begin{bmatrix}
\frac{dG_{h1-re}}{dU_{ld-h1-re}} & \frac{dG_{h1-re}}{dU_{ld-h3-re}}  & \frac{dG_{h1-re}}{dU_{ld-h1-im}} & \frac{dG_{h1-re}}{dU_{ld-h3-im}}\\\\
\\\\
\frac{dG_{h3-re}}{dU_{ld-h1-re}} & \frac{dG_{h3-re}}{dU_{ld-h3-re}}  & \frac{dG_{h3-re}}{dU_{ld-h1-im}} & \frac{dG_{h3-re}}{dU_{ld-h3-im}}\\\\
\\\\
\frac{dG_{h1-im}}{dU_{ld-h1-re}} & \frac{dG_{h1-im}}{dU_{ld-h3-re}}  & \frac{dG_{h1-im}}{dU_{ld-h1-im}} & \frac{dG_{h1-im}}{dU_{ld-h3-im}}\\\\
\\\\
\frac{dG_{h3-im}}{dU_{ld-h1-re}} & \frac{dG_{h3-im}}{dU_{ld-h3-re}}  & \frac{dG_{h3-im}}{dU_{ld-h1-im}} & \frac{dG_{h3-im}}{dU_{ld-h3-im}}\\\\
\end{bmatrix} 
\cdot
\begin{bmatrix}
\Delta U_{ld-h1-re} \\\\
\Delta U_{ld-h3-re} \\\\
\Delta U_{ld-h1-im} \\\\
\Delta U_{ld-h3-im} \\\\
\end{bmatrix}
\=
\begin{bmatrix}
G_{h1-re} \\\\
G_{h3-re} \\\\
G_{h1-im} \\\\
G_{h3-im}
\end{bmatrix}
$$

To add 5<sup>th</sup> harmonic, we just extend the jacobian with two more rows and columns. Matlab function $G_h$ is:

```
function G_h=NR_knowns(U_ld_h,hs,w,t,dt,U_g,Y)
  I_ld_h=nonlinear_load_current_harmonics(U_ld_h,hs,w,t,dt);
  G_h=zeros(1,length(hs));
  U_g_h=[U_g 0 0];
  for hh=1:length(hs)
    G_h(hh)=(U_g_h(hh)-U_ld_h(hh))*Y/hs(hh)-I_ld_h(hh);
  endfor
end
```

Matlab Jacobian function is:

```
function Jac=NR_Jacobian(U_ld_h,hs,w,t,dt,U_g,Y)
  dm=0.001*eye(length(hs));
  dG_h_re=zeros(length(hs),length(hs));
  dG_h_im=zeros(length(hs),length(hs));
  for hh=1:length(hs)
    dG_h_re(hh,:)=(NR_knowns(U_ld_h+dm(hh,:),hs,w,t,dt,U_g,Y)-NR_knowns(U_ld_h-dm(hh,:),hs,w,t,dt,U_g,Y))/(2*dm);
    dG_h_im(hh,:)=(NR_knowns(U_ld_h+1i*dm(hh,:),hs,w,t,dt,U_g,Y)-NR_knowns(U_ld_h-1i*dm(hh,:),hs,w,t,dt,U_g,Y))/(2*dm);
  endfor
  Jac=zeros(length(hs)*2,length(hs)*2);
  Jac(1:length(hs),1:length(hs))=real(dG_h_re)';
  Jac(1:length(hs),length(hs)+1:end)=real(dG_h_im)';
  Jac(length(hs)+1:end,1:length(hs))=imag(dG_h_re)';
  Jac(length(hs)+1:end,length(hs)+1:end)=imag(dG_h_im)';
end
```

To find $I_{ld-h}=f(U_{ld-h})$, we need to take into account all harmonics to reconstruct voltage in time domain (DC offset is not present):

$$
u_{ld}=\sum_{h}|U_{ld-h}|\cdot\sin(\omega t\cdot h+\angle U_{ld-h})
$$

and Matlab function for this:

```
function I_ld_h=nonlinear_load_current_harmonics(U_ld_h,hs,w,t,dt)
  N=length(t);
  u_ld=0;
  for hh=1:length(hs)
    u_ld=u_ld+sqrt(2)*abs(U_ld_h(hh))*sin(w*t*hs(hh)+angle(U_ld_h(hh)));
  endfor
  lambda=zeros(1,N);
  for ii=2:N
    lambda(ii)=lambda(ii-1)+0.5*(u_ld(ii-1)+u_ld(ii))*dt;
  endfor
  i_ld=sign(lambda).*(9.6539e+06*abs(lambda).^3+2.0879e+04*abs(lambda).^2);
  I_ld_h=zeros(1,length(hs));
  for hh=1:length(hs)
    I_ld_h(hh)=mean(i_ld.*exp(-1i*2*pi/N*hs(hh)*[1:N]))*sqrt(2)*exp(1i*pi/2);
  endfor
end
```

Finally, Matlab solver with 5 iterations and initial guesses $U_{ld-h1}=j$, $U_{ld-h3}=0$ and $U_{ld-h5}=0$ looks like this:

```
U_ld_h=[1i 0 0]; ##initial guess
for iter=1:5
  Jac=NR_Jacobian(U_ld_h,hs,w,t,dt,U_g,Y);
  G_h=NR_knowns(U_ld_h,hs,w,t,dt,U_g,Y);
  delta_u=-Jac\[real(G_h)';imag(G_h)'];
  for hh=1:length(hs)
    U_ld_h(hh)=U_ld_h(hh)+delta_u(hh)+1i*delta_u(hh+length(hs));
  endfor
endfor
```

The outcome is all required harmonic voltages at once. The result is presented in the figure below and the magnitudes are very close to the reference values. 

![Harmonic power flow B](/images/harmonic-power-flow/harmonic_pf_B.jpg)  
*Figure 6: Coupled approach: comparison of load voltage harmonic spectrum for time and frequency domain solutions.*

## Conclusion

The results highlight the fundamental difference between the two approaches in terms of accuracy and computational cost. The table below shows a comparison of harmonic magnitudes (relative errors with respect to the reference values) obtained by the two described approaches. It shows that the coupled approach is much more precise, however it requires more  computation efforts. The decoupled approach is much simpler, but it produces noticeable errors. 

Order | Decoupled       |  Coupled
--- | ---       |   ---
1<sup>st</sup> | 0.22% | 0.00035%
3<sup>rd</sup> | 20.5% | 0.0016%
5<sup>th</sup> | 60.7% | 0.014%

It is worth mentioning that precision of the decoupled approach can be improved by doing several iterations in finding harmonic voltages. However, the final result depends on system complexity and might be still relatively poor.