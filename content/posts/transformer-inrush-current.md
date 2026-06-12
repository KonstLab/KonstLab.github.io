+++
title = 'Transformer sympathetic inrush current: analysis in frequency domain'
date = 2026-06-12
draft = false
+++

## Introduction

Transformer sympathetic inrush current occurs in a transformer that is already energized when a parallel transformer is energized. Usually, such a situation is studied using EMTP (electromagnetic transient program) or time domain simulations. However, it is of interest to develop an approach for analysis of this phenomenon in the frequency domain and investigate how it can be integrated into RMS-type simulators with a harmonic power flow (HPF) solver (where a strong second harmonic is typically present). 

This article describes such an approach in detail and investigates it for a simple test system. 

## Reference test system

The figure below demonstrates the test system in PSCAD that will be used for analysis. It has two parallel transformers (with enabled saturation of the magnetic core) loaded with resistors and a common grid with resistive internal impedance. Grid impedance is required to obtain voltage distortion and mutual interaction between transformers. All values are chosen such that simulated values will be equal to per-unit values for simplicity. 

![Test system](/images/transformer-inrush-current/test_system.jpg)  
*Figure 1: Test system in PSCAD.*

The figure below demonstrates DC offset, 1<sup>st</sup> and 2<sup>nd</sup> harmonic RMS, as well as true RMS current of both transformers. 

What is noticeable is an increase in current of the first transformer after energizing the second (simulation time after 0.51 seconds), whereas a decrease would normally be expected due to the voltage drop. It is also seen that the DC offset of the second transformer is negative and it is highly dependent on switching time (for example, there is no sympathetic effect at 0.5 seconds). Finally, as it can be observed, the DC offset (and the 2<sup>nd</sup> harmonic RMS) attenuate which is **the main challenge for representation in frequency domain analysis**. 

![Simulated data](/images/transformer-inrush-current/sim_data.jpg)  
*Figure 2: Time domain simulated data.*

## An approach for transformer inrush current representation in frequency domain

The key idea of the proposed method is to combine frequency-domain network solving with a time-domain representation of transformer magnetization. While the electrical network is solved using harmonic power flow, the nonlinear magnetizing branch is evaluated in time domain over one period. This allows accurate harmonic injection while maintaining compatibility with RMS-type simulation frameworks. To achieve this, transformer magnetization current - flux linkage function is required, i.e. $i_m=f(\lambda)$. In a simple way, it can be modeled as a two-piece linear characteristic shown in the figure below. 

![Transformer characteristics](/images/transformer-inrush-current/magn_current.jpg)  
*Figure 3: Transformer magnetization current - flux linkage characteristics (left) and time domain magnetization current with applied terminal voltage exceeding saturation limit (right).*

The first slope of the characteristic is defined by magnetization current available from a transformer data sheet (typically 1-5%). The second slope represents the saturation and it is typically approximately twice the transformer leakage inductance. Here, $\lambda=\int{u}$ and $\lambda_s$ is a saturation point. 

We will use the following Matlab function for calculations:

```
function im=NonLinearInductor(w,Lm,La,Vs,lambda)
  k1=1/Lm;
  k2=1/La;
  lambda_s=Vs*sqrt(2)/w;
  b=(k1-k2)*lambda_s;
  if abs(lambda)<lambda_s
    im=k1*lambda;
  else
    im=sign(lambda)*(k2*abs(lambda)+b);
  endif
end
```

When we know spectrum of transformer terminal voltage in frequency domain $U_{t-h}$ (where $h$ is harmonic number) and its DC offset $U_{t-DC}$, we can construct time domain curve as 

$$
u_{t}=U_{t-DC} + \sum_{h}\sqrt{2}\cdot|U_{t-h}|\cdot\sin(\omega t\cdot h+\angle U_{t-h})
$$

then integrate it numerically (e.g., using the trapezoidal rule) to get flux linkage $\lambda$ and finally get magnetization current in time domain $i_m$. Here $\omega$ is angular frequency and $t$ is a vector of time samples $[0, 1, 2, ..., N-1]\cdot\Delta t$ with fixed time step $\Delta t$ and number of samples in one period $N=2\pi/(\omega\Delta t)$.

The last step will be to extract harmonics $I_{m-h}$ using Discrete Fourier Transform (DFT):

$$
I_{m-h}=\frac{\sqrt{2}\exp(j\frac{\pi}{2})}{N}\sum_{n=0}^{N-1} i_m(n)\exp(-j2\pi\frac{n}{N}h)
$$

where $I_{m-h}$ is the $h$-order sin-based RMS harmonic phasor. DC offset $I_{m-DC}$ is just mean value of all samples. These quantities will be used further in HPF analysis. 

If, for example, $u_t$ is a pure sine function, then its time domain integration (on a sample-by-sample basis) for flux linkage will give a DC offset. This offset will be present in HPF solution, whereas in reality (as it is seen from Figure 2) it decays. This means that HPF analysis is performed over a time sequence updating transformer voltages $U_{t-h}$, $U_{t-DC}$ at each new time step. This allows introducing a decaying mechanism by adding resistive losses $R$ for transformer terminal voltages:

$$
U_{t-h/DC,\text{ }i}=U_{t-h/DC,\text{ }i-1}-RI_{m-h,\text{ }i-1}
$$

where $i$ and $i-1$ denote the current time step and the previous. The introduced resistive term should not be interpreted as a physical winding resistance. Instead, it acts as an equivalent damping mechanism representing losses in the magnetic core and surrounding network, enabling realistic decay of the DC component.

Note that $\lambda$ is also a vector of length $N$ and at each time step it depends on remanence, i.e. $\lambda=\int{u}+\lambda_{r}$. Thus, for each new step, $\lambda_{r}$ is the second sample of $\lambda$ - vector calculated at the previous step. It ensures continuity between successive time steps and it effectively propagates the magnetic state of the core, mimicking remanent flux evolution between HPF iterations. The Matlab implementation of the function calculating magnetization current harmonics is as follows:

```
function [Im_hs,lambda_r]=MagnetizingCurrentFD(Ut,hs,dt,w,Lm,La,Vs,Im_hs,lambda_r)
  Ut=Ut-0.005*Im_hs;#decaying
  Nw=2*pi/w/dt;
  t=[0:Nw-1]*dt;
  u=Ut(1);#first element is DC
  for hh=1:length(hs)
    u=u+sqrt(2)*abs(Ut(hh+1))*sin(w*t*hs(hh)+angle(Ut(hh+1)));
  endfor
  lambda=cumsum([lambda_r 0.5*(u(1:end-1)+u(2:end))*dt]);
  lambda_r=lambda(2);
  im=zeros(1,length(lambda));
  for ii=1:length(lambda)
    im(ii)=NonLinearInductor(w,Lm,La,Vs,lambda(ii));
  endfor
  Im_hs(1)=mean(im);#first element is DC
  for hh=1:length(hs)
    Im_hs(hh+1)=HarmonicDFT(im,hs(hh));
  endfor
end
```

Here, hardcoded value for resistive losses is $0.005$ (as will be shown later, it gives better precision with respect to the reference simulations). This function also uses DFT function that is defined as:

```
function X=HarmonicDFT(x,h)
  N=length(x);
  X=mean(x.*exp(-1i*2*pi/N*h*[0:N-1]))*sqrt(2)*exp(1i*pi/2);
end
```

### Selection of time step

Time step $\Delta t$ should be small enough to be able to get higher harmonics with sufficient precision. We are going to include the first 5 harmonics at fundamental frequency 50 Hz. According to the Nyquist - Shannon sampling theorem, sampling frequency should be higher than $2\cdot5\cdot50=500$ Hz or $\Delta t<2$ ms. For analysis, we will use 1 ms.

 ## Harmonic Power Flow solver

For HPF analysis, we are going to use the coupled approach described in detail here: [harmonic power flow analysis](https://konstlab.github.io/posts/harmonic-power-flow/). 

The nonlinear function $G$ when only the first transformer is energized can be written as: 

$$
G=\frac{U_g-U_t}{R_g}-\frac{U_t}{X_t\cdot h+R_{ld}}-I_{m1}
$$

Here all quantities are complex and $U_g$ is grid voltage with internal resistance $R_g$, $X_t$ is transformer leakage reactance (no winding losses included for simplicity), $R_{ld}$ is load resistance and $I_{m1}$ is magnetization current of the first transformer.

Note that this function represents a set of sub-functions defined for each harmonic separately as well as for DC offset (or $h=0$). Thus, $U_g=0$ for $h>1$ and $h=0$.   

When the second transformer is energized (with the same parameters and load), this function is changed as:

$$
G=\frac{U_g-U_t}{R_g}-2\cdot\frac{U_t}{X_t\cdot h+R_{ld}}-I_{m1}-I_{m2}
$$

where $I_{m2}$ is magnetization current of the second transformer.

The Matlab implementation looks like this:

```
function G=NR_knowns(Ut,hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged)
  G=zeros(1,length(hs)+1);
  if (transformer2Engaged)
    G(1)=(Ug*0-Ut(1))/Rg-2*Ut(1)/(Xt*0+Rld)-Im1(1)-Im2(1);#first element is DC
  else
    G(1)=(Ug*0-Ut(1))/Rg-Ut(1)/(Xt*0+Rld)-Im1(1);
  endif
  for hh=1:length(hs)
    if hh>1
      Ug=0;
    endif
    if (transformer2Engaged)
      G(hh+1)=(Ug-Ut(hh+1))/Rg-2*Ut(hh+1)/(Xt*hs(hh)+Rld)-Im1(hh+1)-Im2(hh+1);
    else
      G(hh+1)=(Ug-Ut(hh+1))/Rg-Ut(hh+1)/(Xt*hs(hh)+Rld)-Im1(hh+1);
    endif
  endfor
end
```

The Jacobian formation is fully described in the reference mentioned above: it is required to calculate derivatives of $G$ function over all harmonic voltages (split into real and imaginary parts) as well as DC offset. In this article, we demonstrate only the function as implemented in Matlab:

```
function Jac=NR_Jacobian(Ut,hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged)
  dm=0.001*eye(length(hs)+1);
  dG_h_re=zeros(length(hs)+1,length(hs)+1);
  dG_h_im=zeros(length(hs),length(hs)+1);
  for hh=1:length(hs)+1
    dG_h_re(hh,:)=(NR_knowns(Ut+dm(hh,:),hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged)-NR_knowns(Ut-dm(hh,:),hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged))/(2*dm);
    if hh>1
      dG_h_im(hh-1,:)=(NR_knowns(Ut+1i*dm(hh,:),hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged)-NR_knowns(Ut-1i*dm(hh,:),hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged))/(2*dm);
    end
  endfor
  Jac=zeros(length(hs)*2+1,length(hs)*2+1);
  Jac(1:length(hs)+1,1:length(hs)+1)=real(dG_h_re)';
  Jac(1:length(hs)+1,length(hs)+2:end)=real(dG_h_im)';
  Jac(length(hs)+2:end,1:length(hs)+1)=imag(dG_h_re(:,2:end))';
  Jac(length(hs)+2:end,length(hs)+2:end)=imag(dG_h_im(:,2:end))';
end
```

In the output matrix, the first column is derivatives for voltage DC offset, the second - real part of the first voltage harmonic, the third - real part of the second harmonic and so on, the last columns are for imaginary parts. The first row is for $G$ calculated for $h=0$, the second is the real part calculated for $h=1$ and so on for higher harmonics, the last rows are for imaginary parts.

The next Matlab function is HPF solver itself for one specific moment of time:

```
function [It1_out,It2_out,Ut,lambda1,lambda2,Im1,Im2]=CalculateTransformerCurrents(Niter,Ut,hs,dt,w,Lm,La,Vs,lambda1,lambda2,Im1,Im2,Ug,Rg,Xt,Rld,transformer2Engaged)
  for iter=1:Niter
    [Im1,lambda1_iter]=MagnetizingCurrentFD(Ut,hs,dt,w,Lm,La,Vs,Im1,lambda1);
    if transformer2Engaged
      [Im2,lambda2_iter]=MagnetizingCurrentFD(Ut,hs,dt,w,Lm,La,Vs,Im2,lambda2);
    else
      lambda2_iter=lambda2;
    end
    G=NR_knowns(Ut,hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged);
    Jac=NR_Jacobian(Ut,hs,Ug,Rg,Xt,Rld,Im1,Im2,transformer2Engaged);
    delta_u=-Jac\[real(G)';imag(G(2:end))'];
    Ut(1)=Ut(1)+delta_u(1);#first element is DC
    for hh=1:length(hs)
      Ut(hh+1)=Ut(hh+1)+delta_u(hh+1)+1i*delta_u(hh+1+length(hs));
    endfor
  end
  lambda1=lambda1_iter;
  lambda2=lambda2_iter;
  It1=zeros(1,length(hs)+1);
  It2=zeros(1,length(hs)+1);
  It1(1)=Ut(1)/(Xt*0+Rld)+Im1(1);#first element is DC
  for hh=1:length(hs)
    It1(hh+1)=Ut(hh+1)/(Xt*hs(hh)+Rld)+Im1(hh+1);
  endfor
  if transformer2Engaged
    It2(1)=Ut(1)/(Xt*0+Rld)+Im2(1);#first element is DC
    for hh=1:length(hs)
      It2(hh+1)=Ut(hh+1)/(Xt*hs(hh)+Rld)+Im2(hh+1);
    endfor
  end
  It1_out=[It1(1) abs(It1(2:3)) sqrt(sum(abs(It1).^2))]';
  It2_out=[It2(1) abs(It2(2:3)) sqrt(sum(abs(It2).^2))]';
end
```

Its outputs are transformer currents in frequency domain (namely DC offset, 1<sup>st</sup> and 2<sup>nd</sup> harmonic, as well as true RMS) and the following information to be passed to the next time step: terminal voltage in frequency domain, remanent flux linkage and magnetization current in frequency domain for both transformers. 

The core part of the function is Newton - Raphson solver with number of iterations *Niter*. It should be noted that remanent flux linkage is not updated inside these iterations, but only at HPF solver time step. 

Finally, the time running HPF solver with all input parameters is implemented in Matlab as follows:

```
Niter=5;#Newton-Raphson iterations
w=2*pi*50;
Xt=0.1i;
Lt=0.1/w;
Lm=100/w;
La=0.2/w;
Vs=1.17;
Rg=0.01;
Rld=1;
dt=0.001;
t=[0:dt:1];
lambda1=0;lambda2=0;Im1=0;Im2=0;#zero initial conditions
hs=[1:5];#harmonics
It1_out=zeros(4,length(t));
It2_out=zeros(4,length(t));
Ut=zeros(1,length(hs)+1)+0.01;Ut(2)=1;#initial guess
transformer2Engaged=false;
for ii=1:length(t)
  Ug=1*exp(1i*w*(ii-1)*dt);
  if t(ii)==0.51
    transformer2Engaged=true;
  endif
  if t(ii)>=0.1
    [It1_out(:,ii),It2_out(:,ii),Ut,lambda1,lambda2,Im1,Im2]=CalculateTransformerCurrents(Niter,Ut,hs,dt,w,Lm,La,Vs,lambda1,lambda2,Im1,Im2,Ug,Rg,Xt,Rld,transformer2Engaged);
  endif
endfor
```

So, as in the PSCAD tests, the first transformer is energized at 0.1 seconds and the second at 0.51 seconds. Only the first 5 harmonics are considered because it reduces size of the Jacobian which is $11\times11$. It is also worth mentioning that the grid voltage should change its phase (amplitude is fixed) over the time.

## Modeling results

The figure below demonstrates the modeling results as comparison of simulated transformer currents in time domain (PSCAD) and in frequency domain. As can be seen, they are quite similar and even coincide during the transient periods which indicates good accuracy of the described approach. Deviations at the moments of energizing are mainly explained by difference in approaches how RMS values are calculated in the time domain simulations: RMS values increase gradually, whereas in the frequency domain they are instantly determined. The time domain DC curves also have fast vanishing peaks that cannot be reproduced by the frequency domain simulations.

![Modeling results](/images/transformer-inrush-current/modeling_results.jpg)  
*Figure 4: Frequency domain modeling results.*

## Conclusion

The article demonstrates a practical approach for transformer inrush current calculations in frequency domain and compares it against time-domain simulations, demonstrating its validity and applicability in RMS - type simulators. The main challenge in the presented approach is assumed value for active losses introduced for modeling flux linkage decay. However, with a better approach to modeling transformer magnetization processes, this randomness can be excluded.