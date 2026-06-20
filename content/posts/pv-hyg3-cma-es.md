+++
title = 'Finding model parameters of PV panel and hydro governor HYG3 using Covariance Matrix Adaptation Evolution Strategy'
date = 2026-06-19
draft = false
+++

## Introduction

This article investigates how to:
1. determine the five-parameter model of a PV panel using its data sheet,
2. find parameters of the hydro governor HYG3 from the known dynamic response curve.

utilizing the optimization algorithm CMA-ES (Covariance Matrix Adaptation Evolution Strategy) [CMA-ES home page](https://cma-es.github.io/). CMA-ES is a derivative-free optimization algorithm particularly well-suited for non-linear and non-convex problems. It adapts the covariance matrix of the search distribution, allowing efficient exploration of complex parameter spaces without requiring gradient information.

These cases also demonstrate how to model a PV panel and the dynamics of a complex system.  

## PV panel case

The five-parameter model of a PV panel is demonstrated in the figure below and it includes the following five parameters:
1. photo-generated current $I_{ph}$,
2. series resistance $R_s$,
3. shunt resistance $R_{sh}$,
4. thermal voltage $V_t$,
5. diode saturation current $I_0$. 

![PV model](/images/pv-hyg3-cma-es/PVmodel.jpg)  
*Figure 1: The five-parameter model of a PV panel.*

Here, the thermal voltage is calculated using the diode ideality factor $A$ and the number of cells $N_s$ as:

$$
V_t=N_sA\frac{kT_{STC}}{q}
$$

where $k$ is the Boltzmann constant, $q$ is the electron charge, and $T_{STC}=298.15$ K is the cell temperature at Standard Test Conditions. 

The model voltage-current characteristic $I_{pv}=f(V_{pv})$ is described by the following equation:

$$
I_{pv}=I_{ph}-I_0\exp\Big({\frac{V_{pv}+I_{pv}R_s}{V_t}}\Big)-\frac{V_{pv}+I_{pv}R_s}{R_{sh}}
$$

This non-linear equation is solved using a fixed-point iteration approach. It is used here to avoid directly solving the implicit non-linear equation for current. While simple to implement, the method may become numerically unstable near the open-circuit voltage, where the exponential term dominates and small numerical changes lead to large variations in the solution. The Matlab function looks like this: 

```
function Ipv=PVpanel5paramsModel(param,Vpv,Isc)
  Iph=param(1);
  Rs=param(2);
  Rsh=param(3);
  Vt=param(4)*1.3806E-23*(25+273.15)/1.6022E-19;
  I0=param(5)*1e-13;
  Ipv=zeros(1,length(Vpv));
  Ipv(1)=Isc;
  for ii=2:length(Vpv)
    Ipv(ii)=Iph-I0*exp((Vpv(ii)+Ipv(ii-1)*Rs)/Vt)-(Vpv(ii)+Ipv(ii-1)*Rs)/Rsh;
    Ipv(ii)=max([0 Ipv(ii)]);#avoid negative current
  endfor
end
```

Here, *Isc* is the known short circuit current. Note that the fourth parameter is $N_sA$ and the fifth is scaled by a factor $10^{-13}$. The introduced five parameters will be extracted from a PV panel data sheet. Below, the data for Siemens solar module SM55 is presented:

> MPP current 3.15 A ($I_{MPP}$)  
> MPP voltage 17.4 V ($V_{MPP}$)  
> Short circuit current 3.45 A ($I_{SC}$)  
> Open circuit voltage 21.7 V ($V_{OC}$)

MPP denotes a maximum power point. 

The CMA-ES algorithm requires a fitness function for its minimization and it will be based on comparison of three known voltage - current points (except $I_{SC}$) from the data sheet, and the actual values calculated using the model.  

$$
\max\Big(\Big|1-\frac{I_{MPP,calc}}{I_{MPP}}\Big|\text{,   }\Big|1-\frac{V_{MPP,calc}}{V_{MPP}}\Big|\text{,   }\Big|1-\frac{V_{OC,calc}}{V_{OC}}\Big|\text{,   }\min(I_{pv})\Big)
$$

The minimum PV current is added to make sure that the model achieves zero current at open circuit voltage. The Matlab implementation is:

```
function err=FitFunc(param,Impp,Vmpp,Isc,Voc)
  Vpv=linspace(0,Voc,200);
  Ipv=PVpanel5paramsModel(param,Vpv,Isc);
  Ppv=Vpv.*Ipv;
  [~,indx]=max(Ppv);
  Vmpp_actual=Vpv(indx);
  Impp_actual=Ipv(indx);
  Voc_actual=max(Vpv(Ipv>0));
  err=max([abs(1-[Impp_actual/Impp,Vmpp_actual/Vmpp,Voc_actual/Voc]) min(Ipv)]);
end
```

Here, the calculated MPPs are determined from the element-wise vector multiplication with the maximum search - $\max(I_{MPP}\cdot V_{MPP})$. Current and voltage are vectors with 200 samples. It was observed that too many samples can lead to numerical instability around open circuit voltage due to using fixed-point iteration approach as it was described above. Under-sampling has negative impact on fitness precision.

Finally, the CMA-ES algorithm activation code in Matlab is:

```
options = struct();
options.MaxIter = 1000;
options.SaveVariables = 'off';
options.LBounds = [0.9*Isc 1e-3 1 10 1]';
options.UBounds = [1.1*Isc 1 1e3 100 100]';
param=[Isc 0.1 100 10 10]';
tic;
xmin = cmaes('FitFunc', param, 0.3, options,Impp,Vmpp,Isc,Voc)
toc;
```

Function *cmaes* is a separate .m file that can be downloaded from the algorithm home page. The following options are defined: maximum number of iterations is 1000, variables are not saved, lower and upper boundaries for parameters (typical values). The initial guess is $I_{ph}=I_{SC}$, $R_s=0.1$ Ohm, $R_{sh}=100$ Ohm, $N_sA=10$ and $I_0=10$ A. Note that the algorithm performs much better if parameters are of the same order; therefore, $I_0$ is scaled by factor $10^{-13}$ inside function *PVpanel5paramsModel*. 

Running of the algorithm gives the following output:

```
  n=5: (4,8)-CMA-ES(w=[53 29 14 4]%, mu_eff=2.6) on function FitFunc
Iterat, #Fevals:   Function Value    (median,worst) |Axis Ratio|idx:Min SD idx:Max SD
    1 ,      9 : 6.3819095477387e-01 +(3e-02,6e-02) | 1.14e+00 |  1:3.4e-01  4:3.6e-01
    2 ,     17 : 6.2814070351759e-01 +(2e-02,3e-02) | 1.20e+00 |  1:3.4e-01  4:3.7e-01
   87 ,    697 : 5.0251256281407e-03 +(1e-04,5e-03) | 6.21e+00 |  1:5.9e-03  5:3.3e-02

#Fevals:   f(returned x)   |    bestever.f     | stopflag
    698: 5.02512562814e-03 | 3.55224397851e-03 | equalfunvals
mean solution: +3.5e+00 +6.7e-01 +9.9e+01 +2.8e+01 +3.4e+00
std deviation:  5.9e-03  1.3e-02  2.3e-02  2.0e-02  3.3e-02

xmin =
    3.4830
    0.6686
   99.0362
   28.3229
    3.3649

Elapsed time is 4.12997 seconds.
```

So, it finishes in about 4 seconds with 87 iterations and fitness function error about 0.5% with the following parameters: $I_{ph}=3.483$ A, $R_s=0.6686$ Ohm, $R_{sh}=99.0362$ Ohm, $N_sA=28.3229$ and $I_0=3.3649$ A. The PV model voltage-current and voltage-power characteristics with these parameters are demonstrated in the figure below:

![PV characteristics](/images/pv-hyg3-cma-es/PVcharacteristics.jpg)  
*Figure 2: PV panel voltage-current (left) and voltage-power (right) characteristics.*

It can be observed that the fitted parameters reproduce the key operating points from the data sheet with high accuracy. In particular, the maximum power point and open-circuit voltage closely match the reference values, validating the effectiveness of the optimization approach.

## Hydro governor HYG3 case

In contrast to the PV example, this case involves identifying parameters of a dynamic system. The hydro governor model includes multiple interconnected blocks with internal states, making the parameter estimation problem significantly more complex due to time dependency and feedback loops. The hydro governor HYG3 block-diagram is shown in the figure below:

![HYG3](/images/pv-hyg3-cma-es/HYG3.jpg)  
*Figure 3: The hydro governor HYG3 block-diagram [Power World page](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20HYG3.htm).*

Its dynamic response to a reference $P_{ref}$ step increase 0-1 pu and decrease 1-0.01 pu with stable frequency 1 pu (or $\Delta\omega=0$) is demonstrated in the figure below:

![HYG3 response](/images/pv-hyg3-cma-es/HYG3response.jpg)  
*Figure 4: Dynamic response of the HYG3 on input step of $P_{ref}$.*

It is known that the time step is 100 ms, $cflag>0$, $K_1=0$ (no derivative part), $R_{elec}=0$, $P_{aux}=0$, all dead-bands are zero. The parameters of the control part and the guide vanes part are required to be determined; in other words, all blocks up to signal $GV$. Parameters of the blocks representing water and turbine dynamics are known.

Firstly, we will establish the state-space model of the governor where each block is described as:

$$
\frac{dx}{dt}=Ax+Bu \\\\
y=Cx+Du
$$

where $x$ is the block state, $u$ - the input, $y$ - the output and $ABCD$ are parameters to be defined. So, the block output can be calculated from the known last state (if any) and the input. Thus, algebraic equations describing dynamics of the governor in the order of their execution are presented below (with reference to the diagram notations):

$$
y_\text{input filter}=x_1 \\\\
y_\text{pi integrator}=\min(\max(x_3\text{,   }P_{MIN})\text{,   }P_{MAX}) \\\\
y_\text{gate servo filter}=x_4 \\\\
GV=\min(\max(x_5\text{,   }P_{MIN})\text{,   }P_{MAX}) \\\\
q=x_6 \\\\
P_{GV}=N_{GV}(GV) \\\\
H=\Big(\frac{q}{P_{GV}}\Big)^2 \\\\
P_{mech}=(q-qNL)HA_t-GV\Delta\omega D_{turb} \\\\
CV=\min(\max(K_2y_\text{input filter}+y_\text{pi integrator}\text{,   }P_{MIN})\text{,   }P_{MAX}) \\\\
\frac{dx_1}{dt}=-\frac{1}{T_D}x_1+\frac{1}{T_D}(P_{ref}-R_{gate}CV) \\\\
\frac{dx_3}{dt}=K_Iy_\text{input filter} \\\\
\frac{dx_4}{dt}=\min(\max(-\frac{1}{T_P}x_4+\frac{K_G}{T_P}(CV-GV)\text{,   }VEL_{CLOSE})\text{,   }VEL_{OPEN}) \\\\
\frac{dx_5}{dt}=y_\text{gate servo filter} \\\\
\frac{dx_6}{dt}=\frac{1}{T_W}(H_0-H)
$$

The Matlab function based on these equations is:

```
function Pmech=HYG3dynamic(param,dt,Pref,delta_w)
  Rgate=param(1);
  TD=param(2);
  K2=param(3);
  KI=param(4);
  KG=param(5);
  TP=param(6);
  VEL_OPEN=param(7);
  VEL_CLOSE=param(8);

  #turbine and water parameters below are known
  N_GV=@(GV)max([0 min([0.9 0.9/0.7*GV])]);
  TW=1.5;
  H0=1;
  qNL=0;
  At=1;
  Dturb=0.5;
  P_MIN=0;
  P_MAX=1;

  Pmech=zeros(1,length(Pref));
  x=zeros(1,6);#second state is not in use
  x(5)=dt;#avoid division by 0
  x_der=zeros(1,length(x));
  for ii=1:length(Pref)
    y_input_filter=x(1);
    y_pi_integrator=min([max([x(3),P_MIN]),P_MAX]);
    y_gate_servo_filter=x(4);
    GV=min([max([x(5), P_MIN]), P_MAX]);
    q=x(6);
    P_GV=N_GV(GV);
    H=(q/P_GV)^2;
    Pmech(ii)=(q-qNL)*H*At-GV*delta_w(ii)*Dturb;
    CV=min([max([K2*y_input_filter+y_pi_integrator, P_MIN]), P_MAX]);
    x_der(1)=-1/TD*x(1)+1/TD*(Pref(ii)-Rgate*CV);
    x_der(3)=KI*y_input_filter;
    x_der(4)=min([max([-1/TP*x(4)+KG/TP*(CV-GV), VEL_CLOSE]), VEL_OPEN]);
    x_der(5)=y_gate_servo_filter;
    x_der(6)=1/TW*(H0-H);
    #update states using forward Euler method
    x=x+dt*x_der;
    #limit states after update
    x(3)=min([max([x(3),P_MIN]),P_MAX]);
    x(5)=min([max([x(5), P_MIN]), P_MAX]);
  endfor
end
```

This function has eight unknown parameters associated with the control system and guide vanes, $P_{MIN/MAX}$ are guessed based on the reference curve and given $N_{GV}$. For state update, the forward Euler method is used with the following limitation on the relevant states. 

Secondly, we construct the fitness function for minimization - it is based on Root Mean Square error between the reference curve and actual. The Matlab function implementing this is:

```
function err=FitFunc(param,dt,Pref,delta_w,Pmech_ref)
  Pmech=HYG3dynamic(param,dt,Pref,delta_w);
  err=sqrt(sum((Pmech_ref-Pmech).^2)/length(Pmech_ref));
end
```

Finally, the CMA-ES algorithm activation code in Matlab is:

```
options = struct();
options.MaxIter = 2000;
options.SaveVariables = 'off';
options.LBounds = [0.001 0.2 0 0 0 0.2 0 -0.1]';
options.UBounds = [0.1 10 10 10 10 10 0.1 0]';
param=[0.02 1 1 1 1 1 0.1 -0.1]';#initial guess
tic;
xmin = cmaes('FitFunc', param, 0.3, options,dt,Pref,delta_w,Pmech)
toc;
```

Here, maximum number of iterations is 2000. Note that due to 100 ms time step, time constants cannot be less than 200 ms to guarantee stable solution for the governor dynamics. The output of the optimization process looks like this:

```
  n=8: (5,10)-CMA-ES(w=[46 27 16 9 3]%, mu_eff=3.2) on function FitFunc
Iterat, #Fevals:   Function Value    (median,worst) |Axis Ratio|idx:Min SD idx:Max SD
    1 ,     11 : 2.0293511187010e-01 +(3e-01,3e-01) | 1.10e+00 |  2:4.6e-02  1:4.8e-02
    2 ,     21 : 1.6646381996029e-01 +(9e-02,3e-01) | 1.14e+00 |  3:4.2e-02  1:4.4e-02
  100 ,   1001 : 1.2608329722537e-01 +(4e-03,4e-02) | 1.52e+01 |  7:8.0e-04  2:9.1e-03
  200 ,   2001 : 7.0351723823643e-02 +(1e-02,2e-02) | 1.44e+02 |  7:7.6e-04  4:5.9e-02
  300 ,   3001 : 3.4777564701159e-02 +(4e-03,2e-02) | 2.30e+02 |  7:3.2e-04  3:4.4e-02
  400 ,   4001 : 2.4814078961975e-02 +(2e-03,5e-03) | 1.53e+03 |  7:1.4e-04  3:1.0e-01
  500 ,   5001 : 1.7052784820206e-02 +(1e-03,2e-02) | 2.52e+03 |  7:1.0e-04  3:1.5e-01
  600 ,   6001 : 1.5264976150185e-02 +(4e-05,2e-04) | 1.76e+03 |  7:1.8e-05  3:1.8e-02
  700 ,   7001 : 1.3580544205909e-02 +(3e-04,2e-03) | 4.10e+03 |  7:4.8e-05  3:1.3e-01
  800 ,   8001 : 1.0432600793715e-02 +(1e-04,4e-04) | 7.14e+03 |  7:2.9e-05  3:9.2e-02
  900 ,   9001 : 1.0105741043067e-02 +(4e-05,1e-04) | 7.91e+03 |  7:6.1e-06  3:2.4e-02
Iterat, #Fevals:   Function Value    (median,worst) |Axis Ratio|idx:Min SD idx:Max SD
 1000 ,  10001 : 9.7047740178444e-03 +(2e-04,1e-03) | 1.19e+04 |  7:2.8e-05  3:1.6e-01
 1100 ,  11001 : 9.4490323452818e-03 +(4e-06,1e-05) | 5.82e+03 |  7:2.6e-06  3:5.5e-03
 1200 ,  12001 : 8.6126481634275e-03 +(1e-04,3e-04) | 6.01e+03 |  7:1.5e-05  3:2.7e-02
 1300 ,  13001 : 4.3525779224783e-03 +(3e-04,2e-03) | 6.93e+03 |  7:1.8e-05  3:5.6e-02
 1400 ,  14001 : 1.7106519867031e-03 +(4e-05,2e-04) | 8.35e+03 |  7:5.0e-06  3:2.1e-02
 1500 ,  15001 : 1.3499189162254e-03 +(6e-05,1e-04) | 1.56e+04 |  7:7.6e-06  3:6.2e-02
 1600 ,  16001 : 4.0243451572931e-04 +(3e-05,1e-04) | 2.10e+04 |  7:1.7e-06  3:2.3e-02
 1700 ,  17001 : 9.9749166519724e-05 +(2e-05,1e-04) | 3.11e+04 |  7:1.8e-06  3:2.7e-02
 1800 ,  18001 : 3.2664359770340e-07 +(2e-07,6e-07) | 3.26e+04 |  7:1.0e-08  3:1.4e-04
 1900 ,  19001 : 1.1139800071692e-10 +(7e-11,6e-10) | 4.65e+04 |  7:3.3e-12  3:5.9e-08
 1997 ,  19971 : 7.6208640082190e-14 +(4e-14,6e-13) | 6.57e+04 |  7:1.7e-15  3:4.1e-11

#Fevals:   f(returned x)   |    bestever.f     | stopflag
  19972: 5.28872802764e-14 | 5.28872802764e-14 | tolfun
mean solution: +5.0e-02 +2.0e-01 +6.0e+00 +1.5e+00 +1.0e+00 +2.0e-01 +1.0e-02 -3.0e-02
std deviation:  4.2e-14  2.7e-12  4.1e-11  2.2e-12  3.9e-12  2.2e-11  1.7e-15  1.5e-14

xmin =

   5.0000e-02
   2.0134e-01
   6.0085e+00
   1.5004e+00
   1.0000e+00
   2.0000e-01
   1.0000e-02
  -3.0000e-02

Elapsed time is 1719.32 seconds.
```

So, it took almost half an hour and around 2000 iterations to achieve a fitness function value on the order of $10^{-14}$. The dynamic response with the found parameters is demonstrated in the figure below:

![Fitted response](/images/pv-hyg3-cma-es/FittedResponse.jpg)  
*Figure 5: Dynamic response of the HYG3 with the found parameters.*

and it is possible to observe that they are identical to the reference curve. The very low final error indicates that the model accurately reproduces the reference dynamic response. However, it should be noted that such near-perfect fitting may also indicate potential overfitting to the specific dataset, and validation against additional scenarios is recommended.

## Conclusion

The article demonstrates how to model a PV panel and the dynamics of a complex system, as well as how the modern optimization algorithm can be used to parametrize these models. The two selected use cases are common engineering tasks appearing in power system simulations and analysis. Additionally, the examples highlight the importance of proper parameter scaling and constraint selection for achieving robust and efficient optimization performance. One interesting outcome is the fact that optimization can be applied to tune the low-fidelity model to maximize its performance from the precision point of view. 
