+++
title = 'Induction motor direct torque control model in PSCAD'
date = 2026-05-22
draft = false
+++

## Introduction

This article briefly describes how to model induction motor direct torque control (DTC) in PSCAD.  
More about the control principles can be found here
[ABB drives technical guide](https://library.e.abb.com/public/3dc7259dcd684a55841165ddefe7eb3f/Technical_guide_No_1_3AFE58056685_RevD_EN.pdf?x-sign=CV811x+ZJnmbioTCCeJ4s9KDTQPJQxMEy55Ku2MVknKvQK2vJRH5bUL/VaPnXUvt)
and here 
[Matlab Direct Torque Control](https://www.mathworks.com/help/mcb/gs/direct-torque-control-dtc.html).  
The main focus is implementation of the inverter control (space vector modulator) in the most simple way.  
It is also of interest to analyze the simulation results using pure induction motor theory.

In this article you will learn:
- How to build a DTC-controlled induction motor model in PSCAD.
- How torque and flux are estimated.
- How the switching table is implemented.
- Typical simulation results and their interpretation.

## PSCAD model breakdown

The figure below shows the main power system part of the model setup consisting of:
 1. DC voltage source representing an AC supply system with a rectifier.
 2. Three phase IGBT inverter.
 3. Induction motor in torque control mode.
 4. Machine load (mechanical) torque is set as $T_m=\omega_m^2+0.3$, where $T_m$ and $\omega_m$ are mechanical torque and speed in per units.

![Motor and inverter part](/images/dtc-motor-control/inverter_motor.jpg)  
*Figure 1: Motor and inverter part of the PSCAD model.*

The next model part is dedicated to calculation (actually estimation) of the machine electrical torque $T_e$ and stator flux $\psi$ (peak) in per units according to these equations:
$$
\psi_d=\omega_n\int{(v_\alpha-i_\alpha R_s)dt}
$$
$$
\psi_q=\omega_n\int{(v_\beta-i_\beta R_s)dt}
$$
$$
T_e=0.5(\psi_d i_\beta - \psi_q i_\alpha)
$$

where $\psi_{d,q}$ is d,q-component of the flux in per units, $v_{\alpha,\beta}$ is $\alpha,\beta$-component of the measured motor voltage in per units, $i_{\alpha,\beta}$ is motor current $\alpha,\beta$-components in per units, $R_s$ is stator resistance in per units and $\omega_n=314.16$ is nominal radial frequency in rad/s.  

![Torque and flux calculations](/images/dtc-motor-control/torque_flux_calc.jpg)  
*Figure 2: Torque and flux calculation part of the PSCAD model.*

Note that input phase voltage in kV and phase current in kA are scaled to per units before transformation.  
Before voltage is integrated, $\gamma$-component is subtracted from $\alpha,\beta$-components to remove DC offset. 

The final part of the model is shown in figure below and it is used to control IGBT switches in the inverter. 

![Inverter control](/images/dtc-motor-control/inverter_control.jpg)  
*Figure 3: Inverter control part of the PSCAD model.*

Basically, it implements the following switching table:

 Sector 	|	  $\Delta\psi$ 	|	  $\Delta T_e$ 	|	 Voltage Vector 	|	 (Sa, Sb, Sc)      
---	|	---	|	---	|	---	|	---
1	|	1	|	1	|	 V2             	|	 (1,1,0)           
1	|	1	|	-1	|	 V6             	|	 (1,0,1)           
1	|	-1	|	1	|	 V3             	|	 (0,1,0)           
1	|	-1	|	-1	|	 V5             	|	 (0,0,1)           
1	|	1	|	0	|	 V0        	|	 (0,0,0)
1	|	-1	|	0	|	 V7        	|	 (1,1,1) 
2	|	1	|	1	|	 V3             	|	 (0,1,0)           
2	|	1	|	-1	|	 V1             	|	 (1,0,0)           
2	|	-1	|	1	|	 V4             	|	 (0,1,1)           
2	|	-1	|	-1	|	 V6             	|	 (1,0,1)           
2	|	1	|	0	|	 V0        	|	 (0,0,0)
2	|	-1	|	0	|	 V7        	|	 (1,1,1) 
3	|	1	|	1	|	 V4             	|	 (0,1,1)           
3	|	1	|	-1	|	 V2             	|	 (1,1,0)           
3	|	-1	|	1	|	 V5             	|	 (0,0,1)           
3	|	-1	|	-1	|	 V1             	|	 (1,0,0)           
3	|	1	|	0	|	 V0        	|	 (0,0,0)
3	|	-1	|	0	|	 V7        	|	 (1,1,1) 
4	|	1	|	1	|	 V5             	|	 (0,0,1)           
4	|	1	|	-1	|	 V3             	|	 (0,1,0)           
4	|	-1	|	1	|	 V6             	|	 (1,0,1)           
4	|	-1	|	-1	|	 V2             	|	 (1,1,0)           
4	|	1	|	0	|	 V0        	|	 (0,0,0)
4	|	-1	|	0	|	 V7        	|	 (1,1,1) 
5	|	1	|	1	|	 V6             	|	 (1,0,1)           
5	|	1	|	-1	|	 V4             	|	 (0,1,1)           
5	|	-1	|	1	|	 V1             	|	 (1,0,0)           
5	|	-1	|	-1	|	 V3             	|	 (0,1,0)           
5	|	1	|	0	|	 V0        	|	 (0,0,0)
5	|	-1	|	0	|	 V7        	|	 (1,1,1) 
6	|	1	|	1	|	 V1             	|	 (1,0,0)           
6	|	1	|	-1	|	 V5             	|	 (0,0,1)           
6	|	-1	|	1	|	 V2             	|	 (1,1,0)           
6	|	-1	|	-1	|	 V4             	|	 (0,1,1)           
6	|	1	|	0	|	 V0        	|	 (0,0,0)
6	|	-1	|	0	|	 V7        	|	 (1,1,1) 

here, depending on flux angle (which *Sector* it belongs to), flux error $\Delta\psi$ sign (reference minus calculated represented as 1-positive, -1-negative and 0-no error) and electrical torque error $\Delta T_e$ sign (the same representation), a specific voltage vector (*V0-V7*) is set using switches in phase a, b and c (e.g. *(Sa, Sb, Sc)=(1,1,0)* indicates states of upper switches with 1-open and 0-close, then lower have opposite states *(0,0,1)*).  
Each row gets a unique indicator calculated as $100 \cdot Sector + 10 \cdot \Delta\psi + \Delta T_e$. This indicator is used in the block *X-Y table* to determine voltage vector number (*V0-V7*). Mapping is stored in the .txt file loaded by the model during simulations. Finally, voltage vector number (*V0-V7*) is used in the next block *X-Y table* to determine switching states. The figure below show settings of this block: *x* is voltage vector number (*V0-V7*), *y1-y6* are output switching states. 

![Switching states](/images/dtc-motor-control/table2.jpg)  
*Figure 4: Table for mapping voltage vector number to switching states.*

Note, the model does not have dead-band for flux error that leads to noticeable flux ripples. Dead-band for torque error is set to 0.02. 

### Flux and torque references

Torque reference is directly set as a value in per units. Flux reference is discrete:
- it is kept 1 per units if mechanical speed is below 1 per units. Since $\psi \sim v/f$, then voltage will be equal to electrical speed in per units.
- when $\omega_m>1$, voltage is kept 1 per units, then flux reference is $1/\omega_m$ (field weakening). 

However, more precisely would be to use electrical speed. 

## Simulation results

The tests are run with torque reference 1 (0-10 seconds in the plot below), 0.8 (10-20 seconds) and 2 (20-30 seconds) per units. The latter value is of special interest to investigate what happens with the control when requested torque is much higher than available. The following is concluded:
 - Before 20 seconds: adequate and correct electrical torque response with respect to references 1 and 0.8 per units. Mechanical speed and RMS voltage are about the same that indicates flux 1 per units. Some moderate ripples are present. 
 - After 20 seconds: the control cannot provide the requested torque reference and actual steady-state maximum possible is settled (transient torque is high though). Voltage jumps to its maximum value limited by the DC voltage source. Considerable ripples are observed. 

![Simulation results](/images/dtc-motor-control/simData.jpg)  
*Figure 5: Simulated electrical torque, mechanical speed and RMS voltage.*

It is also interesting to look at distorted voltage and current waveforms presented in the figure below. 

![Simulation results electro](/images/dtc-motor-control/simData1.jpg)  
*Figure 6: Simulated instantaneous phase voltage and current.*

### Theoretical analysis of the results

Mechanical speed can be directly found from load torque expression $T_m=\omega_m^2+0.3$ assuming that it is equal to electrical torque reference. Thus, for references 1 and 0.8 per units, we get $\omega_m$ 0.84 and 0.7 respectively. Similar values are observed in Figure 5. For electrical torque reference 2, maximum should be found on torque-slip curve. To plot such curve, equivalent scheme of induction motor is used shown in the figure below.

<img src="/images/dtc-motor-control/motor_scheme.jpg" width="400">

*Figure 7: Induction motor equivalent scheme.*

Using Matlab function for calculation of electrical torque from the given motor parameters and mechanical speed:
```
function Te=tout(R1,L1,Lm,R2,L2,wm,s)
  Te=zeros(1,length(s));
  for ii=1:length(s)
    we=wm/(1-s(ii));
    V=min([1,we]);
    X1=we*L1;
    X2=we*L2;
    Xm=we*Lm;
    R2der=R2/s(ii);
    Rsub=R1*X2+R2der*X1;
    Xsub=X1*X2-R1*R2der;
    Rt=R1+R2der+Rsub/Xm;
    Xt=X1+X2+Xsub/Xm;
    Ir_abs_sqr=V^2/(Rt^2+Xt^2);
    Te(ii)=Ir_abs_sqr*(R2der-R2)/wm;
  endfor
end
```
we can now plot motor characteristics showing operation points for three torque references. 
```
R1=0.0119;
L1=0.1101;
Lm=3.6955;
R2=0.0307;
L2=0.1067;
s=linspace(0.0001,0.9999,1000);

Tref1=1;
wm=sqrt(Tref1-0.3);
Te1=tout(R1,L1,Lm,R2,L2,wm,s);

Tref2=0.8;
wm=sqrt(Tref2-0.3);
Te2=tout(R1,L1,Lm,R2,L2,wm,s);

Tref3=1.456;
wm=sqrt(Tref3-0.3);
Te3=tout(R1,L1,Lm,R2,L2,wm,s);
```
<img src="/images/dtc-motor-control/ts_curves.jpg" width="500">

*Figure 8: Motor torque-slip characteristics.*

For the case with electrical torque reference 2, actual maximum possible reference should be found such that operation point lands exactly on the curve maximum. As it is possible to observe, it is achieved with reference 1.456 (yellow curve). In such case, $\omega_m$ is 1.075. These numbers are similar what is observed in Figure 5. 


## Conclusion

- DTC can be implemented in PSCAD with simple logic blocks.
- It gives fast torque response but introduces ripples.
- The model matches well with theoretical motor characteristics.

## PSCAD files

[Download simulation example](/files/pscad-dtc-motor-control.zip)