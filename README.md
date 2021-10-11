# OpenVectorField
A general ode simulation tool

# Delpoyment
Create an empty folder, for example, YOUR_MATLAB_PATH\toolbox\ovf, and run file install_OVM.m.

# Quick start

  This library provides a structure-oriented ode simulation framework. The traditional simulation approach requires the variables of an equation to be a vector, which is inefficient and error-prone. To better reduce the mental burden, I have created a new simulation framework that uses structures to describe the variables of the vector field and encapsulates some of the associated methods.
 
Example:
 
## step1:model description
Suppose that we want to simulate a control system in 2D system:
x'=v,      (1)
v'=u,      (2)
to track a desired trajectory:
r = [cos(w*t),sin(w*t)].
Here, a simple solution to the control system is to design the control input u as
u = -k1 * (x-r) - k2 * (v-r') + r'',
where constants k1>0, k2>0.

## step2: variables and functions
W first define a structural variable as
```matlab
state0.x=unifrnd(-5,5,[2,1]);
state0.v=unifrnd(-5,5,[2,1]);
```
Note that it is not necessary to write r into the structure "state0",  because we can get r by defining a function using time t:
```matlab
para.w=0.2;
function [r,dr,ddr]=ref(t,para)
         w = para.w;
         r = [cos(w*t);sin(w*t)];
         dt= [-w*sin(w*t);w*cos(w*t)];
         ddr = [-w^2*cos(w*t);-w^2*sin(w*t)];
end
```
Next, we need to give a function to represent the dynamics (1) and (2). That is
```matlab
para.k1=1;
para.k2=1;
mySystem = vectorField(state0,@(t,s)myModel(t,s,para));
function dot = myModel(t,state,para)
    x=state.x;
    v=state.v;
    [r,dr,ddr]=ref(t,para);
    u = -para.k1*(x-r)-para.k2*(v-dr)+ddr;
    dot.x=v;
    dot.v=u;
end
```
Here, remeber to put down the functions "myModel" and "ref" at the bottom   of the .m file if "myModel" or "ref" is not in a single m file.

## step3:simulation
Thus, we can run the simulation by
```matlab
stime=30;
tspans=[0,stime];
state_end=mySystem.solve('ode45',tspans,state0);
```
Here, method "solve" returns the final result in the form of a structure,  which allows you to take state_end as an argument for the next time period of the simulation.

## step4:add signals
You may wonder why "solve" does not give you all the simulation results. This is because in most cases, we need to process the simulation data before visualizing it. In this case, what we need may be a tracking trajectory as well as the error curve of position errors. However, calculating the error using the raw simulation data is tedious. Here, we can use method "addSignals" to handle this issue.
Define a function as
```matlab
function signals = sig(t,state,para)
    [r,~] = ref(t,para);
    er = state.x - r;
    signals.r = r;
    signals.er = er;  
end
```
Next, we add the desired signal to the final result using the following method. 
```matlab
mySystem.addSignals(@(t,s)sig(t,s,para));
```
which means function "sig" will be executed at all time instants in the raw simulation data. Note that the structure returned by the equation "sig" will overwrite the structure of the same name in the final result, so you'd better make sure that it doesn't have the same name unless you know what you're doing. Then, the final result is collected by
```matlab
[t,data] = mySystem.result()
% or
% [t,data] = mySystem.result(2000)
```
where the parameter 2000 means that the final result consisted of 2000 pieces of data sampled from the raw data.
 
## step5:visualization
Finally, we can plot result as
```matlab
figure(1);clf;hold on
plot(data.x(:,1),data.x(:,2));
plot(data.r(:,1),data.r(:,2),'--');
plot(data.x(1,1),data.x(1,2),'x');
plot(data.r(1,1),data.r(1,2),'*');
plot(data.x(end,1),data.x(end,2),'o');
plot(data.r(end,1),data.r(end,2),'p');
title('tracking trajectory')
figure(2)
plot(t,data.er)
xlabel('t(s)')
title('position errors')
```
Full source code for the above example is in the example.m.
