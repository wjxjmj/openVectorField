para.w=0.2;
stime=30;
tspans=[0,stime];
para.k1=1;
para.k2=1;
state0.x=unifrnd(-5,5,[2,1]);
state0.v=unifrnd(-5,5,[2,1]);
mySystem = vectorField(state0,@(t,s)myModel(t,s,para));
state_end=mySystem.solve('ode45',tspans,state0);
mySystem.addSignals(@(t,s)sig(t,s,para));
[t,data] = mySystem.result();

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

function dot = myModel(t,state,para)
         x=state.x;
         v=state.v;
         [r,dr,ddr]=ref(t,para);
         u = -para.k1*(x-r)-para.k2*(v-dr)+ddr;
         dot.x=v;
         dot.v=u;
end

function [r,dr,ddr]=ref(t,para)
         w = para.w;
         r = [cos(w*t);sin(w*t)];
         dr= [-w*sin(w*t);w*cos(w*t)];
         ddr = [-w^2*cos(w*t);-w^2*sin(w*t)];
end
function signals = sig(t,state,para)
         [r,~] = ref(t,para);
         er = state.x - r;
         signals.r = r;
         signals.er = er;  
end