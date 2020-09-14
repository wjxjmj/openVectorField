% this file is a part of the openVectorField project.
% Get the full toolbox at https://github.com/wjxjmj/openVectorField
clc
dim=2;
state0=[];
state0.x=unifrnd(-1,1,[dim,1]);
state0.v=unifrnd(-1,1,[dim,1]);
k=-1;
sys1 = vectorField(state0,@(t,state)myModel(t,state,k));

stime=30;
tspan=[0,stime];
[t,data]=sys1.ode45(tspan,stime,state0);

plot(data.x(:,1),data.x(:,2));
xlabel('time')
ylabel('x')

function dot = myModel(~,state,k)
x=state.x;
v=state.v;
dot.x=v; 
dot.v=-x-k*v; 
end

