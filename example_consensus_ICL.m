% This file is part of the opensource simulation platform 'openVectorField'
% Author: wjxjmj(github)

% system model:
% dot_x_i = -a*sin(x_i)+d_i
% leader dynamic:
% dot_xd = -xd

% here we define the structual parameter variable to contain all parameters
para.n=3; % number of agents
para.dim=2; % dimension of space
para.ka=1; % gain of consensus control
para.kb=1; % gain of leader tracking control
para.ICL_Size=10; % define how many data to collect in ICL algorithm
para.resultSize=500; % define how many data need to collect in result
para.stime=60; % define the simulation length
para.a=1; % unknown parameters in the system model


% here we define the structual state variable to contain all states
state.x     =   unifrnd(-3,3,[para.dim,para.n]);
state.th    =   zeros([1,para.n]);
state.xd    =   unifrnd(-3,3,[para.dim,1]);
state.Y     =   zeros([para.dim,para.n]);
state.U     =   zeros([para.dim,para.n]);

% here we define the structual memory variable to hold memories
data=[];
data.xt=state.x;
data.delta_x='';
data.Y='';
data.U='';

% create a memory instance
mem = Memory(data);

% create a simulation instance
sim = vectorField(state,@(t,state)model(t,state,mem,para)); % create instance using class 'vectorField'
sim.setFunctionAfterEachTspan(@(t,state)recording(t,state,mem,para)); % set callback function at the end of each simulation time span
[t,result]=sim.solve('ode15s',linspace(0,para.stime,para.ICL_Size+1),state,para.resultSize); % solve system using selected solver

% plot result
figure(1)
plot(result.xd(:,1),result.xd(:,2),'r--');
hold on
plot(result.xd(end,1),result.xd(end,2),'rp');
plot(result.x(:,1:para.dim:end),result.x(:,2:para.dim:end),'b-');
plot(result.x(end,1:para.dim:end),result.x(end,2:para.dim:end),'bo');
hold off

figure(2)
plot(t,result.th);
hold on
plot(t,linspace(para.k,para.k,length(t)));
hold off


function grad = model(t,state,mem,para)
x   = state.x;
xd  = state.xd;
th  = state.th;
Y   = state.Y;
U   = state.U;

u=zeros(size(x));
dot_x = zeros(size(x));
dot_xd=-xd;
dot_th=zeros(size(th));
dot_Y=zeros(size(Y));
dot_U=zeros(size(U));

for i=1:para.n
    xi=x(:,i);
    thi=th(:,i);
    for j=1:para.n
        if i==j
            continue
        else
            xj=x(:,j);
            u(:,i)=u(:,i)+para.ka*(xj-xi);
        end
    end
    u(:,i)=u(:,i)-thi.*f(xi)+para.kb*(xd-xi);
    dot_th(:,i)=(xi-xd)'*xi;
    for l=1:mem.n
        Yi = mem.data.Y{l}(:,i);
        Ui = mem.data.U{l}(:,i);
        delta_xi = mem.data.delta_x{l}(:,i);
        dot_th(:,i)=dot_th(:,i)+Yi'*(delta_xi-Ui-Yi*thi);
    end
    dot_x(:,i)=u(:,i)+para.a*f(xi);
    dot_Y(:,i)=f(xi);
    dot_U(:,i)=u(:,i);
end

grad.x  = dot_x;
grad.xd = dot_xd;
grad.th = dot_th;
grad.Y  = dot_Y;
grad.U  = dot_U;
end

function y=f(x)
y=sin(x);
end

function new_state = recording(t,state,mem,para)
new_state = state;
index = mem.n+1;
mem.data.delta_x{index}=state.x-mem.data.xt;
mem.data.Y{index}=state.Y;
mem.data.U{index}=state.U;
mem.data.xt=state.x;
mem.n=mem.n+1;
new_state.U=new_state.U.*0;
new_state.Y=new_state.Y.*0;
end
