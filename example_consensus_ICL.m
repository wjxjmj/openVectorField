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
state=[];
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
sim.solve('ode15s',linspace(0,para.stime,para.ICL_Size+1),state); % solve system using selected solver
[t,result]=sim.result(para.resultSize); % collect simulation result.

% plot result
% all result are stacked in the structual variable result in the fashion of
% time*data matrices.
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
plot(t,linspace(para.a,para.a,length(t)),'--');
hold off

% model function
% The first two arguments must be the time and a structual state variable,
% and you can use any more arguments as you want after that.
% Remeber to return variable with the same structure as the state variable.

function grad = model(t,state,mem,para)

% unpack all states from structual variable 'state'
x   = state.x;
xd  = state.xd;
th  = state.th;
Y   = state.Y;
U   = state.U;

% this is the control input variable
u=zeros(size(x));

% define differential relationship
dot_x = zeros(size(x));
dot_xd=-xd;
dot_th=zeros(size(th));
dot_Y=zeros(size(Y));
dot_U=zeros(size(U));

% implemention of the algorithm
for i=1:para.n
    xi=x(:,i);
    thi=th(:,i);
    for j=1:para.n
        if i==j
            continue
        else
            xj=x(:,j);
            u(:,i)=u(:,i)+para.ka*(xj-xi); % interactions between neighboring agents
        end
    end
    u(:,i)=u(:,i)+para.kb*(xd-xi); % leader information feedback
    u(:,i)=u(:,i)-thi.*f(xi); % adaptive term
     
    dot_th(:,i)=(xi-xd)'*xi; % classic adaptive law
    for l=1:mem.n
        Yi = mem.data.Y{l}(:,i); % get the l-th Y from memory
        Ui = mem.data.U{l}(:,i); % get the l-th U from memory
        delta_xi = mem.data.delta_x{l}(:,i); % get x(t+delta_t)-x(t) from memory
        dot_th(:,i)=dot_th(:,i)+Yi'*(delta_xi-Ui-Yi*thi); % integral concurrent learning terms
    end
    dot_x(:,i)=u(:,i)+para.a*f(xi); % calculate the gradient of x
    dot_Y(:,i)=f(xi); % integral operation needed in the integral concurrenting learning adaptive control
    dot_U(:,i)=u(:,i); % same as to the previous line, note that what we need is the integral of the control input, not the gradient of the state x.
end

% set all gradients
grad.x  = dot_x;
grad.xd = dot_xd;
grad.th = dot_th;
grad.Y  = dot_Y;
grad.U  = dot_U;
end

function y=f(x)
y=sin(x);
end

% callback function for ICL.
% This function will be executed automatically at the end of each
% simulation time span.
% The control inputs arguments is defined to your need, but remeber to
% return a new structual state variable which is used as the initial
% value for the next simulation time span.
function new_state = recording(t,state,mem,para)
new_state = state; % we can directly copy from state.
index = mem.n+1; % indexing in matlab is started from 1, not 0
mem.data.delta_x{index}=state.x-mem.data.xt;
mem.data.Y{index}=state.Y;
mem.data.U{index}=state.U;
mem.data.xt=state.x;
mem.n=mem.n+1;% currently, we must maintain the counter manually.
new_state.U=new_state.U.*0; % reset the integral variable U
new_state.Y=new_state.Y.*0; % reset the integral variable Y
end
