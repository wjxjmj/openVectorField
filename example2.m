dim=2;
n=10;
th_size=3;
state0=[];
state0.q    =[unifrnd(-15,0,[1,n]);unifrnd(-7.5,7.5,[1,n])];
state0.dot_q=unifrnd(-1,1,[dim,n]);
state0.theta=zeros(th_size,n);
state0.s    =zeros(dim,n);
state0.d    =unifrnd(-1,1,[dim,n]);
state0.angle=pi/2;

a=unifrnd(1,2,[n,1]);
b=unifrnd(1,2,[n,1]);
c=unifrnd(1,2,[n,1]);
ks=unifrnd(2,5,[n,1]);
h=unifrnd(-1,1,[d,n]);
l=1;
theta_star=[1;2;3];
rd=10;

sys = vectorField(state0,@(t,state)consensus(t,state,dim,n,a,b,c,ks,h,l,rd,theta_star));

tspan=[0,50];
[t,data]=sys.ode45(tspan,state0);
fig=figure(1);iter=1;
skip=floor(length(t)/200);
r=[rd*2^0.5.*cos(data.angle)./(1+sin(data.angle).^2).*1,...
   rd*2^0.5.*cos(data.angle)./(1+sin(data.angle).^2).*sin(data.angle)];
for i=1:skip:length(t)
    plot(data.q(1,1:2:dim*n),data.q(1,2:2:dim*n),'bx');
    hold on
    plot(data.q(1:i,1:2:dim*n),data.q(1:i,2:2:dim*n));
    plot(r(1:i,1),r(1:i,2),'r--')
    plot(r(i,1),r(i,2),'rp')
    axis equal
    hold off
    title([num2str(round(i/length(t)*100)),'%'])
    drawnow
%     frame = getframe(fig);im{iter}=frame2im(frame);iter=iter+1;
end
title('100%')

% frame = getframe(fig); 
% im{iter}=frame2im(frame);
% close
% filename = 'infty.gif'; % Specify the output file name
% for i = 1:iter
%     [A,map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
%     end
% end





function grad = consensus(t,state,dim,n,a,b,c,ks,h,l,rd,theta_star)

% trajectory
dot_angle = 0.038*pi;
r =     rd*2^0.5*cos(state.angle)/(1+sin(state.angle)^2).*[1;sin(state.angle)];
dot_r = [-((2*sqrt(2)*rd*sin(state.angle)*(cos(2*state.angle) + 5))/(cos(2*state.angle) - 3)^2);
        (rd*(3*cos(2*state.angle) - 1))/(2^0.5*(sin(state.angle)^2 + 1)^2)].*dot_angle;

% define tmp variables
ddot_q=zeros(dim,n);
u=zeros(dim,n);
us=zeros(dim,n);
ud=zeros(dim,n);
u_theta=zeros(size(state.theta));

%calculate control inputs
for i=1:n
    qi = state.q(:,i);
    vi = state.dot_q(:,i);
    ei = a(i)*(qi-r)+(vi-dot_r);
    u(:,i)=u(:,i)-f(qi,vi)*state.theta(:,i)-a(i)*(vi-dot_r)-(1+ks(i))*ei+state.s(:,i);
    u_theta(:,i)=ei'*f(qi,vi);
    us(:,i)=-(1+ks(i))*ei-b(i)*sign(ei);
    H=[0,-h(i);h(i),0];
    ud(:,i)=H*state.d(:,i);
    for j=1:n
        if i==j
            continue
        else
            qj = state.q(:,j);
            vj = state.dot_q(:,j);
            qij = qi-qj;
            dij2 =  qij'*qij;
            if dij2<l^2
                u(:,i)=u(:,i)-c(i)*(vi-vj);
            end
        end
    end
    ddot_q(:,i)=f(qi,vi)*theta_star+state.d(:,i)+u(:,i);
end

% assign values for return structural variable
grad.q=state.dot_q;
grad.dot_q=ddot_q;
grad.theta=u_theta;
grad.s=us;
grad.d=ud;
grad.angle=dot_angle;

end

function y=f(q,dot_q)
y=[ cos(q(1)),  -dot_q(2), -sin(q(2));
    sin(q(2)),  -dot_q(1),  cos(q(1))];
end