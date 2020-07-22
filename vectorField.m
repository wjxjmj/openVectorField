% The class vectorFieldis designed to easy simulation of ode systems by
% developing a wrapper of the ode45 function.
% Auther: wjxjmj(github) ÐÇÃ¢(zhihu.com)
% E-mail: wjxjmj@126.com
% Initial uploading time: 2020-07-22
% To start with, the vectorField has two arguments including a structure
% that gives the initial values of all variablesand a function handle by
% which we describe the differential relationship of the variables.
% Example:
% Suppose we want to sumilate the following system in 2D system:
% ------x'=v                (*)
% ------v'=-x-k*v;      (**)
% with a random initial position x0 and a still initial velocity v0.
% Then we first define a structural variable as
% ------state.x=unifrnd(-5,5,[2,1]);
% ------state.v=zeros([2,1]);
% Next, we need to give a function to represent the dynamics (*,**). That
% is
% ------k=1;
% ------mySystem = vectorField(state0,@(t,state)myModel(t,state,k));
% ------function dot = myModel(~,state,k)
%         ©®---x=state.x;
%         ©®---v=state.v;
%         ©®---dot.x=v; 
%         ©®---dot.v=-x-k*v; 
% ------end
% Here, remeber to put down the function "myModel" to the tail of the m
% file if "myModel" is not seperated in a single m file.
% Thus, we can run the simulation by
% ------stime=30;
% ------tspan=[0,stime];
% ------[t,data]=mySystem.ode45(tspan,stime,state0);
% Finally, we can plot result as
% ------plot(data.x(:,1),data.x(:,2));
% ------xlabel('time')
% ------ylabel('x')
% Full source code for the above example is in the example.m. Just try it.

classdef vectorField < handle
   
    properties
        size_state='';
        index_a=[];
        index_b=[];
        size_vector=[];
        progress_bar=[];
        functionHandle=[];
        stime=[];
    end
    methods
        function self = vectorField(state,handle)
            if nargin==0
                self.clear();
            else
                iter=1;
                names=fieldnames(state);
                k = length(names);
                for i=1:k
                    var_i = state.(names{i});
                    self.size_state.(names{i})=size(var_i);
                    l = numel(var_i);
                    self.index_a.(names{i})=iter;
                    self.index_b.(names{i})=iter+l-1;
                    iter=iter+l;
                end
                self.size_vector=iter-1;
                self.functionHandle=@(t,vector)self.model(t,vector,handle);
            end
        end
        
        function [t,data,last]=ode45(self,tspan,state0,stime)
            if nargin==3
                stime=tspan(2);
            end
            self.stime=stime;
            self.initWaitBar('');
            [t,stream]=ode45(self.functionHandle,tspan,self.toVector(state0));
            data=self.fromVectorStream(stream);
            last=stream(end,:);
            self.closeWaitBar();
        end
        
        function dot_vector=model(self,t,vector,fun)
            state = self.fromVector(vector);
            grad  = fun(t,state);
            dot_vector=self.toVector(grad);
            self.updateWaitBar(t/self.stime);
        end
        
        function self = initWaitBar(self,str)
            self.progress_bar.handle=waitbar(0,str);
            self.progress_bar.initTime=tic;
        end
        
        function updateWaitBar(self,per)
            nt=toc(self.progress_bar.initTime);
            totalTime=nt/per*(1-per);
            timeStr=self.sec2time(totalTime);
            waitbar(per,self.progress_bar.handle,[timeStr,' left']);
        end
        
        function closeWaitBar(self)
            close(self.progress_bar.handle);
        end
        
        function vector = toVector(self,state)
            vector = zeros(self.size_vector,1);
            names = fieldnames(self.size_state);
            fieldname = inputname(2);
            k = length(names);
            for i=1:k
                assert(isfield(state,names{i}),['There is no variable named "',names{i},'" in "',fieldname,'", please check it.']);
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                s1 = self.size_state.(names{i});
                s2 = size(state.(names{i}));
                assert(isequal(s1,s2),['You are trying to change the size of "',names{i},'", please check it.'...
                                      ,' The size of "',names{i},'" should be ',mat2str(s1),', but the size of "'...
                                      ,fieldname,'.',names{i},'" is ',mat2str(s2)]);
                vector(a:b)=state.(names{i})(:);
            end
        end
 
        function state = fromVector(self,vector)
            state=[];
            vectorname = inputname(2);
            s = length(vector);
            assert(s==self.size_vector,['The size of input vector should to be ',mat2str(self.size_vector),', but the length of "'...
                                       ,vectorname,'" is ',num2str(s),', please check it']);
            names = fieldnames(self.size_state);
            n = length(names);
            for i=1:n
                size_i = self.size_state.(names{i});
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                state.(names{i})=reshape(vector(a:b),size_i);
            end
        end
        function result = fromVectorStream(self,vectorStream)
            result=[];
            streamname = inputname(2);
            s = size(vectorStream,2);
            assert(s==self.size_vector,['The column size of input vector stream should to be ',mat2str(self.size_vector),', but the column size of "'...
                                       ,streamname,'" is ',num2str(s),', please check it.']);
            names = fieldnames(self.size_state);
            n = length(names);
            for i=1:n
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                result.(names{i})=vectorStream(:,a:b);
            end
        end
        
        function self=clear(self)
        self.size_state='';
        self.index_a=[];
        self.index_b=[];
        self.size_vector=[];
        end
    end
    
    methods(Static)
        function str=sec2time(sec)
            sec=floor(sec);
            if sec<60
                str=sprintf('%d second(s)',sec);
            elseif sec<60^2
                a=mod(sec,60);
                b=(sec-a)/60;
                str=sprintf('%d minute(s)',b);
            elseif sec<60^2*24
                a=mod(sec,60);
                b=mod((sec-a)/60,60);
                c=(sec-b*60-a)/60/60;
                str=sprintf('%d hour(s)',c);
            elseif sec<60^2*24*30
                d=floor(sec/60/60/24);
                str=sprintf('%d day(s)',d);
            elseif sec<60^2*24*365
                d=floor(sec/60/60/24/30);
                str=sprintf('%d month(s)',d);
            else
                d=floor(sec/60/60/24/365);
                str=sprintf('%d year(s)',d);
            end
        end
    end

end

