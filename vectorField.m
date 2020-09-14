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
        name_state='';
        size_state='';
        index_a=[];
        index_b=[];
        size_vector=[];
        progress_bar=[];
        functionHandle=[];
        functionName='';
        functionAfterEachTspanHandle=[];
        functionAfterEachTspanName='';
        inputName=[];
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
                    self.name_state{i}=names{i};
                    l = numel(var_i);
                    self.index_a.(names{i})=iter;
                    self.index_b.(names{i})=iter+l-1;
                    iter=iter+l;
                end
                self.size_vector=iter-1;
                self.functionName= self.getFunctionName(handle);
                self.functionHandle=@(t,vector)self.model(t,vector,handle);
            end
        end
        
        function strs = translate(self,fileName,solverName,functionName,tspanName,structureName,vectorName,retValName,timeName,resultName)
            idt=[];
            iter=1;
            line1=[vectorName,'=['];
            for i=1:length(self.name_state)
                line1=[line1,structureName,'.',self.name_state{i},'(:);'];
            end
            line1=[line1,'];'];
            strs{iter}=line1;iter=iter+1;idt=[idt,0];
            strs{iter}=['[',timeName,',',resultName,']=',solverName,'(@(t,vec)',functionName,'(t,vec),',tspanName,',',vectorName,');'];iter=iter+1;idt=[idt,0];
            for i=1:length(self.name_state)
                a = self.index_a.(self.name_state{i});
                b = self.index_b.(self.name_state{i});
                strs{iter}=[resultName,'_',self.name_state{i},'=',resultName,'(:,',mat2str(a),':',mat2str(b),');'];iter=iter+1;idt=[idt,0];
            end
            strs{iter}='%%visualization...';iter=iter+1;idt=[idt,0];
            strs{iter}='';iter=iter+1;idt=[idt,0];
            strs{iter}=['function ',retValName,'=',functionName,'(t,',vectorName,')'];iter=iter+1;idt=[idt,0];
            strs{iter}='%%reverse vectorlization...';iter=iter+1;idt=[idt,1];
            for i=1:length(self.name_state)
                a = self.index_a.(self.name_state{i});
                b = self.index_b.(self.name_state{i});
                size_i=self.size_state.(self.name_state{i});
                strs{iter}=[self.name_state{i},'=','reshape(',vectorName,'(',mat2str(a),':',mat2str(b),'),',mat2str(size_i),');'];iter=iter+1;idt=[idt,1];
            end
            strs{iter}='';iter=iter+1;idt=[idt,1];
            strs{iter}='%%define derivative variables...';iter=iter+1;idt=[idt,1];
            for i=1:length(self.name_state)
                size_i=self.size_state.(self.name_state{i});
                strs{iter}=['dot_',self.name_state{i},'=zeros(',mat2str(size_i),');'];iter=iter+1;idt=[idt,1];
            end
            strs{iter}='';iter=iter+1;idt=[idt,1];
            strs{iter}='%%write your algorithm code here...';iter=iter+1;idt=[idt,1];
            strs{iter}='';iter=iter+1;idt=[idt,1];
            strs{iter}='';iter=iter+1;idt=[idt,1];
            strs{iter}='';iter=iter+1;idt=[idt,1];
            strs{iter}='%%vectorlization...';iter=iter+1;idt=[idt,1];
            strs{iter}=[retValName,'=zeros([',mat2str(self.index_b.(self.name_state{end})),',1]);'];iter=iter+1;idt=[idt,1];
            for j=1:length(self.name_state)
                a = self.index_a.(self.name_state{j});
                b = self.index_b.(self.name_state{j});
                strs{iter}=[retValName,'(',mat2str(a),':',mat2str(b),')=dot_',self.name_state{j},'(:);'];iter=iter+1;idt=[idt,1];
            end
            strs{iter}='end';iter=iter+1;idt=[idt,0];
            fid = fopen(fileName,'w');
            for i=1:length(strs)
                for j=1:idt(i)
                    fprintf('    ');
                    fprintf(fid,'%s','    ');
                end
                fprintf(strs{i});
                fprintf(fid,'%s',strs{i});
                fprintf('\n');
                fprintf(fid,'\n');
            end
            fclose(fid);
            
%             for i=1:size(self.name)
%                 strs{i}=[inputName,'=',inputName,'()'];
%                 strs{1}=[name_state'hello';
%                 strs{2}='world';
%             end
        end
        

        function [t,data,last]=solve(self,solver,tspans,state0,N)
            if nargin<=4
                N=-1;
            end
            self.stime=tspans(end);
            self.inputName=inputname(3);
            self.initWaitBar('');
            vector = self.fromStateToVector(state0);
            t_size=length(tspans);
            stream=[];
            t=[];
            for i=1:t_size-1
                tspan_i=tspans(i:i+1);
                switch lower(solver)
                    case {'ode113'}
                        [t_i,stream_i]=ode113(self.functionHandle,tspan_i,vector);
                    case {'ode15s'}
                        [t_i,stream_i]=ode15s(self.functionHandle,tspan_i,vector);
                    case {'rk4'}
                        [t_i,stream_i]=self.rk4(self.functionHandle,tspan_i,vector,0.001);
                    otherwise
                        [t_i,stream_i]=ode45(self.functionHandle,tspan_i,vector);
                end
                vector=stream_i(end,:)';
                if ~isempty(self.functionAfterEachTspanName)
                    state_i=self.fromVector(vector);
                    state_tspan=self.functionAfterEachTspanHandle(t_i(end),state_i);
                    vector=self.toVector(state_tspan);
                end
                t=vertcat(t,t_i);
                stream=vertcat(stream,stream_i);
            end
            if N>0 && N<length(t)
                indexing = self.sampling(t,N);
                t= t(indexing);
                stream=stream(indexing,:);
            end
            data=self.fromVectorStream(stream);
            last=stream(end,:);
            self.closeWaitBar();
        end
        
        function [t,data,last]=ode45c(self,tspans,state0,N)
            if nargin<4
                N=-1;
            end
            self.stime=tspans(end);
            self.inputName=inputname(3);
            self.initWaitBar('');
            vector = self.fromStateToVector(state0);
            t_size=length(tspans);
            stream=[];
            t=[];
            for i=1:t_size-1
                tspan_i=tspans(i:i+1);
                [t_i,stream_i]=ode45(self.functionHandle,tspan_i,vector);
                vector=stream_i(end,:)';
                if ~isempty(self.functionAfterEachTspanName)
                    state_i=self.fromVector(vector);
                    state_tspan=self.functionAfterEachTspanHandle(t_i(end),state_i);
                    vector=self.toVector(state_tspan);
                end
                t=vertcat(t,t_i);
                stream=vertcat(stream,stream_i);
            end
            if N>0 && N<length(t)
                indexing = self.sampling(t,N);
                t= t(indexing);
                stream=stream(indexing,:);
            end
            data=self.fromVectorStream(stream);
            last=stream(end,:);
            self.closeWaitBar();
        end
        
        %         function [t,data,last]=ode45c(self,tspans,state0,stime)
        %             self.stime=tspans(end);
        %             self.inputName=inputname(3);
        %             self.initWaitBar('');
        %             vector = self.fromStateToVector(state0);
        %             [t,stream]=ode45(self.functionHandle,tspans,vector);
        %             data=self.fromVectorStream(stream);
        %             last=stream(end,:);
        %             self.closeWaitBar();
        %         end
        
        function [t,data,last]=ode45(self,tspan,state0,stime)
            if nargin==3
                stime=tspan(2);
            end
            self.inputName=inputname(3);
            self.stime=stime;
            self.initWaitBar('');
            vector = self.fromStateToVector(state0);
            [t,stream]=ode45(self.functionHandle,tspan,vector);
            data=self.fromVectorStream(stream);
            last=stream(end,:);
            self.closeWaitBar();
        end
        %
        %         function [t,data,last]=ode45(self,tspan,state0,stime)
        %             if nargin==3
        %                 stime=tspan(2);
        %             end
        % %             inputname(3)
        %             self.stime=stime;
        %             self.initWaitBar('');
        %             [t,stream]=ode45(self.functionHandle,tspan,self.toVector(state0));
        %             data=self.fromVectorStream(stream);
        %             last=stream(end,:);
        %             self.closeWaitBar();
        %         end
        
        %         function [t,data,last]=ode45(self,varargin)
        %             tspan=varargin{1};
        %             state0=varargin{2};
        %             if nargin==3
        %                 stime=tspan(3);
        %             end
        %             self.stime=varargin{3};
        %             self.initWaitBar('');
        %             [t,stream]=ode45(self.functionHandle,tspan,self.toVector(varargin{2}));
        %             data=self.fromVectorStream(stream);
        %             last=stream(end,:);
        %             self.closeWaitBar();
        %         end
        function self=setFunctionAfterEachTspan(self,handle)
            self.functionAfterEachTspanName= self.getFunctionName(handle);
            self.functionAfterEachTspanHandle=handle;
        end
        
        function dot_vector=model(self,t,vector,fun)
            %             global debug_i
            %             if debug_i==2
            %                 debug_i;
            %             end
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
            waitbar(per,self.progress_bar.handle,[timeStr,' left',' for "',self.functionName,'"']);
        end
        
        function closeWaitBar(self)
            close(self.progress_bar.handle);
        end
        
        function vector = fromStateToVector(self,state)
            vector = zeros(self.size_vector,1);
            names = fieldnames(self.size_state);
            k = length(names);
            for i=1:k
                assert(isfield(state,names{i}),['No variable named "',names{i}...
                    ,'" found. Please check if "',self.inputName,'" contains variable "',names{i},'".']);
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                s1 = self.size_state.(names{i});
                s2 = size(state.(names{i}));
                assert(isequal(s1,s2),['You are trying to change the size of "',names{i},'", please check it.'...
                    ,' The size of "',names{i},'" should be ',mat2str(s1),', but the size of "'...
                    ,self.inputName,'.',names{i},'" is ',mat2str(s2),'.']);
                vector(a:b)=state.(names{i})(:);
            end
        end
        
        function vector = toVector(self,state)
            vector = zeros(self.size_vector,1);
            names = fieldnames(self.size_state);
            k = length(names);
            for i=1:k
                assert(isfield(state,names{i}),['No variable named "',names{i}...
                    ,'" found. Please check if the structural variable returned from function "'...
                    ,self.functionName,'" contains variable "',names{i},'".']);
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                s1 = self.size_state.(names{i});
                s2 = size(state.(names{i}));
                assert(isequal(s1,s2),['You are trying to change the size of "',names{i},'", please check it.'...
                    ,' The size of "',names{i},'" should be ',mat2str(s1),', but the size of "'...
                    ,names{i},'" returned from function "',self.functionName,'" is ',mat2str(s2),'.']);
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
        
        
        function indexing=sampling(t,N)
            expect = linspace(min(t),max(t),N);
            indexing=zeros(N,1);
            index=1;
            for i=1:length(expect)
                [~,index_i]=min(abs(t(index:end)-expect(i)));
                indexing(i)=index+index_i-1;
                index=index_i;
            end
        end
        
        function str = getFunctionName(handle)
            handle_name = func2str(handle);
            strSet=regexp(handle_name,'(?<=\)).+?(?=\()','match');
            if isempty(strSet)
                str=handle_name;
            else
                str=strSet{1};
            end
        end
        function [t,stream]=rk4(ufunc,tspan,state0,dt)
            a=tspan(1);
            b=tspan(2);
            loop=floor((b-a)/dt);
            t=zeros(loop+1,1);
            stream=zeros(loop+1,length(state0));
            t(1)=a;
            stream(1,:)=state0;
            for iter=1:loop
                t(iter+1)=t(iter)+dt;
                k1=ufunc(t(iter),stream(iter,:));
                k2=ufunc(t(iter)+dt/2,stream(iter,:)+dt*k1/2);
                k3=ufunc(t(iter)+dt/2,stream(iter,:)+dt*k2/2);
                k4=ufunc(t(iter)+dt,stream(iter,:)+dt*k3);
                stream(iter+1,:)=stream(iter,:)+(dt*(k1+2*k2+2*k3+k4)/6)';
            end
        end
    end
    
end

