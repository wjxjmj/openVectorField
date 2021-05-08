% The class vectorFieldis designed to easy simulation of ode systems by
% developing a wrapper of the ode45 function.
% Auther: wjxjmj(github) Mito.W(zhihu.com)
% E-mail: wjxjmj@126.com
% Initial uploading time: 2020-07-22
% To start with, the vectorField has two arguments including a structure
% that gives the initial values of all variables and a function handle by
% which we describe the differential relationship of the variables.
%
% Example:
%
% step1:model description
% Suppose that we want to simulate a control system in 2D system:
% ------x'=v,               (*)
% ------v'=u,      (**)
% to track a desired trajectory:
% ------r = [cos(w*t),sin(w*t)].
% Here, a simple solution to the control system is to design the control
% input u as
% ------u = -k1 * (x-r) - k2 * (v-r') + r'',
% where constants k1>0, k2>0.
%
% step2: variables and functions
% W first define a structural variable as
% ------state0.x=unifrnd(-5,5,[2,1]);
% ------state0.v=unifrnd(-5,5,[2,1]);
% Note that it is not necessary % to write r into the structure "state0",
% because we can get r by defining a function using time
% ------para.w=0.2;
% ------function [r,dr,ddr]=ref(t,para)
%          |---w = para.w;
%          |---r = [cos(w*t);sin(w*t)];
%          |-- dt= [-w*sin(w*t);w*cos(w*t)];
%          |---ddr = [-w^2*cos(w*t);-w^2*sin(w*t)];
% ------end
% Next, we need to give a function to represent the dynamics (*,**). That
% is
% ------para.k1=1;
% ------para.k2=1;
% ------mySystem = vectorField(state0,@(t,s)myModel(t,s,para));
% ------function dot = myModel(t,state,para)
%          |---x=state.x;
%          |---v=state.v;
%          |---[r,dr,ddr]=ref(t,para);
%          |---u = -para.k1*(x-r)-para.k2*(v-dr)+ddr;
%          |---dot.x=v;
%          |---dot.v=u;
% ------end
% Here, remeber to put down the functions "myModel" and "ref" at the bottom
% of the .m file if "myModel" or "ref" is not in a single m file.
%
% step3:simulation
% Thus, we can run the simulation by
% ------stime=30;
% ------tspans=[0,stime];
% ------state_end=mySystem.solve('ode45',tspans,state0);
% Here, method "solve" returns the final result in the form of a structure,
% which allows you to take state_end as an argument for the next time 
% period of the simulation.
%
% step4:visulation
% You may wonder why "solve" does not give you all the simulation results. 
% This is because in most cases, we need to process the simulation data 
% before visualizing it.
% In this case, what we need may be a tracking trajectory as well as the 
% error curve of position errors. However, calculating the error using the 
% raw simulation data is tedious.
% Here, we can use method "addSignals" to handle this issue.
% Define a function as
% ------function signals = sig(t,state,para)
%          |---[r,~] = ref(t,para);
%          |---er = state.x - r;
%          |---signals.r = r;
%          |---signals.er = er;  
% ------end
% Next, we add the desired signal to the final result using the following 
% method. 
% -----mySystem.addSignals(@(t,s)sig(t,s,para));
% which means function "sig" will be executed at all time instants in the
% raw simulation data.
% Note that the structure returned by the equation "sig" will overwrite 
% the structure of the same name in the final result, so you'd better make
% sure that it doesn't have the same name unless you know what you're doing.
% Then, the final result is collected by
% ------[t,data] = mySystem.result()
% ------% or
% ------% [t,data] = mySystem.result(2000)
% where the parameter 2000 means that the final result consisted of 2000
% pieces of data sampled from the raw data.
%
% step5:visualization
% Finally, we can plot result as
% ------figure(1);clf;hold on
% ------plot(data.x(:,1),data.x(:,2));
% ------plot(data.r(:,1),data.r(:,2),'--');
% ------plot(data.x(1,1),data.x(1,2),'x');
% ------plot(data.r(1,1),data.r(1,2),'*');
% ------plot(data.x(end,1),data.x(end,2),'o');
% ------plot(data.r(end,1),data.r(end,2),'p');
% ------title('tracking trajectory')
% ------figure(2)
% ------plot(t,data.er)
% ------xlabel('t(s)')
% ------title('position errors')
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
        delayedFunctionHandle=[];
        pastDelayedFunctionHandle=[];
        functionName='';
        functionAfterEachTspanHandle=[];
        functionAfterEachTspanName='';
        inputName=[];
        stime=[];
        delayedResult=[];
        simulation_t=[];
        simulation_data=[];
        simulation_result=[];
        dt=[];
        options=[];
        enableSizeCheck=[];
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
                self.delayedFunctionHandle=@(t,vector,delayedVector)self.delayedModel(t,vector,delayedVector,handle);
                self.pastDelayedFunctionHandle=@(t,vector,delayedVector)self.pastDelayedModel(t,vector,delayedVector,handle);
                self.delayedResult=[];
                self.simulation_t=[];
                self.simulation_data=[];
                self.simulation_result=[];
                self.dt=0.0001;
                self.options=[];
                self.enableSizeCheck=1;
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
        
        function state_t = stateAt(self,at,stream)
            state_vec = stream(at,:);
            state_t = self.fromVector(state_vec);
        end
        
        function [t,result] = result(self,N)
            if nargin<=1
                N=-1;
            end
            t = self.simulation_t;
            if N>0 && N<length(t)
                indexing = self.sampling(t,N);
                t= t(indexing);
                sampled_data = self.simulation_data(indexing,:);
                result=self.fromVectorStream(sampled_data);
            else
                result=self.simulation_result;
            end
        end
        
        function set_dt(self,dt)
            self.dt=dt;
        end
        
        function set_size_check(self,str)
            if str==true
                self.enableSizeCheck=true;
            else
                self.enableSizeCheck=false;
            end
        end
        
        function set_options(self,op)
            self.options = op;
        end
        
        function final_state = solve(self,solver,tspans,state0)
            self.simulation_data=[];
            self.simulation_result=[];
            self.stime=tspans(end);
            self.inputName=inputname(3);
            self.initWaitBar('');
            vector = self.fromStateToVector(state0);
            t_size=length(tspans);
            stream=[];
            t=[];
            for i=1:t_size-1
                tspan_i=tspans(i:i+1);
                if isstruct(self.options)
                    op = odeset('reltol',1e-3,'abstol',1e-6);
                else
                    op = self.options;
                end
                switch lower(solver)
                    case {'ode113'}
                        [t_i,stream_i]=ode113(self.functionHandle,tspan_i,vector,op);
                    case {'ode15s'}
                        [t_i,stream_i]=ode15s(self.functionHandle,tspan_i,vector,op);
                    case {'rk4'}
                        [t_i,stream_i]=self.rk4(self.functionHandle,tspan_i,vector,self.dt);
                    case {'euler'}
                        [t_i,stream_i]=self.euler(self.functionHandle,tspan_i,vector,self.dt);
                    case {'ode45'}
                        [t_i,stream_i]=ode15s(self.functionHandle,tspan_i,vector,op);
                    otherwise
                        fprintf('unknown solver %s.\nUsing ode45 as default.\n',solver);
                        [t_i,stream_i]=ode45(self.functionHandle,tspan_i,vector,op);
                end
                vector=stream_i(end,:)';
                final_state = self.fromVector(vector);
                if ~isempty(self.functionAfterEachTspanName)
                    state_i=self.fromVector(vector);
                    state_tspan=self.functionAfterEachTspanHandle(t_i(end),state_i);
                    vector=self.toVector(state_tspan);
                end
                t=vertcat(t(1:end-1),t_i);
                stream=vertcat(stream(1:end-1,:),stream_i);
            end
            result=self.fromVectorStream(stream);
%             data=stream(end,:);
            self.simulation_t = t;
            self.simulation_data=stream;
            self.simulation_result=result;
            self.closeWaitBar();
            
        end
        
        function addSignal(self,name,functionHandle)
            state_vec_i = self.simulation_data(1,:);
            state_i = self.fromVector(state_vec_i);
            ti = self.simulation_t(1);
            x = functionHandle(ti,state_i);
            d = zeros(length(self.simulation_t),numel(x));
            d(1,:)=x(:)';
            for i=2:length(self.simulation_t)
                state_vec_i = self.simulation_data(i,:);
                state_i = self.fromVector(state_vec_i);
                ti = self.simulation_t(i);
                x = functionHandle(ti,state_i);
                d(i,:)=x(:)';
            end
            data = self.simulation_result;
            data.(name)=d;
            self.simulation_result=data;
        end
        
        function addSignals(self,functionHandle)
            state_vec_i = self.simulation_data(1,:);
            state_i = self.fromVector(state_vec_i);
            ti = self.simulation_t(1);
            signals = functionHandle(ti,state_i);
            signal_names = fieldnames(signals);
            s = length(signal_names);
            data = self.simulation_result;
            for i=1:s
                xi = signals.(signal_names{i});
                data.(signal_names{i})=zeros(length(self.simulation_t),numel(xi));
            end
            h1 = waitbar(0,'adding singnals...');
            for i=1:length(self.simulation_t)
                state_vec_i = self.simulation_data(i,:);
                state_i = self.fromVector(state_vec_i);
                ti = self.simulation_t(i);
                signals = functionHandle(ti,state_i);
                for j=1:s
                    sj = signals.(signal_names{j});
                    data.(signal_names{j})(i,:) = sj(:)';
                end
                waitbar(i/length(self.simulation_t),h1);
            end
            self.simulation_result=data;
            close(h1);
        end
        
%         function addSignalsWithDerivatives(self,functionHandle)
%             state_vec_i = self.simulation_data(1,:);
%             state_i = self.fromVector(state_vec_i);
%             ti = self.simulation_t(1);
%             signals = functionHandle(ti,state_i);
%             signal_names = fieldnames(signals);
%             s = length(signal_names);
%             data = self.simulation_result;
%             for i=1:s
%                 xi = signals.(signal_names{i});
%                 data.(signal_names{i})=zeros(length(self.simulation_t),numel(xi));
%             end
%             h1 = waitbar(0,'adding singnals...');
%             for i=1:length(self.simulation_t)
%                 state_vec_i = self.simulation_data(i,:);
%                 state_i = self.fromVector(state_vec_i);
%                 ti = self.simulation_t(i);
%                 signals = functionHandle(ti,state_i);
%                 for j=1:s
%                     sj = signals.(signal_names{j});
%                     data.(signal_names{j})(i,:) = sj(:)';
%                 end
%                 waitbar(i/length(self.simulation_t),h1);
%             end
%             self.simulation_result=data;
%             close(h1);
%         end
        
%         function data = addSignals(self,names,functionHandle)
%             state_vec_i = self.simulation_data(1,:);
%             state_i = self.fromVector(state_vec_i);
%             ti = self.simulation_t(1);
%             x = functionHandle(ti,state_i);
%             d = zeros(length(self.simulation_t),numel(x));
%             d(1,:)=x(:)';
%             for i=2:length(self.simulation_t)
%                 state_vec_i = self.simulation_data(i,:);
%                 state_i = self.fromVector(state_vec_i);
%                 ti = self.simulation_t(i);
%                 x = functionHandle(ti,state_i);
%                 d(i,:)=x(:)';
%             end
%             data = self.simulation_result;
%             data.(name)=d;
%             self.simulation_result=data;
%         end
        
        function [t,data,last]=dde23(self,tspans,state0,delay,N)
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
                if isempty(self.delayedResult)==1
                    sol_i = dde23(self.delayedFunctionHandle,delay,vector,tspan_i);
                else
                    sol_i = dde23(self.pastDelayedFunctionHandle,delay,self.delayedResult,tspan_i);
                end
                self.delayedResult = sol_i;
                t_i=sol_i.x';
                stream_i=sol_i.y';
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
        
        function dot_vector=delayedModel(self,t,vector,delayedVectors,fun)
            %             global debug_i
            %             if debug_i==2
            %                 debug_i;
            %             end
            state = self.fromVector(vector);
            delayedStates = self.fromVector(delayedVectors);
            
            grad  = fun(t,state,delayedStates);
            dot_vector=self.toVector(grad);
            self.updateWaitBar(t/self.stime);
        end
        
        function dot_vector = pastDelayedModel(self,t,vector,delayedVectors,fun)
            state = self.fromVector(vector);
            delayedStates = self.fromVector(self.delayedResult.y(:,end)');
            grad  = fun(t,state,delayedStates);
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
                if self.enableSizeCheck
                assert(isfield(state,names{i}),['No variable named "',names{i}...
                    ,'" found. Please check if "',self.inputName,'" contains variable "',names{i},'".']);
                end
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                s1 = self.size_state.(names{i});
                s2 = size(state.(names{i}));
                if self.enableSizeCheck
                assert(isequal(s1,s2),['You are trying to change the size of "',names{i},'", please check it.'...
                    ,' The size of "',names{i},'" should be ',mat2str(s1),', but the size of "'...
                    ,self.inputName,'.',names{i},'" is ',mat2str(s2),'.']);
                end
                vector(a:b)=state.(names{i})(:);
            end
        end
        
        function vector = toVector(self,state)
            vector = zeros(self.size_vector,1);
            names = fieldnames(self.size_state);
            k = length(names);
            for i=1:k
                if self.enableSizeCheck
                assert(isfield(state,names{i}),['No variable named "',names{i}...
                    ,'" found. Please check if the structural variable returned from function "'...
                    ,self.functionName,'" contains variable "',names{i},'".']);
                end
                a = self.index_a.(names{i});
                b = self.index_b.(names{i});
                s1 = self.size_state.(names{i});
                s2 = size(state.(names{i}));
                if self.enableSizeCheck
                assert(isequal(s1,s2),['You are trying to change the size of "',names{i},'", please check it.'...
                    ,' The size of "',names{i},'" should be ',mat2str(s1),', but the size of "'...
                    ,names{i},'" returned from function "',self.functionName,'" is ',mat2str(s2),'.']);
                end
                vector(a:b)=state.(names{i})(:);
            end
        end
        
        function state = fromVector(self,vector)
            state=[];
            vectorname = inputname(2);
            s = length(vector);
            if self.enableSizeCheck
            assert(s==self.size_vector,['The size of input vector should to be ',mat2str(self.size_vector),', but the length of "'...
                ,vectorname,'" is ',num2str(s),', please check it']);
            end
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
            if self.enableSizeCheck
            assert(s==self.size_vector,['The column size of input vector stream should to be ',mat2str(self.size_vector),', but the column size of "'...
                ,streamname,'" is ',num2str(s),', please check it.']);
            end
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
                k2=ufunc(t(iter)+dt/2,stream(iter,:)+dt/2*k1);
                k3=ufunc(t(iter)+dt/2,stream(iter,:)+dt/2*k2);
                k4=ufunc(t(iter)+dt,stream(iter,:)+dt*k3);
                stream(iter+1,:)=stream(iter,:)+(dt*(k1+2.*k2+2.*k3+k4)/6)';
            end
        end
        
        function [t,stream]=euler(ufunc,tspan,state0,dt)
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
                stream(iter+1,:)=stream(iter,:)+dt*(k1)';
            end
        end
    end
    
end

