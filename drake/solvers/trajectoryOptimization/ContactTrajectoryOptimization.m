classdef ContactTrajectoryOptimization < DirectTrajectoryOptimization
    % Todorov 2011 Trajectory Optimization
    
%     properties(SetAccess = immutable)
%         TSRBM
%     end
    
    methods
        function obj = ContactTrajectoryOptimization(TSRBM,N,duration)
            if ~strcmp(class(TSRBM),'TimeSteppingRigidBodyManipulator')
                error('Error. \n Inverse Dynamics only defined for TimeSteppingRigidBodyManipulator, not for %s.',class(TSRBM))
            end
            
            %I probably have to live with the u_inds and the general setup.
            % What I don't understand is how I'm going to perform the
            % optimization without gradient descent on the u_inds.  I
            % have not added any dynamics or u dependency.  So really the
            % solver should just return a correct sequence of x's.  Then I
            % need to plug this into reconstruct trajectory.
            %I hope it doesn't take anymore time to include all the
            %variables that don't do anything.
            
            if nargin < 4
                options = struct();
            end
            obj = obj@DirectTrajectoryOptimization(TSRBM,N,duration,options);
            %obj.TSRBM = manipulator;
            %I need to add in the dynamics after I add in the manipulator.
            %I need to add in the manipulator after the call to direct
            %trajectory optimization by matlab syntax.
        end
        
        
        %TODO overwrite get initial vars, can't include the u's.
        
        function obj = addDynamicConstraints(obj)
            %TODO
            %You must add the forward dynamics.  Which depends on the
            %inverse dynamics.
            %called in constructor of DirectTrajectoryOptimization
            manip = obj.plant;
            
            nX = manip.getNumStates();
            nU = manip.getNumInputs();
            N = obj.N;
            %num_pos = manip.num_positions;
            num_pos = nX/2;
            
            %constraints = cell(N-1,1);
            %dyn_inds = cell(N-1,1);
            constraints = cell(N-1,1);
            dyn_inds = cell(N-1,1);
            
            %second state constraint is prespecified by the input velocity
            %x0_dot = (x1 - x0) /h0
            m = size(obj.x_inds);
            n_vars = 1 + 3/2*nX;
            cnstr = FunctionHandleConstraint(zeros(nX/2,1),zeros(nX/2,1),n_vars,@obj.finite_difference);
            %once derivatives are implemented, perhaps remove this line
            cnstr.grad_method = 'numerical';
            dyn_inds{1} = {obj.h_inds(1);obj.x_inds(1:num_pos,1);obj.x_inds(num_pos+1:m,1);obj.x_inds(1:num_pos,2)};
            constraints{1} = cnstr;
            obj = obj.addConstraint(constraints{1},dyn_inds{1});
            
            %you'll have to change these lines 
            n_vars = 2 + 3/2*nX;
            %cnstr = FunctionHandleConstraint(zeros(nX/2,1),zeros(nX/2,1),n_vars,@obj.forward_constraint_fun);
            cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
            %once derivatives are implemented, perhaps remove this line
            cnstr.grad_method = 'numerical';
            
            %n_vars = 2 + 3/2*nX;
            %cnstr = FunctionHandleConstraint(zeros(nX/2,1),zeros(nX/2,1),n_vars,@obj.forward_constraint_fun);
            
            %we're only minimizing over states not over
            %velocities/derivatives 
            %for i=1:obj.N-3,
            for i = 1:obj.N - 3
                dyn_inds{i+1} = {obj.h_inds(i); obj.h_inds(i+1); obj.x_inds(1:num_pos,i);obj.x_inds(1:num_pos,i+1); obj.x_inds(1:num_pos,i+2)};
                
                constraints{i+1} = cnstr;
                
                obj = obj.addConstraint(constraints{i+1}, dyn_inds{i+1});
            end
        end
        
        %TODO
        function [xtraj,utraj,z,F,info,infeasible_constraint_name] = solveTraj(obj,t_init,traj_init)
            % Solve the nonlinear program and return resulting trajectory
            % @param t_init initial timespan for solution.  can be a vector of
            % length obj.N specifying the times of each segment, or a scalar
            % indicating the final time.
            % @param traj_init (optional) a structure containing Trajectory
            % objects specifying the initial guess for the system inputs
            %    traj_init.u , traj_init.x, ...
            % @default small random numbers
            
            if nargin<3, traj_init = struct(); end
            
            z0 = obj.getInitialVars(t_init,traj_init);
            
            %since neither costs nor constraint depend on u_inds, I can
            %only assume this is safe.  Especially with Scott's comment
            %about how number of constraints should roughly scale linearly.
            % Seems a little opaque.  Maybe write a note.  
            [z,F,info,infeasible_constraint_name] = obj.solve(z0);
            xtraj = reconstructStateTrajectory(obj,z);
            if nargout>1, utraj = reconstructInputTrajectory(obj,z); end
        end
        
        %MUST BE IMPLEMENTED
        %Here we take running cost function and take advantage of inverse
        %dynamics.
        %So what do you want from a running cost.  The running cost
        %function I am passed has u's in it, which I must replace with
        %inverse dynamics.  I'm adding this running cost function for which
        %time intervals?
        function obj = addRunningCost(obj,running_cost_function)
            % Adds an integrated cost to all time steps, which is
            % numerical implementation specific (thus abstract)
            % this cost is assumed to be time-invariant
            % @param running_cost_function a function handle
            %  of the form running_cost_function(dt,x,u)
            % This implementation assumes a ZOH, but where the values of
            % x(i),u(i) are held over an interval spanned by .5(dt(i-1) + dt(i))
            
            manip = obj.TSRBM.getManipulator();
            
            nX = manip.getNumStates();
            %nU = obj.plant.getNumInputs();
            
            %             q = x(1:obj.num_positions);
            %             v = x(obj.num_positions+1:end);
            
            
            %running_cost_end = FunctionHandleObjective(1+nX,@(h,x) obj.running_fun_end(running_cost_function,h,x,u));
            %nX/2 + nX/2 = nX
            %it was 2 + nX but that's just not correct.
            cost_handle = FunctionHandleObjective(2+3/2*nX,@(h0,h1,x0,x1,x2) obj.cost_inverse_dynamics(running_cost_function,h0,h1,x0,x1,x2));
            
            %obj = obj.addCost(running_cost_end,{obj.h_inds(1);obj.x_inds(:,1)});
            
            %adding cost will have to accomodate the lack of u's.  how?
            %       for i=2:obj.N-1,
            %         obj = obj.addCost(running_cost_mid,{obj.h_inds(i-1);obj.h_inds(i);obj.x_inds(:,i)});
            %       end
            
            %TODO consider what to do with first iteration
            num_pos = manip.num_positions;
            
            for i=1:obj.N-3
                inds_i = {obj.h_inds(i);obj.h_inds(i+1); obj.x_inds(1:num_pos,i);obj.x_inds(1:num_pos,i+1); obj.x_inds(1:num_pos,i+2)};
                obj = obj.addCost(cost_handle,inds_i);
            end
            %Address last case
        end
        
        function xtraj = reconstructStateTrajectory(obj,z)
            % default behavior is to use first order holds, but this can be
            % re-implemented by a subclass.
            manip = obj.plant.getManipulator();
            num_pos = manip.num_positions;
            
            t = [0; cumsum(z(obj.h_inds))];
            
            x = reshape(z(obj.x_inds),[],obj.N);
            [m,n] = size(x);
            %finite differencing for velocity. What about the last one?
            for i = 1:size(x,2)-1
                try
                    v = (x(1:num_pos,i+1) - x(1:num_pos,i))/(obj.h_inds(i));
                catch
                    v = (x(1:num_pos,i+1) - x(1:num_pos,i))/(obj.h_inds(i));
                end
                x(num_pos+1:m,i) = v;
            end
            x(num_pos+1:m,n) = v;
            xtraj = PPTrajectory(foh(t,x));
            xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
        end
        
        %methods(Access = protected)
        function f = cost_inverse_dynamics(obj,running_handle,h0,h1,x0,x1,x2)
            
            %I could try to fix inverse dynamics right now.  No doubt it is
            %bugridden.  But unecessary just to get the brick to fall.
            %Ignore it for now.
            %             try
            %                 u = obj.TSRBM.inverseDynamics(h0,h1,x0,x1,x2);
            %             catch
            %                disp('VALIDATION')
            %                u = obj.TSRBM.inverseDynamics(h0,h1,x0,x1,x2);
            %             end
            u = 0;
            f = running_handle(h0,x0,u);
        end
        
        function f = forward_constraint_fun(obj,h0,h1,x0,x1,x2)
            [m,n] = size(x0);
            nX = 2*m;
            
            x0_dot = (x1 - x0)/h0;
            x1_dot = (x2 - x1)/h1;
            xmdot = .5*(x0_dot + x1_dot);
            
            x0_full = [x0;x0_dot];
            x1_full = [x1;x1_dot];
            u = zeros(0,1);
            
            
            %[u,~,~,~,~,~] = obj.plant.inverseDynamics(h0,h1,x0,x1,x2);
            
            [x_next,dxdot] = obj.plant.update(0,.5*(x0_full + x1_full),u);
            xmdotdot = (x_next(m+1:2*m) - xmdot)/obj.plant.timestep;
            %f = (x1_dot - x0_dot) - h0*(x_next(m+1:2*m) - .5*(x0_dot + x1_dot))/obj.plant.timestep;
            xdot = [xmdot;xmdotdot];
            f = (x1_full - x0_full) - h0*xdot;
            
            
            
            
            %We are doing this like Dirtran because we currently don't know
            %how the inverse dynamics really works.  
            %df = [-xdot (-eye(nX) - .5*h0*dxdot(:,2:1+nX)) (eye(nX)- .5*h0*dxdot(:,2:1+nX)) -.5*h0*dxdot(:,nX+2:end) -.5*h0*dxdot(:,nX+2:end)];
            
        end
        
        function f = finite_difference(obj,h0,x0,x0_dot,x1)
            %f = (x1 - x0)/h0 - x0_dot;
            f = (x1 - x0) - h0*x0_dot;
            %f = x1 - x1;
        end
    end
    %end
end



