classdef iLQRTrajectoryOptimization < IndirectTrajectoryOptimization
  %ILQRTRAJECTORYOPTIMIZATION Performs trajectory optimization with the
  % Itterative LQR Algorithm.
  %
  % Generally considers cost functions of the form:
  % int(f(x(t),u(t)) + g(xf)
  %
  % This class assumes that there are a fixed number (N) time steps of
  % equal length (h).
  %
  % This implementation assumes that all constraints and costs are
  % time-invariant.
  
    properties (Constant)
        EULER = 1;
        MIDPOINT = 2;
        RK3 = 3; % DEFAULT
        RK4 = 4;
    end
    
    properties (Access = private)
        integrator
        input_lower_bound
        input_upper_bound
    end
    
    methods
        function obj = iLQRTrajectoryOptimization(plant,N,duration,options)
            
            defaults = struct(...
                'tolFun',         1e-7,...  minimum cost reduction
                'tolGrad',        1e-5,...  minimum gradient size
                'maxIter',        500,...   maximum iterations
                'lambdaFactor',   1.6,...   lambda scaling factor
                'lambdaMax',      1e8,...   maximum lambda value
                'lambdaMin',      1e-6,...  minimum lambda value
                'alphaMin',       1e-3,...  minimum line search parameter
                'integration_method', iLQRTrajectoryOptimization.RK3);
            
            if nargin < 4
                options = defaults;
            end
            
            obj = obj@IndirectTrajectoryOptimization(plant,N,duration,options);
            
            %Figure out which integrator to use
            switch obj.options.integration_method
                case iLQRTrajectoryOptimization.EULER
                    obj.integrator = @obj.euler;
                case iLQRTrajectoryOptimization.MIDPOINT
                    obj.integrator = @obj.midpoint;
                case iLQRTrajectoryOptimization.RK3
                    obj.integrator = @obj.rk3;
                case iLQRTrajectoryOptimization.RK4
                    obj.integrator = @obj.rk4;
                otherwise
                    error('Drake:iLQRTrajectoryOptimization:InvalidArgument','Unknown integration method');
            end
        end
        
        function obj = addInputLimits(obj,lower,upper)
            % Add upper and lower bounds to control inputs
            obj.input_lower_bound = lower;
            obj.input_upper_bound = upper;
        end
        
        function [xtraj,utraj,J,info] = solveTraj(obj,tf,traj_init)
            
            n = obj.plant.getNumContStates();
            m = obj.plant.getNumInputs();
            N = obj.getN();
            
            %Allocate matrices
            J = 0;
            A = zeros(n, n, N-1);
            B = zeros(n, m, N-1);
            Q = zeros(n+m, n+m, N);
            q = zeros(n+m, N);
            x = zeros(n, N);
            u = zeros(m, N-1);
            t = zeros(N, 1);
            h = (obj.duration/(N-1))*ones(N-1, 1);
            
            %Set Initial Values
            x(:,1) = traj_init.x.eval(0);
            lambda = 1;
            dlambda = 1;
            
            for iter = 1:obj.options.maxIter
                %Simulate forward and evaluate derivatives of dynamics and cost
                J = 0;
                for k = 1:(N-1)
                    [l, Q(:,:,k), q(:,k)] = obj.running_cost(h(k),x(:,k),u(:,k));
                    J = J + l;
                    t(k+1) = t(k) + h(k);
                    [x(:,k+1), A(:,:,k), B(:,:,k)] = obj.integrator(x(:,k),u(:,k),h(k));
                end
                [l, Q(1:n,1:n,N), q(1:n,N)] = obj.final_cost(t(N),x(:,N));
                J = J + l;
                
                %Do backward pass to calculate du and K
                done = 0;
                while ~done
                    [du,K,cholFail] = obj.backwardPass(t,x,u,A,B,Q,q,lambda);
                    if cholFail
                        %Increase lambda
                        dlambda = max(dlambda*obj.options.lambdaFactor, obj.options.lambdaFactor);
                        lambda = max(lambda*dlambda, obj.options.lambdaMin);
                        if lambda > obj.options.lambdaMax
                            break;
                        end
                    else
                        done = 1;
                    end
                end
                
                %Check gradient norm to see if we should stop
                g_norm = mean(max(abs(du)./(abs(u)+1),[],1));
                if g_norm < obj.options.tolGrad && lambda < 1e-5
                    break;
                end
                
                %Do forward pass with line search to find new J, x, and u
                alpha = 1;
                done = 0;
                while ~done
                    [xnew,unew,Jnew] = obj.forwardPass(t,x,u,h,du,K,alpha);
                    deltaJ = J - Jnew;
                    
                    if deltaJ > 0 
                        % Success
                        x = xnew;
                        u = unew;
                        J = Jnew;
                        
                        %Try to decrease lambda
                        dlambda = min(dlambda/obj.options.lambdaFactor, 1/obj.options.lambdaFactor);
                        lambda = lambda*dlambda*(lambda > obj.options.lambdaMin);
                        break;
                    end
                    
                    if alpha > obj.options.alphaMin
                        alpha = alpha/2;
                    else
                        % Failed - increase lambda
                        deltaJ = obj.options.tolFun;
                        dlambda = max(dlambda*obj.options.lambdaFactor, obj.options.lambdaFactor);
                        lambda = max(lambda*dlambda, obj.options.lambdaMin);
                        break;
                    end
                end
                
                %Check deltaJ to see if we should stop
                if deltaJ < obj.options.tolFun
                    break;
                end
                
                %Check regularization parameter to see if we should stop
                if lambda > obj.options.lambdaMax
                    break;
                end
            end
            
            xtraj = obj.reconstructStateTrajectory(x);
            if nargout>1
                utraj = obj.reconstructInputTrajectory(u);
            end
            
            
        end
    end
    
    methods (Access = private)
        function [du,K,error] = backwardPass(obj,t,x,u,A,B,Q,q,lambda)
            
            %Set up backwards LQR pass
            n = size(x,1);
            m = size(u,1);
            N = size(x,2);
            
            P = zeros(n, n, N);
            P(:,:,N) = Q(1:n,1:n,N);
            
            p = zeros(n, N);
            p(:,N) = q(1:n,N);
            
            K = zeros(m,n,N-1);
            du = zeros(m,N-1);
            
            for k = (N-1):-1:1
                
                % Propagate cost-to-go function through linearized dynamics
                H = [A(:,:,k) B(:,:,k)]'*P(:,:,k+1)*[A(:,:,k) B(:,:,k)];
                g = [A(:,:,k) B(:,:,k)]'*p(:,k+1);
                
                % Add on running cost for this time step
                H = H + Q(:,:,k);
                g = g + q(:,k);
                
                % Break cost to go matrix/vector into convenient blocks
                Hxx = H(1:n,1:n);
                Huu = H((n+1):end,(n+1):end);
                Hxu = H(1:n,(n+1):end);
                Hux = H((n+1):end,1:n);
                gx = g(1:n);
                gu = g((n+1):end);
                
                if isempty(obj.input_lower_bound) && isempty(obj.input_upper_bound)
                    % No bounds - solve for du and K in closed form
                    
                    % Add regularization and check for positive definiteness
                    [Suu,error] = chol(Huu + lambda*eye(m));
                    if error
                        break;
                    end
                
                    duK = Suu\(Suu'\[gu Hux]);
                    du(:,k) = -duK(:,1);
                    K(:,:,k) = duK(:,2:end);
                else
                    % Solve QP to find du and K
                    lower = obj.input_lower_bound - u(:,k);
                    upper = obj.input_upper_bound - u(:,k);
                    
                    [du(:,k),free,Suu,error] = obj.boxQP(Huu+lambda*eye(m),gu,lower,upper,du(:,k));
                    if error
                        break;
                    end
                    
                    K(:,:,k) = zeros(m,n);
                    if any(free)
                        Kfree = Suu\(Suu'\Hux(free,:));
                        K(free,:,k) = Kfree;
                    end
                end
                
                %Calculate new cost-to-go function
                p(:,k) = Hxu*du(:,k) - K(:,:,k)'*Huu*du(:,k) + gx - K(:,:,k)'*gu;
                P(:,:,k) = Hxx + K(:,:,k)'*Huu*K(:,:,k) - Hxu*K(:,:,k) - K(:,:,k)'*Hux;
                
            end
        end
        
        function [xnew,unew,Jnew] = forwardPass(obj,t,x,u,h,du,K,alpha)
            unew = zeros(size(u));
            xnew = zeros(size(x));
            xnew(:,1) = x(:,1);
            N = size(x,2);
            Jnew = 0;
            for k = 1:N-1
                Jnew = Jnew + obj.running_cost(h(k),xnew(:,k),u(:,k));
                unew(:,k) = u(:,k) + alpha*du(:,k) - K(:,:,k)*(xnew(:,k)-x(:,k));
                xnew(:,k+1) = obj.integrator(xnew(:,k),unew(:,k),h(k));
            end
            Jnew = Jnew + obj.final_cost(t(N), xnew(:,N));
        end
        
        function [x,free,Hfree,error] = boxQP(obj,H,g,lower,upper,x0)
            % Minimize 0.5*x'*H*x + x'*g  s.t. lower<=x<=upper
            %
            %  inputs:
            %     H            - positive definite matrix   (n * n)
            %     g            - bias vector                (n)
            %     lower        - lower bounds               (n)
            %     upper        - upper bounds               (n)
            %
            %   optional inputs:
            %     x0           - initial state              (n)
            %
            %  outputs:
            %     x            - solution                   (n)
            %     Hfree        - subspace cholesky factor   (n_free * n_free)
            %     free         - set of free dimensions     (n)
            
            n        = size(H,1);
            clamped  = false(n,1);
            free     = true(n,1);
            nfactor  = 0;
            Hfree    = zeros(n);
            clamp    = @(x) max(lower, min(upper, x));
            error = 0;
            
            % initial state
            if nargin > 5 && numel(x0) == n
                x = clamp(x0(:));
            else
                x = clamp(zeros(n,1));
            end
            x(~isfinite(x)) = 0;
            
            % Algorithm parameters
            maxIter        = 100;       % maximum number of iterations
            minGrad        = 1e-8;      % minimum norm of non-fixed gradient
            minRelImprove  = 1e-8;      % minimum relative improvement
            stepDec        = 0.6;       % factor for decreasing stepsize
            minStep        = 1e-16;     % minimal stepsize for linesearch
            Armijo         = 0.1;   	% Armijo parameter (fraction of linear improvement required)
            
            % initial objective value
            value    = x'*g + 0.5*x'*H*x;
            
            % main loop
            for iter = 1:maxIter
                
                % get gradient
                grad = g + H*x;
                
                % find clamped dimensions
                old_clamped                     = clamped;
                clamped                         = false(n,1);
                clamped((x == lower)&(grad>0))  = true;
                clamped((x == upper)&(grad<0))  = true;
                free                            = ~clamped;
                
                % check for all clamped
                if all(clamped)
                    break;
                end
                
                % factorize if clamped has changed
                if iter == 1
                    factorize = true;
                else
                    factorize = any(old_clamped ~= clamped);
                end
                
                if factorize
                    [Hfree, indef]  = chol(H(free,free));
                    if indef
                        error = 1;
                        break
                    end
                    nfactor = nfactor + 1;
                end
                
                % check gradient norm
                gnorm  = norm(grad(free));
                if gnorm < minGrad
                    break;
                end
                
                % get search direction
                grad_clamped   = g  + H*(x.*clamped);
                search         = zeros(n,1);
                search(free)   = -Hfree\(Hfree'\grad_clamped(free)) - x(free);
                
                % check for descent direction
                sdotg          = sum(search.*grad);
                if sdotg >= 0 % (should not happen)
                    error = 1;
                    break
                end
                
                % armijo linesearch
                step  = 1;
                nstep = 0;
                xc    = clamp(x+step*search);
                vc    = xc'*g + 0.5*xc'*H*xc;
                while (vc - value)/(step*sdotg) < Armijo
                    step  = step*stepDec;
                    nstep = nstep+1;
                    xc    = clamp(x+step*search);
                    vc    = xc'*g + 0.5*xc'*H*xc;
                    if step<minStep
                        break
                    end
                end
                
                % check relative improvement
                if(value - vc < minRelImprove*abs(value) )
                    x = xc;
                    break;
                else
                    value = vc;
                    x = xc;
                end
            end
        end
        
        % All of these integrators integrate Jacobians as well as states
        function [xn,A,B] = euler(obj,xk,uk,hk)
            if nargout == 1
                xdot = obj.plant.dynamics(0,xk,uk);
                xn = xk + hk*xdot;
            else
                [xdot, dxdot] = obj.plant.dynamics(0,xk,uk);
                xn = xk + hk*xdot;
                n = length(xk);
                A = eye(n) + full(dxdot(:,2:(n+1)));
                B = full(dxdot(:,(n+2):end));
            end
        end
        
        function [xn,A,B] = midpoint(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xn = xk + hk*xdot2;
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xn = xk + hk*xdot2;
                
                A1 = eye(n) + (hk/2)*full(dxdot1(:,2:(n+1)));
                B1 = (hk/2)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/2)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/2)*full(dxdot2(:,(n+2):end));
                
                A = A2*A1;
                B = A2*B1 + B2;
            end
        end
        
        function [xn,A,B] = rk3(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xk+.75*hk*xdot2,uk);
                xn = xk + (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xk+.75*hk*xdot2,uk);
                xn = xk + (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
                
                A1 = eye(n) + (2*hk/9)*full(dxdot1(:,2:(n+1)));
                B1 = (2*hk/9)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) + (4*hk/9)*full(dxdot3(:,2:(n+1)));
                B3 = (4*hk/9)*full(dxdot3(:,(n+2):end));
                
                A = A3*A2*A1;
                B = A3*A2*B1 + A3*B2 + B3;
            end
        end
        
        function [xn,A,B] = rk4(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xk+.5*hk*xdot2,uk);
                xdot4 = obj.plant.dynamics(0,xk+hk*xdot3,uk);
                
                xn = xk + (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xk+.5*hk*xdot2,uk);
                [xdot4, dxdot4] = obj.plant.dynamics(0,xk+hk*xdot3,uk);
                
                xn = xk + (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
                
                A1 = eye(n) + (hk/6)*full(dxdot1(:,2:(n+1)));
                B1 = (hk/6)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) + (hk/3)*full(dxdot3(:,2:(n+1)));
                B3 = (hk/3)*full(dxdot3(:,(n+2):end));
                
                A4 = eye(n) + (hk/6)*full(dxdot4(:,2:(n+1)));
                B4 = (hk/6)*full(dxdot4(:,(n+2):end));
                
                A = A4*A3*A2*A1;
                B = A4*A3*A2*B1 + A4*A3*B2 + A4*B3 + B4;
            end
        end
    end
    
end

