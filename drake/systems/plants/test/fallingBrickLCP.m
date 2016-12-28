function fallingBrickLCP

options.floating = true;
%comment out for no ground no ground
%options.terrain = RigidBodyFlatTerrain();
options.update_convex = true;
options.ignore_self_collisions = true;
% options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator(s,.01,options);
p = p.addRobotFromURDF(s,[],[],options);
% x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];
%for i = 1:100
     x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
     x0 = p.resolveConstraints(x0);
%end
start = x0;
if 0
    v = p.constructVisualizer();
    sys = cascade(p,v);
    sys.simulate([0 8],x0);
    return;
end

v = p.constructVisualizer();
v.drawWrapper(0,x0);


%You should support this functionality
N = 20;
cto = ContactTrajectoryOptimization(p,N,[0.1 0.3]);
cto = cto.addStateConstraint(ConstantConstraint(start(1:24)),1);
% cto = cto.setSolverOptions('snopt','iterationslimit',1);
% cto = cto.setSolverOptions('snopt','majoriterationslimit',1);
% cto = cto.setSolverOptions('snopt','minoriterationslimit',1);
xtraj = p.simulate([0.0 0.2],x0);
x_old = xtraj;

v.playback(xtraj);

%just remove one layer of comments
[m,n] = size(xtraj.xx);
xf = xtraj.xx(:,n);
%cto = cto.addStateConstraint(ConstantConstraint(xf),N);
%cto = cto.addRunningCost(@cost);

    function [g,dg] = cost(dt,x,u)
        R = 1;
        g = sum((R*u).*u,1);
        dg = [zeros(1,1+size(x,1)),2*u'*R];
    end
tf0 = 0.2;
traj_init.x = PPTrajectory(foh([0,tf0],[double(start),double(xf)]));

 for attempts=1:10
     tic
     [xtraj,utraj,z,F,info] = cto.solveTraj(tf0,traj_init);
     toc
     break;
     %if info==1 || info == 3 || info == 41, break; end
 end


%  for i = 1:N-2
%      a = 24*(i-1);
%      b = a + 24;
%      c = b + 24;
%      x0 = xtraj.pp.coefs(a+1:a+12,2);
%      x1 = xtraj.pp.coefs(b+1:b+12,2);
%      x2 = xtraj.pp.coefs(c+1:c+12,2);
%      h0 = xtraj.pp.breaks(2);
%      h1 = xtraj.pp.breaks(2);
%      f = cto.forward_constraint_fun(h0,h1,x0,x1,x2);
%  end

 for i = 1:N-2
     a = 24*(i-1);
     b = a + 24;
     x0 = xtraj.pp.coefs(a+1:a+24,2);
     x1 = xtraj.pp.coefs(b+1:b+24,2);
     h = xtraj.pp.breaks(2);
     %f = cto.forward_constraint_fun(h,x0,x1,zeros(0,1));
 end

%xtraj = p.simulate([0 0.5],x0);



% for i = 1:size(x_old.xx,2) - 2
%     [u,f,J,H,B,C] = p.inverseDynamics(x_old.tt(i), x_old.tt(i+1), x_old.xx(1:12,i), x_old.xx(1:12,i+1), x_old.xx(1:12,i+2));
% end

v.playback(xtraj);
end