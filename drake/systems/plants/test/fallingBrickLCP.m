function fallingBrickLCP

options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.update_convex = true;
options.ignore_self_collisions = true;
% options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator(s,.01,options);
p = p.addRobotFromURDF(s,[],[],options);
% x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];
x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
start = x0;
x0 = p.resolveConstraints(x0);

if 0
    v = p.constructVisualizer();
    sys = cascade(p,v);
    sys.simulate([0 8],x0);
    return;
end

v = p.constructVisualizer();
v.drawWrapper(0,x0);


%You should support this functionality
% N = 20;
% cto = ContactTrajectoryOptimization(p,N,[0 4])
% cto = cto.addStateConstraint(ConstantConstraint(start(1:24)),1);
% cto = cto.addRunningCost(@cost);
% 
%     function [g,dg] = cost(dt,x,u)
%         R = 1;
%         g = sum((R*u).*u,1);
%         dg = [zeros(1,1+size(x,1)),2*u'*R];
%     end


% cto = cto.addFinalCost(@finalCost);
% traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
%
% for attempts=1:10
%     tic
%     [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
%     toc
%     if info==1, break; end
% end



xtraj = p.simulate([0 0.5],x0);
for i = 1:size(xtraj.xx,2) - 2
    p.inverseDynamics(xtraj.tt(i), xtraj.tt(i+1), xtraj.tt(i+2), xtraj.xx(1:12,i), xtraj.xx(1:12,i+1), xtraj.xx(1:12,i+2));
end

v.playback(xtraj);
end