function runILQRSwingUp()
% runs trajectory optimization and animates 

pd = PendulumPlant;
pv = PendulumVisualizer();
[utraj,xtraj] = iLQRswingUpTrajectory(pd);

uhist = utraj.eval(0:.1:4);

figure(1);
plot(uhist);

pv.playback(xtraj);

end
