init

delV = 0;
delR = 0;
delA = 0;
dt = 1000;
t = 0;
timespan = linspace(t,t+dt,dt);

% numerical integration
for i = 0:length(timespan)
    y[] = [delr;delV]
    f[] = [delV, delA]

    find 
        delR(t+dt) and delV(t+dt

        compute rosc, vosc (3.4)
        compute ri = roscc + delr(i) (10.6)
        compute vi = vosc(i) +delv(i)
        compute dela (10.7)
        dvi , delr to next iteration

end

