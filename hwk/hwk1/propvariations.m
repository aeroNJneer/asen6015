init;

dt = 1;
COE = [Rearth+200, 0, 0, 0, 0, 0]';


pert0 = [0, 2.3357/(Rearth+200), 0];

% find orbital evolution through time

timespan = 0:dt:10000;
pertini = zeros(length(timespan),3);
pert = [pert0; pertini];
dCOE = [];
for i = 1:length(timespan)
    COE(:,i+1) = propogate_orbital_elements(COE(:,i),Gearth,pert(i,:),dt);
end
