close all
clear
init
mu = Gearth;
af = 10000;
incf = 0;
eccf = 0.0001;
wf = 35;
hf = HfromAe(mu,af,eccf);
Raanf = 0;

r0 = [8276; 5612; 5];
v0 = [-3.142; 4.672; 0];

ai = AfromVR(mu,norm(v0),norm(r0));
Hi = cross(r0,v0);

Ecci = EccfromHRV(mu,r0,v0,Hi);
ecci = norm(Ecci);
inci = IfromH(Hi);
LON = cross([0 0 1],Hi);
n_mag = norm(LON);  
Raani = RaanfromLON(LON) * raan_sign(LON);
%Raani = calculate_RAAN(ai,ecci,inci);

wi = WfromLONe(LON,Ecci)*w_sign(Ecci);
TAi = TAfromRe(r0,Ecci);

orbit_3Di = plot_orbit(ai,ecci,inci,Raani,wi,TAi);
orbit_3Df = plot_orbit(af,eccf,incf,Raanf,wf,0);

figure(1) 
    % Plot the orbit
    plot3(orbit_3Di(1, :), orbit_3Di(2, :), orbit_3Di(3, :), 'b','DisplayName','Initial Orbit');
    hold on;
    plot3(orbit_3Df(1, :), orbit_3Df(2, :), orbit_3Df(3, :), 'g','DisplayName','Final Orbit');
    plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r','DisplayName','Earth'); % Central body
    plot3(r0(1),r0(2),r0(3),'bo','DisplayName','s/c at t_0')
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Orbit Plot');
    legend('Location','bestoutside')
    zlim([-50 50]);


function orbit_3D = plot_orbit(a, e, i, RAAN, w, theta)
    % a: semi-major axis (km)
    % e: eccentricity
    % i: inclination (degrees)
    % RAAN: right ascension of ascending node (degrees)
    % w: argument of periapsis (degrees)
    % theta: true anomaly (degrees)

    % Convert angles from degrees to radians
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    w = deg2rad(w);
    theta = deg2rad(theta);

    % Define the number of points in the orbit
    num_points = 1000;
    theta_vals = linspace(0, 2*pi, num_points);

    % Calculate the radius for each point in the orbit
    r = a * (1 - e^2) ./ (1 + e * cos(theta_vals));

    % Calculate the position in the orbital plane
    x_orbit = r .* cos(theta_vals);
    y_orbit = r .* sin(theta_vals);

    % Rotation matrices
    R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];

    % Combined rotation matrix
    R = R3_W * R1_i * R3_w;

    % Initialize the orbit in 3D space
    orbit_3D = zeros(3, num_points);

    for j = 1:num_points
        % Position in the orbital plane
        pos_orbit = [x_orbit(j); y_orbit(j); 0];
        % Rotate to 3D space
        pos_3D = R * pos_orbit;
        orbit_3D(:, j) = pos_3D;
    end
end
