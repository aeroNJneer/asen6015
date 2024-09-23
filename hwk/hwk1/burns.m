dv1 = 2.454586542816101;
dv2 = 1.477271410678098;
acceleration = 0.03;

burnangle = deg2rad(45);
        Dv1 = [dv1*sind(burnangle) dv1*cosd(burnangle) 0];
        Dv2 = [dv2*sind(burnangle) dv2*cosd(burnangle) 0];

duration1 = dv1/acceleration;
duration2 = dv2/acceleration;
compare_burns(Dv1,duration1,acceleration)
