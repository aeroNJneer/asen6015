        % Helper function to compute RTN to ECI rotation matrix
        function R = compute_RTN_to_ECI_rotation(r, v)
            % R (radial) axis is along position vector
            R_axis = r / norm(r);
            
            % H is angular momentum vector
            h_vec = cross(r, v);
            H_axis = h_vec / norm(h_vec);
            
            % T (transverse) is perpendicular to R and H
            T_axis = cross(H_axis, R_axis);
            
            % Rotation matrix
            R = [R_axis, T_axis, H_axis]';
        end
