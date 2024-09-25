function I = midpoint_rule(f, a, b, n)
    % f: function handle of the integrand
    % a: lower limit of integration
    % b: upper limit of integration
    % n: number of subintervals

    % Calculate the width of each subinterval
    dx = (b - a) / n;

    % Initialize the sum
    midpoint_sum = 0;

    % Loop over each subinterval
    for i = 1:n
        % Calculate the midpoint of the i-th subinterval
        xi = a + (i - 0.5) * dx;
        
        % Add the function value at the midpoint to the sum
        midpoint_sum = midpoint_sum + f(xi);
    end

    % Multiply the sum by the width of the subintervals to get the integral
    I = dx * midpoint_sum;
end
