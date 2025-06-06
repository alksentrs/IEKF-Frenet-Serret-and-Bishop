function S = skew(v)
    % Compute skew-symmetric matrix from vector
    % Input:
    %   v - 3x1 vector
    % Output:
    %   S - 3x3 skew-symmetric matrix
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end 