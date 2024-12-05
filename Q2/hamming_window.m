function w = hamming_window(N)
    % Validate input
    if N <= 0
        error('N must be a positive integer.');
    end
    
    % Compute Hamming window
    n = 0:N-1; % Indices
    w = 0.54 - 0.46 * cos(2 * pi * n / (N-1));
end
