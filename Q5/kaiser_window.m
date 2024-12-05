function w = kaiser_window(N, alpha)
    % Validate inputs
    if N <= 0
        error('N must be a positive integer.');
    end
    
    % Compute Kaiser window
    n = 0:N-1; % Indices
    beta = alpha * sqrt(1 - ((2 * n / (N-1)) - 1).^2); % Argument for the Bessel function
    I0_alpha = besseli(0, alpha); % Denominator term (constant for the entire window)
    I0_beta = besseli(0, beta);   % Numerator term (varies with n)
    
    % Calculate Kaiser window
    w = I0_beta / I0_alpha;
end
