function a = dolph_chebyshev(N, R)
    % DOLPH_CHEBYSHEV Generate Dolph-Chebyshev window coefficients
    %
    % Inputs:
    %   N - Number of coefficients (window length)
    %   R - Desired sidelobe level in dB
    %
    % Output:
    %   a - Coefficients of the Dolph-Chebyshev window

    % Number of zeros
    N1 = N - 1;

    % Convert sidelobe level to absolute units
    Ra = db2mag(R); % Sidelobe level in absolute units

    % Compute the scaling factor
    x0 = cosh(acosh(Ra) / N1);

    % Compute zeros of the Chebyshev polynomial in cosine space
    i = 1:N1;
    xi = cos(pi * (i - 0.5) / N1);

    % Transform zeros into psi-space and z-space
    psi = 2 * acos(xi / x0);  % Zeros in psi-space
    zi = exp(1j * psi);       % Zeros in z-domain (complex exponential)

    % Convert zeros to polynomial coefficients
    a = real(poly2(zi));      % Real coefficients of the polynomial
end
