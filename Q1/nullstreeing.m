function w = nullstreeing(D, thetanull, thetabeam)
% Calculate the steering vector for null directions
wn = steervec(D,thetanull);

% Calculate the steering vectors for lookout directions
wd = steervec(D,thetabeam);

% Compute the response of desired steering at null direction
rn = wn'*wd/(wn'*wn);

% Sidelobe canceler - remove the response at null direction
w = wd-wn*rn;