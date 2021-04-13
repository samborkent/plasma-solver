function curveGauss = fitGaussian( x, y, yThreshold )
%FITGAUSSIAN Output a normalized Gaussian curve
%   Use for plotting a Gaussian curve over the non-collided particles

% Default smoothing value to multiply with the Gaussian curve width
smoothValue = 10;

% Remove first and last points to prevent the curve from sloping up at the
%   limits
xTemp = x(2 : end-1);
yTemp = y(2 : end-1);

% Only data points above threshold
xTemp = xTemp(yTemp > yThreshold);
yTemp = yTemp(yTemp > yThreshold);

% Fit a Gaussian curve
fitGauss = fit( xTemp', yTemp', 'gauss2' );

% Evaluate the Gaussian fit
curveGauss = feval( fitGauss, xTemp );

% Interpolate the data to match the original x-axis
curveGauss = interp1( xTemp, curveGauss, x );

% Smooth data based on width of Gaussian curve
curveGauss = smooth( x, curveGauss, smoothValue*fitGauss.c1 );

% Normalize
curveGauss = curveGauss ./ sum(curveGauss);

end

