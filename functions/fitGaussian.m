function curveGauss = fitGaussian( x, y, yThreshold )
%FITGAUSSIAN Output a normalized Gaussian curve
%   Use for plotting a Gaussian curve over the non-collided particles

% Default smoothing value to multiply with the Gaussian curve width
smoothValue = 10;

% Remove first and last points to prevent the curve from sloping up at the
%   limits
x([1 end]) = 0;
y([1 end]) = 0;

% Only data points above threshold
xTemp = x(y > yThreshold);
yTemp = y(y > yThreshold);

% At least 6 points are required to fit a gaussian curve
if numel(y) > 6
    % Fit a Gaussian curve
    fitGauss = fit( xTemp.', yTemp.', 'gauss2' );

    % Evaluate the Gaussian fit
    curveGauss = feval( fitGauss, xTemp );

    % Interpolate the data to match the original x-axis
    curveGauss = interp1( xTemp, curveGauss, x );

    % Smooth data based on width of Gaussian curve
    curveGauss = smooth( x, curveGauss, smoothValue*fitGauss.c1 );

    % Remove data below the threshold
    curveGauss(curveGauss < yThreshold) = 0;

    % Normalize
    curveGauss = curveGauss ./ sum(curveGauss);
else
    curveGauss = y;
end

end

