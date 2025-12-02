function [best_params, resnorm, yfit] = fit_travis_prof_from_data(xdata, ydata, varargin)
%FIT_TRAVIS_PROF_FROM_DATA  Fit the Travis‐profile model to experimental data.
%
%   [params, resnorm] = FIT_TRAVIS_PROF_FROM_DATA(xdata, ydata) fits the
%   model
%
%       N(p,x) = p1 - p2 + (1-p1+p2)*(1 - x.^p3).^p4 + ...
%                p2*(1 - exp(-x.^2./p5.^2))
%
%   to the supplied vectors xdata and ydata using a non‑linear least‑squares
%   algorithm (lsqcurvefit).  The function returns the best‑fit parameter
%   vector `params` (1×5) and the residual norm `resnorm`.
%
%   [...] = FIT_TRAVIS_PROF_FROM_DATA(..., 'InitialGuess', guess) lets you
%   specify a custom 5‑element initial guess.  The default is
%       guess = [0.0, 6, 10, 0.20, 0.25];
%
%   [...] = FIT_TRAVIS_PROF_FROM_DATA(..., 'Options', opts) lets you pass a
%   custom options structure created with optimoptions('lsqcurvefit',...).
%   If omitted, the function uses
%       opts = optimoptions('lsqcurvefit','Display','iter');
%
%   [...] = FIT_TRAVIS_PROF_FROM_DATA(..., 'Plot', true) will produce a
%   figure showing the raw data and the fitted curve.  The default is false.
%
%   [...] = FIT_TRAVIS_PROF_FROM_DATA(..., 'ReturnFit', true) also returns
%   the fitted y‑values evaluated on a dense grid (100 points) in the
%   variable yfit.  If false (default) yfit is returned as an empty array.
%
%   Example
%   -------
%       % Suppose you have a structure array `prof_data` and you want to fit
%       % the profile of the i‑th element:
%       i = 3;
%       x = prof_data(i).ne.reff ./ prof_data(i).ne.reff(end);
%       y = prof_data(i).ne.value ./ prof_data(i).ne.value(1);
%
%       [p, rss, yfit] = fit_travis_prof_from_data(x, y, ...
%                               'Plot', true, 'ReturnFit', true);
%
%   See also lsqcurvefit, optimoptions.
%
%   Author:  David Kulla
%   Date:    2025‑11‑14
%   Version: 1.0

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
p.CaseSensitive = false;
p.FunctionName   = mfilename;

% Required arguments
addRequired(p, 'xdata', @(v) isnumeric(v) && isvector(v));
addRequired(p, 'ydata', @(v) isnumeric(v) && isvector(v) && numel(v)==numel(xdata));

% Optional name‑value pairs
defaultGuess   = [0.0, 6, 10, 0.20, 0.25];
defaultOpts    = optimoptions('lsqcurvefit','Display','iter');
defaultPlot    = false;
defaultReturnFit = false;

addParameter(p, 'InitialGuess', defaultGuess, @(v) isnumeric(v) && numel(v)==5);
addParameter(p, 'Options',      defaultOpts,  @(v) isstruct(v));
addParameter(p, 'Plot',         defaultPlot,  @(v) islogical(v) || isnumeric(v));
addParameter(p, 'ReturnFit',    defaultReturnFit, @(v) islogical(v) || isnumeric(v));

parse(p, xdata, ydata, varargin{:});

xdata       = p.Results.xdata(:);   % force column vectors
ydata       = p.Results.ydata(:);
initialGuess = p.Results.InitialGuess;
options      = p.Results.Options;
doPlot       = logical(p.Results.Plot);
doReturnFit  = logical(p.Results.ReturnFit);

% -------------------------------------------------------------------------
% Model definition (anonymous function)
% -------------------------------------------------------------------------
% params = [g, p, q, hole_width, sigma]
N = @(params, x) params(1) - params(2) + ...
                 (1 - params(1) + params(2)) .* (1 - x.^params(3)).^params(4) + ...
                 params(2) .* (1 - exp(-x.^2 ./ params(5).^2));

% -------------------------------------------------------------------------
% Perform the fit
% -------------------------------------------------------------------------
% No bounds are supplied (empty matrices) – you can modify this if you need
% constraints.
lb = [];   % lower bounds
ub = [];   % upper bounds

[best_params, resnorm] = lsqcurvefit(N, initialGuess, xdata, ydata, lb, ub, options);

% -------------------------------------------------------------------------
% Optional return of a smooth fitted curve
% -------------------------------------------------------------------------
if doReturnFit
    xfit = linspace(min(xdata), max(xdata), 100);
    yfit = N(best_params, xfit);
else
    yfit = [];   % keep output consistent
end

% -------------------------------------------------------------------------
% Optional plotting
% -------------------------------------------------------------------------
if doPlot
    figure('Name','Travis Profile Fit','NumberTitle','off');
    hold on;
    plot(xdata, ydata, 'ro', 'MarkerFaceColor','r', 'DisplayName','Experimental Data');
    if isempty(yfit)   % compute on‑the‑fly if not already returned
        xfit = linspace(min(xdata), max(xdata), 200);
        yfit = N(best_params, xfit);
    end
    plot(xfit, yfit, 'b-', 'LineWidth',1.5, 'DisplayName','Fitted Curve');
    legend('show','Location','best');
    xlabel('x');
    ylabel('N(x)');
    title('Curve Fitting Result');
    grid on;
    hold off;
end

% -------------------------------------------------------------------------
% Display results (optional – can be suppressed by redirecting output)
% -------------------------------------------------------------------------
fprintf('Best‑fitted parameters:\n');
disp(best_params);
fprintf('Residual norm (sum of squares): %g\n', resnorm);

end   % <--- end of function