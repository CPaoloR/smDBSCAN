function [f,xbin,dbin] = hist_fd(X,varargin)
%written by
%C.P.Richter
%Division of Biophysics / Group J.Piehler
%University of Osnabrueck

%modified 21.03.2015
%modified 16.04.2015

ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'X',@(x)isvector(x))
addParamValue(ip,'verbose', false, @(x)islogical(x))
parse(ip,X,varargin{:});

verbose = ip.Results.verbose;

%%
N = numel(X);
dbin = 2*iqr(X)/nthroot(N,3); % freedman diaconis rule
nbin = ceil(range(X)/dbin); %make # bins discrete
[f,xbin] = hist(X,nbin);
f = f/sum(f)/dbin; % normalize integral to 1
dbin = xbin(2)-xbin(1);

%%
if verbose
    figure('color','w'); hold on
    bar(xbin,f,'hist')
    axis tight
    box on
end
end %fun