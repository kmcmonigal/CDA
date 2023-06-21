function [EOF_X, EOF_Y, pcs_X, pcs_Y, CDAs, CDA_ts_X, CDA_ts_Y] = cda(X, Y)
% Based upon LFCA code by Dr. Robb Jnglin Wills (https://github.com/rcjwills/lfca)

%% CDA  Truncated Covariant Discriminant Analysis
%     performs covariance discriminant analysis on the data in 
%     matrices X and Y based on a the ratio of variance(X)/variance(Y).

%% INPUT
%     X and Y are 2D data matrices with time variations along the first dimension
%     and spatial variations along the second  dimension
%     The number of EOFs is dictated by the length of dimension 2. Practically
%     speaking, this is likely the number of EOFi calculated and projected onto X/Y
%

%% OUTPUT
%     EOF_X are the EOFs of time series X. pcs_X is the associated PC time series
%     EOF_Y are the EOFs of time series Y. pcs_Y is the associated PC time series
%     pcx_X is the percentage of variance explained by each EOF of X
%     pcx_Y is the percentage of variance explained by each EOF of Y
%     CDAs are the CDA patterns 
%     CDA_ts_X is the associated pattern time series for data X
%     CDA_ts_Y is the associated pattern time series for data Y
%     CDA patterns and timeseries are ordered such that the first pattern/timeseries maximizes variance(X)/variance(Y) and the last pattern/timeseries minimizes variance(X)/varianc(Y)
%     r_eofs is 

%% NOTES

%      scaleX and scaleY currently set as all ones, could be implemented to weight by grid cell size

%  narginchk(2)          % check number of input arguments 
  if ndims(X) ~= 2,  error('Data matrix X must be 2-D.'); end % check dimensions
  if ndims(Y) ~= 2,  error('Data matrix Y must be 2-D.'); end
  
  [n,p]         = size(X);
  N=p; % take as many EOFs as we had EOFi

  % center data 
  if any(any(isnan(X)))               % there are missing values in x
    Xm  = nanmean(X);
  else                                % no missing values
    Xm  = mean(X);
  end
  X    = X - repmat(Xm, n, 1);  
      
  %% compute covariance matrix
  Covtot               = cov(X);
  cov_X=Covtot;
  %% scale vector (e.g. square root of normalized grid-cell area)
    scaleX       = ones(1,p);
    Xs          = X;
  clear X
  
  %% eigendecomposition of covariance matrix
  Covtot      = repmat(scaleX',1,p) .* Covtot .* repmat(scaleX,p,1);
  [pcvec,evl,rest] = peigs(Covtot, min(n-1, p));
  trCovtot    = trace(Covtot);
  
  % percent of total sample variation accounted for by each EOF
  pvar          = evl./trCovtot .* 100;
  % principal component time series
  pcs_X           = Xs*pcvec;
  % return EOFs in original scaling as patterns (row vectors)
  EOF_X           = pcvec' ./ repmat(scaleX,rest,1);

  % whitening - this S comes in to the equation later - need it  
  f             = sqrt(evl(1:N));
  
  % get transformation matrices that transform original variables to
  % whitened variables and back
  S		= pcvec(:, 1:N) * diag(1./f); 
  Sadj	        = diag(f) * pcvec(:, 1:N)';

  %%%%%%%%%%%%%%%%% same calculation but on Y

  [n,p]         = size(Y);

  % center data 
  if any(any(isnan(Y)))               % there are missing values in Y
    Ym  = nanmean(Y);
  else                                % no missing values
    Ym  = mean(Y);
  end
  Y    = Y - repmat(Ym, n, 1);  
      
  Covtot               = cov(Y);
  cov_Y=Covtot;

  if any(size(Covtot) ~= [p, p])
    error('Covariance matrix must have same dimension as data.')
  end
  
  %% scale vector (e.g. square root of normalized grid-cell area)
    scaleY       = ones(1,p);
    Ys          = Y;
  %clear Y
  
  %% eigendecomposition of covariance matrix
  Covtot      = repmat(scaleY',1,p) .* Covtot .* repmat(scaleY,p,1);
  [pcvec_Y,evl,rest] = peigs(Covtot, min(n-1, p));
  trCovtot    = trace(Covtot);
  
  % percent of total sample variation accounted for by each EOF
  pvar_Y          = evl./trCovtot .* 100;
  % principal component time series
  pcs_Y           = Ys*pcvec_Y;
  % return EOFs in original scaling as patterns (row vectors)
  EOF_Y           = pcvec_Y' ./ repmat(scaleY,rest,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Now we minimize/maximize ratio of X/Y
  %% whiten variables
  Z		= Y * S;

  %  covariance matrix of whitened variables
  % (i.e. covariance matrix of filtered and whitened principal components)
  Gamma = cov(Z);

  %% SVD of covariance matrix (such that r are eigenvalues and V are eigenvectors)
  [~, r, V]	= csvd(Gamma);

  %% weight vectors (canonical vectors) and patterns (CDAs) in original scaling
  weights	= repmat(scaleX', 1, N) .* (S * V);       % weights are columns
  CDAs	= (V' * Sadj) ./ repmat(scaleX, N, 1);    % patterns are rows 
  
%% get associated timeseries
  
if nargin > 3
    Xs = Xs./repmat(scaleX,n,1);
end

CDA_ts_X = Xs * weights;

if nargin > 3
    Ys = Ys./repmat(scaleY,n,1);
end

CDA_ts_Y = Ys * weights;

% get the discriminant ratios
X_var_eofs = diag(pcvec'*cov_X*pcvec)./diag(pcvec'*pcvec);
Y_var_eofs = diag(pcvec_Y'*cov_Y*pcvec_Y)./diag(pcvec_Y'*pcvec_Y);

