% calculate EOFs, inverse eofs (EOFi)
% when calculating CDA, it is standard to project the EOFi from a common set (such as observations) onto both X and Y prior to calculating the CDA patterns

function [EOF,PC,FEXPVAR,EOFi]=eofi(data,neof,lat,lon)

%%% EOF calculates the first neof EOFs and EOFi of data matrix "data"

%%% INPUT
% data is a 2d or 3d matrix with time variations along the first dimension
% neof is the number of EOF and EOFi patterns to be calculated
% lat, lon are X by 1
%%% OUTPUT
% EOF are the first neof EOFs
% PC are the associated principal component timeseries
% FEXPVAR is the percent variance explained by each EOF
% EOFi are the inverse EOFs

% get the size of the data, put time as last dimension
ntime=size(data,1);
d=ndims(data);
if d==2
    ns=size(data,2);
    data=permute(data,[2 1]);
elseif d==3
    ns=size(data,2)*size(data,3);
    ni=size(data,2);
    nj=size(data,3);
    data=permute(data,[2 3 1]);
end

%weight by cosine of latitude
weight=repmat(abs(cosd(lat)),1,length(lon));
if 1
data_w=data.*weight.';
end

if 0 % could alter to not apply any weighting 
weight_ones=ones(size(weight));
weight=weight_ones;

data_w=data;
end

% if 3d, need to reshape
if d==3
    x=reshape(data_w,ns,ntime);
else
    x=data_w;
end

% calculate EOFs
x=x.';
x(isnan(x))=0;
[P,L,U]=svd(x,'econ'); % U are our EOFs
A=P*L; % PCs 
G=ctranspose(L)*L/(ntime-1);

% reshape and apply weighting 
if d==2
    Q=reshape(U,ns,ntime);
    S=reshape(L,ns,ntime); % for weighting
else
    Q=reshape(U,ni,nj,ntime);
end
if d==2
    EOF=zeros(ns,neof);
    EOFi=zeros(ns,neof);
else
    EOF=zeros(ni,nj,neof);
    EOFi=zeros(ni,nj,neof); 
end
PC=zeros(ntime,neof);

% calculate first neof patterns and pc timeseries
for i=1:neof
    if d==3
        EOF(:,:,i)=(Q(:,:,i)*G(i,i))./(weight.');
        EOFi(:,:,i)=(Q(:,:,i).*weight.')./(G(i,i));
    else
        EOF(:,i)=(Q(:,i)*G(i,i))./(weight.');
        EOFi(:,i)=(Q(:,i).*weight.')./(G(i,i));
    end
    PC(:,i)=A(:,i)/std(A(:,i)); 
    FEXPVAR(i)=G(i,i)/sum(sum(G)); % percent variance explained
end
EOF(EOF==0)=nan; % mask anything that is exactly equal to zero 

end

% weighting by cos(lat) -> think about
