function projected_gamma = projection_cone2(gamma, r, recovStruct, varargin)
%projection_cone projects a vector onto a second-order cone
%constraint as well as a support constraint captured by the indx variable.
%->---
%input
%->---
%gamma:             array(3N,n_f) - molecular parameter estimates
%N:             scalar        - number of grid points
%r:             scalar        - constraint parameter associated with the secon-order cone
%n_f:           scalar        - number of frames
%optinal input parameters:
%indx:          logical array(3N,n_f) - indices of grid points with no
%molecules associated with them
%---->-
%output
%---->-
%projected_gamma:  array(3N,n_f) - projected molecular parameter estimates

numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:projectiong2:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

N = recovStruct.n_grid_p;

% set defaults for optional inputs
optargs = {[], [], []};
optargs(1:numvarargs) = varargin;

indx_xx = optargs{1};

%rearrange molecular parameters
%------------------------------
xx = gamma(1:N, :);
xxdx = gamma(1*N+1:2*N, :);
xxdy = gamma(2*N+1:3*N, :);

% x1=gamma(1:N,:);
% x2=gamma(N+1:2*N,:);
% x3=gamma(2*N+1:3*N,:);
% x23=sqrt(x2.^2+x3.^2);

xxjoint = sqrt(xxdx.^2+xxdy.^2);

%compute the conditions of the projection operator
%-------------------------------------------------

%XX
cxx = (xx + r * xxjoint) / (1 + r^2);
cxx1 = xxjoint <= xx * r;
cxx2 = xxjoint <= -xx / r;
cxx3 = xxjoint > abs(xx*r);
cxx_t1 = (1 - cxx2) .* cxx1;
cxx_t2 = (1 - cxx2) .* cxx3 .* cxx;

projected_gammaxx_t = cxx_t1 .* xx + cxx_t2;
projected_gammaxx_tx = cxx_t1 .* xxdx + cxx_t2 .* r .* (xxdx) ./ (eps + xxjoint);
projected_gammaxx_ty = cxx_t1 .* xxdy + cxx_t2 .* r .* xxdy ./ (eps + xxjoint);
projected_gammaxx_t = projected_gammaxx_t .* (projected_gammaxx_t > 0);


% c=(x1+r*x23)/(1+r^2);
% c1=x23<= x1*r;
% c2=x23<= -x1/r;
% c3=x23> abs(x1*r);
% c_t1=(1-c2).*c1;
% c_t2=(1-c2).*c3.*c;

% projected_gamma_t=c_t1.*x1+c_t2;
% projected_gamma_tx=c_t1.*x2+c_t2.*r.*(x2)./(eps+x23);
% projected_gamma_ty=c_t1.*x3+c_t2.*r.*x3./(eps+x23);
% projected_gamma_t=projected_gamma_t.*(projected_gamma_t>0);

%apply support constraints. Set indices in indx variable to zero
%---------------------------------------------------------------
% projected_gamma_t(indx)=0;
% projected_gamma_tx(indx)=0;
% projected_gamma_ty(indx)=0;

projected_gammaxx_t(indx_xx) = 0;
projected_gammaxx_tx(indx_xx) = 0;
projected_gammaxx_ty(indx_xx) = 0;


projected_gamma(1:N, :) = projected_gammaxx_t;
projected_gamma(1*N+1:2*N, :) = projected_gammaxx_tx;
projected_gamma(2*N+1:3*N, :) = projected_gammaxx_ty;
