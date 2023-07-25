function rhoa = dcfwdf(res,thk,ab)
% DCFWDF Fast Schlumberger Resistivity Sounding
%
% REMARKS
%   The underlying problem is the solution of the 1D Laplace equation using
%   the Wait-Algorithm.
%   As the evaluation point lies at the surface, no additional 
%   propagation from the surface into the earth is required.
%   As the desired observation is the specific electical resistivity, the
%   electrical fields are calculated, rather then the potetial. 
%   Therefore, the filter coefficitents represent the approximated
%   first order Bessel functions, resulting from the derivative of the
%   zeroth order Bessel functions (required for the potential).
%   In priciple the algorithm below represents the following idea,
%   the solution of 1D Laplace is something like:
%       Int_inf Bessel_1(a) * wait_kernel da
%   this can be solved approximately by quadrature or filterapproach
%   providing:
%       wait_kernel * (Sum_j Bessel_1(a_j) * w_j da)
%   as summation is always the same thing, this can be summarized as 
%   tabalized constants (so called filter coeffitients)
%       wait_kernel * fc
%   The wait_kernel results from propagating the behavior at top of the 
%   infinite layer across the layer boundaries exploiting the continuity
%   conditions of the potential.
%
%   Usage: RHOA=DCFWDF(RES,THK,AB)
%   Input arrays:
%   RES: vector of layer resistivities
%   THK: vector of LENGTH(RES)-1 layer thicknesses
%   AB:  vector of electrode spacings
%   Output array:
%   RHOA: vector of apparent resistivities of LENGTH(AB)
nd = length(ab);
rhoa=zeros(nd,1);
nl = length(res);
fc=[-.22247786e-4 .51184989e-4 -.66575186e-4 .86592875e-4 ...
        -.11262944e-3 .14649463e-3 -.19054233e-3 .24783420e-3 ...
        -.32235248e-3 .41927675e-3 -.54534402e-3 .70931693e-3 ...
        -.92259288e-3 .11999962e-2 -.15608086e-2 .20301093e-2 ...
        -.26405183e-2 .34344639e-2 -.44671314e-2 .58102992e-2 ...
        -.75573279e-2 .98296496e-2 -.12785208e-1 .16629439e-1 ...
        -.21629544e-1 .28133073e-1 -.36592072e-1 .47594515e-1 ...
        -.61905179e-1 .80518827e-1 -.10472943e+0 .13622036e+0 ...
        -.17718202e+0 .23046613e+0 -.29978969e+0 .39001009e+0 ...
        -.50751079e+0 .66077997e+0 -.86135847e+0 .11254667e+1 ...
        -.14762271e+1 .19413057e+1 -.25117825e+1 .29397638e+1 ...
        -.22862253e+1 -.71362115e+0  .41491251e+1 -.23169602e+1 ...
        -.16867419e+1 -.32170199e+0  .68963453e+0 .69150845e+0 ...
        .54204064e+0 .32222510e+0  .19033795e+0 .99724447e-1 ...
        .54063095e-1 .27109364e-1 .14239291e-1 .70240599e-2 ...
        .36435998e-2 .17863940e-2 .92183691e-3 .45100602e-3 ...
        .23218367e-3 .11353351e-3 .58376418e-4 .28547132e-4 ...
        .14666868e-4 .14501929e-4];
nc=length(fc);
f=10^0.1;
y1=12.664218;
x = zeros(nd,nc);
x = exp(y1) ./ (ab(:) * f.^[0 : nc-1]);
tt = res(nl)*ones(nd,nc);
for k = nl-1 : -1 : 1
	r = (res(k) - tt) ./ (res(k) + tt) .* exp(-2*thk(k)*x);
	tt = res(k) * (1 - r) ./ (1 + r);
end
rhoa = tt * fc(:);