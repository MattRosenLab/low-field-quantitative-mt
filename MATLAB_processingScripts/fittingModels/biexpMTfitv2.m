function out=biexpMTfitv2(R1a,R1b,Mb0,R,Ma_start,T2b,Tpulse,t,lineshape)
% biexpMTfit.m
% Computes a biexponential recovery curve for [selective inversion]-
% recovery data from a sample containing a semisolid pool 
% Uses the equations from Gochberg & Gore, MRM 2003, as well as from
% Gochberg et al, MRM 1999 to estimate the semisolid Mb saturation caused
% by the inversion pulse
%
% INPUTS:
%   R1a         -   Spin-lattice relaxation rate of water spins [s^-1] 
%   R1b         -   Spin-lattice relaxation rate of semisolid spins [s^-1]
%   Mb0         -   Proton volume fraction of semisolid spins (= Msemisolid/Mwater)
%   R           -   Magnetization transfer rate from semisolid to water spins [s^-1]
%   Ma_start    -   Normalized inverted magnetization of water at time t=0
%   T2b         -   Spin-spin relaxation rate of semisolid spins [s^-1]
%   Tpulse      -   Duration of inversion pulse used [s] - used to
%                   calculate the pulse nutation frequency
%   t           -   Vector containing the recovery times between the inversion 
%                   and excitation pulses [s]
%   lineshape   -   String indicating the semisolid lineshape. Options are:
%                   {'lorentzian','gaussian','superlorentzian','kubo-tomita'}
%
% MAJOR CHANGE: Instead of using a value of Mb_start, we will estimate it
% using Eq. A15 in Gochberg et al, MRM 1999:
%       Mb_start = exp(-w1_pulse^2 * t_pulse * T2b)
%

T2b=T2b*1e-6;   %If t2b is too small (i.e. on order of 1e-5), fitting Jacobian may be singular!

% Set forward and backward exchange rates
kab=Mb0*R;
kba=R;

% Calculate w1 of pulse (in rad/s) from duration
w1pulse=pi/Tpulse;

% Calculate the exponential rates R1p and R1n
R1p=0.5*(R1a + R1b + kab + kba + sqrt((R1a - R1b + kab - kba)^2 ...
    + 4*kab*kba));
R1n=0.5*(R1a + R1b + kab + kba - sqrt((R1a - R1b + kab - kba)^2 ...
    + 4*kab*kba));

% Estimate Mb_start, the partially saturated semisolid magnetization due to
% the selective inversion pulse 
switch lineshape %NOTE: there are probably better expressions for 
    % superlorentzian and kubo-tomita
    case {'lorentzian', 'superlorentzian'}
        adjExpFac=1;            %for Lorentzian semisolid pool
    case {'gaussian', 'kubo-tomita'}
        adjExpFac=sqrt(pi/2);   %for Gaussian semisolid pool
end
Mb_start=exp(-w1pulse^2 * Tpulse * T2b * adjExpFac);

% Calculate the exponential weightings Bp and Bn
Bp=((Ma_start - 1)*(R1a - R1n) + (Ma_start - Mb_start)*kab)/(R1p - R1n);
Bn=-((Ma_start - 1)*(R1a - R1p) + (Ma_start - Mb_start)*kab)/(R1p - R1n);

% Calculate the full curve
out=Bp*exp(-R1p*t) + Bn*exp(-R1n*t) + 1;
end