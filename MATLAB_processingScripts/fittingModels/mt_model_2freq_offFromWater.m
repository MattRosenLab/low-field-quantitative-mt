function out1=mt_model_2freq_offFromWater(t2a,ra,rb,mb0,r,t2b,sigtau,delta_MT,delta,w1,Txoff,lineshape)

% function out1=mtmodel(t2a,ra,rb,mb0,r,t2b,delta,w1)
% Note: may need to fix  t2a , ra, and rb to get good fit;
% Scott Swanson August 13 2024 (c)
% Modified by DK on 9/24/25 for when two offsets symmetric about the MT
% line are used

% DK notes: 
%   Most of these equations are found in Henkelman et al, MRM 1993
%   delta is the ±offset FROM A SPECIFIED OFFSET
%   delta and delta_MT should be in Hz
%   w1 is the nutation frequency PER OFFSET! So if the RF WITHOUT cosine
%       modulation were to have nutation 2*w1, each offset would have half
%       of that, i.e. w1 <-- THIS would be your input to this function!
%   w1 should be in rad/s
%   sigtau [unitless] is for Kubo-Tomita lineshape
%   Txoff [Hz] is the built-in offset about which delta is symmetric
%   lineshape options: {'lorentzian','gaussian','superlorentzian','kubo-tomita'}


if nargin<10 %shift over variables if sigtau not included
    w1=delta;
    delta=delta_MT;
    delta_MT=sigtau;
    sigtau=1;
end

if ~exist('lineshape','var')
    warning('Lineshape function not specified! Setting to Gaussian by default...')
    lineshape='gaussian';
end

t2b=t2b*1e-6;   %If t2b is too small (i.e. on order of 1e-5), fitting Jacobian may be singular!

% Fixing values for white matter fitting
% t2a=0.065;
% ra=0.4;
% rb=1;

delta_pos=Txoff+delta;
delta_neg=Txoff-delta;

rrfa1=t2a.*w1.^2./(1+(t2a.*2*pi*delta_pos).^2);  
rrfa2=t2a.*w1.^2./(1+(t2a.*2*pi*delta_neg).^2);    %eq10, Henkelman paper
                                            %NOTE: typo in paper! Missing
                                            %square of parentheses term

% Switch calculation of parameters based upon variable lineshape
switch lineshape
    case 'lorentzian'   % Lorentzian semisolid lineshape

        rrfb1=2*t2b.*w1.^2./(1+(t2b.*2*pi*(delta-delta_MT)).^2);
        rrfb2=2*t2b.*w1.^2./(1+(t2b.*2*pi*(-delta-delta_MT)).^2);

    case 'gaussian'     % Gaussian semisolid lineshape -- DK note: makes curve "lumpier"

        rrfb1=t2b*sqrt(pi/2)*w1.^2.*exp(-((2*pi*(delta_pos-delta_MT)*t2b).^2)/2);
        rrfb2=t2b*sqrt(pi/2)*w1.^2.*exp(-((2*pi*(delta_neg-delta_MT)*t2b).^2)/2);

    case 'superlorentzian' % Super-Lorentzian semisolid lineshape

        rrfb1=RF_superlorentzian(t2b,delta_MT,w1,delta);    % DK note: set chem shift to 0 for ULF!
        rrfb2=RF_superlorentzian(t2b,delta_MT,w1,-delta);    % DK note: set chem shift to 0 for ULF!

    case 'kubo-tomita'  % Kubo-Tomita semisolid lineshape

        rrfb1=RF_KuboTomita(t2b,sigtau,delta_MT,w1,delta_pos);    % DK note: set chem shift to 0 for ULF!
        rrfb2=RF_KuboTomita(t2b,sigtau,delta_MT,w1,delta_neg);    % DK note: set chem shift to 0 for ULF!

end                                           
                

%These three below come from eq9, Henkelman paper by factoring the
%numerator and multiplying out the denominator
numerator = r.*(ra+mb0*rb)+ra.*(rb+rrfb1+rrfb2); 
denom1 = (ra+rrfa1+rrfa2).*(rb+rrfb1+rrfb2);
denom2 = r.*(ra+rrfa1+rrfa2+mb0.*(rb+rrfb1+rrfb2));   

mtspec=numerator./(denom1 + denom2);

% % Equations from Scott Swanson; from Provotorov eq'ns I think
% % (NOTE: I think this only works if w1 is sqrt(2) larger than it would be
% % with single-offset irradiation)
% num2 = r .* (ra + mb0 .* rb) + ra .* (rb + rrfb);
% den2 = (ra + rrfa) .* (rb + rrfb) + r .* (ra + rrfa + mb0 .* (rb + rrfb));
% mtspec = num2 ./ den2;

% plot(numerator,'linewidth',2)
% hold on
% % plot(denom1,'linewidth',2)
% % plot(denom2,'linewidth',2)
% % plot(denom3,'linewidth',2)
% % legend('num','den1','den2','den3')
% hold off

%plot(delta,mtspec,'o')
out1 = mtspec;
