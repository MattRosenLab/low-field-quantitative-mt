function out1=mt_model_w_dipolar(t2a,ra,rb,mb0,r,t2b,sigtau,delta_MT,td,delta,w1,lineshape)

% Note: may need to fix  t2a , ra, and rb to get good fit;
% Scott Swanson August 13 2024 (c)

% DK notes: 
%   These equations are from Morrison et al, J Magn Reson B 1995
%   delta should be in Hz
%   w1 should be in rad/s
%   sigtau [unitless] is for Kubo-Tomita lineshape
%   lineshape options: {'lorentzian','gaussian','superlorentzian','kubo-tomita'}

if nargin<11 %shift over variables if sigtau not included
    w1=delta;
    delta=td;
    td=delta_MT;
    delta_MT=sigtau;
    sigtau=1;
end

if ~exist('lineshape','var')
    warning('Lineshape function not specified! Setting to Gaussian by default...')
    lineshape='gaussian';
end

% Fixing values for white matter fitting
% t2a=0.065;
% ra=0.4;
% rb=1;

rrfa=t2a.*(w1.^2)./(1+(t2a.*2*pi*delta).^2);  %eq10, Henkelman paper
                                            %NOTE: typo in paper! Missing
                                            %square of parentheses term
% rrfa=((w1./(2*pi*delta)).^2)./t2a;              %in-line eq between eqs 2 and 3,
%                                                 %Morrison et al JMRB 1995,
%                                                 %assumes (2*pi*delta*T2a)^2 >> 1

% Switch calculation of parameters based upon variable lineshape
switch lineshape
    case 'lorentzian'   % Lorentzian semisolid lineshape

        rrfb=t2b.*w1.^2./(1+(t2b.*2*pi*(delta-delta_MT)).^2);
        d=1/sqrt(3)/t2b;    % value of D based upon Td and Gaussian lineshape 
                            % (I don't know what Lorentzian would be!)

    case 'gaussian'     % Gaussian semisolid lineshape -- DK note: makes curve "lumpier"

        rrfb=t2b*sqrt(pi/2)*w1.^2.*exp(-((2*pi*(delta-delta_MT)*t2b).^2)/2);
        d=1/sqrt(3)/t2b;    % value of D based upon Td and Gaussian lineshape

    case 'superlorentzian' % Super-Lorentzian semisolid lineshape

        rrfb=RF_superlorentzian(t2b,delta_MT,w1,delta);    % DK note: set chem shift to 0 for ULF!
        d=1/sqrt(15)/t2b; % value of D based upon Td and super-Lorentzian lineshape

    case 'kubo-tomita'  % Kubo-Tomita semisolid lineshape

        rrfb=RF_KuboTomita(t2b,sigtau,delta_MT,w1,delta);    % DK note: set chem shift to 0 for ULF!
        d=1/sqrt(3)/t2b;    % value of D based upon Td and Gaussian lineshape 
                            % (I don't know what Lorentzian would be!)

end


% Parameters S, T, and U, from Morrison et al paper
% (DK note: I rewrote to include rrfa when appropriate: 
% rrfa = (w1 / (2*pi*delta) )^2 * 1/t2a
S = r*mb0/ra*rb + rb + r;
T = S + rrfa./ra.*(rb + r);
U = r*mb0/ra + 1 + rrfa./ra;

%These three below come from eq9, Henkelman paper by factoring the
%numerator and multiplying out the denominator
numerator = S*(1 + ((2*pi*delta./d).^2).*rrfb*td) + rrfb; 
denom = T.*(1 + ((2*pi*delta./d).^2).*rrfb*td) + rrfb.*U;

mtspec=numerator./denom;

% plot(numerator,'linewidth',2)
% hold on
% % plot(denom1,'linewidth',2)
% % plot(denom2,'linewidth',2)
% % plot(denom3,'linewidth',2)
% % legend('num','den1','den2','den3')
% hold off

%plot(delta,mtspec,'o')
out1 = mtspec;
