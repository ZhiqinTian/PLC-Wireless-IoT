function [z0, gama] = line_parameter(f, type)
% calculate the paramters of the backbone and branch respectively, which
% include the inter-coupling effect.

% constant defination
mu_c = 12.5664e-7; % permeability of the conductor(brone); H/m mu0=4*pi*10e-7
sigma_c = 1.57e7; % conductivity of the conductor(brone); ohm/m
epslon0 = 8.85e-12; % permittivity of the dielectric material(air) F/m
% epslonr = 3.6;     % relative dielectic constant  (PVC)
epslonr = -3.3e-9.*f+3.7;
epslon = epslon0*epslonr;
sigma_d = 1e-4; % conductivity of the dielectirc material (aire + pvc); s/m
tdelta = -5.7e-10.*f+0.075;
switch type
    case 1     % for cable type of NYY 3*1.5
        r = 0.8921e-3;      % radius of the conductor; m
        D = 2*r+1e-3;       % distance between conductors; m
    case 2     % for cable type of NYY 3*2.5
        r = 0.8921e-3;      % radius of the conductor; m
        D = 2*r+1e-3;       % distance between conductors; m
    case 3     % for cable type of NYY 3****
        r = 0.8921e-3;      % radius of the conductor; m
        D = 2*r+1e-3;       % distance between conductors; m
    otherwise
        error('wrong cable type')
end

% calculate the distributed parameters-------resistance R
delta = sqrt(1./pi./f./mu_c./sigma_c);
R_solid = 1./pi./r./delta./sigma_c;
% XR = (acos((r_wire-delta)./r_wire).*r_wire^2-(r_wire-delta).*sqrt(r_wire^2-(r_wire-delta).^2))./(2.*r_wire.*delta);
R = (D./2./r)./sqrt((D./2./r).^2-1).*R_solid;       % ohm/m
% R = R_solid;

% calculate the inductance L
Ls = R./2./pi./f; % self-inductance of single conductor
Lm = mu_c/pi*acosh(D/2/r); % mutual inductance
L = Ls+Lm;

% calculate the capacitance C
C = pi*epslon/acosh(D/2/r);

% calculate the conductance G
% G = sigma_d*C/epslon;
G = 2.*pi.*f.*C.*tdelta;

% calculate the characteristic impedance and transmission constant
temp_a = complex(R, 2.*pi.*f.*L);
temp_b = complex(G, 2.*pi.*f.*C);
zc = sqrt(temp_a./temp_b);
z0 = abs(zc);       % consinder the characteristic impedance as real resistor
% z0 = sqrt(L./C) 

alfa = R./2./z0+G.*z0./2;
beta = 2.*pi.*f.*sqrt(L.*C);
gama = alfa+beta.*1i;
gama = sqrt(temp_a.*temp_b);

% c=1.5e8;%3e8./sqrt(-3.3e-9.*f+3.7);
% % beta=2.*pi.*f./c;
% beta = f.*0.32e-7;
% % c=1.5e8
% a0=5.8e-3;
% a1=6.88e-10;
% alfa=a1.*f+a0;
% gama=alfa+beta.*i;
% z0=50;

