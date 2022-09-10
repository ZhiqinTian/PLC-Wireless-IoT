function [Hf, Zin] = unit_transfer(f, bone1_len, bone1_type, bone2_len, bone2_type, zin_branch, zl)
%  this function computes the transfer function for each unit according to
%  the ratio of voltage
%  input arguments:
%  f----frequency
%  bone1_len----length of the latter section of the backbone
%  bone2_len----length of the former section of the backbone
%  bone1_type----type of the latter section of the backbone
%  bone2_type----type of the former section of the backbone
%  zin_branch----total input impedance of the branch
%  zl----the impedance connected to the latter section of the backbone

%  output arguments:
%  Hf----the transfer function of each unit presented by the voltage ratio
%  Zin----the input impedance of each unit looking from the transmitting
%  port

% compute the characteristic impedance and trasmission constant of backbone
[zc1,gama1]=line_parameter(f,bone1_type);        
[zc2,gama2]=line_parameter(f,bone2_type);

zinl = zin_compute(zc1,gama1,bone1_len,zl); % input impedance of the latter backbone
zinb = zin_branch;      % input impedance of the branch

zp = zin_parallel(zinl,zinb);        

Zin = zin_compute(zc2,gama2,bone2_len,zp);

Tl = (zl-zc1)./(zl+zc1);
Tp = (zp-zc2)./(zp+zc2);

v1 = (1+Tl)./(exp(gama1.*bone1_len)+Tl.*exp(-gama1.*bone1_len));
v2 = (1+Tp)./(exp(gama2.*bone2_len)+Tp.*exp(-gama2.*bone2_len));
% v2 = exp(-gama1.*bb1_len).*(1+Tp)./(1+Tp.*exp(-2.*gama1.*bb1_len))

Hf = v1.*v2;
