function zin = zin_compute(zc,gama,l,zl)
% this function computes the input impedance of the transmission line
% connecting to a load
zin = zc.*(zl+zc.*tanh(gama.*l))./(zc+zl.*tanh(gama.*l));