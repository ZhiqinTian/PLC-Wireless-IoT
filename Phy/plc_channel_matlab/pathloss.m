function ploss = pathloss(z,z_rx,H)
% this function computes the path loss from transfer function

ploss = (H).^2.*real(z_rx)./real(z).*(abs(z).^2)./(abs(z_rx)).^2;