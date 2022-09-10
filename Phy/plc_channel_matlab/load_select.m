function [zload] = load_select(f,type)
switch type
    case 1
        zload = time_variant_load(); % Note: the frequency point may not be the same as that in the initilization in Main.m
    case 2 
        zload = 1e8;
    case 3
        zload = 1e8;
%     case 4 
%         zload = 10.*cos(2.*pi.*50.*t);
%     case 5
%         zload = 1e8;
%     case 6 
%         zload = 50-76*i+10.*sin(2.*pi.*50.*t);
%     case 7
%         zload = 200-2.*pi.*f.*0.0000000009.*i+10.*sin(2.*pi.*50.*t)%130;
    otherwise
        error('no such load type')
end
