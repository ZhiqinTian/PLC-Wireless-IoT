function zin_branch = br_impedance( f, level_cable_num, level_load_type, level_cable_type, level_len, next_connect, level_node_load,root_load,bflag)
[zin_branch] = root_load;
if bflag == 0
    level_num = length(level_cable_num);
    max_level_line_num = max(level_cable_num);
    %level_connect = zeros(level_num,max_value);
    zin = 1e8.*ones(level_num, max_level_line_num,length(f));
    
    for i = 1:level_num
        for j = 1:max_level_line_num
            zin(i,j,:) = level_node_load(i,j);
        end
    end
 
    for l = level_num:-1:2   % from top level to level 1
        for b = 1:level_cable_num(l)
            if level_load_type(l,b) > 0
                z_load = load_select(f,level_load_type(l,b));
            elseif level_load_type(l,b) == -1
                z_load = zin(l,b,1);
            end
            [zc, gama] = line_parameter(f,level_cable_type(l,b));
            z = zin_compute(zc,gama,level_len(l,b),z_load);
            zp = squeeze(zin(l-1,next_connect(l,b),:));
            zin(l-1,next_connect(l,b),:) = zin_parallel(zp',z);
        end
    end
    
    for b1 = 1:level_cable_num(1)
        if level_load_type(1,b1) > 0
            z_load1 = load_select(f,level_load_type(1,b1));
        elseif level_load_type(1,b1) == -1
            z_load1 = squeeze(zin(1,b1,:))';
        end
        [zc, gama] = line_parameter(f,level_cable_type(1,b1));
        z1 = zin_compute(zc,gama,level_len(1,b1),z_load1);
        zin_branch = zin_parallel(zin_branch,z1);
    end
end