function Hf_a_dB=transfer_calc(f,N,units_structure)
% input arguments
tn = 1;
z_rx = 50;  % port imepdance of the receiver
% end of input argument
for n = 1:N     % for each T unit
    %bone structure input
    bone1_len = units_structure{n}.bone1_len;
    bone1_type = units_structure{n}.bone1_type;
    bone2_len = units_structure{n}.bone2_len;
    bone2_type = units_structure{n}.bone2_type;
    
    %branch structure input
    level_cable_num = units_structure{n}.level_cable_num;
    level_load_type = units_structure{n}.level_load_type;
    level_cable_type = units_structure{n}.level_cable_type;
    level_len = units_structure{n}.level_len;
    next_connect = units_structure{n}.next_connect;
    level_node_load = units_structure{n}.level_node_load;
    root_load = units_structure{n}.root_load;
    bflag = units_structure{n}.bflag;
    for t = 1:tn  % t is the time index
%         zin_branch = br_impedance(tm(t), f, level_cable_num, level_load_type, level_cable_type, level_len, next_connect, bflag)
        zin_branch = br_impedance( f, level_cable_num, level_load_type, level_cable_type, level_len, next_connect,level_node_load,root_load, bflag);
        % get the transfer function
        if n == 1      % for the last T unit which is the receiver
            zl = z_rx;
        else           % for other units
            zl = squeeze(z(t,n-1,:))';
        end
        [Hf, Zin] = unit_transfer(f, bone1_len, bone1_type, bone2_len, bone2_type, zin_branch, zl);
        Hv(t,n,:) = Hf;
        z(t,n,:) = Zin;     % get the input impedance 
    end
end
% get the total transfer function
for m = 1:tn
    total = 1;
    for k = 1:N
        temp = squeeze(Hv(m,k,:))';
        Hf_total(m,:) = temp.*total;
%         Hf_total(m,k,:) = temp.*total
        total = squeeze(Hf_total(m,:));
    end
end

tx = 0:1/100e6:length(f)/100e6;
for j = 1:tn
    % get the amplitude response in dB
    Hf_a(j,:) = abs(Hf_total(j,:));
    Hf_a_dB(j,:) = 10.*log10(Hf_a(j,:));

    % get the phase response in rad
    %Hf_p_rad(j,:) = angle(Hf_total(j,:)');
    
    % get the windowed impulse response
    %w = hamming(length(f));
    %wh = [Hf_total(j,:)];
    %hf(j,:) = ifft(wh'.*w);

    % get the path loss from the total transfer function
    %zpass = squeeze(z(j,N,:))';
    %PL = pathloss(zpass,z_rx,Hf_a(j,:));
    %PLdb(j,:) = 10.*log10(PL);
end
% plot the amplitude response
% figure
% subplot(2,1,1);
% plot(f,Hf_a_dB,'r');
% grid on
% xlabel('Frequency(Hz)');
% ylabel('Frequency Response(dB)');
% 
% % plot the impulse response
% subplot(2,1,2);
% plot(tx(1:100),real(hf(1:100)),'r');
% grid on
% xlabel('Time(s)');
% ylabel('Impulse Response');