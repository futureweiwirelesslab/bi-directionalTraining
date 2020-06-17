% Bi-directional Training for 5G
% An approach for interference probing and mitigation to achieve high spectrum efficiency in 5G Massive MIMO
% FUTUREWEI 2019, MIT license
%%
clear all;
close all;
warning('off','all');
% NOTEs: The terms BiT, Bi-directional Training, and interference probing are used interchangeably
% The terms UE, CPE, and MS are used interchangeably throughout; 
% the terms BS, TRP, sector, and gNB are used interchangeably throughout; 
% a cell site (numbered 1~7) has 3 sectors
%% variables

MAX_RANK_CPE = 2; %rank per ms

TXP_BS = 47;      % bs transmisison power in dbm
TXP_MS = 23;      % ms transmisison power in dbm
NF_BS = 3;        % bs noise figure in dB
NF_MS = 9;        % ms noise figure in dB

N_RB = 50;              % number of RBs in a carrier with 10MHz bandwidth
TONE_RB = 12;           % number of tones per RB
N_TONE = N_RB*TONE_RB;
SC_F = 15e3;            % subcarrier spacing
N_RBG = 8;        % number of RBGs
TONE_RBG = 1;     % number of tones per RBG modeled in this simulation

dl_txp_re = 10^((TXP_BS-30)/10)/N_TONE; % downlink transmission power per RE
ul_txp_re = 10^((TXP_MS-30)/10)/N_TONE; % uplink transmission power per RE

dl_np_re = 10^((-174+10*log10(SC_F)+NF_MS-30)/10);  % downlink noise power per RE
ul_np_re = 10^((-174+10*log10(SC_F)+NF_BS-30)/10);  % uplink noise power per RE

% input file
load('.\antenna_64x4_70ue.mat');

% output file
for i=1:length(linkInfo)
  linkInfo_cp(i).dh_f = linkInfo(i).dh_f;   % freq domain DL channel matrices per RBG
  linkInfo_cp(i).uh_f = permute(linkInfo(i).dh_f,[2 1 3]);  % freq domain UL channel matrices
  linkInfo(i).ul_txp_re = ul_txp_re;
  linkInfo(i).dl_txp_re = dl_txp_re;
  linkInfo(i).rate = [];
  linkInfo(i).sinr_zf = [];         % for zero forcing
  linkInfo(i).sinr_intf_prob = [];  % for BiT
end

N_RX_U = size(linkInfo_cp(1).uh_f,1);   % uplink number of receive antennas
N_TX_U = size(linkInfo_cp(1).uh_f,2);   % uplink number of transmit antennas

N_RX_D = size(linkInfo_cp(1).dh_f,1);   % downlink number of receive antennas
N_TX_D = size(linkInfo_cp(1).dh_f,2);   % downlink number of transmit antennas

%% per RBG simulation
for i_rbg=1:N_RBG
  disp(['Processing RBG : ' num2str(i_rbg) '/' num2str(N_RBG)]);
  idx_s = (i_rbg-1)*TONE_RBG+1;
  idx_e = i_rbg*TONE_RBG;  

  for i=1:length(linkInfo)
    linkInfo(i).dv_svd = [];
    linkInfo(i).uv_svd = [];
    linkInfo(i).dh_f = linkInfo_cp(i).dh_f(:,:,idx_s:idx_e);
    linkInfo(i).uh_f = linkInfo_cp(i).uh_f(:,:,idx_s:idx_e);
    
    linkInfo(i).on = 1;
    linkInfo(i).ri = MAX_RANK_CPE;
    linkInfo(i).uv_tx = [];
  end

  %calculate SVD precoder for each served CPE
  for i=1:length(linkInfo)
    if linkInfo(i).servcell==1
      cov_h = linkInfo(i).dh_f(:,:,1)'*linkInfo(i).dh_f(:,:,1);
      for i_tone=2:TONE_RBG
        cov_h = cov_h+linkInfo(i).dh_f(:,:,i_tone)'*linkInfo(i).dh_f(:,:,i_tone);
      end
      cov_h = cov_h/TONE_RBG;
      [~,~,v] = svd(cov_h);
      linkInfo(i).dv_svd = v(:,1:MAX_RANK_CPE); % downlink precoder based on SVD

      cov_h = linkInfo(i).uh_f(:,:,1)'*linkInfo(i).uh_f(:,:,1);
      for i_tone=2:TONE_RBG
        cov_h = cov_h+linkInfo(i).uh_f(:,:,i_tone)'*linkInfo(i).uh_f(:,:,i_tone);
      end
      cov_h = cov_h/TONE_RBG;
      [~,~,v] = svd(cov_h);
      linkInfo(i).uv_svd = v(:,1:MAX_RANK_CPE); % uplink precoder based on SVD
    end
  end
  
  %% uplink simulation: generate signal from served cpe
  bs_info = [];
  for i_c = 1:7     % 7 cell sites
    for i_s=1:3     % 3 sectors per cell site
      bs_info(i_c,i_s).h_cov_s = zeros(N_RX_U,N_RX_U);
      idx_serv = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
      if isempty(idx_serv); continue; end
      
      for i_cpe=idx_serv    % for each served cpe
        tx_v = linkInfo(i_cpe).uv_svd(:,1:linkInfo(i_cpe).ri);
        tx_v = tx_v/norm(tx_v,'fro')*sqrt(linkInfo(i_cpe).ul_txp_re);
        linkInfo(i_cpe).uv_tx = tx_v;
        
        eff_uhf = [];
        for j=1:TONE_RBG
          eff_uhf(:,:,j) = linkInfo(i_cpe).uh_f(:,:,j)*linkInfo(i_cpe).uv_tx;   % uplink channel matrix x uplink precoder
          bs_info(i_c,i_s).h_cov_s = bs_info(i_c,i_s).h_cov_s+eff_uhf(:,:,j)*eff_uhf(:,:,j)';
        end
        linkInfo(i_cpe).H_nb = eff_uhf;
      end
    end
  end
  
  %% uplink simulation: generate signal from interfering cpe
  for i_c = 1:7
    for i_s=1:3
      bs_info(i_c,i_s).h_cov_i = zeros(N_RX_U,N_RX_U);
      
      i_intf = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==0));
      if isempty(i_intf); continue; end
      
      for i_cpe=i_intf      % for each interfering cpe
        cpe_id = linkInfo(i_cpe).cpe_ID;
        i_serv = find(([linkInfo(:).cpe_ID]==cpe_id)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
        if isempty(i_serv); continue; end
        
        intf_eff_uhf = [];
        for i=1:TONE_RBG
          intf_eff_uhf(:,:,i) = linkInfo(i_cpe).uh_f(:,:,i)*linkInfo(i_serv).uv_tx; % uplink interfering channel matrix x uplink precoder
          bs_info(i_c,i_s).h_cov_i = bs_info(i_c,i_s).h_cov_i+intf_eff_uhf(:,:,i)*intf_eff_uhf(:,:,i)';
        end
      end
      
      bs_info(i_c,i_s).h_cov_s = bs_info(i_c,i_s).h_cov_s/TONE_RBG;
      bs_info(i_c,i_s).h_cov_i = bs_info(i_c,i_s).h_cov_i/TONE_RBG;
      bs_info(i_c,i_s).h_cov = bs_info(i_c,i_s).h_cov_s+bs_info(i_c,i_s).h_cov_i;
      bs_info(i_c,i_s).h_cov = bs_info(i_c,i_s).h_cov+ul_np_re*eye(size(bs_info(i_c,i_s).h_cov));   % emulate uplink receive noise      
      
      h_cov_intf_prob = bs_info(i_c,i_s).h_cov;
      idx_serv = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
      if isempty(idx_serv); continue; end
      
      %bit precoder
      for i_cpe=idx_serv
        linkInfo(i_cpe).v_bit = [];
        
        for i=1:linkInfo(i_cpe).ri
          linkInfo(i_cpe).v_bit(:,i) = inv(h_cov_intf_prob)*linkInfo(i_cpe).H_nb(:,i);
          linkInfo(i_cpe).v_bit(:,i) = linkInfo(i_cpe).v_bit(:,i)/norm(linkInfo(i_cpe).v_bit(:,i))*sqrt(linkInfo(i_cpe).dl_txp_re/linkInfo(i_cpe).ri)/sqrt(length(idx_serv));   % equal power allocation per layer
          linkInfo(i_cpe).v_bit(:,i) = conj(linkInfo(i_cpe).v_bit(:,i));    % conjugation needed for downlink transmission
        end
      end
    end
  end  
  
  %% downlink simulation: transmit with ZF precoder for ZF vs with BiT precoder for BiT; simulate the data rx at CPEs
  %ZF precoder
  for i_c = 1:7
    for i_s=1:3
      idx_serv = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
      if isempty(idx_serv); continue; end
      
      V = [];
      for i_cpe=idx_serv
        v = linkInfo(i_cpe).dv_svd(:,1:linkInfo(i_cpe).ri);
        V = [V v];
      end
      V = V*inv(V'*V);
      for i=1:size(V,2)
        V(:,i) = V(:,i)/norm(V(:,i));
      end
      
      i_V = 1;
      for i_cpe=idx_serv
        linkInfo(i_cpe).v_zf = V(:,i_V:i_V+linkInfo(i_cpe).ri-1);
        linkInfo(i_cpe).v_zf = linkInfo(i_cpe).v_zf/norm(linkInfo(i_cpe).v_zf,'fro')*sqrt(linkInfo(i_cpe).dl_txp_re)/sqrt(length(idx_serv));
        i_V = i_V+linkInfo(i_cpe).ri;
      end
    end
  end
  
  %calcuate rate
  for i_cpe=1:length(linkInfo)
    if linkInfo(i_cpe).servcell==0 || linkInfo(i_cpe).on==0 || linkInfo(i_cpe).ri<=0; continue; end
    
    cpe_id = linkInfo(i_cpe).cpe_ID;
    
    %zf
    r_cov = zeros(N_RX_D,N_RX_D,TONE_RBG);
    for i_c = 1:7
      for i_s=1:3
        idx_serv = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
        if isempty(idx_serv); continue; end
        
        for i_serv=idx_serv
          tx_v = linkInfo(i_serv).v_zf;
          i_link = find(([linkInfo(:).cpe_ID]==cpe_id)&([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s));
          if isempty(i_link); continue; end
          
          for i=1:TONE_RBG
            for j=1:linkInfo(i_serv).ri
              eff_dhf = linkInfo(i_link).dh_f(:,:,i)*tx_v(:,j);
              r_cov(:,:,i) = r_cov(:,:,i)+eff_dhf*eff_dhf';
            end
          end
        end
      end
    end
    
    linkInfo(i_cpe).rate(end+1).zf = 0;
    tx_v = linkInfo(i_cpe).v_zf;
    for i=1:TONE_RBG
      for j=1:linkInfo(i_cpe).ri
        eff_dhf = linkInfo(i_cpe).dh_f(:,:,i)*tx_v(:,j);
        cov = r_cov(:,:,i)-eff_dhf*eff_dhf';
        cov = cov+dl_np_re*eye(size(cov));
        linkInfo(i_cpe).rate(end).zf = linkInfo(i_cpe).rate(end).zf+real(log2(1+eff_dhf'*inv(cov)*eff_dhf));
        linkInfo(i_cpe).sinr_zf(end+1) = 10*log10(eff_dhf'*inv(cov)*eff_dhf);
      end
    end
    
    %generate precoder based on interference probing, i.e., BiT
    r_cov = zeros(N_RX_D,N_RX_D,TONE_RBG);
    for i_c = 1:7
      for i_s=1:3
        idx_serv = find(([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s)&([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
        if isempty(idx_serv); continue; end
        
        for i_serv=idx_serv
          tx_v = linkInfo(i_serv).v_bit;
          i_link = find(([linkInfo(:).cpe_ID]==cpe_id)&([linkInfo(:).cell_ID]==i_c)&([linkInfo(:).sector_ID]==i_s));
          if isempty(i_link); continue; end
          
          for i=1:TONE_RBG
            for j=1:linkInfo(i_serv).ri
              eff_dhf = linkInfo(i_link).dh_f(:,:,i)*tx_v(:,j);
              r_cov(:,:,i) = r_cov(:,:,i)+eff_dhf*eff_dhf';
            end
          end
        end
      end
    end
    
    linkInfo(i_cpe).rate(end).intf_prob = 0;
    tx_v = linkInfo(i_cpe).v_bit;
    for i=1:TONE_RBG
      for j=1:linkInfo(i_cpe).ri
        eff_dhf = linkInfo(i_cpe).dh_f(:,:,i)*tx_v(:,j);
        cov = r_cov(:,:,i)-eff_dhf*eff_dhf';
        cov = cov+dl_np_re*eye(size(cov));
        linkInfo(i_cpe).rate(end).intf_prob = linkInfo(i_cpe).rate(end).intf_prob+real(log2(1+eff_dhf'*inv(cov)*eff_dhf));
        linkInfo(i_cpe).sinr_intf_prob(end+1) = 10*log10(eff_dhf'*inv(cov)*eff_dhf);
      end
    end
  end
end

%% Simulation results: ZF vs BiT
idx = find(([linkInfo(:).servcell]==1)&([linkInfo(:).on]==1)&([linkInfo(:).ri]>0));
link_result = linkInfo(idx);

%rate
cpe_rate = [link_result(:).rate];
se_zf = [cpe_rate(:).zf]/TONE_RBG;
se_intf_prob = [cpe_rate(:).intf_prob]/TONE_RBG;

[x_se_zf,cdf_se_zf] = myCDF(se_zf);
[x_se_intf_prob,cdf_se_intf_prob] = myCDF(se_intf_prob);

figure
set(gcf,'Position',[100 100 100+700 100+500]);
plot(x_se_zf,cdf_se_zf, x_se_intf_prob,cdf_se_intf_prob);
grid;
legend('zf','intf-probe','location','best');
xlabel('(b/s/Hz)');
ylabel('CDF')
title(['SE (rank=' num2str(MAX_RANK_CPE) ')'])

SE_zf = sum(double(se_zf))/N_RBG/21
SE_intf_prob = sum(double(se_intf_prob))/N_RBG/21

%sinr
sinr_zf = real([link_result(:).sinr_zf]);
sinr_intf_prob = real([link_result(:).sinr_intf_prob]);

[x_sinr_zf,cdf_sinr_zf] = myCDF(sinr_zf);
[x_sinr_intf_prob,cdf_sinr_intf_prob] = myCDF(sinr_intf_prob);

figure
set(gcf,'Position',[100 100 100+700 100+500]);
plot(x_sinr_zf,cdf_sinr_zf, x_sinr_intf_prob,cdf_sinr_intf_prob);
grid;
legend('zf','intf-probe','location','best');
xlabel('(dB)');
ylabel('CDF')
title(['SINR (rank=' num2str(MAX_RANK_CPE) ')'])