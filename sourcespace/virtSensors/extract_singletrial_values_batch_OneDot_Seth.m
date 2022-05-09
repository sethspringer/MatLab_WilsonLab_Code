function extract_singletrial_values_OneDot_Seth

[Files1,PathName1,~] = uigetfile('*.tfcs','Please Select OneDot TFCSs','MultiSelect','on');

%%%%%PARADIGM%%%%%
%This paradigm is OneDot with an epoch of -600 to 3200, baseline of -600 to
%0, and active/entrainment period of 200 to 2000%


baseline_times = [-600 0];
baseline_start_time = -600;
baseline_time_res = 100;
%Index the locations of the where the baseline starts (1,1) and ends (1,2)%
baseline_time_index = ((baseline_times-baseline_start_time)./baseline_time_res);
baseline_time_index(1,1) = baseline_time_index(1,1)+1;



entrain_freqs = [14.5 15.5];
entrain_freq_res = .5;
entrain_start_freq = 4;
entrain_times = [200 2000];  
entrain_time_res = 100;
entrain_start_time = -600;
%Index the locations of the where the active freq starts (1,1) and ends (1,2)%
entrain_freq_index = ((entrain_freqs-entrain_start_freq)./entrain_freq_res);
entrain_freq_index(1,1) = entrain_freq_index(1,1)+1;
%Index the locations of the where the active time starts (1,1) and ends (1,2)%
entrain_time_index = ((entrain_times-entrain_start_time)./entrain_time_res);
entrain_time_index(1,1) = entrain_time_index(1,1)+1;







% alpha_freqs = [10 16];
% alpha_freq_res = 2;
% alpha_start_freq = 4;
% alpha_times = [350 1950];
% alpha_time_res = 25;
% alpha_start_time = -1800;
% alpha_freq_index = ((alpha_freqs-alpha_start_freq)./alpha_freq_res);
% alpha_freq_index(1,1) = alpha_freq_index(1,1)+1;
% alpha_time_index = ((alpha_times-alpha_start_time)./alpha_time_res);
% alpha_time_index(1,1) = alpha_time_index(1,1)+1;
% 
% beta_freqs = [18 24];
% beta_freq_res = 2;
% beta_start_freq = 4;
% beta_times = [350 1950];
% beta_time_res = 25;
% beta_start_time = -1800;
% beta_freq_index = ((beta_freqs-beta_start_freq)./beta_freq_res);
% beta_freq_index(1,1) = beta_freq_index(1,1)+1;
% beta_time_index = ((beta_times-beta_start_time)./beta_time_res);
% beta_time_index(1,1) = beta_time_index(1,1)+1;
% 
% gamma_freqs = [62 92];
% gamma_freq_res = 2;
% gamma_start_freq = 4;
% gamma_times = [200 1800];
% gamma_time_res = 25;
% gamma_start_time = -1800;
% gamma_freq_index = ((gamma_freqs-gamma_start_freq)./gamma_freq_res);
% gamma_freq_index(1,1) = gamma_freq_index(1,1)+1;
% gamma_time_index = ((gamma_times-gamma_start_time)./gamma_time_res);
% gamma_time_index(1,1) = gamma_time_index(1,1)+1;






%%%%ENTRAIN%%%%


cd(PathName1)
for i = 1:length(Files1)
    tfc_data = readBESAtfcs(Files1{i});
    for ii = 1:size(tfc_data.trials,2)
        tfc_data.trials{ii} = abs(tfc_data.trials{ii});    %Converts the complex numbers into single values per time/freq%
        %VectorSumming the two orientations and then getting rid of the extra dimension%
        tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2))); 
        %Grab the time/frequency bins of interest for the ACTIVE window%
        target_data = tfc_data.trials{ii}(entrain_freq_index(1,1):entrain_freq_index(1,2),entrain_time_index(1,1):entrain_time_index(1,2));
        %Grab the time/frequency bins of interest for the BASELINE window%
        baseline_data = tfc_data.trials{ii}(entrain_freq_index(1,1):entrain_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
        %Average over the time/frequency windows of interest for the active and baseline windows%
        %old way
        target_data_mean = mean(target_data(:));
        baseline_data_mean = mean(baseline_data(:));
        
        trialvalues(i,ii) = (target_data_mean-baseline_data_mean)/baseline_data_mean;

               
        %new appropriate way
        trialvalues_freqWiseCorrection(i,ii) = mean((target_data-mean(baseline_data,2))./(mean(baseline_data,2)),'all');
        trialvalues_raw(i,ii) = mean(target_data,'all');
    end
end


fid = fopen('Entrain_OneDot_singletrial_baselinecorrected.txt','wt');
fprintf(fid,'File');
for i = 1:size(trialvalues,1)
    fprintf(fid,'\n%s\t',Files1{i});
    for ii = 1:size(trialvalues,2)
        fprintf(fid,'%.4f\t',trialvalues(i,ii));
    end
end
fclose(fid);

fid = fopen('Entrain_OneDot_singletrial_baselinecorrected_FreqWise.txt','wt');
fprintf(fid,'File');
for i = 1:size(trialvalues_freqWiseCorrection,1)
    fprintf(fid,'\n%s\t',Files1{i});
    for ii = 1:size(trialvalues_freqWiseCorrection,2)
        fprintf(fid,'%.4f\t',trialvalues_freqWiseCorrection(i,ii));
    end
end
fclose(fid);


fid = fopen('Entrain_OneDot_singletrial_raw.txt','wt');
fprintf(fid,'File');
for i = 1:size(trialvalues_raw,1)
    fprintf(fid,'\n%s\t',Files1{i});
    for ii = 1:size(trialvalues_raw,2)
        fprintf(fid,'%.4f\t',trialvalues_raw(i,ii));
    end
end
fclose(fid);


clear trialvalues trialvalues_raw trialvalues_freqWiseCorrection tfc_data i ii target_data






















% %%%%ALPHA%%%%
% 
% cd(PathName1)
% for i = 1:length(Files1)
%     tfc_data = readBESAtfcs(Files1{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),alpha_time_index(1,1):alpha_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('AlphaAVE-Lmotor_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files1{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName2)
% for i = 1:length(Files2)
%     tfc_data = readBESAtfcs(Files2{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),alpha_time_index(1,1):alpha_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% fid = fopen('AlphaAVE-Locc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files2{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName3)
% for i = 1:length(Files3)
%     tfc_data = readBESAtfcs(Files3{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),alpha_time_index(1,1):alpha_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('AlphaAVE-Lparietal_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files3{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName4)
% for i = 1:length(Files4)
%     tfc_data = readBESAtfcs(Files4{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),alpha_time_index(1,1):alpha_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(alpha_freq_index(1,1):alpha_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('AlphaAVE-Rocc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files4{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% 
% 
% 
% %%%%BETA%%%%
% 
% 
% cd(PathName5)
% for i = 1:length(Files5)
%     tfc_data = readBESAtfcs(Files5{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),beta_time_index(1,1):beta_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('BetaAVE-Locc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files5{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName6)
% for i = 1:length(Files6)
%     tfc_data = readBESAtfcs(Files6{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),beta_time_index(1,1):beta_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('BetaAVE-Rocc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files6{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName10)
% for i = 1:length(Files10)
%     tfc_data = readBESAtfcs(Files10{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),beta_time_index(1,1):beta_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(beta_freq_index(1,1):beta_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% fid = fopen('BetaPairedT-Lparietocc-main_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files10{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%GAMMA%%%%
% 
% 
% cd(PathName7)
% for i = 1:length(Files7)
%     tfc_data = readBESAtfcs(Files7{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('GammaAVE-Locc-lateral_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files7{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName8)
% for i = 1:length(Files8)
%     tfc_data = readBESAtfcs(Files8{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('GammaAVE-Locc-medial_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files8{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName9)
% for i = 1:length(Files9)
%     tfc_data = readBESAtfcs(Files9{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% fid = fopen('GammaAVE-Rocc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files9{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% 
% cd(PathName11)
% for i = 1:length(Files11)
%     tfc_data = readBESAtfcs(Files11{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('GammaPairedT-Locc_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files11{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% cd(PathName12)
% for i = 1:length(Files12)
%     tfc_data = readBESAtfcs(Files12{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% 
% fid = fopen('GammaPairedT-Rocc-inferior_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files12{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% 
% cd(PathName13)
% for i = 1:length(Files13)
%     tfc_data = readBESAtfcs(Files13{i});
%     for ii = 1:size(tfc_data.trials,2)
%         tfc_data.trials{ii} = abs(tfc_data.trials{ii});
%         tfc_data.trials{ii} = squeeze(sqrt((tfc_data.trials{ii}(1,:,:).^2)+(tfc_data.trials{ii}(2,:,:).^2)));
%         target_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),gamma_time_index(1,1):gamma_time_index(1,2));
%         baseline_data = tfc_data.trials{ii}(gamma_freq_index(1,1):gamma_freq_index(1,2),baseline_time_index(1,1):baseline_time_index(1,2));
%         target_data = mean(target_data(:));
%         baseline_data = mean(baseline_data(:));
%         trialvalues(i,ii) = (target_data-baseline_data)/baseline_data;
%     end
% end
% fid = fopen('GammaPairedT-Rocc-superior_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(trialvalues,1)
%     fprintf(fid,'\n%s\t',Files13{i});
%     for ii = 1:size(trialvalues,2)
%         fprintf(fid,'%.4f\t',trialvalues(i,ii));
%     end
% end
% fclose(fid);
% clear trialvalues tfc_data i ii target_data
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%EXTRACT BEHAVIOR%%%%
% function extract_accepted_RTs_MRT
% 
% cd('E:\Mental Rotation Task\fifs_check\Stats\evts\evt_postrecode_reextracted')
% files = dir('*.evt');
% for i = 1:length(files)
% full_evts{i} = readBESAevt(files(i).name);
% end
% 
% cd('E:\Mental Rotation Task\fifs_check\Stats\evts\evt_accepted')
% files = dir('*.evt');
% for i = 1:length(files)
% accepted_evts{i} = readBESAevt(files(i).name);
% end
% 
% for i = 1:length(files)
% for ii = 1:length(accepted_evts{i}(:,1))
% accepted_evt_indices{i}(ii,1) = find(full_evts{i}(:,1) == accepted_evts{i}(ii,1));
% end
% end
% 
% for i = 1:length(files)
% for ii = 1:length(accepted_evt_indices{i})
% RTs{i}(ii,1) = (full_evts{i}(accepted_evt_indices{i}(ii)+1,1)-full_evts{i}(accepted_evt_indices{i}(ii),1))/1000;
% if full_evts{i}(accepted_evt_indices{i}(ii),3) == 12290 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12291 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12292 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12293
% Condition{i,ii} = 'LowRotation'; 
% elseif full_evts{i}(accepted_evt_indices{i}(ii),3) == 12294 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12295 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12296 | full_evts{i}(accepted_evt_indices{i}(ii),3) == 12297
% Condition{i,ii} = 'HighRotation'; 
% end
% end
% end
% 
% 
% %%%% pull timestamps for single trial triggers:
% for i = 1:length(files)
% LowRotation_counter = 1;
% HighRotation_counter = 1;
% for ii = 1:length(RTs{i})
% if strcmp(Condition{i,ii},'LowRotation')
% LowRotation_RT{i}(LowRotation_counter) = RTs{i}(ii);
% LowRotation_counter = LowRotation_counter+1;
% elseif strcmp(Condition{i,ii},'HighRotation')
% HighRotation_RT{i}(HighRotation_counter) = RTs{i}(ii);
% HighRotation_counter = HighRotation_counter+1;
% end
% end
% end
% 
% 
% %% save condition-wise accepted RTs
% 
% fid = fopen('RTs_LowRotation_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(LowRotation_RT,2)
%     fprintf(fid,'\n%s\t',files(i).name(1:end-4));
%     for ii = 1:size(LowRotation_RT{i},2)
%         fprintf(fid,'%.4f\t',LowRotation_RT{i}(1,ii));
%     end
% end
% fclose(fid);
% 
% fid = fopen('RTs_HighRotation_singletrial_baselinecorrected.txt','wt');
% fprintf(fid,'File');
% for i = 1:size(HighRotation_RT,2)
%     fprintf(fid,'\n%s\t',files(i).name(1:end-4));
%     for ii = 1:size(HighRotation_RT{i},2)
%         fprintf(fid,'%.4f\t',HighRotation_RT{i}(1,ii));
%     end
% end
% fclose(fid);
