clear all
close all

PLX_Path='D:\0922';
PLX_FileName=[{'WCJ092102'}];
Select_CH=[5];
fmin=0;fmax=100;  
GuangXian_Path='D:\0922';
GuangXian_Name=[{'2021_09_22-09_22_49'}];
Select_GuangXian_Channel=2;
GuangXian_start_time=0.9156;  
GuangXian_time=10.76; 
Delete_GuangXian_record_number=[];
Cell_filename=[];
PLX_rate=1000;
GuangXian_rate=15; 
No_Channel=8; 
Output_Path='D:\Original_Data';
addpath('D:\Work\RWDfiber_elecrode')
addpath(genpath('D:\My Documents\MATLAB\Work\chronux_2_12'))

Read_GuangXian_Channel=Select_GuangXian_Channel+1;
Stim_Count=0;
Stimuli_Time2=[];
Stimuli_Time1=[];
cw=[]; 
Cell=[];
Wave_Data=[];

for i=1:length(PLX_FileName)
    
    FileName=strcat(PLX_Path,'\',PLX_FileName{i},'.PLX');
    new_cw=[]; 
    for CH_ID=1:No_Channel
       qqw=CH_ID;
     
       if CH_ID<10
        CH_ID=strcat('FP0',num2str(CH_ID));
       else
         CH_ID=strcat('FP',num2str(CH_ID));   
       end
       [~, ~, ~, ~, ad] = plx_ad_v(FileName, CH_ID);
       new_cw(:,qqw)=ad;
    end
    
    Ev_CH=2;
    [b_Stim_Count, b_Stimuli_Time, ~] = plx_event_ts(FileName, Ev_CH);
    Stim_Count=Stim_Count+b_Stim_Count;
    Stimuli_Time2=[Stimuli_Time2;length(cw)/PLX_rate+b_Stimuli_Time];
    
    Ev_CH=1;
    [a_Stim_Count, a_Stimuli_Time, sv] = plx_event_ts(FileName, Ev_CH);
    Stimuli_Time1=[Stimuli_Time1;length(cw)/PLX_rate+a_Stimuli_Time];
    
    cw=[cw;new_cw];
    
  
    if length(Cell_filename)>=1
        CellName=strcat(GuangXian_Path,'\',Cell_filename{i},'.mat');
        load(CellName,'Spike_Continous','Wave_Epoch_continue','LFP_array');
        Cells_name=fieldnames(Spike_Continous);
        Cells= struct2cell(Spike_Continous);
        for ii=1:length(Cells_name)
           Cell_name=Cells_name{ii};
           Cell_number=str2num(Cell_name(5));
           if Cell_number==Select_Cell_CH
               new_Cell=Cells{ii};
               if isempty(Cell)==0
                   new_Cell=Cell(length(Cell))+new_Cell;
               end
           end
        end
        new_Wave_Data=[];
        Cell=[Cell;new_Cell];
        EEG_Ch=find(LFP_array==Select_Cell_CH);
        new_Wave_Data(1,:)=Wave_Epoch_continue(EEG_Ch,:); % as EEG
        new_Wave_Data(2,:)=mean([Wave_Epoch_continue(1,:);Wave_Epoch_continue(end,:)]);% as EMG
        Wave_Data=[Wave_Data new_Wave_Data];
    end
end
Stim_Count=Stim_Count-length(Delete_GuangXian_record_number);
Stimuli_Time2(Delete_GuangXian_record_number)=[];

GuangXian_Name=GuangXian_Name(setdiff(1:length(GuangXian_Name),Delete_GuangXian_record_number));
if Stim_Count==length(GuangXian_Name)
else
error('ERROR:The number of files does not match the number of EEG records.')
end




window=10;
movingwin= [window,1];   
params.Fs=1000; 
params.fpass=[fmin fmax];
params.tapers=[7,9];
params.err=[1,0];
Min_X=0; Max_X=1200;Min_Color=-0.5;Max_Color=0.5;


total_Guangxian_Data=[];
total_CH_data=[];
figure(1)
sub_line=ceil(length(GuangXian_Name)/2);
for file_number=1:length(GuangXian_Name)
    File_Name=[];
    File_Name=strcat(GuangXian_Path,'\',GuangXian_Name{file_number},'\Fluorescence_DFF.csv');
    GuangXian_Data=[];
    GuangXian_Data=xlsread(File_Name);
    
    GuangXian_Data=GuangXian_Data(GuangXian_start_time*60*GuangXian_rate:(GuangXian_start_time*60+GuangXian_time*60)*GuangXian_rate-1,Read_GuangXian_Channel);
   
    [min_GuangXian_Data, I] = sort(GuangXian_Data,1);
    min_GuangXian_Data=min_GuangXian_Data(1:GuangXian_time*2);
    GuangXian_Data=GuangXian_Data-mean(min_GuangXian_Data);
    
    total_Guangxian_Data=[total_Guangxian_Data;GuangXian_Data];

    CH_start_time=(Stimuli_Time2(file_number)+GuangXian_start_time*60-window/2)*PLX_rate;
    CH_data=cw(CH_start_time:CH_start_time+GuangXian_time*60*PLX_rate+window*PLX_rate-1,Select_CH);
    if exist('CH_threshold','var')==1
        for ii=1:length(CH_data)
            if CH_data(ii)>CH_threshold
                CH_data(ii)=CH_threshold;
            elseif CH_data(ii)<-1*CH_threshold
                CH_data(ii)=-1*CH_threshold;
            end
        end
    end
    [S,t,f,Serr]=mtspecgramc(CH_data,movingwin,params);
    if length(GuangXian_Name)>1
        if file_number==1
            CH_start_time_a=(Stimuli_Time2(file_number)+GuangXian_start_time*60-window/2)*PLX_rate;
            CH_data_a=cw(CH_start_time:CH_start_time+GuangXian_time*60*PLX_rate-1,Select_CH);
        elseif file_number==length(GuangXian_Name)
            CH_start_time_a=(Stimuli_Time2(file_number)+GuangXian_start_time*60)*PLX_rate;
            CH_data_a=cw(CH_start_time:CH_start_time+GuangXian_time*60*PLX_rate+window*PLX_rate-1,Select_CH);
        else
            CH_start_time_a=(Stimuli_Time2(file_number)+GuangXian_start_time*60)*PLX_rate;
            CH_data_a=cw(CH_start_time:CH_start_time+GuangXian_time*60*PLX_rate-1,Select_CH);
        end
        total_CH_data=[total_CH_data;CH_data_a];
    end



    subplot(sub_line,2,file_number)
    Mean=mean(S);  
    for a=1:length(Mean)
        S(:,a)=log10(S(:,a)/Mean(a));
    end;
    image(t,f,S','CDataMapping','scaled');
    set(gca,'CLim',[Min_Color Max_Color],'TickDir','out','YDir','normal');
    hold on;
    GuangXian_Data=resample(GuangXian_Data,length(t),length(GuangXian_Data));
    GuangXian_Data=(GuangXian_Data-min(GuangXian_Data))/(max(GuangXian_Data)-min(GuangXian_Data))*fmax;
    plot(t,GuangXian_Data,'w');
    
end


if length(GuangXian_Name)==1
    figure(2)
    
    k=total_Guangxian_Data;
    m=mean(total_Guangxian_Data);
    s=std(total_Guangxian_Data,0);
    plot(k);
    ylabel('deltaF/F');
    xlabel('Time(s)');
    axis([0,length(k),min(k)-1,max(k)+2]);
    hold on;
    plot([0,length(k)],[m,m]);
    hold on;
    t=m+s;
    plot([0,length(k)],[t,t],'--');
    hold on;
    [amp,dec,firerate]=func_calculate_spike(k,m,s,GuangXian_time);
    save(strcat(GuangXian_Path,'\spikedata.mat'),'amp','dec','firerate');
else
    
    figure(2)
    subplot(2,1,2)
    if exist('CH_threshold','var')==1
        for ii=1:length(total_CH_data)
            if total_CH_data(ii)>CH_threshold
                total_CH_data(ii)=CH_threshold;
            elseif total_CH_data(ii)<-1*CH_threshold
                total_CH_data(ii)=-1*CH_threshold;
            end
        end
    end
    [S,t,f,Serr]=mtspecgramc(total_CH_data,movingwin,params);
    Mean=mean(S);  
    for a=1:length(Mean)
        S(:,a)=log10(S(:,a)/Mean(a));
    end;
    image(t,f,S','CDataMapping','scaled');
    set(gca,'CLim',[Min_Color Max_Color],'TickDir','out','YDir','normal');
    hold on;
    k=resample(total_Guangxian_Data,length(t),length(total_Guangxian_Data));
    k=(k-min(k))/(max(k)-min(k))*fmax;
    m=mean(k);
    s=std(k,0);
    plot(k,'w','LineWidth',2);
    ylabel('deltaF/F');
    xlabel('Time(s)');
    axis([0,length(k),min(k)-1,max(k)+2]);
    hold on;
    plot([0,length(k)],[m,m],'w');
    hold on;
    t=m+s;
    plot([0,length(k)],[t,t],'w--');
    hold on;
    [amp,dec,firerate]=func_calculate_spike(k,m,s,GuangXian_time);
    hold on;
    
    subplot(2,1,1)
    fullrecord_CH_data=cw(:,Select_CH);
    if exist('CH_threshold','var')==1
        for ii=1:length(fullrecord_CH_data)
            if fullrecord_CH_data(ii)>CH_threshold
                fullrecord_CH_data(ii)=CH_threshold;
            elseif fullrecord_CH_data(ii)<-1*CH_threshold
                fullrecord_CH_data(ii)=-1*CH_threshold-1*CH_threshold;
            end
        end
    end
    [S,t,f,Serr]=mtspecgramc(fullrecord_CH_data,movingwin,params);
    Mean=mean(S);  
    for a=1:length(Mean)
        S(:,a)=log10(S(:,a)/Mean(a));
    end;
    image(t,f,S','CDataMapping','scaled');
    set(gca,'CLim',[Min_Color Max_Color],'TickDir','out','YDir','normal');
    hold on;
    if length(Cell_filename)>=1
        fs=1000;
        DurEpo=1;
        [stgS,time,PSDx,fx,timex] = SleepAnalysisRV1(Wave_Data,DurEpo,FileName,fs);
        hh=hist(Cell,timex*60);
        standard_hh=(hh-min(hh))/(max(hh)-min(hh));
        hh=100*standard_hh;
        hh=smooth(hh,0.01,'loess');
        plot(hh,'w','LineWidth',2);
    end

    for file_number=1:length(GuangXian_Name)
        if file_number < length(GuangXian_Name)
           figure(2)
           subplot(2,1,2)
           line_GuangXian=file_number*GuangXian_time*60;
           plot([line_GuangXian,line_GuangXian],[min(k)-1,max(k)+1],'r','LineWidth',2);
           hold on;
        end
        figure(2)
        subplot(2,1,1)
        line_PLX_1=Stimuli_Time2(file_number)+GuangXian_start_time*60-window/2;
        line_PLX_2=line_PLX_1+GuangXian_time*60+window;
        plot([line_PLX_1,line_PLX_2],[-0.5,-0.5],'b','LineWidth',5);
        axis([0,length(t),-1,100]);
    end

if exist('GuangXian_analyze_mode','var')==1
    switch(GuangXian_analyze_mode)
        case 'group'
           total_control_amp=[];
           total_control_dec=[];
           total_control_firerate=[];
           total_exp_amp=[];
           total_exp_dec=[];
           total_exp_firerate=[];
           m=mean(total_Guangxian_Data);
           s=std(total_Guangxian_Data,0);
           for file_number=1:length(GuangXian_Name)
               kk=total_Guangxian_Data((file_number-1)*GuangXian_time*60*GuangXian_rate+1:file_number*GuangXian_time*60*GuangXian_rate);
               figure(length(GuangXian_Name)+2)
               [amp,dec,firerate]=func_calculate_spike(kk,m,s,GuangXian_time);
               close(figure(length(GuangXian_Name)+2));
               if file_number<=2
                   total_control_amp=[total_control_amp amp];
                   total_control_dec=[total_control_dec dec];
                   total_control_firerate=[total_control_firerate firerate];
               elseif file_number>=5
                   total_exp_amp=[total_exp_amp amp];
                   total_exp_dec=[total_exp_dec dec];
                   total_exp_firerate=[total_exp_firerate firerate];
               end
           end
           total_control_amp=total_control_amp';
           total_control_dec=total_control_dec';
           total_control_firerate=total_control_firerate';
           total_exp_amp=total_exp_amp';
           total_exp_dec=total_exp_dec';
           total_exp_firerate=total_exp_firerate';
           save(strcat(GuangXian_Path,'\spikedata',num2str(Select_GuangXian_Channel),'.mat'),'total_control_amp','total_control_dec','total_control_firerate','total_exp_amp','total_exp_dec','total_exp_firerate');

        case 'combine'
           total_amp=[];
           total_dec=[];
           total_firerate=[];
           m=mean(total_Guangxian_Data);
           s=std(total_Guangxian_Data,0);
           for file_number=1:length(GuangXian_Name)
               kk=total_Guangxian_Data((file_number-1)*GuangXian_time*60*GuangXian_rate+1:file_number*GuangXian_time*60*GuangXian_rate);
               figure(length(GuangXian_Name)+2)
               [amp,dec,firerate]=func_calculate_spike(kk,m,s,GuangXian_time);
               close(figure(length(GuangXian_Name)+2));
               total_amp=[total_amp amp];
               total_dec=[total_dec dec];
               total_firerate=[total_firerate firerate];
           end
           save(strcat(GuangXian_Path,'\spikedata',num2str(Select_GuangXian_Channel),'.mat'),'total_amp','total_dec','total_firerate');

        case 'single'
            m=mean(total_Guangxian_Data);
            s=std(total_Guangxian_Data,0);
            aaa=0;
            save(strcat(GuangXian_Path,'\spikedata',num2str(Select_GuangXian_Channel),'.mat'),'aaa');
            for file_number=1:length(GuangXian_Name)
               kk=total_Guangxian_Data((file_number-1)*GuangXian_time*60*GuangXian_rate+1:file_number*GuangXian_time*60*GuangXian_rate);
               figure(length(GuangXian_Name)+2)
               [amp,dec,firerate]=func_calculate_spike(kk,m,s,GuangXian_time);
               close(figure(length(GuangXian_Name)+2));
               str_amp=strcat('single_',num2str(file_number),'_amp');
               eval([str_amp,'=amp;']);
               str_dec=strcat('single_',num2str(file_number),'_dec');
               eval([str_dec,'=dec;']);
               str_firerate=strcat('single_',num2str(file_number),'_firerate');
               eval([str_firerate,'=firerate;']);
               save(strcat(GuangXian_Path,'\spikedata',num2str(Select_GuangXian_Channel),'.mat'),...
                   strcat('single_',num2str(file_number),'_amp'),strcat('single_',num2str(file_number),'_dec'),...
                   strcat('single_',num2str(file_number),'_firerate'),'-append');
            end
    end
end

end