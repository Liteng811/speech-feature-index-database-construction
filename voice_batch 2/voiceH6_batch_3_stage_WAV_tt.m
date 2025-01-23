

%%%%%%%%%%%%%%第一部分.循环样本样本分析开始%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%一.0.0前置测试%%%%%%%%
thisfolder=pwd
%{
audiopath='H:\音频开发\音频继续研发20220606集成\V_diagnosis\2023.11.01\杜首燕\ZOOM0013'
cd(audiopath)
[speech,fs]=audioread('ZOOM0013_LR.WAV');
cd(thisfolder)
size(speech)
speech_1=speech(:,1);
[idx,C]=kmeans(abs(speech_1),3);
%}
%{
basic_sec_num=3  %H6实验过程中要调整到13
for i=2:7
    idx=[1 1 2 2 1 1 1 2 2 1 2 1 1 1 1 2 2 3 3 3 3 2 2 1 1 1 2 3 2 1 1 2 2 1 1 1 1 1 1 2 2 3];
    size(idx)
    loc_1=(find(idx==1))%每一段等于特定数字的连续数的位置
    loc_1_dis=diff(loc_1)%每一段等于特定数字的连续数的位置后面减去前面，得到差值，大于一表明其中连续；
    loc_1_dis_all=cat(2,loc_1(1),loc_1_dis) %与位置号码loc_1一一对应的间隔数字
    %find(loc_1_dis>1,1)%第一个
    breakplace_1=[find(loc_1_dis>1,1) diff(find(loc_1_dis>1))] %每一段等于特定数字的连续数的独立长度；
    loc_2=(find(idx==2));%每一段等于特定数字的连续数的位置
    loc_2_dis=diff(loc_2);%每一段等于特定数字的连续数的位置后面减去前面，得到差值，大于一表明其中连续；
    breakplace_2=[find(loc_2_dis>1,1) diff(find(loc_2_dis>1))] %每一段等于特定数字的连续数的独立长度；
    
    if numel(breakplace)<basic_sec_num
        continue
    else
        %if loc_1(1)<loc_2(1)  %如果噪音在前方
            breakplace_1_sorted=sort(breakplace_1,'descend')
            break_length_threathhold=breakplace_1_sorted(basic_sec_num);
            loc=zeros(length(idx),1);
            line=ones(1,break_length_threathhold+1)
            for i=1:length(idx)-break_length_threathhold-1
                if (idx(i:i+break_length_threathhold)==line)
                    loc(i:i+break_length_threathhold+1)=1;
                end
            end
            loc'        
    end
    ssssss
end
                       
     sssssss
%}
%%%如果小于某一个长度限度，就连起来
%idx(1:100:end)%1,2 两类，
%C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%一.1.0参数设置%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
uichoose=1 %如果选择ui设定参数的方式，就选1，全自动就选2
I_want_choose_sample=1;
inGUI=0;%放到GUI中，就变成1，独立运行就是0
repeattime=1;%循环次数，用于写标题,需要标题就写1，不需要就写2以上
    overlapmethod=2;%3级以上文件夹方法，1或2
    txtfileclean=0; %分析文件中，如果没有同名的txt文件选0，有选1
    iwanthearit=0; %不想听选0，想听有选1
    iwantcontent=0;%不想识别语义选0，想识别语义就选1
        %%%%循环起点设置%%%%
    wavstart_i_1=5; %第一层文件夹起点
    wavstart_i_2=1;
    wavstart_i_3=1;
    basic_sec_num=13;  %样本中最小的识别字符数，基于样本设定设定
    basic_sec_num_method=2; %1为kmeans,2为findpeaks
    filtermethod=4; % 1.小波滤波；%2.带通滤波%3.平均滤波；4.平滑
    f_band_low = 0.1; % 低通频率  
    f_band_high = 5; % 高通频率
    lvbochidu=0.01;  %滤波尺度
    figureshow=1;  %显示图形
    unfinished=0  %1等效于观察未完成部分的效果
    namethedatabase='datafileport-2024_1.txt';
    houzuiming_audio_WAV='.WAV';
    houzuiming_audio_m4a='.m4a';
    houzuiming_audio_mp3='.mp3';
    houzuiming_audio_pcm='.pcm';
    houzuiming_audio_txt='.txt';
    uichoose %如果选择ui设定参数的方式，就选1
    if uichoose==1    
        prompt = {'\fontsize{11} 1.循环次数，用于写标题填写1:','\fontsize{11}2.三级以上文件夹方法，1或2:',...
            '\fontsize{11} 3.如果没有同名的txt文件选0，有选1:','\fontsize{11}4.不想听选0，想听有选1:',...
            '\fontsize{11}5.不想识别语义选0，想识别语义就选1:',...
            '\fontsize{11}6.一层文件夹中断继续，参考text显示，默认选1:',...
            '\fontsize{11}7.二层文件夹中断继续，参考text显示，默认选1:',...
            '\fontsize{11}8.三层文件夹中断继续，参考text显示，默认选1:',...
            '\fontsize{11}9.写入数据的文件名称（后缀请用txt格式:'};
        dlg_title = 'Input your information';
        num_lines = 1;
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        %def = {'1','2','0','0','1'};
        def = {num2str(repeattime),num2str(overlapmethod),num2str(txtfileclean),...
            num2str(iwanthearit),num2str(iwantcontent),...
            num2str(wavstart_i_1),num2str(wavstart_i_2),num2str(wavstart_i_3),namethedatabase};
        %def = {P{1}{line},P{2}{line}, P{3}{line}};
        %def = {'阿宝','19860722', '男'};
        mainparameterscell = inputdlg(prompt,dlg_title,num_lines,def,options);
        repeattime=str2double(mainparameterscell{1});
        overlapmethod=str2double(mainparameterscell{2});
        txtfileclean=str2double(mainparameterscell{3});
        iwanthearit=str2double(mainparameterscell{4});
        iwantcontent=str2double(mainparameterscell{5});
        wavstart_i_1=str2double(mainparameterscell{6});
        wavstart_i_2=str2double(mainparameterscell{7});
        wavstart_i_3=str2double(mainparameterscell{8});
        namethedatabase=(mainparameterscell{9});
    end    
    mainparameters=getappdata(0,'mainparameters');
    mainparameters.repeattime=repeattime;
    mainparameters.overlapmethod=overlapmethod;
    mainparameters.txtfileclean=txtfileclean;
    mainparameters.iwanthearit=iwanthearit;
    mainparameters.iwantcontent=iwantcontent;
    setappdata(0,'mainparameters',mainparameters);
    %
    %%%%参数设置完成%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%一.1.1读取五音文件%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fivetonedata='five_tone_data.xlsx';
    [data,~,~]=xlsread(fivetonedata);
    fivetonemel20=data(:, 10:29);   %五音原始mel频谱20分量
    fivetonerank5=data(:, 1:3);
    max_five20=max(fivetonemel20,[],2);  %行向量最大值
    min_five20=min(fivetonemel20,[],2);  %行向量最小值
    M=size(fivetonemel20,2);  %m=20
    maxrank=repmat(max_five20,1,M);  
    minrank=repmat(min_five20,1,M);
    max_five20=(fivetonemel20-minrank)./(maxrank-minrank);    
    %%%%读取五音文件结束%%%%
    clear speech time Yabsvoicefft Yfft Ynoisefft Yvoicefft markspeechtimeline speechenergyrank
    clear MFCCs logFBEs    
    %%%%一.1.2清理之前留下的内存%%%%    
    %%%%一.1.3加载必要的文件%%%%
    this_path = mfilename('fullpath'); %本m文件所在位置
    [thisfolder,~,~]=fileparts(this_path);
    cd(thisfolder);
    cd('..')    
    originalfolder=pwd;
    originalfolder=cat(2,originalfolder,'\');
    VADpath=strcat(originalfolder,'rVAD2.0now\');
    mfccpath=strcat(originalfolder,'mfcc\');
    pythpath=strcat(originalfolder,'pyth\');
    wenetpath=strcat(originalfolder,'wenet\');
    addpath(VADpath); 
    addpath(mfccpath);
    addpath(pythpath); 
    addpath(wenetpath);
    cd(thisfolder);
    [TxtId, ~] = fopen(namethedatabase, 'a+');%在for外面打开文档，避免重复开关%cd('..')
    %%%%加载必要的文件完成%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%一.1.3读取音频文件%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%选择待计算的音频文档%%%%
    if I_want_choose_sample==1
        path_basic_0_stage=uigetdir();%'G:\儿研所论文\听诊样本集合\';
    else
        path_basic_0_stage='H:\音频开发\音频继续研发20220606集成\V_diagnosis\三附院呼吸科数据 未识别\2023.11.20';
    end
    path_1_stage  = dir( path_basic_0_stage );    
    %%%%选择待计算的音频文档完成%%%%
%%%%加载必要的文件完成%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%一.2.0解析音频文件%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%选择待计算的音频文档%%%%
    %%%%开始解析音频文件%%%%
for i_1_stage=wavstart_i_1:length(path_1_stage)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%一.2.1分层级完成音频文件在文件夹中的查询%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(isequal(path_1_stage(i_1_stage).name,'.')||isequal(path_1_stage(i_1_stage).name, '..'))        
        continue;
    elseif ~(path_1_stage(i_1_stage).isdir)               % 如果不是文件夹则跳过
        continue;
    end
    if overlapmethod==1
        if ~isempty(strfind(path_1_stage(i_1_stage).name,houzuiming_audio_WAV))
            wav_file=path_1_stage(i_1_stage).name;
            audiopath=fullfile(path_basic_0_stage,'\');
            i_2_stage=0;i_3_stage=0;i_4_stage=0;
        elseif  path_1_stage(i_1_stage).isdir
            path_1_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\' );  %,{'*.m4a';'*.wav'}
            path_2_stage = dir(path_1_file);
            for i_2_stage = 1 : length(path_2_stage) 
                if ~isempty(strfind(path_2_stage(i_2_stage).name,houzuiming_audio_WAV))
                    wav_file=path_2_stage(i_2_stage).name;
                    audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\');
                    i_3_stage=0;i_4_stage=0;
                elseif  path_2_stage(i_2_stage).isdir
                    path_2_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,'\');
                    path_3_stage = dir(path_2_file);  
                    for i_3_stage = 1 : length(path_3_stage) 
                        if ~isempty(strfind(path_3_stage(i_3_stage).name,houzuiming_audio_WAV))
                            wav_file=path_3_stage(i_3_stage).name;
                            audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,'\');
                            i_4_stage=0;                            
                        elseif  path_3_stage(i_3_stage).isdir
                            path_3_stage(i_3_stage).name
                            path_2_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,path_3_stage(i_3_stage).name,'\');
                            path_4_stage = dir(path_2_file);
                            for i_4_stage = 1 : length(path_4_stage) 
                                if ~isempty(strfind(path_4_stage(i_4_stage).name,houzuiming_audio_WAV))
                                    wav_file=path_4_stage(i_4_stage).name;
                                    audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage( i_2_stage ).name,'\' );
                                else
                                    continue
                                end
                            end                        
                        end                                        
                    end
                end
            end
        end
    elseif overlapmethod==2
        path_1_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\' );  %,{'*.m4a';'*.wav'}
        %addpath(path_1_file);
        path_2_stage = dir(path_1_file);
        for i_2_stage = wavstart_i_2:length(path_2_stage) 
            if( isequal(path_2_stage(i_2_stage).name,'.')||isequal(path_2_stage(i_2_stage).name, '..'))        
                continue;
            elseif ~(path_2_stage(i_2_stage).isdir)               % 如果不是文件夹则跳过
                continue;
            end
            if ~isempty(strfind(path_2_stage(i_2_stage).name,houzuiming_audio_WAV)) || ~isempty(strfind(path_2_stage(i_2_stage).name,houzuiming_audio_m4a))
                wav_file=audiopath(i_2_stage).name;
                audio.recorddate=path_2_stage(i_2_stage).date;
                audio.datpath=audiopath;
                audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\');
                i_3_stage=0;
            else                            
                path_2_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,'\');
                path_3_stage = dir(path_2_file);                     
                for i_3_stage = wavstart_i_3:length(path_3_stage) 
                    if( isequal( path_3_stage(i_3_stage).name,'.')||isequal(path_3_stage(i_3_stage).name, '..'))
                        continue;
                    elseif (path_3_stage(i_3_stage).isdir)               % 如果不是文件夹则跳过
                        continue;
                    %elseif isempty(strfind(path_3_stage(i_3_stage).name,houzuiming_audio_WAV)) && isempty(strfind(path_3_stage(i_3_stage).name,houzuiming_audio_m4a))
                    %    continue
                    end
                    path_basic_0_stage;
                    path_1_stage(i_1_stage).name;
                    path_2_stage(i_2_stage).name;
                    path_3_stage(i_3_stage).name;
                    %isempty(strfind(path_3_stage( i_3_stage ).name,'.WAV'))
                    audiopath = fullfile(path_basic_0_stage,path_1_stage( i_1_stage ).name,path_2_stage( i_2_stage ).name,'\');
                    if ~isempty(strfind(path_3_stage( i_3_stage ).name,houzuiming_audio_WAV))  ||    ~isempty(strfind(path_3_stage( i_3_stage ).name,houzuiming_audio_m4a))         
                        wav_file = path_3_stage(i_3_stage).name;
                        audio.recorddate=path_3_stage(i_3_stage).date;
                        audio.datpath=audiopath;                    
                    elseif ~isempty(strfind(path_3_stage( i_3_stage ).name,houzuiming_audio_txt))    
                        txt_file = path_3_stage(i_3_stage).name;
                        audio_str.recorddate=path_3_stage(i_3_stage).date;
                        audio_str.datpath=audiopath;
                    else
                        continue
                    end     
                    %{
                    if( isequal(wav_file,'.' )||isequal(wav_file,'..')|| path_3_stage(i_3_stage).isdir)  % 如果不是目录则跳过
                    continue         
                    end        
                    if isempty(strfind(wav_file,'.WAV')) && isempty(strfind(wav_file,'.m4a')) && isempty(strfind(wav_file,'.mp3') )
                    continue         
                    end
                    %}
                end
            end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%一.2.1分层级完成音频文件在文件夹中的查询完成%%%%%%%%
    %%%%%%%一.2.2音频文件计算解析开始 %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            % 一.2.2.1此处添加同名txt文件的查阅与读写操作 %
            %%%%path %将文件纳入搜索路径完成，包含音频与文件两种格式的信息%%%%
            %%%如果有同名文档，则写入文档%%%   
            if txtfileclean==1
                %%%%文档目录读取%%%%
                file_list=dir([audiopath,'*.txt']);
                filenum=size(file_list,1);%总共有多少个txt文件
                txt_file = file_list(iii).name; 
                %%%%判断语音文件与信息文件是否一致%%%%    
                txt_filestr=txt_file(1:strfind(txt_file,'.')-1);
                wav_filestr=wav_file(1:strfind(wav_file,'.')-1);
                coupled=strcmp(wav_filestr,txt_filestr);
                textorder=[];
                audioorder=[];   
                if coupled==0  
                    lookingforthistxtfile=strrep(wav_file,houzuiming_audio_WAV,'.txt');
                    lookingforthiswavfile=strrep(txt_file,'.txt',houzuiming_audio_WAV);
                    fileloc=find(ismember({file_list.name},lookingforthistxtfile ),1);
                    wavloc=find(ismember({audio_list.name},lookingforthistxtfile ),1);
                    if sum(fileloc)>0            
                        txt_file = file_list(fileloc).name; 
                    end
                    if sum(wavloc)>0            
                        wav_file = audio_list(wavloc).name; 
                    end      
                    if and (sum(fileloc)==0,sum(wavloc)==0)           
                        continue
                    end
                end
                %%%%判断语音文件与信息文件是否一致完成%%%%  
                %%%%先清空上次的记录%%%%
                [Txt_Id,message] = fopen(namethedatabase, 'w+');
                vain=[];
                fprintf(Txt_Id,'%s',vain);
                fclose(Txt_Id);
                %%%%清空结束%%%%                
            end
            %%%%文档目录读取完成%%%%
            % 一.2.2.1添加同名txt文件的查阅与读写操作完成 %
            %%%%%一.2.2.2显示当前所在的目录阶梯位置%%%%%%%%%%%
            currentsage=[i_1_stage,i_2_stage,i_3_stage] %
            if inGUI==1
                thecondition=sprintf('当前样本数进展到%d,%d,%d。',currentsage);
                set(handles.text_show,'string',thecondition)
                pause(0.05)
            end            
            %wav格式%%%%%%%%%%%%%%%%
            %%%%%一.2.2.2显示当前所在的目录阶梯位置完成%%%%%%%%%%%
            %%%%%一.2.2.3解析音频文件%%%%%%%%%%%
            wav_file   
            audiopath
            isempty(strfind(wav_file,houzuiming_audio_WAV))
            cd(audiopath)
            if ~isempty(strfind(wav_file,houzuiming_audio_WAV))
                [speech_all,fs_ori] =audioread(wav_file,'double');%H6记录仪器中，fs=96000Hz,超高频频段分析
                fs_ori
            elseif ~isempty(strfind(wav_file,houzuiming_audio_pcm))
                %pcm格式%%%%%%%%%%%%%%%%%%%%  
                WavId = fopen(wav_file,'r');               
                speech = fread(WavId,inf,'int16');
                speech=speech/(65536*0.5);  %2的16次方65536
                size(speech)                    
                %if n>=2
                %    speech=sqrt(speech(:,1).*speech(:,1)+speech(:,2).*speech(:,2));
                %end
            end        
            cd(thisfolder)
            audio.wavname=wav_file;
            %%%%%%%确保竖列%%%%%%%%%
            [m,n]=size(speech_all);
            if m<n
                speech_all=speech_all';%换成竖列
            end
            [m,n]=size(speech_all)  %音频大小与音轨数量显示
            if (m/fs_ori)>60
                %continue %腾腾呼吸科样本有一些很大的需要排除
            end
            numel(find((speech_all(:,1)<0)))
            numel(find((speech_all(:,1)>0)))
            if n>=2
            numel(find((speech_all(:,2)<0)))
            numel(find((speech_all(:,2)>0)))
            end
            speech_1=speech_all(:,1); %形成独立纵列数据
            time_end=length(speech_1)/fs_ori;
            timeline = linspace(0,time_end,numel(speech_1)); % 时间向量
            %%%%%%%%%只取一列完成%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%一.2.2.3解析音频文件并单列化完成%%%%%%%%%%%
            %%%%%一.2.2.4压缩音频文件，128 or perbit数量合并在一起，用于压缩%%%%%%%%%%%
            %%%%%%%%%为了防止vad识别的语音片段过小，现设计第二种补偿方式，即按照幅值划分。
            %            
            Ts=10; %ms 压缩点比特对应时长     
            perbit=floor(fs_ori*Ts/1000);
            tail=rem(size(speech_1,1),perbit);
            %size(speech_1,1)
            speech_1=speech_1(1:end-tail);
            %size(speech_1,1)
            speech_1_brifes=abs(reshape(speech_1,perbit,[]));%全部转化为正数
            speech_1_brife=sum(speech_1_brifes,1);%纵列压缩数据
            [speech_1_brife_m,speech_1_brife_n]=size(speech_1_brife)
            %time_brief=interp1(0:length(speech_1)/fs_ori,speech_1_brife,speech_1_brife_n)
            timeline_brife=linspace(0,length(speech_1)/fs_ori,speech_1_brife_n);            
            %%%%%一.2.2.4压缩音频文件，128 or perbit数量合并在一起，用于压缩完成%%%%%%%%%%%
            %%%%%一.2.2.5压缩信号的平滑滤波方法启动%%%%%%%%%%%
            if filtermethod==1
                % 小波函数降噪%%%%
                % 设置小波函数和变换阶数
                wname = 'db4';  % 选用 Daubechies 4 小波
                level = 5;      % 小波变换的阶数                
                [C, L] = wavedec(speech_1_brife, level, wname); % 进行小波变换               
                D = detcoef(C, L, level); % 提取细节系数
                % 对细节系数进行阈值处理
                sigma = median(abs(D)) / 0.6745;  % 计算阈值
                D = wthresh(D, 'h', sigma);       % 硬阈值处理
                % 重构滤波信号
                speech_1_brife_denoised = wrcoef('a', C, L, wname, level);
            elseif filtermethod==2 % 带通滤波  
                %fs_ori = 1000; % 采样频率  
                timeline_brife = linspace(0,time_end,numel(speech_1_brife)); % 时间向量  
                f_band_low = 0.1; % 低通频率  
                f_band_high = 10; % 高通频率  
                % 生成一个包含不同频率成分的信号  
                %speech_1_brife %= sin(2*pi*f1*t) + sin(2*pi*f2*t) + 0.5*sin(2*pi*50*t) + 0.5*sin(2*pi*300*t);  
                % 设计带通滤波器参数  
                Wn = [f_band_low/(fs_ori/2) f_band_high/(fs_ori/2)]; % 归一化频率  
                Rp = 1; % 通带最大衰减（以分贝为单位）  
                Rs = 30; % 阻带最小衰减（以分贝为单位） 
                % 使用butterworth滤波器设计带通滤波器  
                [b,a] = butter(5,Wn,'bandpass')
                % 应用滤波器  
                speech_1_brife_denoised = filter(b,a,speech_1_brife);  
                % 绘制原始信号和滤波后的信号
            elseif filtermethod==3 % 平均滤波
                lvbochidu=0.01;
                windowSize = ceil(numel(speech_1_brife)*lvbochidu);
                b = (1/windowSize)*ones(1,windowSize);
                a=1;
                speech_1_brife_denoised = filter(b,a,speech_1_brife);  
            elseif filtermethod==4 % 平滑 
                speech_1_brife_denoised = smooth(1:speech_1_brife_n,speech_1_brife,lvbochidu,'rloess');
            end
            % 绘制结果
            if figureshow==1
                figure(1)
                subplot(2,1,1);  
                plot(timeline_brife,speech_1_brife);  
                title('原始简化信号');  
                xlabel('时间 (s)');  
                ylabel('幅度'); 
                subplot(2,1,2);  
                plot(timeline_brife,speech_1_brife_denoised);  
                title('滤波后的简化信号');  
                xlabel('时间 (s)');  
                ylabel('幅度');
            end            
            speech_1_brife_smoothed = smooth(1:speech_1_brife_n,speech_1_brife,0.1,'rloess');
            basic_sec_num=13;  %H6实验过程中要调整到13,孟孟的实验调到6
            %%%%%一.2.2.5压缩信号的平滑滤波方法完成%%%%%%%%%%%
            %%%%%一.2.2.6识别基础语音位置方法选择开始%%%%%%%%%%%
            if basic_sec_num_method==1 %1为
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%下面这个for循环是为了找到合适的差值%%%%%%%%%%
                for i=2:7
                    [idx,C]=kmeans(abs(speech_1_brife_denoised),i);
                    %idx=[1 1 2 2 1 1 1 2 2 1 2 1 1 1 1 2 2 3 3 3 3 2 2 1 1 1 2 3 2 1 1 2 2 1 1 1 1 1 1 2 2 3];
                    size(idx)
                    [val loc]=max(C);
                    loc_voice=(find(idx==loc));%第一步，先定义每一段等于特定数字（最大值即语音部分）的连续数的位置；
                    %然后观察其比例
                    loc_1_dis=diff(loc_voice)%每一段等于特定数字的连续数的位置后面减去前面，得到差值，大于一表明其中连续；
                    loc_1_dis_all=cat(1,loc_voice(1),loc_1_dis) %与位置号码loc_1一一对应的间隔数字
                    %find(loc_1_dis>1,1)%第一个
                    breakplace_1=cat(1,find(loc_1_dis>1,1),diff(find(loc_1_dis>1))) %每一段等于特定数字的连续数的独立长度；
                    %%%看看是否平均%%%%%%
                     [value loc]=max(breakplace_1)
                     breakplace_1_nomaxrate=value/(sum(breakplace_1)-value)*(numel(breakplace_1)-1) %第一个指标，不能过大，防止局部化，设定阈值为2，
                    if numel(breakplace_1)<basic_sec_num
                        continue
                    end
                    if breakplace_1_nomaxrate<3
                        continue
                    else
                        break
                    end                    
                end
            elseif basic_sec_num_method==2 %1为
                baseline=abs((mean(speech_1_brife_denoised)-std(speech_1_brife_denoised)))
                speech_1_brife_denoised_cross=abs(speech_1_brife_denoised-baseline);
                crossplace_up=zeros(basic_sec_num,1);
                crossplace_down=zeros(basic_sec_num,1);
                u=1;uu=1;
                for iop=1:numel(speech_1_brife_denoised)-1
                    if speech_1_brife_denoised(iop)<=baseline && speech_1_brife_denoised(iop+1)>=baseline
                        crossplace_up(u)=iop;
                        u=u+1;
                    end
                    if speech_1_brife_denoised(iop)>=baseline && speech_1_brife_denoised(iop+1)<=baseline
                        crossplace_down(uu)=iop;
                        uu=uu+1;
                    end
                end
                crossplace_up=crossplace_up(crossplace_up>0);
                crossplace_down=crossplace_down(crossplace_down>0);
                %minvalue=min(speech_1_brife_denoised_cross)
                %row=find(speech_1_brife_denoised_cross<minvalue*3)
                %row=find(speech_1_brife_denoised<(mean(speech_1_brife_denoised)-std(speech_1_brife_denoised)));%[row,~,~]
                if crossplace_up(1)<crossplace_down(1)
                    crossplace_up=[1;crossplace_up];
                    crossplace_down=[crossplace_up;numel(speech_1_brife_denoised)];
                end
                crossplace_up
                crossplace_down
                if numel(crossplace_up)>numel(crossplace_down)
                    crossplace_up= crossplace_up(1:numel(crossplace_down));
                elseif numel(crossplace_up)<numel(crossplace_down)
                    crossplace_down= crossplace_down(1:numel(crossplace_up));
                end
                voiceplace_1=crossplace_down-crossplace_up;
                breakplace_1=[crossplace_up(1); crossplace_up-crossplace_down];
                
                
               % [loc value]=findpeaks(-speech_1_brife_denoised);
               % [loc1 value1]=findpeaks(value);
                %loc2=intersect(loc(loc1),row);
            if figureshow==1
                figure(2)
                subplot(2,1,1);  
                plot(timeline_brife,speech_1_brife_denoised);  
                hold on
                plot(timeline_brife(crossplace_up),speech_1_brife_denoised(crossplace_up),'r*');  
                plot(timeline_brife(crossplace_down),speech_1_brife_denoised(crossplace_down),'g*');  
                hold off
                title('原始简化信号');xlabel('时间 (s)');ylabel('幅度'); 
                subplot(2,1,2);  
                plot(timeline_brife,speech_1_brife_denoised);  
                title('滤波后的简化信号');xlabel('时间 (s)');ylabel('幅度');
            end           
                breakplace_1
                %findpeaks
            end             
            %%%%%一.2.2.6识别基础语音位置方法选择完成%%%%%%%%%%%
            %%%%%一.2.2.7不清楚是什么开始%%%%%%%%%%%
            %loc_2=(find(idx==2));%每一段等于特定数字的连续数的位置
            %loc_2_dis=diff(loc_2);%每一段等于特定数字的连续数的位置后面减去前面，得到差值，大于一表明其中连续；
            %breakplace_2=[find(loc_2_dis>1,1) diff(find(loc_2_dis>1))] %每一段等于特定数字的连续数的独立长度；
            %if loc_1(1)<loc_2(1)  %如果噪音在前方
            if unfinished==1
                breakplace_1_sorted=sort(breakplace_1,'descend');
                ori_order=[1:length(breakplace_1)]';
                size(breakplace_1)
                size(ori_order)
                breakplace_1_2=cat(2,breakplace_1,ori_order);        
                breakplace_1_2_sorted=sortrows(breakplace_1_2,-1);  %di
                order=[1:numel(breakplace_1)]';
                size(breakplace_1_2_sorted(1:end,:))
                size(order)
                breakplace_1_2_back=[breakplace_1_2_sorted(1:end,:) order];
                breakplace_1_2_back=[sortrows(breakplace_1_2_back,2) order]  %完工，列1时序排列幅值。2，处于所有频域峰值的位置；3幅值排名；4位置排名
                audio.maxEF_rhythmheight=breakplace_1_2_back(:,1);  %音律峰最大5项高度值时间顺序排列
                audio.maxEF_rhythmrank=breakplace_1_2_back(:,3); %音律峰最大5项高度时间顺序排名
                audio.maxEF_rhythmfre=(breakplace_1_2_back(:,2)); %音律峰最大5项频率时间顺序   
                audio.maxEF_rhythmnum=numel(maxEF_value); %音律共振峰数量
                audio.maxEF_valuemean=mean(maxEF_value);%语音节奏频率平均能量
                lastEF_value=maxEF_value(end);
                lastEF_fre=maxEF_location(end);
                audio.lastEF_value=lastEF_value;%最末梢语调节律幅值
                audio.lastEF_fre=lastEF_fre;%最末梢语调节律频率

                break_length_threathhold=breakplace_1_sorted(basic_sec_num);
                loc=zeros(length(idx),1);
                line=ones(break_length_threathhold+1,1);
                for i0=1:length(idx)-break_length_threathhold-1
                    if (idx(i0:i0+break_length_threathhold)==line)
                        loc(i0:i0+break_length_threathhold+1)=1;
                    end
                end
                loc'
            end
           %%%%%一.2.2.7不清楚是什么完成%%%%%%%%%%%
            %speech_all(1:20:40000,:)
            %abs(speech_all(1:20:4000,:))
            %ssssss
            %%%%%一.2.2.8声纹识别开始%%%%%%
            %=spectrogram(speech_all
            %%%%%一.2.2.8声纹识别完成%%%%%%%%%%%%
            %%%%%%%一.2.2.9聆听与否%%%%%%%%%
            if iwanthearit==1
                sound(speech_1(:,1),fs_ori); 
            end
            %%%%%%%一.2.2.9聆听与否结束%%%%%%%%%
            %}   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%一.2.2音频文件计算解析完成 %%%%%%%%
    %%%%%%%一.2.3音频特征计算开始%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%%一.2.3.1语音信号提纯与时长计算%%%%%
            speech_1(isnan(speech_1))=[];
            time_end=length(speech_1)/fs_ori;
            %%%%归一化%%%%
            audio.time=time_end;  %语音时长
            %speech=speech/(65536*0.5);
            %%%%第一步，重复完成语音发生部位的识别设计marktimeline完成标注%%%%    
            opts=1;
            vadThres=0.8;        
            [markFBEtimeline]=vad(speech_1,fs_ori,opts,vadThres);  %opts模式选择，vadthres，严苛度，0-1越高越严
            markFBEtimeline=reshape(markFBEtimeline,[],1); %纵列矩阵  
            voiceFBEstart_points=find(diff(markFBEtimeline)==1)+1;
            for i=1:9
                numel(voiceFBEstart_points)
                if numel(voiceFBEstart_points)>10
                    continue
                else
                    vadThres=0.95-i*0.1
                    [markFBEtimeline]=vad(speech_1,fs_ori,opts,vadThres);  %opts模式选择，vadthres，严苛度，0-1越高越严
                    markFBEtimeline=reshape(markFBEtimeline,[],1); %纵列矩阵  
                    voiceFBEstart_points=find(diff(markFBEtimeline)==1)+1;                    
                end
            end
            if numel(voiceFBEstart_points)<10
                continue
            end
            %diff(markFBEtimeline);
            voiceFBEstart_points=find(diff(markFBEtimeline)==1)+1;
            voiceFBEstop_points=find(diff(markFBEtimeline)==-1);
            if voiceFBEstart_points(end)>voiceFBEstop_points(end)
                voiceFBEstop_points=cat(1,voiceFBEstop_points,length(uuy));
            end
            if voiceFBEstart_points(1)>voiceFBEstop_points(1)
                voiceFBEstart_points=cat(1,1,voiceFBEstart_points);
            end
            %voiceFBEstart_points
            %voiceFBEstop_points
            for i=1:numel(voiceFBEstart_points)-1
                if voiceFBEstart_points(i+1)-voiceFBEstop_points(i)<3
                    voiceFBEstart_points(i+1)=0;
                    voiceFBEstop_points(i)=0;
                end
            end
            voiceFBEstart_points=voiceFBEstart_points(voiceFBEstart_points>0);
            voiceFBEstop_points=voiceFBEstop_points(voiceFBEstop_points>0);

            %%%%第二步，完成语义识别并形成文本结果%%%%            
            wav_file;
            if iwantcontent==1
                path_and_audio_file=strcat(audiopath,wav_file);
                speech_content=py_method(path_and_audio_file);
                speech_content=(char(speech_content));
            else
                speech_content='0_0';
            end
            loc_dot=strfind(speech_content,',');
            recognizedcharactersnum=length(speech_content)-numel(loc_dot);
            
            audio.recognizedcharactersnum=recognizedcharactersnum;
            %cd(thisfolder)
            %%%%第二步，完成语义识别并形成文本结果完成%%%%
                    %%%%第三步，完成声纹识别并检测文本的说话者%%%%  

            %%%%mfcc计算%%%%
            % Define variables
            Tw = 25;                % analysis frame duration (ms)
            Ts = 10;                % analysis frame shift (ms)
            alpha = 0.97;           % preemphasis coefficient
            M = 35;                 % number of filterbank channels, origially is 20
            C = 12;                 % number of cepstral coefficients
            %C = 19;                  %想让输出的两个数据格式相同
            L = 22;                 % cepstral sine lifter parameter
            LF = 30;               % lower frequency limit (Hz)
            HF = fs_ori/2;              % upper frequency limit (Hz)
            %wav_file = 'sp10.wav';  % input audio filename
            % Feature extraction (feature vectors as columns)
            [MFCCs,FBEs,frames,MAG,nfft,quitHz] = ...
                            mfcc( speech_1, fs_ori, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );
            size(MFCCs); %13,3102
            size(FBEs);  %M 3102
            size(frames); %2400,3102
            size(MAG); % 4096        3102
            size(nfft); %1 1 4096
            size(quitHz);%1    20
            nfft;%4096
            quitHz';%0.0033-4.8000
            

            % Generate data needed for plotting 
            [ Nw, NF ] = size( frames );                % frame length and number of frames
            time_frames = (0:NF-1)*Ts*0.001+0.5*Nw/fs_ori;  % time vector (s) for frames 
            time = (0:length(speech_1)-1)/fs_ori;           % time vector (s) for signal samples 
            logFBEs = 20*log10( FBEs );                 % compute log FBEs for plotting
            logFBEs_floor = max(logFBEs(:))-50 ;       % get logFBE floor 50 dB below max
            %size(logFBEs_floor);
            %size(logFBEs);
            logFBEs( logFBEs<logFBEs_floor ) = logFBEs_floor; % limit logFBE dynamic range

        %
        %%%%%%%%%希望通过对FBEs的分析，评估得出自动识别排除其时间段的方法。
            %%%%原始经典mel频率三图%%%%
            % Generate plots
            figure(1);%('Position', [30 30 800 600], 'PaperPositionMode', 'auto', ... 
                   %   'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 
            subplot( 311 );
            plot( time, speech_1, 'k' );
            xlim( [ min(time_frames) max(time_frames) ] );
            xlabel( 'Time (s)' ); 
            ylabel( 'Amplitude' ); 
            title( 'Speech waveform'); 
            subplot( 312 );
            imagesc( time_frames, (1:M), logFBEs ); 
            axis( 'xy' );
            xlim( [ min(time_frames) max(time_frames) ] );
            xlabel( 'Time (s)' ); 
            ylabel( 'Channel index' ); 
            title( 'Log (mel) filterbank energies'); 
            subplot( 313 );
            imagesc( time_frames, (1:C), MFCCs(2:end,:) ); % HTK's TARGETKIND: MFCC
            %imagesc( time_frames, [1:C+1], MFCCs );       % HTK's TARGETKIND: MFCC_0
            axis( 'xy' );
            xlim( [ min(time_frames) max(time_frames) ] );
            xlabel( 'Time (s)' ); 
            ylabel( 'Cepstrum index' );
            title( 'Mel frequency cepstrum' );
            % Set color map to grayscale
            colormap( 1-colormap('gray') ); 
            %%%%mfcc完成%%%%
            %%%%%%%%%%%2024新增加每一个发音单独分析的策略%%%%%%%%%%
            %%%%前期计算分析语音的能量曲线，过零率，短时能量分布基础性信息%%%%
            %framepointnum=floor(length(speech)/line);
            [~,line]=size(FBEs);
            crosszeroline=zeros(1,line);  %语音信号的过零率行列 
            framepointnum=floor(fs_ori*Ts/1000);
            speech_FBEs_energy_frameline=zeros(1,line);  %音频段能量线
            speechenergyrank=zeros(framepointnum,line);
            %speech_1=speech
            for j=1:line
                %speechenergyrank(j)= sum(speech(framepointnum*(j-1)+1:framepointnum(j)).*speech(framepointnum*(j-1)+1:framepointnum(j)));
                for k=1:framepointnum-1
                   speechenergyrank(:,j)= speech_1(framepointnum*(j-1)+k)*speech_1(framepointnum*(j-1)+k);%声音能量矩阵           
                    if speech_1(framepointnum*(j-1)+k)*speech_1(framepointnum*(j-1)+k+1)<0                
                        crosszeroline(j)=crosszeroline(j)+1;
                    end
                end
                speech_FBEs_energy_frameline(j)=sum(speechenergyrank(:,j));%声音能量线性，即时序声压
                %speech_intense_frameline(j)= speech_energy_frameline/2*(j);%声音强度，声强
            end
            %%%%前期计算分析语音的能量曲线，过零率，短时能量分布基础性信息完成%%%%
            %%%%分割音频中的每一次发音节点%%%%
            markspeechtimeline=repmat(markFBEtimeline',floor(fs_ori*Ts/1000),1);
            markspeechtimeline=reshape(markspeechtimeline,[],1);%变成一列数据标注speech。
            if length(markspeechtimeline)>length(speech_1)
                markspeechtimeline=markspeechtimeline(1:length(speech_1));
            elseif length(markspeechtimeline)<length(speech_1)
                speech_1=speech_1(1:length(markspeechtimeline));
            end
            speech_start_points=voiceFBEstart_points*floor(fs_ori*Ts/1000);
            speech_stop_points=voiceFBEstop_points*floor(fs_ori*Ts/1000);
            if speech_stop_points(end)>length(speech_1)
                speech_start_points=speech_start_points(1:end-1);
                speech_stop_points=speech_stop_points(1:end-1);                
            end
            voicepressure_rms=rms(speech_1(markspeechtimeline==1));  %均方根声压均值
            noisepressure_rms=rms(speech_1(markspeechtimeline==0)) ; %Root-mean-square level
            ref_pressure=20*10^(-6);
            voicepressurelevel_rms=20*log10(voicepressure_rms/ref_pressure);
            noisepressurelevel_rms=20*log10(noisepressure_rms/ref_pressure);
            %a_weighting
            speech_sec_cell=cell(numel(voiceFBEstart_points),1);
            speech_sec_mean_energy_20=zeros(numel(voiceFBEstart_points),1);
            speech_sec_std_20=zeros(numel(voiceFBEstart_points),1);
            charanum_upperlimit=20;
            tail=zeros(charanum_upperlimit,1);
            speech_sec_length_list=zeros(numel(voiceFBEstart_points),1);
            for i=1:numel(voiceFBEstart_points)
                speech_sec_cell{i}=speech_1(speech_start_points(i):speech_stop_points(i)); 
                speech_sec_mean_energy_20(i)=mean(abs(speech_sec_cell{i}));%
                speech_sec_std_20(i)=std(abs(speech_sec_cell{i}));
                speech_sec_length_list(i)=voiceFBEstop_points(i)-voiceFBEstart_points(i)+1;
            end
            speech_sec_inter_20=voiceFBEstart_points(2:end)-voiceFBEstop_points(1:end-1)+1;%
            speech_sec_inter_20=cat(1,voiceFBEstart_points(1),speech_sec_inter_20);%长度对齐
            if numel(voiceFBEstart_points)<charanum_upperlimit
                speech_sec_mean_energy_20=cat(1,speech_sec_mean_energy_20,tail);
                speech_sec_std_20=cat(1,speech_sec_std_20,tail);
                speech_sec_length_list_20=cat(1,speech_sec_length_list,tail);
                speech_sec_inter_20=cat(1,speech_sec_inter_20,tail);
                speech_sec_mean_energy_20=speech_sec_mean_energy_20(1:charanum_upperlimit);
                speech_sec_std_20=speech_sec_std_20(1:charanum_upperlimit);
                speech_sec_length_list_20=speech_sec_length_list_20(1:charanum_upperlimit);
                speech_sec_inter_20=speech_sec_inter_20(1:charanum_upperlimit);
            else
                speech_sec_mean_energy_20=speech_sec_mean_energy_20(1:charanum_upperlimit);
                speech_sec_std_20=speech_sec_std_20(1:charanum_upperlimit);
                speech_sec_length_list_20=speech_sec_length_list(1:charanum_upperlimit);
                speech_sec_inter_20=speech_sec_inter_20(1:charanum_upperlimit);
            end
            %%%%%%%%识别每一个声点的长度，少的做好空，多的做删减。%%%%%%%            
            %speech_sec_length_list
            chara_count=numel(speech_sec_length_list);%间隔长度的差异。
            [value,loc1_9doubt]=max(diff(speech_sec_length_list)); %根据其中发音延时长度的变化，找到数字到三元音之间的间隙
            %loc1_9doubt是间隔数的位数，表达其所处的位置
            %numel(chara_count)>9    
                speech_content                
                loc_jiu=strfind(speech_content,'九');
                if ~isempty(loc_jiu)
                    speech_content_1_9 = speech_content(1:loc_jiu(1));
                    loc_dot=strfind(speech_content_1_9,',');
                    recognized_the_nine=1;
                else
                    loc_jiu=strfind(speech_content,'八');                    
                    recognized_the_nine=0;
                    if ~isempty(loc_jiu)
                        speech_content_1_9 = speech_content(1:loc_jiu(1));
                        loc_dot=strfind(speech_content_1_9,',');
                    else
                        loc_jiu(1)=round(length(speech_content)/2);
                        speech_content_1_9 = speech_content(1:loc_jiu(1));     
                        loc_dot=strfind(speech_content_1_9,',');
                    end                    
                end
                if isempty(loc_dot)
                    loc_dot=0;
                end
                audio.recognized_the_nine=recognized_the_nine;
                if loc_jiu(1)>numel(loc_dot)
                    loc1_9doubt=loc_jiu(1)-numel(loc_dot);
                else
                    loc1_9doubt=9;
                    if chara_count<9
                        loc1_9doubt=chara_count-3;
                    end
                end
                count1_9loss=9-loc1_9doubt;
                count1_9_num=loc1_9doubt;
            if  loc1_9doubt<=9                
                count1_9loss=9-loc1_9doubt;
            elseif loc1_9doubt>9
                count1_9loss=0;
            end
            if chara_count-loc1_9doubt<=4
                count_aaaongmmm_num=chara_count-loc1_9doubt-1;
                count_cough_num=1;
            else
                count_aaaongmmm_num=3;
                count_cough_num=chara_count-loc1_9doubt-3;
            end
            %chara_count-loc1_9doubt
            %count_aaaongmmm_num
            %count_cough_num
            threecount=[loc1_9doubt count_aaaongmmm_num count_cough_num]
            %speech_stop_points
            aaaongmmm_stopnum=chara_count-count_cough_num;
            speech_1_1_9=speech_1(1:speech_stop_points(loc1_9doubt));
            speech_1_2_aaaongmmm=speech_1(speech_start_points(loc1_9doubt+1):speech_stop_points(aaaongmmm_stopnum));
            speech_1_3_cough=speech_1(speech_start_points(aaaongmmm_stopnum)+1:speech_stop_points(chara_count));
            if size(speech_1_3_cough,1)<3
                speech_1_2_aaaongmmm=speech_1(speech_start_points(loc1_9doubt+1):speech_stop_points(chara_count)-3);
                speech_1_3_cough=speech_1(speech_start_points(chara_count)-3+1:speech_stop_points(chara_count));
            end
            [length_1_9,~]=size(speech_1_1_9)
            [length_aaaongmmm,~]=size(speech_1_2_aaaongmmm)
            [length_cough,~]=size(speech_1_3_cough)
            line=size(FBEs,2)
            framepointnum=floor(fs_ori*Ts/1000);
            aongm_num=floor(length_aaaongmmm/framepointnum);
            cough_num=floor(length_cough/framepointnum);
            %fs_ori*Ts/1000
            FBEs_1_9=FBEs(:,1:floor(length_1_9/(fs_ori*Ts/1000)));
            frames_1_9=frames(:,1:floor(length_1_9/(fs_ori*Ts/1000)));
            MAG_1_9=MAG(:,1:floor(length_1_9/(fs_ori*Ts/1000)));
            markFBEtimeline_1_9=markFBEtimeline(1:floor(length_1_9/(fs_ori*Ts/1000)));
            %%%%%%%%%%%%%%%%%%
            FBEs_aongm=FBEs(:,floor(length_1_9/(fs_ori*Ts/1000))+1:floor(length_1_9/(fs_ori*Ts/1000))+aongm_num);
            frames_aongm=frames(:,floor(length_1_9/(fs_ori*Ts/1000))+1:floor(length_1_9/(fs_ori*Ts/1000))+aongm_num);
            MAG_aongm=MAG(:,floor(length_1_9/(fs_ori*Ts/1000))+1:floor(length_1_9/(fs_ori*Ts/1000))+aongm_num);
            markFBEtimeline_aongm=markFBEtimeline(floor(length_1_9/(fs_ori*Ts/1000))+1:floor(length_1_9/(fs_ori*Ts/1000))+aongm_num);
            
            FBEs_cough=FBEs(:,line-cough_num+1:end);
            frames_cough=frames(:,line-cough_num+1:end);
            MAG_cough=MAG(:,line-cough_num+1:end);
            markFBEtimeline_cough=markFBEtimeline(line-cough_num+1:end);
            FBEs_sec_length=[size(FBEs_1_9,2) size(FBEs_aongm,2) size(FBEs_cough,2)]
            if size(FBEs_1_9,2)<50 || size(FBEs_aongm,2)<50 ||size(FBEs_cough,2)<50
                continue
            end
            %%%%%%%%分样本audio计算分析开始%%%%%%%%
            %markFBEtimeline_1_9size=size(markFBEtimeline_1_9)
            %markFBEtimeline_aongmsize=size(markFBEtimeline_aongm)
            %markFBEtimeline_coughsize=size(markFBEtimeline_cough)
            [audio_1_9,audio_t_1_9,audio_dat_1_9]=data_generaterH6(speech_1_1_9,audio,fs_ori,Ts,FBEs_1_9,markFBEtimeline_1_9,MAG_1_9,quitHz);  
            [audio_aongm,audio_t_aongm,audio_dat_aongm]=data_generaterH6(speech_1_2_aaaongmmm,audio,fs_ori,Ts,FBEs_aongm,markFBEtimeline_aongm,MAG_aongm,quitHz);  
            [audio_cough,audio_t_cough,audio_dat_cough]=data_generaterH6(speech_1_3_cough,audio,fs_ori,Ts,FBEs_cough,markFBEtimeline_cough,MAG_cough,quitHz);  
            %audio_1_9
            %audio_t_1_9
            %audio_dat_1_9
            
            
            %%%%%%%%audio计算分析结束%%%%%%%%
            
            
            audio.voicepressure_rms=voicepressure_rms;%语音声压均值均方根
            audio.noisepressure_rms=noisepressure_rms;%噪音声压均值均方根
            audio.voicepressurelevel_rms=voicepressurelevel_rms;%语音声压级均值均方根
            audio.noisepressurelevel_rms=noisepressurelevel_rms;%噪音声压级均值均方根
            
            audio.count1_9loss=count1_9loss;%1-9识别疏漏
            audio.count1_9_num=count1_9_num;%1-9过度识别点位
            audio.count_aaaongmmm_num=count_aaaongmmm_num;%aaaongmmm识别疏漏
            audio.count_cough_num=count_cough_num;%咳嗽音数量
            audio.charanumall=numel(voiceFBEstart_points);
            audio.speech_sec_mean_energy_20=speech_sec_mean_energy_20;%单字均能序列
            audio.speech_sec_std_20=speech_sec_std_20;%单字音长序列
            audio.speech_sec_length_list_20=speech_sec_length_list_20;
            audio.speech_sec_inter_20=speech_sec_inter_20;
                
            %p=ones(1,10)*(1/10);
            %speech_FBEs_energy_frameline=filter(p,1,speech_FBEs_energy_frameline);%第一步，平滑完成；
               % voicevolume_average=max(abs(speech_FBEs_energy_frameline(markFBEtimeline==1)));%/(noisetimeratio*length(speech)/fs); %平均语音能量

            
           %%%%%%%%%%%2024新增加每一个发音单独分析的策略完成%%%%%%%%%%

            %%%%%%%%audio整体计算分析开始%%%%%%%%
            [audio,audio_t,audio_dat]=data_generaterH6(speech_1,audio,fs_ori,Ts,FBEs,markFBEtimeline,MAG,quitHz);  
            %%%%%%%%audio整体计算分析结束%%%%%%%%
            audio;
            audio_t;
            %%%%%%%%新增整体与部分数据%%%%%%%%
            audio_str.audio_str_3={audio.recorddate,audio.datpath,audio.wavname};
            audio_t.audio_str_3={'记录日期','路径','文件名'};
            before_sep_para_89={audio.voicepressure_rms;audio.noisepressure_rms;...
                audio.voicepressurelevel_rms;audio.noisepressurelevel_rms;...
                audio.count1_9loss;audio.count1_9_num;...
                audio.count_aaaongmmm_num;audio.count_cough_num;audio.charanumall;...
                audio.speech_sec_mean_energy_20;...    %单字均能序列
                audio.speech_sec_std_20;...%单字均能标准差序列
                audio.speech_sec_length_list_20;...%单字音长序列
                audio.speech_sec_inter_20;...%单字间隔序列
                };
            audio_dat.before_sep_para_89=before_sep_para_89;
            %%%%%%%%%%%
            addpara_10={'语音声压','噪音声压','语音声压级','噪音声压级','1-9丢失位',...
                '1-9识别数量','三元音识别数量','咳嗽音识别数量','总共识别数量'};
            numm=1:charanum_upperlimit;    
                %str1=repmat('mel20频谱峰和',numel(numm),1);
                str2=num2str(reshape(numm,[],1));                
                str3=repmat('s',numel(numm),1);
                titles=strcat('识别',num2str(charanum_upperlimit),'音节均能序列',str2,str3);
            title_speech_sec_mean_energy_20=mat2cell(titles,ones(1,numel(numm)))';                       
                titles=strcat('识别',num2str(charanum_upperlimit),'音节均能标准差序列',str2,str3);
            title_speech_sec_std_20=mat2cell(titles,ones(1,numel(numm)))';
                titles=strcat('识别',num2str(charanum_upperlimit),'音节音长序列',str2,str3);
            title_speech_sec_length_list_20=mat2cell(titles,ones(1,numel(numm)))';                       
                titles=strcat('识别',num2str(charanum_upperlimit),'音节间隔序列',str2,str3);
            title_speech_sec_inter_20=mat2cell(titles,ones(1,numel(numm)))';
            title_before_sep_para_89={addpara_10{1:end},...
                title_speech_sec_mean_energy_20{1:end},...
                title_speech_sec_std_20{1:end},...
                title_speech_sec_length_list_20{1:end},...
                title_speech_sec_inter_20{1:end},...
                };    
            audio_t.title_before_sep_para_89=title_before_sep_para_89;
            
            mainparameters=getappdata(0,'mainparameters');
            repeattime=mainparameters.repeattime;
            %mainparameters.overlapmethod=overlapmethod;
            %mainparameters.txtfiletosamewav=txtfiletosamewav;
            %mainparameters.overlapmethod=overlapmethod;
            %mainparameters.txtfiletosamewav=txtfiletosamewav;
            audio_t
            if TxtId > 0           
                if repeattime==1
                    fprintf(TxtId,'%s\t','一序号','二序号','三序号','识别汉字数量','是否识别了那个-九','文件夹姓名','音频名称','语义');
                    audio_t.audio_str_3={'记录日期','路径','文件名'};  %9
                    fprintf(TxtId,'%s\t',audio_t.audio_str_3{1:end});%9
                    fprintf(TxtId,'%s\t',audio_t.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice4_100{1:end});  
                    %%%%%%%%%%%%新加部分%%%%%%%%%%
                    
                    fprintf(TxtId,'%s\t',audio_t.title_before_sep_para_89{1:end});
                    fprintf(TxtId,'%s\t','1-9自动分割后的片段数据分析');
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice4_100{1:end});
                    fprintf(TxtId,'%s\t','aaa-ong-mmm自动分割后的片段数据分析');
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice4_100{1:end});
                    fprintf(TxtId,'%s\t','cough音自动分割后的片段数据分析');
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice4_100{1:end});
                end
                mainparameters.repeattime=repeattime+1;
                setappdata(0,'mainparameters',mainparameters);
                % Write the output data开始写入数据   
                fprintf(TxtId,'\n');
                fprintf(TxtId,'%5.0f\t',i_1_stage,i_2_stage,i_3_stage,audio.recognizedcharactersnum,audio.recognized_the_nine);%4 num
                fprintf(TxtId,'%s\t',path_1_stage(i_1_stage).name,wav_file,speech_content);%3 str
                fprintf(TxtId,'%s\t',audio_str.audio_str_3{1:end});%audio_dat.recorddate,audio_dat.datpath,audio_dat.wavname);%audio_str.audio_str_3
                       %11+2=13;对应第一条标题 
                fprintf(TxtId,'%10.4f\t',audio_dat.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice4_100{1:end});%  
                %%%%%%%%%%%%%%%%%%%%%                
                fprintf(TxtId,'%10.4f\t',audio_dat.before_sep_para_89{1:end});%  
                fprintf(TxtId,'%s\t','1-9自动分割后的片段数据分析');
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice4_100{1:end});%
                fprintf(TxtId,'%s\t','aaa-ong-mmm自动分割后的片段数据分析');
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice4_100{1:end});%
                fprintf(TxtId,'%s\t','cough音自动分割后的片段数据分析');
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice4_100{1:end});%
            else 
                h = msgbox('Record audio characters in files failed!');
            end
            %显示框的名字text_show
            if inGUI==1
                pause(0.05)
                thecondition=sprintf('当前样本数进展到%d,%d,%d。',currentsage)
                set(handles.text_show,'string',thecondition)
            end
        end
    end
end
fclose(TxtId);