function varargout = Voice_batch_GUI(varargin)
% VOICE_BATCH_GUI MATLAB code for Voice_batch_GUI.fig
%      VOICE_BATCH_GUI, by itself, creates a new VOICE_BATCH_GUI or raises the existing
%      singleton*.
%
%      H = VOICE_BATCH_GUI returns the handle to a new VOICE_BATCH_GUI or the handle to
%      the existing singleton*.
%
%      VOICE_BATCH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOICE_BATCH_GUI.M with the given input arguments.
%
%      VOICE_BATCH_GUI('Property','Value',...) creates a new VOICE_BATCH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Voice_batch_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Voice_batch_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Voice_batch_GUI

% Last Modified by GUIDE v2.5 23-Jan-2024 19:11:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Voice_batch_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Voice_batch_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Voice_batch_GUI is made visible.
function Voice_batch_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Voice_batch_GUI (see VARARGIN)

% Choose default command line output for Voice_batch_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Voice_batch_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
thecondition=sprintf('Voice analysis in standard method');
            set(handles.text_title,'string',thecondition,'FontSize',18);

% --- Outputs from this function are returned to the command line.
function varargout = Voice_batch_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%�����������ݣ�Ӧ���Լ�������%%%%
    repeattime=1;%ѭ������������д����
    overlapmethod=2;%3�������ļ��з�����1��2
    txtfiletosamewav=0; %�����ļ��У����û��ͬ����txt�ļ�ѡ0����ѡ1
    iwanthearit=0; %������ѡ0��������ѡ1
    iwantcontent=0;%����ʶ������ѡ0����ʶ�������ѡ1
        %%%%ѭ���������%%%%
    wavstart_i_1=1;
    wavstart_i_2=1;
    wavstart_i_3=1;
    namethedatabase='datafileport-2024_1.txt';
    uichoose=1
    if uichoose==1    
        prompt = {'\fontsize{11} 1.ѭ������������д������д1:','\fontsize{11}2.���������ļ��з�����1��2:',...
            '\fontsize{11} 3.���û��ͬ����txt�ļ�ѡ0����ѡ1:','\fontsize{11}4.������ѡ0��������ѡ1:',...
            '\fontsize{11}5.����ʶ������ѡ0����ʶ�������ѡ1:',...
            '\fontsize{11}6.һ���ļ����жϼ������ο�text��ʾ��Ĭ��ѡ1:',...
            '\fontsize{11}7.�����ļ����жϼ������ο�text��ʾ��Ĭ��ѡ1:',...
            '\fontsize{11}8.�����ļ����жϼ������ο�text��ʾ��Ĭ��ѡ1:',...
            '\fontsize{11}9.д�����ݵ��ļ����ƣ���׺����txt��ʽ:'};
        dlg_title = 'Input your information';
        num_lines = 1;
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        %def = {'1','2','0','0','1'};
        def = {num2str(repeattime),num2str(overlapmethod),num2str(txtfiletosamewav),...
            num2str(iwanthearit),num2str(iwantcontent),...
            num2str(wavstart_i_1),num2str(wavstart_i_2),num2str(wavstart_i_3),namethedatabase};
        %def = {P{1}{line},P{2}{line}, P{3}{line}};
        %def = {'����','19860722', '��'};
        mainparameterscell = inputdlg(prompt,dlg_title,num_lines,def,options);
        repeattime=str2double(mainparameterscell{1});
        overlapmethod=str2double(mainparameterscell{2});
        txtfiletosamewav=str2double(mainparameterscell{3});
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
    mainparameters.txtfiletosamewav=txtfiletosamewav;
    mainparameters.iwanthearit=iwanthearit;
    mainparameters.iwantcontent=iwantcontent;
    setappdata(0,'mainparameters',mainparameters);

    %
    %%%%path %���ļ���������·��%%%%
    %cd('H:\��Ƶ����\'); %��������ǰ����Ŀ¼
    %%%%��ȡ�����ļ�%%%%
    fivetonedata='five_tone_data.xlsx';
    [data,~,~]=xlsread(fivetonedata);
    fivetonemel20=data(:, 10:29);   %����ԭʼmelƵ��20����
    fivetonerank5=data(:, 1:3);
    max_five20=max(fivetonemel20,[],2);  %���������ֵ
    min_five20=min(fivetonemel20,[],2);  %��������Сֵ
    M=size(fivetonemel20,2);  %m=20
    maxrank=repmat(max_five20,1,M);  
    minrank=repmat(min_five20,1,M);
    max_five20=(fivetonemel20-minrank)./(maxrank-minrank);    
    %%%%��ȡ�����ļ�����%%%%
    clear speech time Yabsvoicefft Yfft Ynoisefft Yvoicefft markspeechtimeline speechenergyrank
    clear MFCCs logFBEs
    
    %%%%����֮ǰ���µ��ڴ�%%%%
    %%%%���ر�Ҫ���ļ�%%%%
    this_path = mfilename('fullpath'); %��m�ļ�����λ��
    [thisfolder,~,~]=fileparts(this_path);
    cd(thisfolder);   
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
    [TxtId, ~] = fopen(namethedatabase, 'a+');%��for������ĵ��������ظ�����%cd('..')
    %%%%���ر�Ҫ���ļ����%%%%
    %%%%ѡ����������Ƶ�ĵ�%%%%
    path_basic_0_stage=uigetdir();%'G:\����������\������������\';
    path_1_stage  = dir( path_basic_0_stage );    
    %%%%ѡ����������Ƶ�ĵ����%%%%

    %%%%��ʼ������Ƶ�ļ�%%%%
for i_1_stage=wavstart_i_1:length(path_1_stage)
    if(isequal(path_1_stage(i_1_stage).name,'.')||isequal(path_1_stage(i_1_stage).name, '..'))        
        continue;
    elseif ~(path_1_stage(i_1_stage).isdir)               % ��������ļ���������
        continue;
    end
    if overlapmethod==1
        if ~isempty(strfind(path_1_stage(i_1_stage).name,'.WAV'))
            wav_file=path_1_stage(i_1_stage).name;
            audiopath=fullfile(path_basic_0_stage,'\');
            i_2_stage=0;i_3_stage=0;i_4_stage=0;
        elseif  path_1_stage(i_1_stage).isdir
            path_1_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\' );  %,{'*.m4a';'*.wav'}
            path_2_stage = dir(path_1_file);
            for i_2_stage = 1 : length(path_2_stage) 
                if ~isempty(strfind(path_2_stage(i_2_stage).name,'.WAV'))
                    wav_file=path_2_stage(i_2_stage).name;
                    audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,'\');
                    i_3_stage=0;i_4_stage=0;
                elseif  path_2_stage(i_2_stage).isdir
                    path_2_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,'\');
                    path_3_stage = dir(path_2_file);  
                    for i_3_stage = 1 : length(path_3_stage) 
                        if ~isempty(strfind(path_3_stage(i_3_stage).name,'.WAV'))
                            wav_file=path_3_stage(i_3_stage).name;
                            audiopath=fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,'\');
                            i_4_stage=0;                            
                        elseif  path_3_stage(i_3_stage).isdir
                            path_3_stage(i_3_stage).name
                            path_2_file = fullfile(path_basic_0_stage,path_1_stage(i_1_stage).name,path_2_stage(i_2_stage).name,path_3_stage(i_3_stage).name,'\');
                            path_4_stage = dir(path_2_file);
                            for i_4_stage = 1 : length(path_4_stage) 
                                if ~isempty(strfind(path_4_stage(i_4_stage).name,'.WAV'))
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
            elseif ~(path_2_stage(i_2_stage).isdir)               % ��������ļ���������
                continue;
            end
            if ~isempty(strfind(path_2_stage(i_2_stage).name,'.WAV'))
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
                    elseif (path_3_stage(i_3_stage).isdir)               % ��������ļ���������
                        continue;
                    elseif isempty(strfind(path_3_stage(i_3_stage).name,'.WAV'))
                        continue
                    end
                    path_basic_0_stage;
                    path_1_stage(i_1_stage).name;
                    path_2_stage(i_2_stage).name;
                    path_3_stage(i_3_stage).name;
                    %isempty(strfind(path_3_stage( i_3_stage ).name,'.WAV'))
                    audiopath = fullfile(path_basic_0_stage,path_1_stage( i_1_stage ).name,path_2_stage( i_2_stage ).name,'\');
                    %addpath(audiopath);              
                    wav_file = path_3_stage(i_3_stage).name;
                    audio.recorddate=path_3_stage(i_3_stage).date;
                    audio.datpath=audiopath;
                    if( isequal(wav_file,'.' )||isequal(wav_file,'..')|| path_3_stage(i_3_stage).isdir)  % �������Ŀ¼������
                    continue         
                    end        
                    if isempty(strfind(wav_file,'.WAV')) && isempty(strfind(wav_file,'.m4a')) && isempty(strfind(wav_file,'.mp3') )
                    continue         
                    end
                end
            end                
                % �˴������Ķ��ļ���д���� %
                %%%%path %���ļ���������·����ɣ�������Ƶ���ļ����ָ�ʽ����Ϣ%%%%
                %%%�����ͬ���ĵ�����д���ĵ�%%%   
            if txtfiletosamewav==1
            %%%%�ĵ�Ŀ¼��ȡ%%%%
            file_list=dir([path_basic_0_stage,'*.txt']);
                %%%%�����%%%%
                [Txt_Id,message] = fopen('datafileport.txt', 'w+');
                vain=[];
                fprintf(Txt_Id,'%s',vain);
                fclose(Txt_Id);
                %%%%��ս���%%%%
                filenum=size(file_list,1);
            end
            %%%%�ĵ�Ŀ¼��ȡ���%%%%
            %%%%%��ǰ���ڵ�Ŀ¼����λ����ʾ%%%%%%%%%%%
            currentsage=[i_1_stage,i_2_stage,i_3_stage];%
            thecondition=sprintf('��ǰ��������չ��%d,%d,%d��',currentsage);
            set(handles.text_show,'string',thecondition)
            pause(0.05)
            %wav��ʽ%%%%%%%%%%%%%%%%
            %��Ƶ�ļ�������ʾ
            wav_file    
            cd(audiopath)
            [speech_all,fs_ori] =audioread(wav_file,'double');%H6��¼�����У�fs=96000Hz,����ƵƵ�η���
            cd(thisfolder)
            audio.wavname=wav_file;
            [m,n]=size(speech_all);
            if m<n
                speech_all=speech_all';%��������
            end
            [m,n]=size(speech_all)  %��Ƶ��С������������ʾ
            if (m/fs_ori)>60
                continue
            end
            numel(find((speech_all(:,1)<0)));
            numel(find((speech_all(:,1)>0)));
            numel(find((speech_all(:,2)<0)));
            numel(find((speech_all(:,2)>0)));
            speech_1=speech_all(:,1); 
            %%%%%%����ʶ��%%%%%%
            %=spectrogram(speech_all
            %%%%%%����ʶ�����%%%%%%
            %%%%%%%�������%%%%%%%%%
            if iwanthearit==1
                sound(speech_1(:,1),fs_ori); 
            end
            %%%%%%%����������%%%%%%%%%
            %}
            %pcm��ʽ%%%%%%%%%%%%%%%%%%%%
            %{
            WavId = fopen(wav_file,'r')
            audio.wavname=wav_file;   
            speech = fread(WavId,inf,'int16');
            speech=speech/(65536*0.5);  %2��16�η�65536
            size(speech)
            %sound(speech(:,1),fs);
            fs=16000;   
             %}
            %if n>=2
            %    speech=sqrt(speech(:,1).*speech(:,1)+speech(:,2).*speech(:,2));
            %end
            %WavId = fopen(wav_file,'r');
            %%%%�ı��ļ����ƶ�ȡ��Ϣ%%%%
            %{
            txt_file = file_list(iii).name; 
            %%%%�ж������ļ�����Ϣ�ļ��Ƿ�һ��%%%%    
            txt_filestr=txt_file(1:strfind(txt_file,'.')-1);
            wav_filestr=wav_file(1:strfind(wav_file,'.')-1);
            coupled=strcmp(wav_filestr,txt_filestr);
            textorder=[];
            audioorder=[];   
            if coupled==0  
                lookingforthistxtfile=strrep(wav_file,'.pcm','.txt');
                lookingforthiswavfile=strrep(txt_file,'.txt','.pcm');
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
            %%%%�ж������ļ�����Ϣ�ļ��Ƿ�һ�����%%%%  
            %}
            %%%%%%%%�����ź��ᴿ��ʱ������%%%%%
            speech_1(isnan(speech_1))=[];
            timespan=length(speech_1)/fs_ori;
            %%%%��һ��%%%%
            audio.time=timespan;  %����ʱ��
            %speech=speech/(65536*0.5);
            %%%%��һ�����ظ��������������λ��ʶ�����marktimeline��ɱ�ע%%%%    
            opts=1;
            vadThres=0.8;        
            [markFBEtimeline]=vad(speech_1,fs_ori,opts,vadThres);  %optsģʽѡ��vadthres���Ͽ��ȣ�0-1Խ��Խ��
            markFBEtimeline=reshape(markFBEtimeline,[],1); %���о���  
            voiceFBEstart_points=find(diff(markFBEtimeline)==1)+1;
            for i=1:9
                numel(voiceFBEstart_points)
                if numel(voiceFBEstart_points)>10
                    continue
                else
                    vadThres=0.95-i*0.1
                    [markFBEtimeline]=vad(speech_1,fs_ori,opts,vadThres);  %optsģʽѡ��vadthres���Ͽ��ȣ�0-1Խ��Խ��
                    markFBEtimeline=reshape(markFBEtimeline,[],1); %���о���  
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

            %%%%�ڶ������������ʶ���γ��ı����%%%%            
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
            %%%%�ڶ������������ʶ���γ��ı�������%%%%
                    %%%%���������������ʶ�𲢼���ı���˵����%%%%  

            %%%%mfcc����%%%%
            % Define variables
            Tw = 25;                % analysis frame duration (ms)
            Ts = 10;                % analysis frame shift (ms)
            alpha = 0.97;           % preemphasis coefficient
            M = 35;                 % number of filterbank channels, origially is 20
            C = 12;                 % number of cepstral coefficients
            %C = 19;                  %����������������ݸ�ʽ��ͬ
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
        %%%%%%%%%ϣ��ͨ����FBEs�ķ����������ó��Զ�ʶ���ų���ʱ��εķ�����
            %%%%ԭʼ����melƵ����ͼ%%%%
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
            %%%%mfcc���%%%%
            %%%%%%%%%%%2024������ÿһ���������������Ĳ���%%%%%%%%%%
            %%%%ǰ�ڼ�������������������ߣ������ʣ���ʱ�����ֲ���������Ϣ%%%%
            %framepointnum=floor(length(speech)/line);
            [~,line]=size(FBEs);
            crosszeroline=zeros(1,line);  %�����źŵĹ��������� 
            framepointnum=floor(fs_ori*Ts/1000);
            speech_FBEs_energy_frameline=zeros(1,line);  %��Ƶ��������
            speechenergyrank=zeros(framepointnum,line);
            %speech_1=speech
            for j=1:line
                %speechenergyrank(j)= sum(speech(framepointnum*(j-1)+1:framepointnum(j)).*speech(framepointnum*(j-1)+1:framepointnum(j)));
                for k=1:framepointnum-1
                   speechenergyrank(:,j)= speech_1(framepointnum*(j-1)+k)*speech_1(framepointnum*(j-1)+k);%������������           
                    if speech_1(framepointnum*(j-1)+k)*speech_1(framepointnum*(j-1)+k+1)<0                
                        crosszeroline(j)=crosszeroline(j)+1;
                    end
                end
                speech_FBEs_energy_frameline(j)=sum(speechenergyrank(:,j));%�����������ԣ���ʱ����ѹ
                %speech_intense_frameline(j)= speech_energy_frameline/2*(j);%����ǿ�ȣ���ǿ
            end
            %%%%ǰ�ڼ�������������������ߣ������ʣ���ʱ�����ֲ���������Ϣ���%%%%
            %%%%�ָ���Ƶ�е�ÿһ�η����ڵ�%%%%
            markspeechtimeline=repmat(markFBEtimeline',floor(fs_ori*Ts/1000),1);
            markspeechtimeline=reshape(markspeechtimeline,[],1);%���һ�����ݱ�עspeech��
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
            voicepressure_rms=rms(speech_1(markspeechtimeline==1));  %��������ѹ��ֵ
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
            speech_sec_inter_20=cat(1,voiceFBEstart_points(1),speech_sec_inter_20);%���ȶ���
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
            %%%%%%%%ʶ��ÿһ������ĳ��ȣ��ٵ����ÿգ������ɾ����%%%%%%%            
            %speech_sec_length_list
            chara_count=numel(speech_sec_length_list);%������ȵĲ��졣
            [value,loc1_9doubt]=max(diff(speech_sec_length_list)); %�������з�����ʱ���ȵı仯���ҵ����ֵ���Ԫ��֮��ļ�϶
            %loc1_9doubt�Ǽ������λ���������������λ��+
            %numel(chara_count)>9    
                speech_content                
                loc_jiu=strfind(speech_content,'��');
                if ~isempty(loc_jiu)
                    speech_content_1_9 = speech_content(1:loc_jiu(1));
                    loc_dot=strfind(speech_content_1_9,',');
                    recognized_the_nine=1;
                else
                    loc_jiu=strfind(speech_content,'��');                    
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
            %%%%%%%%������audio���������ʼ%%%%%%%%
            %markFBEtimeline_1_9size=size(markFBEtimeline_1_9)
            %markFBEtimeline_aongmsize=size(markFBEtimeline_aongm)
            %markFBEtimeline_coughsize=size(markFBEtimeline_cough)
            [audio_1_9,audio_t_1_9,audio_dat_1_9]=data_generaterH6(speech_1_1_9,audio,fs_ori,Ts,FBEs_1_9,markFBEtimeline_1_9,MAG_1_9,quitHz);  
            [audio_aongm,audio_t_aongm,audio_dat_aongm]=data_generaterH6(speech_1_2_aaaongmmm,audio,fs_ori,Ts,FBEs_aongm,markFBEtimeline_aongm,MAG_aongm,quitHz);  
            [audio_cough,audio_t_cough,audio_dat_cough]=data_generaterH6(speech_1_3_cough,audio,fs_ori,Ts,FBEs_cough,markFBEtimeline_cough,MAG_cough,quitHz);  
            %audio_1_9
            %audio_t_1_9
            %audio_dat_1_9
            
            
            %%%%%%%%audio�����������%%%%%%%%
            
            
            audio.voicepressure_rms=voicepressure_rms;%������ѹ��ֵ������
            audio.noisepressure_rms=noisepressure_rms;%������ѹ��ֵ������
            audio.voicepressurelevel_rms=voicepressurelevel_rms;%������ѹ����ֵ������
            audio.noisepressurelevel_rms=noisepressurelevel_rms;%������ѹ����ֵ������
            
            audio.count1_9loss=count1_9loss;%1-9ʶ����©
            audio.count1_9_num=count1_9_num;%1-9����ʶ���λ
            audio.count_aaaongmmm_num=count_aaaongmmm_num;%aaaongmmmʶ����©
            audio.count_cough_num=count_cough_num;%����������
            audio.charanumall=numel(voiceFBEstart_points);
            audio.speech_sec_mean_energy_20=speech_sec_mean_energy_20;%���־�������
            audio.speech_sec_std_20=speech_sec_std_20;%������������
            audio.speech_sec_length_list_20=speech_sec_length_list_20;
            audio.speech_sec_inter_20=speech_sec_inter_20;
                
            %p=ones(1,10)*(1/10);
            %speech_FBEs_energy_frameline=filter(p,1,speech_FBEs_energy_frameline);%��һ����ƽ����ɣ�
               % voicevolume_average=max(abs(speech_FBEs_energy_frameline(markFBEtimeline==1)));%/(noisetimeratio*length(speech)/fs); %ƽ����������

            
           %%%%%%%%%%%2024������ÿһ���������������Ĳ������%%%%%%%%%%

            %%%%%%%%audio������������ʼ%%%%%%%%
            [audio,audio_t,audio_dat]=data_generaterH6(speech_1,audio,fs_ori,Ts,FBEs,markFBEtimeline,MAG,quitHz);  
            %%%%%%%%audio��������������%%%%%%%%
            audio;
            audio_t;
            pause(0.05)
            %%%%%%%%���������벿������%%%%%%%%
            audio_str.audio_str_3={audio.recorddate,audio.datpath,audio.wavname};
            audio_t.audio_str_3={'��¼����','·��','�ļ���'};
            before_sep_para_89={audio.voicepressure_rms;audio.noisepressure_rms;...
                audio.voicepressurelevel_rms;audio.noisepressurelevel_rms;...
                audio.count1_9loss;audio.count1_9_num;...
                audio.count_aaaongmmm_num;audio.count_cough_num;audio.charanumall;...
                audio.speech_sec_mean_energy_20;...    %���־�������
                audio.speech_sec_std_20;...%���־��ܱ�׼������
                audio.speech_sec_length_list_20;...%������������
                audio.speech_sec_inter_20;...%���ּ������
                };
            audio_dat.before_sep_para_89=before_sep_para_89;
            %%%%%%%%%%%
            addpara_10={'������ѹ','������ѹ','������ѹ��','������ѹ��','1-9��ʧλ',...
                '1-9ʶ������','��Ԫ��ʶ������','������ʶ������','�ܹ�ʶ������'};
            numm=1:charanum_upperlimit;    
                %str1=repmat('mel20Ƶ�׷��',numel(numm),1);
                str2=num2str(reshape(numm,[],1));                
                str3=repmat('s',numel(numm),1);
                titles=strcat('ʶ��',num2str(charanum_upperlimit),'���ھ�������',str2,str3);
            title_speech_sec_mean_energy_20=mat2cell(titles,ones(1,numel(numm)))';                       
                titles=strcat('ʶ��',num2str(charanum_upperlimit),'���ھ��ܱ�׼������',str2,str3);
            title_speech_sec_std_20=mat2cell(titles,ones(1,numel(numm)))';
                titles=strcat('ʶ��',num2str(charanum_upperlimit),'������������',str2,str3);
            title_speech_sec_length_list_20=mat2cell(titles,ones(1,numel(numm)))';                       
                titles=strcat('ʶ��',num2str(charanum_upperlimit),'���ڼ������',str2,str3);
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
                    fprintf(TxtId,'%s\t','һ���','�����','�����','ʶ��������','�Ƿ�ʶ�����Ǹ�-��','�ļ�������','��Ƶ����','����');
                    audio_t.audio_str_3={'��¼����','·��','�ļ���'};  %9
                    fprintf(TxtId,'%s\t',audio_t.audio_str_3{1:end});%9
                    fprintf(TxtId,'%s\t',audio_t.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t.title_voice4_100{1:end});  
                    %%%%%%%%%%%%�¼Ӳ���%%%%%%%%%%
                    
                    fprintf(TxtId,'%s\t',audio_t.title_before_sep_para_89{1:end});
                    fprintf(TxtId,'%s\t','1-9�Զ��ָ���Ƭ�����ݷ���');
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_1_9.title_voice4_100{1:end});
                    fprintf(TxtId,'%s\t','aaa-ong-mmm�Զ��ָ���Ƭ�����ݷ���');
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_aongm.title_voice4_100{1:end});
                    fprintf(TxtId,'%s\t','cough���Զ��ָ���Ƭ�����ݷ���');
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice1_236{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice2_116{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice3_150{1:end});
                    fprintf(TxtId,'%s\t',audio_t_cough.title_voice4_100{1:end});
                end
                mainparameters.repeattime=repeattime+1;
                setappdata(0,'mainparameters',mainparameters);
                % Write the output data��ʼд������   
                fprintf(TxtId,'\n');
                fprintf(TxtId,'%5.0f\t',i_1_stage,i_2_stage,i_3_stage,audio.recognizedcharactersnum,audio.recognized_the_nine);%4 num
                fprintf(TxtId,'%s\t',path_1_stage(i_1_stage).name,wav_file,speech_content);%3 str
                fprintf(TxtId,'%s\t',audio_str.audio_str_3{1:end});%audio_dat.recorddate,audio_dat.datpath,audio_dat.wavname);%audio_str.audio_str_3
                       %11+2=13;��Ӧ��һ������ 
                fprintf(TxtId,'%10.4f\t',audio_dat.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat.voice4_100{1:end});%  
                %%%%%%%%%%%%%%%%%%%%%                
                fprintf(TxtId,'%10.4f\t',audio_dat.before_sep_para_89{1:end});%  
                fprintf(TxtId,'%s\t','1-9�Զ��ָ���Ƭ�����ݷ���');
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_1_9.voice4_100{1:end});%
                fprintf(TxtId,'%s\t','aaa-ong-mmm�Զ��ָ���Ƭ�����ݷ���');
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_aongm.voice4_100{1:end});%
                fprintf(TxtId,'%s\t','cough���Զ��ָ���Ƭ�����ݷ���');
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice1_236{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice2_116{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice3_150{1:end});%
                fprintf(TxtId,'%10.4f\t',audio_dat_cough.voice4_100{1:end});%
            else 
                h = msgbox('Record audio characters in files failed!');
            end
            %��ʾ�������text_show
            pause(0.05)
            thecondition=sprintf('��ǰ��������չ��%d,%d,%d��',currentsage)
            set(handles.text_show,'string',thecondition)
        end
    end
end
fclose(TxtId);
