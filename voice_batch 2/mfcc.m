function [ CC,FBE,frames,MAG,nfft,quitHz ] = mfcc( speech, fs, Tw, Ts, alpha, window, R, M, N, L )
% MFCC Mel frequency cepstral coefficient feature extraction.
% fs.采样频率, Tw.25, Ts.10, alpha.0.97, window.@hamming, R[300,3700], M.20, N.12+1,
% L.22
%   MFCC(S,FS,TW,TS,ALPHA,WINDOW,R,M,N,L) returns mel frequency 
%   cepstral coefficients (MFCCs) computed from speech signal given 
%   in vector S and sampled at FS (Hz). The speech signal is first 
%   preemphasised using a first order FIR filter with preemphasis 
%   coefficient ALPHA. The preemphasised speech signal is subjected 
%   to the short-time Fourier transform analysis with frame durations 
%   of TW (ms), frame shifts of TS (ms) and analysis window function 
%   given as a function handle in WINDOW. This is followed by magnitude 
%   spectrum computation followed by filterbank design with M triangular 
%   filters uniformly spaced on the mel scale between lower and upper 
%   frequency limits given in R (Hz). The filterbank is applied to 
%   the magnitude spectrum values to produce filterbank energies (FBEs) 
%   (M per frame). Log-compressed FBEs are then decorrelated using the 
%   discrete cosine transform to produce cepstral coefficients. Final
%   step applies sinusoidal lifter to produce liftered MFCCs that 
%   closely match those produced by HTK [1].
%
%   [CC,FBE,FRAMES]=MFCC(...) also returns FBEs and windowed frames,
%   with feature vectors and frames as columns.
%
%   This framework is based on Dan Ellis' rastamat routines [2]. The 
%   emphasis is placed on closely matching MFCCs produced by HTK [1]
%   (refer to p.337 of [1] for HTK's defaults) with simplicity and 
%   compactness as main considerations, but at a cost of reduced 
%   flexibility. This routine is meant to be easy to extend, and as 
%   a starting point for work with cepstral coefficients in MATLAB.
%   The triangular filterbank equations are given in [3].
%
%   Inputs
%           S is the input speech signal (as vector)
%
%           FS is the sampling frequency (Hz) 
%
%           TW is the analysis frame duration (ms) 
% 
%           TS is the analysis frame shift (ms)
%
%           ALPHA is the preemphasis coefficient
%
%           WINDOW is a analysis window function handle
% 
%           R is the frequency range (Hz) for filterbank analysis
%
%           M is the number of filterbank channels
%
%           N is the number of cepstral coefficients 
%             (including the 0th coefficient)
%
%           L is the liftering parameter
%
%   Outputs
%           CC is a matrix of mel frequency cepstral coefficients
%              (MFCCs) with feature vectors as columns
%
%           FBE is a matrix of filterbank energies
%               with feature vectors as columns
%
%           FRAMES is a matrix of windowed frames
%                  (one frame per column)
%
%   Example
%           Tw = 25;           % analysis frame duration (ms)
%           Ts = 10;           % analysis frame shift (ms)
%           alpha = 0.97;      % preemphasis coefficient
%           R = [ 300 3700 ];  % frequency range to consider
%           M = 20;            % number of filterbank channels 
%           C = 13;            % number of cepstral coefficients
%           L = 22;            % cepstral sine lifter parameter
%       
%           % hamming window (see Eq. (5.2) on p.73 of [1])
%           hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
%       
%           % Read speech samples, sampling rate and precision from file
%           [ speech, fs, nbits ] = wavread( 'sp10.wav' );
%       
%           % Feature extraction (feature vectors as columns)
%           [ MFCCs, FBEs, frames ] = ...
%                           mfcc( speech, fs, Tw, Ts, alpha, hamming, R, M, C, L );
%       
%           % Plot cepstrum over time
%           figure('Position', [30 100 800 200], 'PaperPositionMode', 'auto', ... 
%                  'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 
%       
%           imagesc( [1:size(MFCCs,2)], [0:C-1], MFCCs ); 
%           axis( 'xy' );
%           xlabel( 'Frame index' ); 
%           ylabel( 'Cepstrum index' );
%           title( 'Mel frequency cepstrum' );
%
%   References
%
%           [1] Young, S., Evermann, G., Gales, M., Hain, T., Kershaw, D., 
%               Liu, X., Moore, G., Odell, J., Ollason, D., Povey, D., 
%               Valtchev, V., Woodland, P., 2006. The HTK Book (for HTK 
%               Version 3.4.1). Engineering Department, Cambridge University.
%               (see also: http://htk.eng.cam.ac.uk)
%
%           [2] Ellis, D., 2005. Reproducing the feature outputs of 
%               common programs using Matlab and melfcc.m. url: 
%               http://labrosa.ee.columbia.edu/matlab/rastamat/mfccs.html
%
%           [3] Huang, X., Acero, A., Hon, H., 2001. Spoken Language 
%               Processing: A guide to theory, algorithm, and system 
%               development. Prentice Hall, Upper Saddle River, NJ, 
%               USA (pp. 314-315).
%
%   See also EXAMPLE, COMPARE, FRAMES2VEC, TRIFBANK.

%   Author: Kamil Wojcicki, September 2011


    %% PRELIMINARIES 

    % Ensure correct number of inputs
    if( nargin~= 10 ), help mfcc; return; end; 

    % Explode samples to the range of 16 bit shorts 扩大原始音频信号2的16次方，平移一个长整形
    if( max(abs(speech))<=1 ), speech = speech * 2^15; end;

    Nw = round( 1E-3*Tw*fs );    % frame duration (samples)1e-3*25*12500==round(312.5)=313帧每样本
    Ns = round( 1E-3*Ts*fs );    % frame shift (samples)==313 帧与帧的间隔

    nfft = 2^nextpow2( Nw ) ;    % length of FFT analysis 取二进制指数值大的，应是256
    K = nfft/2+1;                % length of the unique part of the FFT 中间值，应是129


    %% HANDY INLINE FUNCTION HANDLES

    % Forward and backward mel frequency warping (see Eq. (5.13) on p.76 of [1]) 
    % Note that base 10 is used in [1], while base e is used here and in HTK code
    hz2mel = @( hz )( 1127*log(1+hz/700) );     % Hertz to mel warping function
    mel2hz = @( mel )( 700*exp(mel/1127)-700 ); % mel to Hertz warping function

    % Type III DCT matrix routine (see Eq. (5.14) on p.77 of [1])
    dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );%0.32*cos(repmat([0:13-1].',1,20)N12行复制M20列.*M20列复制N行，构成矩阵

    % Cepstral lifter routine (see Eq. (5.12) on p.75 of [1])
    ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );


    %% FEATURE EXTRACTION 

    % Preemphasis filtering (see Eq. (5.1) on p.73 of [1])预加重滤波
    speech = filter( [1 -alpha], 1, speech ); % fvtool( [1 -alpha], 1 );%[1  -0.97]，对原始音频信号进行预先滤波
    size(speech)
    % Framing and windowing (frames as columns)%分割相互重叠的帧信息
    frames = vec2frames( speech, Nw, Ns, 'cols', window, false );% Nw.313, Ns313, 'cols'每帧按列排, window.@hamming, false

    % Magnitude spectrum computation (as column vectors)%每一帧的频率列向量
    MAG = abs( fft(frames,nfft,1) ); 

    % Triangular filterbank with uniformly spaced filters on mel scale
    [H1,f,c] = trifbank( M, K, R, fs, hz2mel, mel2hz ); % size of H is M x K 三角滤波，依据梅尔尺度统一滤波频段空间 
   %H1
   diff(f);%频率被均分，2048项，nfft的一半，,23.4375Hz间隔，合起来就是48000Hz
   c;
   quitHz=c(1:end-1)*c(end)/c(end-1);
    %采样率48000时，频率上限是24000Hz，上下限为300-3700Hz，采样率4000，频率上限为2000Hz。上下限等比为
    %%%自加频率截止点%%%    
    numel(quitHz);
    M+1;
    quitHz(end);
    floor(fs*0.5);
    isequal(quitHz(end),floor(fs*0.5)  );
    numel(quitHz);
    M+1;
    isequal(numel(quitHz),(M+1));
    
    %
    if  ~isequal(numel(quitHz),(M))  
        p=zeros(1,size(H1,1)+1);quitHz=zeros(1,size(H1,1));pall=0;
        for i=1:size(H1,1)
            p(i)=length(find(H1(i,:)));
            pall=pall+p(i);
            quitHz(i)=pall;               
        end           
        %p=p*(8000/pall);
        quitHz=quitHz*( floor(fs*0.5)/quitHz(end));
    end
    quitHz;
    
    %}    
     %%%自加频率截止点结束%%%      
           
    % Filterbank application to unique part of the magnitude spectrum
    FBE = H1 * MAG(1:K,:); % FBE( FBE<1.0 ) = 1.0; % apply mel floor,筛选截取前K段的频率分段，选择H过滤需要的频段,转换成人的听觉所匹配的频率与幅度。

    % DCT matrix computation
    DCT = dctm( N, M );

    % Conversion of logFBEs to cepstral coefficients through DCT
    CC =  DCT * log( FBE );

    % Cepstral lifter computation
    lifter = ceplifter( N, L );  %N=12+1=13or 19+1=20， L=22

    % Cepstral liftering gives liftered cepstral coefficients
    CC = diag( lifter ) * CC; % ~ HTK's MFCCs


% EOF
