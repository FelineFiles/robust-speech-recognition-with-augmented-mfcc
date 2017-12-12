function feature=extract_features(input,fsamp)
% -------------------------------------------------------------------------
winlen      = 25;             % window length in ms
winshft     = 2.5;            % dense window shift in ms
nfft        = 512;            % fft size 
cepnum      = 13;             % number of cepstral coefficients
liftercoe   = 22;             % liftering coefficient (e.g. 22)
numchan     = 26;             % number of channels of the MEL filter bank (e.g. 26)
preemcoeff  = 0.97;           % coefficient for pre-emphasis
doDisp=0;                     % Display option- '1' outputs a graphic display, '0' does not
% -------------------------------------------------------------------------

Fs = fsamp;

winlen=round(winlen*10^(-3)*fsamp);
winshft=winshft*10^(-3)*fsamp;
FrameNo=ceil((length(input)-winlen)/winshft);

% initialize MEL filter bank
fbank = initfiltb(winlen, numchan, fsamp, nfft);

% initialize lifter coefficients
lifter = (1 + (liftercoe/2)*sin((pi/liftercoe)*(0:cepnum)) );

% pre-emphasis
am = [1 0];  % denominator polynomial
bm = [1 -preemcoeff];            % numerator   polynomial
preem = filter(bm, am, input);

% change signal (a vector) into frame (a matrix), where each collum is a frame
frmwin = sig2fm(input, winlen, winshft, FrameNo);    
[winlen, framenum]=size(frmwin); 

% Hamming window each frame
frmwin = frmwin .* (hamming(winlen) * ones(1, framenum));

% FFT
ffto=abs(fft(frmwin,nfft)).^2;

% Pitch Detection for VTLN
vad = voiceActivityDetect(ffto, Fs, nfft);
[avgpitch, pitch_arr] = harmprodspec(ffto(1 : (nfft/2), :), vad);
pitch = avgpitch/nfft * Fs;

% Peforms harmonic demodulation
width = 15;
pad = (width-1)/2;
idx = bsxfun(@plus, (-pad:pad)', (0:nfft-1));
idx= idx'+1;
idx_temp = idx;
idx_temp(idx<1 | idx > nfft) = 1;
ffto_seg = ffto(idx_temp, :);
ffto_seg(idx<1 | idx > nfft) = 0;
ffto_seg = reshape(ffto_seg, [nfft, width, framenum]);
h = cos((-pad:pad)*pi/(width+1));
ffto = bsxfun(@times, h, ffto_seg);
ffto = max(ffto, [], 2);

% VTLN
femamax = 195;
femamin = 175;
vtlnmax1 = 1.22;
alpha = 1.30;
vtlnmax2 = alpha;
vtlnmin = 1.14;
tbeta = 1200;
fedges = tbeta / Fs * nfft;
if pitch >= femamin
    if pitch > femamax
        pitch = femamax;
    end
    scalefactor1 = vtlnmin + (vtlnmax1-vtlnmin)/(femamax-femamin) * (pitch - femamin);
    scalefactor2 = vtlnmin + (vtlnmax2-vtlnmin)/(femamax-femamin) * (pitch - femamin);

    ffto = freqshift(ffto(1 : (nfft/2), :), scalefactor1,scalefactor2, fedges);
end


% Then perform noise flooring
% noise_floor_threshold = 0.4;
% energy = mean(ffto,1);
% for fr=1:framenum
%    ffto(ffto(:,fr) < noise_floor_threshold*energy(fr)) = noise_floor_threshold*energy(fr);
% end

% Get VAD
vad = voiceActivityDetect(ffto, Fs, nfft);

% SNR Noise Surpression
snr = calcSNR(ffto(1:nfft/2, :), vad);
ffto = noiseSurpress(ffto, snr, vad);
ffto = squeeze(ffto);

% Increase silence noise variance when training on clean data
if snr > 14.4
   ffto(:, vad==0) = ffto(:, vad==0) * 1.2;
else
   ffto(:, vad==0) = ffto(:, vad==0) * 0.1;
end
% subplot(2,1,2);
% imagesc(ffto(1:100,:));

% MEL filtering 
% fbank = zeros(numchan, nfft/2);
% [fbank_b, fbank_a] = gamma_bank(winlen, numchan, fsamp, nfft);
% for i=1:size(fbank_b,1)
%     fbank(i,:) = abs(freqz(fbank_b(i,:), fbank_a(i,:), nfft/2));
%     fbank(i,:) = fbank(i,:)/sqrt(sum(abs(fbank(i,:)).^2));
% end

fb=fbank.^2*ffto(1 : (nfft/2), :);

% take logarithm of MEL filter output
fbfloor=mean(mean(fb))*0.00001;  
logfb=log(max(fb, fbfloor*rand(size(fb))));

% 
if doDisp
figure(1);
subplot(2,1,1);
imagesc([0, 2.5*size(ffto,2)], [0, Fs/2], 10*log10(ffto(1:nfft/2,:)));
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Power Spectrum of Speech');
subplot(2,1,2);
plot(linspace(0,2.5*size(ffto,2),numel(vad)),vad);
xlabel('Time (ms)');
ylabel('Activity Level ()');
title('Speech Activity Level');
axis tight;
end


% take DCT
mfcco=dct(logfb);
mfcco=mfcco(1:cepnum+1,:);

% lifer MFCCs
mfcco=mfcco.*(lifter'*ones(1,framenum));

% determine frame entropies
ENT_TEMP=find_entropy(input,fsamp,numchan,nfft);

keep=ones(1,framenum);
ENT=[];
ctr=1;
rate=8;
ho=1;
for fr=7:6:size(fb,2)-6
    curr_VAD=vad(:,fr-3:fr+2);
    if(sum(curr_VAD)==0)
     	start=fr-3-ho-1;
        finish=fr+2;
        keep(start:rate:finish)=0;
        ho=mod(finish-start,rate);
    else
        ho=0;
        ENT=[ENT ENT_TEMP(ctr)];
    end;
    ctr=ctr+1;
end;

if(length(ENT)>0)
    Mx=max(ENT);
    Md=median(ENT);
    Mn=min(ENT);
else
    Mx=1;
    Md=1;
    Mn=1;
end;
 
% perform VFR analysis, according to Section 3 of (You et al., 2004)
ctr=1;
ent_ctr=1;
w1=0.7;
w2=0.8;
w3=0.5;
T1=w1*Mx+(1-w1)*Md;
T2=(1-w2)*Mx+w2*Md;
T3=(1-w3)*Md+w3*Mn;
ho=1;
for fr=7:6:size(fb,2)-6
    curr_VAD=vad(:,fr-3:fr+2);
    if(sum(curr_VAD)>0)
        curr_frames=mfcco(:,fr-3:fr+2);
        curr_ent=ENT(ent_ctr);
        ent_ctr=ent_ctr+1;
        if(curr_ent>T1)
            rate=2;
        elseif(curr_ent>T2)
            rate=3;
        elseif(curr_ent>T3)
            rate=4;
        else
            rate=5;
        end;
        start=fr-3-ho-1;
        finish=fr+2;
        keep(start:rate:finish)=0;
        ho=mod(finish-start,rate);
        plot_ent(ent_ctr-1)=curr_ent;
    else
        ho=0;
        rate=8;
        curr_ent=0;
    end;
    plot_rate(ctr)=rate;
    ctr=ctr+1;
end;
S1=size(mfcco,2);
mfcco=mfcco(:,find(keep==0));
vad = vad(:,find(keep==0));

if(doDisp)
    figure;
    subplot(3,1,1);plot(input);axis tight;
    title('Speech Signal');
    subplot(3,1,2);plot(plot_ent);axis tight;
    title('Frame Entropy Values');
    subplot(3,1,3);imagesc(ones(10,1)*keep);colormap(gray);
    title('Selected Frames');
end;


%% Noise robust stuff

% Peak Isolation
C0 = mfcco(1,:);
%mfcco(1,:) = 0;
mfcco = idct(mfcco);
mfcco(mfcco<0) = 0;

if doDisp
figure(2);
subplot(2,1,1);
imagesc(mfcco(2:end,:));
xlabel('Frame Number ()');
ylabel('Cepstral Number ()');
title('Before Ratio Locking');
colorbar;
end

% Ratio locking
% if ~strcmp(type, 'train')
% p_max = 45;
% for n = 1:size(mfcco, 2)
%     if max(mfcco(2:end,n)) > 0 && vad(n) == 1
%         mfcco(2:end,n) = mfcco(2:end,n)/ max(mfcco(2:end,n)) * p_max;
%     end
% end
% end 

if doDisp

subplot(2,1,2);
imagesc(mfcco(2:end,:));
xlabel('Frame Number ()');
ylabel('Cepstral Number ()');
title('After Ratio Locking');
colorbar;
end

if doDisp
    a=0;
end

mfcco = dct(mfcco);


% perform cepstral mean subtraction (CMS)
mfcco = mfcco - mean(mfcco, 2) * ones(1, size(mfcco, 2));

% determine deltas and double-deltas, and concatenate into a matrix of
% feature vectors
dt1=deltacc(mfcco);
dt2=deltacc(dt1);
mfcco=[mfcco;dt1;dt2];
mfcco=mfcco(1:end,:);
feature = mfcco';


% ---------------------------------------------------------------
% ---------------------------------------------------------------
function ENT=find_entropy(speech,fsamp,numchan,nfft)
% Function approximates the feature space entropy of an input speech signal
% on a frame-by-frame basis, according to Eq. 5 of (You et al., 2004).
% 'speech' is the raw speech signal
% 'fsamp' is the sampling rate of 'speech'
% 'numchan' is the number of Mel-channels used during short-time analysis
% 'nfft' is the size of the DFT used during short-time analysis
% The function returns 'ENT', a vector of frame entropies

winlen=round(25*10^(-3)*fsamp);
winshft=2.5*10^(-3)*fsamp;
FrameNo = ceil((length(speech) - winlen) / winshft);

% initialize MEL filter bank
fbank=initfiltb(winlen,numchan,fsamp,nfft);

% change signal (a vector) into frame (a matrix), where each collum is a frame
frmwin=sig2fm(speech,winlen,winshft,FrameNo);    
[winlen,framenum]=size(frmwin); 

% Hamming window each frame
frmwin=frmwin.*(hamming(winlen)*ones(1,framenum));

% FFT
ffto=abs(fft(frmwin,nfft));

% MEL filtering 
fb=fbank*ffto(1:(nfft/2),:);

% Approximate entropy of each MFCC frame according to 
% Eq. 5 of (You et al., 2004)
ctr=1;
for fr=7:6:framenum-6
    TR=0;
    curr_frames=fb(:,fr-6:fr+5);
    MU=mean(curr_frames')';
    MU=MU*ones(1,size(curr_frames,2));
    SIGMA=(curr_frames-MU)*(curr_frames-MU)';
    for k=1:numchan
        TR=TR+SIGMA(k,k);
    end;
    TR=log(TR)/log(exp(1));
    ENT(ctr)=numchan*log((2*pi)^0.5)/log(exp(1))+TR;
    ctr=ctr+1;
end;

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function mels=mel(freq)
% change frequency from Hz to mel
mels=1127*log(1+(freq/700));

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function wins=sig2fm(input,winlen,winshft,frameno)
% put vector into matrix, each column is a frame. 
% The rest of signal that is less than one frame is discarded
% winlen, winshft are in number of sample, notice winshft is not limited to
% integer
input=input(:);     
wins=zeros(winlen, frameno);
for i=1:frameno
    b=round((i-1)*winshft);
    c=min(winlen,length(input)-b);
    wins(1:c,i)=input(b+1:min(length(input),b+winlen));
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function fbank=initfiltb(framelen,numchan,fsamp,nfft)
% triangle shape melfilter initialization

fftfreqs=((0:(nfft/2-1))/nfft)*fsamp;  % frequency of each fft point (1-fsamp/2)
melfft=mel(fftfreqs);   % mel of each fft point
mel0=0;                  
mel1=mel(fsamp/2);       % highest mel 
melmid=((1:numchan)/(numchan+1))*(mel1-mel0)+mel0; % middle mel of each filter
fbank=zeros(numchan,nfft/2);

% non overlaping triangle window is used to form the mel filter
for k=2:(nfft/2)  % for each fft point, to all the filters,do this:
    chan=max([0 find(melfft(k)>melmid)]); % the highest index of melfft that is larger than the middle mel of all channels
    if(chan==0)  % only the first filter cover here
        fbank(1,k)=(melfft(k)-mel0)/(melmid(1)-mel0);
    elseif(chan==numchan)  % only the last filter covered here
        fbank(numchan,k)=(mel1-melfft(k))/(mel1-melmid(chan));
    else                   % for any other part, there will be two filter cover that frequency, in the complementary manner
        fbank(chan,k)=(melmid(chan+1)-melfft(k))/(melmid(chan+1)-melmid(chan));
        fbank(chan+1,k)=1-fbank(chan,k);  % complementary
 	end
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function dt=deltacc(input)
% calculates derivatives of a matrix, whose columns are feature vectors

tmp=0;
for cnt=1:2
    tmp=tmp+cnt*cnt;
end
nrm=1/(2*tmp);
dt=zeros(size(input));
rows=size(input,1);
cols=size(input,2);
for col=1:cols
    for cnt=1:2
        inx1=col-cnt; 
        inx2=col+cnt;
        if(inx1<1)
            inx1 = 1;     
        end;
        if(inx2>cols)  
            inx2 = cols;  
        end;
        dt(:,col)=dt(:,col)+(input(:,inx2)-input(:,inx1))*cnt;
    end
end
dt=dt*nrm;
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function vad = voiceActivityDetect(ffto, Fs, nfft)
    low_freq = 62.5;
    high_freq = 1015.6;
    low_idx = round(low_freq*nfft/Fs)+1;
    high_idx = round(high_freq*nfft/Fs)+1;
    vad = sum(ffto(low_idx:high_idx,:), 1);
    sil_sc = 2.0;
    [N,edges] = histcounts(vad);
%     figure;
%     histogram(vad);
%     ylabel('Counts');
%     xlabel('Energy');
%     title('VAD Histogram');

    [mx, lc] = max(N);
    sil_min = edges(lc+1);
    sil_u = mean(vad(vad<sil_min));
    sil_std = std(vad(vad<sil_min));
    vad(vad<sil_u + sil_std * sil_sc) = 0;
  
    vad(vad>sil_u + sil_std * sil_sc) = 1;
    vad = medfilt1(vad,30); % smooth irregularities
    vad(vad>0) = 1;
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function [pitch, locs] = harmprodspec(spec, vad)
  nFrame = size(spec,2);
  bpw2 = 68;
  mx =  zeros(1, nFrame);
  locs = zeros(1, nFrame);
  peri = (4:1:18);
  maxi = zeros(size(peri));
  
  for m = 1:nFrame
     
    if vad(m)
      specs = spec(1:bpw2,m);

      for n = 1:length(peri)
        spec_sh = circshift(specs, peri(n), 1);
        spec_sh(1:peri(n)) = 0;
        maxi(n) =  sum(specs .* spec_sh);
      end

      [mxx, loc]  = findpeaks(maxi);
      if isempty(loc)
          locs(m) = 0;
          mx(m) = 0;
      else
        locs(m) = peri(loc(1));
        mx(m) = mxx(1);
      end

    end
  end
  vsig = locs(locs>0);
  pitch = mean(vsig);
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function snr = calcSNR(ffto, vad)
  nFrame = size(ffto,2);
  %histogram binning
  %SNR Calc
  n_pow = 0;
  for n = 1:nFrame
     n_pow = n_pow + sum(ffto(:, n)) * (1 - vad(n));
  end
  n_pow = n_pow / nFrame;
  s_pow = 0;
  for n = 1:nFrame
    s_pow = s_pow + (sum(ffto(:, n)) - n_pow) * vad(n);
  end
  s_pow = s_pow / nFrame;
  snr = 10*log10(s_pow/n_pow);
% ---------------------------------------------------------------
  
% ---------------------------------------------------------------
% ---------------------------------------------------------------
function ffto = noiseSurpress(ffto, snr, vad)
    if snr > 14.4
        return;
    else
      spLevel = 0.01;
      nFrames = length(vad);
      sp_fil = vad;
      x = [-2, -1, 0, 1, 2];
      y = fliplr(sigmf(x,[1 0]) * (1-spLevel) + spLevel);
      vc = 1;
      % Forward Populate
      for n = 1:nFrames
          if  vad(n) == 1
              sp_fil(n) = 1;
              vc = 1;
          elseif vad(n) == 0 && vc <= length(y);
              sp_fil(n) = y(vc);
              vc = vc + 1;
          else 
              sp_fil(n) = spLevel;
          end
      end
      for n = nFrames:-1:1
          if  vad(n) == 1
              sp_fil(n) = 1;
              vc = 1;
          elseif vad(n) == 0 && vc <= length(y);
              sp_fil(n) = max(y(vc),sp_fil(n));
              vc = vc + 1;
          else 
              sp_fil(n) = max(spLevel, sp_fil(n));
          end
      end
      sp_fil(1:5) = spLevel;
%       imagesc(ffto);
      for n = 1:nFrames
          ffto(:,n) = ffto(:,n) * sp_fil(n);
      end
%       figure;
%       imagesc(ffto);
    end
% --------------------------------------------------------------- 

function [b,a] = gamma_bank(framelen,numchan,fsamp,nfft)
    fftfreqs=((0:(nfft/2-1))/nfft)*fsamp;  % frequency of each fft point (1-fsamp/2)
    melfft=frq2mel(fftfreqs);   % mel of each fft point
    mel0=0;                  
    mel1=frq2mel(fsamp/2);       % highest mel 
    mel_mid=((1:numchan)/(numchan+1))*(mel1-mel0)+mel0; % middle mel of each filter
    mel_range = [mel0, mel_mid, mel1]; 
    frq_range = mel2frq(mel_range);
    frq_bw = circshift(frq_range,-2,2) - frq_range;
    frq_bw = frq_bw(1:end-2);

    [b,a] =  gammabank(0, fsamp, [], frq_range(2:end-1), frq_bw);

function ffto = freqshift(spec, sf1, sf2, f1)
    speclen = size(spec,1);
    nFrame = size(spec,2);
    for m = 1:nFrame
      temp = zeros(speclen, 1);
      % Linear interpolation
      for n = 1:speclen

          if n < f1
              oldn = n * sf1;
          else %if n >= f1 %&& n <= f2
              oldn = f1*sf1 + (n-f1) * sf2;
          end
          n1 = ceil(oldn);
          n2 = max(floor(oldn), 1);
          if n1 >= speclen || n2 >= speclen
              break;
          end
          if n1 == n2
            temp(n) = spec(n1, m);
          else
            % interpolate
            temp(n) = spec(n1, m) + (spec(n2, m) - spec(n1, m))/(n2-n1) * (oldn - n1);
          end
      end
      spec(:,m) = temp;
    end
ffto = spec;