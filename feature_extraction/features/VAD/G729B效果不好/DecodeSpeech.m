function DecodeSpeech

FID = fopen('CodedData.cod', 'r');

CActive = 1;
CSilence = 0;
CSID = 2;

Fs = 8000;
NsFrame = 80;        % Frame size (10ms)

% Set up the quantizer tables
[Yq, Xq, Code, ICode] = QuantMuLawTables;

k = 0;
while (true)
  [CCode, N] = fread(FID, 1, 'uint8=>double');
  if (N < 1)        % Jump out at the end-of-file
    break
  end
  ist = k*NsFrame+1;
  ifn = ist+NsFrame-1;
  if (CCode == CActive)
    [xc, N] = fread(FID, NsFrame, 'uint8=>double');
    if (N < NsFrame)
      error('DecodeSpeechNoise: Incomplete Frame');
    end
    % Reconstruct a speech sample from the received code
    Qindex = ICode(xc+1);
    x(ist:ifn) = Yq(Qindex+1);
  else
    % Insert silence
    x(ist:ifn) = zeros(NsFrame, 1);
  end
  k = k + 1;
end

fclose(FID);

% Write the data to an audio file
wavwrite(x, Fs, 16, 'CodedSpeech.wav');

return
