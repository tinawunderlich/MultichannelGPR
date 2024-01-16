function writesegy3D(filename,data,varargin)

% writesegy3D(filename,data,varargin)
% Writes 3D-sgy-file for Kingdom Suite
%
% Revision 1, 4-byte IEEE floating point, big endian
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% filename: Complete filename and path
% data: matrix with traces in columns ([ns,ntraces]=size(data))
% varargin have to be combinations of names and numbers:
% 'dt', 'ns', 'cdpX', 'cdpY', 'Inline', 'Crossline','topo'


    [ns,ntraces]=size(data);

    % evaluate varargin names&numbers
    for i=1:2:length(varargin)
        var=varargin{i};
        val=varargin{i+1};
        eval([var,'=[',num2str(val),'];']);
    end
    
    if ~exist('dt')
        dt=0.1;
        disp('No dt found. Setting dt=0.1');
    end
    if ~exist('ns')
        ns=1024;
        disp('No ns found. Setting =1024');
    end
    


    % h structure
    h.TextualFileHeader=sprintf('%3200s','created with MultichannelGPR: writesegy.m');
    h.Job=0;
    h.Line=0;
    h.Reel=0;
    h.DataTracePerEnsemble=0;
    h.AuxiliaryTracePerEnsemble=0;
    h.dt=round(dt.*1e4);
    h.ns=ns;
    h.SegyFormatRevisionNumber=100; % Revision 1 
    h.DataSampleFormat=5;  % datasample format: 4-byte IEEE floating point
    h.dtOrig=h.dt;
    h.nsOrig=h.ns;
    h.EnsembleFold=1;
    h.TraceSorting=0;
    h.VerticalSumCode=0;
    h.SweepFrequencyStart=0;
    h.SweepFrequencyEnd=0;
    h.SweepLength=0;
    h.SweepType=0;
    h.SweepChannel=0;
    h.SweepTaperlengthStart=0;
    h.SweepTaperLengthEnd=0;
    h.TaperType=0;
    h.CorrelatedDataTraces=0;
    h.BinaryGain=0;
    h.AmplitudeRecoveryMethod=0;
    h.MeasurementSystem=1;
    h.ImpulseSignalPolarity=0;
    h.VibratoryPolarityCode=0;
    % 3261-3500 UNASSIGNED
    h.Unassigned1=zeros(1,120);
    h.FixedLengthTraceFlag=1;
    h.NumberOfExtTextualHeaders=0;
    % 3506-3600 UNASSIGNED
    h.Unassigned2=zeros(1,47);


    % create file handle:
    segyfh = fopen(filename,'w','b'); % big endian

    % write header:
    fseek(segyfh,0,'bof');
    fwrite(segyfh,h.TextualFileHeader(1:3200),'uchar');     % 1-3200
    fwrite(segyfh,h.Job,'int32');                           % 3204
    fwrite(segyfh,h.Line,'int32');                          % 3208
    fwrite(segyfh,h.Reel,'int32');                          % 3212
    fwrite(segyfh,h.DataTracePerEnsemble,'int16');          % 3214
    fwrite(segyfh,h.AuxiliaryTracePerEnsemble,'int16');     % 3216
    fwrite(segyfh,h.dt,'uint16');                           % 3218
    fwrite(segyfh,h.dtOrig,'uint16');                       % 3220
    fwrite(segyfh,h.ns,'uint16');                           % 3222
    fwrite(segyfh,h.nsOrig,'uint16');                       % 3224
    fwrite(segyfh,h.DataSampleFormat,'int16');              % 3226
    fwrite(segyfh,h.EnsembleFold,'int16');                  % 3228
    fwrite(segyfh,h.TraceSorting,'int16');                  % 3230
    fwrite(segyfh,h.VerticalSumCode,'int16');               % 3232
    fwrite(segyfh,h.SweepFrequencyStart,'int16');           % 3234
    fwrite(segyfh,h.SweepFrequencyEnd,'int16');             % 3236
    fwrite(segyfh,h.SweepLength,'int16');                   % 3238
    fwrite(segyfh,h.SweepType,'int16');                     % 3240
    fwrite(segyfh,h.SweepChannel,'int16');                  % 3242
    fwrite(segyfh,h.SweepTaperlengthStart,'int16');         % 3244
    fwrite(segyfh,h.SweepTaperLengthEnd,'int16');           % 3246
    fwrite(segyfh,h.TaperType,'int16');                     % 3248
    fwrite(segyfh,h.CorrelatedDataTraces,'int16');          % 3250
    fwrite(segyfh,h.BinaryGain,'int16');                    % 3252
    fwrite(segyfh,h.AmplitudeRecoveryMethod,'int16');       % 3254
    fwrite(segyfh,h.MeasurementSystem,'int16');             % 3256
    fwrite(segyfh,h.ImpulseSignalPolarity,'int16');         % 3258
    fwrite(segyfh,h.VibratoryPolarityCode,'int16');         % 3260
    % 3261-3500 UNASSIGNED1 => (120 int32 = 240 bytes)
    fwrite(segyfh,h.Unassigned1,'int16');                   % 3500
    fwrite(segyfh,h.SegyFormatRevisionNumber,'uint16');     % 3502
    fwrite(segyfh,h.FixedLengthTraceFlag,'int16');          % 3504
    fwrite(segyfh,h.NumberOfExtTextualHeaders,'uint16');    % 3506
    % 3506-3600 UNASSIGNED2 => 94/2=47 int16
    fwrite(segyfh,h.Unassigned2,'int16');                   


    % create trace header
    trh=createTraceheader(ns,dt);

    % write every trace:
    for i=1:ntraces
        if mod(i,100)==100
            disp(['Writing trace ',num2str(i),' of ',num2str(ntraces)]);
        end

        % set trace header data for current trace
        trh.TraceNumber=i;
        trh.TraceSequenceFile=i;

        trh.cdpX = cdpX(i)*100;
        trh.cdpY = cdpY(i)*100;
        trh.SourceX=cdpX(i)*100;
        trh.SourceY = cdpY(i)*100;
        trh.GroupX=cdpX(i)*100;
        trh.GroupY = cdpY(i)*100;
        trh.ReceiverGroupElevation=topo(i)*100;
        trh.SourceSurfaceElevation=topo(i)*100;
        trh.Inline3D=Inline(i);
        trh.Crossline3D=Crossline(i);

        % Write the current trace
        WriteTrace(segyfh,data(:,i),trh,h);
    end

    fclose(segyfh);
end


function trh=createTraceheader(ns,dt)

    da=clock;
    DayOfYear=datenum(0,da(2),da(3));

    trh.TraceSequenceLine=0;
    trh.TraceSequenceFile=0;
    trh.FieldRecord=0;
    trh.TraceNumber=0;
    trh.EnergySourcePoint=0;
    trh.cdp=0;
    trh.cdpTrace=0;
    trh.TraceIdentificationCode=1;
    trh.NSummedTraces=0;
    trh.NStackedTraces=1;
    trh.DataUse=0;
    trh.offset=0;
    trh.ReceiverGroupElevation=0;
    trh.SourceSurfaceElevation=0;
    trh.SourceDepth=0;
    trh.ReceiverDatumElevation=0;
    trh.SourceDatumElevation=0;
    trh.SourceWaterDepth=0;
    trh.GroupWaterDepth=0;
    trh.ElevationScalar=-100;
    trh.SourceGroupScalar=-100; % divide coordinates by 100 to give meters!
    trh.SourceX=0;
    trh.SourceY=0;
    trh.GroupX=0;
    trh.GroupY=0;
    trh.CoordinateUnits=0;
    trh.WeatheringVelocity=0;
    trh.SubWeatheringVelocity=0;
    trh.SourceUpholeTime=0;
    trh.GroupUpholeTime=0;
    trh.SourceStaticCorrection=0;
    trh.GroupStaticCorrection=0;
    trh.TotalStaticApplied=0;
    trh.LagTimeA=0;
    trh.LagTimeB=0;
    trh.DelayRecordingTime=0;
    trh.MuteTimeStart=0;
    trh.MuteTimeEND=0;
    trh.ns=ns;
    trh.dt=round(dt*1e4);
    trh.GainType=0;
    trh.InstrumentGainConstant=0;
    trh.InstrumentInitialGain=0;
    trh.Correlated=0;
    trh.SweepFrequenceStart=0;
    trh.SweepFrequenceEnd=0;
    trh.SweepLength=0;
    trh.SweepType=0;
    trh.SweepTraceTaperLengthStart=0;
    trh.SweepTraceTaperLengthEnd=0;
    trh.TaperType=0;
    trh.AliasFilterFrequency=0;
    trh.AliasFilterSlope=0;
    trh.NotchFilterFrequency=0;
    trh.NotchFilterSlope=0;
    trh.LowCutFrequency=0;
    trh.HighCutFrequency=0;
    trh.LowCutSlope=0;
    trh.HighCutSlope=0;

    trh.YearDataRecorded=da(1);
    trh.DayOfYear=DayOfYear;
    trh.HourOfDay=da(4);
    trh.MinuteOfHour=da(5);
    trh.SecondOfMinute=round(da(6));
    trh.TimeBaseCode=0;
    trh.TraceWeightningFactor=0; 
    trh.GeophoneGroupNumberRoll1=0;
    trh.GeophoneGroupNumberFirstTraceOrigField=0;
    trh.GeophoneGroupNumberLastTraceOrigField=0;
    trh.GapSize=0;
    trh.OverTravel=0;
    trh.cdpX=0;
    trh.cdpY=0;
    trh.Inline3D=0;
    trh.Crossline3D=0;
    trh.ShotPoint=0;
    trh.ShotPointScalar=0;
    trh.TraceValueMeasurementUnit=0;
    trh.TransductionConstantMantissa=0;
    trh.TransductionConstantPower=0;
    trh.TransductionUnit=0;
    trh.TraceIdentifier=0;
    trh.ScalarTraceHeader=0;

    trh.SourceType=0;
    trh.SourceEnergyDirectionMantissa=0;
    trh.SourceEnergyDirectionExponent=0;
    trh.SourceMeasurementMantissa=0;
    trh.SourceMeasurementExponent=0;
    trh.SourceMeasurementUnit=0;

    trh.UnassignedInt1=0;
    trh.UnassignedInt2=0;
end


function WriteTrace(segyfh,tracedata,trh,h)

    TraceStart=ftell(segyfh);

    fseek(segyfh,TraceStart,'bof');
    fwrite(segyfh,trh.TraceSequenceLine,'int32');    % 0
    fwrite(segyfh,trh.TraceSequenceFile,'int32');    % 4
    fwrite(segyfh,trh.FieldRecord,'int32');          % 8
    fwrite(segyfh,trh.TraceNumber,'int32');          % 12
    fwrite(segyfh,trh.EnergySourcePoint,'int32');    % 16
    fwrite(segyfh,trh.cdp,'int32');                  % 20
    fwrite(segyfh,trh.cdpTrace,'int32');             % 24
    fwrite(segyfh,trh.TraceIdentificationCode,'int16'); % 28
    fwrite(segyfh,trh.NSummedTraces,'int16'); % 30
    fwrite(segyfh,trh.NStackedTraces,'int16'); % 32
    fwrite(segyfh,trh.DataUse,'int16'); % 34
    fwrite(segyfh,trh.offset,'int32');             %36
    fwrite(segyfh,trh.ReceiverGroupElevation,'int32');             %40
    fwrite(segyfh,trh.SourceSurfaceElevation,'int32');             %44
    fwrite(segyfh,trh.SourceDepth,'int32');             %48
    fwrite(segyfh,trh.ReceiverDatumElevation,'int32');             %52
    fwrite(segyfh,trh.SourceDatumElevation,'int32');             %56
    fwrite(segyfh,trh.SourceWaterDepth,'int32');  %60
    fwrite(segyfh,trh.GroupWaterDepth,'int32');  %64
    fwrite(segyfh,trh.ElevationScalar,'int16');  %68

    fwrite(segyfh,trh.SourceGroupScalar,'int16');  %70 % Multiply/divide next number for following 4 values
    fwrite(segyfh,trh.SourceX,'int32');  %72
    fwrite(segyfh,trh.SourceY,'int32');  %76
    fwrite(segyfh,trh.GroupX,'int32');  %80
    fwrite(segyfh,trh.GroupY,'int32');  %84

    fwrite(segyfh,trh.CoordinateUnits,'int16');  %88
    fwrite(segyfh,trh.WeatheringVelocity,'int16');  %90
    fwrite(segyfh,trh.SubWeatheringVelocity,'int16');  %92
    fwrite(segyfh,trh.SourceUpholeTime,'int16');  %94
    fwrite(segyfh,trh.GroupUpholeTime,'int16');  %96
    fwrite(segyfh,trh.SourceStaticCorrection,'int16');  %98
    fwrite(segyfh,trh.GroupStaticCorrection,'int16');  %100
    fwrite(segyfh,trh.TotalStaticApplied,'int16');  %102
    fwrite(segyfh,trh.LagTimeA,'int16');  %104
    fwrite(segyfh,trh.LagTimeB,'int16');  %106
    fwrite(segyfh,trh.DelayRecordingTime,'int16');  %108
    fwrite(segyfh,trh.MuteTimeStart,'int16');  %110
    fwrite(segyfh,trh.MuteTimeEND,'int16');  %112
    fwrite(segyfh,trh.ns,'uint16');  %114
    fwrite(segyfh,trh.dt,'uint16');  %116
    fwrite(segyfh,trh.GainType,'int16');  %118
    fwrite(segyfh,trh.InstrumentGainConstant,'int16');  %120
    fwrite(segyfh,trh.InstrumentInitialGain,'int16');  %%122
    fwrite(segyfh,trh.Correlated,'int16');  %124
    fwrite(segyfh,trh.SweepFrequenceStart,'int16');  %126
    fwrite(segyfh,trh.SweepFrequenceEnd,'int16');  %128
    fwrite(segyfh,trh.SweepLength,'int16');  %130
    fwrite(segyfh,trh.SweepType,'int16');  %132
    fwrite(segyfh,trh.SweepTraceTaperLengthStart,'int16');  %134
    fwrite(segyfh,trh.SweepTraceTaperLengthEnd,'int16');  %136
    fwrite(segyfh,trh.TaperType,'int16');  %138
    fwrite(segyfh,trh.AliasFilterFrequency,'int16');  %140
    fwrite(segyfh,trh.AliasFilterSlope,'int16');  %142
    fwrite(segyfh,trh.NotchFilterFrequency,'int16');  %144
    fwrite(segyfh,trh.NotchFilterSlope,'int16');  %146
    fwrite(segyfh,trh.LowCutFrequency,'int16');  %148
    fwrite(segyfh,trh.HighCutFrequency,'int16');  %150
    fwrite(segyfh,trh.LowCutSlope,'int16');  %152
    fwrite(segyfh,trh.HighCutSlope,'int16');  %154
    fwrite(segyfh,trh.YearDataRecorded,'int16');  %156
    fwrite(segyfh,trh.DayOfYear,'int16');  %158
    fwrite(segyfh,trh.HourOfDay,'int16');  %160
    fwrite(segyfh,trh.MinuteOfHour,'int16');  %162
    fwrite(segyfh,trh.SecondOfMinute,'int16');  %164
    fwrite(segyfh,trh.TimeBaseCode,'int16');  %166
    fwrite(segyfh,trh.TraceWeightningFactor,'int16');  %170
    fwrite(segyfh,trh.GeophoneGroupNumberRoll1,'int16');  %172
    fwrite(segyfh,trh.GeophoneGroupNumberFirstTraceOrigField,'int16');  %174
    fwrite(segyfh,trh.GeophoneGroupNumberLastTraceOrigField,'int16');  %176
    fwrite(segyfh,trh.GapSize,'int16');  %178
    fwrite(segyfh,trh.OverTravel,'int16');  %178
    fwrite(segyfh,trh.cdpX,'int32');  %180
    fwrite(segyfh,trh.cdpY,'int32');  %184
    fwrite(segyfh,trh.Inline3D,'int32');  %188
    fwrite(segyfh,trh.Crossline3D,'int32');  %192
    fwrite(segyfh,trh.ShotPoint,'int32');  %196
    fwrite(segyfh,trh.ShotPointScalar,'int16');  %200
    fwrite(segyfh,trh.TraceValueMeasurementUnit,'int16');  %202
    fwrite(segyfh,trh.TransductionConstantMantissa,'int32');  %204
    fwrite(segyfh,trh.TransductionConstantPower,'int16'); %208
    fwrite(segyfh,trh.TransductionUnit,'int16');  %210
    fwrite(segyfh,trh.TraceIdentifier,'int16');  %212
    fwrite(segyfh,trh.ScalarTraceHeader,'int16');  %214
    fwrite(segyfh,trh.SourceType,'int16');  %216
    fwrite(segyfh,trh.SourceEnergyDirectionMantissa,'int32');  %218
    fwrite(segyfh,trh.SourceEnergyDirectionExponent,'int16');  %222
    fwrite(segyfh,trh.SourceMeasurementMantissa,'int32');  %224
    fwrite(segyfh,trh.SourceMeasurementExponent,'int16');  %228
    fwrite(segyfh,trh.SourceMeasurementUnit,'int16');  %230
    % WRITE UNASSIGNED CHARACTERS FOR THE REST
    fwrite(segyfh,trh.UnassignedInt1,'int32');  %232
    fwrite(segyfh,trh.UnassignedInt2,'int32');  %236
    % 217-240 Unassigned

    fseek(segyfh,TraceStart+240,'bof');
    % write trace data:
    fwrite(segyfh,tracedata,'float32');

end