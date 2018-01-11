function phaseout = Wrap4musdynContinous_FtildeState_vM1(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time) 
    splinestruct = SplineInputData_vM(input.phase.time,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_FtildeState(input);