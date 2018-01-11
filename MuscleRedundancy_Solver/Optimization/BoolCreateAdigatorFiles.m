function [BoolCreate] = BoolCreateAdigatorFiles(filenames,ContName,EndPointName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

BoolCreate=0;
bool_exist=1;
for i=1:length(filenames)
    if ~exist(filenames{i},'file');
        bool_exist=0;
    end
end
if bool_exist==0
    % check if the adigator files exist
    BoolCreate=1;
else    
    % check if continuous function or endpointfunction is changed since the
    % aditator files are created.
    filename=filenames{1};
    AdigatorFile=dir(filename);
    Time_Adigator=AdigatorFile.datenum;  
    ContFile=dir(ContName);
    if isempty(ContFile)
        ContName=which(ContName);   % when the path is not specified by the user
        ContFile=dir(ContName);
    end
    Time_Cont=ContFile.datenum;
    EndPointFile=dir(EndPointName);
    if isempty(EndPointFile)
        EndPointName=which(EndPointName);
        EndPointFile=dir(EndPointName);
    end
    Time_Endpoint=EndPointFile.datenum;    
    if Time_Adigator<Time_Cont || Time_Adigator<Time_Endpoint
        BoolCreate=1;
    end    
end   
end

