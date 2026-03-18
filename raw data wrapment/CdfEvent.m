classdef CdfEvent < handle
    % CdfEvent contains CDF event data

    properties(SetAccess=?CdfDataFile)
        runNumber    % integer:  CDF run number
        eventNumber  % integer:  CDF event number within run
        vertex       % 2x1 float:  transverse location (xy) of primary interaction vertex
        tracks       % array of CdfTrack:  track parameters
    end

    methods

        % constructor:  create CdfEvent in invalid state
        function obj = CdfEvent()
            obj.runNumber = -1;
            obj.eventNumber = -1;
            obj.vertex = [0; 0];
            obj.tracks = [];
        end

        % test if this event is valid
        function valid = isValid(obj)
            valid = (obj.runNumber >= 0);
        end

    end % methods

end % classdef
