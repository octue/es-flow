classdef DiskVar < handle
    %DISKVAR Dummy variable allowing seamless interfacing to 'big' data stored
    %on the hard drive
    
    properties (Access = protected)
        
        Class
        Size
        MatFilesHdl
        Bytes
        FieldRef
        
    end
    
    
    methods
        
        function [dv] = DiskVar(mfhdl, fieldRef, size, class, bytes)
            dv.Class = class;
            dv.Size = size;
            dv.Bytes = bytes;
            dv.MatFilesHdl = mfhdl;
            dv.FieldRef = fieldRef;
        end
        
        function B = subsref(dv, S)
            
            % TODO CHECKS ON SIZE!
            
            % Get data from the database and return it as a numeric array 
            B = ReadData(dv.MatFilesHdl, dv.FieldRef, S);
            
            % Check that the class of the outgoing variable matches the Class
            % stored in the dummy variable
            if ~strcmp(class(B), dv.Class)
                error('MatFiles:IncorrectDataType','Incorrect data type from database. Mismatched fieldname?')
            end
            
        end
        
        function sz = size(dv,varargin)
            if nargin == 1
                sz = dv.Size;
            else
                sz = dv.Size(varargin{1});
            end
        end
        
        function tf = isa(dv,clsnm)
            
            if strcmpi(dv.Class,clsnm)
                tf = true;
            else
                tf = false;
            end
        end    
        
    end
    
end

