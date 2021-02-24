classdef SystemSetUp
%% Setting system environment (Static methods)
    methods(Static, Hidden=true)
        function output = ROOTPATH()            
            [~, computer_name] = system('hostname');
            switch computer_name(1:end-1)
                case 'yuncong-Precision-WorkStation-T7500'
                    output = '/data/Vessel';
                case 'xiang-ThinkPad-T450s'
                    output = '/home/xiang/Data';
                case 'dklab-workstation'
                    output = '/data/Vessel';
                case 'XiangJi-PC'
                    output = 'D:\Data\Vessel';
            end 
        end
        
        function output = ExtLibRootPath()
            [~, computer_name] = system('hostname');
            switch computer_name(1:end-1)
                case 'dklab-workstation'
                    output = '/home/dklab/Documents/Github';
                case 'XiangJi-PC'
                    output = '/home/xiang/Documents/Github';
            end 
        end
        
        function output = Scratch_Folder_Path()
            [~, computer_name] = system('hostname');
            switch computer_name(1:end-1)
                case 'dklab-workstation'
                    output = '/scratch/Vessel';
                case 'XiangJi-PC'
                    output = '/tmp/scratch';
            end
        end
        
        function output = AnalysisResultFolder()
            [~, computer_name] = system('hostname');
            switch computer_name(1:end-1)
                case 'dklab-workstation'
                    output = '/home/dklab/GoogleDrive/Vessel/Analysis';
            end
        end
        
    end
    

    
end
