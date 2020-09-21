%%this script is to take all laser position data from where it is and
%%writing it as a NPY file to the subjects folder this works for the file
%%format from the first laser scanning experiments through to 2019-XX-XX
clear all
if ispc
    baseDir= "F:\laserPosititionData";
    saveBase = "F:\Subjects";
else
    baseDir = "/Users/ckrasnia/Desktop/matlab_scanData";
    saveBase = "/Users/ckrasnia/Desktop/Zador_Lab/scanData/Subjects/";
end
cd(baseDir)
subs = dir;
for sub = 3:length(subs)
    cd(baseDir)
    subject = subs(sub).name;
    if strcmp(subject,'.DS_Store')
        continue
    end
    cd(fullfile(baseDir,subject))
    days = dir;
    for day = 4:length(days) 
        if days(day).isdir
            if strcmp(days(day).name,'.DS_Store')
                day=day+1;
            end
            exptDay = days(day).name;
            cd(fullfile(baseDir,subject,exptDay))
            files = dir;
        end
        for file = 3:length(files)
            if files(file).name(end-2:end) == 'mat'
                load(files(file).name)
                fileparts = split(files(file).name,"_");suffix = split(fileparts(end),".");
                exptNum = strcat("00",suffix(1));
                f1 = fullfile(saveBase,subject,exptDay,exptNum);
                try
                    cd(f1)
                catch
                    sprintf('Warning, no folder %s, creating folder', f1)
                    mkdir(fullfile(saveBase,subject,exptDay),exptNum);
                    cd(f1)
                end
                if isfield(data_struct, 'laserLoc')                   
                    tab = array2table(data_struct.laserLoc);
                    writetable(tab,"laserLoc.csv")      
                end
                cd(fullfile(baseDir,subject,exptDay))
            end
        end
        
    end
end