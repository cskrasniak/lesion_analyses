%%this script is to take all laser position data from where it is and
%%writing it as a NPY file to the subjects folder this works for the file
%%format from the first laser scanning experiments through to 2019-XX-XX

baseDir= "/Users/ckrasnia/Desktop/matlab_scanData";
saveBase = "/Users/ckrasnia/Desktop/Zador_Lab/scanData/Subjects/";
cd(baseDir)
subs = dir;
for sub = 4:length(subs)
    cd(baseDir)
    subject = subs(sub).name;
    cd(fullfile(baseDir,subject))
    days = dir;
    for day = 3:length(days) 
        if days(day).isdir
            exptDay = days(day).name;
            cd(fullfile(baseDir,subject,exptDay))
            files = dir;
        end
        for file = 3:length(files)
            if strcmp(files(file).name,'.DS_Store')
                continue
            end
            load(files(file).name)
            fileparts = split(files(file).name,"_");suffix = split(fileparts(end),".");
            exptNum = strcat("00",suffix(1));
            f1 = fullfile(saveBase,subject,exptDay,exptNum);
            try
                cd(f1)
            catch
                sprintf('Warning, no folder %s, creating folder', f1)
                mkdir(fullfile(saveBase,subject,exptDay),exptNum)
                cd(f1)
            end
            writeNPY(data_struct.XY_list,"laserData")
            cd(fullfile(baseDir,subject,exptDay))
        end
        
    end
end