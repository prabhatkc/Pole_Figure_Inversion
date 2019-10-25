function [] = prompt(init)
    if strcmp(init, 'F')
        for k =1:1
            close all;
            continue;
        end
    else 
    	what = 'enter [r]un to continue or [q]uit to exit out: ';
        str = input(what, 's');
        for k = 1:1
            if strcmp(str(1), 'r')
                close all;
                continue;
            elseif strcmp(str(1), 'q')
                close all;
                error('quitting this program ...');
            end 
        end
    end
                                   
end
