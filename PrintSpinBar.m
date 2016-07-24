%
%
% Prints a spinning bar
%


function PrintSpinBar(i, IMax)

spin = ['-','-','-','\\','\\','\\','|','|','|','\/','\/','\/'];


if(i<IMax)
    if(mod(i,13))
       fprintf('\b\b\b[%s]',spin(mod(i,13)))
       drawnow
    end

else
   
fprintf('\b\b\b\b..done\n')
end


    
    