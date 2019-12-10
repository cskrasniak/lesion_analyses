function newTrialCheck(src,event)
global s2

       if any(event.Data >= 5)
            s2.startForeground
%             saved = event.Data;
            src.stop()
       end
end