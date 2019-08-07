function newTrialCheckTargetted(src,event)
global s2

       if any(event.Data(:,1) >= 4)
            s2.startBackground
%             saved = event.Data;
            src.stop()
       end
end