function stillRunningCheck(src,event)
global s2
    if any(event.Data(:,1) >= 4) && s2.IsLogging
        s2.stop()
        src.stop()
    end   
end