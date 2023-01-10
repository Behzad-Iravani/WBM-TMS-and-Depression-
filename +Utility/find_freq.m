function Estimated_freq = find_freq(TCs)
% This function estimates the frequency with maximum power

for node=1:size(TCs,2)
    clear x f Pxx loc
    x       = TCs(:,node)';
    Nfft    = 128;
    fsamp   = 1/2;
    [Pxx,f] = pwelch(x,gausswin(Nfft),Nfft/2,Nfft,fsamp);
    [~,loc] = max(Pxx);
    Estimated_freq(node) = f(loc);
end
end