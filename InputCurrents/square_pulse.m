function I = square_pulse(Iamp, width, IP, NC, len, PWindow)
% square_pulse -> creates a train of square pulses
% Inputs:
% Iamp: [A] amplitude of the pulses
% stim_onset: [s] stimulus onset
% width: width of the square pulse
% IP_: Inter-Pulse Interval
% NC: Number of cycles
% PWindow = [P0,PF]
% P0 starting of the Train (points)
% PF ending of the Train (points)

if NC*(width + IP)>len

    I = Iamp.*repmat([zeros(1,IP) ones(1,width)], 1, NC);
    I = I(1:len);
else
    I = Iamp.*repmat([zeros(1,IP) ones(1,width)], 1, NC);
    I = cat(2, I, zeros(1,(len-NC*(width + IP))));
end

I(1:PWindow(1))=zeros(1,PWindow(1));
I(PWindow(2):end)=zeros(1,length(I)-PWindow(2)+1);
end