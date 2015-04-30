function pdf=pdfscaling(x,tlag)
% pdf=pdfscaling(x,tlag)
% This function is to calculate the pdf scaling of velocity increments
% Input
% x is the time series to be analyzed
% tlag is the maximum separation scale
% Output
% pdf is the pdf scaling
%      pdf.M is the maxima pdf of increments
%      pdf.Z is the zero-crossing pdf
%      pdf.tau is the separation scales
% See also: pdfscalingc
% 
% Written by Yongxiang HUANG 28/03/2010
% 

%   References:
%   HUANG Y., ZHOU Q. QIU X. SCHMITT F.G., SHANG X., LU Z. LIU Y. Scaling 
%   Scaling of probability density functions of velocity increments in turbulent system
%   Physics of Fluid (submitted)




if min(size(x))>1
    error('x should be one dimensional time series!');
end

if nargin==1
    tlag=fix(length(x)/4);
end

if length(tlag)==1
    tlag=[2 4 6 8 fix(10.^[1:.1:log10(tlag)])];
end

x=x/std(x);


pd=pdfscalingc(x,tlag);
pdf.M=pd(1,:);
pdf.Z=pd(2,:);
pdf.tau=tlag;