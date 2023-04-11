function ecog_prf_plotScatter(opts,prf_params1,prf_params2)

% Description: 
%
% ecog_prf_plotScatter(opts,prf_params1,prf_params2)
%
% Input
% - modeldata          	= Nx1 cell-array of time-series structure with the following fields:
%   - subject
%   - channels
%   - datats            = Mx1 cell-array of time-series data
%   - stimulus          = Mx1 cell-array of stimulus data
% - opts
%   - plotsavedir
%   - plot
%   ---------------------
%   - 
%
% Example
%   opts = [];
%   opts.plotparam      = {'ecc','rfsize'};
%   opts.plot.legend    = {'eccentricity', 'pRF size'};
%   opts.plot.XLabel    = 'broadband';
%   opts.plot.YLabel    = 'alpha';

%   ecog_prf_plotScatter(opts,prf_params1,prf_params2)

% Dependency: ECoG_utils, SetDefault

% 20200324 - Yuasa

narginchk(3,3);

SetDefault('opts.plotparam',{'R2','ecc','ang','rfsize'});
SetDefault('opts.plot.XLim',[]);
SetDefault('opts.plot.YLim',[]);
SetDefault('opts.plot.XLabel','');
SetDefault('opts.plot.YLabel','');
SetDefault('opts.plot.legend','');

lgndflg = ~isempty(opts.plot.legend);

ndats = length(prf_params1);

[ix,iy] = arrangeinrect(ndats,1.5,[1.8 1]);
for iparam = 1:length(opts.plotparam)
  figure('Name',opts.plotparam{iparam},'Position',[150 100 1000 800]);
  for ii = 1:ndats
      subplot_er(iy,ix,ii);
      dat1 = prf_params1{ii}.xval;
      dat2 = prf_params2{ii}.xval;
      scatter(dat1,dat2);
      limmin = floor(min([0;dat1;dat2]));
      limmax = ceil(max([0;dat1;dat2]));
      hold on;
      selelec = prf_params1{ii}.xval>threshbb;
      scatter(dat1(selelec),dat2(selelec));
      hold off;
      set(gca,'FontSize',18);
      xlim([limmin limmax]);      ylim([limmin limmax]);
      straightline(0,'h','k:');   straightline(0,'v','k:');
      xlabel('xR^2 (broadband)');    ylabel('xR^2 (alpha)');
      title(subjectList{ii});
  end
end
