close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'
             'rxnFluxes'};


pdfnames = {};
for a = 1:length(datatypes)
  f000_predictPerturbs_pairwiseCorr(datatypes{a})
  export_fig(['s002_' num2str(gcf)],'-pdf')
  pdfnames{end+1} = ['s002_' num2str(gcf) '.pdf'];
end

append_pdfs('s002_output.pdf',pdfnames{:});
for a = 1:length(pdfnames)
  delete(pdfnames{a})
end