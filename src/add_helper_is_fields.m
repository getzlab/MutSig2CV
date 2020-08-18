function x = add_helper_is_fields(x)

if ~all(isfield(x,{'type','classification','ref_allele','newbase'})), x = add_simple_fieldnames(x); end

% ANNOTATE FLANK, CODING, SILENT, INDEL

x.is_coding = false(slength(x),1);
x.is_coding(grepi('indel|read.?thr|de.?novo|start|stop|missense|nonsense|silent|splice|synon|frame|shift|non.?stop',...
                  x.type,1)) = true;
x.is_coding(grepi('IGR|flank|UTR|intron|promoter|unknown|non.?coding|mismatch|^RNA',x.type,1)) = false;

x.is_flank = ~x.is_coding;

x.is_indel = false(slength(x),1);
x.is_indel(grepi('INS|DEL|INDEL',x.classification,1)) = true;
x.is_indel(grepi('Ins|Del|Frame|Shift|Indel',x.type,1)) = true;
x.is_indel(strcmp('-',x.ref_allele)) = true;
x.is_indel(strcmp('-',x.newbase)) = true;
x.is_ins = x.is_indel & (grepmi('ins',x.classification)|grepmi('ins',x.type)|strcmp('-',x.ref_allele));
x.is_del = x.is_indel & (grepmi('del',x.classification)|grepmi('del',x.type)|strcmp('-',x.newbase));

x.is_missense = false(slength(x),1); x.is_missense(grepi('missense',x.type,1)) = true;
x.is_nonsense = false(slength(x),1); x.is_nonsense(grepi('nonsense',x.type,1)) = true;
x.is_splice = false(slength(x),1); x.is_splice(grepi('splice',x.type,1)) = true;

x.is_silent = false(slength(x),1); x.is_silent(grepi('silent|synon',x.type,1)) = true;

