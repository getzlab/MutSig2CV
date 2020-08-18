function names = find_good_names_for_mutation_categories(autonames)

if ~iscell(autonames)
  noncellflag = true;
  autonames = {autonames};
else
  noncellflag = false;
end

z.code = {'At','Af','As','Atf','Afs','Ats','Atfs',...
          'Ct','Cf','Cs','Ctf','Cfs','Cts','Ctfs',...
          'ACt','ACf','ACs','ACtf','ACfs','ACts','ACtfs'};
z.name = {'A->G','A->T','A->C','A->(G/T)','A->(T/C)','A->(C/G)','A->mut',...
          'C->T','C->G','C->A','C->(T/G)','C->(G/A)','C->(A/T)','C->mut',...
          'N->transit','N->flip','N->skew','N->nonskew','N->transver','N->nonflip','N->mut'};

p = parse(autonames,'^([ACGT]+)\[([AC])+->([tfs]+)\]([ACGT]+)$',{'before','at','change','after'});
p.muttype = mapacross(stringsplice([p.at p.change]),z.code,z.name);

for i=1:length(autonames)
  names{i,1} = regexprep(p.muttype{i},'^(.*)(->.*)$',[p.before{i} 'p*$1p' p.after{i} '$2']);
end

names = regexprep(names,'ACGTp','');
names = regexprep(names,'pACGT','');
names = regexprep(names,'^([ACGT])([ACGT])p','($1/$2)p');
names = regexprep(names,'^([ACGT])([ACGT])([ACGT])p','($1/$2/$3)p');
names = regexprep(names,'p([ACGT])([ACGT])->','p($1/$2)->');
names = regexprep(names,'p([ACGT])([ACGT])([ACGT])->','p($1/$2/$3)->');
names = regexprep(names,'^\*([ACN])->','$1->');
names = regexprep(names,'^N->','');

if ~noncellflag
  autonames = autonames{1};
end
