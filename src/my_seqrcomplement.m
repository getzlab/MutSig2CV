function rc = my_seqrcomplement(seq)

if isempty(seq)
  rc = seq;
else
  if ischar(seq)
    rc = tr(fliplr(seq),'ACGTacgt','TGCAtgca');
  else
    for a=1:size(seq,1)
      for b=1:size(seq,2)
        for c=1:size(seq,3)
          for d=1:size(seq,4)
            for e=1:size(seq,5)
              rc{a,b,c,d,e} = tr(fliplr(seq{a,b,c,d,e}),'ACGTacgt','TGCAtgca');
            end,end,end,end,end
  end
end
