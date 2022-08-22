function [u]= extract_mesg(c,rearranged_cols)
%u= extract_mesg(c,rearranged_cols)
%For examples and more details, please refer to the LDPC toolkit tutorial at
%http://arun-10.tripod.com/ldpc/ldpc.htm 
rows=length(rearranged_cols);
for i=1:rows
   if rearranged_cols(i)~=0
      temp=c(:,i);
      c(:,i)=c(:,rearranged_cols(i));   %size c [K,cols]
      c(:,rearranged_cols(i))=temp;
   end
end
cols=size(c,2);
u=c(:,rows+1:cols);   %size c [K,cols-rows]
