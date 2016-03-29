function newci = reorder_membership(ci)

n = length(ci);
group_map = -ones(1,n);

curnum=0;

for i=1:n
    if group_map(ci(i))<0
        group_map(ci(i)) = curnum;
        curnum = curnum + 1;
    end
end

newci = -ones(1,n);

for i=1:n
   val = group_map(ci(i));
   if val<0
       val = group_map(end);
   end
   newci(i)=val;
end

