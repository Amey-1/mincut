function O = augmentpath(augpath,bottleneck)
%augmentpath Push 'bottleneck' flow through 'augpath' from s to t to saturate edges. Return orphans created
global parent;
global tree;
global G;
global source2px;
global px2sink;
O = []; %n rows 2 cols(r,c)
%O = zeros(size(augpath,1),3);WILL NEED SEPARATE ROW INDEX COUNTER SINCE end CANNOT BE USED %n rows 2 cols(r,c)

%update edge from source to first node in path
source2px(augpath(1,1), augpath(1,2)) = source2px(augpath(1,1), augpath(1,2)) - bottleneck;
if source2px(augpath(1,1), augpath(1,2)) == 0 %saturated
    parent(augpath(1,1), augpath(1,2),1) = -2; %null
    O(end+1,1)= augpath(1,1); O(end,2) = augpath(1,2);
end

for i = 1:(size(augpath,1)-1)
    %push flow from i to i+1 [augpath third column is how to reach cur node from parent]
    G(augpath(i,1),augpath(i,2),augpath(i+1,3)) = G(augpath(i,1),augpath(i,2),augpath(i+1,3)) - bottleneck;
    
    switch augpath(i+1,3)
        case 1
            G(augpath(i+1,1),augpath(i+1,2),4) = G(augpath(i+1,1),augpath(i+1,2),4) + bottleneck;
        case 2
            G(augpath(i+1,1),augpath(i+1,2),3) = G(augpath(i+1,1),augpath(i+1,2),3) + bottleneck;
        case 3
            G(augpath(i+1,1),augpath(i+1,2),2) = G(augpath(i+1,1),augpath(i+1,2),2) + bottleneck;
        case 4
            G(augpath(i+1,1),augpath(i+1,2),1) = G(augpath(i+1,1),augpath(i+1,2),1) + bottleneck;
    end
    
    if G(augpath(i,1),augpath(i,2),augpath(i+1,3)) == 0 %saturated
        if (tree(augpath(i,1),augpath(i,2))==tree(augpath(i+1,1),augpath(i+1,2))) && (tree(augpath(i,1),augpath(i,2))== 1) %source
            parent(augpath(i+1,1), augpath(i+1,2),1) = -2; %null
            O(end+1,1)= augpath(i+1,1); O(end,2) = augpath(i+1,2);
        elseif (tree(augpath(i,1),augpath(i,2))==tree(augpath(i+1,1),augpath(i+1,2))) && (tree(augpath(i,1),augpath(i,2))== 0) %sink
            parent(augpath(i,1), augpath(i,2),1) = -2; %null
            O(end+1,1)= augpath(i,1); O(end,2) = augpath(i,2);
        end
    end
end

%update edge from last node in path to sink
px2sink(augpath(end,1),augpath(end,2)) = px2sink(augpath(end,1),augpath(end,2)) - bottleneck;
if px2sink(augpath(end,1),augpath(end,2))==0 %saturated
    parent(augpath(end,1), augpath(end,2)) = -2; %null
    O(end+1,1)= augpath(end,1); O(end,2) = augpath(end,2);
end
end

