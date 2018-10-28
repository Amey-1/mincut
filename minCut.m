function newlabelling = minCut(img,labelling,prob)
%minCut - calculate new labelling of image based on min cut
    global nrows;
    nrows = size(img,1);
    global ncols;
    ncols = size(img,2);
    global G;
    G = zeros(nrows,ncols,4); %edge weights stored here
    %Neighbour order is UP(1) RIGHT(2) LEFT(3) DOWN(4)
    global source2px;
    source2px = zeros(nrows,ncols);
    global px2sink;
    px2sink = zeros(nrows,ncols);
    global TIME;
    TIME = 1;
    
    global dist;
    dist = ones(nrows,ncols);
    global ts;
    ts = ones(nrows,ncols);
    
    beta = rand*15+5;
    
    for row = 1:nrows
        for col = 1:ncols
            z1 = img(row,col,:);z1 = squeeze(z1)';
            %update G(n) pixel edge weights
            if row ~= 1
                z2 = img(row-1,col,:);z2 = squeeze(z2)';
                G(row,col,1) = beta*((norm(z1-z2)).^2);%*beta=1 Max possible value of this equation is 195075
            end
            if col ~= ncols
                z2 = img(row,col+1,:);z2 = squeeze(z2)';
                G(row,col,2) = beta*((norm(z1-z2)).^2);
            end
            if col ~= 1
                z2 = img(row,col-1,:);z2 = squeeze(z2)';
                G(row,col,3) = beta*((norm(z1-z2)).^2);
            end
            if row ~= nrows
                z2 = img(row+1,col,:);z2 = squeeze(z2)';
                G(row,col,4) = beta*((norm(z1-z2)).^2);
            end
            if labelling(row,col)==1
                %foreground
                source2px(row,col) = -log(prob(row,col));
                px2sink(row,col) = 0;
            else
                px2sink(row,col) = -log(prob(row,col));
                source2px(row,col) = 0;
            end
        end
    end
    
    global tree;
    tree = labelling; %0-sink 1-source 2-empty
    global parent;
    parent = -labelling; %(-2)-empty 0-sink (-1)-source
    parent(:,:,2)=0; %direction of current node w.r.t. parent UP(1) RIGHT(2) LEFT(3) DOWN(4) SOURCE/SINK(0)
    global A;
    A = (1:nrows*ncols); %source : -1, sink: 0
  
    %main algo
    count = 1;tic;
    while true
        [augpath,bottleneck] = growthstage(A);
        fprintf('Growth stage %d completed.\n',count);toc;
        if isempty(augpath)
            break
        end
        TIME = TIME+1;
        orphanset = augmentpath(augpath,bottleneck);
        fprintf('Augment stage %d completed.\n',count);toc;
        adoptorphans(orphanset);
        fprintf('Adoption stage %d completed.\n',count);toc;
        count = count+1;
    end
    
    %assign new labels
    newlabelling = tree==1;
end

function [augpath,bottleneck] = growthstage(A)
global parent;
global tree;
global G;
global nrows;
global ncols;
global dist;
global ts;

while ~isempty(A)
    p = A(1);
    [r,c] = l2maddr(p,ncols);
    if tree(r,c)==1 %part of S
        %for every neighbour q
        if r~=1 %UP
            if G(r,c,1)>0
                if tree(r-1,c)==2
                    tree(r-1,c) = 1;
                    parent(r-1,c,1) = p;
                    parent(r-1,c,2) = 1;
                    A(end+1) = (r-2)*ncols + c;
                    dist(r-1,c) = dist(r,c)+1;
                    ts(r-1,c) = ts(r,c);
                else
                    if tree(r-1,c) ~= 1
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r-1,c,1);return;
                    elseif (dist(r,c)<dist(r-1,c))&&(ts(r,c)>=ts(r-1,c)) %reassign parent
                        parent(r-1,c,1) = p;parent(r-1,c,2) = 1;
                        dist(r-1,c) = dist(r,c)+1;
                        ts(r-1,c) = ts(r,c);
                    end
                end
            end
        end
        if c~=ncols %RIGHT
            if G(r,c,2)>0
                if tree(r,c+1)==2
                    tree(r,c+1) = 1;
                    parent(r,c+1,1) = p;
                    parent(r,c+1,2) = 2;
                    A(end+1)=r*ncols + c+1;%right;
                    dist(r,c+1) = dist(r,c)+1;
                    ts(r,c+1) = ts(r,c);
                else
                    if tree(r,c+1) ~= 1
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r,c+1,2);return;
                    elseif (dist(r,c)<dist(r,c+1))&&(ts(r,c)>=ts(r,c+1)) %reassign parent
                        parent(r,c+1,1) = p;parent(r,c+1,2) = 2;
                        dist(r,c+1) = dist(r,c)+1;
                        ts(r,c+1) = ts(r,c);
                    end
                end
            end
        end
        if c~=1 %LEFT
            if G(r,c,3)>0
                if tree(r,c-1)==2
                    tree(r,c-1)=1;
                    parent(r,c-1,1)=p;
                    parent(r,c-1,2)=3;
                    A(end+1)=r*ncols + c-1;%left;
                    dist(r,c-1) = dist(r,c)+1;
                    ts(r,c-1) = ts(r,c);
                else
                    if tree(r,c-1) ~= 1
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r,c-1,3);return;
                    elseif (dist(r,c)<dist(r,c-1))&&(ts(r,c)>=ts(r,c-1)) %reassign parent
                        parent(r,c-1,1) = p;parent(r,c-1,2)=3;
                        dist(r,c-1) = dist(r,c)+1;
                        ts(r,c-1) = ts(r,c);
                    end
                end
            end
        end
        if r~=nrows %DOWN
            if G(r,c,4)>0
                if tree(r+1,c)==2
                    tree(r+1,c)=1;
                    parent(r+1,c,1)=p;
                    parent(r+1,c,2)=4;
                    A(end+1)=(r+1)*ncols + c;%down;
                    dist(r+1,c) = dist(r,c)+1;
                    ts(r+1,c) = ts(r,c);
                else
                    if tree(r+1,c) ~= 1
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r+1,c,4);return;
                    elseif (dist(r,c)<dist(r+1,c))&&(ts(r,c)>=ts(r+1,c)) %reassign parent
                        parent(r+1,c,1) = p;parent(r+1,c,2)=4;
                        dist(r+1,c) = dist(r,c)+1;
                        ts(r+1,c) = ts(r,c);
                    end
                end
            end
        end
    else %part of T
        %for every neighbour q
        if r~=1 %UP
            if G(r-1,c,4)>0
                if tree(r-1,c)==2
                    tree(r-1,c) = 1;
                    parent(r-1,c,1) = p;
                    parent(r-1,c,2) = 1;
                    A(end+1)=(r-2)*ncols + c;%up;
                    dist(r-1,c) = dist(r,c)+1;
                    ts(r-1,c) = ts(r,c);
                else
                    if tree(r-1,c) ~= 0
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r-1,c,1);return;
                    elseif (dist(r,c)<dist(r-1,c))&&(ts(r,c)>=ts(r-1,c)) %reassign parent
                        parent(r-1,c,1) = p;parent(r-1,c,2) = 1;
                        dist(r-1,c) = dist(r,c)+1;
                        ts(r-1,c) = ts(r,c);
                    end
                end
            end
        end
        if c~=ncols %RIGHT
            if G(r,c+1,3)>0
                if tree(r,c+1)==2
                    tree(r,c+1) = 1;
                    parent(r,c+1,1) = p;
                    parent(r,c+1,2) = 2;
                    A(end+1)=r*ncols + c+1;%right;
                    dist(r,c+1) = dist(r,c)+1;
                    ts(r,c+1) = ts(r,c);
                else
                    if tree(r,c+1) ~= 0
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r,c+1,2);return;
                    elseif (dist(r,c)<dist(r,c+1))&&(ts(r,c)>=ts(r,c+1)) %reassign parent
                        parent(r,c+1,1) = p;parent(r,c+1,2) = 2;
                        dist(r,c+1) = dist(r,c)+1;
                        ts(r,c+1) = ts(r,c);
                    end
                end
            end
        end
        if c~=1 %LEFT
            if G(r,c-1,2)>0
                if tree(r,c-1)==2
                    tree(r,c-1)=1;
                    parent(r,c-1,1)=p;
                    parent(r,c-1,2)=3;
                    A(end+1)=r*ncols + c-1;%left;
                    dist(r,c-1) = dist(r,c) + 1;
                    ts(r,c-1) = ts(r,c);
                else
                    if tree(r,c-1) ~= 0
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r,c-1,3);return;
                    elseif (dist(r,c)<dist(r,c-1))&&(ts(r,c)>=ts(r,c-1)) %reassign parent
                        parent(r,c-1,1) = p;parent(r,c-1,2)=3;
                        dist(r,c-1) = dist(r,c)+1;
                        ts(r,c-1) = ts(r,c);
                    end
                end
            end
        end
        if r~=nrows %DOWN
            if G(r+1,c,1)>0
                if tree(r+1,c)==2
                    tree(r+1,c)=1;
                    parent(r+1,c,1)=p;
                    parent(r+1,c,2)=4;
                    A(end+1)=(r+1)*ncols + c;%down;
                    dist(r+1,c) = dist(r,c)+1;
                    ts(r+1,c) = ts(r,c);
                else
                    if tree(r+1,c) ~= 0
                        %find path
                        [augpath,bottleneck]=findpath(r,c,r+1,c,4);return;
                    elseif (dist(r,c)<dist(r+1,c))&&(ts(r,c)<ts(r+1,c)) %reassign parent
                        parent(r+1,c,1) = p;parent(r+1,c,2)=4;
                        dist(r+1,c) = dist(r,c)+1;
                        ts(r+1,c) = ts(r,c);
                    end
                end
            end
        end
    end
    A = A(2:end);
end
augpath = [];
bottleneck = -1;
end

function [augpath,bottleneck] = growthstage_tbd(A)
global ncols;
global tree;
global parent;
global ts;
global dist;

augpath = [];
bottleneck = 0;

while ~isempty(A)
    p = A(1);
    [pr,pc] = l2maddr(p,ncols);
    %FOR every neighbour of p
    if (pr~=1) && (treecap(pr,pc,1,pr-1,pc,4)>0)%examine UP
        qr = pr-1;qc=pc;
        if tree(qr,qc)==2 %TREE(q) is empty
            tree(qr,qc)=tree(pr,pc);
            parent(qr,qc,1)=p;parent(qr,qc,2)=1;
            A(end+1) = m2laddr(qr,qc,ncols);
        elseif tree(qr,qc)==tree(pr,pc) %TREE(q)=TREE(p)
            if (ts(qr,qc)<=ts(pr,pc)) && (dist(qr,qc)>dist(pr,pc)) %isCloser(p,q)
                parent(qr,qc,1) = p; parent(qr,qc,2)=1;
            end
        end
        
    end
    %ENDFOR every neighbour of p
    A = A(2:end);
end
end

function [paths,bottleneck]=findpath(pr,pc,qr,qc,dir)
global G;
global parent;
global source2px;
global px2sink;
global ncols;
global tree;
minval = 200000000;
paths = []; %n rows 3 cols
if tree(pr,pc)==1 %p in S, q in T
    curr = pr;curc=pc;
    while parent(curr,curc,1)~=-1
        parcur = parent(curr,curc,1);
        [parcurr,parcurc] = l2maddr(parcur,ncols);
        paths(end+1,1)=parr;paths(end,2)=parc;paths(end,3)= parent(curr,curc,2);
        minval = min(minval, G(parcurr,parcurc,parent(curr,curc,2)));
        cur = parcur;
        [curr,curc]=l2maddr(cur,ncols); 
    end
    minval = min(minval, source2px(curr,curc));
    paths = flipud(paths);
    paths(end+1,1)=pr; paths(end,2)=pc; paths(end,3)=parent(pr,pc,2);
    
    minval = min(minval,G(pr,pc,dir)); %edge p->q
    paths(end+1,1)=qr; paths(end,2)=qc; paths(end,3)=dir;
    
    curr = qr;curc=qc;
    while parent(curr,curc,1)~=0
        parcur = parent(curr,curc,1);
        [parcurr,parcurc] = l2maddr(parcur,ncols);
        paths(end+1,1)=parcurr; paths(end,2)=parcurc;
        switch parent(curr,curc,2)
            case 1
                minval = min(minval,G(curr,curc,4));paths(end,3)=4;
            case 2
                minval = min(minval,G(curr,curc,3));paths(end,3)=3;
            case 3
                minval = min(minval,G(curr,curc,2));paths(end,3)=2;
            case 4
                minval = min(minval,G(curr,curc,1));paths(end,3)=1;
        end
        curr = parcurr;curc = parcurc;
    end
    minval = min(minval,px2sink(curr,curc));
    
else %p in T, q in S
    curr=qr;curc=qc;
    while parent(curr,curc,1)~=-1
        parcur = parent(curr,curc,1);
        [parcurr,parcurc] = l2maddr(parcur,ncols);
        paths(end+1,1)=parcurr;paths(end,2)=parcurc;paths(end,3)= parent(curr,curc,2);
        minval = min(minval, G(parcurr,parcurc,parent(curr,curc,2)));
        cur = parcur;
        [curr,curc]=l2maddr(cur,ncols); 
    end
    minval = min(minval, source2px(curr,curc));
    paths = flipud(paths);
    paths(end+1,1)=qr; paths(end,2)=qc; paths(end,3)=parent(qr,qc,2);
    
    %edge q->p
    paths(end+1,1)=pr; paths(end,2)=pc;
    switch dir
        case 1
            minval = min(minval,G(qr,qc,4));paths(end,3)=4;
        case 2
            minval = min(minval,G(qr,qc,3));paths(end,3)=3;
        case 3
            minval = min(minval,G(qr,qc,2));paths(end,3)=2;
        case 4
            minval = min(minval,G(qr,qc,1));paths(end,3)=1;
    end
    
    curr = pr;curc=pc;
    while parent(curr,curc,1)~=0
        parcur = parent(curr,curc,1);
        [parcurr,parcurc] = l2maddr(parcur,ncols);
        paths(end+1,1)=parcurr; paths(end,2)=parcurc;
        switch parent(curr,curc,2)
            case 1
                minval = min(minval,G(curr,curc,4));paths(end,3)=4;
            case 2
                minval = min(minval,G(curr,curc,3));paths(end,3)=3;
            case 3
                minval = min(minval,G(curr,curc,2));paths(end,3)=2;
            case 4
                minval = min(minval,G(curr,curc,1));paths(end,3)=1;
        end
        curr = parcurr;curc = parcurc;
    end
    minval = min(minval,px2sink(curr,curc));
end
bottleneck = minval;
end
