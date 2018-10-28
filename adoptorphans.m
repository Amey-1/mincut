function adoptorphans(orphanset)
%adoptorphans reclaim tree structure
global tree;
global ncols;
global nrows;
global A;
global parent;
global TIME;
global ts;
global dist;

while ~isempty(orphanset)
    p = orphanset(1,:);
    pr = p(1,1);pc = p(1,2);
    orphanset = orphanset(2:end,:);
    par = findparent(pr,pc);
    if isempty(par)
        %scan all neighbours q of p
        if (pr>1)&&(tree(pr-1,pc) == tree(pr,pc) )
            qr = pr-1;qc=pc;
            orphanset = supportfunc1(pr,pc,qr,qc,1,4,orphanset);
        end
        if (pc<ncols)&&(tree(pr,pc+1) == tree(pr,pc) )
            qr = pr;qc=pc+1;
            orphanset = supportfunc1(pr,pc,qr,qc,2,3,orphanset);
        end
        if (pc>1)&&(tree(pr,pc-1) == tree(pr,pc) )
            qr = pr;qc=pc-1;
            orphanset = supportfunc1(pr,pc,qr,qc,3,2,orphanset);
        end
        if (pr<nrows)&&(tree(pr+1,pc) == tree(pr,pc) )
            qr = pr+1;qc=pc;
            orphanset = supportfunc1(pr,pc,qr,qc,4,1,orphanset);
        end
        
        %set tree of p as null and remove p from active set A
        tree(pr,pc)=2;
        A = A(A~=m2laddr(pr,pc,ncols));
    
    else
        %suitable parent found
        parent(pr,pc) = par;
        [parr,parc] = l2maddr(par,ncols);
        dist(pr,pc) = dist(parr,parc)+1;
        ts(pr,pc) = TIME;
    end
end
end

function par = findparent(pr,pc)
%findparent tries to find suitable parent 'par' for node (pr,pc)
global tree;
global G;
global TIME;
global ts;
global dist;
global ncols;
global nrows;

%check for parent q among neighbours of p
if pr>1 % check UP
    qr = pr-1;qc=pc;
    if supportfunc2(pr,pc,qr,qc,1,4)
        par = m2laddr(qr,qc,ncols);
        return;
    end
end
if pc<ncols %check RIGHT
    qr = pr;qc=pc+1;
    if supportfunc2(pr,pc,qr,qc,2,3)
        par = m2laddr(qr,qc,ncols);
        return;
    end
end
if pc>1 %check LEFT
    qr = pr;qc=pc-1;
    if supportfunc2(pr,pc,qr,qc,3,2)
        par = m2laddr(qr,qc,ncols);
        return;
    end
end
if pr<nrows %check DOWN
    qr = pr+1;qc=pc;
    if supportfunc2(pr,pc,qr,qc,2,3)
        par = m2laddr(qr,qc,ncols);
        return;
    end
end
par = [];
end

function truthvalue = validorigin(qr,qc)
%validorigin checks to see if origin of (qr,qc) is either source or sink
global ts;
global TIME;
global parent;
global ncols;

truthvalue = false;
%trace path to terminal
count=1;
curr = qr;curc=qc;
while true %not source 1 or sink 0
    par = parent(curr,curc,1);
    if par==-2
        truthvalue = false;
        return;
    elseif par<0 %source or sink
        updatedistinfo(qr,qc,count);
        truthvalue = true;
        return;
    else %par is some valid parent
        [parr,parc]=l2maddr(par,ncols);
        if ts(parr,parc)==TIME
            %currently checked node found
            updatedistinfo(qr,qc,count,dist(parr,parc));
        else
            count = count+1;
            curr=parr;curc=parc;
        end
    end
end
end

function updatedistinfo(qr,qc,count,dist_r)
%update the distances along path from q to terminal
global ts;
global TIME;
global dist;
global parent;

if nargin<4
    dist_r = 0;
end

curr=qr;curc=qc;
while count>0
    dist(curr,curc)=dist_r+count;
    ts(curr,curc)=TIME;
    count = count-1;
    [curr,curc] = l2maddr(parent(curr,curc));
end
end

function orphanset = supportfunc1(pr,pc,qr,qc,dirp2q,dirq2p,orphanset)
global ncols;
global A;
global parent;
global G;
global tree;

if (tree(qr,qc)==1) && (G(qr,qc,dirq2p)>0)
    A(end+1) = m2laddr(qr,qc,ncols); %add q to active set A
elseif (tree(qr,qc)==0) && (G(pr,pc,dirp2q)>0)
    A(end+1) = m2laddr(qr,qc,ncols); %add q to active set A
end
if (parent(qr,qc,1)==m2laddr(pr,pc,ncols))
    %add q to set of orphans and set parent is null
    orphanset(end+1,1)=qr;orphanset(end,2)=qc;
    parent(qr,qc,1) = -2;
end
end

function truthvalue = supportfunc2(pr,pc,qr,qc,dirp2q,dirq2p)
%supportunc2 tries to find suitable parent in specified direction
global tree;
global G;
truthvalue = false;
if (tree(qr,qc)==tree(pr,pc)) && (tree(qr,qc)==1) && (G(qr,qc,dirq2p)>0) && validorigin(qr,qc)
    truthvalue = true;
elseif (tree(qr,qc)==tree(pr,pc)) && (tree(qr,qc)==2) && (G(pr,pc,dirp2q)>0) && validorigin(qr,qc)
    truthvalue = true;
end
end