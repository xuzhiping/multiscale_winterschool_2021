function [rnnew,linksnew]=cleanupnodes(rn,links);
%this subroutine begins removes all of the bad nodes in the system 
% and cleans up the data structures to include good nodes only.
%Note: this function fails to clean up (-1) nodes surrounded by valid nodes
%       it only removes a continuous chunk of (-1) nodes in the end
%      it would be good if it can remove all (-1) nodes
rnnew=rn;
linksnew=links;
e=1;
while e<2
   lrn=length(rnnew(:,4));
   if lrn>0
        if rnnew(lrn,4)==-1
            rnnew(lrn,:)=[];
        else
            e=2;
        end
    else
        e=2;
    end
end
e=1;
while e<2
    llinks=length(linksnew(:,1));
    if llinks>0
        if linksnew(llinks,1)==0
            linksnew(llinks,:)=[];
        else
            e=2;
        end
    else
        e=2;
    end
end